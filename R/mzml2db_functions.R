

library(XML)
library(DBI)

#' Convert mzML files to database
#'
#' This function is the core of the mzml2db package, converting mzML files into
#' a relational database specified by `db_engine`. It parses the mzML files one
#' at a time using R's XML library and writes the MS1 and MS2 data out to the
#' database specified in `db_name`. This allows for rapid SQL-based queries of
#' the associated data using R-native syntax and minimal memory.
#'
#' The database constructed consists (currently) of two tables, one named MS1
#' and the other named MS2. The MS1 table has columns for filename, retention
#' time (rt), m/z ratio (mz), and intensity (int). The MS2 table has columns
#' for filename, retention time, precursor m/z (premz), fragment m/z (fragmz),
#' and intensity (if the files have MS2 information in them). I'm planning to
#' add additional tables with file and scan data but haven't gotten to that yet.
#'
#' @param ms_files A list of the files to write into the database, including
#' their paths (although only the basename will be written into the filename
#' column)
#' @param db_engine The database connector to use, e.g. RSQLite::SQLite(),
#' duckdb::duckdb(), or RPostgreSQL::PostgreSQL()
#' @param db_name The name of the database you'd like to write out to
#' @param overwrite_ok Ok to overwrite an existing database with the same name
#' if it exists? Defaults to FALSE for safety reasons
#' @param scan_batch_size Some files are too large to read into memory all at
#' once (looking at you, Thermo Astral data). This argument controls the number
#' of scans that should be read into memory before writing them into the
#' database. If the function begins to hang or consume too much memory, reduce
#' this value.
#'
#' @returns db_name, if successful
#' @export
#'
#' @examples
#' # Use demo files available in the RaMS package
#' library(RaMS)
#' ms_files <- list.files(system.file("extdata", package="RaMS"),
#'                        pattern = "mzML", full.names = TRUE)[2:4]
#' # DuckDB are nice and small, pretty speedy
#' mzml2db(ms_files, db_engine = duckdb::duckdb(), db_name = "minidata.duckdb")
#' # SQLite is well supported but can get to be large (~3x mzML size)
#' mzml2db(ms_files, db_engine = RSQLite::SQLite, db_name = "minidata.sqlite")
#'
#' conn <- dbConnect(duckdb::duckdb(), "minidata.duckdb")
#'
#' # Write raw SQL to extract a chromatogram
#' bet_chrom <- dbGetQuery(conn, "SELECT * FROM MS1 WHERE mz BETWEEN 118.0 AND 118.1")
#' plot(bet_chrom$rt, bet_chrom$int, type="l")
#'
#' # Or calculate a BPC/TIC and use ggplot2 to graph it
#' library(ggplot2)
#' bpc_query <- "SELECT filename, rt, MAX(int) AS int FROM MS1 GROUP BY rt, filename"
#' manual_bpc <- dbGetQuery(conn, bpc_query)
#' ggplot(manual_bpc) + geom_line(aes(x=rt, y=int, group=filename))
#'
#' # Alternatively, use the dbplyr interface and treat it like a tibble
#' library(dplyr)
#' library(dbplyr)
#' conn %>%
#'   tbl("MS1") %>%
#'   filter(between(mz, 118.0, 118.1)) %>%
#'   ggplot() +
#'   geom_line(aes(x=rt, y=int, group=filename))
#'
#' # Remember to clean up afterward!
#' dbDisconnect(conn)
#' unlink("minidata.duckdb")
#' unlink("minidata.sqlite")
mzml2db <- function(ms_files, db_engine=duckdb::duckdb(), db_name, verbosity=NULL,
                    scan_batch_size=10000, sort_by=NULL, ...){
  sax_env <- new.env()

  sax_env$engine <- DBI::dbConnect(db_engine, db_name)
  sax_env$scan_batch_size <- scan_batch_size
  empty_MS1 <- data.frame(filename=character(), scan_idx=numeric(), rt=numeric(), mz=numeric(), int=numeric())
  DBI::dbWriteTable(sax_env$engine, "MS1", empty_MS1, ...)
  empty_MS2 <- data.frame(filename=character(), scan_idx=numeric(), rt=numeric(), premz=numeric(),
                          fragmz=numeric(), int=numeric(), voltage=numeric())
  DBI::dbWriteTable(sax_env$engine, "MS2", empty_MS2, ...)
  on.exit(DBI::dbDisconnect(sax_env$engine))

  sax_env$MS1_scan_data <- list()
  sax_env$MS2_scan_data <- list()
  sax_env$parsing_binary <- FALSE
  if(!is.null(sort_by)){
    if(sort_by%in%c("mz", "rt", "int", "scan_idx")){
      sax_env$sort_by <- sort_by
    } else {
      warning(paste("Sorting by", sort_by, "not supported, ignoring."))
      sax_env$sort_by <- NULL
    }
  }
  if(is.null(verbosity)){
    verbosity <- ifelse(length(ms_files)==1, 2, 1)
  }
  sax_env$verbosity <- verbosity

  sax_handlers <- list(
    startElement = function(name, attrs) startElemParser(name, attrs, sax_env),
    text = function(content) textElemParser(content, sax_env),
    endElement = function(name, attrs) endElemParser(name, attrs, sax_env)
  )

  if(verbosity>0){
    if(length(ms_files)>=2){
      pb <- txtProgressBar(min = 0, max = length(ms_files), style = 3)
    }
    start_time <- Sys.time()
  }
  for(i in seq_along(ms_files)){
    ms_file_i <- ms_files[i]
    sax_env$filename <- basename(ms_file_i)
    XML::xmlEventParse(ms_file_i, handlers = sax_handlers)
    if(verbosity>0 & length(ms_files)>=2){
      setTxtProgressBar(pb, i)
    }
  }
  if(verbosity>0){
    time_total <- round(difftime(Sys.time(), start_time), digits = 2)
    cat("\nTotal time:", time_total, units(time_total), "\n")
  }
  db_name
}
startElemParser <- function(name, attrs, sax_env){
  if(name == "cvParam"){
    if(attrs["name"] == "ms level"){
      sax_env$scan_ms_level <- as.numeric(attrs["value"])
    }
    if(attrs["name"] == "scan start time"){
      sax_env$scan_rt <- as.numeric(attrs["value"])/60
    }
    if(attrs["name"] == "isolation window target m/z"){
      sax_env$scan_premz <- as.numeric(attrs["value"])
    }
    if(attrs["accession"] %in% c("MS:1000523", "MS:1000521")){
      sax_env$bin_type <- ifelse(attrs["accession"] == "MS:1000523", "m/z", "int")
      bin_prec <- sub("-bit float", "", attrs["name"])
      sax_env$bin_prec <- as.numeric(bin_prec) / 8
    }
    if(attrs["accession"] %in% c("MS:1000574", "MS:1000576")){
      sax_env$bin_compr <- switch(attrs["name"], `zlib` = "gzip", `zlib compression` = "gzip", `no compression` = "none", `none` = "none")
    }
    if(attrs["name"]=="collision energy"){
      sax_env$scan_voltage <- as.numeric(attrs["value"])
    }
  }
  if(name == "binary"){
    sax_env$parsing_binary <- TRUE
    sax_env$binary_bits <- character()
  }
  if(name == "spectrum"){
    sax_env$scan_idx <- as.numeric(attrs["index"])
    # print(paste0("Starting scan #", sax_env$scan_idx))
    if(sax_env$scan_idx %% 1000 == 0){
      if(sax_env$verbosity>1){
        print(paste0("Starting scan #", sax_env$scan_idx))
      }
    }
  }
}
textElemParser <- function(content, sax_env){
  if(sax_env$parsing_binary){
    sax_env$binary_bits <- c(sax_env$binary_bits, content)
  }
}
endElemParser <- function(name, attrs, sax_env){
  if(name == "binary"){
    sax_env$parsing_binary <- FALSE
    bin_vals <- paste(sax_env$binary_bits, collapse = "")
    bin_vals <- base64enc::base64decode(bin_vals)
    if(length(bin_vals)==0){
      return(NULL)
    }
    bin_vals <- memDecompress(as.raw(bin_vals), type = sax_env$bin_compr)
    bin_vals <- readBin(bin_vals, what = "double", size = sax_env$bin_prec, n = length(bin_vals) / sax_env$bin_prec)
    if(sax_env$bin_type == "m/z"){
      sax_env$scan_mzs <- bin_vals
    } else {
      sax_env$scan_ints <- bin_vals
    }
  }
  if(name == "spectrum"){
    if(sax_env$scan_ms_level == 1){
      scan_data <- data.frame(scan_idx = sax_env$scan_idx, rt = sax_env$scan_rt,
                              mz = sax_env$scan_mzs, int = sax_env$scan_ints)
      sax_env$MS1_scan_data <- c(sax_env$MS1_scan_data, list(scan_data))
    }
    if(sax_env$scan_ms_level == 2){
      scan_data <- data.frame(scan_idx = sax_env$scan_idx, rt = sax_env$scan_rt,
                              premz = sax_env$scan_premz, fragmz = sax_env$scan_mzs,
                              int = sax_env$scan_ints,
                              voltage = sax_env$scan_voltage)
      sax_env$MS2_scan_data <- c(sax_env$MS2_scan_data, list(scan_data))
    }
    if((length(sax_env$MS1_scan_data) + length(sax_env$MS2_scan_data)) > sax_env$scan_batch_size){
      if(sax_env$verbosity>1){
        print("Writing batch to database")
      }
      if(length(sax_env$MS1_scan_data)>0){
        new_MS1_data <- do.call(what = rbind, args = sax_env$MS1_scan_data)
        if(!is.null(sax_env$sort_by)){
          new_MS1_data <- new_MS1_data[order(new_MS1_data[[sax_env$sort_by]], )]
        }
        new_MS1_data$filename <- sax_env$filename
        DBI::dbWriteTable(sax_env$engine, "MS1", new_MS1_data, append = TRUE)
        sax_env$MS1_scan_data <- list()
      }
      if(length(sax_env$MS2_scan_data)>0){
        new_MS2_data <- do.call(what = rbind, args = sax_env$MS2_scan_data)
        if(!is.null(sax_env$sort_by)){
          if(sax_env$sort_by=="mz"){
            # Assume sorting by mz for MS2 data means premz, allow specifying fragmz explicitly
            new_MS1_data <- new_MS1_data[order(new_MS1_data[["premz"]], )]
          } else {
            new_MS1_data <- new_MS1_data[order(new_MS1_data[[sax_env$sort_by]], )]
          }
        }
        new_MS2_data$filename <- sax_env$filename
        DBI::dbWriteTable(sax_env$engine, "MS2", new_MS2_data, append = TRUE)
        sax_env$MS2_scan_data <- list()
      }
    }
  }
  if(name == "mzML"){
    if(sax_env$verbosity>1){
      print("Writing final data to database")
    }
    if(length(sax_env$MS1_scan_data)>0){
      new_MS1_data <- do.call(what = rbind, args = sax_env$MS1_scan_data)
      if(!is.null(sax_env$sort_by)){
        new_MS1_data <- new_MS1_data[order(new_MS1_data[[sax_env$sort_by]]), ]
      }
      new_MS1_data$filename <- sax_env$filename
      DBI::dbWriteTable(sax_env$engine, "MS1", new_MS1_data, append = TRUE)
      sax_env$MS1_scan_data <- list()
    }
    if(length(sax_env$MS2_scan_data)>0){
      new_MS2_data <- do.call(what = rbind, args = sax_env$MS2_scan_data)
      if(!is.null(sax_env$sort_by)){
        if(sax_env$sort_by=="mz"){
          # Assume sorting by mz for MS2 data means premz, allow specifying fragmz explicitly
          new_MS2_data <- new_MS2_data[order(new_MS2_data[["premz"]]), ]
        } else {
          new_MS2_data <- new_MS2_data[order(new_MS2_data[[sax_env$sort_by]]), ]
        }
      }
      new_MS2_data$filename <- sax_env$filename
      DBI::dbWriteTable(sax_env$engine, "MS2", new_MS2_data, append = TRUE)
      sax_env$MS2_scan_data <- list()
    }
  }
}

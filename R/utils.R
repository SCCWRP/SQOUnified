# Internal, per-directory counter used to name the small RDS display snapshots that
# writelog() persists (see the data = branch below). Keyed by the log directory so each
# log file gets its own 0001, 0002, ... sequence. Lives in the package namespace so it
# survives across the many writelog() calls of a single run; init.log() resets it when it
# first creates a log file.
.sqo_logdata_counter <- new.env(parent = emptyenv())

#' Function to write to a log file
#'
#' @param content Content to write to the logfile - text or a dataframe/tibble-like object
#' @param logfile Path to the file where the content will be written
#' @param filetype 'txt' or 'csv' - the type of file
#' @param append Append to existing file? or overwrite existing contents
#' @param verbose This is what it will be called in parent functions - verbose tells whether to actually write to the logfile or not
#'                This will allow us to write the code in the other functions and files without a bunch of if statements .
#'
#' @examples
#' writelog(df, 'logs/test.csv', filetype = 'csv', append = T)
#' writelog("# --- Test Log Statement --- #", 'logs/log.txt')

#' @export
writelog <- function(content, logfile, filetype = 'text', append = T, verbose = F, prefix = NULL, include.row.names = F, code = NULL, data = NULL, echo.code = TRUE, pageLength = 10) {

  if (!verbose) return(NULL)

  # Only text and csv (for now - ...... don't see any reason why we would support any other type at this point)
  if ( !(filetype %in% c('text','csv')) ) {
    stop("In logging utility function - filetype must be text or csv")
  }

  if (filetype == 'csv') {
    write.csv(content, logfile, row.names = include.row.names)
  } else {
    if (!file.exists(logfile)) {
      stop(paste0("Error in ",  sys.call(), " - File '", logfile, "' not found!"))
    }
    if (!is.null(prefix)) {
      write(paste(prefix, content), logfile, append = append)
    } else {
      write(content, logfile, append = append)
    }

    # Code snippets are written as display-only (eval = FALSE) chunks. The log documents the
    # code that was run for auditing, but knitting the report must NOT re-execute the pipeline
    # (it would be slow and would require all upstream inputs to be present). This is what lets
    # several logs be merged into one report and knit once - see SQOUnified() / append_log_section().
    if (!is.null(code)) {
      if (echo.code) {
        write("```{r eval=FALSE}\n", file = logfile, append = append)
      } else {
        write("```{r eval=FALSE, echo=FALSE}\n", file = logfile, append = append)
      }
      write(code, file = logfile, append = append)
      write("```\n", file = logfile, append = append)
    }
    if (!is.null(data)) {
      # Persist a snapshot of the table and display it via readRDS() so the chunk is
      # self-contained at knit time: it does not depend on the originating variable being in
      # scope. Snapshots are tiny (callers pass an already-subset preview, e.g. head(25)) and
      # live in a _logdata/ subfolder so the visible log directory stays as the intermediate
      # CSVs. The full data remains downloadable via create_download_link().
      datadir <- file.path(dirname(logfile), '_logdata')
      if (!dir.exists(datadir)) dir.create(datadir, recursive = TRUE)

      key <- dirname(logfile)
      n <- get0(key, envir = .sqo_logdata_counter, ifnotfound = 0L) + 1L
      assign(key, n, envir = .sqo_logdata_counter)
      snapshot <- sprintf("logdata_%04d.rds", n)
      saveRDS(data, file.path(datadir, snapshot))

      write("```{r echo=FALSE}\n", file = logfile, append = append)
      write(paste0("datatable(readRDS('./_logdata/", snapshot, "'), options = list(pageLength = ", pageLength, ", autoWidth = TRUE))\n"), file = logfile, append = append)
      write("```\n", file = logfile, append = append)
    }
  }
}



#' Function to write to a log file
#'
#' @param logfile Path to the file where the content will be written
#' @param base.func.name the function name where the call originated - to be written at the top of the log file
#' @param current.time time the function got called
#' @param is.base.func log file shouldnt get initialized if the function gets called within a function that is not the originating one
#'                    Just tells us whether to actually create the log file or not
#' @param verbose This is what it will be called in parent functions - verbose tells whether to actually write to the logfile or not
#'                This will allow us to write the code in the other functions and files without a bunch of if statements .
#'
#' @examples
#' init.log(logfile, sys.call())
#' @export
init.log <- function(
    logfile,
    base.func.name,
    type = 'text',
    current.time = Sys.time(),
    is.base.func = T,
    verbose = F,
    title = 'Log',
    libraries = c('rmarkdown','knitr','DT','dplyr','tidyr','stringr','SQOUnified')
) {

  logfile = path.expand(logfile)
  logdir = dirname(logfile)

  if ( !(type %in% c('RMarkdown','text') )) {
    stop("Logfile type must be RMarkdown or text")
  }

  if (!verbose) return(NULL);

  # Create the parent directories if they don't exist
  if (!dir.exists(logdir)) {
    dir.create(logdir, recursive = TRUE)
  }
  # Create the file if it doesn't exist
  if (!file.exists(logfile)) {
    file.create(logfile)
    if (type == 'text') {
      write(
        paste(
          "----------------------------------------------------------- BEGIN LOG -",
          current.time,
          "----------------------------------------------------------- "
        ),
        logfile
      )
    } else {
      write( paste0("---\ntitle: \"", title ,"\"\noutput: html_document\n---", collapse = ''), file = logfile)
      write('```{r setup, include=FALSE, message=FALSE, warning=FALSE}\n', file = logfile, append = TRUE)
      for (lib in libraries) {
        write(paste0('library(', lib, ')'), file = logfile, append = TRUE)
      }
      # Tolerate a failing chunk instead of aborting the whole render - a single bad table
      # should not prevent the rest of an audit report from being produced.
      write('knitr::opts_chunk$set(error = TRUE)', file = logfile, append = TRUE)
      write('```', file = logfile, append = TRUE)
    }
    # Fresh log file -> restart the writelog() snapshot numbering for this directory.
    assign(dirname(logfile), 0L, envir = .sqo_logdata_counter)
  }

  return(NULL)

}



#' Create download link within an RMarkdown log file
#' @param data data to be downloaded
#' @param logfile the main log file filepath
#' @param filename name of the csv file to be downloaded
#' @param linktext text on the download link
#' @param include.rownames should rownames of the dataframe be included in the CSV file?
#'
#' @export
create_download_link <- function(data, logfile, filename, linktext = 'Download the data', include.row.names = FALSE, verbose = F){
  if(verbose){
    csvfile <- file.path(dirname(logfile), paste0(filename))
    # write.csv() cannot serialize list-columns (e.g. tox 'result'/'result_control' hold
    # per-row vectors of replicate values) and errors with "unimplemented type 'list' in
    # 'EncodeElement'". Collapse any list-columns to a comma-separated string for the CSV
    # artifact only - the returned data frame the caller holds is untouched.
    list.cols <- names(data)[vapply(data, is.list, logical(1))]
    for (col in list.cols) {
      data[[col]] <- vapply(
        data[[col]],
        function(x) paste(format(x, trim = TRUE), collapse = ", "),
        character(1)
      )
    }
    write.csv(data, csvfile, row.names = include.row.names)
    # Trailing blank line after <br> is required: a raw-HTML block runs until a
    # blank line, so without it the next markdown header gets absorbed into the
    # <br> paragraph and rendered literally (e.g. "#### BRI Step 2 ...").
    write(paste0("\n[", linktext, "](./", basename(csvfile), ")\n\n<br>\n"), file = logfile, append = TRUE)
  }
}



#' Merge one line-of-evidence log into a consolidated report (internal)
#'
#' Appends the body of a section's R Markdown log (e.g. the Benthic or Chemistry log) into a
#' master report under a single top-level heading. Because writelog() writes self-contained
#' chunks (tables read from RDS snapshots, code shown but not executed) the merged report can be
#' knit once without re-running any pipeline. The section's own YAML header and setup chunk are
#' dropped (the master supplies those), its headings are demoted one level so the supplied title
#' is the section's top heading, and relative file references (download links and RDS snapshots)
#' are rewritten to point into the section's subfolder so they resolve from the master's location.
#'
#' @param master.logfile path to the consolidated report being assembled
#' @param section.logfile path to the section log to merge in
#' @param subdir the section's subfolder name relative to the master (e.g. 'Benthic')
#' @param title the heading to give the section in the report (e.g. 'Benthic')
#' @param verbose if FALSE, do nothing
#'
#' @export
append_log_section <- function(master.logfile, section.logfile, subdir, title, verbose = F) {
  if (!verbose) return(invisible(NULL))
  if (!file.exists(section.logfile)) return(invisible(NULL))

  lines <- readLines(section.logfile, warn = FALSE)

  # Drop the YAML front matter (--- ... --- at the very top)
  fence.idx <- which(lines == '---')
  if (length(fence.idx) >= 2 && fence.idx[1] == 1) {
    lines <- lines[-(fence.idx[1]:fence.idx[2])]
  }

  out <- character(0)
  in.chunk <- FALSE
  drop.chunk <- FALSE
  for (ln in lines) {
    is.fence <- grepl("^```", ln)

    if (is.fence && !in.chunk) {
      # opening a code chunk
      in.chunk <- TRUE
      drop.chunk <- grepl("^```\\{r setup", ln)   # the section's setup chunk is dropped
      if (!drop.chunk) out <- c(out, ln)
      next
    }
    if (is.fence && in.chunk) {
      # closing a code chunk
      if (!drop.chunk) out <- c(out, ln)
      in.chunk <- FALSE
      drop.chunk <- FALSE
      next
    }
    if (in.chunk) {
      if (drop.chunk) next
      # rewrite relative file refs inside chunks (e.g. readRDS('./_logdata/x.rds'))
      ln <- gsub("(['\"])\\./", paste0("\\1./", subdir, "/"), ln)
      out <- c(out, ln)
      next
    }

    # outside any code chunk:
    if (grepl("^#{1,6}[[:space:]]", ln)) {
      ln <- paste0("#", ln)                 # demote heading one level
      ln <- sub("^#######+", "######", ln)  # markdown supports at most 6 levels
    }
    # rewrite relative download links: [text](./file.csv) -> [text](./subdir/file.csv)
    ln <- gsub("\\]\\(\\./", paste0("](./", subdir, "/"), ln)
    out <- c(out, ln)
  }

  write(paste0("\n\n# ", title, "\n"), file = master.logfile, append = TRUE)
  write(out, file = master.logfile, append = TRUE)
  invisible(NULL)
}



#' Knit an R Markdown log file to HTML (internal)
#'
#' No-op unless both verbose and knitlog are TRUE and the logfile is a .Rmd. Renders to a
#' sibling .html. Centralizes the render step used by the SQO logging functions.
#'
#' @param logfile path to the .Rmd log to render
#' @param verbose only render when TRUE
#' @param knitlog only render when TRUE
#'
#' @export
knit.log <- function(logfile, verbose = F, knitlog = F) {
  if (!(verbose && knitlog)) return(invisible(NULL))

  if (tolower(tools::file_ext(logfile)) != 'rmd') {
    warning("knitlog = TRUE but the logfile is not an R Markdown (.Rmd) file. Skipping knitting.")
    return(invisible(NULL))
  }

  html_file <- sub("\\.Rmd$", ".html", logfile, ignore.case = TRUE)
  print(paste0("Rendering ", logfile, " to ", html_file))
  rmarkdown::render(
    input = logfile,
    output_file = html_file,
    output_format = "html_document",
    quiet = TRUE
  )
  print("Done")
  invisible(html_file)
}



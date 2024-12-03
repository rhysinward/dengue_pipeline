# Install missing packages and load necessary packages
options(repos = c(CRAN = "http://cran.us.r-project.org"))
install.packages("xfun")
xfun::pkg_attach2(required_packages, message = FALSE)


# Universal safe file reading function
safe_read_file_param <- function(fpath, read_fn, ..., required = FALSE) {
  if (!is.null(fpath)) {
    safe_read <- purrr::safely(read_fn)
    read_dat <- safe_read(fpath, ...)

    if (is.null(read_dat$error)) {
      return(read_dat$result)
    } else {
      message(sprintf("\n(ERRO) `%s` is read with error(s) below:", fpath))
      message("####### ERROR(S)")
      print(read_dat$error)
      stop("####### END OF ERROR(S)")
    }
  } else {
    if (required) {
      stop(sprintf("\n (ERRO) `%s` is required but not found. Exiting.", fpath))
    }
  }
}


# Functions for message logging
info_msg <- function(msg, ...) {
  message("[", Sys.time(), "] ", "(INFO) ", msg, ...)
}

warn_msg <- function(msg, ...) {
  warning("[", Sys.time(), "] ", "(WARN) ", msg, ...)
}


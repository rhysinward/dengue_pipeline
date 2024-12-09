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


# Metadata cleaning functions
calcDecimalDate_fromTxt <- function(dateTxt, sep = "/", namedMonths = FALSE, dayFirst = FALSE) {
  els <- strsplit(dateTxt, sep)[[1]]
  if (dayFirst) {
    if (length(els) > 1) {
      els <- els[length(els):1]
    }
  }

  year <- as.integer(els[1])

  if (length(els) == 1) {
    month <- 6 # 7
    day <- 15 # 2
    decDate <- year + 0.5
  } else {
    if (length(els) == 2) {
      if (nchar(els[2]) > 0) {
        if (namedMonths) {
          month <- match(els[2], c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
        } else {
          month <- as.integer(els[2])
        }
        day <- 15
        decDate <- calcDecimalDate(day, month, year)
      } else {
        month <- 6 # 7
        day <- 15 # 2
        decDate <- year + 0.5
      }
    } else {
      if (namedMonths) {
        month <- match(els[2], c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
      } else {
        month <- as.integer(els[2])
      }

      if (nchar(els[3]) > 0) {
        day <- as.integer(els[3])
      } else {
        day <- 15
      }
      decDate <- calcDecimalDate(day, month, year)
    }
  }


  return(decDate)
}

calcDecimalDate_from_yymmdd <- function(dateTxt, sep = "/", ycutoff = 15, defaultMonth = 6, defaultDay = 15) {
  els <- strsplit(dateTxt, sep)[[1]]
  yy <- as.integer(els[1])
  mm <- as.integer(els[2])
  dd <- as.integer(els[3])

  if (!is.finite(yy)) {
    return(-1)
  } else {
    if (yy <= ycutoff) {
      yy <- yy + 2000
    }
    if ((yy > ycutoff) & (yy < 99)) {
      yy <- yy + 1900
    }

    if (!is.finite(mm)) {
      mm <- 0
    }
    if (!is.finite(dd)) {
      dd <- 0
    }
    return(calcDecimalDate(dd, mm, yy, defaultMonth = defaultMonth, defaultDay = defaultDay))
  }
}

getEl <- function(line, sep = ",", ind = -1, final = FALSE, reconstruct = FALSE, ex = -1, fromEnd = FALSE) {
  els <- strsplit(line, sep)[[1]]

  if (ind[1] != -1) {
    if (fromEnd) {
      ind <- length(els) - (ind - 1)
    }
  }

  if (final) {
    return(els[length(els)])
  } else {
    if (reconstruct) {
      if (ex[1] > 0) {
        if (fromEnd) {
          ex <- length(els) - (ex - 1)
        }
        ind <- setdiff((1:length(els)), ex)
      }

      newLine <- els[ind[1]]
      if (length(ind) > 1) {
        for (i in 2:length(ind)) {
          newLine <- paste(newLine, els[ind[i]], sep = sep)
        }
      }
      return(newLine)
    } else {
      if (ind[1] == -1) {
        return(els)
      } else {
        return(els[ind])
      }
    }
  }
}

calcDecimalDate <- function(day, month, year, defaultMonth = 6, defaultDay = 15) {
  cd <- c(0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334)

  if (month == 0) {
    if (defaultMonth >= 1) {
      month <- defaultMonth
    } else {
      month <- ceiling(runif(1) * 12)
    }
  }

  if (day == 0) {
    if (defaultDay >= 1) {
      day <- defaultDay
    } else {
      day <- ceiling(runif(1) * 30)
    }
  }

  dd <- cd[month] + day - 1

  decDate <- year + (dd / 365)

  return(decDate)
}

invertDecimalDate <- function(decDate, formatAsTxt = FALSE, ddmmyy = FALSE) {
  cd <- c(0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365)
  fractD <- cd / 365

  year <- floor(decDate)
  fractYear <- decDate - year
  month <- which(fractD >= fractYear)[1] - 1

  if (month > 0) {
    fractMonth <- fractYear - fractD[month]
    day <- round((fractMonth * 365) + 1)
  } else {
    month <- 1
    day <- 1
  }

  res <- c(day, month, year)

  if (formatAsTxt) {
    if (month < 10) {
      mm <- paste("0", month, sep = "")
    } else {
      mm <- month
    }
    res <- paste(year, mm, day, sep = "-")
  }

  if (ddmmyy) {
    months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
    if (day < 10) {
      dd <- paste("0", day, sep = "")
    } else {
      dd <- day
    }
    res <- paste(dd, months[month], year, sep = "-")
  }
  return(res)
}

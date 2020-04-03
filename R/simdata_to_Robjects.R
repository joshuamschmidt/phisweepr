# R functions for converting simdata objects produced by
# coalescent simulators such as ms, msms, discoal.
# #' lines are roxygen comments for documentation purposes.

#' Read in msms data,
#'
#' @param ms.file A string
#' @param software A string
#' @return A list containing the simulation command 
#' physical positions and genotype sparse matrix from \code{ms.file}
#' @examples
#' get.ms.output(ms.file, software="ms")
get.ms.output <- function(ms.file = NA, software = msms) {
  if (!is.na(ms.file)) {
    raw <- scan(
      file = ms.file,
      what = character(0),
      sep = "\n",
      blank.lines.skip = FALSE,
      quiet = TRUE
    )
    com = raw[1]
    # get locus length
    split.com <- strsplit(com, " ")[[1]]
    locusLength <- as.numeric(split.com[match("-r", split.com) + 2])
    pos = as.numeric(strsplit(substr(raw[6], 
                                     12, 
                                     nchar(raw[6]) - 1), 
                              split =" ")[[1]])
    pos = round(locusLength * pos)
    for (i in 1:length(pos)) {
      if (i != 1 && pos[i - 1] >= pos[i]) {
        pos[i] <- pos[i - 1] + 1
      }
    }
    if (tail(raw, 1) == "") {
      raw = head(raw, -1)
    }
    if (split.com[match("-I", split.com) + 1] == 2) {
      offset <- as.numeric(split.com[match("-I", split.com) + 3])
    } else if (is.na(split.com[match("-I", split.com)])) {
      offset <- 0
    }
    geno = Matrix(t(sapply(raw[7:(length(raw) - offset)], function(x)
      as.numeric(strsplit(x, split = "")[[1]]), USE.NAMES = F)),
      sparse = TRUE)
    ret_list = list(
      command = com,
      pos = pos,
      geno = geno
    )
    rm(raw)
    return(ret_list)
  }
}
# 
#   if( is.na(txt[1]) ){
#     print("Usage: read.ms.output(txt), or read.ms.output(file=filename)")
#     return()
#   }



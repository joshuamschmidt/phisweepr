# R functions for converting simdata objects produced by
# coalescent simulators such as ms, msms, discoal.
# #' lines are roxygen comments for documentation purposes.

#' Read in simulated data,
#'
#' @param sim.file A string giving file stroing output from a coalescent simulator in ms compatible format
#' @param software A string giving the software used. Currently only msms/ms and discoal supported.
#' @return A list containing the simulation command 
#' physical positions and genotype sparse matrix from \code{ms.file}
#' We return only the first deme, in multi deme models!
#' @examples
#' get.simulated.data(sim.file, software)
get.simulated.data <- function(sim.file = NA, software = NA, fixed_derived = FALSE) {
  if (!is.na(sim.file) && !is.na(software)) {
    if (software == "discoal") {
      returnlist <- get.discoal.output(sim.file)
    }
    if (software == "msms" || software == "ms") {
      returnlist <- get.ms.output(sim.file)
    } else {} 
  }
  # add SFS information
  one_dimensional_sfs <- get_sfs(returnlist$genotypes, fixed_derived)
  mutation_rate <- estimate_mu(one_dimensional_sfs[["derived_allele_counts"]],returnlist)
  returnlist[["derived_allele_counts"]] <- one_dimensional_sfs$derived_allele_counts
  returnlist[["one_dimensional_sfs"]] <- one_dimensional_sfs$sfs
  returnlist[["mutation_rate"]] <- mutation_rate
  return(returnlist)
}


#' Read in msms data,
#'
#' @param ms.file A string
#' @return A list containing the simulation command 
#' physical positions and genotype sparse matrix from \code{ms.file}
#' @examples
#' get.ms.output(ms.file)
get.ms.output <- function(ms.file = NA) {
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
    rho <- as.numeric(split.com[match("-r", split.com) + 1])
    locusLength <- as.integer(split.com[match("-r", split.com) + 2])
    pos = as.numeric(strsplit(substr(raw[6], 
                                     12, 
                                     nchar(raw[6]) - 1), 
                              split =" ")[[1]])
    gen_pos <- pos * rho
    pos = round(locusLength * pos)
    for (i in 1:length(pos)) {
      if (i != 1 && pos[i - 1] >= pos[i]) {
        pos[i] <- pos[i - 1] + 1
      }
    }
    
    if (tail(raw, 1) == "") {
      raw = head(raw, -1)
    }
    if (is.na(split.com[match("-I", split.com)])) {
      offset <- 0
      sample.size <- as.integer(split.com[[2]])
    } else if (!is.na(split.com[match("-I", split.com)])) {
        if (split.com[match("-I", split.com) + 1] == 2) {
          offset <- as.integer(split.com[match("-I", split.com) + 3])
          sample.size <- as.integer(split.com[match("-I", split.com) + 2])
        }
    }
    geno = Matrix::Matrix(t(sapply(raw[7:(length(raw) - offset)], function(x)
      as.numeric(strsplit(x, split = "")[[1]]), USE.NAMES = F)),
      sparse = TRUE)
    ret_list = list(
      sim_command = com,
      locusLength = locusLength,
      sample.size = sample.size,
      phys_pos = pos,
      gen_pos=gen_pos,
      genotypes = geno
    )
    rm(raw)
    return(ret_list)
  }
}

#' Read in discoal data,
#'
#' @param discoal.file A string
#' @return A list containing the simulation command 
#' physical positions and genotype sparse matrix from \code{discoal.file}
#' @examples
#' get.discoal.output(discoal.file)
get.discoal.output <- function(discoal.file = NA) {
  if (!is.na(discoal.file)) {
    raw <- scan(
      file = discoal.file,
      what = character(0),
      sep = "\n",
      blank.lines.skip = FALSE,
      quiet = TRUE
    )
    com = raw[1]
    # get locus length
    split.com <- strsplit(com, " ")[[1]]
    rho <- as.numeric(split.com[match("-r", split.com) + 1])
    locusLength <- as.integer(split.com[match("-t", split.com) -1])
    pos = as.numeric(strsplit(substr(raw[6], 
                                     12, 
                                     nchar(raw[6]) - 1), 
                              split =" ")[[1]])
    gen_pos <- pos * rho
    pos = round(locusLength * pos)
    for (i in 1:length(pos)) {
      if (i != 1 && pos[i - 1] >= pos[i]) {
        pos[i] <- pos[i - 1] + 1
      }
    }
    
    if (tail(raw, 1) == "") {
      raw = head(raw, -1)
    }
    if (split.com[match("-p", split.com) + 1] == 2) {
      offset <- as.integer(split.com[match("-p", split.com) + 3])
      sample.size <- as.integer(split.com[match("-p", split.com) + 2])
      } else if (is.na(split.com[match("-p", split.com)])) {
        offset <- 0
        sample.size <- as.integer(split.com[[2]])
    }
    geno = Matrix::Matrix(t(sapply(raw[7:(length(raw) - offset)], function(x)
      as.numeric(strsplit(x, split = "")[[1]]), USE.NAMES = F)),
      sparse = TRUE)
    ret_list = list(
      sim_command = com,
      locusLength = locusLength,
      sample.size = sample.size,
      phys_pos = pos,
      gen_pos=gen_pos,
      genotypes = geno
    )
    rm(raw)
    return(ret_list)
  }
}

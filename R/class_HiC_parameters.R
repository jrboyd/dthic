#class for storing and comparing parameter settings used in TAD analysis of HiC_matrix objects
#comparable HiC_matrixes should have identical HiC_parameters




setMethod("initialize", "HiC_parameters", function(.Object,
                                                   bin_size = 40000,
                                                   diagonal_removed = F,
                                                   canonical_chr_only = T,
                                                   depth_normalization = T,
                                                   quantile_normalization = F,
                                                   log2_over_mean_normalization = T,
                                                   n_insulation_bins = "auto",
                                                   min_insulation_distance = 0,
                                                   min_insulation_coverage = .3,
                                                   n_delta_bins = "auto"){
  .Object@bin_size = as.integer(bin_size)
  .Object@diagonal_removed = diagonal_removed
  .Object@canonical_chr_only = canonical_chr_only
  .Object@depth_normalization = depth_normalization
  .Object@quantile_normalization = quantile_normalization
  .Object@log2_over_mean_normalization = log2_over_mean_normalization
  if(n_insulation_bins == "auto"){
    #auto n_insulation_bins is the integer * bin_size >= 500kb
    n_insulation_bins = ceiling(5*10^5 / .Object@bin_size)
  }
  .Object@n_insulation_bins = as.integer(n_insulation_bins)
  .Object@min_insulation_distance = min_insulation_distance
  .Object@min_insulation_coverage = min_insulation_coverage
  if(n_delta_bins == "auto"){
    #auto n_insulation_bins is the integer * bin_size >= 100kb
    n_delta_bins = ceiling(1*10^5 / .Object@bin_size)
  }
  .Object@n_delta_bins = as.integer(n_delta_bins)
  return(.Object)
})

#' HiC_parameters
#'
#' @param bin_size
#' @param diagonal_removed
#' @param canonical_chr_only
#' @param depth_normalization
#' @param quantile_normalization
#' @param log2_over_mean_normalization
#' @param n_insulation_bins
#' @param min_insulation_distance
#' @param min_insulation_coverage
#' @param n_delta_bins
#'
#' @return
#' @export
#'
#' @examples
HiC_parameters = function(bin_size = 40000,
                          diagonal_removed = T,
                          canonical_chr_only = T,
                          depth_normalization = T,
                          quantile_normalization = F,
                          log2_over_mean_normalization = T,
                          n_insulation_bins = "auto",
                          min_insulation_distance = 0,
                          min_insulation_coverage = .3,
                          n_delta_bins = "auto"){
    new("HiC_parameters",
        bin_size = bin_size,
        diagonal_removed = diagonal_removed,
        canonical_chr_only = canonical_chr_only,
        depth_normalization = depth_normalization,
        quantile_normalization = quantile_normalization,
        log2_over_mean_normalization = log2_over_mean_normalization,
        n_insulation_bins = n_insulation_bins,
        min_insulation_distance = min_insulation_distance,
        min_insulation_coverage = min_insulation_coverage,
        n_delta_bins = n_delta_bins)
}

setMethod("show", "HiC_parameters", function(object){
  slots = names(getSlots("HiC_parameters"))
  print("HiC parameters:")
  w = max(nchar(slots))
  hidden = sapply(slots, function(s){
    print(paste(format(s, width = w), "=", slot(object, s)))
  })
  return()
})

setMethod("==", c("HiC_parameters", "HiC_parameters"), function(e1, e2){
  slots = names(getSlots("HiC_parameters"))
  pass = all(sapply(slots, function(s){
    slot(e1, s) == slot(e2, s)
  }))
  return(pass)
})

setMethod("!=", c("HiC_parameters", "HiC_parameters"), function(e1, e2){
  return(!e1 == e2)
})

# HiC_parameters() == HiC_parameters()
# HiC_parameters() != HiC_parameters()
# HiC_parameters() == HiC_parameters(canonical_chr_only = F)
# HiC_parameters() != HiC_parameters(canonical_chr_only = F)
# HiC_parameters() == HiC_parameters(bin_size = T)
# HiC_parameters() != HiC_parameters(bin_size = T)

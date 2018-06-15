

HiC_matrix_wInsulation = setClass(Class = "HiC_matrix_wInsulation", contains = "HiC_matrix")

setMethod("initialize", "HiC_matrix_wInsulation", function(.Object, base_HiC_matrix) {
  if(missing(base_HiC_matrix)){
    # warning("please supply base_HiC_matrix to generate functional HiC_matrix_wInsulation")
    return(.Object)
  }
  to_copy = names(getSlots("HiC_matrix"))
  for(s in to_copy){
    slot(.Object, s) = slot(base_HiC_matrix, s)
  }

  .Object@hic_1d = .Object@hic_1d[, .(index, seqnames, start, end)]
  .Object@hic_2d = .Object@hic_2d[i != j]
  ###score insulation across each chromosome
  full_dt = score_insulation(.Object, .Object@parameters@n_insulation_bins)

  ###conditionally apply log2 / mean normalization
  if(.Object@parameters@log2_over_mean_normalization){
    print("applying log2 mean normalization...")
    #floor value at lowest measured
    #add pseudocount

    full_dt = apply_log2mean_norm(full_dt)
    # meanv = full_dt[value > 0, mean(value, na.rm = T)]
    # pseudo = meanv / 2^10
    # full_dt[value < pseudo, value := pseudo]
    # full_dt[, value := log2(value / meanv)]
  }

  .Object@hic_1d = score_delta(hic_1d = full_dt, n_delta_bins = .Object@parameters@n_delta_bins)

  validObject(.Object)
  .Object
})

score_insulation = function(hic_mat, n_insulation_bins){
  all_chrms = unique(hic_mat@hic_1d$seqnames)
  names(all_chrms) = all_chrms
  n_bins = n_insulation_bins
  print(paste("scoring insulation for ", length(all_chrms), "chromosomes..."))
  if(exists("scores_dt"))remove(scores_dt, pos = ".GlobalEnv")
  hidden = pblapply(all_chrms, function(chr){
    s = hic_mat@hic_1d[seqnames == chr][1, start]
    e = hic_mat@hic_1d[seqnames == chr][.N, end]
    ins = insulation_of_chrRange(hic_mat, chr = chr, start = s, end = e)
    ins = ins[!is.na(names(ins))]
    dt = data.table(value = ins, index = as.integer(names(ins)))
    if(!exists("scores_dt")){
      scores_dt <<- dt
    }else{
      scores_dt <<- rbind(scores_dt, dt)
    }
  })
  setkey(scores_dt, cols = "index")

  full_dt = merge(hic_mat@hic_1d, scores_dt, by = "index", all.x = T)
  return(full_dt)
}

score_delta = function(hic_1d, n_delta_bins){
  new_1d = hic_1d
  new_1d$delta = 0
  new_1d$minmax = 0
  hidden = pbsapply(unique(new_1d$seqnames), function(chr){
    # print(chr)
    fmm_res = find_maxima_minima(vals = new_1d[seqnames == chr]$value, n_delta_bins = n_delta_bins)
    # print(paste("d-", class(fmm_res$d)))
    new_1d[seqnames == chr]$delta <<- as.numeric(fmm_res$d)
    minmax = rep(NA, length(fmm_res$vals))
    minmax[fmm_res$maxima_i] = 1
    minmax[fmm_res$minima_i] = -1
    # print(paste("mm-", class(minmax)))
    new_1d[seqnames == chr]$minmax <<- as.numeric(minmax)
  })
  return(new_1d)
}

apply_log2mean_norm = function(full_dt){
  meanv = full_dt[value > 0, mean(value, na.rm = T)]
  pseudo = meanv / 2^10
  full_dt[value < pseudo, value := pseudo]
  full_dt[, value := log2(value / meanv)]
  return(full_dt)
}

#takes an object of class HiC_matrix and returns HiC_matrix_wInsulation
#if input is already HiC_matrix_wInsulation, will overwrite
recalculate_insulation = function(hic_matrix, insulation_distance = 5*10^5, delta_distance = 1*10^5){
  hic_matrix@parameters@n_insulation_bins =
    as.integer(ceiling(insulation_distance / hic_matrix@parameters@bin_size))
  hic_matrix@parameters@n_delta_bins =
    as.integer(ceiling(delta_distance / hic_matrix@parameters@bin_size))
  HiC_matrix_wInsulation(hic_matrix)
}

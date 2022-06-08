find_maxima_minima = function(vals, n_delta_bins = 3, shift_type = c("extreme_val", "delta_zero", "none")[1]){
  #delta vector for values
  # no_delta
  d = sapply(1:length(vals),  function(i){
    if(is.na(vals[i])) return(NA)
    # if(any(i == no_delta)) return(NA)
    leftside = mean(sapply(i - 1:n_delta_bins, function(j){
      vals[j] - vals[i]
    }))
    rightside = mean(sapply(i + 1:n_delta_bins, function(j){
      vals[j] - vals[i]
    }))
    return(rightside - leftside)
  })


  #find where deltas cross zero - minima and maxima
  #return logical picking up the leftside of all such point pairs
  crosses_zero = apply(cbind(d[-1], d[-length(d)]), 1,
                       function(x){
                         if(any(is.na(x))) return(NA)
                         ifelse(x[1]*x[2] < 0, T, F)
                       })
  #indexes for left point of each pair that crosses zero
  cz = which(crosses_zero)
  #maxima go from increasing to decreasing, ie. leftside above zero
  is_maxima = d[cz] > 0
  is_minima = !is_maxima

  if(shift_type == "delta_zero"){
    closest_to_zero = apply(cbind(d[cz], d[cz + 1]), 1, function(x){
      #lowest absolute value is closest to zero, first is selected in case of tie
      (which(abs(x) == min(abs(x))) - 1)[1]
    })
    cz_closest_to_zero = cz + closest_to_zero
  }

  maxima_i = cz[is_maxima]
  minima_i = cz[is_minima]

  if(shift_type == "extreme_val"){
    # maxima_shift = apply(cbind(vals[maxima_i], vals[maxima_i + 1]), 1, function(x){
    #   (which(abs(x) == max(abs(x))) - 1)[1]
    # })
    #
    # minima_shift = apply(cbind(vals[minima_i], vals[minima_i + 1]), 1, function(x){
    #   (which(abs(x) == min(abs(x))) - 1)[1]
    # })
    if(length(maxima_i) > 0){
      maxima_shift = sapply(maxima_i, function(x){
        ir = x + c(-n_delta_bins:-1, 0, 1:n_delta_bins)
        (which(vals[ir] == max(vals[ir])))[1] - n_delta_bins - 1
      })
      maxima_i = maxima_i + maxima_shift
    }
    if(length(minima_i) > 0){
      minima_shift = sapply(minima_i, function(x){
        ir = x + c(-n_delta_bins:-1, 0, 1:n_delta_bins)
        (which(vals[ir] == min(vals[ir])))[1] - n_delta_bins - 1
      })
      minima_i = minima_i + minima_shift


    }
  }
  return(list("vals" = vals, "d" = d, "maxima_i" = maxima_i, "minima_i" = minima_i))
}

score_minima = function(fmm){

}


plot_fmm_res = function(fmm, main = "", start_f = 0, end_f = 1, maxima_col = "green", minima_col = "blue", delta_col = "red"){
  xlim = round(length(fmm$vals)*c(start_f, end_f))
  xs = 1:length(fmm$vals)
  if(!is.null(fmm$xs)){
    xlim = fmm$xs[xlim]
    xs = fmm$xs
  }


  plot(x = xs, y = fmm$vals, xlim = xlim, type = "l", ylim = c(-2,2), ylab = "value")
  title(main)
  lines(x = xs, y = fmm$d, col = delta_col)
  lines(xlim, c(0,0), col = delta_col, lty  = 2)

  points(xs[fmm$maxima_i], fmm$vals[fmm$maxima_i], col = maxima_col)
  for(i in fmm$maxima_i){
    lines(xs[rep(i, 2)], c(0, fmm$vals[i]), col = maxima_col)
  }
  points(xs[fmm$minima_i], fmm$vals[fmm$minima_i], col = minima_col)
  for(i in fmm$minima_i){
    lines(xs[rep(i, 2)], c(0, fmm$vals[i]), col = minima_col)
  }
  legend(x = "topright", legend = c('value', "delta"), col = c("black", delta_col), lty = 1)
  legend(x = "bottomright", legend = c('maxima', "minima"), col = c(maxima_col, minima_col), pch = 1)
}

ggplot_hic_delta = function(dt, chr, start, end){
  start = start + 1
  gr = GRanges(dt)
  start(gr) = start(gr) + 1
  q_gr = GRanges(chr, IRanges(start, end))
  q_indexes = dt[subjectHits(findOverlaps(query = q_gr, subject = gr, ignore.strand = TRUE))]$index
  dt = dt[q_indexes]
  dt[,xs := (start + end) / 2]
  max_k = dt$minmax == 1
  max_k[is.na(max_k)] = F
  min_k = dt$minmax == -1
  min_k[is.na(min_k)] = F

  # grps = sapply(as.character(dt$minmax), function(x)switch(x, "1" = "max", "-1" = "min", "."))
  # names(grps) = NULL
  # dt[, minmax_grp := switch(as.character(minmax), "1" = "max", "-1" = "min")]
  # dt$max_ys =
  p = ggplot(dt, aes(x = xs, y = value)) + geom_line() +
    # geom_line(mapping = aes(y = delta), col = "red") +
    annotate(geom = "point", x = dt$xs[max_k], y = dt$value[max_k], col = "red") +
    annotate(geom = "point", x = dt$xs[min_k], y = dt$value[min_k], col = "blue")
  return(p)
}

ggplot_hic_delta.df = function(dt, chr, start, end){
  start = start + 1
  gr = GRanges(dt)
  start(gr) = start(gr) + 1
  q_gr = GRanges(chr, IRanges(start, end))
  q_indexes = dt[subjectHits(findOverlaps(query = q_gr, subject = gr, ignore.strand = TRUE))]$index
  dt = dt[q_indexes]
  dt[,xs := (start + end) / 2]
  max_k = dt$minmax == 1
  max_k[is.na(max_k)] = F
  min_k = dt$minmax == -1
  min_k[is.na(min_k)] = F
  insulation_df = data.frame(x = dt$xs, y = dt$value, id = 0, type = "insulation")
  if(any(max_k)){
    insulation_df = rbind(insulation_df,
                          data.frame(x = dt$xs[max_k], y = dt$value[max_k], id = 1, type = "insulation"))
  }
  if(any(min_k)){
    insulation_df = rbind(insulation_df,
                          data.frame(x = dt$xs[min_k], y = dt$value[min_k], id = -1, type = "insulation"))
  }
  return(insulation_df)
}

# chr_vals = my_hic@hic_1d[seqnames == "chr6"]$value
# dt = my_hic@hic_1d[seqnames == "chr6"]
# xs = dt[, (end + start) / 2]
# todo = c(0,5,10)
# layout(1:length(todo))
# for(ndb in todo){
#   print(ndb)
#   fmm_res = find_maxima_minima(chr_vals, n_delta_bins = ndb)
#   fmm_res$xs = xs
#   plot_fmm_res(fmm_res, start_f = .36, end_f = .4, main = paste("n delta bins =", ndb))
# }

#calculate delta for each chromosome
#dt is data.table 1d used by HiC_matrix
#ndb is the number of delta bins, specified in HiC_parameters
calc_delta_for_hic = function(dt, ndb){
  dt$delta = 0
  dt$minmax = 0
  hidden = pbsapply(unique(dt$seqnames), function(chr){
    fmm_res = find_maxima_minima(vals = dt[seqnames == chr]$value, n_delta_bins = ndb)
    dt[seqnames == chr]$delta <<- fmm_res$d
    minmax = rep(NA, length(fmm_res$vals))
    minmax[fmm_res$maxima_i] = 1
    minmax[fmm_res$minima_i] = -1
    dt[seqnames == chr]$minmax <<- minmax
  })
  return(dt)
}
# my_hic@hic_1d = calc_delta_for_hic(my_hic@hic_1d, my_hic@parameters@n_delta_bins)


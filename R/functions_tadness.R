# library(data.table) #super fast data handling
# library(GenomicRanges) #genomic coord intersection
# library(pbapply) #progress bars
# library(preprocessCore) #quantile normalization
# 
# # #calculates the average value of square sub-matrix extended from diagonal at pos_index position by n_bins in both directions
# # insulation_of_index = function(hic_matrix, pos_index, n_bins){
# #   rng = insulation_range(pos_index, n_bins)
# #   matrix_subset = hic_matrix[rng[[1]], rng[[2]]]
# #   ins = sum(matrix_subset$val, na.rm = T) / n_bins^2
# #   return(ins)
# # }
# 
# # #runs insulation_of_index on all indexes in pos_indexes
# # insulation_of_indexes = function(hic_matrix, pos_indexes, n_bins){
# #   ins = sapply(pos_indexes, function(pos_index){
# #     insulation_of_index(hic_matrix, pos_index, n_bins)
# #   })
# #   names(ins) = pos_indexes
# #   return(ins)
# # }
# 
# chr_list2gr = function(chr_list, bin_size){
#   i = 1
#   x = chr_list[[i]]
#   chr = names(chr_list)[i]
#   idx = 1:length(x)
#   GRanges(seqnames = chr, ranges = IRanges(idx * bin_size - bin_size + 1, idx * bin_size))
# }
# 
# #returns centers of bins that intersect the chromosom range
# #most useful for getting genomic coordinates for plotting
# bin_centers_of_chrRange = function(hic_matrix, chr, start, end){
#   m_gr = GRanges(hic_matrix@regions)
#   start(m_gr) = start(m_gr) + 1
#   q_gr = GRanges(seqnames = chr, IRanges(start + 1, end))
#   pos_indexes = subjectHits(findOverlaps(query = q_gr, subject = m_gr))
#   start(m_gr[pos_indexes]) + width(m_gr[pos_indexes])/2
# }
# 
# plot_insulation_score_by_set = function(HiC_set, chr, s, e, score_colors = NA){
#   if(is.na(score_colors)){
#     score_colors = RColorBrewer::brewer.pal(length(HiC_matrix_list), "Paired")
#     names(score_colors) = names(HiC_matrix_list)
#   }
#   inde
#   scores_list = pblapply(HiC_matrix_list, function(x)insulation_of_chrRange(x, chr = chr, start = s, end = e, n_bins = n_bins))
#   ylab = "insulation score"
#   
#   if(apply_normalize.log2_ratio){
#     scores_list = lapply(scores_list, function(x){
#       x = ifelse(x <= 0, min(x[x > 0]), x)
#       x = x / mean(x)
#       log2(x)
#     })
#   }
#   
#   if(apply_normalize.quantiles){
#     #apply quantile normalization
#     all_mat = matrix(unlist(scores_list), ncol = length(scores_list))
#     all_mat_normq = normalize.quantiles(all_mat)
#     normq_list = lapply(1:ncol(all_mat_normq), function(i)all_mat_normq[,i])
#     names(normq_list) = names(scores_list)
#     scores_list = normq_list
#     ylab = "quantile normalized insulation score"
#   }
#   
#   delta_list = lapply(scores_list, function(x){
#     no_delta = c(1:n_delta_bins, (length(x)-n_delta_bins+1):length(x))
#     sapply(1:length(x), function(i){
#       if(any(i == no_delta)){
#         return(NA)
#       }
#       x[i + n_delta_bins] - x[i - n_delta_bins]
#     })
#   })
#   
#   p_max = max(sapply(scores_list, max))
#   p_min = min(sapply(scores_list, min), 0)
#   
#   plot(scores_list[[1]], type = "n", xlim = c(s, e), ylim = c(p_min, p_max), ylab = ylab, xlab = paste0(chr, ": position (bp)"))
#   xs = bin_centers_of_chrRange(all_HiC_matrixes[[1]], chr, s, e)
#   for(i in 1:length(scores_list)){
#     lines(xs, scores_list[[i]], col = score_colors[names(scores_list)[i]])
#   }
#   legend("top", legend = names(scores_list), fill = score_colors[names(scores_list)], ncol = 3, xpd = NA, cex = .7)
#   
#   # for(i in 1:length(scores_list)){
#   #   lines(xs, delta_list[[i]], col = "red")
#   # }
# }
# 
# 
# plot_insulation_score = function(HiC_set, chr, s, e, n_bins = 15, n_delta_bins = 2, score_colors = NA, show_delta = T, apply_normalize.quantiles = F, apply_normalize.log2_ratio = T){
#   if(is.na(score_colors)){
#     score_colors = RColorBrewer::brewer.pal(length(HiC_matrix_list), "Paired")
#     names(score_colors) = names(HiC_matrix_list)
#   }
#   scores_list = pblapply(HiC_matrix_list, function(x)insulation_of_chrRange(x, chr = chr, start = s, end = e, n_bins = n_bins))
#   ylab = "insulation score"
#   
#   if(apply_normalize.log2_ratio){
#     scores_list = lapply(scores_list, function(x){
#       x = ifelse(x <= 0, min(x[x > 0]), x)
#       x = x / mean(x)
#       log2(x)
#     })
#   }
#   
#   if(apply_normalize.quantiles){
#     #apply quantile normalization
#     all_mat = matrix(unlist(scores_list), ncol = length(scores_list))
#     all_mat_normq = normalize.quantiles(all_mat)
#     normq_list = lapply(1:ncol(all_mat_normq), function(i)all_mat_normq[,i])
#     names(normq_list) = names(scores_list)
#     scores_list = normq_list
#     ylab = "quantile normalized insulation score"
#   }
#   
#   delta_list = lapply(scores_list, function(x){
#     no_delta = c(1:n_delta_bins, (length(x)-n_delta_bins+1):length(x))
#     sapply(1:length(x), function(i){
#       if(any(i == no_delta)){
#         return(NA)
#       }
#       x[i + n_delta_bins] - x[i - n_delta_bins]
#     })
#   })
#   
#   p_max = max(sapply(scores_list, max))
#   p_min = min(sapply(scores_list, min), 0)
#   
#   plot(scores_list[[1]], type = "n", xlim = c(s, e), ylim = c(p_min, p_max), ylab = ylab, xlab = paste0(chr, ": position (bp)"))
#   xs = bin_centers_of_chrRange(all_HiC_matrixes[[1]], chr, s, e)
#   for(i in 1:length(scores_list)){
#     lines(xs, scores_list[[i]], col = score_colors[names(scores_list)[i]])
#   }
#   legend("top", legend = names(scores_list), fill = score_colors[names(scores_list)], ncol = 3, xpd = NA, cex = .7)
#   
#   # for(i in 1:length(scores_list)){
#   #   lines(xs, delta_list[[i]], col = "red")
#   # }
# }
# 
# plot_delta_score = function(HiC_matrix_list, chr, s, e, n_bins = 15, n_delta_bins = 2, score_colors = NA, show_delta = T, apply_normalize.quantiles = F, apply_normalize.log2_ratio = T){
#   if(is.na(score_colors)){
#     score_colors = RColorBrewer::brewer.pal(length(HiC_matrix_list), "Paired")
#     names(score_colors) = names(HiC_matrix_list)
#   }
#   scores_list = pblapply(HiC_matrix_list, function(x)insulation_of_chrRange(x, chr = chr, start = s, end = e, n_bins = n_bins))
#   ylab = "insulation score"
#   
#   if(apply_normalize.log2_ratio){
#     scores_list = lapply(scores_list, function(x){
#       x = ifelse(x <= 0, min(x[x > 0]), x)
#       x = x / mean(x)
#       log2(x)
#     })
#   }
#   
#   if(apply_normalize.quantiles){
#     #apply quantile normalization
#     all_mat = matrix(unlist(scores_list), ncol = length(scores_list))
#     all_mat_normq = normalize.quantiles(all_mat)
#     normq_list = lapply(1:ncol(all_mat_normq), function(i)all_mat_normq[,i])
#     names(normq_list) = names(scores_list)
#     scores_list = normq_list
#     ylab = "quantile normalized insulation score"
#   }
#   
#   delta_list = lapply(scores_list, function(x){
#     no_delta = c(1:n_delta_bins, (length(x)-n_delta_bins+1):length(x))
#     sapply(1:length(x), function(i){
#       if(any(i == no_delta)){
#         return(NA)
#       }
#       x[i + n_delta_bins] - x[i - n_delta_bins]
#     })
#   })
#   
#   p_max = max(sapply(scores_list, max))
#   p_min = min(sapply(scores_list, min), 0)
#   
#   plot(scores_list[[1]], type = "n", xlim = c(s, e), ylim = c(p_min, p_max), ylab = ylab, xlab = paste0(chr, ": position (bp)"))
#   xs = bin_centers_of_chrRange(all_HiC_matrixes[[1]], chr, s, e)
#   for(i in 1:length(scores_list)){
#     lines(xs, scores_list[[i]], col = score_colors[names(scores_list)[i]])
#   }
#   legend("top", legend = names(scores_list), fill = score_colors[names(scores_list)], ncol = 3, xpd = NA, cex = .7)
#   
#   for(i in 1:length(scores_list)){
#     lines(xs, delta_list[[i]], col = "red")
#   }
# }

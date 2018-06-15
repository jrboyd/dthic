library(dthic)
hparm = HiC_parameters(bin_size = 150000, depth_normalization = FALSE)
hic_mats = list(
    MCF10AT1 = HiC_matrix(matrix_file = "~/HiC-Pro/outputs/MCF10AT1_pooled/hic_results/matrix/pooled/iced/150000/pooled_150000_iced.matrix",
                          regions_file = "~/HiC-Pro/outputs/MCF10AT1_pooled/hic_results/matrix/pooled/raw/150000/pooled_150000_abs.bed",
                          hic_parameters = hparm),
    MCF10CA1a = HiC_matrix(matrix_file = "~/HiC-Pro/outputs/MCF10CA1a_pooled/hic_results/matrix/pooled/iced/150000/pooled_150000_iced.matrix",
                           regions_file = "~/HiC-Pro/outputs/MCF10CA1a_pooled/hic_results/matrix/pooled/raw/150000/pooled_150000_abs.bed",
                           hic_parameters = hparm)
)

chr = chr
s = 18*10^6
e = 35*10^6

plot_upperMatrix(hic_mats[[1]], chr, s, e)

hic_ins = lapply(hic_mats, HiC_matrix_wInsulation)

plot_upperMatrix_with_insulation(hic_ins[[1]], chr, s, e)

z_mats = lapply(hic_mats, apply_diagonal_zscore)

zcol = c("darkblue", "blue", "gray", "red", "darkred")
plot_upperMatrix(z_mats[[1]], chr, s, e, hmap_colors = zcol)

z_ins = lapply(z_mats, HiC_matrix_wInsulation)

plot_upperMatrix_with_insulation(z_ins[[1]], chr, s, e, hmap_colors = zcol, max_dist = e - s)

chr = "chr2"
s = 20*10^6
e = 35*10^6
plot_upperMatrix_with_insulation(hic_ins[[1]], chr, s, e, max_dist = e - s)
plot_upperMatrix_with_insulation(z_ins[[1]], chr, s, e, hmap_colors = zcol, max_dist = e - s)


hic_di = lapply(hic_mats, HiC_matrix_wDirectionality)

hcol = c("white", "slategray",
         rep("blue", 3), "yellow", rep("orange", 3),
         "red", rep("magenta", 3))

chr = "chr6"
s = 25*10^6
e = 51*10^6
plot_upperMatrix_with_insulation(
    hic_ins[[1]], chr, s, e,
    max_dist = e - s, max_fill = 4,
    hmap_colors = hcol)
plot_upperMatrix_with_insulation(z_ins[[1]], chr, s, e, hmap_colors = zcol, max_dist = e - s)

plot_upperMatrix_with_insulation(
    hic_di[[1]], chr, s, e,
    max_fill = 500, show_minmax = FALSE,
    hmap_colors = hcol,
    max_dist = e - s)

plot(directionality_of_chrRange(hic_mats[[1]], chr, s, e), type = "l")
ins = insulation_of_chrRange(hic_mats[[1]], chr, s, e)
ins = ins - min(ins)
ins = ins / max(ins)
ins = ins - mean(ins)
ins = ins * 3000
lines(ins, type = "l", col = "blue")
lines(c(0, 180), y = c(0,0), lty = 2)

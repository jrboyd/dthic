library(dthic)
hparm = HiC_parameters(bin_size = 40000, depth_normalization = FALSE)
hic_mats = list(
    MCF10AT1 = HiC_matrix(matrix_file = "~/HiC-Pro/outputs/MCF10AT1_pooled/hic_results/matrix/pooled/iced/40000/pooled_40000_iced.matrix",
                          regions_file = "~/HiC-Pro/outputs/MCF10AT1_pooled/hic_results/matrix/pooled/raw/40000/pooled_40000_abs.bed",
                          parameters = hparm)#,
    # MCF10CA1a = HiC_matrix(matrix_file = "~/HiC-Pro/outputs/MCF10CA1a_pooled/hic_results/matrix/pooled/iced/40000/pooled_40000_iced.matrix",
    #                        regions_file = "~/HiC-Pro/outputs/MCF10CA1a_pooled/hic_results/matrix/pooled/raw/40000/pooled_40000_abs.bed",
    #                        hic_parameters = hparm),
    # HepG2 =
    #     HiC_matrix(matrix_file = "~/HiC-Pro/outputs/HepG2/hic_results/matrix/rep1/raw/40000/rep1_40000.matrix",
    #                regions_file = "~/HiC-Pro/outputs/HepG2/hic_results/matrix/rep1/raw/40000/rep1_40000_abs.bed",
    #                hic_parameters = hparm) +
    #     HiC_matrix(matrix_file = "~/HiC-Pro/outputs/HepG2/hic_results/matrix/rep2/raw/40000/rep2_40000.matrix",
    #                regions_file = "~/HiC-Pro/outputs/HepG2/hic_results/matrix/rep2/raw/40000/rep2_40000_abs.bed",
    #                hic_parameters = hparm)
    # HepG2 =
    #     HiC_matrix(matrix_file = "~/HiC-Pro/outputs/HepG2/hic_results/matrix/rep1/iced/40000/rep1_40000_iced.matrix",
    #                regions_file = "~/HiC-Pro/outputs/HepG2/hic_results/matrix/rep1/raw/40000/rep1_40000_abs.bed",
    #                hic_parameters = hparm) +
    #     HiC_matrix(matrix_file = "~/HiC-Pro/outputs/HepG2/hic_results/matrix/rep2/iced/40000/rep2_40000_iced.matrix",
    #                regions_file = "~/HiC-Pro/outputs/HepG2/hic_results/matrix/rep2/raw/40000/rep2_40000_abs.bed",
    #                hic_parameters = hparm)
)

hic_ins = lapply(hic_mats, HiC_matrix_wInsulation)
hic_di = lapply(hic_mats, HiC_matrix_wDirectionality)

hcol = c("white", "slategray",
         rep("blue", 3), rep("orange", 3),
         rep("red", 3), "magenta")

chr = "chr6"
s = 15*10^6
e = 51*10^6

s = 18*10^6
e = 32*10^6




pdf("chr6_hic.pdf")
plot_upperMatrix_with_insulation(
    hic_ins[[3]], chr, s, e,
    max_fill = 500, show_minmax = FALSE,
    hmap_colors = hcol,
    max_dist = e - s)

plot_upperMatrix_with_insulation(
    hic_di[[3]], chr, s, e,
    max_fill = 500, show_minmax = FALSE,
    hmap_colors = hcol,
    max_dist = e - s)

plot_upperMatrix_with_insulation(
    hic_ins[[1]], chr, s, e,
    max_fill = 500, show_minmax = FALSE,
    hmap_colors = hcol,
    max_dist = e - s)

plot_upperMatrix_with_insulation(
    hic_di[[1]], chr, s, e,
    max_fill = 500, show_minmax = FALSE,
    hmap_colors = hcol,
    max_dist = e - s)

plot_upperMatrix_with_insulation(
    hic_ins[[2]], chr, s, e,
    max_fill = 500, show_minmax = FALSE,
    hmap_colors = hcol,
    max_dist = e - s)

plot_upperMatrix_with_insulation(
    hic_di[[2]], chr, s, e,
    max_fill = 500, show_minmax = FALSE,
    hmap_colors = hcol,
    max_dist = e - s)
dev.off()

#DNAJB1 chr19:14,514,770-14,518,420
#PRKACA chr19:14,091,688-14,117,744
chr = "chr19"
s = 14091688
e = 14518420

s = s - 10*10^6
e = e + 10*10^6

pdf("FLC_fusion_hic.pdf")
plot_upperMatrix_with_insulation(
    hic_ins[[3]], chr, s, e,
    max_fill = 500, show_minmax = FALSE,
    hmap_colors = hcol,
    max_dist = e - s)

plot_upperMatrix_with_insulation(
    hic_di[[3]], chr, s, e,
    max_fill = 500, show_minmax = FALSE,
    hmap_colors = hcol,
    max_dist = e - s)

plot_upperMatrix_with_insulation(
    hic_ins[[1]], chr, s, e,
    max_fill = 500, show_minmax = FALSE,
    hmap_colors = hcol,
    max_dist = e - s)

plot_upperMatrix_with_insulation(
    hic_di[[1]], chr, s, e,
    max_fill = 500, show_minmax = FALSE,
    hmap_colors = hcol,
    max_dist = e - s)

plot_upperMatrix_with_insulation(
    hic_ins[[2]], chr, s, e,
    max_fill = 500, show_minmax = FALSE,
    hmap_colors = hcol,
    max_dist = e - s)

plot_upperMatrix_with_insulation(
    hic_di[[2]], chr, s, e,
    max_fill = 500, show_minmax = FALSE,
    hmap_colors = hcol,
    max_dist = e - s)
dev.off()

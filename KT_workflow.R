library(dthic)
library(pbapply)
library(data.table)
library(gridExtra)
# source("R/class_HiC_matrix_helpers.R")
# source("R/class_HiC_matrix_wInsulation.R")
source("~/R/jrb_R_scripts/parse_gtf.dt.R")

ref_dt = parse_gtf.dt("~/gencode.v27.annotation.gtf", additional_attribs = "gene_type")
ref_gr = GRanges(ref_dt)
# hic_mat = HiC_matrix(matrix_file = "~/HiC-Pro/outputs/MCF10A_pooled/hic_results/matrix/MCF10A_20split/iced/40000/MCF10A_20split_40000_iced.matrix",
#                      regions_file = "~/HiC-Pro/outputs/MCF10A_pooled/hic_results/matrix/MCF10A_20split/raw/40000/MCF10A_20split_40000_abs.bed",
#                      hic_parameters = HiC_parameters(bin_size = 40000))
hic_mats = list(
    MCF10A = HiC_matrix(matrix_file = "~/HiC-Pro/outputs/MCF10A_pooled/hic_results/matrix/MCF10A_20split/iced/150000/MCF10A_20split_150000_iced.matrix",
                        regions_file = "~/HiC-Pro/outputs/MCF10A_pooled/hic_results/matrix/MCF10A_20split/raw/150000/MCF10A_20split_150000_abs.bed",
                        hic_parameters = HiC_parameters(bin_size = 150000)),
    MCF10AT1 = HiC_matrix(matrix_file = "~/HiC-Pro/outputs/MCF10AT1_pooled/hic_results/matrix/pooled/iced/150000/pooled_150000_iced.matrix",
                          regions_file = "~/HiC-Pro/outputs/MCF10AT1_pooled/hic_results/matrix/pooled/raw/150000/pooled_150000_abs.bed",
                          hic_parameters = HiC_parameters(bin_size = 150000)),
    MCF10CA1a = HiC_matrix(matrix_file = "~/HiC-Pro/outputs/MCF10CA1a_pooled/hic_results/matrix/pooled/iced/150000/pooled_150000_iced.matrix",
                           regions_file = "~/HiC-Pro/outputs/MCF10CA1a_pooled/hic_results/matrix/pooled/raw/150000/pooled_150000_abs.bed",
                           hic_parameters = HiC_parameters(bin_size = 150000))
)

# hic_mats = list(
#     MCF10A = HiC_matrix(matrix_file = "~/HiC-Pro/outputs/MCF10A_pooled/hic_results/matrix/MCF10A_20split/iced/40000/MCF10A_20split_40000_iced.matrix",
#                         regions_file = "~/HiC-Pro/outputs/MCF10A_pooled/hic_results/matrix/MCF10A_20split/raw/40000/MCF10A_20split_40000_abs.bed",
#                         hic_parameters = HiC_parameters(bin_size = 40000)),
#     MCF10AT1 = HiC_matrix(matrix_file = "~/HiC-Pro/outputs/MCF10AT1_pooled/hic_results/matrix/pooled/iced/40000/pooled_40000_iced.matrix",
#                           regions_file = "~/HiC-Pro/outputs/MCF10AT1_pooled/hic_results/matrix/pooled/raw/40000/pooled_40000_abs.bed",
#                           hic_parameters = HiC_parameters(bin_size = 40000)),
#     MCF10CA1a = HiC_matrix(matrix_file = "~/HiC-Pro/outputs/MCF10CA1a_pooled/hic_results/matrix/pooled/iced/40000/pooled_40000_iced.matrix",
#                            regions_file = "~/HiC-Pro/outputs/MCF10CA1a_pooled/hic_results/matrix/pooled/raw/40000/pooled_40000_abs.bed",
#                            hic_parameters = HiC_parameters(bin_size = 40000))
# )

hic_matIns = lapply(hic_mats, HiC_matrix_wInsulation)
hic_matIns_z = lapply(hic_mats, apply_diagonal_zscore)

# qgr = GRanges("chr10:4650185-5458463")
qgr = GRanges("chr10:1-12458463")

pdf(paste0("KT_", "zscore_zoomout.pdf"), width = 18, height = 12)
for(i in 1:length(hic_mats)){
    hic_mat = hic_matIns_z[[i]]
    plots_mat = plot_upperMatrix_with_insulation(hic_mat = hic_mat, hmap_colors = c("blue", "lightgray", "red"), max_dist = 4*10^6,
                                                 chr = as.character(seqnames(qgr)),
                                                 start = start(qgr), end = end(qgr), show_plot = F)

    high_colors = c("tomato", "darkgreen")
    high_colors = rgb(t(col2rgb(high_colors))/255)
    highlight_gene_names = c("LINC00704", "NET1")
    refp = ggplot_ref(ref_gr, qgr = qgr, gene_types = c("protein_coding", "lincRNA"),
                      arrow_override = NULL, top_spacer = 1, text_size = 2,
                      text_angle = 43, text_y_relative = 0,
                      highlight_gene_names = highlight_gene_names, highlight_gene_color = "red")
    plots = append(plots_mat$ggplots, list(ref = refp))

    plots_rects = lapply(plots, function(x){
        gp = ggplot_build(x)
        ylim = gp$layout$panel_ranges[[1]]$y.range
        for(gn in highlight_gene_names){
            x = x + annotate("rect",
                             xmin = ref_dt[gene_name == gn]$start,
                             xmax = ref_dt[gene_name == gn]$end,
                             ymin = min(ylim), ymax = max(ylim),
                             fill = paste0(rgb(t(col2rgb("red")/255)), "44"))
        }
        x
    })

    grobs = ggplotList2grobList(plots_rects)

    grid.arrange(grobs = grobs, heights = c(4, 2, 3), top = names(hic_mats)[i])

}
dev.off()

pdf(paste0("KT_", "normal_zoomout.pdf"), width = 18, height = 12)
for(i in 1:length(hic_mats)){
    hic_mat = hic_matIns[[i]]
    plots_mat = plot_upperMatrix_with_insulation(hic_mat = hic_mat, max_dist = 4*10^6,
                                                 chr = as.character(seqnames(qgr)),
                                                 start = start(qgr), end = end(qgr), show_plot = F)

    high_colors = c("tomato", "darkgreen")
    high_colors = rgb(t(col2rgb(high_colors))/255)
    highlight_gene_names = c("LINC00704", "NET1")
    refp = ggplot_ref(ref_gr, qgr = qgr, gene_types = c("protein_coding", "lincRNA"),
                      arrow_override = NULL, top_spacer = 1, text_size = 2,
                      text_angle = 43, text_y_relative = 0,
                      highlight_gene_names = highlight_gene_names, highlight_gene_color = "red")
    plots = append(plots_mat$ggplots, list(ref = refp))

    plots_rects = lapply(plots, function(x){
        gp = ggplot_build(x)
        ylim = gp$layout$panel_ranges[[1]]$y.range
        for(gn in highlight_gene_names){
            x = x + annotate("rect",
                             xmin = ref_dt[gene_name == gn]$start,
                             xmax = ref_dt[gene_name == gn]$end,
                             ymin = min(ylim), ymax = max(ylim),
                             fill = paste0(rgb(t(col2rgb("red")/255)), "44"))
        }
        x
    })

    grobs = ggplotList2grobList(plots_rects)

    grid.arrange(grobs = grobs, heights = c(4, 2, 3), top = names(hic_mats)[i])

}
dev.off()

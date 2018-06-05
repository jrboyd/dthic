data_source = setClass(Class = "data_source")
setGeneric("plot_from_data_source", def = function(object, gr, ...){stop("please use a class that contains data_source")})

###bigwig data and plot
data_source.bigwig =
  setClass(Class = "data_source.bigwig", contains = "data_source",
           slots = c("bigwig_path" = "character"))
setMethod("initialize", "data_source.bigwig", function(.Object, bigwig_path) {
  if(missing(bigwig_path)){
    return(.Object)
  }
  .Object@bigwig_path = bigwig_path
  validObject(.Object)
  .Object

})
setMethod("plot_from_data_source",
          signature = c("data_source.bigwig", "GenomicRanges"),
          definition = function(object, gr,
                                n_bins = 200,
                                bin_method = c("mean", "max")[2],
                                bigwig_colors = NULL,
                                bigwig_title = "FE",
                                p = NULL,
                                fe_max = NULL,
                                fe_min = NULL,
                                alpha = NULL, ...){
            add_bigwig_plot.tiles(bigwig_file = object@bigwig_path, chr = as.character(seqnames(gr)), start = start(gr), end = end(gr),
                                  n_bins = n_bins,
                                  bin_method = bin_method,
                                  bigwig_title = bigwig_title,
                                  bigwig_colors = bigwig_colors,
                                  p = p,
                                  fe_min = fe_min,
                                  fe_max = fe_max,
                                  alpha = alpha, ...)
          }
)

###insulation data and plot
data_source.hic_insulation =
  setClass(Class = "data_source.hic_insulation", contains = "data_source",
           slots = c("hic_mat_wIns" = "HiC_matrix_wInsulation"))
setMethod("initialize", "data_source.hic_insulation", function(.Object, hic_mat_wIns) {
  if(missing(hic_mat_wIns)){
    return(.Object)
  }
  .Object@hic_mat_wIns = hic_mat_wIns
  validObject(.Object)
  .Object

})
setMethod("plot_from_data_source",
          signature = c("data_source.hic_insulation", "GenomicRanges"),
          definition = function(object, gr){
            df = object@hic_mat_wIns@hic_1d[queryHits(findOverlaps(GRanges(object@hic_mat_wIns@hic_1d), gr))]
            df$value[df$value == -10] = NA
            df[, x := (start + end) / 2]
            df[, y := value]
            df[is.na(minmax), minmax := 0]
            df$fill = sapply(as.character(df$minmax), function(x)switch(x, "-1" = "min", "0" = "-", "1" = "max", "NA"))
            # df[, fill := switch(as.character(minmax), "-1" = "min", "0" = "-", "1" = "max")]
            p = ggplot() + geom_line(data = df, mapping = aes(x = x, y = y)) + scale_y_continuous(breaks = hic_equal_breaks(bin_size = bin_size, max_dist = max_dist)) +
              labs(x = paste(as.character(seqnames(gr)), "position"),
                   y = "log2(insulation / mean)", fill = "") +
              coord_cartesian(xlim = c(start(gr), end(gr))) +
              scale_size_identity() +
              scale_color_identity() +
              scale_shape_identity()
            if(nrow(df[fill != "-"])){
              p = p + geom_point(data = df[fill != "-"],
                                 mapping = aes(shape = 21, x = x, y = y, fill = fill, color = "black", size = 3, stroke = 1.5)) +
                scale_fill_manual(values = c("min" = "orange", "max" = "green"),
                                  guide = guide_legend(override.aes = list(shape = 21, size = 3, stroke = 1.5)))
            }
            return(p)
          })


###interaction arc plot
data_source.hic_arc =
  setClass(Class = "data_source.hic_arc", contains = "data_source",
           slots = c(hic_mat = "HiC_matrix"))
setMethod("initialize", "data_source.hic_arc", function(.Object, hic_mat) {
  if(missing(hic_mat)){
    return(.Object)
  }
  .Object@hic_mat = hic_mat
  validObject(.Object)
  .Object
})
setMethod("plot_from_data_source",
          signature = c("data_source.hic_arc", "GenomicRanges"),
          definition = function(object, gr){
            add_arch_plot(hic_mat = object@hic_mat, qgr = gr)
          })


###matrix data and plot
data_source.hic_matrix =
  setClass(Class = "data_source.hic_matrix", contains = "data_source",
           slots = c(hic_mat = "HiC_matrix"))
setMethod("initialize", "data_source.hic_matrix", function(.Object, hic_mat) {
  if(missing(hic_mat)){
    return(.Object)
  }
  .Object@hic_mat = hic_mat
  validObject(.Object)
  .Object
})
setMethod("plot_from_data_source",
          signature = c("data_source.hic_matrix", "GenomicRanges"),
          definition = function(object, gr,
                                tile_type = c("diamond", "hex")[1],
                                max_dist = 10*10^6,
                                point_size = 3,
                                text_size = 4,
                                max_fill = NULL,
                                hmap_colors = c("lightgray", "steelblue", 'darkblue', "red", "red"),
                                show_insulation_range = T,
                                show_min = T,
                                show_max = F){
            add_matrix_plot(object@hic_mat, qgr = gr,
                            tile_type = tile_type,
                            max_dist = max_dist,
                            point_size = point_size,
                            text_size = text_size,
                            max_fill = max_fill,
                            hmap_colors = hmap_colors,
                            show_insulation_range = show_insulation_range,
                            show_min = show_min,
                            show_max = show_max)
            # # print("plot from data_source.hic_matrix")
            # hic_mat = object@hic_mat
            # bin_size = hic_mat@parameters@bin_size
            # # hicrng = get_chrRange_Matrix(hic_mat@hic_2d, hic_mat@hic_1d, chr, start, end)
            # chr = as.character(seqnames(gr))
            # start = start(gr)
            # end = end(gr)
            # hicrng = hic_mat[chr, c(start,end)]
            # MIN_I = min(c(hicrng$i, hicrng$j))
            # XMIN = subset(hic_mat@hic_1d, index == MIN_I)$start
            # MAX_I = max(c(hicrng$i, hicrng$j))
            # XMAX = subset(hic_mat@hic_1d, index == MAX_I)$end
            # hicrng[, x := (i + j) /2 ]
            # hicrng[, y := abs(i - j)]
            #
            # hicrng = hicrng[j > i] #limit to upper triangle
            # hicrng = na.omit(hicrng) #leave out zero signal regions
            #
            # index2center = function(index){
            #   (index - MIN_I)/(MAX_I - MIN_I) * (XMAX - XMIN - bin_size) + XMIN + bin_size / 2
            # }
            #
            # hicrng[, x_chr := index2center(x)]
            # hicrng[, y_chr := bin_size * y]
            #
            # if(tile_type == "hex"){
            #   df = hex_coord_df(hicrng$x_chr, hicrng$y_chr, width = bin_size, height = bin_size, size = 1)
            #   df = cbind(df, fill=rep(hicrng$val, each=6))
            # }else if(tile_type == "diamond"){
            #   df = diamond_coord_df(hicrng$x_chr, hicrng$y_chr, width = bin_size, height = bin_size, size = 1)
            #   df = cbind(df, fill=rep(hicrng$val, each=4))
            # }else{
            #   stop("tile_type must match hex or diamond")
            # }
            # #pretty up break positions to be multiple of bin_size
            # ylab = sort(unique(df$y))
            # nlab = 6
            # yfac = min(ceiling(max_dist / nlab / bin_size) * bin_size, #when max_dist is less than start to end
            #            ceiling(length(ylab) / nlab) * bin_size) #when start and end are less than max_dist
            # breaks = 0:nlab * yfac
            #
            # fill_lab = "interaction"
            # if(!is.null(max_fill)){
            #   df$fill[df$fill > max_fill] = max_fill
            #   fill_lab = paste(fill_lab, "\ncapped at", max_fill)
            # }
            # df = subset(df, y <= max_dist)
            #
            # mmdf = make_hic_minmax_df(object@hic_mat@hic_1d, gr)
            # mmdf$minmax = sapply(as.character(mmdf$id), function(x)switch(x, "-1" = "min", "0" = "-", "1" = "max"))
            # ann_df = subset(mmdf, id != 0)
            #
            # minmax2col = c("min" = "darkgreen", "max" = "orange")
            #
            # p = ggplot() + geom_polygon(data = df, aes(x=x, y=y, fill = fill, group = id)) +
            #   scale_fill_gradientn(colours = hmap_colors) +
            #   labs(x = "", y = "distance", fill = fill_lab) +
            #   scale_y_continuous(breaks = hic_equal_breaks(bin_size = bin_size, max_dist = max_dist)) +
            #   coord_cartesian(xlim = c(start, end))
            # if(nrow(ann_df) > 0){
            #   if(show_max){
            #     max_df = subset(ann_df, minmax == "max")
            #     p = p + annotate("point",
            #                      x = max_df$x,
            #                      y = 0 - max(df$y)*.05,
            #                      fill = minmax2col[max_df$minmax],
            #                      color = "#00000000",
            #                      size = point_size,
            #                      stroke = 1.5,
            #                      shape = 24)
            #     p = p + annotate("point",
            #                      x = start + (end - start)*.9,
            #                      y = max(df$y)*.88,
            #                      fill = minmax2col["max"],
            #                      color = "#00000000",
            #                      size = point_size,
            #                      stroke = 1.5,
            #                      shape = 24)
            #     p = p + annotate("text",
            #                      label = "maxima",
            #                      x = start + (end - start)*.92,
            #                      y = max(df$y)*.89,
            #                      color = minmax2col["max"],
            #                      size = text_size,
            #                      hjust = 0,
            #                      vjust = .5)
            #   }
            #   if(show_min){
            #     min_df = subset(ann_df, minmax == "min")
            #     p = p + annotate("point",
            #                      x = min_df$x,
            #                      y = 0 - max(df$y)*.03,
            #                      fill = minmax2col[min_df$minmax],
            #                      color = "#00000000",
            #                      size = point_size,
            #                      stroke = 1.5,
            #                      shape = 25)
            #     p = p + annotate("point",
            #                      x = start + (end - start)*.9,
            #                      y = max(df$y)*.95,
            #                      fill = minmax2col["min"],
            #                      color = "#00000000",
            #                      size = point_size,
            #                      stroke = 1.5,
            #                      shape = 25)
            #     p = p + annotate("text",
            #                      label = "minima",
            #                      x = start + (end - start)*.92,
            #                      y = max(df$y)*.95,
            #                      color = minmax2col["min"],
            #                      size = text_size,
            #                      hjust = 0,
            #                      vjust = .5)
            #   }
            # }
            #
            # if(show_insulation_range){
            #   p = p + annotate("line",
            #                    x = c(start, end),
            #                    y = hic_mat@parameters@n_insulation_bins * hic_mat@parameters@bin_size)
            #   p = p + annotate("text",
            #                    label = "insulation\nbin range",
            #                    x = c(end),
            #                    y = hic_mat@parameters@n_insulation_bins * hic_mat@parameters@bin_size,
            #                    size = text_size,
            #                    hjust = 1,
            #                    vjust = -.2)
            # }
            # return(p)
          })

###annotation reference data and plot
data_source.ref =
  setClass(Class = "data_source.ref", contains = "data_source",
           slots = c("ref" = "GRanges"))

setMethod("initialize", "data_source.ref", function(.Object, ref_gr) {
  if(missing(ref_gr)){
    return(.Object)
  }
  .Object@ref = ref_gr
  validObject(.Object)
  .Object
})

setMethod("plot_from_data_source",
          signature = c("data_source.ref", "GRanges"),
          definition = function(object, gr, ...){
            ggplot_ref(object@ref, gr, ...)
          })
#
#
# #class to contain data for repeated plot generation
# #can be queried for genomic regions to generate an integrated_plot
# integrated_data_source = setClass(Class = "integrated_data_source",
#
#                                   slots = c(
#                                     grobs_list = "list",
#                                     ref_dt = "data.table",
#                                     parameters = "HiC_parameters",
#                                     hic_2d = "data.table",
#                                     hic_1d = "data.table"
#
#                                   )
#
#                                   # validity = function(object){
#                                   #   errors <- character()
#                                   #   mat_cnames = c("i", "j", "val")
#                                   #   if (length(intersect(colnames(object@hic_2d), mat_cnames)) != length(mat_cnames)){
#                                   #     msg <- "colnames of hic_2d must be c(i, j, val)"
#                                   #     errors <- c(errors, msg)
#                                   #   }
#                                   #   reg_cnames = c("seqnames", "start", "end", "index")
#                                   #   if (length(intersect(colnames(object@hic_1d), reg_cnames)) != length(reg_cnames)){
#                                   #     msg <- "colnames of hic_1d must be c(seqnames, start, end, index)"
#                                   #     errors <- c(errors, msg)
#                                   #   }
#                                   #   if (length(errors) == 0) TRUE else errors
#                                   # }
# )
#
# integrated_plot = setClass(Class = "integrated_plot",
#
#                            slots = c(
#                              grobs_list = "list",
#                              ref_dt = "data.table",
#                              parameters = "HiC_parameters",
#                              hic_2d = "data.table",
#                              hic_1d = "data.table"
#
#                            )
#
#                            # validity = function(object){
#                            #   errors <- character()
#                            #   mat_cnames = c("i", "j", "val")
#                            #   if (length(intersect(colnames(object@hic_2d), mat_cnames)) != length(mat_cnames)){
#                            #     msg <- "colnames of hic_2d must be c(i, j, val)"
#                            #     errors <- c(errors, msg)
#                            #   }
#                            #   reg_cnames = c("seqnames", "start", "end", "index")
#                            #   if (length(intersect(colnames(object@hic_1d), reg_cnames)) != length(reg_cnames)){
#                            #     msg <- "colnames of hic_1d must be c(seqnames, start, end, index)"
#                            #     errors <- c(errors, msg)
#                            #   }
#                            #   if (length(errors) == 0) TRUE else errors
#                            # }
# )


###helpers

#
# chr = "chr8"
# start = 25*10^6
# end = 27*10^6
# qgr = GRanges(chr, IRanges(start, end))
#
#
# txdb_gff = makeTxDbFromGFF("gencode.v24.exons.level1.11111.gtf")
# txdb_list = as.list(txdb_gff)
# txdb_list$genes$gene_id = gene_dt[txdb_list$genes$gene_id]$gene_name
# txdb_gn = makeTxDb(transcripts = txdb_list$transcripts, splicings = txdb_list$splicings, genes = txdb_list$genes, chrominfo = txdb_list$chrominfo)
#
#
# p.annot = ggplot(txdb_gn) +
#   geom_alignment(which = qgr,
#                  group.selfish = T, range.geom = "arrowrect",  names.expr = "gene_id") +
#   coord_cartesian(xlim = c(start, end))
# hg38_ideogram = biovizBase::getIdeogram(genome = "hg38")
# p.ideo <- Ideogram(obj = hg38_ideogram, subchr = chr, zoom.region = c(start, end), color = "green", fill = "darkgreen")
#
# p.annot
# p.ideo
# plot_from_data_source(object, qgr, show_max = T)

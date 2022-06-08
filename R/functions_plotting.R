#' Title
#'
#' @param hic_list
#' @param runx_bws
#' @param ctcf_bws
#' @param exons_gr
#' @param max_dist
#' @param hic_names
#' @param target_gene
#' @param extension_size_bp
#' @param output_prefix
#' @param n_arch
#' @param min_arch_dist
#' @param max_fill
#' @param hmap_colors
#'
#' @return
#' @export
#'
#' @examples
plot_combined_hic = function(hic_list, runx_bws, ctcf_bws,
                             exons_gr = NULL, max_dist = NULL,
                             hic_names = names(hic_list),
                             target_gene, extension_size_bp, output_prefix = "combined_plots",
                             n_arch = 50, min_arch_dist = 10*10^5,
                             max_fill = NULL,
                             hmap_colors = c("lightgray", "steelblue", 'darkblue', "red", "red")){
    qgr = fetch_region_by_gene_name(target_gene, ext = extension_size_bp)
    start(qgr) = start(qgr)
    chr = as.character(seqnames(qgr))
    start = start(qgr)
    start = ifelse(start < 0, 0, start)
    end = end(qgr)
    bin_size = hic_list[[1]]@parameters@bin_size
    start = floor(start / bin_size) * bin_size
    end = ceiling(end / bin_size) * bin_size
    if(is.null(exons_gr)){
        p.annot = ggplot(txdb_gn) +
            geom_alignment(which = qgr,
                           group.selfish = T, names.expr = "gene_id") +
            coord_cartesian(xlim = c(start, end))
    }else{
        p.annot = ggplot(exons_gr) +
            geom_alignment(which = qgr,
                           group.selfish = T, names.expr = "gene_id") +
            coord_cartesian(xlim = c(start, end))
    }
    hg38_ideogram = biovizBase::getIdeogram(genome = "hg38")
    p.ideo <- Ideogram(obj = hg38_ideogram, subchr = chr, zoom.region = c(start, end), color = "green", fill = "darkgreen")

    # cells = c("MCF10A", "MCF10A-AT1", "MCF10A-CA1a")

    for(i in 1:length(hic_list)){
        main = hic_names[i]
        if(is.null(max_dist)){
            max_dist = end - start
        }
        p.list = plot_upperMatrix_with_insulation(hic_mat = hic_list[[i]], hmap_colors = hmap_colors, max_fill = max_fill,
                                                  chr, start, end, show_plot = F, max_dist = max_dist)
        # p.list = plot_upperMatrix_with_insulation(hic_list[[i]], chr, start, end, show_plot = F)
        arch_dt = hic_list[[i]][chr, start:end][!is.na(val)]
        arch_start = hic_list[[i]]@hic_1d[arch_dt$i, (start + end) / 2]
        arch_end = hic_list[[i]]@hic_1d[arch_dt$j, (start + end) / 2]
        arch_gr = GRanges(chr, IRanges(arch_start, arch_end))
        arch_gr$value = arch_dt$val
        k = width(arch_gr) > min_arch_dist
        arch_gr = arch_gr[k]
        arch_gr = arch_gr[order(arch_gr$value, decreasing = T)][1:min(n_arch, length(arch_gr))]
        p.arch = ggplot(arch_gr) + geom_arch(aes(height = value)) + coord_cartesian(xlim = c(start, end))
        p.runx = add_bigwig_plot.tiles(runx_bws[i], chr, start, end, bigwig_title = "RUNX1 FE", fe_min = 1, fe_max = 40)
        p.ctcf = add_bigwig_plot.tiles(ctcf_bws[i], chr, start, end, bigwig_title = "CTCF FE", fe_min = 1, fe_max = 80)
        #manual assignment
        p.list = list(ggplotGrob(p.ideo),
                      p.list[[1]],
                      ggplotGrob(p.arch),
                      ggplotGrob(p.annot),
                      p.list[[2]],
                      ggplotGrob(p.runx),
                      ggplotGrob(p.ctcf))
        p.list = lapply(p.list, function(gt){
            gt$layout$clip[gt$layout$name == "panel"] = "off"
            gt
        })

        #get max widths to assign to all
        maxWidth = as.list(grid::unit.pmax(p.list[[1]]$widths, p.list[[2]]$widths,
                                     p.list[[3]]$widths, p.list[[4]]$widths,
                                     p.list[[5]]$widths, p.list[[6]]$widths, p.list[[7]]$widths))
        for(j in 1:length(p.list)){
            p.list[[j]]$widths = maxWidth
        }
        pdf(paste0(output_prefix, "_", target_gene, "_ext", extension_size_bp, "_", main, "_bin", bin_size, ".pdf"), height = 12)
        grid.arrange(arrangeGrob(grobs = p.list, ncol=1,heights=c(.2,.6,.4,.2,.3, .3, .3)), top = main)
        dev.off()
    }

}

#' Title
#'
#' @param hic_mat
#' @param chr
#' @param start
#' @param end
#' @param hmap_colors
#'
#' @return
#' @export
#'
#' @examples
plot_upperMatrix = function(hic_mat, chr, start, end, hmap_colors = c("lightgray", "steelblue", 'darkblue', "red", "red")){
    bin_size = hic_mat@parameters@bin_size
    # hicrng = get_chrRange_Matrix(hic_mat@hic_2d, hic_mat@hic_1d, chr, start, end)
    hicrng = hic_mat[chr, start:end]
    MIN_I = min(c(hicrng$i, hicrng$j))
    XMIN = hic_mat@hic_1d[index == MIN_I]$start
    MAX_I = max(c(hicrng$i, hicrng$j))
    XMAX = hic_mat@hic_1d[index == MAX_I]$end
    hicrng[, x := (i + j) /2]
    hicrng[, y := abs(i - j)]

    hicrng = hicrng[j > i] #limit to upper triangle
    hicrng = na.omit(hicrng) #leave out zero signal regions

    index2center = function(index){
        (index - MIN_I)/(MAX_I - MIN_I) * (XMAX - XMIN - bin_size) + XMIN + bin_size / 2
    }

    hicrng[, x_chr := index2center(x)]
    hicrng[, y_chr := bin_size * y]

    # hicrng[, i_chr := (i - MIN_I)/(MAX_I - MIN_I) * (XMAX - XMIN - bin_size) + XMIN + bin_size / 2]
    # hicrng[, j_chr := (j - MIN_I)/(MAX_I - MIN_I) * (XMAX - XMIN - bin_size) + XMIN + bin_size / 2]

    hex_df = hex_coord_df(hicrng$x_chr, hicrng$y_chr, width = bin_size, height = bin_size, size = 1)
    hex_df = cbind(hex_df, fill=rep(hicrng$val, each=6))

    #pretty up break positions to be multiple of bin_size
    ylab = sort(unique(hex_df$y))
    nlab = 6
    breaks = 0:nlab * ceiling(length(ylab) / nlab) * bin_size
    p = ggplot(hex_df, aes(x=x, y=y)) +
        geom_polygon(aes(group=id, fill=fill)) + scale_fill_gradientn(colours = hmap_colors) +
        labs(x = paste(chr, "position"), y = "distance") + scale_y_continuous(breaks = breaks)
    # p = p + coord_fixed()
    return(p)
    #different plot types facet trick
    #https://statbandit.wordpress.com/2011/07/29/a-ggplot-trick-to-plot-different-plot-types-in-facets/

    # plot_df = hex_df
    # plot_df$type = "hex"
    # f1 <- ggplot(dat2, aes(x=Date,y=value,ymin=0,ymax=value))+facet_grid(variable~., scales='free')
    # ins_gr = GRanges(hic_mat@hic_1d)
    # start(ins_gr) = rowMeans(cbind(start(ins_gr), end(ins_gr)))
    # end(ins_gr)  = start(ins_gr)
    # ins_vals = ins_gr[queryHits(findOverlaps(ins_gr, GRanges(chr, IRanges(start, end))))]$value
}

#' Title
#'
#' @param hic_mat
#' @param chr
#' @param start
#' @param end
#' @param tile_type
#' @param show_plot
#' @param main_title
#' @param max_dist
#' @param point_size
#' @param max_fill
#' @param hmap_colors
#' @param show_insulation_range
#' @param show_minmax
#'
#' @return
#' @export
#'
#' @examples
plot_upperMatrix_with_insulation = function(hic_mat,
                                            chr, start, end,
                                            tile_type = c("diamond", "hex")[1],
                                            show_plot = T, main_title = NULL,
                                            max_dist = 10*10^6,  point_size = 3.5,
                                            max_fill = NULL,
                                            hmap_colors = c("lightgray", "steelblue", 'darkblue', "red", "red"),
                                            show_insulation_range = T, show_minmax = T){
    bin_size = hic_mat@parameters@bin_size
    # hicrng = get_chrRange_Matrix(hic_mat@hic_2d, hic_mat@hic_1d, chr, start, end)
    hicrng = hic_mat[chr, c(start,end)]
    MIN_I = min(c(hicrng$i, hicrng$j))
    XMIN = hic_mat@hic_1d[index == MIN_I]$start
    MAX_I = max(c(hicrng$i, hicrng$j))
    XMAX = hic_mat@hic_1d[index == MAX_I]$end
    # hicrng[, data.table::`:=`(x, (i + j) /2) ]
    # hicrng[, data.table::`:=`(y, abs(i - j))]
    # hicrng[, data.table::`:=`(x = (i + j) /2) ]
    # hicrng[, data.table::`:=`(y = abs(i - j))]
    hicrng[, x := (i + j) /2 ]
    hicrng[, y := abs(i - j)]

    hicrng = hicrng[j > i] #limit to upper triangle
    hicrng = na.omit(hicrng) #leave out zero signal regions

    index2center = function(index){
        (index - MIN_I)/(MAX_I - MIN_I) * (XMAX - XMIN - bin_size) + XMIN + bin_size / 2
    }

    # hicrng[, data.table::`:=`(x_chr = index2center(x))]
    # hicrng[, data.table::`:=`(y_chr = bin_size * y)]
    hicrng[, x_chr := index2center(x)]
    hicrng[, y_chr := bin_size * y]

    # hicrng[, i_chr := (i - MIN_I)/(MAX_I - MIN_I) * (XMAX - XMIN - bin_size) + XMIN + bin_size / 2]
    # hicrng[, j_chr := (j - MIN_I)/(MAX_I - MIN_I) * (XMAX - XMIN - bin_size) + XMIN + bin_size / 2]

    if(tile_type == "hex"){
        hex_df = hex_coord_df(hicrng$x_chr, hicrng$y_chr, width = bin_size, height = bin_size, size = 1)
        hex_df = cbind(hex_df, fill=rep(hicrng$val, each=6))
    }else if(tile_type == "diamond"){
        hex_df = diamond_coord_df(hicrng$x_chr, hicrng$y_chr, width = bin_size, height = bin_size, size = 1)
        hex_df = cbind(hex_df, fill=rep(hicrng$val, each=4))
    }else{
        stop("tile_type must match hex or diamond")
    }
    #pretty up break positions to be multiple of bin_size
    ylab = sort(unique(hex_df$y))
    nlab = 6
    yfac = min(ceiling(max_dist / nlab / bin_size) * bin_size, #when max_dist is less than start to end
               ceiling(length(ylab) / nlab) * bin_size) #when start and end are less than max_dist
    breaks = 0:nlab * yfac
    # p = ggplot(cbind(hex_df, fill=rep(hicrng$val, each=6)), aes(x=x, y=y)) +
    #   geom_polygon(aes(group=id, fill=fill)) + scale_fill_gradientn(colours = hmap_colors) +
    #   labs(x = paste(chr, "position"), y = "distance") + scale_y_continuous(breaks = breaks)
    # p = p + coord_fixed()
    # print(p)
    #different plot types facet trick
    #https://statbandit.wordpress.com/2011/07/29/a-ggplot-trick-to-plot-different-plot-types-in-facets/
    hex_df$type = "hex"
    ins_df = ggplot_hic_delta.df(dt = hic_mat@hic_1d, chr, start, end)


    # plot_df = rbind(hex_df,
    #                 ins_df)
    # plot_df$type = factor(plot_df$type, levels = c("hex", "insulation"))



    add_hex_plot = function(p, ptype = "hex"){
        df = hex_df

        fill_lab = "interaction"
        if(!is.null(max_fill)){
            df$fill[df$fill > max_fill] = max_fill
            fill_lab = paste(fill_lab, "\ncapped at", max_fill)
        }
        df = subset(df, y <= max_dist)

        mmdf = ins_df
        # if(ptype != astype) mmdf$type = astype
        mmdf$minmax = sapply(as.character(mmdf$id), function(x)switch(x, "-1" = "min", "0" = "-", "1" = "max"))
        ann_df = subset(mmdf, id != 0)

        minmax2col = c("min" = "darkgreen", "max" = "orange")
        p = p + geom_polygon(data = df, aes(x=x, y=y, fill = fill, group = id)) +
            scale_fill_gradientn(colours = hmap_colors) +
            labs(x = "", y = "distance", fill = fill_lab) +
            scale_y_continuous(breaks = hic_equal_breaks(bin_size = bin_size, max_dist = max_dist)) +
            coord_cartesian(xlim = c(start, end))
        # if(nrow(ann_df) > 0){
        # ann_df$type = "hex"
        # p = p + annotate("point", x = ann_df$x, y = 0 - max(df$y)*.01, fill = minmax2col[ann_df$minmax], color = "black", size = point_size, stroke = 1.5, shape = 21)
        # }

        if(show_insulation_range){
            p = p + annotate("line", x = c(start, end), y = hic_mat@parameters@n_insulation_bins * hic_mat@parameters@bin_size) +
                annotate("text", x = c(end), y = hic_mat@parameters@n_insulation_bins * hic_mat@parameters@bin_size, label = "insulation\nbin range", hjust = 0)
        }
        return(p)
    }
    add_ins_plot = function(p, ptype = "insulation", astype = ptype){
        # ins_df = ggplot_hic_delta.df(hic_mat@hic_1d, chr, start, end)
        df = subset(ins_df, type == ptype)
        if(ptype != astype) df$type = astype
        df$minmax = sapply(as.character(df$id), function(x)switch(x, "-1" = "min", "0" = "-", "1" = "max"))
        ann_df = subset(df, id != 0)
        # ins_min = -10 #-10 was used for no signal
        df$y[df$y == -10] = NA
        p = p + geom_line(data = subset(df, id == 0), mapping = aes(x = x, y = y)) +
            scale_y_continuous(breaks = hic_equal_breaks(bin_size = bin_size, max_dist = max_dist)) +
            # geom_point(data = ann_df, mapping = aes(shape = rep(21, nrow(ann_df)), x = x, y = y, fill = minmax, color = "black", size = point_size, stroke = 1.5)) +
            scale_fill_manual(values = c("min" = "darkgreen", "max" = "orange")) +
            labs(x = paste(chr, "position"), y = "log2(insulation / mean)") +
            coord_cartesian(xlim = c(start, end)) +
            scale_size_identity() +
            scale_color_identity() +
            scale_shape_identity()
        if(nrow(ann_df) > 0 & show_minmax){
            p = p + geom_point(data = ann_df, mapping = aes(shape = 21, x = x, y = y, fill = minmax, color = "black", size = point_size, stroke = 1.5))
        }
        return(p)
    }
    ###funky facet way
    # p = ggplot()
    # p = add_hex_plot(p)
    # p = add_ins_plot(p)
    # a = p + facet_grid(type~.,scales="free_y")
    # gta = ggplot_gtable(ggplot_build(a))
    # gta$heights[[6]] <- unit(5, "cm")
    # gta$heights[[8]] <- unit(1, "cm")
    # grid.newpage()
    # grid.draw(gta)
    pA = add_hex_plot(ggplot())
    pB = add_ins_plot(ggplot())
    gA = ggplotGrob(pA)
    gB = ggplotGrob(pB)
    maxWidth = grid::unit.pmax(gA$widths, gB$widths)
    gA$widths <- as.list(maxWidth)
    gB$widths <- as.list(maxWidth)
    if(show_plot){
        grid.arrange(arrangeGrob(gA,gB,nrow=2,heights=c(.8,.3)), top = main_title)
    }

    # grid.newpage()
    # grid.draw(rbind(ggplotGrob(a), ggplotGrob(b), size = "last"))
    #
    # arrangeGrob(ggplotGrob(),gB,nrow=2,heights=c(.8,.3))

    invisible(list(ggplots = list(upperMatrix = pA, insulation = pB),
                   grobs = list(upperMatrix = gA, insulation = gB)))
}

equal_breaks <- function(n = 3, s = 0.05, ...){
    function(x){
        if(max(x, na.rm = T) > 1000){
            0:n * min(ceiling(max_dist / n / bin_size) * bin_size, #when max_dist is less than start to end
                      ceiling(max(x, na.rm = T) / n / bin_size) * bin_size, na.rm = T)
        }else{
            c(min(x),0, max(x))
        }
        # rescaling
        # d <- s * diff(range(x)) / (1+2*s)
        # seq(min(x)+d, max(x)-d, length=n)
    }
}

hic_equal_breaks <- function(bin_size, max_dist, n = 3, s = 0.05, ...){
    function(x){
        if(max(x, na.rm = T) > 1000){
            0:n * min(ceiling(max_dist / n / bin_size) * bin_size, #when max_dist is less than start to end
                      ceiling(max(x, na.rm = T) / n / bin_size) * bin_size, na.rm = T)
        }else{
            c(min(x),0, max(x))
        }
        # rescaling
        # d <- s * diff(range(x)) / (1+2*s)
        # seq(min(x)+d, max(x)-d, length=n)
    }
}

ggplotList2grobList = function(ggplotList){
    grobList = lapply(ggplotList, function(x){
        ggplotGrob(x)
    })
    maxWidth = grid::unit.pmax(grobList[[1]]$widths, grobList[[2]]$widths)
    i = 3
    while(i <= length(grobList)){
        maxWidth = grid::unit.pmax(maxWidth, grobList[[3]]$widths)
        i = i + 1
    }

    grobList = lapply(grobList, function(x){
        x$widths = maxWidth
        x
    })
    grobList
}


#' Title
#'
#' @param gene_gr
#' @param qgr
#' @param label_method
#' @param gene_types
#' @param highlight_gene_names
#' @param highlight_gene_color
#' @param olap_extension
#' @param top_spacer
#' @param strand_colors
#' @param arrow_override
#' @param line_size
#' @param line_end
#' @param text_hjust
#' @param text_vjust
#' @param text_angle
#' @param text_size
#' @param text_x_relative
#' @param text_y_relative
#'
#' @return
#' @export
#'
#' @examples
ggplot_ref = function(gene_gr, qgr, label_method = c("text", "label", "none")[1],
                      gene_types = "protein_coding", highlight_gene_names = "",
                      highlight_gene_color = "red",
                      olap_extension = .03,top_spacer = .8, strand_colors = c("+" = "gray0", "-" = "gray40"),
                      arrow_override = arrow(ends = "last", length = unit(.045, "npc")),
                      line_size = 2, line_end = c("round", "butt", "square")[2],
                      text_hjust = 0, text_vjust = -1.5, text_angle = 30, text_size = 3.5,
                      text_x_relative = .5, text_y_relative = .1
){
    pdt = as.data.table(subsetByOverlaps(gene_gr[gene_gr$gene_type %in% gene_types],
                                         qgr,
                                         ignore.strand = TRUE))
    pdt[, y := -1]
    pgr = GRanges(pdt)
    ext = (end(qgr) - start(qgr)) * olap_extension
    egr = pgr
    start(egr) = start(egr) - ext
    end(egr) = end(egr) + ext
    to_decide = 1:nrow(pdt)
    i = 0
    while(length(to_decide) > 0){
        to_keep = numeric()
        to_move = numeric()
        olaps =  as.data.table(findOverlaps(egr[to_decide], egr[to_decide], ignore.strand = T))
        olaps = olaps[queryHits < subjectHits]
        olaps$queryHits = to_decide[olaps$queryHits]
        olaps$subjectHits = to_decide[olaps$subjectHits]
        setkey(olaps, queryHits, subjectHits)
        to_keep = setdiff(to_decide, union(olaps$queryHits, olaps$subjectHits))
        sapply(unique(olaps$queryHits), function(x){
            if(any(x == c(to_keep, to_move))) return()
            to_keep <<- c(to_keep, x)
            to_move <<- union(to_move, olaps[.(x)]$subjectHits)
            olaps <<- olaps[!.(c(to_move, to_keep))]
        })
        pdt[to_keep]$y = i
        i = i + 1
        to_decide = to_move
    }
    pdt[strand == "-", start := end]
    pdt[strand == "-", end := start]
    pdt[y == -1, y := 0]
    pdt[, txt_x := min(start, end) + abs(end - start)*text_x_relative, by = gene_id ]
    pdt[, txt_y := y + text_y_relative, by = gene_id ]
    p = ggplot() +
        geom_segment(data = pdt,
                     aes(x = start, xend = end, y = y, yend = y, col = strand, size = line_size),
                     lineend = line_end,
                     arrow = arrow_override) +
        labs(x = "", y = "") +
        scale_size_identity() +
        scale_color_manual(values = strand_colors) +
        coord_cartesian(xlim = c(start(qgr), end(qgr)),
                        ylim = c(0, max(pdt$y) + top_spacer)) +
        scale_y_continuous(breaks = NULL)
    if(label_method == "text"){
        p = p + annotate("text", label = pdt$gene_name, x = pdt$txt_x, y = pdt$txt_y,
                         color = ifelse(pdt$gene_name %in% highlight_gene_names, highlight_gene_color[1], "black"),
                         hjust = text_hjust, vjust = text_vjust, angle = text_angle, size = text_size)
    }else if(label_method == "label"){
        xpos = pdt[, (start + end) / 2]
        p = p + annotate("label", label = pdt$gene_name, x = pdt$txt_x, y = pdt$txt_y,
                         color = ifelse(pdt$gene_name %in% highlight_gene_names, highlight_gene_color[1], "black"),
                         hjust = .5, vjust = text_vjust, angle = text_angle, size = text_size)
    }
    return(p)
}

ggplot_list = function(my_plots, top_text = "", bottom_text = "position", heights = rep(1, length(my_plots))){
    options(warn = -1)
    my_plots[[length(my_plots)]] = my_plots[[length(my_plots)]] + scale_x_continuous(labels = function(x)paste(x/10^6, "Mb"))
    my_plots[[length(my_plots)]] = my_plots[[length(my_plots)]] + labs(x = bottom_text)

    for(i in 1:(length(my_plots) - 1)){
        if(i < 1) next
        my_plots[[i]] = my_plots[[i]] + scale_x_continuous(labels = NULL)
        my_plots[[i]] = my_plots[[i]] + labs(x = "")
    }
    options(warn = 1)

    my_grobs = lapply(my_plots, function(x){
        ggplotGrob(x)
    })

    #removes clipping like par(xpd = NA or maybe xpd = T)
    my_grobs = lapply(my_grobs, function(gt){
        gt$layout$clip[gt$layout$name == "panel"] = "off"
        gt
    })

    my_widths = lapply(my_grobs, function(gt){
        gt$widths
    })

    # if(exists("maxWidth")) remove(maxWidth)
    maxWidth = my_widths[[1]]
    if(length(my_widths) > 1){
        for(i in 2:length(my_widths)){
            maxWidth = grid::unit.pmax(maxWidth, my_widths[[i]])
        }
    }
    for(j in 1:length(my_grobs)){
        my_grobs[[j]]$widths = maxWidth
    }
    grid.arrange(arrangeGrob(grobs = my_grobs, ncol=1,heights=heights), top = top_text, newpage = F)
}

add_bigwig_plot = function(bigwig_file, chr, start, end, bigwig_title = "FE", p = NULL, fe_max = NULL){
    qgr = GRanges(chr, IRanges(start, end))
    if(is.null(p)) p = ggplot()
    runx_gr = import(con = bigwig_file,
                     format = "BigWig", which = qgr)

    MIN_SCORE = 1
    runx_gr.sub =  subset(runx_gr, score > MIN_SCORE)
    if(is.null(fe_max)) fe_max = max(runx_gr$score)
    runx_gr.bg = setdiff(qgr, runx_gr.sub)
    runx_gr.bg$score = MIN_SCORE

    tmp = as.data.table(c(runx_gr.sub, runx_gr.bg))
    tmp = tmp[order(start)]
    runx_dt = data.table(x = sort(c(tmp$start, tmp$end)),
                         y = rep(tmp$score, each = 2))
    runx_dt = rbind(data.table(x = start(qgr), y = MIN_SCORE), runx_dt, data.table(x = end(qgr), y = MIN_SCORE))

    p2 = ggplot() + labs(x = "", y = bigwig_title) +
        geom_line(data = runx_dt, aes(x = x, y = y)) + coord_cartesian(ylim = c(0, fe_max)) +
        annotate("rect", xmin = start(qgr), xmax = end(qgr), ymin = .9, ymax = MIN_SCORE)
    return(p2)
}

#' add_arch_plot
#'
#' @param hic_mat
#' @param qgr
#' @param src_gr
#' @param min_arch_dist
#' @param n_arch
#' @param p
#'
#' @return
#' @export
#' @importFrom ggbio geom_arch
#' @examples
add_arch_plot = function(hic_mat, qgr, src_gr = qgr, min_arch_dist = 1*10^6, n_arch = 50, p = NULL){
    qidx = subsetByOverlaps(GRanges(hic_mat@hic_1d),
                            qgr,
                            ignore.strand = TRUE)$index
    srcidx = subsetByOverlaps(GRanges(hic_mat@hic_1d),
                              src_gr,
                              ignore.strand = TRUE)$index

    chr = as.character(seqnames(qgr))
    # start = start(qgr) + 1
    # end = end(qgr)
    # arch_dt = hic_mat[chr, start:end][!is.na(val)]
    arch_dt = hic_mat@hic_2d[(i %in% qidx & j %in% srcidx) | (j %in% qidx & i %in% srcidx)]
    arch_dt = arch_dt[val > 0]
    arch_start = hic_mat@hic_1d[arch_dt$i, (start + end) / 2]
    arch_end = hic_mat@hic_1d[arch_dt$j, (start + end) / 2]
    arch_gr = GRanges(chr, IRanges(arch_start, arch_end))
    arch_gr$value = arch_dt$val
    k = width(arch_gr) > min_arch_dist
    arch_gr = arch_gr[k]
    arch_gr = arch_gr[order(arch_gr$value, decreasing = T)][1:min(n_arch, length(arch_gr))]
    if(is.null(p)) p = ggplot()
    p + ggbio::geom_arch(data = arch_gr, aes(height = value)) + coord_cartesian(xlim = c(start(qgr), end(qgr))) +
        labs(y = "interaction")
}

#' Title
#'
#' @param hic_mat
#' @param qgr
#' @param tile_type
#' @param max_dist
#' @param point_size
#' @param text_size
#' @param max_fill
#' @param hmap_colors
#' @param minmax2col
#' @param show_insulation_range
#' @param fill_limits
#' @param show_min
#' @param show_max
#' @param p
#'
#' @return
#' @export
#'
#' @examples
add_matrix_plot = function(hic_mat, qgr,
                           tile_type = c("diamond", "hex")[1],
                           max_dist = 10*10^6,
                           point_size = 3,
                           text_size = 4,
                           max_fill = NULL,
                           hmap_colors = c("lightgray", "steelblue", 'darkblue', "red", "red"),
                           minmax2col = c("min" = "orange", "max" = "green"),
                           show_insulation_range = T,
                           fill_limits = NULL,
                           show_min = T,
                           show_max = F,
                           p = NULL){
    # print("plot from data_source.hic_matrix")
    bin_size = hic_mat@parameters@bin_size
    # hicrng = get_chrRange_Matrix(hic_mat@hic_2d, hic_mat@hic_1d, chr, start, end)
    chr = as.character(seqnames(qgr))
    start = start(qgr)
    end = end(qgr)
    hicrng = hic_mat[chr, c(start,end)]
    MIN_I = min(c(hicrng$i, hicrng$j))
    XMIN = subset(hic_mat@hic_1d, index == MIN_I)$start
    MAX_I = max(c(hicrng$i, hicrng$j))
    XMAX = subset(hic_mat@hic_1d, index == MAX_I)$end
    hicrng[, x := (i + j) /2 ]
    hicrng[, y := abs(i - j)]

    hicrng = hicrng[j > i] #limit to upper triangle
    hicrng = na.omit(hicrng) #leave out zero signal regions

    index2center = function(index){
        (index - MIN_I)/(MAX_I - MIN_I) * (XMAX - XMIN - bin_size) + XMIN + bin_size / 2
    }

    hicrng[, x_chr := index2center(x)]
    hicrng[, y_chr := bin_size * y]

    if(tile_type == "hex"){
        df = hex_coord_df(hicrng$x_chr, hicrng$y_chr, width = bin_size, height = bin_size, size = 1)
        df = cbind(df, fill=rep(hicrng$val, each=6))
    }else if(tile_type == "diamond"){
        df = diamond_coord_df(hicrng$x_chr, hicrng$y_chr, width = bin_size, height = bin_size, size = 1)
        df = cbind(df, fill=rep(hicrng$val, each=4))
    }else{
        stop("tile_type must match hex or diamond")
    }
    #pretty up break positions to be multiple of bin_size
    ylab = sort(unique(df$y))
    nlab = 6
    yfac = min(ceiling(max_dist / nlab / bin_size) * bin_size, #when max_dist is less than start to end
               ceiling(length(ylab) / nlab) * bin_size) #when start and end are less than max_dist
    breaks = 0:nlab * yfac

    fill_lab = "interaction"
    if(!is.null(max_fill)){
        df$fill[df$fill > max_fill] = max_fill
        fill_lab = paste(fill_lab, "\ncapped at", max_fill)
    }
    df = subset(df, y <= max_dist)
    if(is.null(hic_mat@hic_1d$value)){
        show_insulation_range = F
        show_max = F
        show_min = F
        warning("no insulation data supplied")
    }else{
        mmdf = make_hic_minmax_df(dt_1d = hic_mat@hic_1d, qgr)
        mmdf$minmax = sapply(as.character(mmdf$id), function(x)switch(x, "-1" = "min", "0" = "-", "1" = "max"))
        ann_df = subset(mmdf, id != 0)
    }
    if(length(minmax2col) != 2 & !all(names(minmax2col) == c("min", "max"))){
        warning("supplied minmax2col not valid, must be length 2 and have names = c('min', 'max')")
        minmax2col = c("min" = "orange", "max" = "green")
    }

    if(is.null(fill_limits)){
        fill_limits = range(df$fill)
    }else{
        if(length(fill_limits) == 1){
            fill_limits = sort(c(-fill_limits, fill_limits))
        }
        fill_limits = range(fill_limits)
        df$fill[df$fill > max(fill_limits)] = max(fill_limits)
        df$fill[df$fill < min(fill_limits)] = min(fill_limits)
    }

    if(is.null(p)) p = ggplot()
    p = p + geom_polygon(data = df, aes(x=x, y=y, fill = fill, group = id)) +
        scale_fill_gradientn(colours = hmap_colors, limits = fill_limits) +
        labs(x = "", y = "distance", fill = fill_lab) +
        scale_y_continuous(breaks = hic_equal_breaks(bin_size = ifelse(max_dist > 4e6, 1e6, bin_size), max_dist = max_dist)) +
        coord_cartesian(xlim = c(start, end))

    if(exists("ann_df")) if(nrow(ann_df) > 0){
        if(show_max){
            max_df = subset(ann_df, minmax == "max")
            #display triangles along bottom of plot indicating maxima positions
            p = p + annotate("point",
                             x = max_df$x,
                             y = 0 - max(df$y)*.05,
                             fill = minmax2col[max_df$minmax],
                             color = "#00000000",
                             size = point_size,
                             stroke = 1.5,
                             shape = 24)
            #triangle for in plot key in top-right position
            p = p + annotate("point",
                             x = start + (end - start)*.9,
                             y = max(df$y)*.88,
                             fill = minmax2col["max"],
                             color = "#00000000",
                             size = point_size,
                             stroke = 1.5,
                             shape = 24)
            #text for in plot key
            p = p + annotate("text",
                             label = "maxima",
                             x = start + (end - start)*.92,
                             y = max(df$y)*.89,
                             color = minmax2col["max"],
                             size = text_size,
                             hjust = 0,
                             vjust = .5)
        }
        if(show_min){
            min_df = subset(ann_df, minmax == "min")
            #plot triangles along bottom indicating minima positions
            p = p + annotate("point",
                             x = min_df$x,
                             y = 0 - max(df$y)*.03,
                             fill = minmax2col[min_df$minmax],
                             color = "#00000000",
                             size = point_size,
                             stroke = 1.5,
                             shape = 25)
            #triangle for in plot key in top right
            p = p + annotate("point",
                             x = start + (end - start)*.9,
                             y = max(df$y)*.95,
                             fill = minmax2col["min"],
                             color = "#00000000",
                             size = point_size,
                             stroke = 1.5,
                             shape = 25)
            #text for in plot key in top right
            p = p + annotate("text",
                             label = "minima",
                             x = start + (end - start)*.92,
                             y = max(df$y)*.95,
                             color = minmax2col["min"],
                             size = text_size,
                             hjust = 0,
                             vjust = .5)
        }
    }

    if(show_insulation_range){
        #horizontal black line indicating maximum distance considered when calcing insulation
        p = p + annotate("line",
                         x = c(start, end),
                         y = hic_mat@parameters@n_insulation_bins * hic_mat@parameters@bin_size)
        #text label for insulation range line
        p = p + annotate("text",
                         label = "insulation\nbin range",
                         x = c(end),
                         y = hic_mat@parameters@n_insulation_bins * hic_mat@parameters@bin_size,
                         size = text_size,
                         hjust = 1,
                         vjust = -.2)
    }
    return(p)
}

#' sync ggplot widths
#'
#' @param my_plots a list of ggplots
#'
#' @return a list of grobs wtih x margins all equal
#' @import ggplot2
#' @import grid
#' @export
#'
#' @examples
sync_width = function(my_plots){
    stopifnot(class(my_plots) == "list")
    stopifnot(all(sapply(my_plots, function(x)"ggplot" %in% class(x))))
    my_grobs = lapply(my_plots, function(x){
        ggplotGrob(x)
    })

    my_widths = lapply(my_grobs, function(gt){
        gt$widths
    })
    maxWidth = my_widths[[1]]
    if(length(my_widths) > 1){
        for(i in 2:length(my_widths)){
            maxWidth = grid::unit.pmax(maxWidth, my_widths[[i]])
        }
    }
    for(j in 1:length(my_grobs)){
        my_grobs[[j]]$widths = maxWidth
    }
    my_grobs
}

#' sync ggplot heights
#'
#' @param my_plots a list of ggplots
#'
#' @return a list of grobs wtih x margins all equal
#' @import ggplot2
#' @import grid
#' @export
#'
#' @examples
sync_height = function(my_plots){
    stopifnot(class(my_plots) == "list")
    stopifnot(all(sapply(my_plots, function(x)"ggplot" %in% class(x))))
    my_grobs = lapply(my_plots, function(x){
        ggplotGrob(x)
    })

    my_widths = lapply(my_grobs, function(gt){
        gt$heights
    })
    maxWidth = my_widths[[1]]
    if(length(my_widths) > 1){
        for(i in 2:length(my_widths)){
            maxWidth = grid::unit.pmax(maxWidth, my_widths[[i]])
        }
    }
    for(j in 1:length(my_grobs)){
        my_grobs[[j]]$heights = maxWidth
    }
    my_grobs
}


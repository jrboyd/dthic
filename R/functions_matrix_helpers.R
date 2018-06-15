.datatable.aware=TRUE
#these are helper function used by HiC_matrix internally but with possible general usefulness
#these function "do the work" on basic classes.
#HiC_matrix class functions feed appropriate slots to these helpers.

#' Subset data.table by i,j indexes
#'
#'efficiently subset data.table as one would a subrange of a matrix
#'
#'
#' @param dt a data.table containing i, j, and val
#' @param i_range integer, range of i indexes to retrieve
#' @param j_range integer, range of j indexes to retrieve
#' @param fill_square logical, should region spanning below diagonal be filled in?
#'
#' @return subset of dt within ranges
#' @export
#'
#' @examples
subset_dt = function(dt, i_range, j_range, fill_square = FALSE){
    if(!any(class(dt) == "data.table")){
        stop("stop: dt must have class data.table")
    }
    if(!all(key(dt) == c("i", "j"))){
        stop("stop: data.table must be keyed by i and j")
    }

    qi = rep(i_range, length(j_range))
    qj = as.vector(sapply(j_range, function(x){
        rep(x, length(i_range))
    }))
    if(fill_square){
        out = unique(rbind(dt[.(qi, qj)][!is.na(val)],
                           dt[.(qj, qi), .(i = j, j = i, val)][!is.na(val)]))
        out = out[i %in% i_range & j %in% j_range]

    }else{
        out = dt[.(qi, qj)]
    }
    out
}

#' fetch_hic_byGRanges
#'
#' Retrieves rectangular region of a HiC_matrix.
#'
#' @param hic_mat a valid HiC_matrix object
#' @param gr1 GRanges, specifies range of i indexes
#' @param gr2 Granges, specifies range of j indexes, default is same as gr1
#' @param ext integer, extends gr1 and gr2 by set amount, default is 0
#'
#' @return data.table correspoding to gr1xgr2 region of hic_mat
#' @export
#'
#' @examples
fetch_hic_byGRanges = function(hic_mat, gr1, gr2 = gr1, ext = 0){
    stopifnot(any(class(hic_mat) == "HiC_matrix"))
    stopifnot(class(gr1) == "GRanges")
    stopifnot(class(gr2) == "GRanges")
    idx1 = get_chrRange_indexes(hic_mat@hic_1d,
                                chr = as.character(seqnames(gr1)),
                                start(gr1) - ext,
                                end(gr1) + ext)
    idx2 = get_chrRange_indexes(hic_mat@hic_1d,
                                chr = as.character(seqnames(gr2)),
                                start(gr2) - ext,
                                end(gr2) + ext)
    res_dt = subset_dt(dt = hic_mat@hic_2d, j_range = idx1, i_range = idx2, fill_square = TRUE)

    res_dt[, i_start := hic_mat@hic_1d[.(i)]$start]
    res_dt[, i_end := hic_mat@hic_1d[.(i)]$end]
    res_dt[, j_start := hic_mat@hic_1d[.(j)]$start]
    res_dt[, j_end := hic_mat@hic_1d[.(j)]$end]

    res_dt
}

#Found via Dave Shirley's hic package, HiC_functions.R
#Adapted from geotheory package, https://gist.github.com/geotheory/5748388 author: Robin Edwards
#Helper function so not imported
hex_coord_df <- function(x, y, width, height, size = 1) {
    # like hex_coord but returns a dataframe of vertices grouped by an id variable
    dx <- size * width / 6
    dy <- size * height

    hex_x <- rbind(x - 2 * dx, x - dx, x + dx, x + 2 * dx, x + dx, x - dx)
    hex_y <- rbind(y, y + dy, y + dy, y, y - dy, y - dy)
    id    <- rep(1:length(x), each=6)

    data.frame(cbind(x=as.vector(hex_x), y=as.vector(hex_y), id))
}

#adapted from hex_coord_df
diamond_coord_df <- function(x, y, width, height, size = 1) {
    # like hex_coord but returns a dataframe of vertices grouped by an id variable
    dx <- size * width / 4
    dy <- size * height
    #for diamonds, delete 3rd and
    hex_x <- rbind(x - 2 * dx, x, x + 2 * dx, x)
    hex_y <- rbind(y, y + dy, y, y - dy)
    id    <- rep(1:length(x), each=4)

    data.frame(cbind(x=as.vector(hex_x), y=as.vector(hex_y), id))
}

#plot min calls



ggplot_tad_scores = function(tad_calls, regions, chr, start, end, bin_dist_range = c(-Inf, Inf), ins_delta_range = c(-Inf, Inf)){
    pos_indexes = get_chrRange_indexes(regions, chr, start, end)
    ggplot(tad_calls[min_index %in% pos_indexes &
                         bin_dist >= bin_dist_range[1] &
                         bin_dist <= bin_dist_range[2] &
                         ins_delta_mednorm >= ins_delta_range[1] &
                         ins_delta_mednorm <= ins_delta_range[2]], aes(x = mid, y = grp)) + geom_point()
}

#convert a list of data.table to a single tall data.table with added column of list names
list_dt2tall_dt = function(list_dt, new_col = NULL, show_progress = T){
    if(is.null(names(list_dt))) names(list_dt) = paste("group", 1:length(list_dt))
    if(exists("tmp_dt")) remove("tmp_dt")
    for(i in 1:length(list_dt)){
        dt = list_dt[[i]]
        if(!is.null(names(list_dt)) & !is.null(new_col)){
            dt[[new_col]] = names(list_dt)[i]
        }
        if(!exists("tmp_dt")){
            tmp_dt = dt
        }else{
            tmp_dt = rbind(tmp_dt, dt)
        }
        if(show_progress){
            cat('\r',round(i / length(list_dt) * 100, digits = 2), "%", rep(" ", 20))
            flush.console()
        }
    }
    return(tmp_dt)
}


#calculates z-score of log2 transformed data
#z-scores are segmented by chromosome and by distance from diagonal
#returns recalculated HiC_matrix_wInsulation
apply_diagonal_zscore = function(hic_mat){
    test_dt = hic_mat@hic_2d
    test_reg = hic_mat@hic_1d

    #merge to include chromosome information
    print("assembling chromosome info...")
    test_dt = merge(merge(test_dt, test_reg, by.x = "i", by.y = "index"), test_reg, by.x = "j", by.y = "index")
    test_dt[, val := log2(val)]

    MAX = test_dt[seqnames.x == seqnames.y, max(j-i)]
    all_chrms = unique(test_reg$seqnames)
    test_dt[, index_diff := abs(i - j)]
    setkey(test_dt, index_diff, seqnames.x, seqnames.y)

    print("calculating z-scores...")
    test_dt = test_dt[seqnames.y == seqnames.x, ]
    test_dt = test_dt[, .(i, j, val, seqnames = seqnames.x, index_diff)]
    #ninja data.table segmented z-score calculation - kachow
    test_z_dt = test_dt[, .(i, j, val = (val - mean(val, na.rm = T)) / sd(val, na.rm = T)), by = c("index_diff", "seqnames")]

    setkey(test_z_dt, i, j)
    test_z_dt = test_z_dt[!is.na(val),]

    hic_mat@hic_2d = test_z_dt
    hic_mat@hic_1d = hic_mat@hic_1d[,1:4]
    hic_mat@parameters@log2_over_mean_normalization = F
    return(HiC_matrix_wInsulation(hic_mat))
}


decrease_resolution = function(hic_mat, pool_factor, calc_insulation = F){
    bin_size = hic_mat@parameters@bin_size
    newb_size = bin_size * pool_factor
    chr_sizes =  hic_mat@hic_1d[, .(size = max(end)),by = seqnames]
    setkey(chr_sizes, seqnames)
    grs = lapply(chr_sizes$seqnames, function(chr){
        ends = 1:ceiling(chr_sizes[chr]$size / newb_size) * newb_size
        starts = ends - newb_size
        ends = ifelse(ends > chr_sizes[chr]$size, chr_sizes[chr]$size, ends)
        GRanges(chr, IRanges(starts, ends))
    })
    grs = unlist(GRangesList(grs))
    grs$pooled_index = 1:length(grs)
    odt = as.data.table(findOverlaps(grs, GRanges(hic_mat@hic_1d), minoverlap = 5))
    setkey(odt, queryHits)
    print("reduce 1d...")
    new_1d = as.data.table(grs)
    new_1d = new_1d[, .(seqnames, start, end, index = pooled_index)]
    new_1d$seqnames = as.character(new_1d$seqnames)

    print("reduce 2d...")
    nbins = pool_factor^2
    setkey(odt, subjectHits)
    new_2d  = hic_mat@hic_2d
    new_2d[ , pooled_i := odt[.(i), queryHits]]
    new_2d[ , pooled_j := odt[.(j), queryHits]]
    agg_func = function(x){
        sum(x) / nbins
    }
    new_2d = new_2d[, .(val = agg_func(val)), by = c("pooled_i", "pooled_j") ]
    new_2d = new_2d[, .(i = pooled_i, j = pooled_j, val = val)]
    setkey(new_2d, i, j)

    print("assemble HiC_matrix...")
    new_hic = hic_mat
    new_hic@hic_2d = new_2d
    new_hic@hic_1d = new_1d
    new_hic@parameters = HiC_parameters(bin_size = hic_mat@parameters@bin_size * pool_factor)
    if(calc_insulation){
        print("calc insulation...")
        new_hic = HiC_matrix_wInsulation(new_hic)
    }
    return(new_hic)
}

#sparse data.table to dense matrix conversions
mat2dt = function(mat){
    dt2 = as.data.table(mat)
    colnames(dt2) = as.character(1:ncol(dt2))
    dt2$i = 1:nrow(dt2)

    dt2m = melt(dt2, id.vars = "i", variable.name = "j")
    dt2m$j = as.integer(dt2m$j)
    return(dt2m)
}
dt2mat = function(dt){
    adj_i = dt[, min(i)] - 1
    adj_j = dt[, min(j)] - 1
    as.matrix(sparseMatrix(dt$i - adj_i, dt$j - adj_j, x = dt$val))
}


gaussian_cluster_chrm = function(hic_mat, chrm, sigma, nclust, output_prefix = "gaussian_blur", sample_description = paste0("matrix_", hic_mat@parameters@bin_size)){
    idxs = hic_mat@hic_1d[seqnames == chrm]$index
    chrm_dt = hic_mat@hic_2d[i %in% idxs & j %in% idxs]
    chrm_dt[, index_diff := abs(i - j)]
    chrm_dt[, seqnames := chrm]
    chrm_dt = rbind(chrm_dt, chrm_dt[, .(index_diff, seqnames, i = j, j = i, val)])
    picture = dt2mat(chrm_dt)
    picture = ifelse(picture < 0, 0, picture)
    picture2 <- as.matrix(blur(as.im(picture), sigma=sigma))

    png(paste0(output_prefix, "_", sample_description, "_sigma", sigma, "_blur_comparison.png"), width = 1000, height = 1000)
    layout(matrix(c(1:4), nrow=2))
    image.plot(picture, col=gray.colors(50), main="original image", asp=1)
    image.plot(picture2, col=gray.colors(50), main=paste0("blurred with sigma = ", sigma), asp=1)
    drape.plot(1:nrow(picture), 1:ncol(picture), picture, border=NA, theta=0, phi=45, main="original spectrogram")
    drape.plot(1:nrow(picture), 1:ncol(picture), picture2, border=NA, theta=0, phi=45, main=paste0("blurred with sigma = ", sigma))
    dev.off()

    # kpic = kmeans(pic, centers = 6)
    hdist = dist(picture2)^2
    hpic = hclust(d = hdist, method = "average")
    hclusters = cutree(hpic, k = nclust)
    side_cols = rainbow(n = nclust)[hclusters]


    png(paste0(output_prefix, "_", sample_description, "_sigma", sigma, "_nclust", nclust, "_clustering.png"), width = 1000, height = 2000)
    layout(1:2)
    par(xpd = NA)
    image.plot(picture2, col=gray.colors(50))
    sc2 = lapply(sort(unique(side_cols)), function(x){
        mtch = which(x == side_cols)
        mis = which(!mtch[-1] - 1 == mtch[-length(mtch)])
        if(length(mis) == 0)return(data.frame(starts = min(mtch), ends = max(mtch), color = x, stringsAsFactors = F))
        starts = c(min(mtch), mtch[mis+1])
        ends = c(mtch[mis], max(mtch))
        data.frame(starts, ends, color = rep(x, length(starts)), stringsAsFactors = F)
    })

    tmp = sc2[[1]]
    for(i in 2:length(sc2)){
        tmp = rbind(tmp, sc2[[i]])
    }
    tmp$starts = tmp$starts - 1
    tmp$starts = tmp$starts / max(tmp$ends)
    tmp$ends = tmp$ends / max(tmp$ends)
    par(xpd = NA)
    for(i in 1:nrow(tmp)){
        # print(tmp[i,])
        rect(xleft = tmp[i,1], xright = tmp[i,2], ybottom = .99, ytop = 1.03, col = tmp[i,3])
        rect(ybottom = tmp[i,1], ytop = tmp[i,2], xleft = .99, xright = 1.03, col = tmp[i,3])
    }


    image.plot(picture2[hpic$order, hpic$order], col=gray.colors(50))


    # sc = side_cols[hpic$order]
    sc = sapply(sort(unique(side_cols[hpic$order])), function(x){
        range(which(x == side_cols[hpic$order]))
    })
    sc[1,]  = sc[1,] - 1
    sc = sc / max(sc)
    par(xpd = NA)
    for(i in 1:ncol(sc)){
        rect(xleft = sc[1,i], xright = sc[2,i], ybottom = .99, ytop = 1.03, col = colnames(sc)[i])
        rect(ybottom = sc[1,i], ytop = sc[2,i], xleft = .99, xright = 1.03, col = colnames(sc)[i])
    }


    # side_cols = cutree(hpic, k = k)
    # side_cols = rainbow(n = k)[hclusters]
    # sc = side_cols

    dev.off()
}

fetch_bigwig_as_dt = function(bigwig_file, chr, start, end, n_bins = 200, bin_method = c("mean", "max")[2], show_progress_bar = T){
    qgr = GRanges(chr, IRanges(start, end))

    runx_gr = import(con = bigwig_file,
                     format = "BigWig", which = qgr)
    runx_dt = as.data.table(runx_gr)
    tiles_gr = tile(qgr, n = n_bins)[[1]]
    tiles_dt = as.data.table(tiles_gr)
    olaps_dt = as.data.table(findOverlaps(runx_gr, tiles_gr))
    setkey(olaps_dt, "subjectHits")

    if(show_progress_bar){
        app_fun = pbsapply
    }else{
        app_fun = sapply
    }

    tiles_mat = t(app_fun(unique(olaps_dt$subjectHits), function(x){
        hit_dt = runx_dt[olaps_dt[.(x)]$queryHits]
        # hit_gr = runx_gr[olaps_dt[.(x)]$queryHits]
        s = tiles_dt[x]$start
        e = tiles_dt[x]$end
        hit_dt[, start := ifelse(start < s, s, start)]
        hit_dt[, end := ifelse(end > e, e, end)]

        # start(hit_gr)[start(hit_gr) < s] = s
        # end(hit_gr)[end(hit_gr) > e] = e
        mean_score = sum(hit_dt[, (end - start + 1) * score]) / tiles_dt[x, end - start + 1]
        # mean_score = sum(width(hit_gr) * hit_gr$score) / width(tiles_gr[x])
        max_score = max(hit_dt$score)
        min_score = min(hit_dt$score)
        return(c(s, e+1, mean_score, min_score, max_score))
    }))
    colnames(tiles_mat) = c("start", "end", "mean", "min", "max")
    tiles_dt = as.data.table(tiles_mat)
    tiles_dt[, xmin := start]
    tiles_dt[, xmax := end]
    tiles_dt[, x := (xmin + xmax) / 2]
    tiles_dt[, y := get(bin_method)]
    tiles_dt[, c("xmin", "xmax", "x", "y")]
}

add_bigwig_plot.tiles = function(bigwig_files, chr, start, end, n_bins = 200,
                                 bin_method = c("mean", "max")[2], bigwig_title = "FE",
                                 bigwig_colors = NULL, p = NULL, fe_max = NULL, fe_min = NULL,
                                 alpha = NULL, sample_desc = NULL, just_data = F, show_progress_bar = F){
    tiles_dtlist = lapply(bigwig_files, function(x){
        fetch_bigwig_as_dt(x, chr, start, end, n_bins = n_bins, bin_method = bin_method, show_progress_bar = show_progress_bar)
    })
    tiles_dt = list_dt2tall_dt(tiles_dtlist, new_col = "group")
    # tiles_dt = fetch_bigwig_as_dt(bigwig_file, chr, start, end, n_bins = n_bins, bin_method = bin_method)

    if(is.null(fe_max)) fe_max = max(tiles_dt$y)
    if(is.null(fe_min)) fe_min = min(tiles_dt$y)
    if(is.null(alpha)){
        if(length(bigwig_files) > 1){
            alpha = .4
        }else{
            alpha = 1
        }
    }
    # p = p + geom_rect(data = tiles_dt, aes(xmin = xmin, xmax = xmax, ymin= 0, ymax = y, fill = group))
    # p = p + geom_line(data = tiles_dt, aes(x = x, y = y, col = group))
    todo = unique(tiles_dt$group)
    names(todo) = todo
    plot_list = lapply(todo, function(grp){
        grp_dt = tiles_dt[group == grp]
        y1 = grp_dt[nrow(grp_dt),y]
        y2 = 0
        y3 = 0
        y4 = grp_dt[1,y]
        x1 = grp_dt[, max(xmax)]
        x2 = x1
        x3 = grp_dt[, min(xmin)]
        x4 = x3
        main_dt = rbind(grp_dt[, .(x = xmax, y, group)],
                        grp_dt[, .(x = xmin, y, group)]
        )
        main_dt = main_dt[order(x)]
        ending_dt = data.table(x = c(x1,x2,x3,x4), y = c(y1,y2,y3,y4), group = grp)
        rbind(main_dt, ending_dt)
    })
    plot_dt = list_dt2tall_dt(plot_list, new_col = NULL)
    if(is.null(p)){
        p = ggplot() +
            labs(x = "", y = bigwig_title) +
            coord_cartesian(xlim = c(start, end), ylim = c(fe_min, fe_max))

    }

    alpha = substr(rgb(.1, .1, .1, alpha), start = 8, stop = 9)
    if(!is.null(bigwig_colors)){
        group_colors = bigwig_colors
        if(substr(group_colors, 0, 1)[1] != "#"){
            group_colors = col2hex(group_colors)
        }
    }else{
        len = length(unique(tiles_dt$group))
        group_colors = RColorBrewer::brewer.pal(max(len, 3), "Dark2")[1:len]
    }

    group_colors = paste0(group_colors, alpha)
    names(group_colors) = unique(tiles_dt$group)
    if(just_data) return(plot_dt)
    p = p + geom_polygon(data = plot_dt, aes(x = x, y = y, col = NULL, fill = group)) +
        labs(fill = sample_desc) +
        scale_fill_manual(values = group_colors) +
        scale_color_manual(values = group_colors)
    return(p)
}

fetch_region_by_gene_name = function(gene_names, ext = 0){
    # if(!exists("gene_ref")){
    #   load("../HiC_histone_genes/hg38ref.save")
    #   ref_dt = as.data.table(ref)
    #   ref_dt = ref_dt[!duplicated(gene_name)]
    #   ref_dt = ref_dt[, .(seqnames = chrm, start, end, strand, symbol = gene_name, ensembl_id = gene_id)]
    #   setkey(ref_dt, symbol)
    #   genesymbol = GRanges(ref_dt)
    #   names(genesymbol) = genesymbol$symbol
    #   gene_ref <<- genesymbol
    # }
    # rng = range(gene_ref[gene_names], ignore.strand = TRUE)
    rng = range(gene_gr[gene_names], ignore.strand = TRUE)
    start(rng) = start(rng) - ext
    end(rng) = end(rng) + ext
    return(rng)
}

fetch_genes_in_region = function(qgr = NULL, chr = NULL, start = NULL, end = NULL, target_gr = NULL){
    if(is.null(qgr)){
        if(is.null(chr) | is.null(start) | is.null(end)){
            stop("need qgr GRanges or chr, start, and end")
        }
        qgr = GRanges(chr, IRanges(start, end))
    }
    if(is.null(target_gr)){
        target_gr = exon_gr
    }
    target_gr[queryHits(findOverlaps(query = target_gr, qgr))]

}

# plot_combined_hic = function(hic_list, hic_names, target_gene, extension_size_bp, file_prefix = "combined_plots"){
#   tp_hic_list = my_p_hicsdz_capped
#
#   qgr = fetch_region_by_gene_name("ZEB1", ext = 2*10^6)
#   start(qgr) = start(qgr)
#   chr = as.character(seqnames(qgr))
#
#   start = start(qgr)
#   end = end(qgr)
#   bin_size = tp_hic_list[[1]]@parameters@bin_size
#   start = floor(start / bin_size) * bin_size
#   end = ceiling(end / bin_size) * bin_size
#
#   p.annot = ggplot(txdb_gn) +
#     geom_alignment(which = qgr,
#                    group.selfish = T, range.geom = "arrowrect",  names.expr = "gene_id") +
#     coord_cartesian(xlim = c(start, end))
#   hg38_ideogram = biovizBase::getIdeogram(genome = "hg38")
#   p.ideo <- Ideogram(obj = hg38_ideogram, subchr = chr, zoom.region = c(start, end), color = "green", fill = "darkgreen")
#
#   cells = c("MCF10A", "MCF10A-AT1", "MCF10A-CA1a")
#   for(i in 1:3){
#     main = cells[i]
#     p.list = plot_upperMatrix_with_insulation(tp_hic_list[[i]], chr, start - 1*10^6, end + 1*10^6, show_plot = F)
#     # p.list = plot_upperMatrix_with_insulation(tp_hic_list[[i]], chr, start, end, show_plot = F)
#     p.runx = add_bigwig_plot(runx_bw[i], chr, start, end, bigwig_title = "RUNX1 FE", fe_max = 30)
#     p.ctcf = add_bigwig_plot(ctcf_bw[i], chr, start, end, bigwig_title = "CTCF FE", fe_max = 30)
#     #manual assignment
#     p.list = list(ggplotGrob(p.ideo),
#                   p.list[[1]],
#                   ggplotGrob(p.annot),
#                   p.list[[2]],
#                   ggplotGrob(p.runx),
#                   ggplotGrob(p.ctcf))
#     #get max widths to assign to all
#     maxWidth = as.list(unit.pmax(p.list[[1]]$widths, p.list[[2]]$widths,
#                                  p.list[[3]]$widths, p.list[[4]]$widths,
#                                  p.list[[5]]$widths, p.list[[6]]$widths))
#     for(j in 1:length(p.list)){
#       p.list[[j]]$widths = maxWidth
#     }
#     pdf(paste0("plots_ZEB1_tss_combined_", main, "_", bin_size, ".pdf"), height = 12)
#     grid.arrange(arrangeGrob(grobs = p.list, ncol=1,heights=c(.2,.6,.6,.3, .3, .3)), top = main)
#     dev.off()
#   }
# }

#combines insulation profiles and performs quantile normalization
quant_norm_insulation = function(my_hicsd){
    my_hicsd_tad_calls = assemble_tad_calls(my_hicsd)
    ggplot_tad_scores(my_hicsd_tad_calls, my_hicsd[[1]]@hic_1d, "chr2", 40*10^6, 55*10^6, bin_dist_range = c(0,20), ins_delta_range = c(2,Inf))

    comb_ins_dt = list_dt2tall_dt(lapply(my_hicsd, function(x)x@hic_1d))
    comb_ins_dt[, mid := (start + end) / 2]

    dc_comb_ins_dt = dcast(comb_ins_dt[, c("index", "value", "group", "mid")],
                           formula = index + mid ~ group,
                           value.var = "value")
    val_cn = colnames(dc_comb_ins_dt)[-1:-2]
    raw_mat = as.matrix(dc_comb_ins_dt[,val_cn, with = F])

    normq_mat = preprocessCore::normalize.quantiles(raw_mat)
    for(i in 1:ncol(normq_mat)){
        dc_comb_ins_dt[,val_cn[i]] = normq_mat[,i]
    }
    return(list(normq_mat, dc_comb_ins_dt))
}





#creates data.frame containing min and max info for insulation profiles
make_hic_minmax_df = function(dt_1d, qgr){
    chr = as.character(seqnames(qgr))
    start = start(qgr) + 1
    end = end(qgr)

    gr = GRanges(dt_1d)
    start(gr) = start(gr) + 1
    q_gr = GRanges(chr, IRanges(start, end))
    q_indexes = dt_1d[subjectHits(findOverlaps(query = q_gr, subject = gr))]$index
    dt_1d = dt_1d[q_indexes]
    dt_1d[, xs := (start + end) / 2]
    max_k = dt_1d$minmax == 1
    max_k[is.na(max_k)] = F
    min_k = dt_1d$minmax == -1
    min_k[is.na(min_k)] = F
    insulation_df = data.frame(x = dt_1d$xs, y = dt_1d$value, id = 0, type = "insulation")
    if(any(max_k)){
        insulation_df = rbind(insulation_df,
                              data.frame(x = dt_1d$xs[max_k], y = dt_1d$value[max_k], id = 1, type = "insulation"))
    }
    if(any(min_k)){
        insulation_df = rbind(insulation_df,
                              data.frame(x = dt_1d$xs[min_k], y = dt_1d$value[min_k], id = -1, type = "insulation"))
    }
    return(insulation_df)
}

library(rtracklayer)
library(pbapply)
#windowed view of bigwig file, filtered by qgr if supplied
fetch_windowed_bw = function(bw_file, win_size = 50, qgr = NULL){
    if(is.null(qgr)){
        bw_gr = import.bw(bw_file)
    }else{
        bw_gr = import.bw(bw_file, which = qgr)
    }

    rng = range(bw_gr)
    start(rng) = start(rng) - start(rng) %% win_size + 1
    end(rng) = end(rng) - end(rng) %% win_size + win_size
    # end(rng) = ceiling(end(rng)/win_size) * win_size
    win = slidingWindows(rng, width = win_size, step = win_size)
    print(object.size(win), units = "GB")
    print(object.size(bw_gr), units = "GB")
    win = unlist(win)
    sn = names(seqlengths(win))
    names(sn) = sn
    suppressWarnings({
        new_seqlengths = pbsapply(sn, function(x){

            max(end(subset(win, seqnames == x)))

        })
        seqlengths(win) = new_seqlengths
        seqlengths(bw_gr) = seqlengths(win)
    })
    mid_gr = function(gr){
        start(gr) + floor((width(gr) - 1)/2)
    }
    mids = mid_gr(win)
    start(win) = mids
    end(win) = mids
    olaps = findOverlaps(win, bw_gr)
    win = win[queryHits(olaps)]
    win$FE = bw_gr[subjectHits(olaps)]$score
    return(win)
}

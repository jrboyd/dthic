score_tadness = function(dt){
    if(exists("score_dt")) remove(score_dt, pos = ".GlobalEnv")
    hidden = pbapply::pblapply(unique(dt$seqnames), function(chr){
        chr_dt = dt[seqnames == chr] # & value > -10] #-10 indicator values removed
        min_dt = chr_dt[minmax == -1 & !is.na(delta) & value < 0] # minima must be below 0 value
        hidden =  lapply(min_dt$index, function(ind){
            min_val = chr_dt[index == ind]$value

            right_ind = chr_dt[ minmax == 1 & index > ind][order(index, decreasing = F)[1]]$index
            right_delta = chr_dt[index == right_ind]$value - min_val
            right_dist = abs(ind - right_ind)

            left_ind = chr_dt[ minmax == 1 & index < ind][order(index, decreasing = T)[1]]$index
            left_delta = chr_dt[index == left_ind]$value - min_val
            left_dist = abs(ind - left_ind)

            k = !is.na(c(left_dist, right_dist))
            if(sum(k) > 0){
                new_dt = data.table(min_index = rep(ind, 2)[k],
                                    max_index = c(left_ind, right_ind)[k],
                                    direction = c("left", "right")[k],
                                    bin_dist = c(left_dist, right_dist)[k],
                                    ins_delta = c(left_delta, right_delta)[k])
                if(!exists("score_dt")){
                    score_dt <<- new_dt
                }else{
                    score_dt <<- rbind(score_dt, new_dt)
                }
            }
            return(NULL)
        })
    })
    return(score_dt)
}

#taking 1 entry as baseline for tad call, assembles data.table of minmax deltas for full hic_list
#for each minmax pair in baseline, measure delta at those genomic coordinates
compare_tadness = function(hic_list, baseline = 1){
    hic_scores = score_tadness(hic_list[[baseline]]@hic_1d)
    # names(hic_list)[baseline]
    for(i in 1:length(hic_list)){
        cn = names(hic_list)[i]
        if(i == baseline){
            hic_scores[[cn]] = hic_scores$ins_delta
        }else{
            ins_dt = hic_list[[i]]@hic_1d
            hic_scores[[cn]] = ins_dt[hic_scores$max_index]$value - ins_dt[hic_scores$min_index]$value
        }
    }
    return(hic_scores)
}

compare_tadness.plot_prep = function(hic_ins_list, baseline = 1){
    print(paste("using", names(hic_ins_list)[baseline], "as reference"))
    ins_dt = hic_ins_list[[baseline]]@hic_1d
    tad_ins_deltas = compare_tadness(hic_ins_list, baseline = baseline)
    # tad_ins_deltas[, 6:11]
    tad_ins_deltas = cbind(tad_ins_deltas, ins_dt[tad_ins_deltas$min_index, .(seqnames = seqnames, mid_min = (start + end) / 2)])
    tad_ins_deltas = cbind(tad_ins_deltas, ins_dt[tad_ins_deltas$max_index, .(mid_max = (start + end) / 2)])

    # all_valid = apply(as.matrix(tad_ins_deltas[,6:11]), 1, function(x)!any(is.na(x)))
    quant_deltas = preprocessCore::normalize.quantiles(as.matrix(tad_ins_deltas[,6:(5+length(hic_ins_list))]))
    colnames(quant_deltas) = paste0(colnames(tad_ins_deltas)[6:(5+length(hic_ins_list))], "_qnorm")
    tad_ins_deltas = cbind(tad_ins_deltas, quant_deltas)
}

plot_tad_schematic = function(hic_ins_list, tad_ins_deltas, chr, start, end, baseline, delta_cutoff = .6){
    pos_indexes = get_chrRange_indexes(hic_ins_list[[1]]@hic_1d, chr, start - 10^6, end + 10^6)
    ins_yrng = range(hic_ins_list[[baseline]]@hic_1d[pos_indexes, ]$value, na.rm = T)
    b = -min(ins_yrng)
    m = 1 / (max(ins_yrng) - min(ins_yrng))
    plot_deltas = tad_ins_deltas[min_index %in% pos_indexes]

    plot(0, xlim = c(start, end), ylim = c(-1, 1))
    tp = which(grepl("qnorm", colnames(plot_deltas)))
    delta_max = max(plot_deltas[, tp, with = F], na.rm = T)
    grps = sapply(strsplit(colnames(plot_deltas)[tp], "_"), function(x)x[1])
    rcols = RColorBrewer::brewer.pal(length(unique(grps)), "Dark2")
    names(rcols) = unique(grps)
    rcols = rcols[grps]

    for(i in 1:nrow(plot_deltas)){

        x_s = plot_deltas[i, ]$mid_min
        x_e = plot_deltas[i, ]$mid_max
        # y_s = 0
        # y_e = plot_deltas[i, ]$ins_delta
        y_s = (hic_ins_list[[baseline]]@hic_1d[plot_deltas[i, ]$min_index]$value + b) * m
        y_e = (hic_ins_list[[baseline]]@hic_1d[plot_deltas[i, ]$max_index]$value + b) * m
        pol_col = rcols[baseline]
        fill_col = ifelse(plot_deltas[i, ]$ins_delta >= delta_cutoff, rcols[baseline], "white")
        polygon(c(x_s, x_e, x_e, x_s), c(y_s, y_s, y_e, y_s),col = fill_col)

        #do not indclude barplot for weak deltas
        if(plot_deltas[i, ]$ins_delta < delta_cutoff) next
        center = mean(c(x_s, x_e))
        yvals = plot_deltas[i, tp, with = F] / delta_max * .6
        bp_w = (end - start) * .025
        box_min = min(c(x_s, x_e)) #center - bp_w
        box_max = max(c(x_s, x_e)) #center + bp_w
        box_win = (box_max - box_min) / length(yvals)

        rect(xleft = box_min,
             ybottom = -.8,
             ytop = -.8 + .6,
             xright = box_min + length(yvals) * box_win, col = "white")

        rect(xleft = box_min + (1:length(yvals) - 1) * box_win,
             ybottom = -.8,
             ytop = -.8 + yvals,
             xright = box_min + 1:length(yvals) * box_win, col = rcols)

    }
}


###assemble tall data.table of tad calls
assemble_tad_calls = function(hic_ins_list){
    tad_scores = lapply(hic_ins_list, function(x) score_tadness(x@hic_1d))
    tad_scores_dt = tad_scores[[1]]
    tad_scores_dt$grp = names(tad_scores)[1]
    for(i in 2:length(tad_scores)){
        new_dt = tad_scores[[i]]
        new_dt$grp = names(tad_scores)[i]
        tad_scores_dt = rbind(tad_scores_dt, new_dt)
    }
    tad_scores_dt = cbind(tad_scores_dt, hic_ins_list[[1]]@hic_1d[tad_scores_dt$min_index, .(seqnames, mid = (start + end) / 2)])
    tad_scores_dt[, ins_delta_mednorm := ins_delta / quantile(ins_delta, .5, na.rm = T), by = grp]
    return(tad_scores_dt)
}

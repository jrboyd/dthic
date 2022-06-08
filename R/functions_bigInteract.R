# library(data.table)
# query_gr = GRanges("chr19:11326847-14093246")
# fibed = "~/R/hic_analysis/data/chicago4/CHi_AF_03122019_MCF10A/fullOutput.ibed"
# fhicpro = "~/HiC-Pro/outputs/AF_dpnII_03112019/MCF10A_full/hic_results/matrix/rep4/raw/50000/rep4_50000.matrix"

cn_interact = c(
    "seqnames",
    "start",
    "end",
    "name",
    "score",
    "value",
    "exp",
    "color",
    "sourceChrom",
    "sourceStart",
    "sourceEnd",
    "sourceName",
    "sourceStrand",
    "targetChrom",
    "targetStart",
    "targetEnd",
    "targetName",
    "targetStrand"
)

base_cmd = "bigBedToBed FETCHBI_FILE /dev/stdout -chrom=FETCHBI_CHR -start=FETCHBI_START -end=FETCHBI_END"

#' fetchBigInteract
#'
#' @param bif
#' @param qgr
#' @param qgr2
#' @param silent if TRUE, no warnings will be shown. default is FALSE.
#'
#' @return
#' @export
#'
#' @examples
fetchBigInteract = function(bif, qgr, qgr2 = NULL, silent = FALSE){
    chr = as.character(seqnames(qgr))
    s = start(qgr)
    e = end(qgr)
    cmd = base_cmd
    cmd = sub("FETCHBI_FILE", bif, cmd)
    cmd = sub("FETCHBI_CHR", chr, cmd)
    cmd = sub("FETCHBI_START", s, cmd)
    cmd = sub("FETCHBI_END", e, cmd)
    dt = data.table(str =  system(cmd, intern = TRUE))
    dt = dt[, tstrsplit(str, "\t")]
    if(ncol(dt) != 18){
        if(!silent) warning(as.character(qgr), " of ", bif, "was empty.")
        dt = as.data.table(matrix("", nrow = 0, ncol = 18))
    }
    colnames(dt) = cn_interact
    dt$start = as.numeric(dt$start)
    dt$end = as.numeric(dt$end)
    dt$score = as.numeric(dt$score)
    dt$value = as.numeric(dt$value)
    dt$sourceStart  = as.numeric(dt$sourceStart)
    dt$sourceEnd = as.numeric(dt$sourceEnd)
    dt$targetStart = as.numeric(dt$targetStart)
    dt$targetEnd = as.numeric(dt$targetEnd)
    if(!is.null(qgr2)){
        dt = dt[targetChrom == as.character(seqnames(qgr2)) & targetEnd >= start(qgr2) & targetStart <= end(qgr2)]
    }
    dt
}

#' writeBigInteract_ibed
#'
#' @param fibed
#' @param group_name
#' @param output_root
#' @param cols
#' @param index_method
#' @param remove_trans
#' @param query_gr
#' @param max_dist
#' @param remove_bait2bait
#' @param remove_bait2non
#' @param min_N
#' @param min_score
#' @param cleanup_interact
#' @param chr_sizes
#' @param autosql
#'
#' @return
#' @export
#'
#' @examples
writeBigInteract_ibed = function(fibed,
                                 #output options
                                 group_name = NULL,
                                 output_root = NULL,
                                 cols = seqsetvis::safeBrew(5),
                                 index_method = c("a", "b", "ab")[1],
                                 #filtering
                                 remove_trans = TRUE,
                                 query_gr = NULL,
                                 max_dist = Inf,
                                 remove_bait2bait = TRUE,
                                 remove_bait2non = FALSE,
                                 min_N = 0,
                                 min_score = 5,
                                 #cleanup
                                 cleanup_interact = TRUE,
                                 #supporting files
                                 chr_sizes = "/slipstream/home/joeboyd/hg38_chrsizes.txt",
                                 autosql = "/slipstream/home/joeboyd/autoSql_bigInteract.txt"

){
    if(remove_bait2bait & remove_bait2non){
        stop("can't remove_bait2bait and remove_bait2non. At least one must be FALSE.")
    }
    if(is.null(group_name)){
        group_name = sub(".ibed", "", basename(fibed))
    }
    if(is.null(output_root)){
        output_root = file.path(dirname(fhicpro), group_name)
    }
    fout = paste0(output_root, ".interact")
    fbig = paste0(output_root, ".bigInteract")
    message("preparing ", fbig)
    idt = fread(fibed)
    nr = nrow(idt)
    bn = unique(idt$bait_name)
    bn_cols = cols[(seq_along(bn) %% length(cols)) +1]
    names(bn_cols) = bn
    message(round(nrow(idr)/1e6, 1), " M pairs")
    if(remove_trans){
        #remove trans
        idt = idt[bait_chr == otherEnd_chr]
        message(round(nrow(idt)/nr*100, 6), "% cis")
    }
    if(!is.null(query_gr)){
        igr = GRanges(idt[, .(seqnames = bait_chr, start = bait_start, end = bait_end)])
        idt = idt[subjectHits(findOverlaps(query_gr, igr, ignore.strand = TRUE)),]
        message(round(nrow(idt)/nr*100, 6), "% filtered for genes")
    }
    if(!is.infinite(max_dist)){
        idt[, dist := round((otherEnd_end + otherEnd_start)/2 - (bait_end  + bait_start)/2), ]
        idt = idt[dist < max_dist]
        message(round(nrow(idt)/nr*100, 2), "% below max_dist ", max_dist)
    }
    if(remove_bait2bait){
        idt = idt[otherEnd_name == "."]
        message(round(nrow(idt)/nr*100, 6), "% bait to non-bait")
    }
    if(remove_bait2non){
        idt = idt[otherEnd_name != "."]
        message(round(nrow(idt)/nr*100, 6), "% bait to bait")
    }
    if(min_N > 0){
        idt = idt[N_reads > min_N]
        message(round(nrow(idt)/nr*100, 6), "% over min_N ", min_N)
    }
    if(min_score > 0){
        idt = idt[score > min_score]
        message(round(nrow(idt)/nr*100, 6), "% over min_score ", min_score)
    }

    switch(index_method,
           a = {
               interact = idt[, .(chrom = bait_chr,
                                  chromStart = bait_start,
                                  chromEnd = bait_end,
                                  name = bait_name,
                                  score = N_reads,
                                  value = score,
                                  exp = group_name,
                                  color = bn_cols[bait_name],
                                  sourceChrom = bait_chr,
                                  sourceStart = bait_start,
                                  sourceEnd = bait_end,
                                  sourceName = bait_name,
                                  sourceStrand = ".",
                                  targetChrom = otherEnd_chr,
                                  targetStart = otherEnd_start,
                                  targetEnd = otherEnd_end,
                                  targetName = paste0(otherEnd_chr, ":", otherEnd_start, "-", otherEnd_end),
                                  targetStrand = "."
               )]
           },
           b = {
               interact = idt[, .(chrom = otherEnd_chr,
                                  chromStart = otherEnd_start,
                                  chromEnd = otherEnd_end,
                                  name = bait_name,
                                  score = N_reads,
                                  value = score,
                                  exp = group_name,
                                  color = bn_cols[bait_name],
                                  sourceChrom = bait_chr,
                                  sourceStart = bait_start,
                                  sourceEnd = bait_end,
                                  sourceName = bait_name,
                                  sourceStrand = ".",
                                  targetChrom = otherEnd_chr,
                                  targetStart = otherEnd_start,
                                  targetEnd = otherEnd_end,
                                  targetName = paste0(otherEnd_chr, ":", otherEnd_start, "-", otherEnd_end),
                                  targetStrand = "."
               )]
           },
           ab = {
               idt_cis = idt[bait_chr == otherEnd_chr]
               idt_trans = idt[bait_chr != otherEnd_chr]
               interact =
                   rbind(
                       idt_cis[, .(chrom = bait_chr,
                                   chromStart = pmin(bait_start, otherEnd_start),
                                   chromEnd = pmax(bait_end, otherEnd_end),
                                   name = bait_name,
                                   score = N_reads,
                                   value = score,
                                   exp = group_name,
                                   color = bn_cols[bait_name],
                                   sourceChrom = bait_chr,
                                   sourceStart = bait_start,
                                   sourceEnd = bait_end,
                                   sourceName = bait_name,
                                   sourceStrand = ".",
                                   targetChrom = otherEnd_chr,
                                   targetStart = otherEnd_start,
                                   targetEnd = otherEnd_end,
                                   targetName = paste0(otherEnd_chr, ":", otherEnd_start, "-", otherEnd_end),
                                   targetStrand = "."
                       )],
                       idt_trans[, .(chrom = bait_chr,
                                     chromStart = bait_start,
                                     chromEnd = bait_end,
                                     name = bait_name,
                                     score = N_reads,
                                     value = score,
                                     exp = group_name,
                                     color = bn_cols[bait_name],
                                     sourceChrom = bait_chr,
                                     sourceStart = bait_start,
                                     sourceEnd = bait_end,
                                     sourceName = bait_name,
                                     sourceStrand = ".",
                                     targetChrom = otherEnd_chr,
                                     targetStart = otherEnd_start,
                                     targetEnd = otherEnd_end,
                                     targetName = paste0(otherEnd_chr, ":", otherEnd_start, "-", otherEnd_end),
                                     targetStrand = "."
                       )],
                       idt_trans[, .(chrom = otherEnd_chr,
                                     chromStart = otherEnd_start,
                                     chromEnd = otherEnd_end,
                                     name = bait_name,
                                     score = N_reads,
                                     value = score,
                                     exp = group_name,
                                     color = bn_cols[bait_name],
                                     sourceChrom = bait_chr,
                                     sourceStart = bait_start,
                                     sourceEnd = bait_end,
                                     sourceName = bait_name,
                                     sourceStrand = ".",
                                     targetChrom = otherEnd_chr,
                                     targetStart = otherEnd_start,
                                     targetEnd = otherEnd_end,
                                     targetName = paste0(otherEnd_chr, ":", otherEnd_start, "-", otherEnd_end),
                                     targetStrand = "."

                       )]
                   )
           }
    )
    interact = interact[order(chromStart)][order(chrom)]

    fwrite(interact, fout, sep = "\t", col.names = FALSE)

    system(paste0("bedToBigBed ", fout, " ", chr_sizes, " ", fbig, " -type=bed5+13 -as=", autosql))

    if(cleanup_interact){
        file.remove(fout)
    }

    fbig
}

#' writeBigInteract_hicpro
#'
#' @param fhicpro
#' @param bed_file
#' @param group_name
#' @param output_root
#' @param cols
#' @param index_method
#' @param remove_trans
#' @param query_gr
#' @param max_dist
#' @param min_value
#' @param cleanup_interact
#' @param chr_sizes
#' @param autosql
#'
#' @return
#' @export
#'
#' @examples
writeBigInteract_hicpro = function(fhicpro,
                                   bed_file = sub(".matrix$", "_abs.bed", fhicpro),
                                   #output options
                                   group_name = NULL,
                                   output_root = NULL,
                                   cols = seqsetvis::safeBrew(1),
                                   index_method = c("a", "b", "ab")[1],
                                   #filtering
                                   remove_trans = TRUE,
                                   query_gr = NULL,
                                   max_dist = Inf,
                                   min_value = 0,
                                   #cleanup
                                   cleanup_interact = TRUE,
                                   #supporting files
                                   chr_sizes = "/slipstream/home/joeboyd/hg38_chrsizes.txt",
                                   autosql = "/slipstream/home/joeboyd/autoSql_bigInteract.txt"

){
    stopifnot(file.exists(fhicpro))
    stopifnot(file.exists(bed_file))
    if(is.null(group_name)){
        group_name = sub(".matrix", "", basename(fhicpro))
    }
    if(is.null(output_root)){
        output_root = file.path(dirname(fhicpro), group_name)
    }
    fout = paste0(output_root, ".interact")
    fbig = paste0(output_root, ".bigInteract")
    message("preparing ", fbig)
    idt = fread(fhicpro)
    colnames(idt) = c("i", "j", "value")
    bdt = fread(bed_file)
    colnames(bdt) = c("seqnames", "start", "end", "index")
    idt = merge(idt, bdt, by.x = "i", by.y = "index")
    idt = merge(idt, bdt, by.x = "j", by.y = "index")
    colnames(idt) = sub("\\.x$", "_i", colnames(idt))
    colnames(idt) = sub("\\.y$", "_j", colnames(idt))
    nr = nrow(idt)
    message(round(nrow(idr)/1e6, 1), " M pairs")
    if(remove_trans){
        #remove trans
        idt = idt[seqnames_i == seqnames_j]
        message(round(nrow(idt)/nr*100, 6), "% cis")
    }
    idt = rbind(idt, idt[i != j, .(j = i, i = j, value,
                                   seqnames_i = seqnames_j,
                                   start_i = start_j, end_i = end_j,
                                   seqnames_j = seqnames_i,
                                   start_j = start_i, end_j = end_i)])
    colnames(idt)[4:6] = c("bait_chr", "bait_start", "bait_end")
    colnames(idt)[7:9] = c("otherEnd_chr", "otherEnd_start", "otherEnd_end")
    nr = nrow(idt)
    if(!is.null(query_gr)){
        igr = GRanges(idt[, .(seqnames = bait_chr, start = bait_start, end = bait_end)])
        idt = idt[subjectHits(findOverlaps(query_gr, igr, ignore.strand = TRUE)),]
        message(round(nrow(idt)/nr*100, 6), "% filtered for genes")
    }
    if(!is.infinite(max_dist)){
        idt[, dist := round((otherEnd_end + otherEnd_start)/2 - (bait_end  + bait_start)/2), ]
        idt = idt[dist < max_dist]
        message(round(nrow(idt)/nr*100, 2), "% below max_dist ", max_dist)
    }
    if(min_value > 0){
        idt = idt[value > min_value]
        message(round(nrow(idt)/nr*100, 6), "% over min_value ", min_value)
    }

    interact = idt[, .(chrom = bait_chr,
                       chromStart = bait_start,
                       chromEnd = bait_end,
                       name = paste(i,j, sep = "_"),
                       score = value,
                       value = value,
                       exp = group_name,
                       color = cols[1],
                       sourceChrom = bait_chr,
                       sourceStart = bait_start,
                       sourceEnd = bait_end,
                       sourceName = paste0(bait_chr, ":", bait_start, "-", bait_end),
                       sourceStrand = ".",
                       targetChrom = otherEnd_chr,
                       targetStart = otherEnd_start,
                       targetEnd = otherEnd_end,
                       targetName = paste0(otherEnd_chr, ":", otherEnd_start, "-", otherEnd_end),
                       targetStrand = "."
    )]
    interact = interact[order(chromStart)][order(chrom)]

    fwrite(interact, fout, sep = "\t", col.names = FALSE)

    system(paste0("bedToBigBed ", fout, " ", chr_sizes, " ", fbig, " -type=bed5+13 -as=", autosql))

    if(cleanup_interact){
        file.remove(fout)
    }

    fbig
}

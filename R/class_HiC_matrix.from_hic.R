#' Load a HiC_matrix from a .hic file.
#'
#' This HiC_matrix is only a small slice of complete matrix therefore stats like insulation are not valid.
#'
#'
#' @param hic_f a .hic file
#' @param query_gr a query GRanges
#' @param bin_size the bin size to retrieve. Use strawr::readHicBpResolutions to find valid resoltions.
#'
#' @return
#' @export
#' @import strawr
#'
#' @examples
HiC_matrix.from_hic = function(hic_f, query_gr, bin_size = NULL){
    if(is.null(bin_size)){
        bin_size.available = strawr::readHicBpResolutions(hic_f)
        bin_size = max(bin_size.available)
        message("default to largest bin size of: ",
                paste(sort(bin_size.available), collapse = ", "))
    }

    chr = as.character(seqnames(query_gr))
    options(scipen = 999)
    s = start(query_gr)
    e = end(query_gr)


    # undebug(fetch_hic)
    hic_dt = fetch_hic(hic_f, chr, s, e, fill_matrix = TRUE, res = bin_size)
    all_idx = c(
        hic_dt$x/bin_size,
        hic_dt$y/bin_size
    ) %>% sort %>% unique
    dt_1d = data.table(
        index = all_idx,
        seqnames = chr,
        start = all_idx*bin_size,
        end = all_idx*bin_size + bin_size)
    setkey(dt_1d, index)
    dt_2d = data.table(
        i = hic_dt$x / bin_size,
        j = hic_dt$y / bin_size,
        value = hic_dt$counts
    )
    setkey(dt_2d, i, j)

    hmat = new("HiC_matrix")
    hmat@hic_2d = dt_2d
    hmat@hic_1d = dt_1d
    hmat@parameters = HiC_parameters(bin_size = bin_size)
    hmat
}


fetch_hic = function(hic_f,
                     chr,
                     s,
                     e,
                     chr2 = NULL,
                     s2 = NULL,
                     e2 = NULL,
                     fill_matrix = F,
                     res = NULL){
    base_cmd = "VC SUB_HIC_FILE SUB_POS_A SUB_POS_B BP SUB_RES"
    pos_a = paste(c(chr, s, e), collapse = ":")
    if(is.null(chr2)) chr2 = chr
    if(is.null(s2)) s2 = s
    if(is.null(e2)) e2 = e
    pos_b = paste(c(chr2, s2, e2), collapse = ":")
    # cmd = base_cmd
    # cmd = sub("SUB_HIC_FILE", hic_f, cmd)
    # cmd = sub("SUB_POS_A", pos_a, cmd)
    # cmd = sub("SUB_POS_B", pos_b, cmd)
    # cmd = sub("SUB_RES", res, cmd)
    # hic_dt <- as.data.table(straw_R(cmd))
    if(is.null(res)){
        res = strawr::readHicBpResolutions(hic_f)[1]
        message("defaulting resolution to ", res)
    }
    strawr::straw(norm = "VC",
                  fname = hic_f,
                  chr1loc = pos_a,
                  chr2loc = pos_b,
                  unit = "BP",
                  binsize = res,
                  matrix = "oe")

    if(fill_matrix){
        hic_dt = rbind(hic_dt, hic_dt[x != y, .(x = y, y = x, counts)])
    }
    return(hic_dt)
}

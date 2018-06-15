directionality_of_chrRange = function(hic_mat, chr, start, end){
    #determines index positions of specified interval and runs directionality_of_indexes
    #determine overlapping indexes
    m_gr = GRanges(hic_mat@hic_1d)
    start(m_gr) = start(m_gr) + 1
    q_gr = GRanges(seqnames = chr, IRanges(start + 1, end))
    pos_indexes = subjectHits(findOverlaps(query = q_gr, subject = m_gr))
    #pos_indexes have n_bins buffer from start or end of chromosome
    chr_indexes = m_gr[seqnames(m_gr) == chr]$index
    chr_len = length(chr_indexes)
    chr_range = range(chr_indexes)
    # n_bins = round(2*10^6 / hic_mat@parameters@bin_size)
    n_bins = hic_mat@parameters@n_insulation_bins
    to_remove = chr_indexes[c(1:n_bins, (chr_len-n_bins+1):chr_len)]
    if(length(intersect(to_remove, pos_indexes)) > 0){
        pos_indexes[pos_indexes %in% to_remove] = NA
        # warning("some bins in range were too close to chromosome ends and were set to NA")
    }

    directionality_of_indexes(hic_mat, pos_indexes, n_bins)
}

#the directionality range is a square region adjacent to diagonal of the matrix.
#pos_index : the position along diagonal, 1 indexed
#n_bins : the size of the directionality square. smaller n_bins may be sensitive to more local features are dependent on read depth
directionality_range = function(pos_index, n_bins){
    list(pos_index - 1:n_bins, pos_index + 1:n_bins)
}

#

#' DI calculation for single index
#'
#' calculates the directionality index of
#' pos_index position using n_bins in both directions
#'
#' @param hic_matrix
#' @param pos_index
#' @param n_bins
#'
#' @return
#' @export
#'
#' @examples
directionality_of_index = function(hic_matrix, pos_index, n_bins){
    rng = directionality_range(pos_index, n_bins)
    #this is a bit fragile but should hold

    max_miss = 1 - hic_matrix@parameters@min_insulation_coverage
    missed_A = hic_matrix@hic_2d[.(rng[[1]], pos_index)][is.na(val), .N] / n_bins > max_miss
    missed_B = hic_matrix@hic_2d[.(pos_index, rng[[2]])][is.na(val), .N] / n_bins > max_miss
    if(missed_A | missed_B){
        DI = NA
    }else{
        A = hic_matrix@hic_2d[.(rng[[1]], pos_index), sum(val, na.rm = TRUE)]
        B = hic_matrix@hic_2d[.(pos_index, rng[[2]]), sum(val, na.rm = TRUE)]
        E = (A + B) / 2
        DI = ( (B - A) / abs(B - A) )*( (A - E)^2 / E + (B - E)^2 / E )
    }

    return(DI)
}

#runs directionality_of_index on all indexes in pos_indexes
directionality_of_indexes = function(hic_matrix, pos_indexes, n_bins){
    DIs = sapply(pos_indexes, function(pos_index){
        directionality_of_index(hic_matrix, pos_index, n_bins)
    })
    names(DIs) = pos_indexes
    return(DIs)
}

insulation_of_chrRange = function(hic_mat, chr, start, end){
    #determines index positions of specified interval and runs insulation_of_indexes
    #determine overlapping indexes
    m_gr = GRanges(hic_mat@hic_1d)
    start(m_gr) = start(m_gr) + 1
    q_gr = GRanges(seqnames = chr, IRanges(start + 1, end))
    pos_indexes = subjectHits(findOverlaps(query = q_gr, subject = m_gr))
    #pos_indexes have n_bins buffer from start or end of chromosome
    chr_indexes = m_gr[seqnames(m_gr) == chr]$index
    chr_len = length(chr_indexes)
    chr_range = range(chr_indexes)
    n_bins = hic_mat@parameters@n_insulation_bins
    to_remove = chr_indexes[c(1:n_bins, (chr_len-n_bins+1):chr_len)]
    if(length(intersect(to_remove, pos_indexes)) > 0){
        pos_indexes[pos_indexes %in% to_remove] = NA
        # warning("some bins in range were too close to chromosome ends and were set to NA")
    }

    insulation_of_indexes(hic_mat, pos_indexes, n_bins)
}

#the insulation range is a square region adjacent to diagonal of the matrix.
#pos_index : the position along diagonal, 1 indexed
#n_bins : the size of the insulation square. smaller n_bins may be sensitive to more local features are dependent on read depth
insulation_range = function(pos_index, n_bins){
    list(pos_index - 1:n_bins, pos_index + 1:n_bins)
}

#calculates the average value of square sub-matrix extended from diagonal at pos_index position by n_bins in both directions
insulation_of_index = function(hic_matrix, pos_index, n_bins){
    rng = insulation_range(pos_index, n_bins)
    matrix_subset = hic_matrix[rng[[1]], rng[[2]]]
    ins = sum(matrix_subset$val, na.rm = T) / n_bins^2
    cov = sum(!is.na(matrix_subset$val)) / n_bins^2
    min_cov = hic_matrix@parameters@min_insulation_coverage
    if(cov < min_cov){
        ins = NA
    }
    return(ins)
}

#runs insulation_of_index on all indexes in pos_indexes
insulation_of_indexes = function(hic_matrix, pos_indexes, n_bins){
    ins = sapply(pos_indexes, function(pos_index){
        insulation_of_index(hic_matrix, pos_index, n_bins)
    })
    names(ins) = pos_indexes
    return(ins)
}

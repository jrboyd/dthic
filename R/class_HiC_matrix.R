#S4 class to load and analyze HiC hic_2d
# load data and validate
# TAD boundaries by insulation as in:
#   Condensin-driven remodelling of X chromosome topology during dosage compensation
#   http://www.nature.com/nature/journal/v523/n7559/full/nature14450.html
#

setMethod("initialize", "HiC_matrix", function(.Object, matrix_file, regions_file, parameters) {
  # print(matrix_file)
  if(missing(matrix_file) & missing(regions_file) & missing(parameters)){
    return(.Object)
  }
  .Object@matrix_file = matrix_file
  .Object@regions_file = regions_file
  .Object@parameters = parameters

  ###load tall hic_2d data
  mat_dt = fread(matrix_file)
  bed_dt = fread(regions_file)
  #colnames are critical
  colnames(mat_dt) = c("i", "j", "val")
  colnames(bed_dt) = c("seqnames", "start", "end", "index")
  ###conditionally filter out noncanonical chromosomes
  if(.Object@parameters@canonical_chr_only){ #warning("filtering to canonical chr")
    chrms = unique(bed_dt$seqnames)
    k = grepl("^chr[0-9XY][0-9]?$", chrms)
    print(paste("removing", sum(!k), "noncanonical chromosomes (of", length(k), "total)..."))
    chrms[k]
    bed_dt = bed_dt[seqnames %in% chrms[k]]
    mat_dt = mat_dt[i %in% bed_dt$index & j %in% bed_dt$index]
  }
  ###conditionally remove diagonal values
  if(.Object@parameters@diagonal_removed){ #warning("removing diagonal")
    print("removing bins on diagonal...")
    mat_dt = mat_dt[j != i]
  }
  ###conditionally apply a counts per million type adjustment
  if(.Object@parameters@depth_normalization){
    print("normalizing for read depth ( x / total * 10^6 )...")
    total = mat_dt[, sum(val)]
    mat_dt[, val := val / total * 10^6]
  }
  #hic_2d must be keyed by i and j for rapid access
  setkey(mat_dt, i, j)
  setkey(bed_dt, index)
  #set final hic_2d and hic_1d
  .Object@hic_2d = mat_dt
  .Object@hic_1d = bed_dt



  validObject(.Object)
  .Object

})

#' HiC_matrix constructor
#'
#' @param matrix_file HiC-Pro matrix file
#' @param regions_file HiC-Pro region bed file matching matrix file
#' @param parameters HiC_parameters object
#'
#' @return
#' @export
#' @import data.table
#' @import pbapply
#'
#' @examples
HiC_matrix = function(matrix_file, regions_file = NULL, parameters = NULL){
    if(is.null(regions_file)){
        #fix ice first
        mfile = matrix_file
        mfile = sub("_iced.matrix$", ".matrix", mfile)
        mfile = sub("/iced/", "/raw/", mfile)
        regions_file = sub("\\.matrix$", "_abs.bed", mfile)
        if(!file.exists(regions_file)){
            stop("could not auto determine valid regions_file, please supply.")
        }
    }
    if(is.null(parameters)){
        tst = read.table(regions_file, nrows = 2)
        bin_size = tst[2, 2] - tst[1, 2]
        parameters = HiC_parameters(bin_size = bin_size)
    }
    new("HiC_matrix",
        matrix_file = matrix_file,
        regions_file = regions_file,
        parameters = parameters)
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

setMethod("show", "HiC_matrix",
          function(object) {
            nr_mat = nrow(object@hic_2d)
            nr_reg = nrow(object@hic_1d)
            covered = nr_mat / ((nr_reg^2 - nr_reg) / 2)
            print(paste("size is", format(object.size(object), units = "GB")))
            print(paste0(round(covered*100, 2), "% of bins have signal"))
          }
)



# setGeneric("insulation_of_chrRange", function(object, chr, start, end){
#   print(paste("insulation_of_chrRange(object) is not defined for input's class:", class(object)))
# })

# setMethod("insulation_of_chrRange", c("HiC_matrix", "character", "numeric", "numeric"), function(object, chr, start, end){
#hic_mat = object


setMethod("[", c("HiC_matrix", "numeric", "numeric", "ANY"),
          function(x, i, j, ..., drop=TRUE)
          {
            subset_dt(x@hic_2d, i, j)
          })

setMethod("[", c("HiC_matrix", "character", "numeric", "ANY"),
          function(x, i, j, ..., drop=TRUE)
          {
            get_chrRange_Matrix(dt = x@hic_2d, gr = x@hic_1d, chr = i, start = min(j), end = max(j))
          })

setMethod("+", c("HiC_matrix", "HiC_matrix"),
          function(e1, e2)
          {
            if(e1@parameters != e2@parameters){
              stop("ERROR paramters must be identical to add.")
            }
            if(!all(e1@hic_1d[,1:4] == e2@hic_1d[,1:4])){
              stop("ERROR regions must be identical to add.")
            }
            print("valid HiC_matrix sum, proceeding.")
            sum_hic = new("HiC_matrix")
            sum_hic@matrix_file = c(e1@matrix_file, e2@matrix_file)
            sum_hic@regions_file = c(e1@regions_file, e2@regions_file)
            sum_hic@parameters = e1@parameters
            print("merging matrices...")
            m2d = merge(e1@hic_2d, e2@hic_2d, all = T)
            print("adding interaction values...")
            m2d[, c("val.x", "val.y") := .(ifelse(is.na(val.x), 0, val.x), ifelse(is.na(val.y), 0, val.y))]
            m2d = m2d[, .(i, j, val = val.x + val.y)]
            ###conditionally reapply a counts per million type adjustment
            if(sum_hic@parameters@depth_normalization){
              print("normalizing for read depth ( x / total * 10^6 )...")
              total = m2d[, sum(val)]
              m2d[, val := val / total * 10^6]
            }
            sum_hic@hic_2d = m2d
            sum_hic@hic_1d = e1@hic_1d[,1:4]
            print("done!")
            return(sum_hic)
          })

get_chrRange_Matrix = function(dt, gr, chr, start, end)
{
  pos_indexes = get_chrRange_indexes(gr, chr, start, end)
  subset_dt(dt, pos_indexes, pos_indexes)

}

#for input chr start an end, find overlaps with input i_gr
#attempts to allow more natural range specification
#  chrX:100000-300000 with 40k bins should return 5 entries with leftmost starting at 100k and rightmost ending at 300k
#uses i_gr$index if present, otherwise rows in i_gr within range
get_chrRange_indexes = function(i_gr, chr, start, end){
  if(is.data.table(i_gr)) i_gr = GRanges(i_gr)
  start = start + 1 #allows more natural range query 1M-2M should start with left most at 1M and end with rightmost at 2M
  q_gr = GRanges(seqnames = chr, IRanges(start, end))
  start(i_gr) = start(i_gr) + 1

  pos_indexes = subjectHits(findOverlaps(query = q_gr, subject = i_gr, ignore.strand = TRUE))
  if(!is.null(i_gr$index)) pos_indexes = i_gr$index[pos_indexes] #allow set index to override row number if present
  return(pos_indexes)
}

# my_hic_10a = HiC_matrix(matrix_file = matrix_files[1], regions_file = region_files[1], hic_parameters = HiC_parameters(min_insulation_coverage = .3))
# my_hic_10a@hic_1d = calc_delta_for_hic(my_hic_10a@hic_1d, my_hic_10a@parameters@n_delta_bins)
# plot_upperMatrix_with_insulation(hic_mat = my_hic_10a, chr = "chr7", start = 32*10^6, end = 52*10^6, max_fill = .2, max_dist = 3*10^6)

# my_hic = HiC_matrix(matrix_file = matrix_files[3], regions_file = region_files[3], hic_parameters = HiC_parameters(min_insulation_coverage = .3))
# my_hic@hic_1d = calc_delta_for_hic(my_hic@hic_1d, my_hic@parameters@n_delta_bins)
# plot_upperMatrix_with_insulation(hic_mat = my_hic, chr = "chr5", start = 45*10^6 , end = 48*10^6, max_fill = .5)
# plot_upperMatrix_with_insulation(hic_mat = my_hic, chr = "chr6", start = 25*10^6, end = 29*10^6, max_fill = .5)
# plot_upperMatrix_with_insulation(hic_mat = my_hic, chr = "chr7", start = 42*10^6, end = 46*10^6, max_fill = .5)


# MAX = .07
# plot(density(my_hic@hic_2d[val < MAX]$val), type = "n")
# lines(density(my_hic@hic_2d[val < MAX]$val))
# lines(density(my_hic_10a@hic_2d[val < MAX]$val), col = 'red')


setClass(Class = "HiC_parameters",

         slots = c(
             bin_size = "integer",
             diagonal_removed = "logical",
             canonical_chr_only = "logical",
             depth_normalization = "logical",
             quantile_normalization = "logical",
             log2_over_mean_normalization = "logical",
             n_insulation_bins = "integer",
             min_insulation_distance = "numeric",
             min_insulation_coverage = "numeric",
             n_delta_bins = "integer"
         )#,
         #
         # validity = function(object){
         #   errors <- character()
         #   if (!all(colnames(object@matrix) == c("i", "j", "val"))){
         #     msg <- "colnames of matrix must be c(i, j, val)"
         #     errors <- c(errors, msg)
         #   }
         #   if (!all(colnames(object@regions) == c("seqnames", "start", "end", "index"))){
         #     msg <- "colnames of regions must be c(seqnames, start, end, index)"
         #     errors <- c(errors, msg)
         #   }
         #   if (length(errors) == 0) TRUE else errors
         # }
)

setClass(Class = "HiC_matrix",

         slots = c(
             matrix_file = "character",
             regions_file = "character",
             parameters = "HiC_parameters",
             hic_2d = "data.table",
             hic_1d = "data.table"

         ),

         validity = function(object){
             errors <- character()
             mat_cnames = c("i", "j", "val")
             if (length(intersect(colnames(object@hic_2d), mat_cnames)) != length(mat_cnames)){
                 msg <- "colnames of hic_2d must be c(i, j, val)"
                 errors <- c(errors, msg)
             }
             reg_cnames = c("seqnames", "start", "end", "index")
             if (length(intersect(colnames(object@hic_1d), reg_cnames)) != length(reg_cnames)){
                 msg <- "colnames of hic_1d must be c(seqnames, start, end, index)"
                 errors <- c(errors, msg)
             }
             if (length(errors) == 0) TRUE else errors
         }
)

HiC_matrix_wDirectionality = setClass(Class = "HiC_matrix_wDirectionality", contains = "HiC_matrix")

HiC_matrix_wInsulation = setClass(Class = "HiC_matrix_wInsulation", contains = "HiC_matrix")

#class to ensure HiC_matrix objects are handled uniformly for analysis
HiC_sample_set = setClass(Class = "HiC_sample_set",

                          slots = c(
                              HiC_matrix_list = "list",
                              parameters = "HiC_parameters"
                          )
)

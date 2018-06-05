#class to ensure HiC_matrix objects are handled uniformly for analysis
HiC_sample_set = setClass(Class = "HiC_sample_set",

                          slots = c(
                            HiC_matrix_list = "list",
                            parameters = "HiC_parameters"
                          )
)

setMethod("initialize", "HiC_sample_set", function(.Object, matrix_files, regions_files, parameters) {
 if(length(matrix_files) != length(regions_files)){
   stop(paste("length of matrix_files and regions_files must match, were ",
              length(matrix_files),
              length(regions_files)))
 }
  len = length(matrix_files)
  print(paste("loading", len, "samples."))
  .Object@HiC_matrix_list = lapply(1:len, function(i){
    print(paste("sample", i, "of", len))
    HiC_matrix(matrix_files[i], regions_files[i], parameters)
  })
  .Object@parameters = parameters

  ###conditionally apply quantile normalization
  if(.Object@parameters@quantile_normalization){
    print("applying quantile normalization...")
    #apply quantile normalization
    all_mat = matrix(unlist(scores_list), ncol = length(scores_list))
    all_mat_normq = normalize.quantiles(all_mat)
    normq_list = lapply(1:ncol(all_mat_normq), function(i)all_mat_normq[,i])
    names(normq_list) = names(scores_list)
    scores_list = normq_list
  }


  .Object
})



# my_set = HiC_sample_set(matrix_files, region_files, HiC_parameters())

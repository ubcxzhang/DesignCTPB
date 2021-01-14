#' Initiate the reticulate environment
#' @description A utility function to initiate the python environment , check CUDA device and source the python function into the R environment and set it as global function. 
#' @return Error in CudaSupport means there is no cuda device or cuda driver lib cannot be found; otherwise a message of successfully set up the environment will be shown
#' @export 
#' @examples 
#'\dontrun{
#' py_ini()
#' }
py_ini <- function(){
  
  reticulate::py_run_string("import numba; numba.cuda.select_device(0)", convert = TRUE)
  return("Reticulate environment has been set up successfully!")
}

#' Initiate the reticulate environment
#' @description A utility function to initiate the python environment and source the python function into the R environment and set it as global function. 
#' @export 
#' @examples 
#'\dontrun{
#' py_ini()
#' }
py_ini <- function(){
  reticulate::py_config()
  reticulate::source_python(system.file("python","power4R.py",package="DesignCTPB"), envir = .GlobalEnv, convert = TRUE) # source python4R.py into the environment
  assign(Power.sampling,  Power_sampling, envir = .GlobalEnv, inherits=TRUE) # set Power.sampling as an global function
}

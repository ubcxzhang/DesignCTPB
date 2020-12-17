#' Setup reticulate environment
#' @description A utility function to create and setup reticulate virtual environment and source the python function into the R environment. 
#'
#' @param reticulate_venv character, the name of the reticulate virtual environment to setup.
#' @param py_path character, python path, the default value is NULL which the default python path will be applied
#' @details set up the python environment, and check whether the required python modules installed or not, if not install it. 
#' @export 
#' @examples 
#' create_venv()
#' 
create_venv <- function(reticulate_venv, py_path=NULL) {
  if(is.element(reticulate_venv, reticulate::virtualenv_list())){
    print(paste("Reticulate environment",reticulate_venv, "is already exists!", sep=' ' ))
  }
  else{
    print(paste("Creating reticulate environment with the name:", reticulate_venv))
    if(is.null(py_path)){
      reticulate::virtualenv_create(reticulate_venv)
    }
    else{
      reticulate::virtualenv_create(reticulate_venv, python = py_path)
    }
  }
  reticulate::virtualenv_install(reticulate_venv,c('scipy','pandas','numba'))
}

#' Initiate the reticulate environment
#' @description A utility function to initiate the python environment created by create_venv(), source the python function into the R environment and set it as global function. 
#' @param reticulate_venv character, the name of the reticulate virtual environment has been setup
#' @export 
#' @examples 
#' py_initial()
#' 
py_initial <- function(reticulate_venv){
  reticulate::use_virtualenv(reticulate_venv, required = TRUE)
  reticulate::py_config()
  reticulate::source_python(system.file("python","power4R.py",package="DesignCTPB")) # source python4R.py into the environment
  Power.sampling <<- Power_sampling # set Power.sampling as an global function
}


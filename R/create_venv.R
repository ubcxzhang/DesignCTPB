#' A utility function to create and setup reticulate virtual environment and source the python function into the R environment. 
#'
#' @export 
#' @examples 
#' create_venv()
#' 
create_venv <- function(reticulate_venv="designctpb_numba",py_path=NULL) {
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

#' A utility function to source the python function into the R environment. 
#'
#' @export 
#' @examples 
#' py_initial()
#' 
py_initial <- function(reticulate_venv="designctpb_numba"){
  reticulate::use_virtualenv(reticulate_venv, required = TRUE)
  reticulate::py_config()
  reticulate::source_python(system.file("python","power4R.py",package="DesignCTPB")) # source python4R.py into the environment
  Power.sampling <<- Power_sampling # set Power.sampling as an global function
}


## This function ensures that all directories we save output to exist
#' Fix Directories
#' 
#' Ensures a set of subdirectories of the current working directory exist;
#' if they do not exist, they are created.
#' 
#' @param directories A character vector giving the subdirectories needed
fix_directories <- function(directories = c("plots", "model-output")) {
    for ( directory in directories ) {
        if ( !dir.exists(directory) ) {
            dir.create(directory)
        }
    }
    return(invisible(NULL))
}

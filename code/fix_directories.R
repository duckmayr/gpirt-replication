## This function ensures that all directories we save output to exist
fix_directories <- function() {
    if ( !dir.exists("plots") ) {
        dir.create("plots")
    }
    if ( !dir.exists("model-output") ) {
        dir.create("model-output")
    }
}

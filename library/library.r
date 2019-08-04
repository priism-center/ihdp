############################################
# Library of functions for propensity analysis work

############################################


#' Compare model formula.
#' 
#' @param spec1 formula specification 1 or model 1
#' @param spec2 formula specification 2 or model 2
#' @param rev BOOLEAN, whether to use spec 2 as the base of comparison
formDiff <- function(spec1, spec2, rev=FALSE){
    if (rev) {
        temp <- spec1
        spec1 <- spec2
        spec2 <- temp
    }
    if (!identical(class(spec1), 'formula')) {
        try(spec1 <- spec1$formula)
    }
    if (!identical(class(spec2), 'formula')) {
        try(spec2 <- spec2$formula)
    }

    spec1 <- sort(all.vars(spec1))
    spec2 <- sort(all.vars(spec2))

    if (identical(spec1, spec2)) {
        print('Specifications are identical')
    } else {
        setdiff(sort(all.vars(spec1)), sort(all.vars(spec2)))
    }
}

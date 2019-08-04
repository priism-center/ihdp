############################################
# Library of functions for propensity analysis work

############################################


#' Compare model formula.
#' 
#' @param spec1 formula specification 1 or model 1
#' @param spec2 formula specification 2 or model 2
formDiff <- function(spec1, spec2){
    if (!identical(class(spec1), 'formula')) {
        spec1 <- spec1$formula
    }
    if (!identical(class(spec2), 'formula')) {
        spec2 <- spec2$formula
    }

    spec1 <- sort(all.vars(spec1))
    spec2 <- sort(all.vars(spec2))

    if (identical(spec1, spec2)) {
        print('Specifications are identical')
    } else {
        print(cat('Spec 1 adds: ', setdiff(spec1, spec2)))
        print(cat('Spec 1 misses: ', setdiff(spec2, spec1)))
    }
}

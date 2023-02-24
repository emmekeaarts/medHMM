#' @keywords internal
# Calculates the probabilities of observing each state at each point in time given
# the observations of all dependent variables, used for the forward probabilities
# Based on Zuchini 2016.
get_all1 <- function(x, emiss, n_dep, data_distr){
    inp <- rep(list(NULL), n_dep)
    if(data_distr == "categorical"){
        for(q in 1:n_dep){
            inp[[q]] <- t(emiss[[q]][,x[,q]])
        }
    } else if (data_distr == "continuous"){
        for(q in 1:n_dep){
            inp[[q]] <- outer(x[,q], Y = emiss[[q]][,1], FUN = dnorm, sd = rep(sqrt(emiss[[q]][,2]), each = dim(x)[1]))
        }
    } else if (data_distr == "poisson") {
        for(q in 1:n_dep){
            inp[[q]] <- outer(x[,q], emiss[[q]][,1], FUN = dpois)
        }
    }
    allprobs <- Reduce("*", inp)
    return(allprobs)
}

#' @keywords internal
# Calculate initial probabilities assuming stationary distribution
get_delta <- function(gamma, m){
    solve(t(diag(m) - gamma + 1), rep(1, m))
}

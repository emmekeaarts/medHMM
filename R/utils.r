#----------------------------------------------------------------------#
#                       General utility functions                      #
#----------------------------------------------------------------------#

#' @keywords internal
#' get d distribution for log normal
get.d.lognorm <- function(run.p, Mx, m){
    input.d <-  matrix(rep(0:(Mx-1), m), nrow = m, byrow = TRUE)
    d <- cbind(rep(0,m), apply(input.d, 2, dlnorm, meanlog = run.p$logmu, sdlog = run.p$logsd))
    return(d)
}

#' @keywords internal
#' computes probabilities from intercepts
int.to.prob <- function(int1) {
    exp.int1 	<- matrix(exp(c(0, int1)), nrow  = 1)
    prob1 		<- exp.int1 / as.vector(exp.int1 %*% c(rep(1,(dim(exp.int1)[2]))))
    return(prob1)
}


#' @keywords internal
#' computes probabilities from intercepts
int_to_prob <- function(int1) {
    if(is.matrix(int1)){
        prob1 <- matrix(nrow = nrow(int1), ncol = ncol(int1) + 1)
        for(r in 1:nrow(int1)){
            exp_int1 	<- matrix(exp(c(0, int1[r,])), nrow  = 1)
            prob1[r,] <- exp_int1 / as.vector(exp_int1 %*% c(rep(1, (dim(exp_int1)[2]))))
        }
    } else {
        exp_int1 	<- matrix(exp(c(0, int1)), nrow  = 1)
        prob1 		<- exp_int1 / as.vector(exp_int1 %*% c(rep(1, (dim(exp_int1)[2]))))
    }
    return(round(prob1,4))
}

#' @keywords internal
#' computes intercepts from probabilities, per row of input matrix
#' first catagory is reference catagory
prob_to_int <- function(prob1){
    prob1 <- prob1 + 0.00001
    b0 <- matrix(NA, nrow(prob1), ncol(prob1)-1)
    sum_exp <- numeric(nrow(prob1))
    for(r in 1:nrow(prob1)){
        sum_exp[r] <- (1/prob1[r,1]) - 1
        for(cr in 2:ncol(prob1)){
            #for every b0 except the first collumn (e.g. b012 <- log(y12/y11-y12))
            b0[r,(cr-1)] <- log(prob1[r,cr]*(1+sum_exp[r]))
        }
    }
    return(round(b0,4))
}

#' @keywords internal
#' Old version: evaluate log-like for intercept only MNL, bassed on p.rossi 2004
llmnl_int_old <- function(int, Obs, n_cat) {
    n_Obs 		<- length(Obs)
    Xint 		<- matrix(c(0, int), byrow = T, ncol = n_cat, nrow = n_Obs)
    ind			<- cbind(c(1:n_Obs),Obs)
    Xby			<- Xint[ind]
    Xint		<- exp(Xint)
    iota		<- c(rep(1,(n_cat)))
    denom		<- log(Xint%*%iota)
    return(sum(Xby-denom))
}

#' @keywords internal
#' Evaluate loglikelihood for intercept only MNL, bassed on P. Rossi 2004
llmnl_int <- function(int, Obs) {
    n_Obs <- length(Obs)
    betas <- c(0, int)
    return(sum(betas[Obs]) - log(sum(exp(betas))) * n_Obs)
}

#' @keywords internal
#' fractional log likelihood for multinomial intercept only model
llmnl_int_frac <- function(int, Obs, pooled_likel, w, wgt){
    return((1-w)*llmnl_int(int = int, Obs = Obs) + w* wgt* pooled_likel)
}

#' @keywords internal
#' Old version: compute mnl -Expected[Hessian]  for intercept only model, bassed on p.rossi 2004
mnlHess_int_old <- function(int, Obs, n_cat){
    k 			<- n_cat - 1
    n_Obs 		<- length(Obs)
    Xint 		<- matrix(c(0, int), byrow = T, ncol = n_cat, nrow = n_Obs)
    Xint		<- exp(Xint)
    iota		<- c(rep(1,n_cat))
    denom		<- Xint%*%iota
    Prob		<- Xint/as.vector(denom)

    Hess=matrix(double(k*k),ncol=k)
    for (hi in 1:n_Obs) {
        p <- as.vector(Prob[hi,])
        A <- diag(p)-outer(p,p)
        Xt <- diag(n_cat)[,-1]
        Hess <- Hess+crossprod(Xt,A)%*%Xt
    }
    return(Hess)
}

#' @keywords internal
#' faster version
mnlHess_int <- function(int, Obs, n_cat){
    n_Obs 	<- length(Obs)
    betas   <- matrix(c(0, int), byrow = TRUE, ncol = n_cat)
    prob    <- exp(betas) / sum(exp(betas))
    Hess    <- (diag(x = prob[-1], nrow = n_cat-1) - outer(prob[-1],prob[-1])) * n_Obs
    return(Hess)
}


#----------------------------------------------------------------#
#                       Metropolis-Hastings                      #
#----------------------------------------------------------------#

#' @keywords internal
#' one run of the random walk metropolis sampler for an intercept only multinomial distribution
mnl_IndepM_once <- function(int1, Obs, n_cat, mu_int_bar1, V_int1, scalar, candcov1, int.star1, nu.M1) {
    # obtain likelihood and transition prob with the parameters sampled in the previous iteration and current sampled state sequence
    oldloglike	 		<- llmnl_int(int = int1, Obs = Obs, n_cat = n_cat)
    oldpostlike	 		<- oldloglike + dmvnorm(int1, mu_int_bar1, V_int1, log = TRUE)
    constant_t 			<- (nu.M1/2)*log(nu.M1)+lgamma((nu.M1+(n_cat-1))/2)-((n_cat-1)/2)*log(pi)-lgamma(nu.M1/2)
    oldlimp				<- dmvt(int1, delta = int.star1, sigma = candcov1, df = nu.M1) - constant_t
    probold				<- int.to.prob(int1)

    # obtain new parameters for gamma from proposal distribution plus new likelihood
    int.new		 		<- rmvt(1, sigma = candcov1, delta = int.star1, df = nu.M1, method = "svd")
    newloglike	 		<- llmnl_int(int = int.new, Obs = Obs, n_cat = n_cat)
    newpostlike	 		<- newloglike + dmvnorm(int.new, mu_int_bar1, V_int1, log = TRUE)
    newlimp				<- dmvt(int.new, delta = int.star1, sigma = candcov1, df = nu.M1) - constant_t
    probnew				<- int.to.prob(int.new)

    # determine to use the updated or current (previous iteration) gamma values of the parameters
    acc 				<- min(log(1), (newpostlike + oldlimp - oldpostlike - newlimp))
    if(acc < log(1))
    {unif = log(runif(1))} else {unif = log(1)}
    if (unif <= acc) {
        draw.int		 	<- int.new
        accept			<- 1
        prob			<- probnew
    } else {
        draw.int			<- int1
        accept			<- 0
        prob			<- probold
    }
    return(list(draw.int = draw.int, accept = accept, prob = prob))
}


#' @keywords internal
#' one run of the random walk metropolis sampler for an intercept only multinomial distribution
mnl_RW_once <- function(int1, Obs, n_cat, mu_int_bar1, V_int1, scalar, candcov1) {
    # obtain likelihood and transition prob with the parameters sampled in the previous iteration and current sampled state sequence
    oldloglike	 		<- llmnl_int(int = int1,
                               # n_cat = n_cat,
                               Obs = Obs)
    oldpostlike	 		<- oldloglike + dmvnorm(int1, mu_int_bar1, V_int1, log = TRUE)
    probold				<- int.to.prob(int1)

    # obtain new parameters for gamma from proposal distribution plus new likelihood
    int.new		 		<- int1 + rmvnorm(1, rep(0,(n_cat-1)), scalar^2*candcov1, method = "svd")
    newloglike	 		<- llmnl_int(int = int.new,
                               # n_cat = n_cat,
                               Obs = Obs)
    newpostlike	 		<- newloglike + dmvnorm(int.new, mu_int_bar1, V_int1, log = TRUE)
    probnew				<- int.to.prob(int.new)

    # determine to use the updated or current (previous iteration) gamma values of the parameters
    acc 				<- min(log(1), (newpostlike - oldpostlike))
    if(acc < log(1))
    {unif = log(runif(1))} else {unif = log(1)}
    if (unif <= acc) {
        draw.int		<- int.new
        accept			<- 1
        prob			<- probnew
    } else {
        draw.int		<- int1
        accept			<- 0
        prob			<- probold
    }
    return(list(draw_int = draw.int, accept = accept, prob = prob))
}



#-----------------------------------------------------------------------#
#                       Utility function mHMMbayes                      #
#-----------------------------------------------------------------------#

#' @keywords internal
#' Whenever you use C++ code in your package, you need to clean up after yourself
#' when your package is unloaded. This function unloads the DLL (H. Wickham(2019). R packages)
.onUnload <- function (libpath) {
    library.dynam.unload("mHMMbayes", libpath)
}

#' @keywords internal
#' simple functions used in mHMM
dif_matrix <- function(rows, cols, data = NA){
    return(matrix(data, ncol = cols, nrow = rows))
}

#' @keywords internal
nested_list <- function(n_dep, m){
    return(rep(list(vector("list", n_dep)),m))
}

#' @keywords internal
dif_vector <- function(x){
    return(numeric(x))
}

#' @keywords internal
is.whole <- function(x) {
    return(is.numeric(x) && floor(x) == x)
}

#' @keywords internal
is.mHMM <- function(x) {
    inherits(x, "mHMM")
}

#' @keywords internal
is.mHMM_cont <- function(x) {
    inherits(x, "mHMM_cont")
}

#' @keywords internal
is.mHMM_gamma <- function(x) {
    inherits(x, "mHMM_gamma")
}

#' @keywords internal
hms <- function(t){
    paste(formatC(t %/% (60*60) %% 24, width = 2, format = "d", flag = "0"),
          formatC(t %/% 60 %% 60, width = 2, format = "d", flag = "0"),
          formatC(t %% 60, width = 2, format = "d", flag = "0"),
          sep = ":")
}

#----------------------------------------------------------------------#
#                       General utility functions                      #
#----------------------------------------------------------------------#

# get d distribution for log normal 
get.d.lognorm <- function(run.p, Mx, m){
    input.d <-  matrix(rep(0:(Mx-1), m), nrow = m, byrow = TRUE)
    d <- cbind(rep(0,m), apply(input.d, 2, dlnorm, meanlog = run.p$logmu, sdlog = run.p$logsd))
    return(d)
}


# computes probabilities from intercepts
int.to.prob <- function(int1) {
    exp.int1 	<- matrix(exp(c(0, int1)), nrow  = 1)
    prob1 		<- exp.int1 / as.vector(exp.int1 %*% c(rep(1,(dim(exp.int1)[2]))))
    return(prob1)
}


# fractional log likelihood for multinomial intercept only model
llmnl_int_frac <- function(int, Obs, n_cat, pooled.likel, w, wgt){
    return((1-w)*llmnl_int(int = int, Obs = Obs, n_cat = n_cat) + w* wgt* pooled.likel)		
}


# evaluate log-like for intercept only MNL, bassed on p.rossi 2004
llmnl_int <- function(int, Obs, n_cat) {
    n_Obs 		<- length(Obs)
    Xint 		<- matrix(c(0, int), byrow = T, ncol = n_cat, nrow = n_Obs)
    ind			<- cbind(c(1:n_Obs),Obs)
    Xby			<- Xint[ind]
    Xint		<- exp(Xint)
    iota		<- c(rep(1,(n_cat)))
    denom		<- log(Xint%*%iota)
    return(sum(Xby-denom))
}


# compute mnl -Expected[Hessian]  for intercept only model, bassed on p.rossi 2004
mnlHess_int <- function(int, Obs, n_cat){
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



#----------------------------------------------------------------#
#                       Metropolis-Hastings                      #
#----------------------------------------------------------------#

# one run of the random walk metropolis sampler for an intercept only multinomial distribution
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


# one run of the random walk metropolis sampler for an intercept only multinomial distribution
mnl_RW_once <- function(int1, Obs, n_cat, mu_int_bar1, V_int1, scalar, candcov1) {
    # obtain likelihood and transition prob with the parameters sampled in the previous iteration and current sampled state sequence
    oldloglike	 		<- llmnl_int(int = int1, Obs = Obs, n_cat = n_cat)
    oldpostlike	 		<- oldloglike + dmvnorm(int1, mu_int_bar1, V_int1, log = TRUE)
    probold				<- int.to.prob(int1)
    
    # obtain new parameters for gamma from proposal distribution plus new likelihood
    int.new		 		<- int1 + rmvnorm(1, rep(0,(n_cat-1)), scalar^2*candcov1, method = "svd")
    newloglike	 		<- llmnl_int(int = int.new, Obs = Obs, n_cat = n_cat)
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
    return(list(draw.int = draw.int, accept = accept, prob = prob))
}

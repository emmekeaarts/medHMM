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

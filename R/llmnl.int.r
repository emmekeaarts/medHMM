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

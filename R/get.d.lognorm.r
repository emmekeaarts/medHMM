# get d distribution for log normal 
get.d.lognorm <- function(run.p, Mx, m){
	input.d <-  matrix(rep(0:(Mx-1), m), nrow = m, byrow = TRUE)
	d <- cbind(rep(0,m), apply(input.d, 2, dlnorm, meanlog = run.p$logmu, sdlog = run.p$logsd))
	return(d)
}

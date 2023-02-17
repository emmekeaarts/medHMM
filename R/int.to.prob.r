# computes probabilities from intercepts
int.to.prob <- function(int1) {
	exp.int1 	<- matrix(exp(c(0, int1)), nrow  = 1)
	prob1 		<- exp.int1 / as.vector(exp.int1 %*% c(rep(1,(dim(exp.int1)[2]))))
	return(prob1)
}

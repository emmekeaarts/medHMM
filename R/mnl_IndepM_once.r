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

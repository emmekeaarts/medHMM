

# Data: list contaning data properties and data; y = observed data, m = number of states, q_pr = number of observed behavioral acts 
# (= cateogries), maximum runlength Mx and maximum runlength per state Mx2
# Start_Val: vector containing start values for the parameters in the model, in the order: gamma, pr, logmu, logsigma. 
# Gamma_Sampler: list containing the start values for the mle algorithm, int.mle0, mu for the proposition distribution, mu_prop 
# and the scalar of the covariance matrix, scalar. 
# Pr_Sampler: list containing the start values for the mle algorithm, pr.int.mle0, mu for the proposition distribution, mu_pr.prop 
# and the scalar of the covariance matrix, pr.scalar.
# Gamma_Prior: list containing values for prior parameters of gamma, mu_int0 and root0
# Pr_Prior: list containing values for prior parameters of pr, mu_pr.int0 and pr.root0
# D_Prior: list containing values for prior parameters of state runlength distribution, logmu0, tau20, alpha.sigma20 and beta.sigma20
# Mcmc: list containing number of iterations J, and number of iterations until burn in, burn.in


HSMM_mnl_Gibbs <- function(Data, Start_Val, Gamma_Prior = NULL, Pr_Prior = NULL, D_Prior, Mcmc, Fast_vers = FALSE, Fast_vers_par = NULL){ 

# starting_values, tau20, logmu0, alpha.sigma20, beta.sigma20, J = 100, inf.prior = FALSE, burn.in = 50, Mx, Mx2)

	# Initialize data	
		y 				<- Data$y
		y2 				<- y - 1
		n 				<- length(y)
		switch 			<- rep(0, n) 
		switch[1] 		<- 1
		for (i in 2:n) {if(y[i] != y[i-1]) {switch[i] <- 1}}
		switch2 		<- c(switch[-1],1)
		m 				<- Data$m
		q_pr 			<- Data$q_pr
		Mx 				<- Data$Mx
		Mx2 			<- Data$Mx2

		if (Fast_vers == TRUE){
			if(is.null(Fast_vers_par)){
				act.sleep <- 6
				lag.sleep <- 60
			}	else {
				act.sleep	<- Fast_vers_par$act.sleep
				lag.sleep 	<- Fast_vers_par$lag.sleep 
			}
			if(sum(y == act.sleep) > 0){		
				sleep.ind <- which(switch == 1 & y == act.sleep)
				sleep.ind2 <- numeric()
				for(k in 1:length(sleep.ind)){
					sleep.ind2 <- c(sleep.ind2, ((sleep.ind[k])-lag.sleep):sleep.ind[k])
				}
				sleep.ind2 <- sleep.ind2[sleep.ind2 > 0]
				sleep.ind3 <- rep(0,n)
				sleep.ind3[sleep.ind2] <- 1
				switch.sl <- as.numeric(switch == 1 & sleep.ind3 == 1)
			} else {
				switch.sl <- sleep.ind3 <- rep(0,n)
			} 
		}



	# Initialize Gamma prior
	if(is.null(Gamma_Prior)){
		alpha.gamma0	<- (diag(m) -1) * -1 # matrix with diagonal elements zero and off diagonal elements 1. 
	} else {
		alpha.gamma0	<- Gamma_Prior$alpha.gamma0
	}

	# Initialize Pr prior
	if(is.null(Pr_Prior)){
		alpha.pr0	<- matrix(rep(1,(q_pr*m)), ncol = q_pr, nrow = m)
		inf.prior	<- FALSE
	} else {
		alpha.pr0 	<- Pr_Prior$alpha.pr0
		inf.prior 	<- Pr_Prior$inf.prior
	}
	
	# Initialize D prior
	logmu0			<- D_Prior$logmu0
	tau20			<- D_Prior$tau20
	alpha.sigma20	<- D_Prior$alpha.sigma20
	beta.sigma20	<- D_Prior$beta.sigma20
	
	# Initialize Mcmc argumetns
	J 				<- Mcmc$J
	burn.in			<- Mcmc$burn.in	
	
	# Define arguments used in Mcmc algorithm
	gamma.counts <- matrix(0,m,m)
	pr.counts <- matrix(0,m,q_pr)
	tau2 <- logmu2 <- numeric(m)
	a1 <- b1 <- numeric(m)
	sampled.states <- matrix(,nrow = n, ncol = J)
	
	# Create outcome matrix of the parameters
	PD 				<- matrix(, nrow = J, ncol = m*q_pr + m*m + m*2 + 1)
	colnames(PD)	<- c(paste("pr",rep(1:q_pr, m),"_S",rep(1:m, each = q_pr), sep = ""), paste("S", rep(1:m, each = m), "to", rep					(1:m, m), sep = ""), paste("logmu", 1:m, sep = ""), paste("logvar", 1:m, sep = ""), "LL")
	PD[1, 1:((m*q_pr + m*m) + m*2)] <- Start_Val
	#write(PDcolnames, file = namefileR, append = F, ncol = m*q_pr + m*m + m*2 + 1, sep = ";")
	
	# Define model parameters by starting values	
	pr 				<- matrix(Start_Val[1:(m*q_pr)], ncol = q_pr, nrow = m, byrow = TRUE)
	gamma 			<- matrix(Start_Val[(m*q_pr + 1) : (m*q_pr + m*m)], byrow = TRUE, ncol = m)
	delta 			<- solve(t(diag(m) - gamma + 1), rep(1, m))
	logmu			<- Start_Val[((m*q_pr + m*m) + 1) : ((m*q_pr + m*m) + m)]
	logsigma2		<- Start_Val[((m*q_pr + m*m) + m + 1) : ((m*q_pr + m*m) + m*2)]
	
	alpha.pr0.sc 	<- alpha.pr0
	pr0.sc.counter  <- matrix(, burn.in, m) 

	# RUN THE GIBBS SAMPLER
	itime=proc.time()[3]
	for (iter in 2 : J){
		# block sample the state sequence (xt) from its conditional distribuition by HSMM message-passing scheme. 

		# first, obtain duration densities and employ the forward backward algorithm
		# outcome FBalgC : 1 = N, 2 = Norm, 3 = Forward, 4 = Backward, 5 = B_star 6 = Occupancy
		d 				<- get.d.lognorm(run.p = list(logmu = logmu, logsd = sqrt(logsigma2)), Mx = Mx, m = m)
		if (Fast_vers == FALSE){
			FB				<- FBalgC(y2 = y2, m = m, pr = pr, Mx = Mx, Mx2 = Mx2, gamma = gamma, d = d, S2 = switch2, S = switch, delta = delta) 
		} else {
			FB				<- FBalgCF(y2 = y2, m = m, pr = pr, Mx = Mx, Mx2 = Mx2, gamma = gamma, d = d, S2 = switch2, S = switch, delta = delta, S3 = switch.sl) 
		}		
		sample.path 	<- numeric()
		Dur 			<- numeric()
		
		B_star 			<- FB[[2]]
		Occupancy 		<- FB[[3]]
		N 				<- FB[[1]]
		
	# first state
		sample.path[1] <- sample(1:m, 1, prob = delta * exp(B_star[,1])) 
		Dur[1] <- sample(1:Mx, 1, prob = (Occupancy[[1]][sample.path[1],] / exp(B_star[sample.path[1],1])))		
		sampled.states[1:Dur[1], iter] <- sample.path[1]
	# following states 	
		gamma.counts[,] <- 0
		s <- 1
		while(sum(Dur) < n){
			s <- s + 1
			Mx.l <- min(Mx, n-sum(Dur))
			sample.path[s] <- sample(1:m, 1, prob = gamma[sample.path[s-1],] * exp(B_star[,sum(Dur)+1])) 
			gamma.counts[sample.path[s-1],sample.path[s]] <- gamma.counts[sample.path[s-1],sample.path[s]] + 1
			Dur[s] <- sample(1:Mx.l, 1, prob = (Occupancy[[sum(Dur)+1]][sample.path[s],] / exp(B_star[sample.path[s], sum(Dur)+1])))		
			sampled.states[sum(Dur[1:s-1],1):sum(Dur), iter] <- sample.path[s]
		}
		n.Dur <- length(Dur)
	
	# scale prior for conditional probabilities to counts if prior is informative
			if(inf.prior == TRUE & iter <= burn.in){
				for (i in 1:m){
					pr0.sc.counter[iter, i] <- sum(sampled.states[, iter] == i) 
					alpha.pr0.sc[i,] <- c(alpha.pr0[i,] * pr0.sc.counter[iter, i] / (n/10)) + 1
				}
			}
			if(inf.prior == TRUE & iter == (burn.in + 1)){
				for (i in 1:m){
					alpha.pr0.sc[i,] <- c(alpha.pr0[i,] * mean(pr0.sc.counter[(burn.in-burn.in/2):burn.in, i]) / (n/10)) + 1
				}
			}
	
	# with the sample path, updata gamma, pr and duration dist
		# update gamma
			for (i in 1:m){
				gamma[i,] <- rdirichlet(1, alpha = (gamma.counts[i,] + alpha.gamma0[i,]))
			}
			PD[iter,(m*q_pr+1):(m*q_pr+m*m)] <- as.vector(t(gamma))
		
		# update pr	
			pr.counts[,] <- 0
			for (i in 2:n){
				for (j in 1:m){
					for (k in 1:q_pr){
						if (sampled.states[i, iter] == j & y[i] == k) {
							pr.counts[j,k] <- pr.counts[j,k] + 1
						}
					}
				}
			}	
			for (i in 1:m){
				pr[i,] <- rdirichlet(1, alpha = (pr.counts[i,] + alpha.pr0.sc[i,]))
			}	
			PD[iter,1:(m*q_pr)] <- as.vector(t(pr))

		# update duration dist
			for (j in 1:m){
				tau2[j] <- 1/ ((1/tau20[i]) + (1/logsigma2[j]) * sum(sample.path[-n.Dur] == j))
				logmu2[j] <- tau2[j] * ((1/tau20[i]) * logmu0[j] + (1/logsigma2[j]) * sum(log(Dur[-n.Dur][sample.path[-n.Dur] == j])))
				logmu[j] <- rnorm(1, mean = logmu2[j], sd = sqrt(tau2[j]))	
			}
			PD[iter,(m*q_pr+m*m+1):(m*q_pr+m*m+m)] <- logmu
			
			for (j in 1:m){
				a1[j] <- alpha.sigma20[i] + sum(sample.path[-n.Dur] == j) / 2
				b1[j] <- beta.sigma20[i] + (1/2) * (sum((log(Dur[-n.Dur][sample.path[-n.Dur] == j]) - logmu[j])^2)) 
				logsigma2[j] <- 1 / rgamma(1, shape = a1[j], rate = b1[j])	
			}
			PD[iter,(m*q_pr+m*m+m+1):(m*q_pr+m*m+m*2)] <- logsigma2
	
	# Save log likelihood of iteration
		print(c(iter, sum(N)))
	}
	ctime=proc.time()[3]
	print(paste("total time elapsed", round((ctime-itime)/60,2)))
	return(list(PD = PD, sampled.states = sampled.states))
}

#######################################
#######################################



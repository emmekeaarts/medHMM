
#############################################################################################
# 				GIBBS SAMPLER WITH METROPOLIS STEP, RW METROPOLIS SAMPLER					#
#############################################################################################



# Model parameters are sampled from distributions based on sampled state sequence. 
# The state sequence is sampled using the outcomes of the forward backward algorithm: first run forward sequence, then backward sequence, 
# then use backward probabilities to do forward sampling of the states. The forward recursion is needed since the likelihood in the backward 
# recursion is scaled using the (pseudo) likelihood obtained in the forward recursion. When not scaling, the probabilities of the likelihood
# become infinitely small after about 500 observations, also when put on log scale. 

# gamma	= matrix that gives the transition probability from state i to j, where the diagonal (so i = j) is fixed to zero
# pr = matrix that gives the conditional distributions of the observations given the state

# Model:
# pi_ (rows i of gamma) 	~(idd) 	Multinomial logit(int) for every row of gamma
									# i.e., categorical distribution modeled by an multinomial intercept only logit model.  
									# i = j is fixed at zero. The first intercept is fixed at zero for identification. 
									# Therefore, we estimate J - 2 intercepts for each row of gamma.  
# int						~		N(mu_int0, root0^-1), i.e., the intercepts have a normal prior	with precision root0				
# Zs (superstate s) 		~ 		probability vector equal to the row of gamma(Zs-1)
# Yt1s:t2s(theta_i) (obs) 	~(idd) 	H(theta_Zs)
									# Conditional distribution for observations given the (super)state
# H(theta_i)				~		Multinomial logit(pr_int) for every row of pr	
									# i.e., categorical distribution modeled by an multinomial intercept only logit model.  
									# the first intercept is always fixed at zero, so q_pr - 1 intercept are estimated
# pr_int					_		N(mu_pr.int0, pr_root0), i.e., the intercepts have a normal prior with precision pr_root0								
# Ds 						~ 		G(w_Zs), 
									# superstate specific runlength distribution G
# G(w_Zsq_pr)					~		LogN(logmu, logsigma)
# logmu						~		N(logmu0, tau20)
# logsima					~		invGamma(alpha.simga20, beta.sigma20)
# Xt1s:t2s 					= 		Zs 
									# (label sequence of the states, so not the superstates, but just each state for each observation) 
# t1s 						= 		sum over all previous durations 
# t2s 						= 		t1s + Ds (current duration) - 1
# (theta_i, w_i) 			~ (idd)	H*G = the distribution of the duration of the states and the conditional distribution of the the observations given the states are independent of each other. 


# For model, see p. 8 of Johnson and Willsky 2012, Bayesian non-parametric hidden semi Markov models
# label sequence xt is defined for convenience, the par Zs, Ds is the runlenght encoding of xt. 

# The outcomes of the Forward backward algorithm are:
# -  B_star: the probability that the state sequence from t1s to t2s is Y given the observations, i.e., p(Y_t+1:d | X_t+1). 
#	 B_star is composed of lists, where the lists corresponds to the probabilities starting at time point t1s + 1, 
#	 so B_star[[1]] corresponds to t1s = 0. B_star denotes the most likely duration of the state that started at t1s.
# - Occupancy: The sum over each list of B_star, denotes the overall most likely state starting at t1s.
# - N: the (pseudo) likelihood for each observation, it is pseudo as it only is obtained when a new act starts instead of every point in time.

# To reduce the computational cost and because it is biological logical that states only start when the animal has switched to another act, 
# switches from one state to another can only occur when the observed act changes. These switches from one act to another is given by 
# the parameter switch (dummy indicator indicating the start of a new act) and switch2 (dummy indicator indicating the end of an act)

# int and pr_int are sampled using a random walk metropolis sampler. 
# This requires a proposal distribution, which is set to be a normal distribution with mu = mu_prop / mu_pr.prop and scalar^2 / pr.scalar^2 
# * covariance matrix. The scalar is by default set to 2.93 / sqrt(number of parameters to be estimated).
# The covariance matrix is obtained from the maximum likelihood estimate (mle) of the multinomial logit model using the sampled states, 
# The start values for the mle algorithm are give by int.mle0 and pr.int.mle0.  

# Data: list containing data properties and data; y = observed data, m = number of states, q_pr = number of observed behavioral acts 
# (= categories), maximum runlength Mx and maximum runlength per state Mx2
# Start_Val: vector containing start values for the parameters in the model, in the order: gamma, pr, logmu, logsigma. 
# Gamma_Sampler: list containing the start values for the mle algorithm, int.mle0, mu for the proposition distribution, mu_prop 
# and the scalar of the covariance matrix, scalar. 
# Pr_Sampler: list containing the start values for the mle algorithm, pr.int.mle0, mu for the proposition distribution, mu_pr.prop 
# and the scalar of the covariance matrix, pr.scalar.
# Gamma_Prior: list containing values for prior parameters of gamma, mu_int0 and root0
# Pr_Prior: list containing values for prior parameters of pr, mu_pr.int0 and pr.root0
# D_Prior: list containing values for prior parameters of state runlength distribution, logmu0, tau20, alpha.sigma20 and beta.sigma20
# Mcmc: list containing number of iterations J, and number of iterations until burn in, burn.in


HSMM_mnl_RW <- function(Data, Start_Val, Gamma_Sampler = NULL, Pr_Sampler = NULL, Gamma_Prior = NULL, Pr_Prior = NULL, D_Prior, Mcmc, Fast_vers = FALSE, Fast_vers_par = NULL){
	
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
	
# Initialize gamma sampler
	if(is.null(Gamma_Sampler)){
		int.mle0	<- mu_prop <- rep(0, m-2)
		scalar 		<- 2.93 / sqrt(m - 2)
	} else {
		int.mle0 	<- Gamma_Sampler$int.mle
		mu_prop 	<- Gamma_Sampler$mu_prop
		scalar 		<- Gamma_Sampler$scalar
	}
	
# Initialize pr sampler
	if(is.null(Pr_Sampler)){
		pr.int.mle0	<- mu_pr.prop <- rep(0, q_pr-1)
		pr.scalar 	<- 2.93 / sqrt(q_pr - 1)
	} else {
		pr.int.mle0	<- Pr_Sampler$pr.int.mle
		mu_pr.prop 	<- Pr_Sampler$mu_pr.prop
		pr.scalar 	<- Pr_Sampler$pr.scalar
	}
	
# Initialize Gamma prior
	if(is.null(Gamma_Prior)){
		mu_int0		<- rep(0, m-2)
		root0 		<- diag(m-2) * 0.01
	} else {
		mu_int0 	<- Gamma_Prior$mu_int0 
		root0 		<- Gamma_Prior$root0
	}
	sigma_int0 <- solve(root0)	
	
# Initialize Pr prior
	if(is.null(Pr_Prior)){
		mu_pr.int0	<- rep(0, q_pr-1)
		pr.root0 	<- diag(q_pr-1) * 0.01
	} else {
		mu_pr.int0 	<- Pr_Prior$mu_pr.int0 
		pr.root0 	<- Pr_Prior$pr.root0
	}
	sigma_pr.int0	<- solve(pr.root0)	

# Initialize D prior
	logmu0			<- D_Prior$logmu0
	tau20			<- D_Prior$tau20
	alpha.sigma20	<- D_Prior$alpha.sigma20
	beta.sigma20	<- D_Prior$beta.sigma20
	
# Initialize Mcmc argumetns
	J 				<- Mcmc$J
	burn.in			<- Mcmc$burn.in	
	
	
# Define arguments used in Mcmc algorithm	
	int <- int.new <- mle <- int.mle <- mhess <- candcov <- vector("list", m)
	pr.int <- pr.int.new <- pr.mle <- pr.int.mle <- pr.mhess <- pr.candcov <- vector("list", m)
	naccept <- rep(0,m)
	pr.naccept <- rep(0,m)
	tau2 <- logmu2 <- a1 <- b1 <- double(m)
	sample.path <- Dur <-numeric()
	sampled.states 	<- matrix(,nrow = n, ncol = J)
	Trans <- vector("list", m)	
		
# Create outcome matrix of the parameters
	PD 				<- matrix(, nrow = J, ncol = m*q_pr + m*m + m*2 + 1)
	colnames(PD) 	<- c(paste("pr",rep(1:q_pr, m),"_S",rep(1:m, each = q_pr), sep = ""), paste("S", rep(1:m, each = m), "to", 
					rep(1:m, m), sep = ""), paste("logmu", 1:m, sep = ""), paste("logvar", 1:m, sep = ""), "LL")
	PD[1, 1:((m*q_pr + m*m) + m*2)] <- Start_Val
	
# Define model parameters by starting values	
	pr 				<- matrix(Start_Val[1:(m*q_pr)], ncol = q_pr, nrow = m, byrow = TRUE)
	gamma 			<- matrix(Start_Val[(m*q_pr + 1) : (m*q_pr + m*m)], byrow = TRUE, ncol = m)
	delta 			<- solve(t(diag(m) - gamma + 1), rep(1, m))
	logmu			<- Start_Val[((m*q_pr + m*m) + 1) : ((m*q_pr + m*m) + m)]
	logsigma2		<- Start_Val[((m*q_pr + m*m) + m + 1) : ((m*q_pr + m*m) + m*2)]
	
	
###################################		
# Run the MCMC algorithm
	itime=proc.time()[3]
	for (iter in 2 : J){	
		
	# Run forward backward algorithm in C++, using the runlength distribution d for each state	
		d 				<- get.d.lognorm(run.p = list(logmu = logmu, logsd = sqrt(logsigma2)), Mx = Mx, m = m)
		if (Fast_vers == FALSE){
			FB				<- FBalgC(y2 = y2, m = m, pr = pr, Mx = Mx, Mx2 = Mx2, gamma = gamma, d = d, S2 = switch2, S = switch, delta = delta) 
		} else {
			FB				<- FBalgCF(y2 = y2, m = m, pr = pr, Mx = Mx, Mx2 = Mx2, gamma = gamma, d = d, S2 = switch2, S = switch, delta = delta, S3 = switch.sl)
		}
		B_star 			<- FB[[2]]
		Occupancy 		<- FB[[3]]
		N 				<- FB[[1]]		
		
	# Using the outcomes of the forward backward algorithm, sample the state sequence	
		sample.path 	<- Dur <- NULL
		Trans			<- vector("list", m)
		sample.path[1] 	<- sample(1:m, 1, prob = delta * exp(B_star[,1]))  
		Dur[1] 			<- sample(1:Mx, 1, prob = (Occupancy[[1]][sample.path[1],] / exp(B_star[sample.path[1],1])))		
		sampled.states[1:Dur[1], iter] <- sample.path[1]
		s <- 1
		while(sum(Dur) < n){
			s 				<- s + 1
			Mx.l 			<- min(Mx, n-sum(Dur))
			sample.path[s] 	<- sample(1:m, 1, prob = gamma[sample.path[s-1],] * exp(B_star[,sum(Dur)+1])) 
			Trans[[sample.path[s-1]]] <- c(Trans[[sample.path[s-1]]], sample.path[s])
			Dur[s]			<- sample(1:Mx.l, 1, prob = (Occupancy[[sum(Dur)+1]][sample.path[s],] / exp(B_star[sample.path[s], sum(Dur)+1])))	
			sampled.states[sum(Dur[1:s-1],1):sum(Dur), iter] <- sample.path[s]
		}
		n.Dur <- length(Dur)
			
	##############################################################
	# with the sample path, updata gamma, pr and duration distribution	

		for(i in 1:m){		
		# Update gamma and pr using Independence metropolis step 
				
			# renumber state sequence from 1 to m-1 such that the number of transition from state i to j = i is cut out  
			Trans.factor 		<- factor(Trans[[i]])	
			Trans.factor2		<- factor(Trans.factor, levels = paste(c(1:m)[-i])) 			
			Trans2 				<- as.numeric(Trans.factor2)
			Trans.factor 		<- factor(c(Trans[[i]], c(1:m)[-i]))	
			Trans.factor3		<- factor(Trans.factor, levels = paste(c(1:m)[-i])) 			
			Trans3				<- as.numeric(Trans.factor3)
			
			# obtain observations that are observed within sampled state i
			cond.y 				<- y[sampled.states[, iter] == i]					
	
			# obtain mle estimates for gamma
			mle[[i]]			<- optim(int.mle0, llmnl_int, Obs =  Trans3, n_cat = (m-1), method="BFGS",hessian=TRUE,control=list(fnscale=-1))
			int.mle[[i]]		<- mle[[i]]$par
			mhess[[i]]			<- mnlHess_int(int = int.mle[[i]], Obs = Trans3, n_cat =  (m-1))
			candcov[[i]]		<- chol2inv(chol(mhess[[i]]))	
			
			# obtain mle estimates for pr
			pr.mle[[i]]			<- optim(pr.int.mle0, llmnl_int, Obs = c(cond.y, c(1:q_pr)), n_cat = q_pr, method="BFGS",hessian=TRUE,control=list(fnscale=-1))
			pr.int.mle[[i]]		<- pr.mle[[i]]$par
			pr.mhess[[i]]		<- mnlHess_int(int = pr.int.mle[[i]], Obs = c(cond.y, c(1:q_pr)), n_cat =  q_pr)
			pr.candcov[[i]]		<- chol2inv(chol(pr.mhess[[i]]))
			# pr.candcov.all[iter, (1 + (i-1) * (q_pr-1)^2):((q_pr-1)^2 + (i-1) * (q_pr-1)^2)] <- as.vector(t(pr.candcov[[i]]))		
			
			
			# if this is the first run of the MCMC sampler, use mle values as start values for int
			if (iter == 2){
				int[[i]] 		<- int.mle[[i]]
				pr.int[[i]] 	<- pr.int.mle[[i]]
			}		
			
			RWout				<- mnl_RW_once(int1 = int[[i]], Obs = Trans2, n_cat = (m-1), mu_int_bar1 = mu_int0, V_int1 = sigma_int0, scalar = scalar, candcov1 = candcov[[i]])
			gamma[i, -i] 		<- RWout$prob
			naccept[i]			<- naccept[i] + RWout$accept
			int[[i]] 			<- RWout$draw.int
					
			pr.RWout			<- mnl_RW_once(int1 = pr.int[[i]], Obs = cond.y, n_cat = q_pr, mu_int_bar1 = mu_pr.int0, V_int1 = sigma_pr.int0, scalar = pr.scalar, candcov1 = pr.candcov[[i]])
			pr[i,]				<- pr.RWout$prob
			pr.naccept[i]		<- pr.naccept[i] + pr.RWout$accept
			pr.int[[i]] 		<- pr.RWout$draw.int
						
		# Update duration distribution using Gibbs step			
			tau2[i] 			<- 1/ ((1/tau20[i]) + (1/logsigma2[i]) * sum(sample.path[-n.Dur] == i))
			logmu2[i] 			<- tau2[i] * ((1/tau20[i]) * logmu0[i] + (1/logsigma2[i]) * sum(log(Dur[-n.Dur][sample.path[-n.Dur] == i])))
			logmu[i] 			<- rnorm(1, mean = logmu2[i], sd = sqrt(tau2[i]))	
			
			a1[i] 				<- alpha.sigma20[i] + sum(sample.path[-n.Dur] == i) / 2
			b1[i] 				<- beta.sigma20[i] + (1/2) * (sum((log(Dur[-n.Dur][sample.path[-n.Dur] == i]) - logmu[i])^2)) 
			logsigma2[i] 		<- 1 / rgamma(1, shape = a1[i], rate = b1[i])			
		}				
		delta 					<- solve(t(diag(m) - gamma + 1), rep(1, m)) 	
			
		# save values			 	
		PD[iter,(m*q_pr+1):(m*q_pr+m*m)] <- as.vector(t(gamma))
		PD[iter,1:(m*q_pr)] <- as.vector(t(pr))
		PD[iter,(m*q_pr+m*m+1):(m*q_pr+m*m+m)] <- logmu
		PD[iter,(m*q_pr+m*m+m+1):(m*q_pr+m*m+m*2)] <- logsigma2	
	
		# Save log likelihood of iteration
		PD[iter, (m*q_pr+m*m+m*2+1)] <- sum(N)
		print(c(iter, sum(N)))
	}
	
	ctime=proc.time()[3]
	print(paste("total time elapsed", round((ctime-itime)/60,2)))
	return(list(PD = PD, sampled.states = sampled.states, naccept = naccept, pr.naccept = pr.naccept))
}

#######################################
#######################################




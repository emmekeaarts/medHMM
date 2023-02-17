# Data : List of list, gen: containing m and q_pr, m.data: mouse data, one list per mouse containing y, Mx, Mx2

# pi_ (rows i of gamma) 	~(idd) 	Multinomial logit(int_i) 
# int_i 					~ N(mu_int_bar, V_int)
# mu_int_bar 				~ N(mu0, 1/K0 * V_int)
# V_int 					~ IW(nu, V)	
# Gamma_hyp_Prior list containing mu0, K0, nu and V

# H(theta_i)				~ Multinomial logit(pr_int_i) 
# pr_int					~ N(mu_pr.int_bar, V_pr.int)
# mu_pr.int_bar 			~ N(pr.mu0, 1/pr.K0 * V_pr.int)
# V_pr.int 					~ IW(pr.nu, pr.V)	
# Pr_hyp_Prior Prior list containing pr.mu0, pr.K0, pr.nu and pr.V

# G(w_Zs)					~ LogN(logmu, logsigma)
# logmu						~ N(mu_d_bar, tau2_d_bar)
# mu_d_bar					~ N(d_mu0, s2_0)
# tau2_d_bar				~ invGamma(alpha.tau20, beta.tau20)
# logsima					~ invGamma(alpha.simga20, beta.sigma20)
# D_hyp_Prior Prior list containing d_mu0, s2_0, alpha.tau20, beta.tau20, alpha.sigma20, beta.sigma20. Nota that logsigma is assumed invariant over mice.

# Gamma_Sampler list containing start values for mle estimates of pooled data for gamma, int.mle0, scalar and weight for the overall ll in the fractional likelihood, w
# Pr_Sampler list containing start values for mle estimates of pooled data for pr, pr.int.mle0, pr.scalar and weight for the overall ll in the fractional likelihood, pr.w
# Startvalues for logmu, logsigma (for first run FB), mu_d_bar, tau2_bar (to start gibbs sampler for duration distribution), gamma and pr (for first run FB), gamma and pr get their start values from the mle at the first run. Note: now programmed to set startvalues logmu = startvalues mu_d_bar, and startvalues logsigma = startvalues tau2_bar.  
# Mcmc: values for the mcmc algorithm, number of iterations J, burn in time, burn.in, and number of cores to be used for the paralel part, n.cores. 


Hier_HSMM_mnl_IndepM_par <- function(Data, Start_Val, Gamma_Sampler = NULL, Pr_Sampler = NULL, Gamma_hyp_Prior = NULL, Pr_hyp_Prior = NULL, D_hyp_Prior, Mcmc, Fast_vers = FALSE, Fast_vers_par = NULL){
	pandterm=function(message) { stop(message,call.=FALSE) }
	
# Initialize data	
	m.data 			<- Data$m.data
	n.mice			<- length(Data$m.data)
	ypooled			<- n <- switch <- switch2 <- NULL
	for(mi in 1:n.mice){
		ypooled 		<- c(ypooled, m.data[[mi]]$y)
		n				<- length(m.data[[mi]]$y)
		switch 			<- rep(0, n) 
		switch[1] 		<- 1
		for (ni in 2:n) {if(m.data[[mi]]$y[ni] != m.data[[mi]]$y[ni-1]) {switch[ni] <- 1}}
		switch2 		<- c(switch[-1],1)
		m.data[[mi]]	<- c(m.data[[mi]], list(y2 = m.data[[mi]]$y - 1, switch = switch, switch2 = switch2, n = n))
		m.data[[mi]]	<- c(m.data[[mi]], list(converge = numeric(m), int.mle = matrix(,m,(m-2)), mhess = matrix(,(m-2)*m, (m-2)), pr.converge = numeric(m), pr.int.mle = matrix(,q_pr,(q_pr-1)), pr.mhess = matrix(,(q_pr-1)*m, (q_pr-1))))	
		}	
	ypooled2		<- ypooled - 1
	n.total 		<- length(ypooled)
	m 				<- Data$gen$m
	q_pr 				<- Data$gen$q_pr

		if (Fast_vers == TRUE){
		if(is.null(Fast_vers_par)){
			act.sleep <- 6
			lag.sleep <- 60
		}	else {
			act.sleep	<- Fast_vers_par$act.sleep
			lag.sleep 	<- Fast_vers_par$lag.sleep 
		}		
		for (mi in 1:n.mice){
			if(sum(m.data[[mi]]$y == act.sleep) > 0){
				sleep.ind <- which(m.data[[mi]]$switch == 1 & m.data[[mi]]$y == act.sleep)
				sleep.ind2 <- numeric()
				for(k in 1:length(sleep.ind)){
					sleep.ind2 <- c(sleep.ind2, ((sleep.ind[k])-lag.sleep):sleep.ind[k])
				}
				sleep.ind2 <- sleep.ind2[sleep.ind2 > 0]
				sleep.ind3 <- rep(0,m.data[[mi]]$n)
				sleep.ind3[sleep.ind2] <- 1
				switch.sl <- as.numeric(m.data[[mi]]$switch == 1 & sleep.ind3 == 1)	
			} else {
				switch.sl <- rep(0,m.data[[mi]]$n)
			}			
			m.data[[mi]]	 <- c(m.data[[mi]], list(switch.sl = switch.sl))
		}			
	}

	
# Initialize gamma sampler
	if(is.null(Gamma_Sampler)){
		int.mle0	<- mu_prop <- rep(0, m-2)
		w			<- .1
		nu.M 		<- 3 + (m-2)
	} else {
		int.mle0 	<- Gamma_Sampler$int.mle0
		nu.M 		<- Gamma_Sampler$nu.M
		w 			<- Gamma_Sampler$w
	}
	
# Initialize pr sampler
	if(is.null(Pr_Sampler)){
		pr.int.mle0	<- mu_pr.prop <- rep(0, q_pr-1)
		pr.w		<- .1
		pr.nu.M 	<- 3 + (q_pr-1)
	} else {
		pr.int.mle0	<- Pr_Sampler$pr.int.mle0
		pr.nu.M 	<- Pr_Sampler$pr.nu.M
		pr.w		<- Pr_Sampler$pr.w
	}
	
# Initialize Gamma hyper prior
	if(is.null(Gamma_hyp_Prior)){
		mu0			<- rep(0, m-2)
		K0			<- 1
		nu			<- 3 + m-2
		V			<- nu * diag(m-2)	
	} else {
		mu0			<- Gamma_hyp_Prior$mu0
		K0			<- Gamma_hyp_Prior$K0
		nu			<- Gamma_hyp_Prior$nu
		V			<- Gamma_hyp_Prior$V
	}	
	
# Initialize Pr hyper prior
	if(is.null(Pr_hyp_Prior)){
		pr.mu0		<- matrix(rep(0, q_pr-1), ncol = q_pr-1, nrow = m)
		pr.K0		<- 1
		pr.nu		<- 3 + q_pr-1
		pr.V		<- pr.nu * diag(q_pr-1)	
	} else {
		pr.mu0		<- Pr_hyp_Prior$pr.mu0
		pr.K0		<- Pr_hyp_Prior$pr.K0
		pr.nu		<- Pr_hyp_Prior$pr.nu
		pr.V		<- Pr_hyp_Prior$pr.V
	}
	
# Initialize D hyper prior
	d_mu0			<- D_hyp_Prior$d_mu0
	s2_0			<- D_hyp_Prior$s2_0
	alpha.sigma20	<- D_hyp_Prior$alpha.sigma20
	beta.sigma20	<- D_hyp_Prior$beta.sigma20
	alpha.tau20		<- D_hyp_Prior$alpha.tau20
	beta.tau20		<- D_hyp_Prior$beta.tau20
	
# Initialize Mcmc argumetns
	if(is.null(Mcmc$n.cores)) {pandterm("Number of cores needs to be specified")}
	J 				<- Mcmc$J
	burn.in			<- Mcmc$burn.in	
	n.cores 		<- Mcmc$n.cores
	
	
# Define arguments used in Mcmc algorithm	
	d <- vector("list", n.mice)
	B_star <- Occupancy <- N <- numeric()
	sample.path <- Dur <- vector("list", n.mice)
	n.Dur <- numeric(n.mice)
	sampled.states 	<- rep(list(matrix(,nrow = n, ncol = J)), n.mice)
	Trans <- rep(list(vector("list", m)), n.mice)
	cond.y <- rep(list(vector("list", m)), n.mice)
	mle.pooled <-int.mle.pooled <- mhess.pooled <- pooled.ll <- vector("list", m)
	pr.mle.pooled <- pr.int.mle.pooled <- pr.mhess.pooled <- pr.pooled.ll <- vector("list", m)
	c_int <- rep(list(matrix(,n.mice, (m-2))), m)
	pr.c_int <- rep(list(matrix(,n.mice, (q_pr-1))), m)
	mu_int_bar <- V_int <- mu_pr.int_bar <- V_pr.int <- vector("list", m)
	mu.gamma.prob_bar <- rep(list(numeric(m)),m)
	mu.pr.prob_bar <- rep(list(numeric(q_pr)),m)
	naccept <- matrix(0, n.mice, m)
	pr.naccept <- matrix(0, n.mice, m)
	
	tau2 <- logmu2 <- a1 <- b1 <- double(m)

		
		
# Create outcome matrix of the parameters
	PD 					<- matrix(, nrow = J, ncol = m*q_pr + m*m + m + 1)
	colnames(PD) 		<- c(paste("pr",rep(1:q_pr, m),"_S",rep(1:m, each = q_pr), sep = ""), paste("S", rep(1:m, each = m), "to", 
					rep(1:m, m), sep = ""), paste("logmu", 1:m, sep = ""), "LL")
	PD[1, 1:((m*q_pr + m*m) + m)] <- Start_Val[1:((m*q_pr + m*m) + m)]
	PD.mice				<- rep(list(PD), n.mice)

	int.mice			<- rep(list(matrix(, nrow = J, ncol = (m-2) * m)), n.mice)
	int_bar				<- matrix(, nrow = J, ncol = (m-2) * m)
	pr.int_bar			<- matrix(, nrow = J, ncol = (q_pr-1) * m)
	pr.int.mice			<- rep(list(matrix(, nrow = J, ncol = (q_pr-1) * m)), n.mice)
	gamma.prob_bar		<- matrix(, nrow = J, ncol = (m * m))
	pr.prob_bar			<- matrix(, nrow = J, ncol = (q_pr * m))
	
	duration_bar		<- matrix(, nrow = J, ncol = m * 3)
	colnames(duration_bar)	<- c(paste("mu_d_bar", 1:m, sep = ""), paste("logsigma2", 1:m, sep = ""), paste("tau2_d_bar", 1:m, sep = ""))
		
# Define model parameters by starting values, start values here are for first FB run 	
	pr 				<- rep(list(matrix(Start_Val[1:(m*q_pr)], ncol = q_pr, nrow = m, byrow = TRUE)), n.mice) 
	gamma 			<- rep(list(matrix(Start_Val[(m*q_pr + 1) : (m*q_pr + m*m)], byrow = TRUE, ncol = m)), n.mice)
	delta 			<- rep(list(solve(t(diag(m) - gamma[[1]] + 1), rep(1, m))), n.mice)
	logmu			<- matrix(Start_Val[((m*q_pr + m*m) + 1) : ((m*q_pr + m*m) + m)], n.mice, m, byrow = TRUE)
	logsigma2		<- Start_Val[((m*q_pr + m*m) + m + 1) : ((m*q_pr + m*m) + m*2)]
	
# start values for duration gibbs sampler
	mu_d_bar		<- Start_Val[((m*q_pr + m*m) + 1) : ((m*q_pr + m*m) + m)]
	tau2_d_bar		<- Start_Val[((m*q_pr + m*m) + m + 1) : ((m*q_pr + m*m) + m*2)]

			
###################################		
# Run the MCMC algorithm
	itime=proc.time()[3]
	for (iter in 2 : J){	
		
		# for each mouse, obtain sampled state sequence with mouse individual parameters
		sample.path <- Dur <- vector("list", n.mice)		
		for(mi in 1:n.mice){
			d[[mi]] 		<- get.d.lognorm(run.p = list(logmu = logmu[mi,], logsd = sqrt(logsigma2)), Mx = m.data[[mi]]$Mx, m = m)
		}
			
		# Run forward backward algorithm in parralel C++, using the runlength distribution d for each state
		cl <- makeCluster(n.cores)
		registerDoParallel(cl)
		if (Fast_vers == FALSE){
			FB <- foreach(mi = 1:n.mice) %dopar% {
			FBalgC(y2 = m.data[[mi]]$y2, m = m, pr = pr[[mi]], Mx = m.data[[mi]]$Mx, Mx2 = m.data[[mi]]$Mx2, gamma = gamma[[mi]], d = d[[mi]], S2 = m.data[[mi]]$switch2, S = m.data[[mi]]$switch, delta = delta[[mi]])} 
		} else {
			FB <- foreach(mi = 1:n.mice) %dopar% {
			FBalgCF(y2 = m.data[[mi]]$y2, m = m, pr = pr[[mi]], Mx = m.data[[mi]]$Mx, Mx2 = m.data[[mi]]$Mx2, gamma = gamma[[mi]], d = d[[mi]], S2 = m.data[[mi]]$switch2, S = m.data[[mi]]$switch, delta = delta[[mi]], S3 = m.data[[mi]]$switch.sl)}
		}
		stopCluster(cl)
		
		for(mi in 1:n.mice){
			B_star					<- FB[[mi]][[2]]
			Occupancy				<- FB[[mi]][[3]]
			N 						<- FB[[mi]][[1]]
			PD.mice[[mi]][iter, m*q_pr + m*m + m + 1] 	<- sum(N)	# adjust index	
		
		# Using the outcomes of the forward backward algorithm, sample the state sequence	
			Trans[[mi]]				<- vector("list", m)
			sample.path[[mi]][1] 	<- sample(1:m, 1, prob = delta[[mi]] * exp(B_star[,1]))  
			Dur[[mi]][1] 			<- sample(1:m.data[[mi]]$Mx, 1, prob = (Occupancy[[1]][sample.path[[mi]][1],] / exp(B_star[sample.path[[mi]][1],1])))		
			sampled.states[[mi]][1:Dur[[mi]][1], iter] <- sample.path[[mi]][1]
			s <- 1
			while(sum(Dur[[mi]]) < m.data[[mi]]$n){
				s 						<- s + 1
				Mx.l 					<- min(m.data[[mi]]$Mx, m.data[[mi]]$n-sum(Dur[[mi]]))
				sample.path[[mi]][s] 	<- sample(1:m, 1, prob = gamma[[mi]][sample.path[[mi]][s-1],] * exp(B_star[,sum(Dur[[mi]])+1])) 
				Trans[[mi]][[sample.path[[mi]][s-1]]] <- c(Trans[[mi]][[sample.path[[mi]][s-1]]], sample.path[[mi]][s])
				Dur[[mi]][s]			<- sample(1:Mx.l, 1, prob = (Occupancy[[sum(Dur[[mi]])+1]][sample.path[[mi]][s],] / exp(B_star[sample.path[[mi]][s], sum(Dur[[mi]])+1])))
				sampled.states[[mi]][sum(Dur[[mi]][1:s-1],1):sum(Dur[[mi]]), iter] <- sample.path[[mi]][s]
			}
			n.Dur[mi] <- length(Dur[[mi]])
			for (i in 1:m){
				cond.y[[mi]][[i]] <- m.data[[mi]]$y[sampled.states[[mi]][, iter] == i]						
			}		
		}

			
	##############################################################
 
		for(i in 1:m){		
			
			#################
			# Obtain hierarchical and mouse specific parameters for gamma and pr using metropolis sampler
			#################	
			
			# First obtain optim on pooled data 
			Trans.pooled.factor		<- factor(c(unlist(sapply(Trans, "[[", i)), c(1:m)[-i]))
			Trans.pooled.factor2	<- factor(Trans.pooled.factor, levels = paste(c(1:m)[-i])) 
			Trans.pooled2			<- as.numeric(Trans.pooled.factor2)				
			mle.pooled[[i]] 		<- optim(int.mle0, llmnl_int, Obs = Trans.pooled2, n_cat = (m-1), method="BFGS",hessian=TRUE,control=list(fnscale=-1))
			int.mle.pooled[[i]] 	<- mle.pooled[[i]]$par
			pooled.ll[[i]]			<- mle.pooled[[i]]$value
	
			cond.y.pooled			<- unlist(sapply(cond.y, "[[", i))	
			pr.mle.pooled[[i]]		<- optim(pr.int.mle0, llmnl_int, Obs = c(cond.y.pooled, c(1:q_pr)), n_cat = q_pr, method="BFGS",hessian=TRUE,control=list(fnscale=-1))
			pr.int.mle.pooled[[i]]	<- pr.mle.pooled[[i]]$par
			pr.pooled.ll[[i]]		<- pr.mle.pooled[[i]]$value

			
			# Obtain optim for each mouse specific	
			for (mi in 1:n.mice){
				# renumber state sequence from 1 to m-1 such that the number of transition from state i to j = i is cut out  
				Trans.factor 		<- factor(c(Trans[[mi]][[i]], c(1:m)[-i]))	
				Trans.factor3		<- factor(Trans.factor, levels = paste(c(1:m)[-i])) 			
				Trans3				<- as.numeric(Trans.factor3)
								
				wgt 				<- m.data[[mi]]$n / n.total
				out					<- optim(int.mle.pooled[[i]], llmnl_int_frac, Obs = Trans3, n_cat = (m-1), pooled.likel = pooled.ll[[i]], w = w, wgt = wgt, method="BFGS",hessian=TRUE,control=list(fnscale=-1))
				if(out$convergence == 0){
					m.data[[mi]]$converge[i]									<- 1											
					m.data[[mi]]$mhess[(1+(i-1)*(m-2)):((m-2)+(i-1)*(m-2)),]	<- mnlHess_int(int = out$par, Obs = Trans3, n_cat =  (m-1))
					m.data[[mi]]$int.mle[i,]									<- out$par
				} else {
					m.data[[mi]]$converge[i]									<- 0
					m.data[[mi]]$mhess[(1+(i-1)*(m-2)):((m-2)+(i-1)*(m-2)),]	<- diag(m-2)
					m.data[[mi]]$int.mle[i,]									<- rep(0, m-2)
				}
				pr.out				<- optim(pr.int.mle.pooled[[i]], llmnl_int_frac, Obs = c(cond.y[[mi]][[i]], c(1:q_pr)), n_cat = q_pr, pooled.likel = pr.pooled.ll[[i]], w = pr.w, wgt = wgt, method="BFGS", hessian=TRUE, control=list(fnscale=-1))
				if(pr.out$convergence == 0){
					m.data[[mi]]$pr.converge[i]									<- 1
					m.data[[mi]]$pr.mhess[(1+(i-1)*(q_pr-1)):((q_pr-1)+(i-1)*(q_pr-1)),]	<- mnlHess_int(int = pr.out$par, Obs = c(cond.y[[mi]][[i]], c(1:q_pr)), n_cat =  q_pr)
					m.data[[mi]]$pr.int.mle[i,]									<- pr.out$par
				} else {
					m.data[[mi]]$pr.converge[i]									<- 0
					m.data[[mi]]$pr.mhess[(1+(i-1)*(m-2)):((m-2)+(i-1)*(m-2)),]	<- diag(q_pr-1)
					m.data[[mi]]$pr.int.mle[i,]									<- rep(0, q_pr-1)
				}
								
				if (iter == 2){
					c_int[[i]][mi,]		<- out$par
					pr.c_int[[i]][mi,]	<- pr.out$par
				}
			}	
			
			# draw hierarchical parameters for for gamma and pr			
			ybar 				<- apply(c_int[[i]], 2, mean)
			mu0.n 				<- n.mice*ybar / (n.mice+K0) + K0*mu0 / (n.mice + K0)
			V.n					<- V + (n.mice - 1)*var(c_int[[i]]) + K0*n.mice*(ybar - mu0)%*%t(ybar - mu0) / (K0 + n.mice)
			V_int[[i]]			<- solve(rwish(S = solve(V.n), v = nu + n.mice))
			mu_int_bar[[i]] 	<- rmvnorm(1, mean = mu0.n, sigma = V_int[[i]] / (K0 + n.mice))
			exp.int				<- matrix(exp(c(0, mu_int_bar[[i]] )), nrow  = 1)
			mu.gamma.prob_bar[[i]][-i] 	<- exp.int / as.vector(exp.int %*% c(rep(1,(m-1))))
			
			pr.ybar 			<- apply(pr.c_int[[i]], 2, mean)
			pr.mu0.n 			<- n.mice*pr.ybar / (n.mice+pr.K0) + pr.K0*pr.mu0[i,] / (n.mice + pr.K0)
			pr.V.n				<- pr.V + (n.mice - 1)*var(pr.c_int[[i]]) + pr.K0*n.mice*(pr.ybar - pr.mu0[i,])%*%t(pr.ybar - pr.mu0[i,]) / (pr.K0 + n.mice)
			V_pr.int[[i]]		<- solve(rwish(S =solve(pr.V.n), v = pr.nu + n.mice))
			mu_pr.int_bar[[i]] 	<- rmvnorm(1, mean = pr.mu0.n, sigma = V_pr.int[[i]] / (pr.K0 + n.mice))
			exp.pr.int			<- matrix(exp(c(0, mu_pr.int_bar[[i]])), nrow  = 1)
			mu.pr.prob_bar[[i]] <- exp.pr.int / as.vector(exp.pr.int %*% c(rep(1,(q_pr))))
			
			# draw mouse individual paramers, based on mouse data and hierarchical distribution, for both gamma and pr. 
			for (mi in 1:n.mice){
				Trans.factor 			<- factor(Trans[[mi]][[i]])	
				Trans.factor2			<- factor(Trans.factor, levels = paste(c(1:m)[-i])) 			
				Trans2 					<- as.numeric(Trans.factor2)
				
				candcov.comb 			<- chol2inv(chol(m.data[[mi]]$mhess[(1+(i-1)*(m-2)):((m-2)+(i-1)*(m-2)),] + chol2inv(chol(V_int[[i]]))))
				int.star.comb			<- as.vector(candcov.comb %*% (m.data[[mi]]$mhess[(1+(i-1)*(m-2)):((m-2)+(i-1)*(m-2)),] %*% matrix(m.data[[mi]]$int.mle[i,], ncol = 1) + chol2inv(chol(V_int[[i]])) %*% matrix(mu_int_bar[[i]], ncol = 1)) )
				IndepMout				<- mnl_IndepM_once(int1 = c_int[[i]][mi,], Obs = Trans2, n_cat = (m-1), mu_int_bar1 = mu_int_bar[[i]], V_int1 = V_int[[i]], nu.M1 = nu.M, candcov1 = candcov.comb, int.star1 = int.star.comb)
				gamma[[mi]][i, -i]  	<- PD.mice[[mi]][iter, c((m*q_pr + 1+(i-1)*m):(m*q_pr + (i-1)*m + m))[-i]] <- IndepMout$prob
				naccept[mi, i]			<- naccept[mi, i] + IndepMout$accept
				c_int[[i]][mi,]		 	<- IndepMout$draw.int

				pr.candcov.comb			<- chol2inv(chol(m.data[[mi]]$pr.mhess[(1+(i-1)*(q_pr-1)):((q_pr-1)+(i-1)*(q_pr-1)),] + chol2inv(chol(V_pr.int[[i]]))))
				pr.int.star.comb		<- as.vector(pr.candcov.comb %*% (m.data[[mi]]$pr.mhess[(1+(i-1)*(q_pr-1)):((q_pr-1)+(i-1)*(q_pr-1)),]	%*% matrix(m.data[[mi]]$pr.int.mle[i,], ncol = 1) + chol2inv(chol(V_pr.int[[i]])) %*% matrix(mu_pr.int_bar[[i]], ncol = 1)))
				pr.IndepMout			<- mnl_IndepM_once(int1 = pr.c_int[[i]][mi,], Obs = cond.y[[mi]][[i]], n_cat = q_pr, mu_int_bar1 = mu_pr.int_bar[[i]], V_int1 = V_pr.int[[i]], nu.M1 = pr.nu.M, candcov1 = pr.candcov.comb, int.star1 = pr.int.star.comb)
				pr[[mi]][i,]			<- PD.mice[[mi]][iter, (1+(i-1)*q_pr):((i-1)*q_pr + q_pr)] <- pr.IndepMout$prob
				pr.naccept[mi, i]		<- pr.naccept[mi, i] + pr.IndepMout$accept
				pr.c_int[[i]][mi,]		<- pr.IndepMout$draw.int				
			
				int.mice[[mi]][iter, (1 + (i-1)*(m-2)) : ((m-2) + (i-1)*(m-2))]		<- c_int[[i]][mi,] #check index 
				pr.int.mice[[mi]][iter, (1 + (i-1)*(q_pr-1)) : ((q_pr-1) + (i-1)*(q_pr-1))]	<- pr.c_int[[i]][mi,] #check index 
				
				if(i == m){
					delta[[mi]] 		<- solve(t(diag(m) - gamma[[mi]] + 1), rep(1, m)) 
				}
			
			}
		
			# ADD TO SCRIPT: log likelihood of data given parameters is the sum of the log likelihoods over all mice, can do seperately for gamma and pr (and distribution), think that complete log likelihood is sum over these, how does this compare with N obtained in forward backward? Sum these over all mice to get the complete log likelihood? 

			#################
			# Obtain hierarchical and mouse specific parameters for duration distribuiton using gibbs sampler
			#################
		
			#draw logmu's
			for(mi in 1:n.mice){		
				tau2 		<- 1/ ((1/tau2_d_bar[i]) + (1/logsigma2[i]) * sum(sample.path[[mi]][-n.Dur[mi]] == i))
				logmu[mi,i]	<- rnorm(1,
					mean = tau2 * ((1/tau2_d_bar[i]) * mu_d_bar[i] + (1/logsigma2[i]) * sum(log(Dur[[mi]][-n.Dur[mi]][sample.path[[mi]][-n.Dur[mi]] == i]))),
					sd = sqrt(tau2))		
				PD.mice[[mi]][iter, m*m + m*q_pr + i]			<- logmu[mi, i] # adjust index
			}
			
			# draw mu_d_bar
			s2				<- 1 / ((1/s2_0[i]) + (1/tau2_d_bar[i]) * n.mice)
			mu_d_bar[i] 	<- rnorm(1, 
				mean = s2 * ((1/s2_0[i]) * d_mu0[i] + (1/tau2_d_bar[i]) * sum(logmu[,i])), 
				sd = sqrt(s2))

			# draw tau2_d_bar
			a1 <- alpha.tau20 + n.mice/2
			b1 <- beta.tau20 + (1/2) * (sum((logmu[,i] - mu_d_bar[i])^2)) 
			tau2_d_bar[i] <- 1 / rgamma(1, shape = a1, rate = b1)			
			
			#draw logsigma
			ss <- numeric(1)
			n.ss <- numeric(1)
			for (mi in 1:n.mice){
				ss <- ss + sum((log(Dur[[mi]][-n.Dur[mi]][sample.path[[mi]][-n.Dur[mi]] == i]) - logmu[mi,i])^2)
				n.ss <- n.ss + sum(sample.path[[mi]][-n.Dur[mi]] == i)	
			}			
			c1 <- alpha.sigma20 + n.ss/2 
			d1 <- beta.sigma20 + (1/2) * ss 
			logsigma2[i] <- 1 / rgamma(1, shape = c1, rate = d1)	
				# compare n.ss with length(unlist(sapply(Trans, "[[", i)))			
					
		}
		###############################################################	
				
		# save global values
		duration_bar[iter, ]			<- c(mu_d_bar, logsigma2, tau2_d_bar) # logsima is fixed over mice, so only need to record general 		
		int_bar[iter, ]					<- unlist(mu_int_bar)
		pr.int_bar[iter, ]				<- unlist(mu_pr.int_bar)
		gamma.prob_bar[iter,]			<- unlist(mu.gamma.prob_bar)	
		pr.prob_bar[iter,]				<- unlist(mu.pr.prob_bar)
		
		print(c(iter))
	}
	
	ctime=proc.time()[3]
	print(paste("total time elapsed", round((ctime-itime)/60,2)))
	return(list(PD.mice = PD.mice, int.mice = int.mice, pr.int.mice = pr.int.mice, duration_bar = duration_bar, int_bar = int_bar, pr.int_bar = pr.int_bar, gamma.prob_bar = gamma.prob_bar, pr.prob_bar = pr.prob_bar, naccept = naccept, pr.naccept = pr.naccept, sampled.states = sampled.states))
}

#######################################
#######################################



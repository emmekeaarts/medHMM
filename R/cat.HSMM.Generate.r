# simulate HSMM sequence (hidden states and observations with rd runlength distribution for state duration)
cat.HSMM.Generate <- function(n, pr, gamma, rd, run.p, Mx = NA, seed = NULL, delta = NA){
  if (!is.null(seed)){ 
    set.seed(seed)
    }	
	m <- dim(gamma)[1]
	if (is.na(Mx)){
    	Mx <- as.integer(max(n, 1000))
    }
    n.pr <- dim(pr)[2]
    state <- y <- c()	
#	if (rd == "poiss"){ # generate runlength distribution
#		RL.gen <- get.d.poiss(run.p = run.p, Mx = Mx, m = m)
#	} else if (rd == "nbinom"){
#		RL.gen <- get.d.nbinom(run.p = run.p, Mx = Mx, m = m)	
	#} else 
	if (rd == "lognorm") {
		RL.gen <- get.d.lognorm(run.p = run.p, Mx = Mx, m = m)
	}
	st <- numeric()
	if (is.na(delta)) { # first state
		st <- sample(1:m, 1, prob = solve(t(diag(m) - gamma + 1), rep(1, m)))
	} else {
		st <- sample(1:m, 1, prob = delta)
	}		
	t <- 0
	z <- st
	dur.z <- numeric()
	while (t < n){
		dur.state <- sample(c(1:length(RL.gen[st, RL.gen[st,]!=0])), prob=RL.gen[st, RL.gen[st,]!=0], 1) # sample duration of the state
		y <- c(y, sample(1:dim(pr)[2], dur.state, prob = pr[st,], replace = TRUE)) # generation of observations with lenght dur.state		
		state <- c(state, rep(st, dur.state)) #put the state with lenght dur.state in state vector 
		st <- sample(1:m, 1, prob = gamma[st,]) # sample the new state
		z <- c(z, st)
		dur.z <- c(dur.z, dur.state)
		t <- t + dur.state
	}
	state <- state[1:n]
	y <- y[1:n]
	out <- list(y = y, state = state, z=z, dur.z = dur.z)
	return(out)
}


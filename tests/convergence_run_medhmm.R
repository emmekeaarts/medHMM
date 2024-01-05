#-------------------------------------------------------------------------------
# Convergence checks script for mhmm
#-------------------------------------------------------------------------------
#library(devtools)
#devtools::install_github("https://github.com/emmekeaarts/medHMM/tree/dev")
library(medHMM)
library(pbmcapply)
library(tidyverse)
library(R.utils)
cores<-1
#------------------------------------------------------------------------------
#Gelman rubin
gelman_rubin <- function(par_matrix, J=J,burn_in=burn_in){
  chains<-2
  samples<-J-burn_in
  # Coerce to matrix
  par_matrix <- as.data.frame(par_matrix)
  out<-c()
  for(param in 2:ncol(par_matrix)){
    # Mean over all samples
    all_mean <- mean(par_matrix[,param])

    # Mean of each chain
    chain_mean <- tapply(par_matrix[,param], par_matrix[,1], mean)

    # Variance of each chain
    chain_var <- tapply(par_matrix[,param], par_matrix[,1], stats::var)
    W <- (1 / chains) * sum(chain_var)
    B <- samples / (chains - 1) * sum((chain_mean - all_mean)^2)
    V <- (1 - 1 / samples) * W + (1 / samples) * B
    out[param-1]<-round(sqrt(V / W), 4)
  }
  return(round(out,1))
}

ge_rub_two_models_emp<-function(model1=NULL, model2=NULL, m=2, J=4000,burn_in=2000){
  if(m>2){
    chains_matrix_gamma<-matrix(ncol = (m^2-(2*m)))
    chains_matrix_gamma<-rbind(chains_matrix_gamma,model1$gamma_int_bar[(burn_in+1):J,1:(m^2-(2*m))],model2$gamma_int_bar[(burn_in+1):J,1:(m^2-(2*m))])
    chains_matrix_gamma<-na.omit(chains_matrix_gamma)
  }
  var_data<-rbind(model1$emiss_var_bar[[1]][(burn_in+1):(J),],model2$emiss_var_bar[[1]][(burn_in+1):(J),])
  mean_data<-rbind(model1$emiss_mu_bar[[1]][(burn_in+1):J,1:m],model2$emiss_mu_bar[[1]][(burn_in+1):J,1:m])
  var_data2<-rbind(model1$emiss_var_bar[[2]][(burn_in+1):(J),],model2$emiss_var_bar[[2]][(burn_in+1):(J),])
  mean_data2<-rbind(model1$emiss_mu_bar[[2]][(burn_in+1):J,1:m],model2$emiss_mu_bar[[2]][(burn_in+1):J,1:m])
  chains_matrix_emiss<-cbind(mean_data,var_data,mean_data2,var_data2)
  chains_matrix_emiss<-na.omit(chains_matrix_emiss)
  mean_dwell<-rbind(model1$dwell_mu_bar[(burn_in+1):J,1:m],model2$dwell_mu_bar[(burn_in+1):J,1:m])
  var_dwell<-rbind(model1$dwell_var_bar[(burn_in+1):(J),],model2$dwell_var_bar[(burn_in+1):(J),])
  chains_matrix_dwell<-cbind(mean_dwell,var_dwell)
  index<-c(rep(1,J-burn_in),rep(2,J-burn_in))
  if(m==2){
    gelman_input<-cbind(index,chains_matrix_emiss,chains_matrix_dwell)
  }else{
    gelman_input<-cbind(index,chains_matrix_emiss,chains_matrix_dwell,chains_matrix_gamma)
  }

  out_scenario_gelman_empirical<-gelman_rubin(gelman_input,J = J,burn_in = burn_in)
  out_scenario_gelman_empirical<-t(as.matrix(out_scenario_gelman_empirical))
  colnames(out_scenario_gelman_empirical)<-colnames(gelman_input)[-1]
  return(out_scenario_gelman_empirical)
}

MAP_medhmm<-function(case_out=NULL,iteration=1,J=NULL,B=NULL,m=NULL){
  #---emiss_mu_bar---------
  emiss_mu_bar_list<-list()
  for(i in 1:length(case_out[["emiss_mu_bar"]])){
    median<-case_out[["emiss_mu_bar"]][[i]][(B+1):J,] %>% apply(.,2,median) %>% t()
    colnames(median)<-paste0(colnames(median),"_median")
    mean<-case_out[["emiss_mu_bar"]][[i]][(B+1):J,] %>% apply(.,2,mean) %>% t()
    colnames(mean)<-paste0(colnames(mean),"_mean")
    sd<-case_out[["emiss_mu_bar"]][[i]][(B+1):J,] %>% apply(.,2,sd) %>% t()
    colnames(sd)<-paste0(colnames(sd),"_sd")
    ci_025<-case_out[["emiss_mu_bar"]][[i]][(B+1):J,] %>% apply(.,2,FUN = function(x){quantile(x,0.025)}) %>% t()
    ci_975<-case_out[["emiss_mu_bar"]][[i]][(B+1):J,] %>% apply(.,2,FUN = function(x){quantile(x,0.975)}) %>% t()
    colnames(ci_025)<-paste0(colnames(ci_025),"_ci_025")
    colnames(ci_975)<-paste0(colnames(ci_025),"_ci_975")
    emiss_mu_bar<-cbind(median,mean,sd,ci_025,ci_975)
    rownames(emiss_mu_bar)<-paste0("emiss_mu_bar",iteration)
    emiss_mu_bar_list[[i]]<-emiss_mu_bar
  }
  #---emiss_var_bar--------
  emiss_var_bar_list<-list()
  for(i in 1:length(case_out[["emiss_var_bar"]])){
    median<-case_out[["emiss_var_bar"]][[i]][(B+1):J,] %>% apply(.,2,median) %>% t()
    colnames(median)<-paste0(colnames(median),"_median")
    mean<-case_out[["emiss_var_bar"]][[i]][(B+1):J,] %>% apply(.,2,mean) %>% t()
    colnames(mean)<-paste0(colnames(mean),"_mean")
    sd<-case_out[["emiss_var_bar"]][[i]][(B+1):J,] %>% apply(.,2,sd) %>% t()
    colnames(sd)<-paste0(colnames(sd),"_sd")
    ci_025<-case_out[["emiss_var_bar"]][[i]][(B+1):J,] %>% apply(.,2,FUN = function(x){quantile(x,0.025)}) %>% t()
    ci_975<-case_out[["emiss_var_bar"]][[i]][(B+1):J,] %>% apply(.,2,FUN = function(x){quantile(x,0.975)}) %>% t()
    colnames(ci_025)<-paste0(colnames(ci_025),"_ci_025")
    colnames(ci_975)<-paste0(colnames(ci_975),"_ci_975")
    emiss_var_bar<-cbind(median,mean,sd,ci_025,ci_975)
    rownames(emiss_var_bar)<-paste0("emiss_var_bar",iteration)
    emiss_var_bar_list[[i]]<-emiss_var_bar
  }
  #---emiss_varmu_bar------
  emiss_varmu_bar_list<-list()
  for(i in 1:length(case_out[["emiss_varmu_bar"]])){
    median<-case_out[["emiss_varmu_bar"]][[i]][(B+1):J,] %>% apply(.,2,median) %>% t()
    colnames(median)<-paste0(colnames(median),"_median")
    mean<-case_out[["emiss_varmu_bar"]][[i]][(B+1):J,] %>% apply(.,2,mean) %>% t()
    colnames(mean)<-paste0(colnames(mean),"_mean")
    sd<-case_out[["emiss_varmu_bar"]][[i]][(B+1):J,] %>% apply(.,2,sd) %>% t()
    colnames(sd)<-paste0(colnames(sd),"_sd")
    ci_025<-case_out[["emiss_varmu_bar"]][[i]][(B+1):J,] %>% apply(.,2,FUN = function(x){quantile(x,0.025)}) %>% t()
    ci_975<-case_out[["emiss_varmu_bar"]][[i]][(B+1):J,] %>% apply(.,2,FUN = function(x){quantile(x,0.975)}) %>% t()
    colnames(ci_025)<-paste0(colnames(ci_025),"_ci_025")
    colnames(ci_975)<-paste0(colnames(ci_975),"_ci_975")
    emiss_varmu_bar<-cbind(median,mean,sd,ci_025,ci_975)
    rownames(emiss_varmu_bar)<-paste0("emiss_varmu_bar",iteration)
    emiss_varmu_bar_list[[i]]<-emiss_varmu_bar
  }
  #---dwell_mu_bar---------
  median<-case_out[["dwell_mu_bar"]][(B+1):J,] %>% apply(.,2,median) %>% t()
  colnames(median)<-paste0(colnames(median),"_median")
  mean<-case_out[["dwell_mu_bar"]][(B+1):J,] %>% apply(.,2,mean) %>% t()
  colnames(mean)<-paste0(colnames(mean),"_mean")
  sd<-case_out[["dwell_mu_bar"]][(B+1):J,] %>% apply(.,2,sd) %>% t()
  colnames(sd)<-paste0(colnames(sd),"_sd")
  ci_025<-case_out[["dwell_mu_bar"]][(B+1):J,] %>% apply(.,2,FUN = function(x){quantile(x,0.025)}) %>% t()
  ci_975<-case_out[["dwell_mu_bar"]][(B+1):J,] %>% apply(.,2,FUN = function(x){quantile(x,0.975)}) %>% t()
  colnames(ci_025)<-paste0(colnames(ci_025),"_ci_025")
  colnames(ci_975)<-paste0(colnames(ci_025),"_ci_975")
  dwell_mu_bar<-cbind(median,mean,sd,ci_025,ci_975)
  rownames(dwell_mu_bar)<-paste0("dwell_mu_bar",iteration)

  #---dwell_var_bar--------
  median<-case_out[["dwell_var_bar"]][(B+1):J,] %>% apply(.,2,median) %>% t()
  colnames(median)<-paste0(colnames(median),"_median")
  mean<-case_out[["dwell_var_bar"]][(B+1):J,] %>% apply(.,2,mean) %>% t()
  colnames(mean)<-paste0(colnames(mean),"_mean")
  sd<-case_out[["dwell_var_bar"]][(B+1):J,] %>% apply(.,2,sd) %>% t()
  colnames(sd)<-paste0(colnames(sd),"_sd")
  ci_025<-case_out[["dwell_var_bar"]][(B+1):J,] %>% apply(.,2,FUN = function(x){quantile(x,0.025)}) %>% t()
  ci_975<-case_out[["dwell_var_bar"]][(B+1):J,] %>% apply(.,2,FUN = function(x){quantile(x,0.975)}) %>% t()
  colnames(ci_025)<-paste0(colnames(ci_025),"_ci_025")
  colnames(ci_975)<-paste0(colnames(ci_975),"_ci_975")
  dwell_var_bar<-cbind(median,mean,sd,ci_025,ci_975)
  rownames(dwell_var_bar)<-paste0("dwell_var_bar",iteration)

  #---dwell_varmu_bar------
  median<-case_out[["dwell_varmu_bar"]][(B+1):J,] %>% apply(.,2,median) %>% t()
  colnames(median)<-paste0(colnames(median),"_median")
  mean<-case_out[["dwell_varmu_bar"]][(B+1):J,] %>% apply(.,2,mean) %>% t()
  colnames(mean)<-paste0(colnames(mean),"_mean")
  sd<-case_out[["dwell_varmu_bar"]][(B+1):J,] %>% apply(.,2,sd) %>% t()
  colnames(sd)<-paste0(colnames(sd),"_sd")
  ci_025<-case_out[["dwell_varmu_bar"]][(B+1):J,] %>% apply(.,2,FUN = function(x){quantile(x,0.025)}) %>% t()
  ci_975<-case_out[["dwell_varmu_bar"]][(B+1):J,] %>% apply(.,2,FUN = function(x){quantile(x,0.975)}) %>% t()
  colnames(ci_025)<-paste0(colnames(ci_025),"_ci_025")
  colnames(ci_975)<-paste0(colnames(ci_975),"_ci_975")
  dwell_varmu_bar<-cbind(median,mean,sd,ci_025,ci_975)
  rownames(dwell_varmu_bar)<-paste0("dwell_varmu_bar",iteration)

  if(m>2){
    #---gamma_prob_bar------
    median<-case_out[["gamma_prob_bar"]][(B+1):J,] %>% apply(.,2,median) %>% t()
    colnames(median)<-paste0(colnames(median),"_median")
    mean<-case_out[["gamma_prob_bar"]][(B+1):J,] %>% apply(.,2,mean) %>% t()
    colnames(mean)<-paste0(colnames(mean),"_mean")
    sd<-case_out[["gamma_prob_bar"]][(B+1):J,] %>% apply(.,2,sd) %>% t()
    colnames(sd)<-paste0(colnames(sd),"_sd")
    ci_025<-case_out[["gamma_prob_bar"]][(B+1):J,] %>% apply(.,2,FUN = function(x){quantile(x,0.025)}) %>% t()
    ci_975<-case_out[["gamma_prob_bar"]][(B+1):J,] %>% apply(.,2,FUN = function(x){quantile(x,0.975)}) %>% t()
    colnames(ci_025)<-paste0(colnames(ci_025),"_ci_025")
    colnames(ci_975)<-paste0(colnames(ci_975),"_ci_975")
    gamma_prob_bar<-cbind(median,mean,sd,ci_025,ci_975)
    rownames(gamma_prob_bar)<-paste0("gamma_prob_bar",iteration)

    #---gamma_int_prob_bar------
    median<-case_out[["gamma_int_bar"]][(B+1):J,] %>% apply(.,2,median) %>% t()
    colnames(median)<-paste0(colnames(median),"_median")
    mean<-case_out[["gamma_int_bar"]][(B+1):J,] %>% apply(.,2,mean) %>% t()
    colnames(mean)<-paste0(colnames(mean),"_mean")
    sd<-case_out[["gamma_int_bar"]][(B+1):J,] %>% apply(.,2,sd) %>% t()
    colnames(sd)<-paste0(colnames(sd),"_sd")
    ci_025<-case_out[["gamma_int_bar"]][(B+1):J,] %>% apply(.,2,FUN = function(x){quantile(x,0.025)}) %>% t()
    ci_975<-case_out[["gamma_int_bar"]][(B+1):J,] %>% apply(.,2,FUN = function(x){quantile(x,0.975)}) %>% t()
    colnames(ci_025)<-paste0(colnames(ci_025),"_ci_025")
    colnames(ci_975)<-paste0(colnames(ci_975),"_ci_975")
    gamma_int_bar<-cbind(median,mean,sd,ci_025,ci_975)
    rownames(gamma_int_bar)<-paste0("gamma_int_bar",iteration)

    #--gamma_V_int_bar-----------
    #subject specific parameters of variance covariance between
    all_st<-c()
    for(st_cur in 1:m){
      to_st<-2:m
      to_st<-to_st[which(to_st!=st_cur)]
      all_st<-append(all_st,paste0("S",st_cur,"toS",to_st))
    }
    names_need<-c()
    for(al in 1:length(all_st)){
      names_need<-append(names_need,paste0("var_int_",all_st[al],"_with_int_",all_st[al]))
    }
    names_need<-names_need[-1]
    case_out_gamma_V_int_bar<-case_out[["gamma_V_int_bar"]][(B+1):J,names_need]
    median<-case_out_gamma_V_int_bar %>% apply(.,2,median) %>% t()
    colnames(median)<-paste0(colnames(median),"_median")
    mean<-case_out_gamma_V_int_bar %>% apply(.,2,mean) %>% t()
    colnames(mean)<-paste0(colnames(mean),"_mean")
    sd<-case_out_gamma_V_int_bar %>% apply(.,2,sd) %>% t()
    colnames(sd)<-paste0(colnames(sd),"_sd")
    ci_025<-case_out_gamma_V_int_bar %>% apply(.,2,FUN = function(x){quantile(x,0.025)}) %>% t()
    ci_975<-case_out_gamma_V_int_bar %>% apply(.,2,FUN = function(x){quantile(x,0.975)}) %>% t()
    colnames(ci_025)<-paste0(colnames(ci_025),"_ci_025")
    colnames(ci_975)<-paste0(colnames(ci_975),"_ci_975")
    gamma_V_int_bar<-cbind(median,mean,sd,ci_025,ci_975)
    rownames(gamma_V_int_bar)<-paste0("gamma_V_int_bar",iteration)

    MAP<-list(emiss_mu_bar_list,emiss_var_bar_list,emiss_varmu_bar_list,dwell_mu_bar,dwell_var_bar,dwell_varmu_bar,gamma_prob_bar,gamma_int_bar,gamma_V_int_bar)
    names(MAP)<-c("emiss_mu_bar","emiss_var_bar","emiss_varmu_bar","dwell_mu_bar","dwell_var_bar","dwell_varmu_bar","gamma_prob_bar","gamma_int_bar","gamma_V_int_bar")
  }else{
    MAP<-list(emiss_mu_bar_list,emiss_var_bar_list,emiss_varmu_bar_list,dwell_mu_bar,dwell_var_bar,dwell_varmu_bar)
    names(MAP)<-c("emiss_mu_bar","emiss_var_bar","emiss_varmu_bar","dwell_mu_bar","dwell_var_bar","dwell_varmu_bar")
  }


  return(MAP)
}
#--------------------------------------------------------------------------------

################################################################################### START ##################################################################################
#-------------------------------------------------------------------------------
#parameter specification
#-------------------------------------------------------------------------------
#----Legend---------------------------------
# There are 3 hidden states scenarios m={1,2,3} (m)
# There are 4 senarios for observation number (o)
# * 1 - 200 obs
# * 2 - 500 obs
# * 3 - 1000 obs
# There are 3 scenarios for dwell time (d)
# * 1 - 1.4 time points (self-transition in mhmm is 0.5)
# * 2 - 3.5 time points (self-transition in mhmm is 0.75)
# * 3 - 19.5 time points (self-transition in mhmm is 0.95)
# * 4 - 99.5 time points (self-transition in mhmm is 0.99)
# J is representation of total number of iterations in MCMC chain (J)
# B is representation of burn_in period for iterations in MCMC chain (B)
# There are 100 simulated datasets from each scenario and each of the datasets
#  will to be fitted to the mhmm model (s)

# Creating vectors
m <- 2
o <- 1:3
d <- 2:4
s <-1:5

# Calling expand.grid() Function to create the parameter dataframe
pars = expand.grid(m=m,o=o,d=d,s=s)
#set global variables
J  <- 20#00
B  <- 10#00
n_dep=2


# Save data files
four_state_chosen_datasets<-readRDS("tests/data/convergence_medhmm_run_chosen_data_four_states.rds")
four_state_chosen_datasets<-four_state_chosen_datasets[[1]]

two_state_chosen_datasets<-readRDS("tests/data/convergence_medhmm_run_chosen_data_two_states.rds")
two_state_chosen_datasets<-two_state_chosen_datasets[[1]]

#==================================================================================================================================================================================================
#Outstart clean
#==================================================================================================================================================================================================
gelman_list_all<-list()

pars
parameter_row <- 1

d <- 1
d <- 4

out_list1 <-pbmcapply::pbmclapply(1:nrow(pars), function(parameter_row) {
  case_out<-NULL
  out1<-NULL
  out2<-NULL

  # Specify simulation parameters
  m  <- pars[parameter_row,1]
  o  <- pars[parameter_row,2]
  d  <- pars[parameter_row,3]
  s  <- pars[parameter_row,4]
  #==================================================================================================================================================================================================
  #load: simulated data, hyper-priors and starting values according (m)
  #==================================================================================================================================================================================================
  #define the expected medians
  exp_mean1=1/(-log(0.5))
  exp_mean2=1/(-log(0.75))
  exp_mean3=1/(-log(0.95))
  exp_mean4=1/(-log(0.99))
  #define standard deviations
  exp_sigma1=1.10093 #such that exp_sigma1=1/3*exp_mean1
  exp_sigma2=1.10093 #such that exp_sigma2=1/3*exp_mean2
  exp_sigma3=1.10093 #such that exp_sigma3=1/3*exp_mean3
  exp_sigma4=1.05902 #such that exp_sigma4=1/4*exp_mean4




  if(m==2){
    #simulated data
    sim_data<-two_state_chosen_datasets
    #---Emission distribution starting values
    emiss_distr <- matrix(c( 30,60*1.2,
                             60, 120*1.2), nrow = m, byrow = TRUE)
    emiss_distr2 <- matrix(c(60,(120+0.25*120)*1.2,
                             30, (60+0.25*60)*1.2), nrow = m, byrow = TRUE)
    #---Emission distribution starting values for second chain
    emiss_distr_2nd_chain <- matrix(c( 30*(base::sample(c(0.8,1.2),1)),60*1.2*(base::sample(c(0.8,1.2),1)),
                                       60*(base::sample(c(0.8,1.2),1)), 120*1.2*(base::sample(c(0.8,1.2),1))), nrow = m, byrow = TRUE)
    emiss_distr2_2nd_chain <- matrix(c(60*(base::sample(c(0.8,1.2),1)),(120+0.25*120)*1.2*(base::sample(c(0.8,1.2),1)),
                                       30*(base::sample(c(0.8,1.2),1)), (60+0.25*60)*1.2*(base::sample(c(0.8,1.2),1))), nrow = m, byrow = TRUE)
    #---Gamma starting values
    gamma <- matrix(c(0, 1,
                      1, 0), ncol = m, byrow = TRUE)
    gamma2 <- matrix(c(0, 1,
                       1, 0), ncol = m, byrow = TRUE)
    if(d==1){
      #-Dwell time starting values
      dwell_mu_var<-log(matrix(
        data = c(exp_mean1,exp_sigma1*1.2,
                 exp_mean1,exp_sigma1*1.2),
        byrow = T,
        ncol = 2,
        nrow = m))

      dwell_mu_var_2nd_chain<-matrix(
        data = c(log(exp_mean1),log(exp_sigma1)*1.1*(base::sample(c(0.9,1.1),1)),
                 log(exp_mean1),log(exp_sigma1)*1.1*(base::sample(c(0.9,1.1),1))),
        byrow = T,
        ncol = 2,
        nrow = m)
      #-Dwell time hyper priors
      dwell_hyp_prior<-list(d_mu0=log(c(exp_mean1,exp_mean1)),
                            s2_0=log(rep(10, m)),
                            alpha.sigma20	= rep(0.01, m),
                            beta.sigma20	= rep(0.01, m),
                            alpha.tau20		= rep(0.01, m),
                            beta.tau20		= rep(0.01, m))


    }else if(d==2){
      #-Dwell time starting values
      dwell_mu_var<-log(matrix(
        data = c(exp_mean2, exp_sigma2*1.2,
                 exp_mean2,  exp_sigma2*1.2),
        byrow = T,
        ncol = 2,
        nrow = m))

      dwell_mu_var_2nd_chain<-matrix(data = c(log(exp_mean2), log(exp_sigma2)*1.1*base::sample(c(0.9,1.1),1),
                                              log(exp_mean2),  log(exp_sigma2)*1.1*base::sample(c(0.9,1.1),1)), byrow = T,ncol = 2,nrow = m)
      #-Dwell time hyper priors
      dwell_hyp_prior<-list(d_mu0=log(c(exp_mean2,exp_mean2)),
                            s2_0=log(rep(10, m)),
                            alpha.sigma20	= rep(0.01, m),
                            beta.sigma20	= rep(0.01, m),
                            alpha.tau20		= rep(0.01, m),
                            beta.tau20		= rep(0.01, m))

    }else if(d==3){
      #-Dwell time starting values
      dwell_mu_var<- log(matrix(
        data = c(exp_mean3,exp_sigma3*1.2,
                 exp_mean3,exp_sigma3*1.2),
        byrow = T,
        ncol = 2,
        nrow = m))
      dwell_mu_var_2nd_chain<- matrix(
        data = c(log(exp_mean3),log(exp_sigma3)*1.1*(base::sample(c(0.9,1.1),1)),
                 log(exp_mean3),log(exp_sigma3)*1.1*(base::sample(c(0.9,1.1),1))),
        byrow = T,
        ncol = 2,
        nrow = m)
      #-Dwell time hyper priors
      dwell_hyp_prior<-list(d_mu0=log(c(exp_mean3,exp_mean3)),
                            s2_0=log(rep(10, m)),
                            alpha.sigma20	= rep(0.01, m),
                            beta.sigma20	= rep(0.01, m),
                            alpha.tau20		= rep(0.01, m),
                            beta.tau20		= rep(0.01, m))

    }else if(d==4){
      #-Dwell time starting values
      dwell_mu_var<- log(matrix(
        data = c(exp_mean4,exp_sigma4*1.2,
                 exp_mean4,exp_sigma4*1.2),
        byrow = T,
        ncol = 2,
        nrow = m))
      dwell_mu_var_2nd_chain<- matrix(
        data = c(log(exp_mean4),log(exp_sigma4)*1.1*(base::sample(c(0.9,1.1),1)),
                 log(exp_mean4),log(exp_sigma4)*1.1*(base::sample(c(0.9,1.1),1))),
        byrow = T,
        ncol = 2,
        nrow = m)
      #-Dwell time hyper priors
      dwell_hyp_prior<-list(d_mu0=log(c(exp_mean4,exp_mean4)),
                            s2_0=log(rep(10, m)),
                            alpha.sigma20	= rep(0.01, m),
                            beta.sigma20	= rep(0.01, m),
                            alpha.tau20		= rep(0.01, m),
                            beta.tau20		= rep(0.01, m))

    }

    #---Specify hyper-prior for the continuous emission distribution----------------
    hyp_pr <- list(
      emiss_mu0 = list(matrix(c(30,60), nrow = 1),
                       matrix(c(60,30), nrow = 1)),
      emiss_K0  = list(1,1),
      emiss_nu  = list(1,1),
      emiss_V   = list(rep(100, m),
                       rep(100, m)),
      emiss_a0  = list(rep(0.001,m),
                       rep(0.001,m)),
      emiss_b0  = list(rep(0.001,m),
                       rep(0.001,m)))

  }


  #Fit 2nd chain MHMM
  # Try statement, to catch errors and avoid interrupting the script
  ti <- Sys.time()

  if(d==1){
    max_dwell=10
  }else if(d==2){
    max_dwell=30
  }else if(d==3){
    max_dwell=100
  }else{
    max_dwell=600
  }
  input_data_mhmm<-sim_data[[o]][[d]][[s]]
  input_data_mhmm<-sim_data[[1]][[4]][[5]]
  out1<- try(medHMM_cont(s_data = input_data_mhmm$observations,
                         gen = list(m = m, n_dep = n_dep),
                         start_val = c(list(gamma), list(emiss_distr,emiss_distr2), list(dwell_mu_var)),
                         emiss_hyp_prior = hyp_pr,
                         dwell_hyp_prior = dwell_hyp_prior,
                         show_progress = TRUE,
                         mcmc = list(J = J, burn_in = B), return_path = TRUE, max_dwell = max_dwell),TRUE)
  print(.Last.value)
  options(show.error.messages = TRUE)
  out1[["time1"]] <- Sys.time() - ti
  if(!is.null(out1)){
    MAP1<-try(MAP_medhmm(case_out =out1,J=J,B=B,m=m,iteration=s),TRUE)
  }

  ti <- Sys.time()
  out2<- try(medHMM_cont(s_data = input_data_mhmm$observations,
                         gen = list(m = m, n_dep = n_dep),
                         start_val = c(list(gamma2), list(emiss_distr_2nd_chain,emiss_distr2_2nd_chain), list(dwell_mu_var_2nd_chain)),
                         emiss_hyp_prior = hyp_pr,
                         dwell_hyp_prior = dwell_hyp_prior,
                         show_progress = TRUE,
                         mcmc = list(J = J, burn_in = B), return_path = TRUE, max_dwell = max_dwell),TRUE)
  print(.Last.value)
  options(show.error.messages = TRUE)
  out2[["time1"]] <- Sys.time() - ti
  if(!is.null(out2)){
    MAP2<-try(MAP_medhmm(case_out =out2,J=J,B=B,m=m,iteration=s),TRUE)
  }
  # Calculate necessary statistics to be saved
  if(!is.null(out1) & !is.null(out2)){
    gelman_rubin_out<-try(ge_rub_two_models_emp(model1 = out1,model2 = out2,m = m,J = J,burn_in = B))

    # Define new output object
    case_out <- list(
      chain1=out1,
      MAP1=MAP1,
      chain2=out2,
      MAP2=MAP2,
      gelman_rubin_stats = gelman_rubin_out
    )
  }else{
    case_out <- list(
      chain1=out1,
      MAP1=MAP1,
      chain2=out2,
      MAP2=MAP2
    )
  }


  # Save output
  try(saveRDS(case_out, paste0("outputs/convergence_run",
                               "_m",m,
                               "_nt",o,
                               "_dwell",d,
                               "dt_set",s,
                               "_it",J,
                               "_burn_in",B,
                               ".rds")),TRUE)



  return(case_out)
}, mc.cores = cores, mc.set.seed = 42L)


# Calling expand.grid() Function to create the parameter dataframe
pars4 = expand.grid(m=4,o=o,d=d,s=s)

gelman_list_all4 <-pbmcapply::pbmclapply(1:nrow(pars4), function(parameter_row) {


  # Specify simulation parameters
  m  <- pars4[parameter_row,1]
  o  <- pars4[parameter_row,2]
  d  <- pars4[parameter_row,3]
  s  <- pars4[parameter_row,4]

  dataset<- readRDS(paste0("outputs/convergence_run",
                           "_m",m,
                           "_nt",o,
                           "_dwell",d,
                           "dt_set",s,
                           "_it",J,
                           "_burn_in",B,
                           ".rds"))
  out_gelman_row<-dataset$gelman_rubin_stats
  rownames(out_gelman_row)<- paste0("convergence_run",
                                    "_m",m,
                                    "_nt",o,
                                    "_dwell",d,
                                    "dt_set",s,
                                    "_it",J,
                                    "_burn_in",B)
  return(out_gelman_row)


}, mc.cores = cores, mc.set.seed = 42L)

gelman_list_all4<-do.call(rbind,gelman_list_all4)


saveRDS(gelman_list_all4,paste0("outputs/all_gelaman_medhmm_m4.rds"))



#Incorporate the trace plots output

library(ggplot2)
library(gridExtra)
library(grid)

traceplots_emiss_mu<-function(model1,model2,dep=1){
  model1_em<-as.data.frame(model1$emiss_mu_bar[[dep]])
  model2_em<-as.data.frame(model2$emiss_mu_bar[[dep]])
  names<-colnames(model2_em)
  data_points<-nrow(model2_em)
  p<-list()
  for(i in 1:ncol(model2_em)){
    p[[i]] <-ggplot(data=model1_em,aes_string(x=1:data_points, y=model1_em[,i]))+geom_line(color="red")+geom_line(data=model2_em, aes_string(x=1:data_points,y=model2_em[,i]),color="blue")+theme_minimal()+theme(plot.title = element_text(size=8,hjust = 0.5))+xlab("iteration")+ylab(colnames(model2_em)[i])+ggtitle(paste0("Emiss mu",dep,":dwell",d," :obs",o," :state",i," :set",s ))

  }
  out<-do.call(grid.arrange,c(p, list(top=textGrob(conv_name))))
  print(out)
}
traceplots_emiss_var<-function(model1,model2,dep=1){
  model1_em<-as.data.frame(model1$emiss_var_bar[[dep]])
  model2_em<-as.data.frame(model2$emiss_var_bar[[dep]])
  names<-colnames(model2_em)
  data_points<-nrow(model2_em)
  p<-list()
  for(i in 1:ncol(model2_em)){
    p[[i]]<-ggplot(data=model1_em,aes_string(x=1:data_points, y=model1_em[,i]))+geom_line(color="red")+geom_line(data=model2_em, aes_string(x=1:data_points,y=model2_em[,i]),color="blue")+theme_minimal()+theme(plot.title = element_text(size=8,hjust = 0.5))+xlab("iteration")+ylab(colnames(model2_em)[i])+ggtitle(paste0("Emiss var",dep,":dwell",d," :obs",o," :state",i," :set",s ))

  }
  out<-do.call(grid.arrange,c(p, list(top=textGrob(conv_name))))
  print(out)
}

traceplots_dwell_mu<-function(model1,model2){
  model1_em<-as.data.frame(model1$dwell_mu_bar)
  model2_em<-as.data.frame(model2$dwell_mu_bar)
  names<-colnames(model2_em)
  data_points<-nrow(model2_em)
  p<-list()
  for(i in 1:ncol(model2_em)){
    p[[i]]<-ggplot(data=model1_em,aes_string(x=1:data_points, y=model1_em[,i]))+geom_line(color="red")+geom_line(data=model2_em, aes_string(x=1:data_points,y=model2_em[,i]),color="blue")+theme_minimal()+theme(plot.title = element_text(size=8,hjust = 0.5))+xlab("iteration")+ylab(colnames(model2_em)[i])+ggtitle(paste0("Dwell mu: dwell",d," :obs",o," :state",i," :set",s ))

  }
  out<-do.call(grid.arrange,c(p, list(top=textGrob(conv_name))))
  print(out)}
traceplots_dwell_var<-function(model1,model2){
  model1_em<-as.data.frame(model1$dwell_var_bar)
  model2_em<-as.data.frame(model2$dwell_var_bar)
  names<-colnames(model2_em)
  data_points<-nrow(model2_em)
  p<-list()
  for(i in 1:ncol(model2_em)){
    p[[i]]<-ggplot(data=model1_em,aes_string(x=1:data_points, y=model1_em[,i]))+geom_line(color="red")+geom_line(data=model2_em, aes_string(x=1:data_points,y=model2_em[,i]),color="blue")+theme_minimal()+theme(plot.title = element_text(size=8,hjust = 0.5))+xlab("iteration")+ylab(colnames(model2_em)[i])+ggtitle(paste0("Dwell var: dwell",d," :obs",o," :state",i," :set",s ))

  }
  out<-do.call(grid.arrange,c(p, list(top=textGrob(conv_name))))
  print(out)
}

traceplots_gamma<-function(model1,model2){
  model1_em<-as.data.frame(model1$gamma_prob_bar)
  model2_em<-as.data.frame(model2$gamma_prob_bar)
  names<-colnames(model2_em)
  data_points<-nrow(model2_em)
  p<-list()
  for(i in 1:ncol(model2_em)){
    p[[i]]<-ggplot(data=model1_em,aes_string(x=1:data_points, y=model1_em[,i]))+geom_line(color="red")+geom_line(data=model2_em, aes_string(x=1:data_points,y=model2_em[,i]),color="blue")+theme_minimal()+theme(plot.title = element_text(size=8,hjust = 0.5))+xlab("iteration")+ylab(colnames(model2_em)[i])+ggtitle(paste0("gamma: dwell",d," :obs",o," :state",i," :set",s))
  }
  out<-do.call(grid.arrange, c(p, list(top=textGrob(conv_name))))
  print(out)
}


# Calling expand.grid() Function to create the parameter dataframe
pars = expand.grid(m=m,o=o,d=d,s=s)

pdf("outputs/trace_plots_medhmm.pdf")
for(parameter in 1:nrow(pars)){
  m  <- pars[parameter,1]
  o  <- pars[parameter,2]
  d  <- pars[parameter,3]
  s  <- pars[parameter,4]


  dataset<- readRDS(paste0("outputs/convergence_run",
                           "_m",m,
                           "_nt",o,
                           "_dwell",d,
                           "dt_set",s,
                           "_it",J,
                           "_burn_in",B,
                           ".rds"))
  conv_name<-paste0("medhmm_convergence_run","_m",m,
                    "_nt",o,
                    "_dwell",d,
                    "dt_set",s,
                    "_it",J,
                    "_burn_in",B)
  model1=try(dataset$chain1)
  model2=try(dataset$chain2)
  if(is.list(model2) && is.list(model1)){
    traceplots_emiss_mu(model1 = model1,model2 = model2,dep=1)
    traceplots_emiss_var(model1 = model1,model2 = model2,dep=1)
    traceplots_emiss_mu(model1 = model1,model2 = model2,dep=2)
    traceplots_emiss_var(model1 = model1,model2 = model2,dep=2)
    traceplots_dwell_mu(model1 = model1,model2 = model2)
    traceplots_dwell_var(model1 = model1,model2 = model2)
    traceplots_gamma(model1 = model1,model2 = model2)
  }

}
dev.off()

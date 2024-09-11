#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

library(tidyverse)
library(mHMMbayes)
library(medHMM)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# Define utility functions:

local_decoding <- function(object, burn_in = NULL){

    # Set up
    n_iter <- object$input$J
    if(is.null(burn_in)){
        burn_in <- object$input$burn_in
    }
    m <- object$input$m
    n_subj <- object$input$n_subj
    n_vary <- object$input$n_vary

    # Find local decoding
    probs <- do.call(rbind, lapply(
        1:n_subj, function(s) {
            cbind(s,
                  t(apply(object$sample_path[[s]][,burn_in:n_iter],
                          1,
                          function(e) {c(which.max(table(factor(e,levels = 1:m))),
                                         table(factor(e,levels = 1:m))/sum(table(e) ) )} ) ),
                  1:n_vary[s])
        }
    ))
    probs <- as.data.frame(probs)
    names(probs) <- c("subj","state",paste0("pr_state",1:m),"occasion")

    return(probs)

}

cross_entropy_loss <- function(y_true, y_pred) {
    # Ensure predictions are probabilities
    epsilon <- 1e-15
    y_pred <- pmin(pmax(y_pred, epsilon, na.rm = TRUE), 1 - epsilon)

    # Calculate cross-entropy loss
    loss <- -sum(y_true * log(y_pred), na.rm = TRUE) / nrow(y_true)
    # loss <- -sum(y_true * log(y_pred))

    return(loss)
}

#==============================================================================#
# Repeat the same outside of the function to play around more easily:

# Parameters
n_t = 250
n = 20
m = 3
n_dep = 2
gamma_medhmm = matrix(c(0, 0.7, 0.3,
                        0.5, 0, 0.5,
                        0.6, 0.4, 0), nrow = m, ncol = m, byrow = TRUE)
gamma_mhmm = matrix(c(0.95, 0.035, 0.015,
                      0.5, 0.9, 0.5,
                      0.006, 0.004, 0.99), nrow = m, ncol = m, byrow = TRUE)
emiss_distr = list(matrix(c(10,2,
                            50,2,
                            2,2), nrow = m, ncol = 2, byrow = TRUE),
                   matrix(c(5,2,
                            20,2,
                            5,2), nrow = m, ncol = 2, byrow = TRUE))
dwell_distr = matrix(log(c(10,
                           5,
                           30)), nrow = m, ncol = 1, byrow = TRUE)
var_gamma = 0.1
var_emiss = c(10,5)
var_dwell = 0.01
hyp_prior_emiss =  prior_emiss_cont(
    gen = list(m = m, n_dep = n_dep),
    emiss_mu0 = list(matrix(c(10, 50, 2), nrow = 1),
                     matrix(c(5, 20, 5), nrow = 1)),
    emiss_K0 = list(1, 1),
    emiss_V =  list(rep(10, m), rep(25, m)),
    emiss_nu = list(1, 1),
    emiss_a0 = list(rep(0.01, m), rep(0.01, m)),
    emiss_b0 = list(rep(0.01, m), rep(0.01, m)))

hyp_prior_dwell = list(
    dwell_mu0 = matrix(log(c(10,5,30)), nrow = 1, ncol = 3), # nrow = number of covariates + 1; ncol = number of hidden states
    dwell_K0  = c(0.01),
    dwell_nu  = c(1),
    dwell_V   = rep(0.1, m)
)

fit_fraction = 0.8
forecast_steps = 50
n_iter = 500
burn_in = 250
Mx = 50


# Run
set.seed(42)
data_cont <- medHMM::sim_medHMM(n_t, n, data_distr = 'continuous', m, n_dep = n_dep,
                                dwell_distr = dwell_distr, dwell_type = 'poisson',
                                gamma = gamma_medhmm, emiss_distr = emiss_distr,
                                var_gamma = var_gamma, var_emiss = var_emiss, var_dwell = var_dwell, return_ind_par = TRUE)

# Extract observations and states
obs <- as.data.frame(data_cont$obs)
true_states <- as.data.frame(data_cont$states)

# Ensure columns are properly named
colnames(obs) <- c("subj", "obs1", "obs2")
colnames(true_states) <- c("subj", "state")

# Add a time column to the observations
obs <- obs %>%
    group_by(subj) %>%
    mutate(time = row_number()) %>%
    ungroup()

# Step 2: Split data into training and forecasting sets
split_data <- function(data, fit_fraction) {
    unique_subj <- unique(data$subj)
    split_data_list <- lapply(unique_subj, function(subj) {
        subj_data <- data %>% filter(subj == !!subj)
        n_fit <- floor(nrow(subj_data) * fit_fraction)
        list(
            fit = subj_data[1:n_fit, ],
            forecast = subj_data[(n_fit + 1):nrow(subj_data), ]
        )
    })
    list(
        fit = bind_rows(lapply(split_data_list, `[[`, "fit")),
        forecast = bind_rows(lapply(split_data_list, `[[`, "forecast"))
    )
}

split_data_abs <- function(data, n_train) {
    unique_subj <- unique(data$subj)
    split_data_list <- lapply(unique_subj, function(subj) {
        subj_data <- data %>% filter(subj == !!subj)
        n_fit <- floor(n_train)
        list(
            fit = subj_data[1:n_fit, ],
            forecast = subj_data[(n_fit + 1):nrow(subj_data), ]
        )
    })
    list(
        fit = bind_rows(lapply(split_data_list, `[[`, "fit")),
        forecast = bind_rows(lapply(split_data_list, `[[`, "forecast"))
    )
}

split_obs <- split_data(obs, fit_fraction)
fit_data <- split_obs$fit
forecast_data <- split_obs$forecast

# Add a time column to the true states
train_states <- split_data(true_states, fit_fraction)[[1]] %>%
    mutate(state = as.integer(state)) %>%
    group_by(subj) %>%
    mutate(time = row_number()) %>%
    ungroup()

forecast_states <- split_data(true_states, fit_fraction)[[2]] %>%
    mutate(state = as.integer(state)) %>%
    group_by(subj) %>%
    mutate(time = row_number()) %>%
    ungroup()


# Step 3: Fit the HMM model
out_medhmm <- medHMM_cont_shiftpois(s_data = as.matrix(fit_data[,-4]),
                                    gen = list(m = m, n_dep = n_dep),
                                    start_val = c(list(gamma_medhmm), emiss_distr, list(exp(dwell_distr))),
                                    emiss_hyp_prior = hyp_prior_emiss,
                                    dwell_hyp_prior = hyp_prior_dwell,
                                    shift = NULL,
                                    show_progress = TRUE,
                                    mcmc = list(J = n_iter, burn_in = burn_in),
                                    return_path = TRUE,
                                    max_dwell = Mx)

out_mhmm <- mHMM(s_data = as.matrix(fit_data[,-4]),
                 data_distr = 'continuous',
                 gen = list(m = m, n_dep = n_dep),
                 start_val = c(list(gamma_mhmm), emiss_distr),
                 emiss_hyp_prior = hyp_prior_emiss,
                 show_progress = TRUE,
                 return_path = TRUE,
                 mcmc = list(J = n_iter, burn_in = burn_in))

#==============================================================================#

# Step 4a: Forecast using the fitted model
# forecast_results <- forecast_mHMM1(object = out, s_data = as.matrix(forecast_data[,-4]), forecast_steps = forecast_steps)
forecast_results_medhmm <- forecast_medHMM1(object = out_medhmm, s_data = as.matrix(fit_data[,-4]), forecast_steps = forecast_steps, Mx = Mx, return_all = TRUE)
forecast_results_mhmm <- forecast_mHMM1(object = out_mhmm, s_data = as.matrix(fit_data[,-4]), forecast_steps = forecast_steps, return_all = TRUE)

# Check decoding (sanity check):
decoding_medhmm <- local_decoding(out_medhmm)
decoding_mhmm <- local_decoding(out_mhmm)

mean(decoding_medhmm$state == train_states$state)
mean(decoding_mhmm$state == train_states$state)


#==============================================================================#

# Cross-entropy: mHMM does better (medHMM over confidetend?)
for (h in c(1, seq(5,50,5))){
    true_h_states <- forecast_states %>%
        as.data.frame() %>%
        filter(time <= h)

    pred_h_states_medhmm <- forecast_results_medhmm %>%
        as.data.frame() %>%
        filter(horizon <= h & horizon > 0)

    pred_h_states_mhmm <- forecast_results_mhmm %>%
        as.data.frame() %>%
        filter(horizon <= h & horizon > 0)

    one_hot_encoded <- matrix(0, nrow = length(true_h_states$state), ncol = m)
    one_hot_encoded[cbind(1:length(true_h_states$state), true_h_states$state)] <- 1

    cat("Average cross-entropy for horizon:",h,"\n",
        "   medHMM:",cross_entropy_loss(y_true = one_hot_encoded, y_pred = pred_h_states_medhmm[,4:6]),"\n",
        "   mHMM:  ",cross_entropy_loss(y_true = one_hot_encoded, y_pred = pred_h_states_mhmm[,4:6]),"\n")
}

# PR AUC: mHMM does better (medHMM over confidetend?)
for (h in c(1, seq(5,50,5))){
    true_h_states <- forecast_states %>%
        as.data.frame() %>%
        filter(time <= h)

    pred_h_states_medhmm <- forecast_results_medhmm %>%
        as.data.frame() %>%
        filter(horizon <= h & horizon > 0)

    pred_h_states_mhmm <- forecast_results_mhmm %>%
        as.data.frame() %>%
        filter(horizon <= h & horizon > 0)

    one_hot_encoded <- matrix(0, nrow = length(true_h_states$state), ncol = m)
    one_hot_encoded[cbind(1:length(true_h_states$state), true_h_states$state)] <- 1

    cat("Average cross-entropy for horizon:",h,"\n",
        "   medHMM:",yardstick::pr_auc_vec(truth = factor(true_h_states$state), estimate = as.matrix(pred_h_states_medhmm[,4:6])),"\n",
        "   mHMM:  ",yardstick::pr_auc_vec(truth = factor(true_h_states$state), estimate = as.matrix(pred_h_states_mhmm[,4:6])),"\n")
}

# Accuracy: medHMM does better
for (h in c(1,seq(5,50,5))){
    true_h_states <- forecast_states %>%
        as.data.frame() %>%
        filter(time <= h)

    pred_h_states_medhmm <- forecast_results_medhmm %>%
        as.data.frame() %>%
        filter(horizon <= h & horizon > 0)

    pred_h_states_mhmm <- forecast_results_mhmm %>%
        as.data.frame() %>%
        filter(horizon <= h & horizon > 0)

    one_hot_encoded <- matrix(0, nrow = length(true_h_states$state), ncol = m)
    one_hot_encoded[cbind(1:length(true_h_states$state), true_h_states$state)] <- 1

    cat("Accuracy for horizon:",h,"\n",
        "   medHMM:",mean(true_h_states$state == pred_h_states_medhmm$state),"\n",
        "   mHMM:  ",mean(true_h_states$state == pred_h_states_mhmm$state),"\n")
}

# F-one score: medHMM does better
for (h in c(1,seq(5,50,5))){
    true_h_states <- forecast_states %>%
        as.data.frame() %>%
        filter(time <= h)

    pred_h_states_medhmm <- forecast_results_medhmm %>%
        as.data.frame() %>%
        filter(horizon <= h & horizon > 0)

    pred_h_states_mhmm <- forecast_results_mhmm %>%
        as.data.frame() %>%
        filter(horizon <= h & horizon > 0)

    one_hot_encoded <- matrix(0, nrow = length(true_h_states$state), ncol = m)
    one_hot_encoded[cbind(1:length(true_h_states$state), true_h_states$state)] <- 1

    cat("Accuracy for horizon:",h,"\n",
        "   medHMM:",yardstick::f_meas_vec(factor(true_h_states$state), factor(pred_h_states_medhmm$state)),"\n",
        "   mHMM:  ",yardstick::f_meas_vec(factor(true_h_states$state), factor(pred_h_states_mhmm$state)),"\n")
}


# In data frame format:

# Initialize an empty data frame to store the results
results_df <- data.frame()

# Define the horizons
# horizons <- c(1, seq(5, 50, 5))

horizons <- seq(1, 50, 1)

# Loop through each horizon and calculate metrics
for (h in horizons) {
    # Filter true and predicted states based on the horizon
    true_h_states <- forecast_states %>%
        as.data.frame() %>%
        filter(time <= h)

    pred_h_states_medhmm <- forecast_results_medhmm %>%
        as.data.frame() %>%
        filter(horizon <= h & horizon > 0)

    pred_h_states_mhmm <- forecast_results_mhmm %>%
        as.data.frame() %>%
        filter(horizon <= h & horizon > 0)

    # One-hot encode true states
    one_hot_encoded <- matrix(0, nrow = length(true_h_states$state), ncol = m)
    one_hot_encoded[cbind(1:length(true_h_states$state), true_h_states$state)] <- 1

    # Calculate metrics
    cross_entropy_medhmm <- cross_entropy_loss(y_true = one_hot_encoded, y_pred = pred_h_states_medhmm[, 4:6])
    cross_entropy_mhmm <- cross_entropy_loss(y_true = one_hot_encoded, y_pred = pred_h_states_mhmm[, 4:6])

    pr_auc_medhmm <- yardstick::pr_auc_vec(truth = factor(true_h_states$state), estimate = as.matrix(pred_h_states_medhmm[, 4:6]))
    pr_auc_mhmm <- yardstick::pr_auc_vec(truth = factor(true_h_states$state), estimate = as.matrix(pred_h_states_mhmm[, 4:6]))

    accuracy_medhmm <- mean(true_h_states$state == pred_h_states_medhmm$state)
    accuracy_mhmm <- mean(true_h_states$state == pred_h_states_mhmm$state)

    f1_score_medhmm <- yardstick::f_meas_vec(factor(true_h_states$state), factor(pred_h_states_medhmm$state, levels = 1:m))
    f1_score_mhmm <- yardstick::f_meas_vec(factor(true_h_states$state), factor(pred_h_states_mhmm$state, levels = 1:m))

    # Append results to the data frame
    results_df <- rbind(results_df, data.frame(
        Horizon = h,
        CrossEntropy_medHMM = cross_entropy_medhmm,
        CrossEntropy_mHMM = cross_entropy_mhmm,
        PRAUC_medHMM = pr_auc_medhmm,
        PRAUC_mHMM = pr_auc_mhmm,
        Accuracy_medHMM = accuracy_medhmm,
        Accuracy_mHMM = accuracy_mhmm,
        F1Score_medHMM = f1_score_medhmm,
        F1Score_mHMM = f1_score_mhmm
    ))
}

# Return the results data frame
results_df

#==============================================================================#

# Check cases:
#   - near perfect decoding both: 5% emiss overlap, 5% dwell overlap, 25 ind, 500 occ, 2 n_dep, 3 dwell times
#   - medHMM better: 50% emiss overlap, 5% dwell overlap, 25 ind, 500 occ, 1 n_dep, 3 dwell times

# Scenarios A: near perfect decoding
#   dur_s = 1, 5, 9
#   dep_s = 6
#   occ_s = 1
#   ind_s = 3

# Scenarios B: medHMM with better decoding
#   dur_s = 1, 5, 9
#   dep_s = 1
#   occ_s = 1
#   ind_s = 3

dur_s = 1
dep_s = 6
occ_s = 1
ind_s = 3

# Load design
design <- read_csv("/Users/a6159737/Documents/Utrecht University/PhD/Projects/Simulation studies/medhmm-sim/seb/main_sim_fat/design.csv")


# Define simulation arguments
cores = 1
J = 500
burn_in = 250
m = 3


# #number of set of mean duration matrices
# dur_s <- crit$dur_s[j]
# #number of scenarios of the number of observations per individual
# ind_s <- crit$ind_s[j]

n_subj=n_distr[[ind_s]]

# #number of scenarios of the number of observations per individual
# occ_s <-crit$occ_s[j]

n_t=n_t_distr[[occ_s]]
n_t = 625

# # #number of dependent variable scenarios
# # dep_s<-crit$dep_s[j]
# #number of samples within the same scenario
# sampl <- 1
# #repetition number
# rep_n <- crit$rep_n[j]

if(dep_s<=3 | dep_s == 7){
    n_dep<-1
}else{
    n_dep=2
}

#=============== MHMM + data simulation ==========================================

#dwell_distr1 is expected dwell time and I calculate the start gamma matrix out of it
dwell_distr1=matrix(unlist(dwell_distr[[dur_s]]),ncol=1)
gam_start<-data.frame(expected=dwell_distr1) %>% reframe(gamma=exp(-1/expected)-0.1) %>% reframe(gamma=gamma,rest=(1-gamma)/2)

# new
gamma_start <- matrix(c(0, 0.7, 0.3,
                        0.5, 0, 0.5,
                        0.6, 0.4, 0), nrow = m, ncol = m, byrow = TRUE)

rest <- gam_start$rest*2

gamma_start <- matrix(rep(rest,3), nrow = 3) * gamma_start
diag(gamma_start) <- gam_start$gamma


if(n_dep==1){
    emiss=matrix(unlist(dep_distr[[dep_s]]), nrow = m, byrow = FALSE)
    emiss_start=cbind(emiss[,1],emiss[,2]*1.5)
    emiss<-list(emiss)

}else{
    emiss=lapply(dep_distr[[dep_s]],function(x){matrix(unlist(x), nrow = m, byrow = FALSE)}  )
    emiss_start=lapply(emiss, function(x){x[,2]<-x[,2]*1.5
    cbind(x[,1],x[,2])})
}

if(n_dep == 1){
    emiss_hyp_pr <- prior_emiss_cont(
        gen = list(m = m, n_dep = n_dep),
        emiss_mu0 = list(matrix(unlist(dep_distr[[dep_s]][,1]), nrow = 1)),
        emiss_K0  = list(1),
        emiss_nu  = list(1),
        emiss_V   = list(rep(10, m)),
        emiss_a0  = list(rep(0.01, m)),
        emiss_b0  = list(rep(0.01, m))
    )
} else {
    emiss_hyp_pr <- prior_emiss_cont(
        gen = list(m = m, n_dep = n_dep),
        emiss_mu0 = list(matrix(unlist(dep_distr[[dep_s]][[1]][,1]), nrow = 1),
                         matrix(unlist(dep_distr[[dep_s]][[2]][,1]), nrow = 1)),
        emiss_K0  = list(1, 1),
        emiss_nu  = list(1, 1),
        emiss_V   = list(rep(10, m), rep(10, m)),
        emiss_a0  = list(rep(0.01, m), rep(0.01, m)),
        emiss_b0  = list(rep(0.01, m), rep(0.01, m))
    )
}

gamma <-matrix(c(0, 0.7, 0.3,
                 0.5, 0, 0.5,
                 0.6, 0.4, 0), nrow = m, ncol = m, byrow = TRUE)

dwell_ss_var<-dwell_var_ss_s[[dur_s]]

emiss_ss_var<-rep(10, length(emiss))

# SIMULATE DATA
sim_data <- sim_medHMM(n_t, n_subj, data_distr = 'continuous', m, n_dep = n_dep,
                       dwell_distr = log(dwell_distr1), dwell_type = 'poisson', shift = 0,
                       start_state = NULL, q_emiss = NULL, gamma = gamma, emiss_distr = emiss, xx_vec = NULL, beta = NULL,
                       var_gamma = gamma_ss_var, var_emiss = emiss_ss_var, var_dwell = as.numeric(dwell_ss_var), return_ind_par = TRUE)

assign(paste0("sim_data_","dur_s_",dur_s,"ind_s_",ind_s,"occ_s_",occ_s,"dep_s_",dep_s), sim_data)

#check simulated data
colnames(sim_data$states)<-c("subject","state")
summary_dwell_ss<-sim_data$states %>%as.data.frame() %>%  group_by(subject)%>%
    reframe(length=rle(state)[[1]],state=rle(state)[[2]])%>%
    group_by(subject,state) %>%reframe(mean_emp_dwell=mean(length),median_emp_dwell=median(length))
summary_dwell_ss$subject<-as.factor(summary_dwell_ss$subject)
summary_dwell_ss$state<-as.factor(summary_dwell_ss$state)
summary_dwell_ss=as.data.frame(summary_dwell_ss)

ob<-data.frame(sim_data[['obs']])
colnames(ob)<-c('subject',paste0('dep',1:n_dep))
dw<-data.frame(sim_data[['states']])
colnames(dw)<-c('subject','state')
data=cbind(ob,state=dw$state) %>% mutate(., state=as.factor(state), subject=as.factor(subject))

# New:

# Extract observations and states
obs <- as.data.frame(sim_data$obs)
true_states <- as.data.frame(sim_data$states)

# Ensure columns are properly named
colnames(obs) <- c("subj", paste0("obs", 1:n_dep))
colnames(true_states) <- c("subj", "state")

# Add a time column to the observations
obs <- obs %>%
    group_by(subj) %>%
    mutate(time = row_number()) %>%
    ungroup()

# Step 2: Split data into training and forecasting sets
split_data <- function(data, fit_fraction) {
    unique_subj <- unique(data$subj)
    split_data_list <- lapply(unique_subj, function(subj) {
        subj_data <- data %>% filter(subj == !!subj)
        n_fit <- floor(nrow(subj_data) * fit_fraction)
        list(
            fit = subj_data[1:n_fit, ],
            forecast = subj_data[(n_fit + 1):nrow(subj_data), ]
        )
    })
    list(
        fit = bind_rows(lapply(split_data_list, `[[`, "fit")),
        forecast = bind_rows(lapply(split_data_list, `[[`, "forecast"))
    )
}

split_obs <- split_data(obs, fit_fraction = 0.8)
fit_data <- split_obs$fit
forecast_data <- split_obs$forecast

# Add a time column to the true states
train_states <- split_data(true_states, fit_fraction = 0.8)[[1]] %>%
    mutate(state = as.integer(state)) %>%
    group_by(subj) %>%
    mutate(time = row_number()) %>%
    ungroup()

forecast_states <- split_data(true_states, fit_fraction = 0.8)[[2]] %>%
    mutate(state = as.integer(state)) %>%
    group_by(subj) %>%
    mutate(time = row_number()) %>%
    ungroup()

# MHMM
mhmm_case_out <- NULL

try({

    mHMM_cont <- NULL

    start_time1 <- Sys.time()


    mHMM_cont <- mHMMbayes::mHMM(s_data = as.matrix(fit_data[,-(n_dep+2)]),
                                 data_distr = 'continuous',
                                 gen = list(m = m, n_dep = n_dep),
                                 start_val = c(list(gamma_start), emiss),
                                 emiss_hyp_prior = emiss_hyp_pr,
                                 show_progress = TRUE,
                                 mcmc = list(J = J, burn_in = burn_in),
                                 return_path = TRUE)

    end_time1 <- Sys.time()

    if(!is.null(mHMM_cont)){
        stored_map_mhmm<-try(MAP_mHMM(case_out = mHMM_cont, iteration = rep_n, J=J, B = burn_in, m = m))
        state_decoding<-try(local_decoding(out=mHMM_cont))
        true_st<-as.data.frame(sim_data$states)
        mhmm_state_decoding_table<-try(data.frame(subject=train_states[,1],true_state_decoding=train_states[,2],vit_state_decoding=state_decoding$state))
    }

    mhmm_case_out<-list(MAP=stored_map_mhmm,state_decoding=mhmm_state_decoding_table,execution_time=end_time1-start_time1)
    # rm(mHMM_cont)

    # saveRDS(mhmm_case_out, paste0("/Users/a6159737/Documents/Utrecht University/PhD/Projects/Simulation studies/medhmm-sim/seb/res/","mhmm_res_","dur_s_",dur_s,"ind_s_",ind_s,"occ_s_",occ_s,"dep_s_",dep_s,"_iter",J,".rds"))

})

#============ MEDHMM ==============================================================

if(dur_s<=4){
    max_dwell=60 # qpois(0.999, 26*1.5)
}else if(dur_s<=8){
    max_dwell=112 # qpois(0.999, 55*1.5)
}else{
    max_dwell=241 # qpois(0.999, 131*1.5)
}

dwell_distr1=matrix(unlist(dwell_distr[[dur_s]]),ncol=1)

dwell_hyp_pr <- list(
    dwell_mu0 = matrix(log(as.numeric(dwell_distr[[dur_s]])), nrow = 1, ncol = 3),
    dwell_K0  = c(1),
    dwell_nu  = c(1),
    dwell_V   = rep(0.1, m)
)
if(n_dep==1){
    emiss=matrix(unlist(dep_distr[[dep_s]]), nrow = m, byrow = FALSE)
    emiss_start=cbind(emiss[,1],emiss[,2]*1.5)
    emiss<-list(emiss)

}else{
    emiss=lapply(dep_distr[[dep_s]],function(x){matrix(unlist(x), nrow = m, byrow = FALSE)}  )
    emiss_start=lapply(emiss, function(x){x[,2]<-x[,2]*1.5
    cbind(x[,1],x[,2])})
}

if(n_dep == 1){
    emiss_hyp_pr <- list(
        emiss_mu0 = list(matrix(unlist(dep_distr[[dep_s]][,1]), nrow = 1)),
        emiss_K0  = list(1),
        emiss_nu  = list(1),
        emiss_V   = list(rep(10, m)),
        emiss_a0  = list(rep(0.01, m)),
        emiss_b0  = list(rep(0.01, m))
    )
}else{
    emiss_hyp_pr <- list(
        emiss_mu0 = list(matrix(unlist(dep_distr[[dep_s]][[1]][,1]), nrow = 1),
                         matrix(unlist(dep_distr[[dep_s]][[2]][,1]), nrow = 1)),
        emiss_K0  = list(1, 1),
        emiss_nu  = list(1, 1),
        emiss_V   = list(rep(10, m), rep(10, m)),
        emiss_a0  = list(rep(0.01, m), rep(0.01, m)),
        emiss_b0  = list(rep(0.01, m), rep(0.01, m))
    )
}

medhmm_case_out <- NULL

try({

    medHMM_cont_shiftpois <- NULL

    start_time2 <- Sys.time()
    medHMM_cont_shiftpois <- try(medHMM_cont_shiftpois(s_data = as.matrix(fit_data[,-(n_dep+2)]),
                                                       shift = 1,
                                                       gen = list(m = m, n_dep = n_dep),
                                                       start_val = c(list(gamma), emiss, list(dwell_distr1)),
                                                       emiss_hyp_prior = emiss_hyp_pr,
                                                       dwell_hyp_prior = dwell_hyp_pr,
                                                       show_progress = TRUE,
                                                       mcmc = list(J = J, burn_in = burn_in),
                                                       return_path = TRUE,
                                                       max_dwell = max_dwell))

    end_time2 <- Sys.time()


    if(!is.null(medHMM_cont_shiftpois)){
        stored_map_medhmm<-try(MAP_medhmm(case_out = medHMM_cont_shiftpois, iteration = rep_n, J=J, B = burn_in, m = m))
        state_decoding<-try(local_decoding(out=medHMM_cont_shiftpois))
        true_st<-as.data.frame(sim_data$states)
        medhmm_state_decoding_table<-try(data.frame(subject=train_states[,1],true_state_decoding=train_states[,2],vit_state_decoding=state_decoding$state))
    }

    medhmm_case_out<-list(MAP=stored_map_medhmm,state_decoding=medhmm_state_decoding_table,execution_time=end_time2-start_time2)

    # rm(medHMM_cont_shiftpois)

    # saveRDS(medhmm_case_out, paste0("/Users/a6159737/Documents/Utrecht University/PhD/Projects/Simulation studies/medhmm-sim/seb/res/","medhmm_res_","dur_s_",dur_s,"ind_s_",ind_s,"occ_s_",occ_s,"dep_s_",dep_s,"_iter",J,".rds"))

})

# Create list of results
mhmm_and_medhhmm_case_out<-list(MHMM=mhmm_case_out, MEDHMM=medhmm_case_out,sim_data=sim_data)

# Save intermediate results
# saveRDS(mhmm_and_medhhmm_case_out, paste0("/gpfs/home3/moragas/main_sim_fat/outputs/job",job,"/","mhmm_and_medhhmm_res_","dur_s_",dur_s,"ind_s_",ind_s,"occ_s_",occ_s,"dep_s_",dep_s,"_iter",J,"_rep",rep_n,".rds"))

# return(mhmm_and_medhhmm_case_out)

forecast_steps = 50
# Step 4a: Forecast using the fitted model
# forecast_results <- forecast_mHMM1(object = out, s_data = as.matrix(forecast_data[,-4]), forecast_steps = forecast_steps)
forecast_results_medhmm <- forecast_medHMM1(object = medHMM_cont_shiftpois, s_data = as.matrix(fit_data[,-(n_dep+2)]), forecast_steps = forecast_steps, Mx = Mx, return_all = TRUE)
forecast_results_mhmm <- forecast_mHMM1(object = mHMM_cont, s_data = as.matrix(fit_data[,-(n_dep+2)]), forecast_steps = forecast_steps, return_all = TRUE)





# Initialize an empty data frame to store the results
results_df <- data.frame()

# Define the horizons
# horizons <- c(1, seq(5, 50, 5))

horizons <- seq(1, 50, 1)

# Loop through each horizon and calculate metrics
for (h in horizons) {

    # Filter true and predicted states based on the horizon
    true_h_states <- forecast_states %>%
        as.data.frame() %>%
        filter(time <= h)

    pred_h_states_medhmm <- forecast_results_medhmm %>%
        as.data.frame() %>%
        filter(horizon <= h & horizon > 0)

    pred_h_states_mhmm <- forecast_results_mhmm %>%
        as.data.frame() %>%
        filter(horizon <= h & horizon > 0)

    # One-hot encode true states
    one_hot_encoded <- matrix(0, nrow = length(true_h_states$state), ncol = m)
    one_hot_encoded[cbind(1:length(true_h_states$state), true_h_states$state)] <- 1

    # Calculate metrics
    cross_entropy_medhmm <- cross_entropy_loss(y_true = one_hot_encoded, y_pred = pred_h_states_medhmm[, 4:6])
    cross_entropy_mhmm <- cross_entropy_loss(y_true = one_hot_encoded, y_pred = pred_h_states_mhmm[, 4:6])

    pr_auc_medhmm <- yardstick::pr_auc_vec(truth = factor(true_h_states$state, levels = 1:m), estimate = as.matrix(pred_h_states_medhmm[, 4:6]))
    pr_auc_mhmm <- yardstick::pr_auc_vec(truth = factor(true_h_states$state, levels = 1:m), estimate = as.matrix(pred_h_states_mhmm[, 4:6]))

    accuracy_medhmm <- mean(true_h_states$state == pred_h_states_medhmm$state)
    accuracy_mhmm <- mean(true_h_states$state == pred_h_states_mhmm$state)

    f1_score_medhmm <- yardstick::f_meas_vec(factor(true_h_states$state, levels = 1:m), factor(pred_h_states_medhmm$state, levels = 1:m))
    f1_score_mhmm <- yardstick::f_meas_vec(factor(true_h_states$state, levels = 1:m), factor(pred_h_states_mhmm$state, levels = 1:m))

    # Append results to the data frame
    results_df <- rbind(results_df, data.frame(
        Horizon = h,
        CrossEntropy_medHMM = cross_entropy_medhmm,
        CrossEntropy_mHMM = cross_entropy_mhmm,
        PRAUC_medHMM = pr_auc_medhmm,
        PRAUC_mHMM = pr_auc_mhmm,
        Accuracy_medHMM = accuracy_medhmm,
        Accuracy_mHMM = accuracy_mhmm,
        F1Score_medHMM = f1_score_medhmm,
        F1Score_mHMM = f1_score_mhmm
    ))
}

# Return the results data frame
results_df

# results_df_dur_s1_dep_s1_occ_s2_ind_s3 <- results_df
# results_df_dur_s5_dep_s1_occ_s2_ind_s3 <- results_df
# results_df_dur_s9_dep_s1_occ_s2_ind_s3 <- results_df



# Checks
mHMM_cont$gamma_prob_bar %>%
    as.data.frame() %>%
    mutate(iter = row_number()) %>%
    filter(iter < n_iter) %>%
    gather(state, value, -iter) %>%
    ggplot(aes(x=iter, y = value)) +
    geom_line() +
    facet_wrap(state~.)

medHMM_cont_shiftpois$gamma_prob_bar %>%
    as.data.frame() %>%
    mutate(iter = row_number()) %>%
    filter(iter < n_iter) %>%
    gather(state, value, -iter) %>%
    ggplot(aes(x=iter, y = value)) +
    geom_line() +
    facet_wrap(state~.)

exp(medHMM_cont_shiftpois$dwell_mu_bar) %>%
    as.data.frame() %>%
    mutate(iter = row_number()) %>%
    filter(iter < n_iter) %>%
    gather(state, value, -iter) %>%
    ggplot(aes(x=iter, y = value)) +
    geom_line() +
    facet_wrap(state~.)

#==============================================================================#

# Test one-step-ahead predictions:

# Step 4a: Forecast using the fitted model
# forecast_results <- forecast_mHMM1(object = out, s_data = as.matrix(forecast_data[,-4]), forecast_steps = forecast_steps)
# forecast_results_medhmm <- forecast_medHMM3(object = medHMM_cont_shiftpois,
#                                             s_data = as.matrix(fit_data[,-(n_dep+2)]),
#                                             forecast_steps = forecast_steps, Mx = Mx, return_all = TRUE)
#
# forecast_results_mhmm <- forecast_mHMM3(object = mHMM_cont, s_data = as.matrix(fit_data[,-(n_dep+2)]), forecast_steps = forecast_steps, return_all = TRUE)

h_step <- 1
i_step <- 0.8*n_t
# i_step <- 0.95*n_t
f_step <- n_t
forecast_steps <- 125

# Extract observations and states
obs <- as.data.frame(sim_data$obs)
states <- as.data.frame(sim_data$states)

# Ensure columns are properly named
colnames(obs) <- c("subj", paste0("obs",1:n_dep))
colnames(states) <- c("subj", "state")

one_step_predictions <- pbapply::pblapply(seq(i_step, f_step, h_step), function(t) {


    h_steps <- min(forecast_steps, f_step - t + 1)

    # Get correspondiong splits
    split_obs <- split_data_abs(obs, t)
    fit_data <- split_obs$fit
    forecast_data <- split_obs$forecast

    # Add a time column to the true states

    train_states <- split_data_abs(states, t)[[1]] %>%
        mutate(state = as.integer(state)) %>%
        group_by(subj) %>%
        mutate(time = row_number()) %>%
        ungroup()

    forecast_states <- split_data_abs(states, t)[[2]] %>%
        mutate(state = as.integer(state)) %>%
        group_by(subj) %>%
        mutate(time = row_number()) %>%
        ungroup() %>%
        filter(time <= forecast_steps)


    pred_medhmm <- forecast_medHMM1(object = medHMM_cont_shiftpois,
                     s_data = as.matrix(fit_data),
                     forecast_steps = h_steps, Mx = Mx, return_all = FALSE)

    pred_mhmm <- forecast_mHMM1(object = mHMM_cont,
                   s_data = as.matrix(fit_data),
                   forecast_steps = h_steps, return_all = FALSE)

    pred_medhmm <- inner_join(as.data.frame(pred_medhmm),
               forecast_states %>%
                   rename("true_state" = "state",
                          "horizon" = "time")) %>%
        mutate(t = t,
               model = "medhmm")
    pred_mhmm <- inner_join(as.data.frame(pred_mhmm),
               forecast_states %>%
                   rename("true_state" = "state",
                          "horizon" = "time")) %>%
        mutate(t = t,
               model = "mhmm")

    return(bind_rows(pred_medhmm, pred_mhmm))

})


one_step_predictions <- do.call(rbind, one_step_predictions)

one_step_metrics <- list("f_meas" = one_step_predictions %>%
    mutate(true_state = factor(true_state, levels = 1:m),
           state = factor(state, levels = 1:m)) %>%
    group_by(model, horizon) %>%
    yardstick::f_meas(truth = true_state,
                            estimate = state,
                            na_rm = TRUE) %>%
    spread(model, .estimate) %>%
    as.data.frame(),
    "bal_accuracy" = one_step_predictions %>%
    mutate(true_state = factor(true_state, levels = 1:m),
           state = factor(state, levels = 1:m)) %>%
    group_by(model, horizon) %>%
    yardstick::bal_accuracy(truth = true_state,
                                           estimate = state,
                                           na_rm = TRUE) %>%
    spread(model, .estimate) %>%
    as.data.frame(),
    "pr_auc" = one_step_predictions %>%
    as.data.frame() %>%
    mutate(true_state = factor(true_state, levels = 1:m)) %>%
    group_by(model, horizon) %>%
    yardstick::pr_auc(true_state,
                       pr_state_1:pr_state_3,
                       na_rm = TRUE) %>%
    spread(model, .estimate) %>%
    as.data.frame(),
    "roc_auc" = one_step_predictions %>%
    as.data.frame() %>%
    mutate(true_state = factor(true_state, levels = 1:m)) %>%
    group_by(model, horizon) %>%
    yardstick::roc_auc(true_state,
                       pr_state_1:pr_state_3,
                       na_rm = TRUE) %>%
    spread(model, .estimate) %>%
    as.data.frame()
)

one_step_metrics

# Lo
inner_join(as.data.frame(pred_medhmm),forecast_states)

inner_join(as.data.frame(pred_medhmm),
           forecast_states %>%
               rename("true_state" = "state",
                      "horizon" = "time"))

forecast_results_medhmm


#==============================================================================#

# Safe version:

library(dplyr)
library(purrr)
library(pbapply)
library(yardstick)

# Set steps
h_step <- 1
i_step <- 0.8 * n_t
# i_step <- 0.95 * n_t
f_step <- n_t
forecast_steps <- 125

# Extract observations and states
obs <- as.data.frame(sim_data$obs)
states <- as.data.frame(sim_data$states)

# Ensure columns are properly named
colnames(obs) <- c("subj", paste0("obs", 1:n_dep))
colnames(states) <- c("subj", "state")

# Safe version of the forecasting functions
safe_forecast_medHMM <- function(...) {
    tryCatch(forecast_medHMM1(...), error = function(e) NULL)
}
safe_forecast_mHMM <- function(...) {
    tryCatch(forecast_mHMM1(...), error = function(e) NULL)
}

# Main loop for one-step predictions
one_step_predictions <- pbapply::pblapply(seq(i_step, f_step, h_step), function(t) {
    tryCatch({
        h_steps <- min(forecast_steps, f_step - t + 1)

        # Get corresponding splits
        split_obs <- split_data_abs(obs, t)
        fit_data <- split_obs$fit
        forecast_data <- split_obs$forecast

        # Add a time column to the true states
        train_states <- split_data_abs(states, t)[[1]] %>%
            mutate(state = as.integer(state)) %>%
            group_by(subj) %>%
            mutate(time = row_number()) %>%
            ungroup()

        forecast_states <- split_data_abs(states, t)[[2]] %>%
            mutate(state = as.integer(state)) %>%
            group_by(subj) %>%
            mutate(time = row_number()) %>%
            ungroup() %>%
            filter(time <= forecast_steps)

        # Forecast using safe functions
        pred_medhmm <- safe_forecast_medHMM(
            object = medHMM_cont_shiftpois,
            s_data = as.matrix(fit_data),
            forecast_steps = h_steps, Mx = Mx, return_all = FALSE
        )

        pred_mhmm <- safe_forecast_mHMM(
            object = mHMM_cont,
            s_data = as.matrix(fit_data),
            forecast_steps = h_steps, return_all = FALSE
        )

        # Check if the predictions are NULL and handle accordingly
        if (is.null(pred_medhmm) || is.null(pred_mhmm)) return(NULL)

        pred_medhmm <- inner_join(as.data.frame(pred_medhmm),
                                  forecast_states %>%
                                      rename("true_state" = "state", "horizon" = "time")
        ) %>%
            mutate(t = t, model = "medhmm")

        pred_mhmm <- inner_join(as.data.frame(pred_mhmm),
                                forecast_states %>%
                                    rename("true_state" = "state", "horizon" = "time")
        ) %>%
            mutate(t = t, model = "mhmm")

        # Return combined results
        return(bind_rows(pred_medhmm, pred_mhmm))

    }, error = function(e) {
        # Handle errors in the main loop without stopping the simulation
        message(paste("Error at time step", t, ":", e$message))
        return(NULL)
    })
})

# Combine predictions
one_step_predictions <- do.call(rbind, one_step_predictions)
if (is.null(one_step_predictions)) stop("No valid predictions generated.")

# Compute metrics with error handling
compute_metrics_safe <- function(predictions, metric_func, ...) {
    tryCatch({
        metric_func(...) %>%
            spread(model, .estimate) %>%
            as.data.frame()
    }, error = function(e) {
        message("Error in metric computation: ", e$message)
        return(NULL)
    })
}

# One-step metrics with safe computation
one_step_metrics <- list(
    "f_meas" = compute_metrics_safe(one_step_predictions,
                                    function() one_step_predictions %>%
                                        mutate(true_state = factor(true_state, levels = 1:m),
                                               state = factor(state, levels = 1:m)) %>%
                                        group_by(model, horizon) %>%
                                        yardstick::f_meas(truth = true_state, estimate = state, na_rm = TRUE)
    ),
    "bal_accuracy" = compute_metrics_safe(one_step_predictions,
                                          function() one_step_predictions %>%
                                              mutate(true_state = factor(true_state, levels = 1:m),
                                                     state = factor(state, levels = 1:m)) %>%
                                              group_by(model, horizon) %>%
                                              yardstick::bal_accuracy(truth = true_state, estimate = state, na_rm = TRUE)
    ),
    "pr_auc" = compute_metrics_safe(one_step_predictions,
                                    function() one_step_predictions %>%
                                        as.data.frame() %>%
                                        mutate(true_state = factor(true_state, levels = 1:m)) %>%
                                        group_by(model, horizon) %>%
                                        yardstick::pr_auc(true_state, pr_state_1:pr_state_3, na_rm = TRUE)
    ),
    "roc_auc" = compute_metrics_safe(one_step_predictions,
                                     function() one_step_predictions %>%
                                         as.data.frame() %>%
                                         mutate(true_state = factor(true_state, levels = 1:m)) %>%
                                         group_by(model, horizon) %>%
                                         yardstick::roc_auc(true_state, pr_state_1:pr_state_3, na_rm = TRUE)
    )
)

# Check if metrics were computed properly
one_step_metrics <- one_step_metrics[!sapply(one_step_metrics, is.null)]

one_step_metrics



#==============================================================================#

# Ensure forecast_results has the necessary columns
colnames(forecast_results_medhmm) <- c("subj", "state", "time", "pr_state_1", "pr_state_2", "pr_state_3")
forecast_results_medhmm <- as.data.frame(forecast_results_medhmm) %>%
    filter(time > 0) %>%
    mutate(state = as.integer(state))

colnames(forecast_results_mhmm) <- c("subj", "state", "time", "pr_state_1", "pr_state_2", "pr_state_3")
forecast_results_mhmm <- as.data.frame(forecast_results_mhmm) %>%
    filter(time > 0) %>%
    mutate(state = as.integer(state))

# Return the evaluation results
evaluation_results <- list("medHMM" = evaluate_forecasts(forecast_results_medhmm, true_states, m, stepsize = 1, max_time = 25),
                           "mHMM" = evaluate_forecasts(forecast_results_mhmm, true_states, m, stepsize = 1, max_time = 25))

# Step 4b: Forecast using the fitted model
# forecast_results <- forecast_mHMM1(object = out, s_data = as.matrix(forecast_data[,-4]), forecast_steps = forecast_steps)
forecast_results_medhmm <- forecast_medHMM1(object = out_medhmm, s_data = as.matrix(fit_data[,-4]), forecast_steps = forecast_steps, Mx = Mx, return_all = TRUE)
forecast_results_mhmm <- forecast_mHMM1(object = out_mhmm, s_data = as.matrix(fit_data[,-4]), forecast_steps = forecast_steps, return_all = TRUE)

# Ensure forecast_results has the necessary columns
colnames(forecast_results_medhmm) <- c("subj", "state", "time", "pr_state_1", "pr_state_2", "pr_state_3")
forecast_results_medhmm <- as.data.frame(forecast_results_medhmm) %>%
    group_by(subj) %>%
    mutate(time = row_number()) %>%
    mutate(state = as.integer(state))

colnames(forecast_results_mhmm) <- c("subj", "state", "time", "pr_state_1", "pr_state_2", "pr_state_3")
forecast_results_mhmm <- as.data.frame(forecast_results_mhmm) %>%
    group_by(subj) %>%
    mutate(time = row_number()) %>%
    mutate(state = as.integer(state))

seen_states <- split_data(as.data.frame(data_cont$states), fit_fraction)[[1]] %>%
    mutate(state = as.integer(state)) %>%
    group_by(subj) %>%
    mutate(time = row_number()) %>%
    ungroup()











#==============================================================================#
# Test simulation function

# Parameters
n_t = 250
n = 25
m = 3
n_dep = 2
gamma_medhmm = matrix(c(0, 0.7, 0.3,
                        0.5, 0, 0.5,
                        0.6, 0.4, 0), nrow = m, ncol = m, byrow = TRUE)
gamma_mhmm = matrix(c(0.95, 0.035, 0.015,
                      0.5, 0.9, 0.5,
                      0.006, 0.004, 0.99), nrow = m, ncol = m, byrow = TRUE)
emiss_distr = list(matrix(c(10,2,
                            50,2,
                            2,2), nrow = m, ncol = 2, byrow = TRUE),
                   matrix(c(5,2,
                            20,2,
                            5,2), nrow = m, ncol = 2, byrow = TRUE))
dwell_distr = matrix(log(c(10,
                           5,
                           30)), nrow = m, ncol = 1, byrow = TRUE)
var_gamma = 0.1
var_emiss = c(10,5)
var_dwell = 0.01
hyp_prior_emiss =  prior_emiss_cont(
    gen = list(m = m, n_dep = n_dep),
    emiss_mu0 = list(matrix(c(10, 50, 2), nrow = 1),
                     matrix(c(5, 20, 5), nrow = 1)),
    emiss_K0 = list(1, 1),
    emiss_V =  list(rep(10, m), rep(25, m)),
    emiss_nu = list(1, 1),
    emiss_a0 = list(rep(0.01, m), rep(0.01, m)),
    emiss_b0 = list(rep(0.01, m), rep(0.01, m)))

hyp_prior_dwell = list(
    dwell_mu0 = matrix(log(c(10,5,30)), nrow = 1, ncol = 3), # nrow = number of covariates + 1; ncol = number of hidden states
    dwell_K0  = c(0.01),
    dwell_nu  = c(1),
    dwell_V   = rep(0.1, m)
)

fit_fraction = 0.8
forecast_steps = 50
n_iter = 500
burn_in = 250
Mx = 50


# Run
set.seed(42)
data_cont <- medHMM::sim_medHMM(n_t, n, data_distr = 'continuous', m, n_dep = n_dep,
                                dwell_distr = dwell_distr, dwell_type = 'poisson',
                                gamma = gamma_medhmm, emiss_distr = emiss_distr,
                                var_gamma = var_gamma, var_emiss = var_emiss, var_dwell = var_dwell, return_ind_par = TRUE)

set.seed(42)
data_cont2 <- medHMM::sim_medHMM2(n_t, n, data_distr = 'continuous', m, n_dep = n_dep,
                                  dwell_distr = dwell_distr, dwell_type = 'poisson',
                                  gamma = gamma_medhmm, emiss_distr = emiss_distr,
                                  var_gamma = var_gamma, var_emiss = var_emiss, var_dwell = var_dwell, return_ind_par = TRUE)


# Overall:
emp_dur_old <- data.frame("model" = "old",
                           "state" = rle(data_cont$states[,2])$values,
                           "duration" = rle(data_cont$states[,2])$lengths)

emp_dur_new <- data.frame("model" = "new",
                           "state" = rle(data_cont2$states[,2])$values,
                           "duration" = rle(data_cont2$states[,2])$lengths)

# Plot distributions
emp_dur <- bind_rows(emp_dur_old,emp_dur_new)

# Density
emp_dur %>%
    # filter(duration > 25) %>%
    ggplot(aes(x = duration, colour = model)) +
    geom_density() +
    # geom_histogram(position = "dodge") +
    scale_colour_viridis_d() +
    facet_grid(state~.) +
    theme_minimal()

# Histogram
emp_dur %>%
    # filter(duration > 25) %>%
    ggplot(aes(x = duration, colour = model, fill = model)) +
    # geom_density() +
    geom_histogram(position = "dodge") +
    scale_colour_viridis_d() +
    scale_fill_viridis_d() +
    facet_grid(state~.) +
    theme_minimal()

# Summary statistics
emp_dur %>%
    # filter(duration > 25) %>%
    group_by(model, state) %>%
    reframe(mean_dur = mean(duration), median_dur = median(duration), sd_dur = sd(duration))



# Step 3: Fit the HMM model
out_medhmm <- medHMM_cont_shiftpois(s_data = data_cont$obs,
                                    gen = list(m = m, n_dep = n_dep),
                                    start_val = c(list(gamma_medhmm), emiss_distr, list(exp(dwell_distr))),
                                    emiss_hyp_prior = hyp_prior_emiss,
                                    dwell_hyp_prior = hyp_prior_dwell,
                                    shift = NULL,
                                    show_progress = TRUE,
                                    mcmc = list(J = n_iter, burn_in = burn_in),
                                    return_path = TRUE,
                                    max_dwell = Mx)

out_mhmm <- mHMM(s_data = data_cont$obs,
                 data_distr = 'continuous',
                 gen = list(m = m, n_dep = n_dep),
                 start_val = c(list(gamma_mhmm), emiss_distr),
                 emiss_hyp_prior = hyp_prior_emiss,
                 show_progress = TRUE,
                 return_path = TRUE,
                 mcmc = list(J = n_iter, burn_in = burn_in))

out_medhmm2 <- medHMM_cont_shiftpois(s_data = data_cont2$obs,
                                    gen = list(m = m, n_dep = n_dep),
                                    start_val = c(list(gamma_medhmm), emiss_distr, list(exp(dwell_distr))),
                                    emiss_hyp_prior = hyp_prior_emiss,
                                    dwell_hyp_prior = hyp_prior_dwell,
                                    shift = NULL,
                                    show_progress = TRUE,
                                    mcmc = list(J = n_iter, burn_in = burn_in),
                                    return_path = TRUE,
                                    max_dwell = Mx)

out_mhmm2 <- mHMM(s_data = data_cont2$obs,
                 data_distr = 'continuous',
                 gen = list(m = m, n_dep = n_dep),
                 start_val = c(list(gamma_mhmm), emiss_distr),
                 emiss_hyp_prior = hyp_prior_emiss,
                 show_progress = TRUE,
                 return_path = TRUE,
                 mcmc = list(J = n_iter, burn_in = burn_in))






# Overall:
emp_dur_1 <- data.frame("model" = "true",
                        "type" = "old",
                        "state" = rle(data_cont$states[,2])$values,
                        "duration" = rle(data_cont$states[,2])$lengths)

emp_dur_2 <- data.frame("model" = "true",
                        "type" = "new",
                        "state" = rle(data_cont2$states[,2])$values,
                        "duration" = rle(data_cont2$states[,2])$lengths)

emp_dur_3 <- data.frame("model" = "medhmm",
                        "type" = "old",
                        "state" = rle(local_decoding(object = out_medhmm)[,2])$values,
                        "duration" = rle(local_decoding(object = out_medhmm)[,2])$lengths)

emp_dur_4 <- data.frame("model" = "medhmm",
                        "type" = "new",
                        "state" = rle(local_decoding(object = out_medhmm2)[,2])$values,
                        "duration" = rle(local_decoding(object = out_medhmm2)[,2])$lengths)

emp_dur_5 <- data.frame("model" = "mhmm",
                        "type" = "old",
                        "state" = rle(local_decoding(object = out_mhmm)[,2])$values,
                        "duration" = rle(local_decoding(object = out_mhmm)[,2])$lengths)

emp_dur_6 <- data.frame("model" = "mhmm",
                        "type" = "new",
                        "state" = rle(local_decoding(object = out_mhmm2)[,2])$values,
                        "duration" = rle(local_decoding(object = out_mhmm2)[,2])$lengths)


# Plot distributions
emp_dur <- bind_rows(emp_dur_1,emp_dur_2,emp_dur_3,emp_dur_4,emp_dur_5,emp_dur_6)

# Density
emp_dur %>%
    # filter(duration > 25) %>%
    ggplot(aes(x = duration, colour = model)) +
    geom_density() +
    # geom_histogram(position = "dodge") +
    scale_colour_viridis_d() +
    facet_grid(type~state) +
    theme_minimal()

# Histogram
emp_dur %>%
    # filter(duration > 25) %>%
    ggplot(aes(x = duration, colour = model, fill = model)) +
    # geom_density() +
    geom_histogram(position = "dodge") +
    scale_colour_viridis_d() +
    scale_fill_viridis_d() +
    facet_grid(type~state) +
    theme_minimal()

# Summary statistics
emp_dur %>%
    group_by(model,type, state) %>%
    reframe(mean_dur = mean(duration), median_dur = median(duration), sd_dur = sd(duration))




# Check chains
out_medhmm$gamma_prob_bar %>%
    as.data.frame() %>%
    mutate(iter = row_number()) %>%
    gather(state, value, -iter) %>%
    ggplot(aes(x = iter, y = value)) +
    geom_line() +
    facet_wrap(state~.)

out_medhmm2$gamma_prob_bar %>%
    as.data.frame() %>%
    mutate(iter = row_number()) %>%
    gather(state, value, -iter) %>%
    ggplot(aes(x = iter, y = value)) +
    geom_line() +
    facet_wrap(state~.)

out_mhmm$gamma_prob_bar %>%
    as.data.frame() %>%
    mutate(iter = row_number()) %>%
    gather(state, value, -iter) %>%
    ggplot(aes(x = iter, y = value)) +
    geom_line() +
    facet_wrap(state~.)

out_mhmm2$gamma_prob_bar %>%
    as.data.frame() %>%
    mutate(iter = row_number()) %>%
    gather(state, value, -iter) %>%
    ggplot(aes(x = iter, y = value)) +
    geom_line() +
    facet_wrap(state~.)




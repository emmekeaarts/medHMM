#-----------------------------------------------------------------------#
#               For shifet Poisson-lognormal dwell times                #
#-----------------------------------------------------------------------#

library(tidyverse)
library(mHMMbayes)

# ## 3 states
# n_t <- 500
# n <- 20
# m <- 3
# n_dep <- 2
#
#
# gamma <- matrix(c(0, 0.7, 0.3,
#                   0.5, 0, 0.5,
#                   0.6, 0.4, 0), nrow = m, ncol = m, byrow = TRUE)
#
# emiss_distr <- list(matrix(c(10,2,
#                              50,2,
#                              2,2), nrow = m, ncol = 2, byrow = TRUE),
#                     matrix(c(-5,2,
#                              -20,2,
#                              5,2), nrow = m, ncol = 2, byrow = TRUE))
#
# dwell_distr <- dwell_start <- matrix(log(c(5,1,
#                                            2,1,
#                                            10,1)), nrow = m, ncol = 2, byrow = TRUE)
#
# # Simulate data
# set.seed(42)
# sim_data_pois <- sim_data <- medHMM::sim_medHMM(n_t, n, data_distr = 'continuous', m, n_dep = n_dep,
#                                                 dwell_distr = matrix(dwell_distr[,1], ncol = 1), dwell_type = 'poisson',
#                                                 start_state = NULL, q_emiss = NULL, gamma = gamma, emiss_distr = emiss_distr, xx_vec = NULL, beta = NULL,
#                                                 var_gamma = 0.1, var_emiss = c(0.1,0.1), var_dwell = 0.1, return_ind_par = TRUE)
#
#
# # Specify hyper-prior for the continuous emission distribution
# emiss_hyp_pr <- list(
#     emiss_mu0 = list(matrix(c(10,50,2), nrow = 1),
#                      matrix(c(-5, -20, 5), nrow = 1)),
#     emiss_K0  = list(1, 1),
#     emiss_nu  = list(1, 1),
#     emiss_V   = list(rep(10, m), rep(10, m)),
#     emiss_a0  = list(rep(0.01, m), rep(0.01, m)),
#     emiss_b0  = list(rep(0.01, m), rep(0.01, m))
# )
#
# ## Define hyper-priors
# dwell_hyp_pr_plnorm <- list(
#     dwell_mu0 = matrix(log(c(5,2,10)), nrow = 1, ncol = 3), # nrow = number of covariates + 1; ncol = number of hidden states
#     dwell_K0  = c(1),
#     dwell_nu  = c(1),
#     dwell_V   = rep(0.1, m)
# )
#
#
# # Run the model on the simulated data:
# out_pois <- medHMM_cont_pois(s_data = sim_data$obs,
#                              gen = list(m = m, n_dep = n_dep),
#                              start_val = c(list(gamma), emiss_distr, list(matrix(dwell_start[,1], ncol = 1))),
#                              emiss_hyp_prior = emiss_hyp_pr,
#                              dwell_hyp_prior = dwell_hyp_pr_plnorm,
#                              show_progress = TRUE,
#                              mcmc = list(J = 500, burn_in = 350), return_path = TRUE, max_dwell = 50)
#
# out_shiftpois <- medHMM_cont_shiftpois(s_data = sim_data$obs,
#                                        gen = list(m = m, n_dep = n_dep),
#                                        start_val = c(list(gamma), emiss_distr, list(matrix(dwell_start[,1], ncol = 1))),
#                                        emiss_hyp_prior = emiss_hyp_pr,
#                                        dwell_hyp_prior = dwell_hyp_pr_plnorm,
#                                        show_progress = TRUE,
#                                        mcmc = list(J = 500, burn_in = 350), return_path = TRUE, max_dwell = 50)
#
# # Check results
# apply(out_pois$dwell_mu_bar[301:500,],2,median)
# apply(out_shiftpois$dwell_mu_bar[301:500,],2,median)
# log(c(5,2,10))
#
# bind_rows(
#     # data.frame(state = paste0("mu_d_bar",1:3),
#     #            value = c(5,2,10),
#     #            model = "true values"),
#     exp(out_pois$dwell_mu_bar[101:500,]) |>
#         as.data.frame() |>
#         gather(state, value) |>
#         mutate(model = "Poisson"),
#     exp(out_shiftpois$dwell_mu_bar[101:500,]) |>
#         as.data.frame() |>
#         gather(state, value) |>
#         mutate(model = "shifted-Poisson")) |>
#     mutate(truth = case_when(state == "mu_d_bar1" ~ 5,
#                              state == "mu_d_bar2" ~ 2,
#                              state == "mu_d_bar3" ~ 10)) |>
#     ggplot(aes(state, value, fill = model)) +
#     geom_boxplot() +
#     geom_hline(aes(yintercept = truth), linetype = "dashed") +
#     facet_wrap(state~., scales = "free_x") +
#     theme_minimal()





source("R/fb_utils.R")
source("R/utils.r")

library(mvtnorm)
library(MCMCpack)

library(Rcpp)

cppFunction('List mult_ed_fb_cpp(int m, int n, NumericVector delta, NumericMatrix allprobs, int Mx, IntegerVector Mx2, NumericMatrix gamma, NumericMatrix d, IntegerVector S, IntegerVector S2) {
                int j, t, i, u, uMax, v, k, Len;

                int zer = 0;
                NumericMatrix d2 = clone(d);

                double x;
                NumericMatrix D2(m, n);
                NumericVector dSum(m);

                NumericVector N(n);
                NumericMatrix Norm(m, n);
                NumericMatrix Forward(m, n);
                NumericMatrix StateIn(m, n);
                double Observ = 0;

                NumericMatrix Backward(m, n);
                NumericMatrix B_star(m, n + 2);
                IntegerVector VarL(Mx - 1);
                for (i = 1; i < Mx; i++) {
                    VarL(i - 1) = Mx - i;
                }
                int occNcol = Mx * (n - Mx + 1) + sum(VarL);
                NumericMatrix Occupancy(m, occNcol);

                IntegerVector lengthID(n);
                for (i = 0; i < n; i++) {
                    if (i < n - Mx + 1) {
                        lengthID(i) = Mx;
                    } else {
                        lengthID(i) = VarL(i - (n - Mx + 1));
                    }
                }

                IntegerVector endID(n + 2);
                for (i = 0; i < n + 2; i++) {
                    if (i == 0) {
                        endID(i) = 0;
                    }
                    if (i > 0) {
                        if (i < n + 1) {
                            endID(i) = endID(i - 1) + lengthID(i - 1);
                        }
                    }
                    if (i == n + 1) {
                        endID(i) = endID(i - 1);
                    }
                }

                IntegerVector::const_iterator first = S.begin() + 0;

                // forward recursion
                for (t = 0; t <= n - 1; t++) {

                    uMax = std::min(t + 1, Mx + 1);

                    IntegerVector::const_iterator last = S.begin() + (t + 1);
                    IntegerVector SSh(first, last);
                    Len = SSh.size();
                    IntegerVector SShrev(Len);
                    for (i = 0; i < Len; i++) {
                        SShrev(i) = SSh(Len - 1 - i);
                    }

                    IntegerVector::const_iterator first2 = SShrev.begin() + 0;
                    IntegerVector::const_iterator last2 = SShrev.begin() + (uMax);
                    IntegerVector SShrev2(first2, last2);

                    d2 = clone(d);
                    for (i = 1; i < uMax; i++) {
                        for (j = 0; j < m; j++) {
                            d2(j, i) *= SShrev2(i - 1);
                        }
                    }

                    for (j = 0; j < m; j++) {
                        dSum(j) = 0;
                        for (i = 0; i <= Mx; i++) {
                            dSum(j) += d2(j, i);
                        }
                    }

                    for (j = 0; j < m; j++) {
                        if (dSum(j) != 0) {
                            for (i = 0; i <= Mx; i++) {
                                d2(j, i) /= dSum(j);
                            }
                        }
                    }

                    if (t == n - 1) {
                        for (j = 0; j < m; j++) {
                            for (u = 1; u <= Mx; u++) {
                                x = 0;
                                for (v = u; v < Mx + 1; v++)
                                    x += d2(j, v);
                                    D2(j, (u - 1)) = x;
                            }
                            for (u = Mx + 1; u <= n; u++) {
                                D2(j, (u - 1)) = 0;
                            }
                        }
                    }

                    N(t) = 0;
                    for (j = 0; j < m; j++) {
                        if (t == 0) {
                            // Calculate initial state probabilities at t = 0
                            if (std::any_of(allprobs(_, 0).cbegin(), allprobs(_, 0).cend(), NumericVector::is_na)) {
                                Norm(j, 0) = log(delta(j));
                            } else {
                                Norm(j, 0) = log(delta(j)) + log(allprobs(j, 0));
                            }
                        } else {
                            // Check for missing values in the row at time t of allprobs
                            if (std::any_of(allprobs(_, t).cbegin(), allprobs(_, t).cend(), NumericVector::is_na)) {
                                Norm(j, t) = log(std::abs(exp(StateIn(j, t)) - exp(Forward(j, (t - 1))) + exp(Norm(j, (t - 1)))));
                            } else {
                                Norm(j, t) = log(allprobs(j, t)) + log(std::abs(exp(StateIn(j, t)) - exp(Forward(j, (t - 1))) + exp(Norm(j, (t - 1)))));
                            }
                        }
                        N(t) += exp(Norm(j, t));
                    }
                    N(t) = log(N(t));
                    for (j = 0; j < m; j++) {
                        Norm(j, t) -= N(t);
                    }

                    for (j = 0; j < m; j++) {
                        Forward(j, t) = 0;
                        Observ = 0;

                        if (t < n - 1) {
                            for (u = 1; u <= std::min(t + 1, Mx2(j)); u++) {
                                // Check for missing values in allprobs
                                if (std::any_of(allprobs(_, t - u + 1).cbegin(), allprobs(_, t - u + 1).cend(), NumericVector::is_na)) {
                                    if (SShrev2(u - 1) == 1) {
                                        if (u < t + 1) {
                                            // Forward(j, t) += exp(log(d2(j, u)) + StateIn(j, (t - u + 1)));
                                             Forward(j, t) += exp(Observ + log(d2(j, u)) + StateIn(j, (t - u + 1)));
                                        } else {
                                            // Forward(j, t) += exp(log(d2(j, t + 1)) + log(delta(j)));
                                            Forward(j, t) += exp(Observ + log(d2(j, t + 1)) + log(delta(j)));
                                        }
                                    }
                                } else {
                                    Observ += log(allprobs(j, t - u + 1)) - N(t - u + 1);
                                    if (SShrev2(u - 1) == 1) {
                                        if (u < t + 1) {
                                            Forward(j, t) += exp(Observ + log(d2(j, u)) + StateIn(j, (t - u + 1)));
                                        } else {
                                            Forward(j, t) += exp(Observ + log(d2(j, t + 1)) + log(delta(j)));
                                        }
                                    }
                                }
                            }
                            Forward(j, t) = log(Forward(j, t));
                        } else {
                            for (u = 1; u <= std::min(n, Mx2(j)); u++) {
                                if (std::any_of(allprobs(_, t - u + 1).cbegin(), allprobs(_, t - u + 1).cend(), NumericVector::is_na)) {
                                    if (SShrev2(u - 1) == 1) {
                                        if (u < n) {
                                            // Forward(j, n - 1) += exp(log(D2(j, u)) + StateIn(j, n - u));
                                            Forward(j, n - 1) += exp(Observ + log(D2(j, u)) + StateIn(j, n - u));
                                        } else {
                                            // Forward(j, n - 1) += exp(log(D2(j, n)) + log(delta(j)));
                                            Forward(j, n - 1) += exp(Observ + log(D2(j, n)) + log(delta(j)));
                                        }
                                    }
                                } else {
                                    Observ += log(allprobs(j, t - u + 1)) - N(t - u + 1);
                                    if (SShrev2(u - 1) == 1) {
                                        if (u < n) {
                                            Forward(j, n - 1) += exp(Observ + log(D2(j, u)) + StateIn(j, n - u));
                                        } else {
                                            Forward(j, n - 1) += exp(Observ + log(D2(j, n)) + log(delta(j)));
                                        }
                                    }
                                }
                            }
                            Forward(j, n - 1) = log(Forward(j, n - 1));
                        }
                    }
                    if (t < n - 1) {
                        for (j = 0; j < m; j++) {
                            StateIn(j, t + 1) = 0;
                            for (i = 0; i < m; i++) {
                                StateIn(j, t + 1) += exp(Forward(i, t) + log(gamma(i, j)));
                            }
                            StateIn(j, t + 1) = log(StateIn(j, t + 1));
                        }
                    }
                }

                // Backward recursion

                for (t = n - 1; t >= 0; t--) {
                    if (S(t) == 1) {

                        uMax = std::min(n - t, Mx);
                        IntegerVector::const_iterator first3 = S2.begin() + t;
                        IntegerVector::const_iterator last3 = S2.begin() + (n);
                        IntegerVector SShB(first3, last3);

                        IntegerVector::const_iterator first4 = SShB.begin() + 0;
                        IntegerVector::const_iterator last4 = SShB.begin() + (uMax);
                        IntegerVector SShB2(first4, last4);

                        d2 = clone(d);
                        for (i = 1; i < uMax + 1; i++) {
                            for (j = 0; j < m; j++) {
                                d2(j, i) *= SShB2(i - 1);
                            }
                        }

                        for (j = 0; j < m; j++) {
                            dSum(j) = 0;
                            for (i = 0; i <= Mx; i++) {
                                dSum(j) += d2(j, i);
                            }
                        }

                        for (j = 0; j < m; j++) {
                            if (dSum(j) != 0) {
                                for (i = 0; i <= Mx; i++) {
                                    d2(j, i) /= dSum(j);
                                }
                            }
                        }

                        for (j = 0; j < m; j++) {
                            for (u = 1; u <= Mx; u++) {
                                x = 0;
                                for (v = u; v < Mx + 1; v++)
                                    x += d2(j, v);
                                    D2(j, (u - 1)) = x;
                            }
                            for (u = Mx + 1; u <= n; u++) {
                                D2(j, (u - 1)) = 0;
                            }
                        }

                        for (j = 0; j < m; j++) {
                            B_star(j, t) = 0;
                            Observ = 0;
                            for (u = 1; u <= std::min(n - t, Mx2(j)); u++) {
                                if (std::any_of(allprobs(_, t + u - 1).cbegin(), allprobs(_, t + u - 1).cend(), NumericVector::is_na)) {
                                    if (SShB2(u - 1) == 1) {
                                        if (u < n - t) {
                                            // Occupancy(j, endID(t) + (u - 1)) = exp(Backward(j, t + u) + log(d2(j, u)));
                                            Occupancy(j, endID(t) + (u - 1)) = exp(Backward(j, t + u) + Observ + log(d2(j, u)));
                                        } else {
                                            // Occupancy(j, endID(t) + (u - 1)) = exp(log(D2(j, n - 1 - t)));
                                            Occupancy(j, endID(t) + (u - 1)) = exp(Observ + log(D2(j, n - 1 - t)));
                                        }
                                        B_star(j, t) += Occupancy(j, endID(t) + (u - 1));
                                    }
                                } else {
                                    Observ += log(allprobs(j, t + u - 1)) - N(t + u - 1);
                                    if (SShB2(u - 1) == 1) {
                                        if (u < n - t) {
                                            Occupancy(j, endID(t) + (u - 1)) = exp(Backward(j, t + u) + Observ + log(d2(j, u)));
                                        } else {
                                            Occupancy(j, endID(t) + (u - 1)) = exp(Observ + log(D2(j, n - 1 - t)));
                                        }
                                        B_star(j, t) += Occupancy(j, endID(t) + (u - 1));
                                    }
                                }
                            }
                            B_star(j, t) = log(B_star(j, t));
                        }
                        for (j = 0; j < m; j++) {
                            Backward(j, t) = 0;
                            for (k = 0; k < m; k++) {
                                Backward(j, t) += exp(B_star(k, t) + log(gamma(j, k)));
                            }
                            Backward(j, t) = log(Backward(j, t));
                        }
                    }
                }

                List H(n);
                for (i = 0; i < n; i++) {
                    if (S(i) == 0) {
                        H(i) = zer;
                    } else {
                        NumericMatrix foo(m, lengthID(i));
                        for (k = 0; k < lengthID(i); k++) {
                            foo(_, k) = Occupancy(_, k + endID(i));
                        }
                        H(i) = foo;
                    }
                }

                return List::create(N, B_star, H);
            }
            ')



# load("/Users/a6159737/Documents/Utrecht University/PhD/Collaborations/Gabriel/eeg_data/tests/tests.RData")
#
# for(i in 1:m){
#     start_gamma[i,i] <- 0
#     start_gamma[i,] <- start_gamma[i,]/sum(start_gamma[i,])
# }



library(tidyverse)
library(mHMMbayes)

## 3 states
n_t <- 500
n <- 20
m <- 3
n_dep <- 2


gamma <- matrix(c(0, 0.7, 0.3,
                  0.5, 0, 0.5,
                  0.6, 0.4, 0), nrow = m, ncol = m, byrow = TRUE)

emiss_distr <- list(matrix(c(10,2,
                             50,2,
                             2,2), nrow = m, ncol = 2, byrow = TRUE),
                    matrix(c(-5,2,
                             -20,2,
                             5,2), nrow = m, ncol = 2, byrow = TRUE))

dwell_distr <- dwell_start <- matrix(log(c(5,1,
                                           2,1,
                                           10,1)), nrow = m, ncol = 2, byrow = TRUE)

# Simulate data
set.seed(42)
sim_data_pois <- sim_data <- medHMM::sim_medHMM(n_t, n, data_distr = 'continuous', m, n_dep = n_dep,
                                                dwell_distr = matrix(dwell_distr[,1], ncol = 1), dwell_type = 'poisson',
                                                start_state = NULL, q_emiss = NULL, gamma = gamma, emiss_distr = emiss_distr, xx_vec = NULL, beta = NULL,
                                                var_gamma = 0.1, var_emiss = c(0.1,0.1), var_dwell = 0.1, return_ind_par = TRUE)


# Specify hyper-prior for the continuous emission distribution
emiss_hyp_pr <- list(
    emiss_mu0 = list(matrix(c(10,50,2), nrow = 1),
                     matrix(c(-5, -20, 5), nrow = 1)),
    emiss_K0  = list(1, 1),
    emiss_nu  = list(1, 1),
    emiss_V   = list(rep(10, m), rep(10, m)),
    emiss_a0  = list(rep(0.01, m), rep(0.01, m)),
    emiss_b0  = list(rep(0.01, m), rep(0.01, m))
)

## Define hyper-priors
dwell_hyp_pr_plnorm <- list(
    dwell_mu0 = matrix(log(c(5,2,10)), nrow = 1, ncol = 3), # nrow = number of covariates + 1; ncol = number of hidden states
    dwell_K0  = c(1),
    dwell_nu  = c(1),
    dwell_V   = rep(0.1, m)
)


set.seed(42)
sim_data$obs[sample(1:nrow(sim_data$obs), nrow(sim_data$obs)*0.1),-1] <- NA

s_data = sim_data$obs
gen = list(m = m, n_dep = n_dep)
start_val = c(list(gamma), emiss_distr, list(dwell_distr))
emiss_hyp_prior = emiss_hyp_pr
dwell_hyp_prior = dwell_hyp_pr_plnorm
show_progress = TRUE
shift = NULL
mcmc = list(J = 500, burn_in = 250)
gamma_hyp_prior = NULL
xx = NULL
gamma_sampler = NULL
dwell_sampler = NULL
return_path = TRUE
max_dwell = 50





medHMM_cont_shiftpois <- function(s_data, gen, xx = NULL, start_val,
                                  # dwell_family = "poisson",
                                  emiss_hyp_prior, dwell_hyp_prior, shift = NULL,
                                  mcmc, return_path = FALSE, show_progress = TRUE,
                                  gamma_hyp_prior = NULL, gamma_sampler = NULL, dwell_sampler = NULL, max_dwell = NULL
){

    # Initialize data -----------------------------------
    # dependent variable(s), sample size, dimensions gamma and conditional distribuiton
    n_dep			 <- gen$n_dep
    dep_labels <- colnames(s_data[,2:(n_dep+1)])
    id         <- unique(s_data[,1])
    n_subj     <- length(id)
    subj_data  <- rep(list(NULL), n_subj)
    if(sum(sapply(s_data, is.factor)) > 0 ){
        stop("Your data contains factorial variables, which cannot be used as input in the function mHMM. All variables have to be numerical.")
    }
    for(s in 1:n_subj){
        subj_data[[s]]$y <- as.matrix(s_data[s_data[,1] == id[s],][,-1], ncol = n_dep)
    }
    ypooled    <- n <- NULL
    n_vary     <- numeric(n_subj)
    m          <- gen$m
    if(is.null(shift)){
        shift <- 1
    }
    # dwell_mle       <- list(NULL)
    # dwell_mhess     <- list(NULL)
    # dwell_converge  <- list(NULL)

    for(s in 1:n_subj){
        ypooled   <- rbind(ypooled, subj_data[[s]]$y)
        n         <- dim(subj_data[[s]]$y)[1]

        switched 			<- rep(0, n)
        switched[1] 		<- 1
        for (t in 2:n) {
            if(!any(is.na(subj_data[[s]]$y[t,]) | is.na(subj_data[[s]]$y[t-1,]))){
                if(any(subj_data[[s]]$y[t,] != subj_data[[s]]$y[t-1,])) {
                    switched[t] <- 1
                }
            }
        }
        switched2 		<- c(switched[-1],1)

        if (is.null(max_dwell)){
            max_dwell_n <- n
        } else if (max_dwell > n) {
            max_dwell_n <- n
        } else {
            max_dwell_n <- max_dwell
        }

        n_vary[s] <- n
        subj_data[[s]]	<- c(
            subj_data[[s]], n = n, list(gamma_converge = numeric(m), gamma_int_mle = matrix(NA_real_, m, (m - 2)),
                                        gamma_mhess = matrix(NA_real_, (m - 2) * m, (m - 2))),
            list(Mx = max_dwell_n),
            list(Mx2 = rep(max_dwell_n, m)),
            list(switch = switched),
            list(switch2 = switched2),
            list(dwell_converge = numeric(m),
                 dwell_mle = numeric(m),
                 dwell_mhess = numeric(m))
        )
    }
    n_total 		<- dim(ypooled)[1]

    # Intialize covariates
    # n_dep1 <- 1 + n_dep
    n_dep2 <- 2 + n_dep
    # nx <- numeric(n_dep1)
    nx <- numeric(n_dep2)
    if (is.null(xx)){
        xx <- rep(list(matrix(1, ncol = 1, nrow = n_subj)), n_dep2)
        nx[] <- 1
    } else {
        if(!is.list(xx) | length(xx) != n_dep2){
            stop("If xx is specified, xx should be a list, with the number of elements equal to the number of dependent variables + 2")
        }
        for(i in 1:n_dep2){
            if (is.null(xx[[i]])){
                xx[[i]] <- matrix(1, ncol = 1, nrow = n_subj)
                nx[i] <- 1
            } else {
                nx[i] <- ncol(xx[[i]])
                if (sum(xx[[i]][,1] != 1)){
                    stop("If xx is specified, the first column in each element of xx has to represent the intercept. That is, a column that only consists of the value 1")
                }
                if(nx[i] > 1){
                    for(j in 2:nx[i]){
                        if(is.factor(xx[[i]][,j])){
                            stop("Factors currently cannot be used as covariates, see help file for alternatives")
                        }
                        if((length(unique(xx[[i]][,j])) == 2) & (sum(xx[[i]][,j] != 0 & xx[[i]][,j] !=1) > 0)){
                            stop("Dichotomous covariates in xx need to be coded as 0 / 1 variables. That is, only conisting of the values 0 and 1")
                        }
                        if(length(unique(xx[[i]][,j])) > 2){
                            xx[[i]][,j] <- xx[[i]][,j] - mean(xx[[i]][,j])
                        }
                    }
                }
            }
        }
    }

    # Initialize mcmc argumetns
    J 				<- mcmc$J
    burn_in			<- mcmc$burn_in


    # Initalize priors and hyper priors --------------------------------
    # Initialize gamma sampler
    if(is.null(gamma_sampler)) {
        gamma_int_mle0  <- rep(0, m - 2)
        gamma_scalar    <- 2.93 / sqrt(m - 2)
        gamma_w         <- .1
    } else {
        gamma_int_mle0  <- gamma_sampler$gamma_int_mle0
        gamma_scalar    <- gamma_sampler$gamma_scalar
        gamma_w         <- gamma_sampler$gamma_w
    }

    # Initialize dwell sampler
    if(is.null(dwell_sampler)) {
        dwell_mle0      <- 0
        dwell_scalar    <- 2.38
        dwell_w         <- .1
    } else {
        dwell_mle0      <- dwell_sampler$emiss_mle0
        dwell_scalar    <- dwell_sampler$emiss_scalar
        dwell_w         <- dwell_sampler$emiss_w
    }



    # Initialize Gamma hyper prior
    if(is.null(gamma_hyp_prior)){
        gamma_mu0	        <- rep(list(matrix(0,nrow = nx[1], ncol = m - 2)), m)
        gamma_K0			<- diag(1, nx[1])
        gamma_nu			<- 3 + m - 2
        gamma_V			    <- gamma_nu * diag(m - 2)
    } else {
        ###### BUILD in a warning / check if gamma_mu0 is a matrix when given, with  nrows equal to the number of covariates
        gamma_mu0			<- gamma_hyp_prior$gamma_mu0
        gamma_K0			<- diag(gamma_hyp_prior$gamma_K0, nx[1])
        gamma_nu			<- gamma_hyp_prior$gamma_nu
        gamma_V			    <- gamma_hyp_prior$gamma_V
    }


    # Initialize emiss hyper prior
    if(missing(emiss_hyp_prior)){
        stop("The hyper-prior values for the Normal emission distribution(s) denoted by emiss_hyp_prior needs to be specified")
    }

    # emiss_mu0: a list containing n_dep matrices with in the first row the hypothesized mean values of the Normal emission
    # distributions in each of the states over the m coloumns. Subsequent rows contain the hypothesised regression
    # coefficients for covariates influencing the state dependent mean value of the normal distribution
    emiss_mu0	  <- rep(list(NULL), n_dep)
    emiss_a0	  <- rep(list(NULL), n_dep)
    emiss_b0	  <- rep(list(NULL), n_dep)
    emiss_V	  <- rep(list(NULL), n_dep)
    emiss_nu	    <- rep(list(NULL), n_dep)
    emiss_K0     <- rep(list(NULL), n_dep)
    for(q in 1:n_dep){
        # emiss_hyp_prior[[q]]$emiss_mu0 has to contain a list with lenght equal to m, and each list contains matrix with number of rows equal to number of covariates for that dep. var.
        # stil build in a CHECK, with warning / stop / switch to default prior
        emiss_mu0[[q]]	<- emiss_hyp_prior$emiss_mu0[[q]]
        emiss_nu[[q]]	<- emiss_hyp_prior$emiss_nu[[q]]
        emiss_V[[q]]	<- emiss_hyp_prior$emiss_V[[q]]
        emiss_K0[[q]]	<- diag(emiss_hyp_prior$emiss_K0, nx[1 + q])
        emiss_a0[[q]] <- emiss_hyp_prior$emiss_a0[[q]]
        emiss_b0[[q]] <- emiss_hyp_prior$emiss_b0[[q]]
    }


    # Initialize dwell hyper prior =====================================================================================================
    # d_mu0			<- dwell_hyp_prior$d_mu0
    # s2_0			<- dwell_hyp_prior$s2_0
    # alpha.sigma20	<- dwell_hyp_prior$alpha.sigma20
    # beta.sigma20	<- dwell_hyp_prior$beta.sigma20
    # alpha.tau20		<- dwell_hyp_prior$alpha.tau20
    # beta.tau20		<- dwell_hyp_prior$beta.tau20

    # Initialize emiss hyper prior
    if(missing(dwell_hyp_prior)){
        stop("The hyper-prior values for the Poisson-lognormal state dwell distribution denoted by dwell_hyp_prior needs to be specified")
    }

    # emiss_mu0: for each dependent variable, emiss_mu0 is a list, with one element for each state.
    # Each element is a matrix, with number of rows equal to the number of covariates (with the intercept being one cov),
    # and the number of columns equal to q_emiss[q] - 1.
    # dwell_mu0	<- list(vector("list", m))
    # dwell_V	    <- list(vector("list", m))
    # dwell_nu	<- list(NULL)
    # dwell_K0    <- list(NULL)
    # emiss_hyp_prior[[q]]$emiss_mu0 has to contain a list with lenght equal to m, and each list contains matrix with number of rows equal to number of covariates for that dep. var.
    # stil build in a CHECK, with warning / stop / switch to default prior
    dwell_mu0	<- dwell_hyp_prior$dwell_mu0
    dwell_nu	<- dwell_hyp_prior$dwell_nu
    dwell_V	    <- dwell_hyp_prior$dwell_V
    dwell_K0	<- diag(dwell_hyp_prior$dwell_K0, nx[2 + n_dep])



    # Define objects used to store data in mcmc algorithm, not returned ----------------------------
    # overall
    c <- llk <- numeric(1)
    sample_path <- lapply(n_vary, dif_matrix, cols = J)
    trans <- rep(list(vector("list", m)), n_subj)

    # gamma
    gamma_int_mle_pooled <- gamma_pooled_ll <- vector("list", m)
    gamma_c_int <- rep(list(matrix(NA_real_, n_subj, (m-2))), m)
    gamma_mu_int_bar <- gamma_V_int <- vector("list", m)
    gamma_mu_prob_bar <- rep(list(numeric(m)), m)
    gamma_naccept <- matrix(0, n_subj, m)

    # emiss
    cond_y <- lapply(rep(n_dep, n_subj), nested_list, m = m)
    emiss_c_mu <- rep(list(rep(list(matrix(NA_real_,ncol = 1, nrow = n_subj)),n_dep)), m)
    for(i in 1:m){
        for(q in 1:n_dep){
            emiss_c_mu[[i]][[q]][,1] <- start_val[[1 + q]][i,1]
        }
    }
    emiss_V_mu <- emiss_c_mu_bar <- emiss_c_V <- rep(list(rep(list(NULL),n_dep)), m)
    ss_subj <- n_cond_y <- numeric(n_subj)
    # label_switch <- matrix(0, ncol = n_dep, nrow = m, dimnames = list(c(paste("mu_S", 1:m, sep = "")), dep_labels))

    # dwell ============================================================================================================================= Added dwell pars
    d <- B_star <- Occupancy <- N <- numeric()
    Dur <- vector("list", n_subj)
    n.Dur <- numeric(n_subj)
    dwell_c_mu <- rep(list(matrix(NA_real_,ncol = 1, nrow = n_subj)), m)
    for(i in 1:m){
        for(q in 1:n_dep){
            dwell_c_mu[[i]][,1] <- start_val[[2 + n_dep]][i,1]
        }
    }
    dwell_V_mu <- dwell_c_mu_bar <- dwell_c_V <- rep(list(NULL), m)
    dwell_mle_pooled <- dwell_pooled_ll <- rep(list(NULL), m)
    dwell_naccept <- matrix(0, n_subj, m)


    # Define objects that are returned from mcmc algorithm ----------------------------
    # Define object for subject specific posterior density, put start values on first row
    if(length(start_val) != n_dep + 2){
        stop("The number of elements in the list start_val should be equal to 2 + the number of dependent variables,
         and should not contain nested lists (i.e., lists within lists)")
    }
    PD 					  <- matrix(NA_real_, nrow = J, ncol = n_dep * m * 2 + m * m + m + 1) # ===================================================== Added dwell time pars and positions
    colnames(PD) 	<- c(paste("dep", rep(1:n_dep, each = m), "_mu", "_S", rep(1:m), sep = ""),
                       paste("dep", rep(1:n_dep, each = m), "_fixvar", "_S", rep(1:m), sep = ""),
                       paste("S", rep(1:m, each = m), "toS", rep(1:m, m), sep = ""),
                       paste("dur_lambda", 1:m, sep = ""),
                       # paste("dur_logsigma2", 1:m, sep = ""),
                       # paste("dur_tau2_bar", 1:m, sep = ""),
                       "LL")

    PD[1, ((n_dep * m * 2 + 1)) :((n_dep * m * 2 + m * m))] <- unlist(sapply(start_val, t))[1:(m*m)]
    for(q in 1:n_dep){
        PD[1, ((q-1) * m + 1):(q * m)] <- start_val[[q + 1]][,1]
        PD[1, (n_dep * m + (q-1) * m + 1):(n_dep * m + q * m)] <- start_val[[q + 1]][,2]
    }
    PD[1, ((n_dep * m * 2 + m * m + 1)) :((n_dep * m * 2 + m * m + m))] <- start_val[[1 + n_dep + 1]][,1] # Check if log or exp
    # PD[1, ((n_dep * m * 2 + m * m + m + 1)) :((n_dep * m * 2 + m * m + m + m))] <- start_val[[1 + n_dep + 1]][,2]
    # PD[1, ((n_dep * m * 2 + m * m + 1)) :((n_dep * m * 2 + m * m + m))] <- unlist(sapply(start_val, t))[(n_dep*m*2+m*m+1):(n_dep*m*2+m*m+m*2)]

    PD_subj				<- rep(list(PD), n_subj)

    # Define object for population posterior density (probabilities and regression coefficients parameterization )
    gamma_prob_bar		<- matrix(NA_real_, nrow = J, ncol = (m * m))
    colnames(gamma_prob_bar) <- paste("S", rep(1:m, each = m), "toS", rep(1:m, m), sep = "")
    gamma_prob_bar[1,] <- PD[1,(n_dep * 2 * m + 1):(n_dep * 2 * m + m * m)]
    emiss_mu_bar			<- rep(list(matrix(NA_real_, ncol = m, nrow = J, dimnames = list(NULL, c(paste("mu_", 1:m, sep = ""))))), n_dep)
    names(emiss_mu_bar) <- dep_labels
    for(q in 1:n_dep){
        emiss_mu_bar[[q]][1,] <- PD[1, ((q-1) * m + 1):(q * m)]
    }

    gamma_int_bar				<- matrix(NA_real_, nrow = J, ncol = ((m-2) * m))
    temp_gamma <- matrix(gamma_prob_bar[1,], byrow = TRUE, ncol = m, nrow = m)
    diag(temp_gamma) <- NA
    gamma_int_bar[1,] <- as.vector(prob_to_int(t(apply(matrix(as.numeric(t(temp_gamma))[!is.na(as.numeric(t(temp_gamma)))], nrow = m, byrow = TRUE),1, function(r) r/sum(r)))))

    temp <- matrix(paste0("int_S", rep(1:m, each = m), "toS", rep(1:m, m)), nrow = m, byrow = TRUE)
    diag(temp) <- NA
    colnames(gamma_int_bar) <- as.vector(t(matrix(t(temp)[!is.na(as.vector(t(temp)))], nrow = m, byrow = TRUE)[,-1]))

    gamma_V_int_bar <- matrix(NA_real_, nrow = J, ncol = ((m-2) * (m-2) * m))
    gamma_V_int_bar[1,] <- unlist(lapply(gamma_V, function(e) as.vector(t(e))))

    int_names <- paste0("int_S", rep(1:m, each = m), "toS", rep(1:m, m))[!(paste0("int_S", rep(1:m, each = m), "toS", rep(1:m, m)) %in% as.vector(t(matrix(t(temp)[!is.na(as.vector(t(temp)))], nrow = m, byrow = TRUE)[,-1])))]
    var_names <- paste0("var_int_S", rep(1:m, each = (m-1)*(m-1)), "toS", rep(2:m, each=m-1), "_with_", "int_S", rep(1:m, each = (m-1)*(m-1)), "toS", rep(2:m, m))
    colnames(gamma_V_int_bar) <- var_names[!grepl(paste(int_names, collapse = "|"), var_names)]

    if(nx[1] > 1){
        gamma_cov_bar				<- matrix(NA_real_, nrow = J, ncol = ((m-2) * m) * (nx[1] - 1))
        gamma_cov_bar[1,] <- 0
    } else{
        gamma_cov_bar <- "No covariates where used to predict the transition probability matrix"
    }
    emiss_varmu_bar			<- rep(list(matrix(NA_real_, ncol = m, nrow = J, dimnames = list(NULL, c(paste("varmu_", 1:m, sep = ""))))), n_dep)
    names(emiss_varmu_bar) <- dep_labels
    emiss_var_bar			<- rep(list(matrix(NA_real_, ncol = m, nrow = J, dimnames = list(NULL, c(paste("var_", 1:m, sep = ""))))), n_dep)
    names(emiss_var_bar) <- dep_labels
    for(q in 1:n_dep){
        emiss_var_bar[[q]][1,] <- PD[1, (n_dep * m + (q-1) * m + 1):(n_dep * m + q * m)]
    }
    if(sum(nx[2:(length(nx)-1)]) > n_dep){
        emiss_cov_bar			<- lapply(m * (nx[2:(length(nx)-1)] - 1 ), dif_matrix, rows = J)
        names(emiss_cov_bar) <- dep_labels
        for(q in 1:n_dep){
            if(nx[1 + q] > 1){
                colnames(emiss_cov_bar[[q]]) <-  paste( paste("cov", 1 : (nx[1+q] - 1), "_", sep = ""), "mu_S", rep(1:m, each = (nx[1 + q] - 1)), sep = "")
                emiss_cov_bar[[q]][1,] <- 0
            } else {
                emiss_cov_bar[[q]] <- "No covariates where used to predict the emission probabilities for this outcome"
            }
        }
    } else{
        emiss_cov_bar <- "No covariates where used to predict the emission probabilities"
    }

    # dwell ============================================================================================================================= Added dwell pars
    # duration_bar		<- matrix(NA_real_, nrow = J, ncol = m * 3)
    # colnames(duration_bar)	<- c(paste("mu_d_bar", 1:m, sep = ""), paste("logsigma2", 1:m, sep = ""), paste("tau2_d_bar", 1:m, sep = ""))
    dwell_mu_bar <- matrix(NA_real_, nrow = J, ncol = m)
    colnames(dwell_mu_bar) <- paste("mu_d_bar", 1:m, sep = "")
    dwell_mu_bar[1,] <- matrix(log(PD[1, ((n_dep * m * 2 + m * m + 1)) :((n_dep * m * 2 + m * m + m))]), nrow = 1, ncol = m, byrow = TRUE)
    # dwell_mu_bar[1,] <- matrix(PD[1, ((n_dep * m * 2 + m * m + 1)) :((n_dep * m * 2 + m * m + m))], nrow = 1, ncol = m, byrow = TRUE)

    dwell_varmu_bar <- matrix(NA_real_, nrow = J, ncol = m)
    # dwell_var_bar <- matrix(NA_real_, nrow = J, ncol = m)
    colnames(dwell_varmu_bar) <- paste("tau2_d_bar", 1:m, sep = "") # Check name
    # colnames(dwell_var_bar) <- paste("logsigma2", 1:m, sep = "")
    # dwell_var_bar[1,] <- dwell_varmu_bar[1,] <- matrix(PD[1, ((n_dep * m * 2 + m * m + m + 1)) :((n_dep * m * 2 + m * m + m + m))], nrow = 1, ncol = m, byrow = TRUE)

    if(sum(nx[length(nx)]) > 1){
        dwell_cov_bar			<- lapply(m * (nx[length(nx)]), dif_matrix, rows = J)
        names(dwell_cov_bar) <- dep_labels

        if(nx[1 + q] > 1){
            colnames(dwell_cov_bar) <-  paste( paste("cov", 1 : (nx[length(nx)]), "_", sep = ""), "mu_S", rep(1:m, each = (nx[length(nx)])), sep = "")
            dwell_cov_bar[1,] <- 0
        } else {
            dwell_cov_bar <- "No covariates where used to predict the emission probabilities for this outcome"
        }

    } else{
        dwell_cov_bar <- "No covariates where used to predict the emission probabilities"
    }

    # Define object for subject specific posterior density (regression coefficients parameterization )
    gamma_int_subj			<- rep(list(gamma_int_bar), n_subj)

    # Put starting values in place for fist run forward algorithm
    emiss_sep 	<- rep(list(matrix(NA_real_, ncol = 2, nrow = m)), n_dep)
    for(q in 1:n_dep){
        emiss_sep[[q]][,1] <- PD[1, ((q-1) * m + 1):(q * m)]
        emiss_sep[[q]][,2] <- PD[1, (n_dep * m + (q-1) * m + 1):(n_dep * m + q * m)]
    }
    emiss	<- rep(list(emiss_sep), n_subj)
    gamma 	<- rep(list(matrix(PD[1,(m * 2 * n_dep + 1):(m * 2 * n_dep + m * m)], byrow = TRUE, ncol = m)), n_subj)
    delta 	<- rep(list(solve(t(diag(m) - gamma[[1]] + 1), rep(1, m))), n_subj)
    dwell   <- rep(list(matrix(start_val[[n_dep + 2]][,1], nrow = m, ncol = 1)), n_subj)
    # dwell   <- rep(list(matrix(exp(start_val[[n_dep + 2]][,1]), nrow = m, ncol = 1)), n_subj)

    # # dwell ============================================================================================================================= Added dwell starting values
    # mu_d_bar        <- matrix(PD[1, ((n_dep * m * 2 + m * m + 1)) :((n_dep * m * 2 + m * m + m))], nrow = 1, ncol = m, byrow = TRUE)
    # logmu           <- matrix(PD[1, ((n_dep * m * 2 + m * m + 1)) :((n_dep * m * 2 + m * m + m))], nrow = n_subj, ncol = m, byrow = TRUE)
    # logsigma2		<- tau2_d_bar <- PD[1, ((n_dep * m * 2 + m * m + m + 1)) :((n_dep * m * 2 + m * m + m + m))]


    # Start analysis --------------------------------------------
    # Run the MCMC algorithm

    if (m==2){cat("A model with only two hidden states has been specified. As a result, all transition probabilities will be fixed to 1, and no variation between individuals will be modelled.\n")}

    itime <- proc.time()[3]
    if(show_progress == TRUE){
        cat("Progress of the Bayesian mHMM algorithm:", "\n")
        pb <- utils::txtProgressBar(min = 2, max = J, style = 3)
    }

    for (iter in 2 : J){

        if (m == 2){

            # # For each subject, obtain sampled state sequence with subject individual parameters ----------
            # sample_path_state <- Dur <- vector("list", n_subj)
            #
            # for(s in 1:n_subj){
            #
            #     # Idea: pre-compute emissions likelihood to pass to FBalgC() here:
            #
            #     # Run forward backward algorithm in C++, using the runlength distribution d for each state ================================================================
            #     d 	<- get.d.lognorm(run.p = list(logmu = logmu[s,], logsd = sqrt(logsigma2)), Mx = subj_data[[s]]$Mx, m = m)
            #
            #     delta[[s]] <- get_delta(gamma[[s]], m)
            #
            #     allprobs <- get_all1(x = subj_data[[s]]$y, emiss = emiss[[s]], n_dep = n_dep, data_distr = "continuous")
            #
            #     FB	<- mult_ed_fb_cpp(
            #         # y2 = subj_data[[s]]$y,
            #         m = m,
            #         n = subj_data[[s]]$n,
            #         allprobs = t(allprobs),
            #         Mx = subj_data[[s]]$Mx,
            #         Mx2 = subj_data[[s]]$Mx2,
            #         gamma = gamma[[s]],
            #         d = d,
            #         S2 = subj_data[[s]]$switch2,
            #         S = subj_data[[s]]$switch,
            #         delta = delta[[s]]
            #     )
            #
            #     B_star				<- FB[[2]]
            #     Occupancy			<- FB[[3]]
            #     N 					<- FB[[1]]
            #     PD_subj[[s]][iter-1, m*n_dep*2 + m*m + m*2 + 1] <- llk <- sum(N)	# adjust index; we may need to do log-sum-ex
            #
            #     # Using the outcomes of the forward backward algorithm, sample the state sequence ==========================================================================
            #     trans[[s]]				                <- vector("list", m)
            #     sample_path_state[[s]][1] 	            <- sample(1:m, 1, prob = delta[[s]] * exp(B_star[,1]))
            #     Dur[[s]][1] 			                <- sample(1:subj_data[[s]]$Mx, 1, prob = (Occupancy[[1]][sample_path_state[[s]][1],] / exp(B_star[sample_path_state[[s]][1],1])))
            #     sample_path[[s]][1:Dur[[s]][1], iter]   <- sample_path_state[[s]][1]
            #
            #     t <- 1
            #     while(sum(Dur[[s]]) < subj_data[[s]]$n){
            #         t 						                                        <- t + 1
            #         Mx.l 					                                        <- min(subj_data[[s]]$Mx, subj_data[[s]]$n-sum(Dur[[s]]))
            #         sample_path_state[[s]][t] 	                                    <- sample(1:m, 1, prob = gamma[[s]][sample_path_state[[s]][t-1],] * exp(B_star[,sum(Dur[[s]])+1]))
            #         trans[[s]][[sample_path_state[[s]][t-1]]]                       <- c(trans[[s]][[sample_path_state[[s]][t-1]]], sample_path_state[[s]][t])
            #         Dur[[s]][t]			                                            <- sample(1:Mx.l, 1, prob = (Occupancy[[sum(Dur[[s]])+1]][sample_path_state[[s]][t],] / exp(B_star[sample_path_state[[s]][t], sum(Dur[[s]])+1])))
            #         sample_path[[s]][sum(Dur[[s]][1:t-1],1):sum(Dur[[s]]), iter]    <- sample_path_state[[s]][t]
            #     }
            #
            #     n.Dur[s] <- length(Dur[[s]])
            #     for (i in 1:m){
            #         # trans[[s]][[i]]         <- c(trans[[s]][[i]], 1:m) # to avoid errors, check if we can drop
            #         for (q in 1:n_dep) {
            #             # cond_y[[s]][[i]]    <- subj_data[[s]]$y[sample_path[[s]][, iter] == i, q]
            #             if(iter == 2){
            #                 cond_y[[s]][[i]][[q]] <- c(subj_data[[s]]$y[sample_path[[s]][, iter] == i, q][!is.na(subj_data[[s]]$y[sample_path[[s]][, iter] == i, q])],emiss_mu0[[q]][1,i])
            #             } else {
            #                 cond_y[[s]][[i]][[q]] <- c(subj_data[[s]]$y[sample_path[[s]][, iter] == i, q][!is.na(subj_data[[s]]$y[sample_path[[s]][, iter] == i, q])],emiss_c_mu_bar[[i]][[q]][1])
            #             }
            #         }
            #
            #     }
            # }
            #
            # The remainder of the mcmc algorithm is state specific
            # for(i in 1:m){
            #
            #     # Sample populaton values for gamma and conditional probabilities using Gibbs sampler -----------
            #     # gamma_mu0_n and gamma_mu_int_bar are matrices, with the number of rows equal to the number of covariates, and ncol equal to number of intercepts estimated
            #     gamma_V_int[[i]]      <- 0
            #     gamma_mu_int_bar[[i]] <- matrix(matrix(c(Inf,
            #                                              -Inf), nrow = 2, byrow = TRUE)[i,], nrow = 1)
            #     gamma_mu_prob_bar[[i]] 	<- as.vector(matrix(c(0,1,
            #                                                   1,0), nrow = 2, byrow = TRUE)[i,])
            #
            #     # sample population mean (and regression parameters if covariates) of the Normal emission distribution, and it's variance (so the variance between the subject specific means)
            #     # note: the posterior is thus one of a Bayesian linear regression because of the optional regression parameters
            #     for(q in 1:n_dep){
            #         emiss_mu0_n                    <- solve(t(xx[[1 + q]]) %*% xx[[1 + q]] + emiss_K0[[q]]) %*% (t(xx[[1 + q]]) %*% emiss_c_mu[[i]][[q]] + emiss_K0[[q]] %*% emiss_mu0[[q]][,i])
            #         emiss_a_mu_n                   <- (emiss_K0[[q]] + n_subj) / 2
            #         emiss_b_mu_n                   <- (emiss_nu[[q]] * emiss_V[[q]][i]) / 2 + (t(emiss_c_mu[[i]][[q]]) %*% emiss_c_mu[[i]][[q]] +
            #                                                                                        t(emiss_mu0[[q]][,i]) %*% emiss_K0[[q]] %*% emiss_mu0[[q]][,i] -
            #                                                                                        t(emiss_mu0_n) %*% (t(xx[[1 + q]]) %*% xx[[1 + q]] + emiss_K0[[q]]) %*% emiss_mu0_n) / 2
            #         emiss_V_mu[[i]][[q]]       <- solve(stats::rgamma(1, shape = emiss_a_mu_n, rate = emiss_b_mu_n))
            #         if(all(dim(emiss_V_mu[[i]][[q]]) == c(1,1))){
            #             emiss_c_mu_bar[[i]][[q]]	  <- emiss_mu0_n + rnorm(1 + nx[1 + q] - 1, mean = 0, sd = sqrt(diag(as.numeric(emiss_V_mu[[i]][[q]]) * solve(t(xx[[1 + q]]) %*% xx[[1 + q]] + emiss_K0[[q]]))))
            #         } else {
            #             emiss_c_mu_bar[[i]][[q]]	  <- emiss_mu0_n + rnorm(1 + nx[1 + q] - 1, mean = 0, sd = sqrt(diag(emiss_V_mu[[i]][[q]] * solve(t(xx[[1 + q]]) %*% xx[[1 + q]] + emiss_K0[[q]]))))
            #         }
            #     }
            #
            #     # Sample subject values  -----------
            #     for (s in 1:n_subj){
            #
            #         gamma[[s]][i,]  	    <- PD_subj[[s]][iter, c((n_dep * 2 * m + 1 + (i - 1) * m):(n_dep * 2 * m + (i - 1) * m + m))] <- as.vector(matrix(c(0,1,
            #                                                                                                                                                  1,0), nrow = 2, byrow = TRUE)[i,])
            #         gamma_naccept[s, i]		<- gamma_naccept[s, i] + 1
            #         gamma_c_int[[i]][s,]	<- matrix(c(Inf,
            #                                          -Inf), nrow = 2, byrow = TRUE)[i,]
            #         # gamma_int_subj[[s]][iter, c((1 + (i - 1) * (m - 2)):((m - 2) + (i - 1) * (m - 2)))[-i]] <- gamma_c_int[[i]][s,] # CHECK
            #         gamma_int_subj[[s]][iter, c((1 + (i - 1) * (m - 2)):((m - 2) + (i - 1) * (m - 2)))] <- gamma_c_int[[i]][s,]
            #
            #         if(i == m){
            #             delta[[s]] 		<- solve(t(diag(m) - gamma[[s]] + 1), rep(1, m))
            #         }
            #     }
            #     # Sample subject values for normal emission distribution using Gibbs sampler   ---------
            #
            #     # population level, conditional probabilities, seperate for each dependent variable
            #     for(q in 1:n_dep){
            #         for (s in 1:n_subj){
            #             ss_subj[s] <- t(matrix(cond_y[[s]][[i]][[q]] - emiss_c_mu[[i]][[q]][s,1], nrow = 1) %*%
            #                                 matrix(cond_y[[s]][[i]][[q]] - emiss_c_mu[[i]][[q]][s,1], ncol = 1))
            #             n_cond_y[s]       <- length(cond_y[[s]][[i]][[q]])
            #         }
            #         emiss_a_resvar_n <- sum(n_cond_y) / 2 + emiss_a0[[q]][i]
            #         emiss_b_resvar_n <- (sum(ss_subj) + 2 * emiss_b0[[q]][i]) / 2
            #         emiss_c_V[[i]][[q]] <- emiss_var_bar[[q]][iter, i] <- solve(stats::rgamma(1, shape = emiss_a_resvar_n, rate = emiss_b_resvar_n))
            #     }
            #
            #     ### sampling subject specific means for the emission distributions, assuming known mean and var, see Lynch p. 244
            #     for(q in 1:n_dep){
            #         emiss_c_V_subj    <- (emiss_V_mu[[i]][[q]] * emiss_c_V[[i]][[q]]) / (2 * emiss_V_mu[[i]][[q]] + emiss_c_V[[i]][[q]])
            #         for (s in 1:n_subj){
            #             emiss_mu0_subj_n  <- (emiss_V_mu[[i]][[q]] * sum(cond_y[[s]][[i]][[q]]) +  emiss_c_V[[i]][[q]] * c(t(emiss_c_mu_bar[[i]][[q]]) %*% xx[[q+1]][s,])) /
            #                 (n_cond_y[s] * emiss_V_mu[[i]][[q]] + emiss_c_V[[i]][[q]])
            #             emiss[[s]][[q]][i,1] <- PD_subj[[s]][iter, ((q - 1) * m + i)] <- emiss_c_mu[[i]][[q]][s,1] <- rnorm(1, emiss_mu0_subj_n, sqrt(emiss_c_V_subj))
            #             emiss[[s]][[q]][i,2] <- PD_subj[[s]][iter, (n_dep * m + (q - 1) * m + i)] <- emiss_c_V[[i]][[q]]
            #         }
            #     }
            #
            #
            #     #################
            #     # Obtain hierarchical and mouse specific parameters for duration distribuiton using gibbs sampler ====================================
            #     #################
            #
            #     #draw logmu's
            #     for(s in 1:n_subj){
            #         tau2 		<- 1/ ((1/tau2_d_bar[i]) + (1/logsigma2[i]) * sum(sample_path_state[[s]][-n.Dur[s]] == i))
            #         logmu[s,i]	<- rnorm(1,
            #                             mean = tau2 * ((1/tau2_d_bar[i]) * mu_d_bar[i] + (1/logsigma2[i]) * sum(log(Dur[[s]][-n.Dur[s]][sample_path_state[[s]][-n.Dur[s]] == i]))),
            #                             sd = sqrt(tau2))
            #         PD_subj[[s]][iter, n_dep*m*2 + m*m + i]			<- logmu[s, i] # adjust index
            #     }
            #
            #     # draw mu_d_bar
            #     s2				<- 1 / ((1/s2_0[i]) + (1/tau2_d_bar[i]) * n_subj)
            #     mu_d_bar[i] 	<- rnorm(1,
            #                           mean = s2 * ((1/s2_0[i]) * d_mu0[i] + (1/tau2_d_bar[i]) * sum(logmu[,i])),
            #                           sd = sqrt(s2))
            #
            #     # draw tau2_d_bar
            #     a1 <- alpha.tau20[i] + n_subj/2
            #     b1 <- beta.tau20[i] + (1/2) * (sum((logmu[,i] - mu_d_bar[i])^2))
            #     tau2_d_bar[i] <- 1 / rgamma(1, shape = a1, rate = b1)
            #
            #     #draw logsigma
            #     ss <- numeric(1)
            #     n.ss <- numeric(1)
            #     for (s in 1:n_subj){
            #         ss <- ss + sum((log(Dur[[s]][-n.Dur[s]][sample_path_state[[s]][-n.Dur[s]] == i]) - logmu[s,i])^2)
            #         n.ss <- n.ss + sum(sample_path_state[[s]][-n.Dur[s]] == i)
            #     }
            #     c1 <- alpha.sigma20[i] + n.ss/2
            #     d1 <- beta.sigma20[i] + (1/2) * ss
            #     logsigma2[i] <- 1 / rgamma(1, shape = c1, rate = d1)
            #     for (s in 1:n_subj) {
            #         PD_subj[[s]][iter, n_dep*m*2 + m*m + m + i] <- logsigma2[i]
            #     }
            #
            # }

        } else if (m >= 3){

            iter = 2

            # For each subject, obtain sampled state sequence with subject individual parameters ----------
            sample_path_state <- Dur <- vector("list", n_subj)

            for(s in 1:n_subj){

                # Idea: pre-compute emissions likelihood to pass to FBalgC() here:

                # Run forward backward algorithm in C++, using the runlength distribution d for each state ================================================================
                # d 	<- get.d.pois(run.p = list(lambda = dwell_c_mu[[i]][s,1]), Mx = subj_data[[s]]$Mx, m = m) # Check
                d 	<- get.d.shiftpois(run.p = list(lambda = t(dwell[[s]]), shift = shift), Mx = subj_data[[s]]$Mx, m = m)

                delta[[s]] <- get_delta(gamma[[s]], m)

                allprobs <- get_all1(x = subj_data[[s]]$y, emiss = emiss[[s]], n_dep = n_dep, data_distr = "continuous")

                FB	<- mult_ed_fb_cpp(
                    # y2 = subj_data[[s]]$y,
                    m = m,
                    n = subj_data[[s]]$n,
                    allprobs = t(allprobs),
                    Mx = subj_data[[s]]$Mx,
                    Mx2 = subj_data[[s]]$Mx2,
                    gamma = gamma[[s]],
                    d = d,
                    # S2 = subj_data[[s]]$switch2,
                    # S = subj_data[[s]]$switch,
                    S2 = rep(1,500),
                    S = rep(1,500),
                    delta = delta[[s]]
                )

                B_star				<- FB[[2]]
                Occupancy			<- FB[[3]]
                N 					<- FB[[1]]
                PD_subj[[s]][iter-1, m*n_dep*2 + m*m + m + 1] <- llk <- sum(N)	# adjust index; we may need to do log-sum-ex

                # Using the outcomes of the forward backward algorithm, sample the state sequence ==========================================================================
                trans[[s]]				                <- vector("list", m)
                sample_path_state[[s]][1] 	            <- sample(1:m, 1, prob = delta[[s]] * exp(B_star[,1]))
                Dur[[s]][1] 			                <- sample(1:subj_data[[s]]$Mx, 1, prob = (Occupancy[[1]][sample_path_state[[s]][1],] / exp(B_star[sample_path_state[[s]][1],1])))
                sample_path[[s]][1:Dur[[s]][1], iter]   <- sample_path_state[[s]][1]

                t <- 1
                while(sum(Dur[[s]]) < subj_data[[s]]$n){
                    t 						                                        <- t + 1
                    Mx.l 					                                        <- min(subj_data[[s]]$Mx, subj_data[[s]]$n-sum(Dur[[s]]))
                    sample_path_state[[s]][t] 	                                    <- sample(1:m, 1, prob = gamma[[s]][sample_path_state[[s]][t-1],] * exp(B_star[,sum(Dur[[s]])+1]))
                    trans[[s]][[sample_path_state[[s]][t-1]]]                       <- c(trans[[s]][[sample_path_state[[s]][t-1]]], sample_path_state[[s]][t])
                    Dur[[s]][t]			                                            <- sample(1:Mx.l, 1, prob = (Occupancy[[sum(Dur[[s]])+1]][sample_path_state[[s]][t],] / exp(B_star[sample_path_state[[s]][t], sum(Dur[[s]])+1])))
                    sample_path[[s]][sum(Dur[[s]][1:t-1],1):sum(Dur[[s]]), iter]    <- sample_path_state[[s]][t]
                }

                n.Dur[s] <- length(Dur[[s]])
                for (i in 1:m){
                    # trans[[s]][[i]]         <- c(trans[[s]][[i]], 1:m) # to avoid errors, check if we can drop
                    for (q in 1:n_dep) {
                        # cond_y[[s]][[i]]    <- subj_data[[s]]$y[sample_path[[s]][, iter] == i, q]
                        if(iter == 2){
                            cond_y[[s]][[i]][[q]] <- c(subj_data[[s]]$y[sample_path[[s]][, iter] == i, q][!is.na(subj_data[[s]]$y[sample_path[[s]][, iter] == i, q])],emiss_mu0[[q]][1,i])
                        } else {
                            cond_y[[s]][[i]][[q]] <- c(subj_data[[s]]$y[sample_path[[s]][, iter] == i, q][!is.na(subj_data[[s]]$y[sample_path[[s]][, iter] == i, q])],emiss_c_mu_bar[[i]][[q]][1])
                        }
                    }

                }
            }

            # The remainder of the mcmc algorithm is state specific
            for(i in 1:m){

                # Obtain MLE of the covariance matrices and log likelihood of gamma and emiss at subject and population level -----------------
                # used to scale the propasal distribution of the RW Metropolis sampler

                # population level, transition matrix ============================================================ Check number of categories in gamma_int_mle0
                trans_pooled			<- factor(c(unlist(sapply(trans, "[[", i)), c(1:m)[-i]))
                trans_pooled            <- factor(trans_pooled, levels = paste(c(1:m)[-i]))
                trans_pooled            <- as.numeric(trans_pooled)
                gamma_mle_pooled		<- optim(gamma_int_mle0, llmnl_int, Obs = trans_pooled,
                                           method = "BFGS", hessian = FALSE,
                                           control = list(fnscale = -1))

                gamma_int_mle_pooled[[i]]   <- gamma_mle_pooled$par
                gamma_pooled_ll[[i]]        <- gamma_mle_pooled$value

                # subject level
                for (s in 1:n_subj){
                    wgt 		        <- subj_data[[s]]$n / n_total

                    trans_subj          <- factor(c(trans[[s]][[i]], c(1:m)[-i]))
                    trans_subj          <- factor(trans_subj, levels = paste(c(1:m)[-i]))
                    trans_subj          <- as.numeric(trans_subj)

                    # subject level, transition matrix ============================================================ Check number of categories in gamma_int_mle_pooled
                    gamma_out			<- optim(gamma_int_mle_pooled[[i]], llmnl_int_frac, Obs = trans_subj,
                                         # n_cat = m,
                                         pooled_likel = gamma_pooled_ll[[i]], w = gamma_w, wgt = wgt,
                                         method="BFGS", hessian = FALSE, control = list(fnscale = -1))
                    if(gamma_out$convergence == 0){
                        subj_data[[s]]$gamma_converge[i] <- 1
                        subj_data[[s]]$gamma_int_mle[i,] <- gamma_out$par
                        subj_data[[s]]$gamma_mhess[(1 + (i - 1) * (m - 2)):((m - 2) + (i - 1) * (m - 2)), ]	<-
                            mnlHess_int(int = gamma_out$par, Obs = c(trans[[s]][[i]], c(1:m)), n_cat =  (m-1) )
                    } else {
                        subj_data[[s]]$gamma_converge[i] <- 0
                        subj_data[[s]]$gamma_int_mle[i,] <- rep(0, m - 2)
                        subj_data[[s]]$gamma_mhess[(1 + (i - 1) * (m - 2)):((m - 2) + (i - 1) * (m - 2)), ]	<- diag(m-2)
                    }

                    # if this is first iteration, use MLE for current values RW metropolis sampler
                    if (iter == 2){
                        gamma_c_int[[i]][s,]		<- gamma_out$par
                    }
                }

                # Sample pouplaton values for gamma and conditional probabilities using Gibbs sampler -----------
                # gamma_mu0_n and gamma_mu_int_bar are matrices, with the number of rows equal to the number of covariates, and ncol equal to number of intercepts estimated
                gamma_mu0_n           <- solve(t(xx[[1]]) %*% xx[[1]] + gamma_K0)  %*% (t(xx[[1]]) %*% gamma_c_int[[i]] + gamma_K0 %*% gamma_mu0[[i]])
                gamma_V_n             <- gamma_V + t(gamma_c_int[[i]] - xx[[1]] %*% gamma_mu0_n) %*% (gamma_c_int[[i]] - xx[[1]] %*% gamma_mu0_n) + t(gamma_mu0_n - gamma_mu0[[i]]) %*% gamma_K0 %*% (gamma_mu0_n - gamma_mu0[[i]])
                gamma_V_int[[i]]      <- solve(rwish(S = solve(gamma_V_n), v = gamma_nu + n_subj))
                gamma_mu_int_bar[[i]] <- gamma_mu0_n + solve(chol(t(xx[[1]]) %*% xx[[1]] + gamma_K0)) %*% matrix(rnorm((m - 2) * nx[1]), nrow = nx[1]) %*% t(solve(chol(solve(gamma_V_int[[i]]))))
                gamma_exp_int				  <- matrix(exp(c(0, gamma_mu_int_bar[[i]][1,] )), nrow  = 1)
                gamma_mu_prob_bar[[i]] 	<- gamma_exp_int / as.vector(gamma_exp_int %*% c(rep(1,(m-1))))

                # sample population mean (and regression parameters if covariates) of the Normal emission distribution, and it's variance (so the variance between the subject specific means)
                # note: the posterior is thus one of a Bayesian linear regression because of the optional regression parameters
                for(q in 1:n_dep){
                    emiss_mu0_n                    <- solve(t(xx[[1 + q]]) %*% xx[[1 + q]] + emiss_K0[[q]]) %*% (t(xx[[1 + q]]) %*% emiss_c_mu[[i]][[q]] + emiss_K0[[q]] %*% emiss_mu0[[q]][,i])
                    emiss_a_mu_n                   <- (emiss_K0[[q]] + n_subj) / 2
                    emiss_b_mu_n                   <- (emiss_nu[[q]] * emiss_V[[q]][i]) / 2 + (t(emiss_c_mu[[i]][[q]]) %*% emiss_c_mu[[i]][[q]] +
                                                                                                   t(emiss_mu0[[q]][,i]) %*% emiss_K0[[q]] %*% emiss_mu0[[q]][,i] -
                                                                                                   t(emiss_mu0_n) %*% (t(xx[[1 + q]]) %*% xx[[1 + q]] + emiss_K0[[q]]) %*% emiss_mu0_n) / 2
                    emiss_V_mu[[i]][[q]]       <- solve(stats::rgamma(1, shape = emiss_a_mu_n, rate = emiss_b_mu_n))
                    if(all(dim(emiss_V_mu[[i]][[q]]) == c(1,1))){
                        emiss_c_mu_bar[[i]][[q]]	  <- emiss_mu0_n + rnorm(1 + nx[1 + q] - 1, mean = 0, sd = sqrt(diag(as.numeric(emiss_V_mu[[i]][[q]]) * solve(t(xx[[1 + q]]) %*% xx[[1 + q]] + emiss_K0[[q]]))))
                    } else {
                        emiss_c_mu_bar[[i]][[q]]	  <- emiss_mu0_n + rnorm(1 + nx[1 + q] - 1, mean = 0, sd = sqrt(diag(emiss_V_mu[[i]][[q]] * solve(t(xx[[1 + q]]) %*% xx[[1 + q]] + emiss_K0[[q]]))))
                    }
                }

                # Sample subject values  -----------
                for (s in 1:n_subj){

                    trans_subj          <- factor(c(trans[[s]][[i]]))
                    trans_subj          <- factor(trans_subj, levels = paste(c(1:m)[-i]))
                    trans_subj          <- as.numeric(trans_subj)

                    # Sample subject values for gamma using RW Metropolis sampler   ---------
                    gamma_candcov_comb 		<- chol2inv(chol(subj_data[[s]]$gamma_mhess[(1 + (i - 1) * (m - 2)):((m - 2) + (i - 1) * (m - 2)), ] + chol2inv(chol(gamma_V_int[[i]]))))
                    gamma_RWout				<- mnl_RW_once(int1 = gamma_c_int[[i]][s,],
                                                  Obs = trans_subj,
                                                  n_cat = m-1,
                                                  mu_int_bar1 = c(t(gamma_mu_int_bar[[i]]) %*% xx[[1]][s,]),
                                                  V_int1 = gamma_V_int[[i]],
                                                  scalar = gamma_scalar,
                                                  candcov1 = gamma_candcov_comb)
                    gamma[[s]][i,-i]  	    <- (gamma_RWout$prob + .00001) / (1 + .00001 * m)
                    PD_subj[[s]][iter, c((n_dep * 2 * m + 1 + (i - 1) * m):(n_dep * 2 * m + (i - 1) * m + m))[-i]] <- gamma_RWout$prob
                    PD_subj[[s]][iter, c((n_dep * 2 * m + 1 + (i - 1) * m):(n_dep * 2 * m + (i - 1) * m + m))[i]] <- 0
                    gamma_naccept[s, i]		<- gamma_naccept[s, i] + gamma_RWout$accept
                    gamma_c_int[[i]][s,]	<- gamma_RWout$draw_int
                    # gamma_int_subj[[s]][iter, c((1 + (i - 1) * (m - 2)):((m - 2) + (i - 1) * (m - 2)))[-i]] <- gamma_c_int[[i]][s,] # CHECK
                    gamma_int_subj[[s]][iter, c((1 + (i - 1) * (m - 2)):((m - 2) + (i - 1) * (m - 2)))] <- gamma_c_int[[i]][s,]

                    if(i == m){
                        delta[[s]] 		<- solve(t(diag(m) - gamma[[s]] + 1), rep(1, m))
                    }
                }
                # Sample subject values for normal emission distribution using Gibbs sampler   ---------

                # population level, conditional probabilities, separate for each dependent variable
                for(q in 1:n_dep){
                    for (s in 1:n_subj){
                        ss_subj[s] <- t(matrix(cond_y[[s]][[i]][[q]] - emiss_c_mu[[i]][[q]][s,1], nrow = 1) %*%
                                            matrix(cond_y[[s]][[i]][[q]] - emiss_c_mu[[i]][[q]][s,1], ncol = 1))
                        n_cond_y[s]       <- length(cond_y[[s]][[i]][[q]])
                    }
                    emiss_a_resvar_n <- sum(n_cond_y) / 2 + emiss_a0[[q]][i]
                    emiss_b_resvar_n <- (sum(ss_subj) + 2 * emiss_b0[[q]][i]) / 2
                    emiss_c_V[[i]][[q]] <- emiss_var_bar[[q]][iter, i] <- solve(stats::rgamma(1, shape = emiss_a_resvar_n, rate = emiss_b_resvar_n))
                }

                ### sampling subject specific means for the emission distributions, assuming known mean and var, see Lynch p. 244
                for(q in 1:n_dep){
                    emiss_c_V_subj    <- (emiss_V_mu[[i]][[q]] * emiss_c_V[[i]][[q]]) / (2 * emiss_V_mu[[i]][[q]] + emiss_c_V[[i]][[q]])
                    for (s in 1:n_subj){
                        emiss_mu0_subj_n  <- (emiss_V_mu[[i]][[q]] * sum(cond_y[[s]][[i]][[q]]) +  emiss_c_V[[i]][[q]] * c(t(emiss_c_mu_bar[[i]][[q]]) %*% xx[[q+1]][s,])) /
                            (n_cond_y[s] * emiss_V_mu[[i]][[q]] + emiss_c_V[[i]][[q]])
                        emiss[[s]][[q]][i,1] <- PD_subj[[s]][iter, ((q - 1) * m + i)] <- emiss_c_mu[[i]][[q]][s,1] <- rnorm(1, emiss_mu0_subj_n, sqrt(emiss_c_V_subj))
                        emiss[[s]][[q]][i,2] <- PD_subj[[s]][iter, (n_dep * m + (q - 1) * m + i)] <- emiss_c_V[[i]][[q]]
                    }
                }


                #################
                # Obtain hierarchical and individual specific parameters for duration distribution ====================================
                #################

                # population level, conditional probabilities, separate for each dependent variable
                Dur_pooled					      <- numeric()
                for(s in 1:n_subj){
                    Dur_pooled             <- c(Dur_pooled, Dur[[s]][-n.Dur[s]][sample_path_state[[s]][-n.Dur[s]] == i])
                }
                if(iter <= 2) {
                    dwell_c_mu_bar[[i]] <- log(start_val[[n_dep+2]][i,1])
                }
                dwell_mle_pooled[[i]]  <- mean(c(Dur_pooled, round(exp(dwell_c_mu_bar[[i]][1]),0)))
                if(dwell_mle_pooled[[i]] == 0){
                    dwell_mle_pooled[[i]] <- 1
                }
                dwell_pooled_ll[[i]]	<- llshiftpois(lambda = dwell_mle_pooled[[i]], Obs = c(Dur_pooled, round(exp(dwell_c_mu_bar[[i]][1]),0)), shift = shift) # Check Dur cond_y pooled?

                # subject level, conditional probabilities, seperate for each dependent variable
                for(s in 1:n_subj){
                    dwell_out	<- optim(log(dwell_mle_pooled[[i]]), llshiftpois_frac_log, Obs = c(Dur[[s]][-n.Dur[s]][sample_path_state[[s]][-n.Dur[s]] == i],
                                                                                                 round(exp(dwell_c_mu_bar[[i]][1]),0)),
                                       pooled_likel = dwell_pooled_ll[[i]],
                                       w = dwell_w, wgt = wgt, shift = shift,
                                       method = "BFGS",
                                       hessian = TRUE,
                                       control = list(fnscale = -1))
                    if(dwell_out$convergence == 0){
                        subj_data[[s]]$dwell_converge[i]   <- 1
                        subj_data[[s]]$dwell_mle[i]        <- exp(dwell_out$par)
                        subj_data[[s]]$dwell_mhess[i]      <- as.numeric(-dwell_out$hessian)
                    } else {
                        subj_data[[s]]$dwell_converge[i]   <- 0
                        subj_data[[s]]$dwell_mle[i]        <- 0
                        subj_data[[s]]$dwell_mhess[i]	   <- 0
                    }

                }

                # sample population mean (and regression parameters if covariates) of the Normal emission distribution, and it's variance (so the variance between the subject specific means)
                # note: the posterior is thus one of a Bayesian linear regression because of the optional regression parameters
                dwell_mu0_n                    <- solve(t(xx[[2 + n_dep]]) %*% xx[[2 + n_dep]] + dwell_K0) %*% (t(xx[[2 + n_dep]]) %*% log(dwell_c_mu[[i]]) + dwell_K0 %*% dwell_mu0[,i]) # CHECK THIS
                dwell_a_mu_n                   <- (dwell_K0 + n_subj) / 2
                dwell_b_mu_n                   <- (dwell_nu * dwell_V[i]) / 2 + (t(log(dwell_c_mu[[i]])) %*% log(dwell_c_mu[[i]]) +
                                                                                     t(dwell_mu0[,i]) %*% dwell_K0 %*% dwell_mu0[,i] -
                                                                                     t(dwell_mu0_n) %*% (t(xx[[2 + n_dep]]) %*% xx[[2 + n_dep]] + dwell_K0) %*% dwell_mu0_n) / 2
                dwell_V_mu[[i]]       <- solve(rgamma(1, shape = dwell_a_mu_n, rate = dwell_b_mu_n))
                if(all(dim(dwell_V_mu[[i]]) == c(1,1))){ # CHECK THIS
                    dwell_c_mu_bar[[i]]	  <- dwell_mu0_n + rnorm(1 + nx[2 + n_dep] - 1, mean = 0, sd = sqrt(diag(as.numeric(dwell_V_mu[[i]]) * solve(t(xx[[2 + n_dep]]) %*% xx[[2 + n_dep]] + dwell_K0))))
                } else {
                    dwell_c_mu_bar[[i]]	  <- dwell_mu0_n + rnorm(1 + nx[2 + n_dep] - 1, mean = 0, sd = sqrt(diag(dwell_V_mu[[i]] * solve(t(xx[[2 + n_dep]]) %*% xx[[2 + n_dep]] + dwell_K0))))
                }

                # Sample subject values for normal emission distribution using Gibbs sampler   ---------

                ### sampling subject specific means for the emission distributions, assuming known mean and var, see Lynch p. 244
                for (s in 1:n_subj){

                    # Sample subject specific lambda with RW-MH
                    dwell_mu0_subj_bar <- c(t(dwell_c_mu_bar[[i]]) %*% xx[[2 + n_dep]][s,]) # Update: add xx for dwell time?
                    dwell_candcov_comb <- (subj_data[[s]]$dwell_mhess[i] + dwell_V_mu[[i]]^-1)^-1
                    dwell_rw_out <- shiftpoisLN_RW_once(lambda = dwell_c_mu[[i]][s,1],
                                                        # Obs = c(cond_y[[s]][[i]], round(exp(dwell_mu0_subj_bar), 0)), # Dur instead of cond_y
                                                        Obs = c(Dur[[s]][-n.Dur[s]][sample_path_state[[s]][-n.Dur[s]] == i], round(exp(dwell_mu0_subj_bar), 0)),
                                                        mu_bar1 = dwell_mu0_subj_bar,
                                                        V_1 = sqrt(dwell_V_mu[[i]]),
                                                        scalar = dwell_scalar,
                                                        candcov1 = dwell_candcov_comb,
                                                        shift = shift)

                    dwell[[s]][i,1]        <- PD_subj[[s]][iter, n_dep*m*2 + m*m + i] <- dwell_c_mu[[i]][s,1] <- dwell_rw_out$draw_lambda # Check PD dwell time indexation
                    dwell_naccept[s, i]    <- dwell_naccept[s, i] + dwell_rw_out$accept

                }

            }

        }


        # End of MCMC iteration, save output values --------
        gamma_int_bar[iter, ]				   	<- unlist(lapply(gamma_mu_int_bar, "[",1,))
        gamma_V_int_bar[iter, ] <- unlist(lapply(gamma_V_int, function(e) as.vector(t(e))))
        if(nx[1] > 1){
            gamma_cov_bar[iter, ]      	<- unlist(lapply(gamma_mu_int_bar, "[",-1,))
        }
        if(m == 2) {
            gamma_prob_bar[iter,]                                   <- unlist(gamma_mu_prob_bar)
        } else {
            gamma_prob_bar[iter,(1:m^2)[-seq(1,m*m, m+1)]]			<- unlist(gamma_mu_prob_bar)
            gamma_prob_bar[iter,seq(1,m*m, m+1)]			        <- 0
        }
        for(q in 1:n_dep){
            emiss_mu_bar[[q]][iter, ]	<- as.vector(unlist(lapply(
                lapply(emiss_c_mu_bar, "[[", q), "[",1,)
            ))
            if(nx[1+q] > 1){
                emiss_cov_bar[[q]][iter, ]  <- as.vector(unlist(lapply(
                    lapply(emiss_c_mu_bar, "[[", q), "[",-1,)
                ))
            }
            emiss_varmu_bar[[q]][iter,]	<- as.vector(unlist(sapply(emiss_V_mu, "[[", q)))
        }

        # dwell time parameters
        dwell_mu_bar[iter, ]	<- as.vector(unlist(lapply(dwell_c_mu_bar, "[",1,)))
        if(nx[2+n_dep] > 1){
            dwell_cov_bar[iter, ]  <- as.vector(unlist(lapply(dwell_c_mu_bar, "[",-1,)))
        }
        # dwell_varmu_bar[iter,]	<- as.vector(unlist(sapply(dwell_V_mu, "[", )))
        dwell_varmu_bar[iter,]	<- unlist(dwell_V_mu)



        # # dwell time parameters
        # # duration_bar[iter,] <- c(mu_d_bar, logsigma2, tau2_d_bar)
        # dwell_mu_bar[iter, ] <- mu_d_bar # IMPROVE memory use
        # dwell_varmu_bar[iter, ] <- tau2_d_bar
        # dwell_var_bar[iter, ] <- logsigma2

        if(show_progress == TRUE){
            utils::setTxtProgressBar(pb, iter)
        }
    }
    if(show_progress == TRUE){
        close(pb)
    }

    # End of function, return output values --------
    ctime = proc.time()[3]
    message(paste("Total time elapsed (hh:mm:ss):", hms(ctime-itime)))
    if(return_path == TRUE){
        out <- list(input = list(m = m, n_dep = n_dep, J = J,
                                 burn_in = burn_in, n_subj = n_subj, n_vary = n_vary, dep_labels = dep_labels),
                    PD_subj = PD_subj,
                    gamma_prob_bar = gamma_prob_bar,
                    gamma_int_bar = gamma_int_bar,
                    gamma_cov_bar = gamma_cov_bar,
                    gamma_V_int_bar = gamma_V_int_bar,
                    gamma_int_subj = gamma_int_subj,
                    gamma_naccept = gamma_naccept,
                    emiss_mu_bar = emiss_mu_bar, emiss_varmu_bar = emiss_varmu_bar, emiss_cov_bar = emiss_cov_bar,
                    emiss_var_bar = emiss_var_bar,
                    dwell_mu_bar = dwell_mu_bar,
                    dwell_varmu_bar = dwell_varmu_bar,
                    dwell_cov_bar = dwell_cov_bar,
                    dwell_naccept = dwell_naccept,
                    sample_path = sample_path)
    } else {
        out <- list(input = list(m = m, n_dep = n_dep, J = J,
                                 burn_in = burn_in, n_subj = n_subj, n_vary = n_vary, dep_labels = dep_labels),
                    PD_subj = PD_subj,
                    gamma_prob_bar = gamma_prob_bar,
                    gamma_int_bar = gamma_int_bar,
                    gamma_cov_bar = gamma_cov_bar,
                    gamma_V_int_bar = gamma_V_int_bar,
                    gamma_int_subj = gamma_int_subj,
                    gamma_naccept = gamma_naccept,
                    emiss_mu_bar = emiss_mu_bar, emiss_varmu_bar = emiss_varmu_bar, emiss_cov_bar = emiss_cov_bar,
                    emiss_var_bar = emiss_var_bar,
                    dwell_mu_bar = dwell_mu_bar, dwell_varmu_bar = dwell_varmu_bar,
                    dwell_cov_bar = dwell_cov_bar,
                    dwell_naccept = dwell_naccept)
    }
    class(out) <- append(class(out), "medHMM_cont")
    return(out)
}



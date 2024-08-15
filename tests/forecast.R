#' Forecasts (forward predictions) using a multilevel hidden Markov model
#'
#' Ideas:
#' (1) use MAPs similarly as done with the Viterbi algorithm:
#'  +++ fast
#'  --- no measure of uncertainty on forecasts
#' (2) obtain forecasts over iterations (after burn-in):
#'  +++ uncertainty on forecasts
#'  --- slow
#'
#' To add:
#'  - option to request simulating (forecasting) observations or at least
#'    predicting the probability of observations for a number of timesteps.
#'
#' @examples
#'
#' ### Example on continuous simulated data
#' library(mHMMbayes)
#'
#' n_t     <- 500
#' n       <- 10
#' m       <- 3
#' n_dep   <- 2
#'
#' gamma   <- matrix(c(0.99, 0.005, 0.005,
#'                     0.08, 0.9, 0.02,
#'                     0.05, 0.15, 0.8), ncol = m, byrow = TRUE)
#'
#' emiss_distr <- list(matrix(c( 50, 10,
#'                               100, 10,
#'                               150, 10), nrow = m, byrow = TRUE),
#'                     matrix(c(5, 2,
#'                              10, 5,
#'                              20, 3), nrow = m, byrow = TRUE))
#'
#' data_cont <- sim_mHMM(n_t = n_t, n = n, data_distr = 'continuous', gen = list(m = m, n_dep = n_dep),
#'                       gamma = gamma, emiss_distr = emiss_distr, var_gamma = .1, var_emiss = c(.5, 0.01))
#'
#' # Specify hyper-prior for the continuous emission distribution
#' manual_prior_emiss <- prior_emiss_cont(
#'   gen = list(m = m, n_dep = n_dep),
#'   emiss_mu0 = list(matrix(c(30, 70, 170), nrow = 1),
#'                    matrix(c(7, 8, 18), nrow = 1)),
#'   emiss_K0 = list(1, 1),
#'   emiss_V =  list(rep(100, m), rep(25, m)),
#'   emiss_nu = list(1, 1),
#'   emiss_a0 = list(rep(1, m), rep(1, m)),
#'   emiss_b0 = list(rep(1, m), rep(1, m)))
#'
#' # Run the model on the simulated data:
#' # Note that for reasons of running time, J is set at a ridiculous low value.
#' # One would typically use a number of iterations J of at least 1000,
#' # and a burn_in of 200.
#' out_3st_cont_sim <- mHMM(s_data = data_cont$obs,
#'                          data_distr = 'continuous',
#'                          gen = list(m = m, n_dep = n_dep),
#'                          start_val = c(list(gamma), emiss_distr),
#'                          emiss_hyp_prior = manual_prior_emiss,
#'                          mcmc = list(J = 500, burn_in = 250))
#'
#' summary(out_3st_cont_sim)
#'
#' t(forecast_mHMM1(object = out_3st_cont_sim, s_data = data_cont$obs, forecast_steps = 50)[[3]])
#'
#' @export

forecast_medHMM1 <- function(s_data, object, forecast_steps = 1, value_range = NULL, burn_in = NULL, Mx = NULL, shift = 1, return_all = FALSE) {

  class(object) <- c(class(object),"mHMM")

  if(is.null(Mx)){
    Mx <- forecast_steps
  }

  # if (!is.mHMM(object)){
  #   stop("The input object used should be from the class mHMM, obtained by using the function mHMM.")
  # }

  input      <- object$input
  # data_distr <- input$data_distr # TBD
  data_distr <- "continuous"
  id         <- unique(s_data[,1])
  n_subj     <- length(id)
  if(length(object$PD_subj) != n_subj){
    stop("s_data used should be from the same subjects used for creating the object in mHMM.
         The number of subjects in the datasets are not the same.")
  }
  n_vary     <- table(s_data[,1])
  max_n      <- max(n_vary)
  state_seq  <- matrix(NA, ncol = n_subj, nrow = max_n)
  probs      <- vector(mode = "list", length = n_subj)
  n_dep      <- input$n_dep
  m          <- input$m
  if(is.null(burn_in)){
    burn_in  <- input$burn_in
  }
  J          <- input$J
  if (burn_in >= (J-1)){
    stop(paste("The specified burn in period should be at least 2 points smaller
               compared to the number of iterations J, J =", J))
  }

  # Get subject-specific emissions
  if(data_distr == "categorical"){
    # TBD
  } else if(data_distr == "continuous"){
    est_emiss  <- rep(list(rep(list(matrix(NA_real_,nrow = m, ncol = 2)),n_dep)), n_subj)
    for(s in 1:n_subj){
      for(q in 1:n_dep){
        est_emiss[[s]][[q]][] <- matrix(c(apply(object$PD_subj[[s]][((burn_in + 1): J),((q - 1) * m + 1):((q - 1) * m + m)], 2, median),
                                                apply(object$PD_subj[[s]][((burn_in + 1): J), (n_dep * m + (q - 1) * m + 1): (n_dep * m + (q - 1) * m + m)], 2, median)),
                                        ncol = 2, nrow = m)
      }
    }
  } else if(data_distr == "count"){
    # TBD
  }

  # Get subject-specific transitions
  est_gamma <- rep(list(matrix(NA_real_,nrow = m, ncol = m)),n_subj)
  for(s in 1:n_subj){
    est_gamma[[s]] <- matrix(c(apply(object$PD_subj[[s]][((burn_in + 1): J),(n_dep*m*2 + 1):(n_dep*m*2 + m*m)], 2, median)),
                               ncol = m, nrow = m, byrow = TRUE)
  }

  # Get subject-specific dwell times
  est_dwell <- rep(list(matrix(NA_real_,nrow = m, ncol = 1)),n_subj)
  for(s in 1:n_subj){
    est_dwell[[s]][] <- matrix(c(apply(object$PD_subj[[s]][((burn_in + 1): J),(n_dep*m*2 + m*m + 1):(n_dep*m*2 + m*m + m)], 2, median)),
                                                         ncol = 1, nrow = m)
  }


  # Obtain the forward probabilities and make forecasts

  for(s in 1:n_subj){
    emiss   <- est_emiss[[s]]
    gamma   <- est_gamma[[s]]
    dwell   <- est_dwell[[s]]
    if(data_distr == "categorical"){
      # TBD
    } else if(data_distr == "continuous"){

      d 	      <- get.d.shiftpois(run.p = list(lambda = t(dwell), shift = shift), Mx = Mx, m = m)
      delta     <- get_delta(gamma, m)
      allprobs  <- get_all1(x = rbind(as.matrix(s_data[s_data[,1] == id[s],][,-1], ncol = n_dep),
                                     matrix(NA_real_, nrow = forecast_steps, ncol = n_dep)),
                           emiss = emiss, n_dep = n_dep, data_distr = "continuous")
      if(return_all){
        probs[[s]]    <- t(apply(exp(t(mult_ed_fb_cpp(
          m = m,
          n = as.numeric(n_vary[s]+forecast_steps),
          allprobs = t(allprobs),
          Mx = Mx,
          Mx2 = rep(Mx, m),
          gamma = gamma,
          d = d,
          S2 = rep(1,n_vary[s]+forecast_steps),
          S = rep(1,n_vary[s]+forecast_steps),
          delta = delta)[[4]])),1,function(r) r/sum(r)))
      } else {
        probs[[s]]    <- t(apply(exp(t(mult_ed_fb_cpp(
          m = m,
          n = as.numeric(n_vary[s]+forecast_steps),
          allprobs = t(allprobs),
          Mx = Mx,
          Mx2 = rep(Mx, m),
          gamma = gamma,
          d = d,
          S2 = rep(1,n_vary[s]+forecast_steps),
          S = rep(1,n_vary[s]+forecast_steps),
          delta = delta)[[4]][,(n_vary[s]+1):(n_vary[s]+forecast_steps)])),1,function(r) r/sum(r)))
      }

    } else if(data_distr == "count"){
      # TBD
    }

    # Sample states
    colnames(probs[[s]]) <- paste0("pr_state_",1:m)
    if(return_all){
      probs[[s]] <- cbind(subj = s, state = apply(probs[[s]], 1, which.max), horizon = c(rep(0, n_vary[s]), 1:forecast_steps), probs[[s]])
    } else {
      probs[[s]] <- cbind(subj = s, state = apply(probs[[s]], 1, which.max), horizon = 1:forecast_steps, probs[[s]])
    }

    # # If a range of observed values is given, obtain the forecasting distribution (forecast probabilities):
    # if(!is.null(x)){
    #   for(q in 1:n_dep){
    #     # obs_probs[[q]][h,] <- emiss
    #
    #     # Forecast probs for scoring values
    #     probs[[1]][,4:6] %*% t(outer(seq(0,100,25), Y = emiss[[q]][,1], FUN = stats::dnorm, sd = emiss[[q]][,2]))
    #
    #     # Forecast prob of scoring over value
    #     probs[[1]][,4:6] %*% t(1-outer(seq(0,100,25), Y = emiss[[q]][,1], FUN = stats::pnorm, sd = emiss[[q]][,2]))
    #
    #     # Forecast prob of scoring under value
    #     probs[[1]][,4:6] %*% t(outer(seq(0,100,25), Y = emiss[[q]][,1], FUN = stats::pnorm, sd = emiss[[q]][,2]))
    #
    #   }
    #   # allprobs <- all1(x = x, emiss = emiss, n_dep = n_dep, data_distr = "categorical")
    #
    # }

  }

  # Bind into a single matrix
  forecast_probs <- do.call(rbind, lapply(probs, function(s) s))

  # Return results
  return(forecast_probs = forecast_probs)

}

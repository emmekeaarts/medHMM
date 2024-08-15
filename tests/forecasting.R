### Example on continuous simulated data
library(tidyverse)
library(mHMMbayes)

n_t     <- 500
n       <- 30
m       <- 3
n_dep   <- 2

gamma   <- matrix(c(0.99, 0.001, 0.009,
                    0.05, 0.9, 0.05,
                    0.005, 0.005, 0.99), ncol = m, byrow = TRUE)

emiss_distr <- list(matrix(c( 50, 10,
                              100, 10,
                              150, 10), nrow = m, byrow = TRUE),
                    matrix(c(5, 2,
                             10, 5,
                             20, 3), nrow = m, byrow = TRUE))

data_cont <- sim_mHMM(n_t = n_t, n = n, data_distr = 'continuous', gen = list(m = m, n_dep = n_dep),
                  gamma = gamma, emiss_distr = emiss_distr, var_gamma = .5, var_emiss = c(10, 5))

# Specify hyper-prior for the continuous emission distribution
manual_prior_emiss <- prior_emiss_cont(
                        gen = list(m = m, n_dep = n_dep),
                        emiss_mu0 = list(matrix(c(30, 70, 170), nrow = 1),
                                         matrix(c(7, 8, 18), nrow = 1)),
                        emiss_K0 = list(1, 1),
                        emiss_V =  list(rep(100, m), rep(25, m)),
                        emiss_nu = list(1, 1),
                        emiss_a0 = list(rep(1, m), rep(1, m)),
                        emiss_b0 = list(rep(1, m), rep(1, m)))

# Run the model on the simulated data:
# Note that for reasons of running time, J is set at a ridiculous low value.
# One would typically use a number of iterations J of at least 1000,
# and a burn_in of 200.
out_3st_cont_sim <- mHMM(s_data = data_cont$obs,
                         data_distr = 'continuous',
                         gen = list(m = m, n_dep = n_dep),
                         start_val = c(list(gamma), emiss_distr),
                         emiss_hyp_prior = manual_prior_emiss,
                         mcmc = list(J = 500, burn_in = 250))

summary(out_3st_cont_sim)

forward_probs <- forecast_mHMM1(object = out_3st_cont_sim, s_data = data_cont$obs, forecast_steps = 50)


vit_mHMM(out_3st_cont_sim, data_cont$obs, return_state_prob = TRUE)

# Plot forward probabilities:
set.seed(50)
forward_probs %>%
  as.data.frame() %>%
  dplyr::select(-state) %>%
  group_by(subj) %>%
  mutate(horizon = row_number()) %>%
  ungroup() %>%
  gather(state, value, -subj, -horizon) %>%
  mutate(subj = factor(subj)) %>%
  filter(subj %in% sample(size = 5, 1:n)) %>%
  ggplot(aes(x = horizon, y = value, group = subj, colour = subj)) +
  geom_line() +
  facet_grid(state~.) +
  theme_minimal()



data_cont <- sim_mHMM(n_t = n_t, n = n, data_distr = 'continuous', gen = list(m = m, n_dep = n_dep),
                      gamma = gamma, emiss_distr = emiss_distr, var_gamma = .5, var_emiss = c(10, 5))

out <- mHMM(s_data = data_cont$obs,
            data_distr = 'continuous',
            gen = list(m = m, n_dep = n_dep),
            start_val = c(list(gamma), emiss_distr),
            emiss_hyp_prior = manual_prior_emiss,
            mcmc = list(J = 500, burn_in = 250))

summary(out)


# hsmm


forward_probs <- forecast_mHMM1(object = out, s_data = data_cont$obs, forecast_steps = 50)

forward_probs2 <- forecast_mHMM2(object = out, s_data = data_cont$obs, forecast_steps = 50, show_progress = TRUE)

solve(t(diag(m) - gamma + 1), rep(1, m))

forward_probs
forward_probs2$forecast_probs

lapply(forward_probs2$predicted_states, function(s) do.call(rbind, lapply(1:nrow(s), function(r) table(s[r,])/length(s[r,]) )))

drawn_obs1 <- do.call(rbind,
        lapply(1:length(forward_probs2$predicted_obs), function(s) cbind(s, 1:nrow(forward_probs2$predicted_obs[[s]][[1]]), forward_probs2$predicted_obs[[s]][[1]])))
drawn_obs1 <- as.data.frame(drawn_obs1)
names(drawn_obs1) <- c("subj","horizon", paste0(1:n_iter))
drawn_obs1 <- drawn_obs1 %>%
  gather(iter, value, -subj, -horizon) %>%
  arrange(desc(subj), desc(horizon))

set.seed(42)

drawn_obs1 %>%
  mutate(subj = factor(subj)) %>%
  filter(subj %in% sample(size = 5, 1:n),
         iter > 250) %>%
  ggplot(aes(x = horizon, y = value, group = subj, colour = subj)) +
  geom_smooth()

# Plot mean line with standard error
drawn_obs1 %>%
  mutate(subj = factor(subj)) %>%
  filter(subj %in% sample(size = 5, 1:n),
         iter > 250) %>%
  ggplot(aes(x = horizon, y = value, group = subj, fill = subj)) +
  stat_summary(fun.data = mean_se, geom = "line") +
  stat_summary(fun.data = mean_se, geom = "ribbon", alpha = 0.2) +
  theme_minimal() +
  labs(title = "Mean Line with Standard Error by Group",
       x = "Horizon",
       y = "Value",
       colour = "Subject")

drawn_obs1 %>%
  mutate(subj = factor(subj)) %>%
  filter(subj %in% sample(size = 1, 1:n),
         # iter > 250,
         iter %in% sample(size = 20, 251:500)) %>%
  ggplot(aes(x = horizon, y = value, group = interaction(subj,iter), colour = subj)) +
  geom_line(alpha = 0.5) +
  theme_minimal() +
  labs(title = "Mean Line with Standard Error by Group",
       x = "Horizon",
       y = "Value",
       colour = "Subject")

# Calculate mean and standard error
summary_data <- drawn_obs1 %>%
  mutate(subj = factor(subj)) %>%
  filter(subj %in% sample(size = 5, 1:n),
         iter > 250) %>%
  group_by(subj, horizon) %>%
  summarise(mean_value = mean(value),
            median_value = median(value),
            se = sd(value) / sqrt(n()))

# Plot mean line with standard error ribbon
ggplot(summary_data, aes(x = horizon, y = mean_value, group = subj, fill = subj)) +
  geom_line() +
  geom_ribbon(aes(ymin = mean_value - se, ymax = mean_value + se), alpha = 0.3) +
  theme_minimal() +
  labs(title = "Mean Line with Standard Error Ribbon",
       x = "Horizon",
       y = "Value")

ggplot(summary_data, aes(x = horizon, y = median_value, group = subj, fill = subj)) +
  geom_line() +
  geom_ribbon(aes(ymin = median_value - se, ymax = median_value + se), alpha = 0.3) +
  theme_minimal() +
  labs(title = "Mean Line with Standard Error Ribbon",
       x = "Horizon",
       y = "Value")


ggplot(summary_data, aes(x = horizon, y = mean_value, group = subj, fill = subj)) +
  geom_ribbon(aes(ymin = as.numeric(subj) - se, ymax = as.numeric(subj) + se), alpha = 0.3) +
  theme_minimal() +
  labs(title = "Mean Line with Standard Error Ribbon",
       x = "Horizon",
       y = "Value")

# warnings()


data_cont$obs
data_cont$states

forward_probs

yardstick::f_meas_vec(truth = factor(data_cont$states[,2], levels = 1:m),
                      estimate = factor(vit_mHMM(out_3st_cont_sim, data_cont$obs, return_state_prob = FALSE)[,2], levels = 1:m))



# Load necessary libraries
library(yardstick)
library(dplyr)
library(purrr)

# Define the wrapping function
simulate_fit_forecast_evaluate <- function(n_t, n, m, n_dep,
                                           gamma, emiss_distr, dwell_distr,
                                           var_gamma, var_emiss, var_dwell,
                                           hyp_prior_emiss, hyp_prior_dwell,
                                           fit_fraction = 0.8, forecast_steps = 50,
                                           n_iter = 500, burn_in = 250) {

  # Step 1: Simulate the data
  data_cont <- sim_mHMM(n_t = n_t, n = n, data_distr = 'continuous', gen = list(m = m, n_dep = n_dep),
                        gamma = gamma, emiss_distr = emiss_distr, var_gamma = var_gamma, var_emiss = var_emiss)

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

  # Add a time column to the true states
  true_states <- true_states %>%
    mutate(state = as.integer(state)) %>%
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

  split_obs <- split_data(obs, fit_fraction)
  fit_data <- split_obs$fit
  forecast_data <- split_obs$forecast

  # Step 3: Fit the HMM model
  out <- mHMM(s_data = as.matrix(fit_data[,-4]),
              data_distr = 'continuous',
              gen = list(m = m, n_dep = n_dep),
              start_val = c(list(gamma), emiss_distr),
              emiss_hyp_prior = hyp_prior_emiss,
              mcmc = list(J = n_iter, burn_in = burn_in))

  # Step 4: Forecast using the fitted model
  forecast_results <- forecast_mHMM1(object = out, s_data = as.matrix(forecast_data[,-4]), forecast_steps = forecast_steps)

  # Ensure forecast_results has the necessary columns
  colnames(forecast_results) <- c("subj", "state", "time", "pr_state_1", "pr_state_2", "pr_state_3")

  forecast_results <- as.data.frame(forecast_results) %>%
    mutate(state = as.integer(state))

  # Step 5: Evaluate the predictions
  evaluate_forecasts <- function(forecast_results, true_states, m) {
    cumulative_f1 <- function(results, true_states, max_time) {
      results %>%
        filter(time <= max_time) %>%
        left_join(true_states %>% filter(time <= max_time), by = c("subj", "time")) %>%
        group_by(subj) %>%
        summarize(f1 = yardstick::f_meas_vec(factor(state.y, levels = 1:m), factor(state.x, levels = 1:m)), .groups = 'drop') %>%
        summarize(mean_f1 = mean(f1, na.rm = TRUE))
    }

    cumulative_f1_by_state <- function(results, true_states, max_time) {
      results %>%
        filter(time <= max_time) %>%
        left_join(true_states %>% filter(time <= max_time), by = c("subj", "time")) %>%
        group_by(subj, state.x) %>%
        summarize(f1 = yardstick::f_meas_vec(factor(state.y, levels = 1:m), factor(state.x, levels = 1:m)), .groups = 'drop') %>%
        group_by(state.x) %>%
        summarize(mean_f1 = mean(f1, na.rm = TRUE)) %>%
        mutate(time = max_time)
    }

    max_time <- max(forecast_results$time)
    f1_scores <- tibble(time = 1:max_time, f1 = map_dbl(1:max_time, ~ cumulative_f1(forecast_results, true_states, .x)$mean_f1))
    f1_scores_by_state <- bind_rows(lapply(1:max_time, function(t) cumulative_f1_by_state(forecast_results, true_states, t)))

    list(overall = f1_scores, by_state = f1_scores_by_state)
  }

  # Return the evaluation results
  evaluate_forecasts(forecast_results, true_states, m)
}

# Example usage
evaluation_results <- simulate_fit_forecast_evaluate(n_t = 500, n = 200, m = 3, n_dep = 2,
                                                     gamma = gamma, emiss_distr = emiss_distr,
                                                     var_gamma = 0.5, var_emiss = c(10, 5),
                                                     manual_prior_emiss = manual_prior_emiss,
                                                     fit_fraction = 0.8, forecast_steps = 50,
                                                     n_iter = 300, burn_in = 200)
evaluation_results$overall %>%
  as.data.frame()

evaluation_results$by_state %>%
  rename("state" = "state.x") %>%
  spread(key = state, value = mean_f1) %>%
  as.data.frame()

evaluation_results$by_state %>%
  rename("state" = "state.x") %>%
  mutate(state = factor(state)) %>%
  ggplot(aes(x = time, y = mean_f1, group =state, colour = state)) +
  geom_line() +
  theme_minimal()






probs[[1]][,4:6] %*% t(outer(1:100, Y = emiss[[q]][,1], FUN = stats::dnorm, sd = emiss[[q]][,2]))


apply(probs[[1]][,4:6] %*% t(outer(0:100, Y = emiss[[q]][,1], FUN = stats::dnorm, sd = emiss[[q]][,2])),1,sum)

# Forecast probs for scoring values
# probs[[1]][,4:6] %*% t(outer(seq(0,100,25), Y = emiss[[q]][,1], FUN = stats::dnorm, sd = emiss[[q]][,2]))
probs[[1]][,4:6] %*% sapply(seq(-10,40,1), function(r) dnorm(r, emiss[[q]][,1],emiss[[q]][,2]))

# Forecast prob of scoring over value
# probs[[1]][,4:6] %*% t(1-outer(seq(0,100,25), Y = emiss[[q]][,1], FUN = stats::pnorm, sd = emiss[[q]][,2]))
probs[[1]][,4:6] %*% (1-sapply(seq(-10,40,1), function(r) pnorm(r, emiss[[q]][,1],emiss[[q]][,2])))

# Forecast prob of scoring under value
# probs[[1]][,4:6] %*% t(outer(seq(0,100,25), Y = emiss[[q]][,1], FUN = stats::pnorm, sd = emiss[[q]][,2]))
probs[[1]][,4:6] %*% (sapply(seq(-10,40,1), function(r) pnorm(r, emiss[[q]][,1],emiss[[q]][,2])))

# Forecast prob
sapply(apply(probs[[1]][,4:6] %*% t(outer(seq(0,100,.1), Y = emiss[[q]][,1], FUN = stats::dnorm, sd = emiss[[q]][,2])),1,which.max), function(h) seq(0,100,.1)[h])

apply(probs[[1]][,4:6] %*% t(outer(seq(0,100,.1), Y = emiss[[q]][,1], FUN = stats::dnorm, sd = emiss[[q]][,2])),1,function(r) seq(0,100,.1)[which(r == quantile(r,0.5))])
apply(probs[[1]][,4:6] %*% t(outer(seq(0,100,.1), Y = emiss[[q]][,1], FUN = stats::dnorm, sd = emiss[[q]][,2])),1,function(r) seq(0,100,.1)[which(r == quantile(r,0.75))])
apply(probs[[1]][,4:6] %*% t(outer(seq(0,100,.1), Y = emiss[[q]][,1], FUN = stats::dnorm, sd = emiss[[q]][,2])),1,function(r) seq(0,100,.1)[which(r == quantile(r,0.25))])


# Forecast distribution
fvals <- probs[[29]][,4:6] %*% sapply(seq(-10,40,1), function(r) dnorm(r, emiss[[q]][,1],emiss[[q]][,2]))

data.frame(value = seq(-10,40,1), h_1 = fvals[1,], h_2 = fvals[2,], h_3 = fvals[3,], h_4 = fvals[4,], h_5 = fvals[5,],h_50 = fvals[5,]) %>%
  gather(horizon, prob, -value) %>%
  ggplot(aes(x = value, y = prob)) +
  geom_col() +
  facet_grid(horizon~.) +
  theme_minimal()


1-pnorm(1, mean = emiss[[q]][,1], sd = emiss[[q]][,2])

1-outer(1:10, Y = emiss[[q]][,1], FUN = stats::pnorm, sd = emiss[[q]][,2])




apply(exp(B_star),2,which.max)
apply(exp(Forward),2,which.max)
apply(exp(StateIn),2,which.max)
sample_path[[s]][,2]

mean(apply(exp(B_star),2,which.max)[1:500] == apply(exp(Forward),2,which.max))
mean(apply(exp(B_star),2,which.max)[1:500] == apply(exp(StateIn),2,which.max))
mean(apply(exp(B_star),2,which.max)[1:500] == sample_path[[s]][,2])

mean(apply(exp(Forward),2,which.max)[1:500] == apply(exp(StateIn),2,which.max))
mean(apply(exp(Forward),2,which.max)[1:500] == sample_path[[s]][,2])

mean(apply(exp(StateIn),2,which.max) == sample_path[[s]][,2])







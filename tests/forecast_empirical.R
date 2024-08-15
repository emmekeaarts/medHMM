
##
library(tidyverse)
library(medHMM)

# Load output:
path_medhmm <- "/Users/a6159737/Documents/Utrecht University/PhD/Projects/Simulation studies/medhmm-sim/adacko/MEDHMM_vs_MHMM/example/outputs/out_medHMM_m4_12dv_it2000_c1.rds"
path_mhmm <- "/Users/a6159737/Documents/Utrecht University/PhD/Projects/Simulation studies/medhmm-sim/adacko/MEDHMM_vs_MHMM/example/outputs/out_mHMM_m4_12dv_it2000_c1.rds"

path_mhmm <- "/Users/a6159737/Documents/Utrecht University/PhD/Projects/UMCG-MHMM/ESM-B HMM/outputs/results/out_cont_gamma_prior_emiss_prior_m4_12dv_it4000_c2.rds"

out <- out_medhmm <- readRDS(path_medhmm)
out_mhmm <- readRDS(path_mhmm)

#------------------------------------------------------------------------------#

# Backup data
data_list <- readRDS("/Users/a6159737/Documents/Utrecht University/PhD/Projects/UMCG-MHMM/inputs/backup.rds")

data <- data_list[[1]]

# Put data in right format and elongate data with "missing" night occasions
extended_data <- data %>%
    group_by(patient_id) %>%
    summarise(maxdays = max(dayno)) %>%
    group_by(patient_id) %>%
    summarise(dayno = rep(1:maxdays,each=8), beepno = rep(1:8, maxdays))

train_df <- left_join(extended_data, data) %>%
    dplyr::select(patient_id, time,
                  bs_diary_5, bs_diary_13, bs_diary_22,
                  bs_diary_15, bs_diary_9, bs_diary_10,
                  bs_diary_7, bs_diary_11, bs_diary_17,
                  bs_diary_8, bs_diary_14, bs_diary_16) %>%
    group_by(patient_id) %>%
    mutate(patient_id = cur_group_id(), time = row_number()) %>%
    ungroup() %>%
    arrange(patient_id, time) %>%
    dplyr::select(-time) %>%
    # drop_na() %>%
    as.matrix()

# Fit model:
m <- out$input$m
n_dep <- out$input$n_dep
burn_in <- out$input$burn_in
n_iter <- out$input$J


## Fixed effects:
# Transitions
gamma_ppc <- out$gamma_prob_bar %>%
    as.data.frame() %>%
    mutate(iter = row_number()) %>%
    filter(iter > burn_in) %>%
    dplyr::select(-iter) %>%
    summarise(across(.cols = everything(), median)) %>%
    as.numeric() %>%
    matrix(., nrow = m, byrow = TRUE)

gamma_ppc[1,-1] <- gamma_ppc[1,-1] + (1-apply(gamma_ppc, 1, sum)[1])/3
gamma_ppc[2,-2] <- gamma_ppc[2,-2] + (1-apply(gamma_ppc, 1, sum)[2])/3
gamma_ppc[3,-3] <- gamma_ppc[3,-3] + (1-apply(gamma_ppc, 1, sum)[3])/3
gamma_ppc[4,-4] <- gamma_ppc[4,-4] + (1-apply(gamma_ppc, 1, sum)[4])/3

gamma_start <- gamma_ppc

# Emissions
emiss_ppc <- lapply(1:length(out_mhmm$emiss_mu_bar), function(s) {

    emiss <- out_mhmm$emiss_mu_bar[[s]] %>%
        as.data.frame() %>%
        mutate(iter = row_number()) %>%
        filter(iter > burn_in) %>%
        dplyr::select(-iter) %>%
        summarise(across(.cols = everything(), median)) %>%
        gather(mu, value) %>%
        pull(value) %>%
        matrix(., nrow = m)

    emiss_var<- out_mhmm$emiss_var_bar[[s]] %>%
        as.data.frame() %>%
        mutate(iter = row_number()) %>%
        filter(iter > burn_in) %>%
        dplyr::select(-iter) %>%
        summarise(across(.cols = everything(), median)) %>%
        gather(mu, value) %>%
        pull(value) %>%
        matrix(., nrow = m)

    cbind(emiss, emiss_var)

})

emiss_start <- emiss_ppc

# Dwell time
self_trans <- readRDS(path_mhmm)$gamma_prob_bar %>%
    as.data.frame() %>%
    mutate(iter = row_number()) %>%
    filter(iter > burn_in) %>%
    dplyr::select(-iter) %>%
    summarise(across(.cols = everything(), median)) %>%
    as.numeric() %>%
    matrix(., nrow = m, byrow = TRUE) %>%
    diag()
expected_mean_dwell_times<-1/-log(self_trans)

# expected_mean_dwell_times <- c(17.1, 9.8, 8.8, 16.3)
dwell_distr <- matrix(log(expected_mean_dwell_times), nrow = m, ncol = 1, byrow = TRUE)
dwell_start <- matrix(expected_mean_dwell_times, nrow = m, ncol = 1, byrow = TRUE)

## Hyper priors
# Specify hyper-prior for the emission distribution
hyp_prior_emiss <- list(
    emiss_mu0 = list(matrix(c(5, 25, 50, 50), nrow = 1), # down
                     matrix(c(5, 25, 50, 50), nrow = 1), # dread rest of day
                     matrix(c(5, 25, 50, 50), nrow = 1), # worry
                     matrix(c(5, 25, 50, 50), nrow = 1), # inadequate
                     matrix(c(25, 25, 50, 50), nrow = 1), # tired
                     matrix(c(50, 25, 5, 5), nrow = 1), # content

                     matrix(c(5, 25, 25, 10), nrow = 1), # agitated
                     matrix(c(5, 25, 25, 10), nrow = 1), # irritated
                     matrix(c(25, 25, 25, 10), nrow = 1), # switch
                     matrix(c(25, 25, 25, 10), nrow = 1), # extremely well
                     matrix(c(25, 25, 25, 10), nrow = 1), # ideas
                     matrix(c(25, 25, 25, 10), nrow = 1)  # thoughts racing
    ),
    emiss_K0  = rep(list(1),n_dep),
    emiss_nu  = rep(list(1),n_dep),
    emiss_V   = rep(list(rep(400, m)),n_dep),
    emiss_a0  = rep(list(rep(0.001, m)),n_dep),
    emiss_b0  = rep(list(rep(0.001, m)),n_dep))

# hyp_prior_dwell = list(
#     dwell_mu0 = matrix(log(expected_mean_dwell_times), nrow = 1, ncol = m), # nrow = number of covariates + 1; ncol = number of hidden states
#     dwell_K0  = c(0.01),
#     dwell_nu  = c(1),
#     dwell_V   = rep(0.01, m)
# )

hyp_prior_dwell = list(
    dwell_mu0 = matrix(log(expected_mean_dwell_times), nrow = 1, ncol = m), # nrow = number of covariates + 1; ncol = number of hidden states
    dwell_K0  = c(0.1),
    dwell_nu  = c(1),
    dwell_V   = rep(10, m)
)


# Fit model:
n_iter <- 1000
burn_in <- 500

# Step 3: Fit the HMM model
set.seed(42)
out_medhmm <- medHMM_cont_shiftpois(s_data = train_df,
                                    gen = list(m = m, n_dep = n_dep),
                                    start_val = c(list(gamma_start), emiss_start, list(dwell_start)),
                                    emiss_hyp_prior = hyp_prior_emiss,
                                    dwell_hyp_prior = hyp_prior_dwell,
                                    shift = 1,
                                    show_progress = TRUE,
                                    mcmc = list(J = n_iter, burn_in = burn_in),
                                    return_path = TRUE,
                                    max_dwell = 56)

train_df

out <- out_medhmm

# saveRDS(out, "/Users/a6159737/Documents/Utrecht University/PhD/Projects/UMCG-MHMM/ESM-B HMM/medHMM/results/out_medHMM_shifpois_m4_12dv_it1000_Mx56_c1.rds")

out_medhmm <- read_rds("/Users/a6159737/Documents/Utrecht University/PhD/Projects/UMCG-MHMM/ESM-B HMM/medHMM/results/out_medHMM_shifpois_m4_12dv_it1000_Mx56_c1.rds")

#==============================================================================#
# Main results (figures and tables)
#==============================================================================#

#------------------------------------------------------------------------------#
# Figure: Group- and patient-level emission means, linked by patient:

# Figure R.1: group-level emissions

# Extract emission means
emiss_mu_bar <-  do.call(rbind, lapply(1:length(out$emiss_mu_bar), function(q){
    out$emiss_mu_bar[[q]] %>%
        as.data.frame() %>%
        mutate(iter = row_number(),
               dep = names(out$emiss_mu_bar)[q]) %>%
        gather(mu, value, -iter, -dep)
}))

# Add labels for the constructs

dep_vars <- c("bs_diary_5", "bs_diary_13", "bs_diary_22",
              "bs_diary_15", "bs_diary_9", "bs_diary_10")
man_vars <- c("bs_diary_7", "bs_diary_11", "bs_diary_17",
              "bs_diary_8", "bs_diary_14", "bs_diary_16")

# Plot

# Patient-specific values instead of group-level MAP
subj_emiss <- do.call(rbind, lapply(1:length(out$PD_subj),function(s){
    out$PD_subj[[s]] %>%
        as.data.frame() %>%
        dplyr::select("dep1_mu_S1":"dep12_mu_S4") %>%
        mutate(patient_id = s, iter = row_number()) %>%
        filter(iter > burn_in) %>%
        dplyr::select(-iter) %>%
        gather(dep, value, -patient_id) %>%
        group_by(patient_id, dep) %>%
        summarise(value = median(value)) %>%
        separate(dep, into = c("dep","mu"), sep = "_mu_") %>%
        mutate(dep = factor(dep, levels = paste0("dep",1:12), labels = c("bs_diary_5", "bs_diary_13", "bs_diary_22",
                                                                         "bs_diary_15", "bs_diary_9", "bs_diary_10",
                                                                         "bs_diary_7", "bs_diary_11", "bs_diary_17",
                                                                         "bs_diary_8", "bs_diary_14", "bs_diary_16")))

} ))

# Subject and Group level parameters
subj_emiss_nice <- subj_emiss %>%
    mutate(construct = case_when(dep %in% dep_vars ~ "Depression items",
                                 dep %in% man_vars ~ "Mania items"),
           construct = factor(construct, levels = c("Depression items","Mania items")),
           dep = factor(dep,
                        levels = c("bs_diary_5", "bs_diary_13", "bs_diary_22",
                                   "bs_diary_15", "bs_diary_9", "bs_diary_10",
                                   "bs_diary_7", "bs_diary_11", "bs_diary_17",
                                   "bs_diary_8", "bs_diary_14", "bs_diary_16"),
                        labels = c("down","dread rest of day","worry",
                                   "inadequate", "tired", "content",
                                   "agitated", "irritated", "switch and focus",
                                   "extremely well", "full of ideas", "racing thoughts")),
           mu = factor(mu, levels = paste0("S",1:4), labels = c("neutral","elevated","mixed","lowered"))) %>%
    rename("state" = "mu")

group_emiss_nice <- emiss_mu_bar %>%
    filter(iter > burn_in) %>%
    group_by(dep, mu) %>%
    summarise(group_mean = mean(value),
              group_median = median(value),
              group_sd = sd(value),
              group_CCIlwr = quantile(value, 0.025),
              group_CCIupr = quantile(value, 0.975)) %>%
    mutate(construct = case_when(dep %in% dep_vars ~ "Depression items",
                                 dep %in% man_vars ~ "Mania items"),
           construct = factor(construct, levels = c("Depression items","Mania items")),
           dep = factor(dep,
                        levels = c("bs_diary_5", "bs_diary_13", "bs_diary_22",
                                   "bs_diary_15", "bs_diary_9", "bs_diary_10",
                                   "bs_diary_7", "bs_diary_11", "bs_diary_17",
                                   "bs_diary_8", "bs_diary_14", "bs_diary_16"),
                        labels = c("down","dread rest of day","worry",
                                   "inadequate", "tired", "content",
                                   "agitated", "irritated", "switch and focus",
                                   "extremely well", "full of ideas", "racing thoughts")),
           mu = factor(mu, levels = paste0("mu_",1:4), labels = c("neutral","elevated","mixed","lowered"))) %>%
    rename("state" = "mu")

p <- ggplot(data = subj_emiss_nice) +
    geom_line(aes(x = state, y = value, group = patient_id), alpha = 0.2) +
    geom_jitter(aes(x = state, y = value, colour = state), width = 0.15, alpha = 0.5) +
    geom_pointrange(data = group_emiss_nice,
                    aes(x = state, y = group_mean,
                        ymin = group_CCIlwr, ymax = group_CCIupr), size = 0.25) +
    scale_color_viridis_d(option = "plasma", direction = -1) +
    geom_hline(yintercept = 0) +
    facet_wrap(construct+dep~., scales = "fixed", nrow = 2) +
    theme_minimal() +
    ylab(label = "EMA score") +
    xlab(label = "") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    theme(legend.position = "none") +
    # ggtitle(label = "Composition of momentary mood states: group-level and\npatient-specific item emission scores by state") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(strip.text.y = element_text(angle = 0))

p

# ggsave(p + theme(text = element_text(size = 14)), filename = "outputs/figures/linked_emiss_a.pdf", dpi = 600)
# ggsave(p, filename = "outputs/figures/linked_emiss_a.jpeg", dpi = 600)
# ggsave(p, filename = "outputs/figures/linked_emiss_a.png", dpi = 600)


# # Reordering variables:
# # Subject and Group level parameters
# subj_emiss_nice <- subj_emiss %>%
#     mutate(construct = case_when(dep %in% dep_vars ~ "Depression items",
#                                  dep %in% man_vars ~ "Mania items"),
#            construct = factor(construct, levels = c("Depression items","Mania items")),
#            dep = factor(dep,
#                         levels = c("bs_diary_5", "bs_diary_13", "bs_diary_22",
#                                    "bs_diary_15", "bs_diary_9", "bs_diary_10",
#                                    "bs_diary_17",
#                                    "bs_diary_8", "bs_diary_14", "bs_diary_16",
#                                    "bs_diary_7", "bs_diary_11"),
#                         labels = c("down","dread rest of day","worry",
#                                    "inadequate", "tired", "content",
#                                    "switch and focus",
#                                    "extremely well", "full of ideas", "racing thoughts",
#                                    "agitated", "irritated")),
#            mu = factor(mu, levels = paste0("S",1:4), labels = c("neutral","elevated","mixed","lowered"))) %>%
#     rename("state" = "mu")
#
# group_emiss_nice <- emiss_mu_bar %>%
#     filter(iter > 1000) %>%
#     group_by(dep, mu) %>%
#     summarise(group_mean = mean(value),
#               group_median = median(value),
#               group_sd = sd(value),
#               group_CCIlwr = quantile(value, 0.025),
#               group_CCIupr = quantile(value, 0.975)) %>%
#     mutate(construct = case_when(dep %in% dep_vars ~ "Depression items",
#                                  dep %in% man_vars ~ "Mania items"),
#            construct = factor(construct, levels = c("Depression items","Mania items")),
#            dep = factor(dep,
#                         levels = c("bs_diary_5", "bs_diary_13", "bs_diary_22",
#                                    "bs_diary_15", "bs_diary_9", "bs_diary_10",
#                                    "bs_diary_17",
#                                    "bs_diary_8", "bs_diary_14", "bs_diary_16",
#                                    "bs_diary_7", "bs_diary_11"),
#                         labels = c("down","dread rest of day","worry",
#                                    "inadequate", "tired", "content",
#                                    "switch and focus",
#                                    "extremely well", "full of ideas", "racing thoughts",
#                                    "agitated", "irritated")),
#            mu = factor(mu, levels = paste0("mu_",1:4), labels = c("neutral","elevated","mixed","lowered"))) %>%
#     rename("state" = "mu")
#
# p <- ggplot(data = subj_emiss_nice) +
#     geom_line(aes(x = state, y = value, group = patient_id), alpha = 0.2) +
#     geom_jitter(aes(x = state, y = value, colour = state), width = 0.15, alpha = 0.5) +
#     geom_pointrange(data = group_emiss_nice,
#                     aes(x = state, y = group_mean,
#                         ymin = group_CCIlwr, ymax = group_CCIupr), size = 0.25) +
#     scale_color_viridis_d(option = "plasma", direction = -1) +
#     geom_hline(yintercept = 0) +
#     facet_wrap(construct+dep~., scales = "fixed", nrow = 2) +
#     theme_minimal() +
#     ylab(label = "EMA score") +
#     xlab(label = "") +
#     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
#     theme(legend.position = "none") +
#     # ggtitle(label = "Composition of momentary mood states: group-level and\npatient-specific item emission scores by state") +
#     theme(plot.title = element_text(hjust = 0.5)) +
#     theme(strip.text.y = element_text(angle = 0))
#
# p
#
# ggsave(p + theme(text = element_text(size = 14)), filename = "outputs/figures/linked_emiss_b.pdf", dpi = 600)
# # ggsave(p, filename = "outputs/figures/linked_emiss_b.jpeg", dpi = 600)
# # ggsave(p + theme(text = element_text(size = 12)), filename = "outputs/figures/linked_emiss_b.png", dpi = 600)


#------------------------------------------------------------------------------#
# Figure R.4: validated scales (decoding)

# Load weekly data
bs_asrm <- data_list[[2]] %>%
    arrange(patient_id)
bs_qids <- data_list[[3]] %>%
    arrange(patient_id)

# Put together
bs_scales <- full_join(bs_asrm, bs_qids) %>%
    dplyr::select(patient_id, bs_asrm_open_from, bs_asrm_tot, bs_qids_tot) %>%
    rename("open_from" = "bs_asrm_open_from") %>%
    group_by(patient_id) %>%
    mutate(from_time_day = format(as.POSIXct(lag(open_from, 1, default = min(open_from)-(3600*24*7))),format = "%Y-%m-%d"),
           to_time_day = format(as.POSIXct(open_from-1),format = "%Y-%m-%d"),
           patient_id = factor(patient_id)) %>%
    ungroup()

# Get decoding
# states <- vit_mHMM_cont_mar(object = out, s_data = train_df, burn_in = 1000) %>%
#     as.data.frame() %>%
#     gather(subject, state) %>%
#     # drop_na() %>%
#     mutate(subject = factor(subject, levels = paste0("Subj_",1:20))) %>%
#     group_by(subject) %>%
#     mutate(subject = cur_group_id()) %>%
#     ungroup() %>%
#     rename("patient_id" = "subject") %>%
#     group_by(patient_id) %>%
#     mutate(occasion = row_number())

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

states <- local_decoding(out_medhmm)

states <- states %>%
    rename("patient_id" = "subj") %>%
    select(-c(paste0("pr_state",1:m)))

data_labelled <- left_join(as.data.frame(train_df) %>% group_by(patient_id) %>%
                               mutate(occasion = row_number()), states) %>%
    group_by(patient_id) %>%
    mutate(occasion = row_number()) %>%
    as.data.frame()

# Add week info
data_labelled <- bind_cols(data_labelled,
                           left_join(extended_data, data) %>%
                               group_by(patient_id) %>%
                               mutate(patient_id = cur_group_id()) %>%
                               ungroup() %>%
                               arrange(patient_id) %>%
                               dplyr::select(open_from, bs_diary_open_from, bs_diary_date)) %>%
    fill(open_from)

for(r in 1:length(data_labelled$bs_diary_open_from)){
    if(is.na(data_labelled$bs_diary_open_from[r])){
        data_labelled$bs_diary_open_from[r] <- data_labelled$bs_diary_open_from[r-1] + 3*60*60
    }
}


# Validated scales
state_data <- data_labelled %>%
    select(patient_id, state, bs_diary_open_from, occasion) %>%
    rename("open_from" = "bs_diary_open_from") %>%
    group_by(patient_id) %>%
    mutate(state = factor(state, levels = 1:4, labels = c("euthymic","manic","mixed","depressive")),
           from_time = lag(open_from, 1, default = min(open_from)-(3600*3)),
           to_time = open_from-1,
           patient_id = factor(patient_id),
           from_occ = occasion-1,
           to_occ = occasion) %>%
    mutate(day = format(as.POSIXct(to_time), format = "%Y-%m-%d")) %>%
    ungroup()

scale_data <- bs_scales %>%
    group_by(patient_id) %>%
    mutate(patient_id = cur_group_id(),
           patient_id = as.factor(patient_id),
           occasion = row_number(),
           day = from_time_day,
           bs_asrm_tot = case_when(is.na(bs_asrm_tot) ~ 99,
                                   !is.na(bs_asrm_tot) ~ bs_asrm_tot),
           bs_qids_tot = case_when(is.na(bs_qids_tot) ~ 99,
                                   !is.na(bs_qids_tot) ~ bs_qids_tot)
    ) %>%
    ungroup() %>%
    select(-occasion, -open_from)

join_data <- full_join(state_data, scale_data) %>%
    tidyr::fill(bs_asrm_tot:bs_qids_tot, .direction = "down") %>%
    gather(variable, value, -patient_id,-state,-open_from,-occasion,-from_time,-to_time,-from_occ,-to_occ,-from_time_day,-to_time_day,-day) %>%
    mutate(value = case_when(value == 99 ~ NA_real_,
                             value != 99 ~ value))

library(viridis)
p <- ggplot() +
    geom_line(data = join_data %>%
                  group_by(patient_id) %>%
                  mutate(variable = factor(variable,
                                           levels = c("bs_asrm_tot","bs_qids_tot"),
                                           labels = c("ASRM","QIDS"))
                         # patient_id = factor(patient_id, levels =  c(15,12,19,17,1,3,
                         #                                             6,5,7,18,16,13,
                         #                                             20,2,14,11,9,8,4,10),
                         #                     labels = paste0("patient ",c(15,12,19,17,1,3,
                         #                                                  6,5,7,18,16,13,
                         #                                                  20,2,14,11,9,8,4,10)))
                         ),
              aes(x = from_occ, y = value, colour = variable, linetype=variable)) +
    geom_rect(data = state_data
              # %>%
              #     mutate(patient_id = factor(patient_id, levels = c(15,12,19,17,1,3,
              #                                                       6,5,7,18,16,13,
              #                                                       20,2,14,11,9,8,4,10),
              #                                labels = paste0("patient ",c(15,12,19,17,1,3,
              #                                                             6,5,7,18,16,13,
              #                                                             20,2,14,11,9,8,4,10))))
              , aes(xmin = from_occ, xmax = to_occ,
                                                                                                       ymin = 25, ymax = 30, fill = state)) +
    scale_fill_viridis(discrete = TRUE, alpha=0.6,direction = -1, option = "plasma") +
    geom_hline(yintercept = 6, linetype = "dashed") +
    facet_wrap(patient_id~., scales = "free", ncol = 4) +
    theme_minimal() +
    labs(colour = "Questionnaire", fill = "Mood state", linetype = "Questionnaire") +
    theme(legend.position = "bottom") +
    # ggtitle(label = "Temporal alignment between mood states and weekly symptom scores") +
    theme(plot.title = element_text(hjust = 0.5)) +
    ylab(label = "Weekly symptom score") +
    xlab(label = "Measurement occasion")

p

# ggsave(plot = p + theme(text = element_text(size = 12)),
#        filename = "outputs/figures/decoding_12.pdf", width = 14, height = 11, units = "in", dpi = 300)
# ggsave(plot = p + theme(text = element_text(size = 12)),
#        filename = "outputs/figures/decoding_12.tiff", width = 14, height = 11, units = "in", dpi = 300)
# ggsave(plot = p + theme(text = element_text(size = 14)),
#        filename = "outputs/figures/decoding_14.pdf", width = 14, height = 11, units = "in", dpi = 300)
# ggsave(plot = p + theme(text = element_text(size = 14)),
#        filename = "outputs/figures/decoding_14.tiff", width = 14, height = 11, units = "in", dpi = 300)
# ggsave(plot = p + theme(text = element_text(size = 16)),
#        filename = "outputs/figures/decoding_16.pdf", width = 14, height = 11, units = "in", dpi = 300)
# ggsave(plot = p + theme(text = element_text(size = 16)),
#        filename = "outputs/figures/decoding_16.tiff", width = 14, height = 11, units = "in", dpi = 300)
# ggsave(plot = p + theme(text = element_text(size = 18)),
#        filename = "outputs/figures/decoding_18.pdf", width = 14, height = 11, units = "in", dpi = 300)
# ggsave(plot = p + theme(text = element_text(size = 18)),
#        filename = "outputs/figures/decoding_18.tiff", width = 14, height = 11, units = "in", dpi = 300)
# ggsave(plot = p + theme(text = element_text(size = 20)),
#        filename = "outputs/figures/decoding_20.pdf", width = 14, height = 11, units = "in", dpi = 300)
# ggsave(plot = p + theme(text = element_text(size = 20)),
#        filename = "outputs/figures/decoding_20.tiff", width = 14, height = 11, units = "in", dpi = 300)


#------------------------------------------------------------------------------#

# Retrain including clinical scales:

# Calculate scale averages by states:
scale_scores <- join_data %>%
    drop_na() %>%
    group_by(variable, state) %>%
    summarise(mean_scale = mean(value), var_scale = var(value)) %>%
    select(-state) %>%
    group_by(variable) %>%
    group_split(.keep = FALSE)

scale_scores <- lapply(scale_scores, as.matrix)

# Prepare data
scales_one_day <- full_join(state_data, scale_data) %>%
    # tidyr::fill(bs_asrm_tot:bs_qids_tot, .direction = "down") %>%
    mutate(bs_asrm_tot = case_when(bs_asrm_tot == 99 ~ NA_real_,
                                   bs_asrm_tot != 99 ~ bs_asrm_tot),
           bs_qids_tot = case_when(bs_qids_tot == 99 ~ NA_real_,
                                   bs_qids_tot != 99 ~ bs_qids_tot),
           patient_id = as.numeric(as.character(patient_id))) %>%
    select(patient_id, bs_asrm_tot, bs_qids_tot, occasion)

scales_full_days <- full_join(state_data, scale_data) %>%
    tidyr::fill(bs_asrm_tot:bs_qids_tot, .direction = "down") %>%
    mutate(bs_asrm_tot = case_when(bs_asrm_tot == 99 ~ NA_real_,
                                   bs_asrm_tot != 99 ~ bs_asrm_tot),
           bs_qids_tot = case_when(bs_qids_tot == 99 ~ NA_real_,
                                   bs_qids_tot != 99 ~ bs_qids_tot),
           patient_id = as.numeric(as.character(patient_id))) %>%
    select(patient_id, bs_asrm_tot, bs_qids_tot, occasion)

# Create new train sets:
train_df_scales_one_day <- inner_join(train_df %>%
               as.data.frame() %>%
               group_by(patient_id) %>%
               mutate(occasion = row_number()), scales_one_day) %>%
    select(-occasion) %>%
    as.matrix()

train_df_scales_full_days <- inner_join(train_df %>%
               as.data.frame() %>%
               group_by(patient_id) %>%
               mutate(occasion = row_number()), scales_full_days) %>%
    select(-occasion) %>%
    as.matrix()



# Define new hyper priors and train:
n_dep <- 14

# Emissions
emiss_ppc <- lapply(1:length(out_mhmm$emiss_mu_bar), function(s) {

    emiss <- out_mhmm$emiss_mu_bar[[s]] %>%
        as.data.frame() %>%
        mutate(iter = row_number()) %>%
        filter(iter > burn_in) %>%
        dplyr::select(-iter) %>%
        summarise(across(.cols = everything(), median)) %>%
        gather(mu, value) %>%
        pull(value) %>%
        matrix(., nrow = m)

    emiss_var<- out_mhmm$emiss_var_bar[[s]] %>%
        as.data.frame() %>%
        mutate(iter = row_number()) %>%
        filter(iter > burn_in) %>%
        dplyr::select(-iter) %>%
        summarise(across(.cols = everything(), median)) %>%
        gather(mu, value) %>%
        pull(value) %>%
        matrix(., nrow = m)

    cbind(emiss, emiss_var)

})

emiss_start <- emiss_ppc

emiss_start <- c(emiss_start, scale_scores)

# Dwell time
self_trans <- readRDS(path_mhmm)$gamma_prob_bar %>%
    as.data.frame() %>%
    mutate(iter = row_number()) %>%
    filter(iter > burn_in) %>%
    dplyr::select(-iter) %>%
    summarise(across(.cols = everything(), median)) %>%
    as.numeric() %>%
    matrix(., nrow = m, byrow = TRUE) %>%
    diag()
expected_mean_dwell_times<-1/-log(self_trans)

# expected_mean_dwell_times <- c(17.1, 9.8, 8.8, 16.3)
dwell_distr <- matrix(log(expected_mean_dwell_times), nrow = m, ncol = 1, byrow = TRUE)
dwell_start <- matrix(expected_mean_dwell_times, nrow = m, ncol = 1, byrow = TRUE)

## Hyper priors
# Specify hyper-prior for the emission distribution
hyp_prior_emiss <- list(
    emiss_mu0 = list(matrix(c(5, 25, 50, 50), nrow = 1), # down
                     matrix(c(5, 25, 50, 50), nrow = 1), # dread rest of day
                     matrix(c(5, 25, 50, 50), nrow = 1), # worry
                     matrix(c(5, 25, 50, 50), nrow = 1), # inadequate
                     matrix(c(25, 25, 50, 50), nrow = 1), # tired
                     matrix(c(50, 25, 5, 5), nrow = 1), # content

                     matrix(c(5, 25, 25, 10), nrow = 1), # agitated
                     matrix(c(5, 25, 25, 10), nrow = 1), # irritated
                     matrix(c(25, 25, 25, 10), nrow = 1), # switch
                     matrix(c(25, 25, 25, 10), nrow = 1), # extremely well
                     matrix(c(25, 25, 25, 10), nrow = 1), # ideas
                     matrix(c(25, 25, 25, 10), nrow = 1),  # thoughts racing

                     matrix(c(1.5, 5, 2.5, 1.5), nrow = 1), # asrm
                     matrix(c(2.5, 10, 10, 10), nrow = 1)  # qids
    ),
    emiss_K0  = rep(list(1),n_dep),
    emiss_nu  = rep(list(1),n_dep),
    emiss_V   = rep(list(rep(400, m)),n_dep),
    emiss_a0  = rep(list(rep(0.001, m)),n_dep),
    emiss_b0  = rep(list(rep(0.001, m)),n_dep))

hyp_prior_dwell = list(
    dwell_mu0 = matrix(log(expected_mean_dwell_times), nrow = 1, ncol = m), # nrow = number of covariates + 1; ncol = number of hidden states
    dwell_K0  = c(0.1),
    dwell_nu  = c(1),
    dwell_V   = rep(10, m)
)


# Fit model:
n_iter <- 100
burn_in <- 50

# Step 3: Fit the HMM model
set.seed(42)
out_medhmm_scales1 <- medHMM_cont_shiftpois(s_data = train_df_scales_one_day,
                                    gen = list(m = m, n_dep = n_dep),
                                    start_val = c(list(gamma_start), emiss_start, list(dwell_start)),
                                    emiss_hyp_prior = hyp_prior_emiss,
                                    dwell_hyp_prior = hyp_prior_dwell,
                                    shift = 1,
                                    show_progress = TRUE,
                                    mcmc = list(J = n_iter, burn_in = burn_in),
                                    return_path = TRUE,
                                    max_dwell = 56)

out_medhmm_scales2 <- medHMM_cont_shiftpois(s_data = train_df_scales_full_days,
                                            gen = list(m = m, n_dep = n_dep),
                                            start_val = c(list(gamma_start), emiss_start, list(dwell_start)),
                                            emiss_hyp_prior = hyp_prior_emiss,
                                            dwell_hyp_prior = hyp_prior_dwell,
                                            shift = 1,
                                            show_progress = TRUE,
                                            mcmc = list(J = n_iter, burn_in = burn_in),
                                            return_path = TRUE,
                                            max_dwell = 500)

n_iter <- 1000
burn_in <- 500
set.seed(42)
out_medhmm_scales3 <- medHMM_cont_shiftpois(s_data = train_df_scales_full_days,
                                            gen = list(m = m, n_dep = n_dep),
                                            start_val = c(list(gamma_start), emiss_start, list(dwell_start)),
                                            emiss_hyp_prior = hyp_prior_emiss,
                                            dwell_hyp_prior = hyp_prior_dwell,
                                            shift = 1,
                                            show_progress = TRUE,
                                            mcmc = list(J = n_iter, burn_in = burn_in),
                                            return_path = TRUE,
                                            max_dwell = 448)


out <- out_medhmm_scales1

out <- out_medhmm_scales2

out_medhmm_scales1$emiss_mu_bar



#
#------------------------------------------------------------------------------#
# Figure S.3: posterior predictive checks

m <- out$input$m
n_dep <- out$input$n_dep
burn_in <- out$input$burn_in
n_iter <- out$input$J


## Fixed effects:
# Transitions
gamma_ppc <- out$gamma_prob_bar %>%
    as.data.frame() %>%
    mutate(iter = row_number()) %>%
    filter(iter > burn_in) %>%
    dplyr::select(-iter) %>%
    summarise(across(.cols = everything(), median)) %>%
    as.numeric() %>%
    matrix(., nrow = m, byrow = TRUE)

gamma_ppc[1,-1] <- gamma_ppc[1,-1] + (1-apply(gamma_ppc, 1, sum)[1])/3
gamma_ppc[2,-2] <- gamma_ppc[2,-2] + (1-apply(gamma_ppc, 1, sum)[2])/3
gamma_ppc[3,-3] <- gamma_ppc[3,-3] + (1-apply(gamma_ppc, 1, sum)[3])/3
gamma_ppc[4,-4] <- gamma_ppc[4,-4] + (1-apply(gamma_ppc, 1, sum)[4])/3

# Emissions
emiss_ppc <- lapply(1:length(out$emiss_mu_bar), function(s) {

    emiss <- out$emiss_mu_bar[[s]] %>%
        as.data.frame() %>%
        mutate(iter = row_number()) %>%
        filter(iter > burn_in) %>%
        dplyr::select(-iter) %>%
        summarise(across(.cols = everything(), median)) %>%
        gather(mu, value) %>%
        pull(value) %>%
        matrix(., nrow = m)

    emiss_var<- out$emiss_var_bar[[s]] %>%
        as.data.frame() %>%
        mutate(iter = row_number()) %>%
        filter(iter > burn_in) %>%
        dplyr::select(-iter) %>%
        summarise(across(.cols = everything(), median)) %>%
        gather(mu, value) %>%
        pull(value) %>%
        matrix(., nrow = m)

    cbind(emiss, emiss_var)

})


## Random effect: between subject variance
# Transitions
gamma_var_ppc <- apply(out_mhmm$gamma_V_int_bar[(burn_in+1):n_iter,] %>%
                           as.data.frame() %>%
                           dplyr::select(paste0("var_int_S",rep(1:4,each=3),"toS",2:4,"_with_int_S",rep(1:4,each=3),"toS",2:4)),
                       2, median) %>%
    as.numeric()
gamma_var_ppc <- matrix(gamma_var_ppc, nrow = m, ncol = m-1, byrow = TRUE)

# Emissions
emiss_varmu_ppc <- lapply(out$emiss_varmu_bar, function(emiss) {
    emiss %>%
        as.data.frame() %>%
        mutate(iter = row_number()) %>%
        filter(iter > burn_in) %>%
        dplyr::select(-iter) %>%
        summarise(across(.cols = everything(), median)) %>%
        gather(mu, value) %>%
        dplyr::select(value) %>%
        as.matrix()
})

# Dwell
out$dwell_mu_bar
out$dwell_var_bar
out$dwell_varmu_bar

# Simulate data
set.seed(42)
sim_data <- pbapply::pblapply(1:500, function(s) mHMMbayes::sim_mHMM_plnorm(n_t = 642, n = 20,
                                                                            data_distr = "continuous",
                                                                            m = out$input$m, n_dep = out$input$n_dep,
                                                                            gamma = gamma_ppc,
                                                                            emiss_distr = emiss_ppc,
                                                                            var_gamma = gamma_var_ppc,
                                                                            var_emiss = emiss_varmu_ppc,
                                                                            return_ind_par = TRUE))

## PPC1: group-level mean

true_data <- train_df %>%
    as.data.frame() %>%
    group_by(patient_id) %>%
    mutate(occasion = row_number()) %>%
    gather(variable, value, -patient_id, -occasion)

ppc_data <- do.call(rbind, pbapply::pblapply(1:length(sim_data), function(s){

    data <- as.data.frame(sim_data[[s]]$obs)
    names(data) <- c("patient_id", out$input$dep_labels)
    data %>%
        group_by(patient_id) %>%
        mutate(occasion = row_number()) %>%
        gather(variable, value, -patient_id, -occasion) %>%
        mutate(rep = s)

}))

ppc1_truth <- true_data %>%
    group_by(variable) %>%
    summarise(mean_value = mean(value, na.rm = TRUE), sd_value = sd(value, na.rm = TRUE),
              min_value = min(value, na.rm = TRUE), max_value = max(value, na.rm = TRUE))

dep_vars <- c("bs_diary_5", "bs_diary_9", "bs_diary_10",
              "bs_diary_13", "bs_diary_15", "bs_diary_22")
man_vars <- c("bs_diary_7", "bs_diary_8", "bs_diary_11",
              "bs_diary_14", "bs_diary_16", "bs_diary_17")

# Mean
p <- ppc_data %>%
    group_by(variable, rep) %>%
    summarise(mean_value = mean(value), sd_value = sd(value)) %>%
    mutate(variable = factor(variable, levels = c("bs_diary_5", "bs_diary_13", "bs_diary_22",
                                                  "bs_diary_15", "bs_diary_9", "bs_diary_10",
                                                  "bs_diary_7", "bs_diary_11", "bs_diary_17",
                                                  "bs_diary_8", "bs_diary_14", "bs_diary_16"),
                             labels = c("down","dread nrest of day","worry",
                                        "inadequate", "tired", "content",
                                        "agitated", "irritated", "switch and focus",
                                        "extremely well", "full of ideas", "thoughts are racing"))) %>%
    ggplot(aes(x = mean_value)) +
    geom_histogram() +
    geom_vline(data = ppc1_truth %>%
                   mutate(variable = factor(variable, levels = c("bs_diary_5", "bs_diary_13", "bs_diary_22",
                                                                 "bs_diary_15", "bs_diary_9", "bs_diary_10",
                                                                 "bs_diary_7", "bs_diary_11", "bs_diary_17",
                                                                 "bs_diary_8", "bs_diary_14", "bs_diary_16"),
                                            labels = c("down","dread nrest of day","worry",
                                                       "inadequate", "tired", "content",
                                                       "agitated", "irritated", "switch and focus",
                                                       "extremely well", "full of ideas", "thoughts are racing"))), aes(xintercept = mean_value), colour = "dark red") +
    facet_wrap(variable~., nrow = 2) +
    theme_minimal() +
    ylab(label = "Counts") +
    xlab(label = "Value")

p

ggsave(plot = p + theme(text = element_text(size = 12)),
       filename = "outputs/figures/group_mean_ppc_12.pdf", width = 9, height = 6, units = "in", dpi = 300)
ggsave(plot = p + theme(text = element_text(size = 12)),
       filename = "outputs/figures/group_mean_ppc_12.tiff", width = 9, height = 6, units = "in", dpi = 300)
ggsave(plot = p + theme(text = element_text(size = 14)),
       filename = "outputs/figures/group_mean_ppc_14.pdf", width = 9, height = 6, units = "in", dpi = 300)
ggsave(plot = p + theme(text = element_text(size = 14)),
       filename = "outputs/figures/group_mean_ppc_14.tiff", width = 9, height = 6, units = "in", dpi = 300)
ggsave(plot = p + theme(text = element_text(size = 16)),
       filename = "outputs/figures/group_mean_ppc_16.pdf", width = 9, height = 6, units = "in", dpi = 300)
ggsave(plot = p + theme(text = element_text(size = 16)),
       filename = "outputs/figures/group_mean_ppc_16.tiff", width = 9, height = 6, units = "in", dpi = 300)
ggsave(plot = p + theme(text = element_text(size = 18)),
       filename = "outputs/figures/group_mean_ppc_18.pdf", width = 9, height = 6, units = "in", dpi = 300)
ggsave(plot = p + theme(text = element_text(size = 18)),
       filename = "outputs/figures/group_mean_ppc_18.tiff", width = 9, height = 6, units = "in", dpi = 300)
ggsave(plot = p + theme(text = element_text(size = 20)),
       filename = "outputs/figures/group_mean_ppc_20.pdf", width = 9, height = 6, units = "in", dpi = 300)
ggsave(plot = p + theme(text = element_text(size = 20)),
       filename = "outputs/figures/group_mean_ppc_20.tiff", width = 9, height = 6, units = "in", dpi = 300)




## PPC2: between-subject variation

# SD
ppc_data %>%
    group_by(variable, rep) %>%
    summarise(mean_value = mean(value), sd_value = sd(value)) %>%
    mutate(variable = factor(variable, levels = c("bs_diary_5", "bs_diary_13", "bs_diary_22",
                                                  "bs_diary_15", "bs_diary_9", "bs_diary_10",
                                                  "bs_diary_7", "bs_diary_11", "bs_diary_17",
                                                  "bs_diary_8", "bs_diary_14", "bs_diary_16"),
                             labels = c("down","dread nrest of day","worry",
                                        "inadequate", "tired", "content",
                                        "agitated", "irritated", "switch and focus",
                                        "extremely well", "full of ideas", "thoughts are racing"))) %>%
    ggplot(aes(x = sd_value)) +
    geom_histogram() +
    geom_vline(data = ppc1_truth %>% mutate(variable = factor(variable, levels = c("bs_diary_5", "bs_diary_13", "bs_diary_22",
                                                                                   "bs_diary_15", "bs_diary_9", "bs_diary_10",
                                                                                   "bs_diary_7", "bs_diary_11", "bs_diary_17",
                                                                                   "bs_diary_8", "bs_diary_14", "bs_diary_16"),
                                                              labels = c("down","dread nrest of day","worry",
                                                                         "inadequate", "tired", "content",
                                                                         "agitated", "irritated", "switch and focus",
                                                                         "extremely well", "full of ideas", "thoughts are racing"))), aes(xintercept = sd_value), colour = "dark red") +
    facet_wrap(variable~., nrow = 2) +
    theme_minimal() +
    ylab(label = "Counts") +
    xlab(label = "Value")


## PPC3: ranked ids plots

# Patient means
p <- ppc_data %>%
    mutate(variable = factor(variable, levels = c("bs_diary_5", "bs_diary_13", "bs_diary_22",
                                                  "bs_diary_15", "bs_diary_9", "bs_diary_10",
                                                  "bs_diary_7", "bs_diary_11", "bs_diary_17",
                                                  "bs_diary_8", "bs_diary_14", "bs_diary_16"),
                             labels = c("down","dread nrest of day","worry",
                                        "inadequate", "tired", "content",
                                        "agitated", "irritated", "switch and focus",
                                        "extremely well", "full of ideas", "thoughts are racing"))) %>%
    group_by(rep, patient_id, variable) %>%
    summarise(mean_value = mean(value)) %>%
    group_by(rep, variable) %>%
    arrange(desc(mean_value)) %>%
    mutate(id = row_number()) %>%
    group_by(variable, rep) %>%
    arrange(id) %>%
    group_by(variable, id) %>%
    summarise(map_value = median(mean_value),
              cci_lwr = quantile(mean_value, 0.025),
              cci_upr = quantile(mean_value, 0.975)) %>%
    ggplot(aes(x = factor(id), y = map_value )) +
    geom_pointrange(aes(ymin = cci_lwr, ymax = cci_upr), alpha = 1, size = 1, fatten = 1, shape = 3) +
    geom_point(data = true_data %>%
                   mutate(variable = factor(variable,
                                            levels = c("bs_diary_5", "bs_diary_13", "bs_diary_22",
                                                       "bs_diary_15", "bs_diary_9", "bs_diary_10",
                                                       "bs_diary_7", "bs_diary_11", "bs_diary_17",
                                                       "bs_diary_8", "bs_diary_14", "bs_diary_16"),
                                            labels = c("down","dread nrest of day","worry",
                                                       "inadequate", "tired", "content",
                                                       "agitated", "irritated", "switch and focus",
                                                       "extremely well", "full of ideas", "thoughts are racing"))
                   ) %>%
                   group_by(patient_id, variable) %>%
                   summarise(mean_value = mean(value, na.rm = TRUE)) %>%
                   group_by(variable) %>%
                   arrange(desc(mean_value)) %>%
                   mutate(id = row_number()), aes(x = factor(id), y = mean_value), colour = "dark red", shape = 4) +
    facet_wrap(variable~., nrow = 2) +
    coord_flip() +
    theme_minimal() +
    ylab(label = "Value") +
    xlab(label = "Patient no.")

p

ggsave(plot = p + theme(text = element_text(size = 12)),
       filename = "outputs/figures/patient_mean_ppc_12.pdf", width = 9, height = 7, units = "in", dpi = 300)
ggsave(plot = p + theme(text = element_text(size = 12)),
       filename = "outputs/figures/patient_mean_ppc_12.tiff", width = 9, height = 7, units = "in", dpi = 300)
ggsave(plot = p + theme(text = element_text(size = 14)),
       filename = "outputs/figures/patient_mean_ppc_14.pdf", width = 9, height = 7, units = "in", dpi = 300)
ggsave(plot = p + theme(text = element_text(size = 14)),
       filename = "outputs/figures/patient_mean_ppc_14.tiff", width = 9, height = 7, units = "in", dpi = 300)
ggsave(plot = p + theme(text = element_text(size = 16)),
       filename = "outputs/figures/patient_mean_ppc_16.pdf", width = 9, height = 7, units = "in", dpi = 300)
ggsave(plot = p + theme(text = element_text(size = 16)),
       filename = "outputs/figures/patient_mean_ppc_16.tiff", width = 9, height = 7, units = "in", dpi = 300)
ggsave(plot = p + theme(text = element_text(size = 18)),
       filename = "outputs/figures/patient_mean_ppc_18.pdf", width = 9, height = 7, units = "in", dpi = 300)
ggsave(plot = p + theme(text = element_text(size = 18)),
       filename = "outputs/figures/patient_mean_ppc_18.tiff", width = 9, height = 7, units = "in", dpi = 300)
ggsave(plot = p + theme(text = element_text(size = 20)),
       filename = "outputs/figures/patient_mean_ppc_20.pdf", width = 9, height = 7, units = "in", dpi = 300)
ggsave(plot = p + theme(text = element_text(size = 20)),
       filename = "outputs/figures/patient_mean_ppc_20.tiff", width = 9, height = 7, units = "in", dpi = 300)



# Patient sd
ppc_data %>%
    mutate(variable = factor(variable, levels = c("bs_diary_5", "bs_diary_13", "bs_diary_22",
                                                  "bs_diary_15", "bs_diary_9", "bs_diary_10",
                                                  "bs_diary_7", "bs_diary_11", "bs_diary_17",
                                                  "bs_diary_8", "bs_diary_14", "bs_diary_16"),
                             labels = c("down","dread nrest of day","worry",
                                        "inadequate", "tired", "content",
                                        "agitated", "irritated", "switch and focus",
                                        "extremely well", "full of ideas", "thoughts are racing"))) %>%
    group_by(rep, patient_id, variable) %>%
    summarise(sd_value = sd(value)) %>%
    group_by(rep, variable) %>%
    arrange(desc(sd_value)) %>%
    mutate(id = row_number()) %>%
    group_by(variable, rep) %>%
    arrange(id) %>%
    group_by(variable, id) %>%
    summarise(map_value = median(sd_value),
              cci_lwr = quantile(sd_value, 0.025),
              cci_upr = quantile(sd_value, 0.975)) %>%
    ggplot(aes(x = factor(id), y = map_value )) +
    geom_pointrange(aes(ymin = cci_lwr, ymax = cci_upr), alpha = 1, size = 1, fatten = 1, shape = 3) +
    geom_point(data = true_data %>%
                   mutate(variable = factor(variable,
                                            levels = c("bs_diary_5", "bs_diary_13", "bs_diary_22",
                                                       "bs_diary_15", "bs_diary_9", "bs_diary_10",
                                                       "bs_diary_7", "bs_diary_11", "bs_diary_17",
                                                       "bs_diary_8", "bs_diary_14", "bs_diary_16"),
                                            labels = c("down","dread nrest of day","worry",
                                                       "inadequate", "tired", "content",
                                                       "agitated", "irritated", "switch and focus",
                                                       "extremely well", "full of ideas", "thoughts are racing"))
                   ) %>%
                   group_by(patient_id, variable) %>%
                   summarise(sd_value = sd(value, na.rm = TRUE)) %>%
                   group_by(variable) %>%
                   arrange(desc(sd_value)) %>%
                   mutate(id = row_number()), aes(x = factor(id), y = sd_value), colour = "dark red", shape = 4) +
    facet_wrap(variable~., nrow = 2) +
    coord_flip() +
    theme_minimal() +
    ylab(label = "Value") +
    xlab(label = "Patient no.")



















# Extract emission means
emiss_mu_bar <-  do.call(rbind, lapply(1:length(out$emiss_mu_bar), function(q){
    out$emiss_mu_bar[[q]] %>%
        as.data.frame() %>%
        mutate(iter = row_number(),
               dep = names(out$emiss_mu_bar)[q]) %>%
        gather(mu, value, -iter, -dep)
}))

# Add labels for the constructs

dep_vars <- c("bs_diary_5", "bs_diary_13", "bs_diary_22",
              "bs_diary_15", "bs_diary_9", "bs_diary_10","bs_qids_tot")
man_vars <- c("bs_diary_7", "bs_diary_11", "bs_diary_17",
              "bs_diary_8", "bs_diary_14", "bs_diary_16","bs_asrm_tot")

# Plot

# Patient-specific values instead of group-level MAP
subj_emiss <- do.call(rbind, lapply(1:length(out$PD_subj),function(s){
    out$PD_subj[[s]] %>%
        as.data.frame() %>%
        dplyr::select("dep1_mu_S1":"dep14_mu_S4") %>%
        mutate(patient_id = s, iter = row_number()) %>%
        filter(iter > burn_in) %>%
        dplyr::select(-iter) %>%
        gather(dep, value, -patient_id) %>%
        group_by(patient_id, dep) %>%
        summarise(value = median(value)) %>%
        separate(dep, into = c("dep","mu"), sep = "_mu_") %>%
        mutate(dep = factor(dep, levels = paste0("dep",1:14), labels = c("bs_diary_5", "bs_diary_13", "bs_diary_22",
                                                                         "bs_diary_15", "bs_diary_9", "bs_diary_10", "bs_qids_tot",
                                                                         "bs_diary_7", "bs_diary_11", "bs_diary_17",
                                                                         "bs_diary_8", "bs_diary_14", "bs_diary_16", "bs_asrm_tot")))

} ))

# Subject and Group level parameters
subj_emiss_nice <- subj_emiss %>%
    mutate(construct = case_when(dep %in% dep_vars ~ "Depression items",
                                 dep %in% man_vars ~ "Mania items"),
           construct = factor(construct, levels = c("Depression items","Mania items")),
           dep = factor(dep,
                        levels = c("bs_diary_5", "bs_diary_13", "bs_diary_22",
                                   "bs_diary_15", "bs_diary_9", "bs_diary_10", "bs_qids_tot",
                                   "bs_diary_7", "bs_diary_11", "bs_diary_17",
                                   "bs_diary_8", "bs_diary_14", "bs_diary_16", "bs_asrm_tot"),
                        labels = c("down","dread rest of day","worry",
                                   "inadequate", "tired", "content","QIDS",
                                   "agitated", "irritated", "switch and focus",
                                   "extremely well", "full of ideas", "racing thoughts","ASRM")),
           mu = factor(mu, levels = paste0("S",1:4), labels = c("neutral","elevated","mixed","lowered"))) %>%
    rename("state" = "mu")

group_emiss_nice <- emiss_mu_bar %>%
    filter(iter > burn_in) %>%
    group_by(dep, mu) %>%
    summarise(group_mean = mean(value),
              group_median = median(value),
              group_sd = sd(value),
              group_CCIlwr = quantile(value, 0.025),
              group_CCIupr = quantile(value, 0.975)) %>%
    mutate(construct = case_when(dep %in% dep_vars ~ "Depression items",
                                 dep %in% man_vars ~ "Mania items"),
           construct = factor(construct, levels = c("Depression items","Mania items")),
           dep = factor(dep,
                        levels = c("bs_diary_5", "bs_diary_13", "bs_diary_22",
                                   "bs_diary_15", "bs_diary_9", "bs_diary_10", "bs_qids_tot",
                                   "bs_diary_7", "bs_diary_11", "bs_diary_17",
                                   "bs_diary_8", "bs_diary_14", "bs_diary_16", "bs_asrm_tot"),
                        labels = c("down","dread rest of day","worry",
                                   "inadequate", "tired", "content","QIDS",
                                   "agitated", "irritated", "switch and focus",
                                   "extremely well", "full of ideas", "racing thoughts","ASRM")),
           mu = factor(mu, levels = paste0("mu_",1:4), labels = c("neutral","elevated","mixed","lowered"))) %>%
    rename("state" = "mu")

p <- ggplot(data = subj_emiss_nice) +
    geom_line(aes(x = state, y = value, group = patient_id), alpha = 0.2) +
    geom_jitter(aes(x = state, y = value, colour = state), width = 0.15, alpha = 0.5) +
    geom_pointrange(data = group_emiss_nice,
                    aes(x = state, y = group_mean,
                        ymin = group_CCIlwr, ymax = group_CCIupr), size = 0.25) +
    scale_color_viridis_d(option = "plasma", direction = -1) +
    geom_hline(yintercept = 0) +
    facet_wrap(construct+dep~., scales = "fixed", nrow = 2) +
    theme_minimal() +
    ylab(label = "EMA score") +
    xlab(label = "") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    theme(legend.position = "none") +
    # ggtitle(label = "Composition of momentary mood states: group-level and\npatient-specific item emission scores by state") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(strip.text.y = element_text(angle = 0))

p




#------------------------------------------------------------------------------#

# Decode states
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

states <- local_decoding(out_medhmm)

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

train_df

data <- as.data.frame(train_df) %>%
    rename("subj" = "patient_id")
split_obs <- split_data(data, 0.8)
train_data <- split_obs$fit
test_data <- split_obs$forecast

# Make predictions
forecast_results_medhmm <- forecast_medHMM1(s_data = train_data, object = out_medhmm, forecast_steps = 500, return_all = TRUE, Mx = 224)

# Ensure forecast_results has the necessary columns
forecast_results_medhmm <- as.data.frame(forecast_results_medhmm) %>%
    group_by(subj) %>%
    mutate(occasion = row_number()) %>%
    mutate(state = as.integer(state))


# Plot state decoding for a few individuals
model_states <- inner_join(states %>%
                               as.data.frame() %>%
                               select(subj, state, occasion) %>%
                               rename("true_state" = "state"),
                           forecast_results_medhmm %>%
                               as.data.frame() %>%
                               rename("pred_state" = "state") %>%
                               gather(state_prob, value, -c(subj, pred_state, horizon, occasion))) %>%
    mutate(pred_state = factor(pred_state),
           true_state = factor(true_state),
           subj = factor(subj),
           pred_state = case_when(horizon == 0 ~ NA,
                                  horizon != 0 ~ pred_state)) %>%
    gather(model, state, -c(subj, horizon, occasion, state_prob, value)) %>%
    group_by(subj) %>%
    mutate(thresh = max(occasion)*0.8)

ggplot(data = model_states,
       aes(x = occasion, y = model, fill = state)) +
    geom_tile(height = 0.9) +
    geom_vline(data = model_states, aes(xintercept = thresh), linetype = "dashed", color = "black") +
    scale_fill_viridis_d(option = "plasma") +
    facet_wrap(subj~., ncol = 4, scales = "free_x") +
    theme_minimal() +
    theme(legend.position = "bottom")


# Add forward probabilities
forecast_results_medhmm <- forecast_medHMM2(s_data = train_data, object = out_medhmm, forecast_steps = 500, return_all = TRUE, Mx = 224)
forecast_results_medhmm


# Look at predicted observations


#------------------------------------------------------------------------------#
forecast_results_medhmm <- forecast_medHMM3(s_data = train_df, object = out_medhmm, initial_window = 600, Mx = 56, shift = 1)

h_step <- 8
forecast_results_medhmm <- forecast_medHMM3(s_data = train_df, object = out_medhmm, initial_window = 600, Mx = 56, shift = 1, h_step = h_step)

# forecast_results_medhmm_h1 <- forecast_results_medhmm

# Plot state decoding for a few individuals
model_states <- left_join(states %>%
                              as.data.frame() %>%
                              select(subj, state, occasion) %>%
                              mutate(baseline_state = ifelse(row_number() == 1, state, lag(state, 1))) %>%
                              mutate(baseline_state = ifelse((row_number()) %% h_step == 1, baseline_state, NA)) %>%
                              fill(baseline_state, .direction = "down") %>%
                              mutate(baseline_state = case_when(
                                  occasion <= 600 ~ NA,
                                  occasion > 600 ~ baseline_state
                              )) %>%
                              rename("true_state" = "state"),
                          forecast_results_medhmm %>%
                              as.data.frame() %>%
                              rename("pred_state" = "state",
                                     "occasion" = "horizon") %>%
                              gather(state_prob, value, -c(subj, pred_state, occasion))) %>%
    mutate(pred_state = factor(pred_state),
           baseline_state = factor(baseline_state),
           true_state = factor(true_state),
           subj = factor(subj),
           pred_state = case_when(occasion == 0 ~ NA,
                                  occasion != 0 ~ pred_state)) %>%
    gather(model, state, -c(subj, occasion, state_prob, value)) %>%
    group_by(subj) %>%
    mutate(thresh = 600)

ggplot(data = model_states,
       aes(x = occasion, y = model, fill = state)) +
    geom_tile(height = 0.9) +
    geom_vline(data = model_states, aes(xintercept = thresh), linetype = "dashed", color = "black") +
    scale_fill_viridis_d(option = "plasma") +
    facet_wrap(subj~., ncol = 4, scales = "free_x") +
    theme_minimal() +
    theme(legend.position = "bottom") +
    ggtitle("One step ahead prediction")

s_data = train_df
object = out_medhmm
initial_window = 600
value_range = NULL
burn_in = NULL
Mx = 56
shift = 1
return_all = FALSE
return_full_sequence = FALSE
h_step = 5









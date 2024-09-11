
library(mHMMbayes)
library(medHMM)


# Load file:
file_path <- "/Users/a6159737/Documents/Utrecht University/PhD/Projects/Simulation studies/medhmm-sim/seb/main_sim_fat/outputs/job8/mhmm_and_medhhmm_res_dur_s_5ind_s_2occ_s_3dep_s_3_iter4000_rep107.rds"

# Load test output:
out <- readRDS(file_path)

# Sanity checks (decoding on train data):
out$MHMM$execution_time
out$MEDHMM$execution_time

mean(out$MHMM$state_decoding[,2] == out$MHMM$state_decoding[,3])
mean(out$MEDHMM$state_decoding[,2] == out$MEDHMM$state_decoding[,3])

out$MEDHMM$MAP$dwell_mu_bar
out$MEDHMM$MAP$dwell_varmu_bar

out$MEDHMM$MAP

out$MEDHMM$MAP$emiss_mu_bar[[1]]
# out$MEDHMM$MAP$emiss_mu_bar[[2]]

# Test function:
out$sim_data






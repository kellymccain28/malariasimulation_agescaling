suppressPackageStartupMessages(library(ggplot2))
library(malariasimulation)
library(malariaEquilibrium)
library(cowplot)
library(reshape2)
library(ggplot2)
library(gridExtra)

year <- 365
month <- 30 # ish
sim_length <- 10 * year
human_population <- 25000
starting_EIR <- 70

#One way to link seasonality to admin units- using input from deterministic package
seas <- read.csv('C:\\Users\\jchallen\\Dropbox\\NC_data_code\\deterministic-malaria-model-master-NC\\data-raw\\admin_units_seasonal.csv')
lst <- as.numeric(seas[seas$NAME_0=='Burkina Faso' & seas$NAME_1=='Sahel',5:11])

simparams <- get_parameters(list(
  human_population = human_population,
  model_seasonality = T, 
  individual_mosquitoes = F,
  phi_bednets = 0.79, #proportion of bites taken in bed. Check with Ellie's review!
  species = 1,
  g0 = lst[1],
  g = c(lst[2], lst[4], lst[6]), 
  h = c(lst[3], lst[5], lst[7]),
  #tbv_md = 0.00028045,
  #tbv_mt = 0.00028045,
  #tbv_ma = 0.00028045,
  #tbv_mu = 0.00028045,
  prevalence_rendering_min_ages = 0 * year,
  prevalence_rendering_max_ages = 100 * year, 
  clinical_incidence_rendering_min_ages = 0,
  clinical_incidence_rendering_max_ages = 100 * year,
  incidence_rendering_min_ages = 0,
  incidence_rendering_max_ages = 100 * year
)
)

simparams <- set_drugs(simparams, list(AL_params, SP_AQ_params))
simparams <- set_clinical_treatment(simparams, 1, 1, 0.5 )

simparams <- set_equilibrium(simparams, starting_EIR)

out1 <- run_simulation(sim_length, simparams)

ggplot(out1) + geom_line(aes(x = timestep, y = n_detect_lm_0_36500/ n_age_0_36500)) + 
  theme_classic()

peak <- peak_season_offset(simparams)
peak

######################### 2. PEV only (R21) ##########################

# Mass PEV
R21params <- set_mass_pev(
  simparams,
  profile = r21_profile,
  coverages = c(0.9),
  min_wait = 0,
  min_ages = 5*year, # The minimum age for the target population to be vaccinated. 
  max_ages = 63 * year, # The maximum age for the target population to be vaccinated.
  timesteps = c(7)*year + peak - 3.5*month,
  adult_scaling = 1,
  booster_spacing = c(1*year), # The booster is given 10 years after?
  booster_coverage = matrix(rep(0.95,1), nrow = 1, ncol= 1),#, # Coverage of the booster dose is 5%.
  booster_profile = list(r21_booster_profile)
)

out2 <- run_simulation(sim_length, R21params)

R21params_A <- set_mass_pev(
  simparams,
  profile = r21_profile,
  coverages = c(0.9),
  min_wait = 0,
  min_ages = 5*year, # The minimum age for the target population to be vaccinated. 
  max_ages = 63 * year, # The maximum age for the target population to be vaccinated.
  timesteps = c(7)*year + peak - 3.5*month,
  adult_scaling = 0.6,
  booster_spacing = c(1*year), # The booster is given 10 years after?
  booster_coverage = matrix(rep(0.95,1), nrow = 1, ncol= 1),#, # Coverage of the booster dose is 5%.
  booster_profile = list(r21_booster_profile)
)

out2A <- run_simulation(sim_length, R21params_A)

R21params_B <- set_mass_pev(
  simparams,
  profile = r21_profile,
  coverages = c(0.9),
  min_wait = 0,
  min_ages = 5*year, # The minimum age for the target population to be vaccinated. 
  max_ages = 63 * year, # The maximum age for the target population to be vaccinated.
  timesteps = c(7)*year + peak - 3.5*month,
  adult_scaling = 0.4,
  booster_spacing = c(1*year), # The booster is given 10 years after?
  booster_coverage = matrix(rep(0.95,1), nrow = 1, ncol= 1),#, # Coverage of the booster dose is 5%.
  booster_profile = list(r21_booster_profile)
)

out2B <- run_simulation(sim_length, R21params_B)

R21params_C <- set_mass_pev(
  simparams,
  profile = r21_profile,
  coverages = c(0.9),
  min_wait = 0,
  min_ages = 5*year, # The minimum age for the target population to be vaccinated. 
  max_ages = 63 * year, # The maximum age for the target population to be vaccinated.
  timesteps = c(7)*year + peak - 3.5*month,
  adult_scaling = 0.2,
  booster_spacing = c(1*year), # The booster is given 10 years after?
  booster_coverage = matrix(rep(0.95,1), nrow = 1, ncol= 1),#, # Coverage of the booster dose is 5%.
  booster_profile = list(r21_booster_profile)
)

out2C <- run_simulation(sim_length, R21params_C)

df2 <- data.frame('t' = out1$timestep/365, 
                  'pr1'=out1$n_detect_lm_0_36500/out1$n_age_0_36500,
                  'pr2'=out2$n_detect_lm_0_36500/out2$n_age_0_36500,
                  'pr2A'=out2A$n_detect_lm_0_36500/out2A$n_age_0_36500,
                  'pr2B'=out2B$n_detect_lm_0_36500/out2B$n_age_0_36500,
                  'pr2C'=out2C$n_detect_lm_0_36500/out2C$n_age_0_36500)
ggplot(df2) + geom_line(aes(x = t, y = pr1)) + theme_classic() + 
  geom_line(aes(x = t, y = pr2), color = 'purple') + xlim(c(7,10)) + 
  geom_line(aes(x = t, y = pr2A), color = 'magenta') +
  geom_line(aes(x = t, y = pr2B), color = 'orange') +
  geom_line(aes(x = t, y = pr2C), color = 'darkgreen') 

df2i <- data.frame('t' = out1$timestep/365, 
                  'inc1' = out1$n_inc_clinical_0_36500/out1$n_age_0_36500,
                  'inc2' = out2$n_inc_clinical_0_36500/out2$n_age_0_36500,
                  'inc2A' = out2A$n_inc_clinical_0_36500/out2A$n_age_0_36500,
                  'inc2B' = out2B$n_inc_clinical_0_36500/out2B$n_age_0_36500,
                  'inc2C' = out2C$n_inc_clinical_0_36500/out2C$n_age_0_36500)

ggplot(df2i) + geom_line(aes(x = t, y = inc1)) + theme_classic() + 
  geom_line(aes(x = t, y = inc2), color = 'purple') + xlim(c(7,10)) + 
  geom_line(aes(x = t, y = inc2A), color = 'magenta') +
  geom_line(aes(x = t, y = inc2B), color = 'orange') +
  geom_line(aes(x = t, y = inc2C), color = 'darkgreen')

#cumulative inc
tt <- R21params$mass_pev_timesteps + 90

ci <- out1$n_inc_clinical_0_36500[(tt):sim_length]
ci2 <- out2$n_inc_clinical_0_36500[(tt):sim_length]
ci2A <- out2A$n_inc_clinical_0_36500[(tt):sim_length]
ci2B <- out2B$n_inc_clinical_0_36500[(tt):sim_length]
ci2C <- out2C$n_inc_clinical_0_36500[(tt):sim_length]

for(i in 2:(length(ci))){
  ci[i] <- ci[i] + ci[i-1]
  ci2[i] <- ci2[i] + ci2[i-1]
  ci2A[i] <- ci2A[i] + ci2A[i-1]
  ci2B[i] <- ci2B[i] + ci2B[i-1]
  ci2C[i] <- ci2C[i] + ci2C[i-1]
}
gh <- data.frame('t' = seq((tt):sim_length)/365, 
                 'ci' = ci, 'ci2' = ci2, 'ci2A' = ci2A, 'ci2B' = ci2B, 'ci2C' = ci2C)#, 'prop' = (1-ci3/ci2))
ggplot(gh) + geom_line(aes(x=t, y=ci)) + 
  geom_line(aes(x=t, y=ci2), color = 'purple') + theme_classic() + 
  geom_line(aes(x=t, y=ci2A),color = 'magenta') +  
  geom_line(aes(x=t, y=ci2B),color = 'orange') +
  geom_line(aes(x=t, y=ci2C),color = 'darkgreen') 

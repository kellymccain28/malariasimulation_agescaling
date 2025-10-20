#suppressPackageStartupMessages(library(ggplot2))
# devtools::install_github('https://github.com/kellymccain28/malariasimulation_agescaling.git')
library(malariasimulation)
library(malariaEquilibrium)
library(cowplot)
library(reshape2)
library(ggplot2)
library(gridExtra)

year <- 365
month <- 30 # ish
sim_length <- 10 * year
human_population <- 7000
starting_EIR <- 70

#One way to link seasonality to admin units- using input from deterministic package
#seas <- read.csv('C:\\Users\\jchallen\\Dropbox\\NC_data_code\\deterministic-malaria-model-master-NC\\data-raw\\admin_units_seasonal.csv')
#lst <- as.numeric(seas[seas$NAME_0=='Burkina Faso' & seas$NAME_1=='Sahel',5:11])

lst <- c(0.2852770,-0.0248801,-0.0529426,-0.0168910,-0.0216681,-0.0242904,-0.0073646)

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

peak <- peak_season_offset(simparams)
peak
# simparams <- set_mass_pev(
#   simparams,
#   profile = r21_profile,
#   coverages = c(0),
#   min_wait = 0,
#   min_ages = 6*month, #
#   max_ages = 63 * year, #
#   timesteps = c(7)*year + peak - 3.5*month,
#   adult_scaling = 1,
#   adolesc_scaling = 1,
#   u5_scaling = 1, 
#   booster_spacing = c(1*year), 
#   booster_coverage = matrix(rep(0.95,1), nrow = 1, ncol= 1),
#   booster_profile = list(r21_booster_profile)
# )
# 
# 
# out1 <- malariasimulation::run_simulation(sim_length, simparams)
# 
# ggplot(out1) + geom_line(aes(x = timestep, y = n_detect_lm_0_36500/ n_age_0_36500)) + 
#   theme_classic()

######################### 2. PEV only (R21) ##########################

# Mass PEV
R21params <- set_mass_pev(
  simparams,
  profile = r21_profile,
  coverages = c(0.9),
  min_wait = 0,
  min_ages = 6*month, #
  max_ages = 63 * year, #
  timesteps = c(3)*year + peak - 3.5*month,
  adult_scaling = 1,
  adolesc_scaling = 1,
  booster_spacing = c(1*year), 
  booster_coverage = matrix(rep(0.95,1), nrow = 1, ncol= 1),
  booster_profile = list(r21_booster_profile)
)

out2 <- run_simulation(sim_length, R21params)

R21params_A <- set_mass_pev(
  simparams,
  profile = r21_profile,
  coverages = c(0.9),
  min_wait = 0,
  min_ages = 6*month, # The minimum age for the target population to be vaccinated. 
  max_ages = 63 * year, # The maximum age for the target population to be vaccinated.
  timesteps = c(3)*year + peak - 3.5*month,
  adult_scaling = 1,
  adolesc_scaling = 0.64,
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
  min_ages = 6*month, # The minimum age for the target population to be vaccinated. 
  max_ages = 63 * year, # The maximum age for the target population to be vaccinated.
  timesteps = c(3)*year + peak - 3.5*month,
  adult_scaling = 1,
  adolesc_scaling = 0.4,
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
  min_ages = 6*month, # The minimum age for the target population to be vaccinated. 
  max_ages = 63 * year, # The maximum age for the target population to be vaccinated.
  timesteps = c(3)*year + peak - 3.5*month,
  adult_scaling = 1,
  adolesc_scaling = 0.2,
  booster_spacing = c(1*year), # The booster is given 10 years after?
  booster_coverage = matrix(rep(0.95,1), nrow = 1, ncol= 1),#, # Coverage of the booster dose is 5%.
  booster_profile = list(r21_booster_profile)
)

out2C <- run_simulation(sim_length, R21params_C)




# Test vaccinating school-aged children with scaling dependent on age group 
R21params_D <- set_mass_pev(
  simparams,
  profile = r21_profile,
  coverages = c(0.9),
  min_wait = 0,
  min_ages = 6 * month, # The minimum age for the target population to be vaccinated. 
  max_ages = 15 * year, # The maximum age for the target population to be vaccinated.
  timesteps = c(3)*year + peak - 3.5*month,
  adult_scaling = 1,
  adolesc_scaling = 1,
  booster_spacing = c(1*year), # The booster is given 1 year after
  booster_coverage = matrix(rep(0.95,1), nrow = 1, ncol= 1),#, # Coverage of the booster dose is 95%.
  booster_profile = list(r21_booster_profile)
)

out2D <- run_simulation(sim_length, R21params_D)

R21params_E <- set_mass_pev(
  simparams,
  profile = r21_profile,
  coverages = c(0.9),
  min_wait = 0,
  min_ages = 6 * month, # The minimum age for the target population to be vaccinated. 
  max_ages = 15 * year, # The maximum age for the target population to be vaccinated.
  timesteps = c(3)*year + peak - 3.5*month,
  adult_scaling = 1,
  adolesc_scaling = 0.64,
  booster_spacing = c(1*year), # The booster is given 1 year after
  booster_coverage = matrix(rep(0.95,1), nrow = 1, ncol= 1),#, # Coverage of the booster dose is 95%.
  booster_profile = list(r21_booster_profile)
)

out2E <- run_simulation(sim_length, R21params_E)

R21params_F <- set_pev_epi(
  simparams,
  profile = r21_profile,
  coverages = c(0.9),
  min_wait = 0,
  age = 6 * month, 
  timesteps = c(3)*year,
  adolesc_scaling = 1,
  booster_spacing = c(7*year), # The booster is given 1 and 5 years after
  booster_coverage = matrix(rep(0.95,1), nrow = 1, ncol= 1),#, # Coverage of the booster dose is 95%.
  booster_profile = list(r21_booster_profile)
)

out2F <- run_simulation(sim_length, R21params_F)

R21params_G <- set_pev_epi(
  simparams,
  profile = r21_profile,
  coverages = c(0.9),
  min_wait = 0,
  age = 6 * month, 
  timesteps = c(3)*year,
  adolesc_scaling = 0.05,
  booster_spacing = c(7*year), # The booster is given 1 and 5 years after
  booster_coverage = matrix(rep(0.95,1), nrow = 1, ncol= 1),#, # Coverage of the booster dose is 95%.
  booster_profile = list(r21_booster_profile)
)

out2G <- run_simulation(sim_length, R21params_G)

R21params_u5A <- set_mass_pev(
  simparams,
  profile = r21_profile,
  coverages = c(0.9),
  min_wait = 0,
  min_ages = 6 * month, # The minimum age for the target population to be vaccinated. 
  max_ages = 5 * year, # The maximum age for the target population to be vaccinated.
  timesteps = c(3)*year + peak - 3.5*month,
  u5_scaling = 1,
  booster_spacing = c(1*year), # The booster is given 1 year after
  booster_coverage = matrix(rep(0.95,1), nrow = 1, ncol= 1),#, # Coverage of the booster dose is 95%.
  booster_profile = list(r21_booster_profile)
)

outu5A <- run_simulation(sim_length, R21params_u5A)

R21params_u5b <- set_mass_pev(
  simparams,
  profile = r21_profile,
  coverages = c(0.9),
  min_wait = 0,
  min_ages = 6 * month, # The minimum age for the target population to be vaccinated. 
  max_ages = 5 * year, # The maximum age for the target population to be vaccinated.
  timesteps = c(3)*year + peak - 3.5*month,
  u5_scaling = 0.1,
  booster_spacing = c(1*year), # The booster is given 1 year after
  booster_coverage = matrix(rep(0.95,1), nrow = 1, ncol= 1),#, # Coverage of the booster dose is 95%.
  booster_profile = list(r21_booster_profile)
)

outu5b <- run_simulation(sim_length, R21params_u5b)


df4 <- data.frame('t' = out1$timestep/365, 
                  'pr1'=out1$n_detect_lm_0_36500/out1$n_age_0_36500, # no mass vaccination
                  'pr2'=out2$n_detect_lm_0_36500/out2$n_age_0_36500, # normal mass vaccination to adults
                  'pr2A'=out2A$n_detect_lm_0_36500/out2A$n_age_0_36500, # adult - 0.6 scaling
                  'pr2B'=out2B$n_detect_lm_0_36500/out2B$n_age_0_36500, # adult - 0.4 scaling
                  'pr2C'=out2C$n_detect_lm_0_36500/out2C$n_age_0_36500, # adult - 0.2 scaling
                  'pr2D'=out2D$n_detect_lm_0_36500/out2D$n_age_0_36500, # no age scaling for ados - mass
                  'pr2E'=out2E$n_detect_lm_0_36500/out2E$n_age_0_36500, # age scaling for ados - mass
                  'pr2F'=out2F$n_detect_lm_0_36500/out2F$n_age_0_36500, # no age scaling for ados - epi
                  'pr2G'=out2G$n_detect_lm_0_36500/out2G$n_age_0_36500, # age scaling for ados - epi
                  'pru5a'=outu5A$n_detect_lm_0_36500/outu5A$n_age_0_36500, # no age scaling for children - mass
                  'pru5b'=outu5b$n_detect_lm_0_36500/outu5b$n_age_0_36500 # age scaling for childre - mass
) 
ggplot(df4) +
  geom_line(aes(x = t, y = pru5a, color = 'no u5 scaling'))+
  geom_line(aes(x = t, y = pru5b, color = 'u5 scaling, 0.1'))+
  geom_line(aes(x=t, y=pr2, color = 'mass adults, no scaling'), linewidth= 1) +  
  geom_line(aes(x=t, y=pr2A, color = 'mass adults, 0.64'), linewidth= 1) 

ggplot(df4) + #geom_line(aes(x=t, y=ci)) + 
  geom_line(aes(x=t, y=pr1, color = 'no vac'), linewidth= 1) + theme_classic() + 
  geom_line(aes(x=t, y=pr2, color = 'mass adults, no scaling'), linewidth= 1) +  
  geom_line(aes(x=t, y=pr2A, color = 'mass adults, 0.64'), linewidth= 1) +
  # geom_line(aes(x=t, y=pr2B),color = 'black') +
  # geom_line(aes(x=t, y=pr2C),color = 'aquamarine') +
  geom_line(aes(x=t, y=pr2D, color = 'mass 5-15y, no scaling'), linewidth= 1) +
  geom_line(aes(x=t, y=pr2E, color = 'mass 5-15y, 0.64'), linewidth= 1) +
  geom_line(aes(x=t, y=pr2F, color = 'EPI, no scaling'), linewidth= 1) +
  geom_line(aes(x=t, y=pr2G, color = 'EPI, 0.05'), linewidth= 1) +
  scale_color_manual(values = c('purple','magenta','brown','darkgoldenrod','darkgreen','darkblue','darkorange')) + labs(color = '')
ggsave("prevalence_withtransformation.png")

df2i <- data.frame('t' = out1$timestep/365, 
                   'inc1' = out1$n_inc_clinical_0_36500/out1$n_age_0_36500,
                   'inc2' = out2$n_inc_clinical_0_36500/out2$n_age_0_36500,
                   'inc2A' = out2A$n_inc_clinical_0_36500/out2A$n_age_0_36500,
                   'inc2B' = out2B$n_inc_clinical_0_36500/out2B$n_age_0_36500,
                   'inc2C' = out2C$n_inc_clinical_0_36500/out2C$n_age_0_36500,
                   'inc2D'=out2D$n_inc_clinical_0_36500/out2D$n_age_0_36500, # no age scaling for ados - mass
                   'inc2E'=out2E$n_inc_clinical_0_36500/out2E$n_age_0_36500, # age scaling for ados - mass
                   'inc2F'=out2F$n_inc_clinical_0_36500/out2F$n_age_0_36500, # no age scaling for ados - epi 
                   'inc2G'=out2G$n_inc_clinical_0_36500/out2G$n_age_0_36500) # age scaling for ados - epi )

ggplot(df2i) + geom_line(aes(x = t, y = inc1)) + theme_classic() + 
  geom_line(aes(x = t, y = inc2), color = 'purple') + xlim(c(7,10)) + 
  geom_line(aes(x = t, y = inc2A), color = 'magenta') +
  # geom_line(aes(x=t, y=inc2A),color = 'brown') +
  # geom_line(aes(x=t, y=inc2B),color = 'black') +
  # geom_line(aes(x=t, y=inc2C),color = 'aquamarine') +
  geom_line(aes(x=t, y=inc2D),color = 'darkgoldenrod1') +
  geom_line(aes(x=t, y=inc2E),color = 'darkgreen') +
  geom_line(aes(x=t, y=inc2F),color = 'darkblue') +
  geom_line(aes(x=t, y=inc2G),color = 'darkorange') 

#cumulative inc
tt <- R21params$mass_pev_timesteps + 90

ci <- out1$n_inc_clinical_0_36500[(tt):sim_length]
ci2 <- out2$n_inc_clinical_0_36500[(tt):sim_length]
ci2A <- out2A$n_inc_clinical_0_36500[(tt):sim_length]
ci2B <- out2B$n_inc_clinical_0_36500[(tt):sim_length]
ci2C <- out2C$n_inc_clinical_0_36500[(tt):sim_length]
ci2D <- out2D$n_inc_clinical_0_36500[(tt):sim_length]
ci2E <- out2E$n_inc_clinical_0_36500[(tt):sim_length]
ci2F <- out2F$n_inc_clinical_0_36500[(tt):sim_length]
ci2G <- out2G$n_inc_clinical_0_36500[(tt):sim_length]

for(i in 2:(length(ci))){
  ci[i] <- ci[i] + ci[i-1]
  ci2[i] <- ci2[i] + ci2[i-1]
  ci2A[i] <- ci2A[i] + ci2A[i-1]
  ci2B[i] <- ci2B[i] + ci2B[i-1]
  ci2C[i] <- ci2C[i] + ci2C[i-1]
  ci2D[i] <- ci2D[i] + ci2D[i-1]
  ci2E[i] <- ci2E[i] + ci2E[i-1]
  ci2F[i] <- ci2F[i] + ci2F[i-1]
  ci2G[i] <- ci2G[i] + ci2G[i-1]
}
gh <- data.frame('t' = seq((tt):sim_length)/365, 
                 'ci' = ci, 'ci2' = ci2, 'ci2A' = ci2A, 'ci2B' = ci2B, 'ci2C' = ci2C, 'ci2D' = ci2D,
                 'ci2E' = ci2E, 'ci2F' = ci2F, 'ci2G' =ci2G)#, 'prop' = (1-ci3/ci2))
ggplot(gh) + geom_line(aes(x=t, y=ci)) + 
  geom_line(aes(x=t, y=ci2), color = 'purple') + theme_classic() + 
  geom_line(aes(x = t, y = ci2A), color = 'magenta') +
  # geom_line(aes(x=t, y=ci2A),color = 'brown') +
  # geom_line(aes(x=t, y=ci2B),color = 'black') +
  # geom_line(aes(x=t, y=ci2C),color = 'aquamarine') +
  geom_line(aes(x=t, y=ci2D),color = 'darkgoldenrod1') +
  geom_line(aes(x=t, y=ci2E),color = 'darkgreen') +
  geom_line(aes(x=t, y=ci2F),color = 'darkblue') +
  geom_line(aes(x=t, y=ci2G),color = 'darkorange') 


# Check on number of vaccine doses per strategy 
df <- dplyr::bind_rows(out1 %>% mutate(strategy = 'no vaccine'),
                       out2 %>% mutate(strategy = 'mass - 6m-63y'),
                       out2A %>% mutate(strategy = 'mass - 6m-63y - 0.6'),
                       out2B %>% mutate(strategy = 'mass - 6m-63y - 0.4'),
                       out2C %>% mutate(strategy = 'mass - 6m-63y - 0.2'),
                       out2D %>% mutate(strategy = 'mass - 6m-15y'),
                       out2E %>% mutate(strategy = 'mass - 6m-15y - 0.64'),
                       out2F %>% mutate(strategy = 'EPI'),
                       out2G %>% mutate(strategy = 'EPI - 0.4'))
ggplot(df) + 
  geom_col(aes(x = strategy, y = n_pev_mass_dose_1))
ggplot(df) + 
  geom_col(aes(x = strategy, y = n_pev_epi_dose_1))


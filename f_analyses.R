library(tidyverse)
library(emmeans)
library(tidybayes)
library(brms)
library(bayesplot)
library(patchwork)
library(corrplot)
library(gtsummary)
library(gt)
library(rstanarm)
library(broom.mixed)

# Set CI to 0.90
emm_options(summary = list(level = 0.90))
theme_set(theme_classic())
color_scheme_set(scheme = 'brightblue')

# brms seed for reproducibility
my_seed <- 2023

d_saliva <- read_rds('data/d_f_saliva.rds')
d_saliva_wide <- read_rds('data/d_f_saliva_wide.rds')
h_saliva <- read_rds('data/h_f_saliva.rds')
h_saliva_wide <- read_rds('data/h_f_saliva_wide.rds')

d_saliva <- ungroup(d_saliva); d_saliva_wide <- ungroup(d_saliva_wide)
h_saliva <- ungroup(h_saliva); h_saliva_wide <- ungroup(h_saliva_wide)

behav <- read_rds('data/behavior_for_analysis.rds')

b_prior <- prior(normal(0,1), class = 'b')

# read in and tweak demographics ----
d_demos <- readxl::read_xlsx('data/demographics.xlsx')
d_demos <- select(d_demos, -(Spayed_Neutered:breed_category))
h_demos <- readxl::read_xlsx('data/demographics.xlsx', 2)
h_demos <- select(h_demos, -(Study:Ethnicity))

h_demos <- rename(h_demos, id = 'ID', sex = 'Child_Sex', age = 'Child_Age_yrs')
h_demos$id <- sprintf("%02d", h_demos$id)
d_demos <- rename(d_demos, id = 'ID', sex = 'Dog_Sex', weight = 'Dog_Approx_Weight_lbs', age = 'Dog_Age_yrs')
d_demos$id <- sprintf("%02d", d_demos$id)

# Make demo variables factor
h_demos$sex <- if_else(h_demos$sex == "F", "female", "male")
h_demos$sex <- as.factor(h_demos$sex)
table(h_demos$sex)

d_demos$sex <- if_else(d_demos$sex == "F", "female", "male")
d_demos$sex <- as.factor(d_demos$sex)
table(d_demos$sex)

# add a row for sisu and do the standardization in the demos file; then pull out sisu's information to plug in later.
add_sisu <- slice_head(d_demos, n = 1)
add_sisu$id <- "00"; add_sisu$sex <- "female"; add_sisu$weight <- 50; add_sisu$age <- 10
d_demos <- bind_rows(d_demos, add_sisu)
d_demos$sex <- as.factor(d_demos$sex)

# standardize age and weight for dog demographics
d_demos <- mutate(d_demos, weight = as.vector(scale(weight)), age = as.vector(scale(age)))
sisu_info <- slice_tail(d_demos, n = 1)
d_demos <- slice(d_demos, -nrow(d_demos)) # remove sisu

# add demographics to data
h_saliva <- left_join(h_saliva, h_demos)
h_saliva_wide <- left_join(h_saliva_wide, h_demos)

# add behavior to human data
h_saliva <- left_join(h_saliva, behav)
h_saliva_wide <- left_join(h_saliva_wide, behav)

# standardize the age variable
h_saliva$age <- as.vector(scale(h_saliva$age))
h_saliva_wide$age <- as.vector(scale(h_saliva_wide$age))

d_saliva <- left_join(d_saliva, d_demos)
d_saliva_wide <- left_join(d_saliva_wide, d_demos)

# add behavior to dog data
d_saliva <- left_join(d_saliva, behav)
d_saliva_wide <- left_join(d_saliva_wide, behav)

# add condition order to human data
orders <- read_rds(file = 'data/condition_order.rds')
orders <- select(orders, id, visit_1)
h_saliva <- left_join(h_saliva, orders)
h_saliva_wide <- left_join(h_saliva_wide, orders)

# Dog saliva ----
## comparison of UD and PD cort levels ----

# Comparison of F concentrations in Sisu vs pet dogs by timepoint
tmp_dat <- d_saliva
tmp_dat$tmp_id <- if_else(tmp_dat$condition == 'PD', tmp_dat$id, 'sisu')

tmp_dat <- mutate(tmp_dat, f_norm = as.vector(scale(log(corr_conc))))

# Fill in Sisu information as needed
tmp_dat$age <- if_else(tmp_dat$condition == "UD", sisu_info$age, tmp_dat$age)
tmp_dat$sex <- if_else(tmp_dat$condition == "UD", "female", as.character(tmp_dat$sex))
tmp_dat$sex <- as.factor(tmp_dat$sex)
tmp_dat$weight <- if_else(tmp_dat$condition == "UD", sisu_info$weight, tmp_dat$weight)

tmp_t1 <- filter(tmp_dat, timepoint == 'T1')
tmp_t3 <- filter(tmp_dat, timepoint == 'T3')

tmp_t1 <- mutate(tmp_t1, f_norm = as.vector(scale(f_norm)))
tmp_t3 <- mutate(tmp_t3, f_norm = as.vector(scale(f_norm)))

t1_dog_compare <- brm(f_norm ~ condition + age + sex + weight, data = tmp_t1, prior = b_prior, seed = my_seed, refresh = 0)
summary(t1_dog_compare, prob = 0.90)

t3_dog_compare <- brm(f_norm ~ condition + age + sex + weight, data = tmp_t3, prior = b_prior, seed = my_seed, refresh = 0)
summary(t3_dog_compare, prob = 0.90)

# z score log_f by group for remaining analyses
d_saliva <- d_saliva %>% group_by(condition) %>% mutate(f_norm = as.vector(scale(log(corr_conc))))
hist(d_saliva$f_norm)

## dog salivary F change models ----

#  __> Model for unfamiliar dog - dog cortisol
ud_d_sal_mod <- brm(f_norm~timepoint + (1|id), data = filter(d_saliva, condition == 'UD'), prior = b_prior, sample_prior = TRUE, refresh = 0, seed = my_seed, control = list(adapt_delta = 0.99)) 
pp_check(ud_d_sal_mod, ndraws = 50)
summary(ud_d_sal_mod, prob = 0.90)

get_variables(ud_d_sal_mod)
ud_dog_sal_betas <- gather_draws(ud_d_sal_mod, b_timepointT3) %>% ungroup() %>% select(ud = .value)

# Data for plot.  This is getting posterior distribution of estimated marginal means. 
ud_d_sal_plot_data <- emmeans(ud_d_sal_mod, pairwise ~ timepoint) %>% gather_emmeans_draws() %>% filter(.grid == 'emmeans') %>% mutate(condition = "unfamiliar dog") 

#  __> Model for pet dog - dog cortisol
pd_d_sal_mod <- brm(f_norm~timepoint + sex + age + weight + (1|id), data = filter(d_saliva, condition == 'PD'), prior = b_prior, sample_prior = TRUE, refresh = 0, seed = my_seed, control = list(adapt_delta = 0.99))

pd_dog_sal_betas <- gather_draws(pd_d_sal_mod, b_timepointT3) %>% ungroup() %>% select(pd = .value)

pp_check(pd_d_sal_mod, ndraws = 50)
summary(pd_d_sal_mod, prob = 0.90)

dog_changes <- bind_cols(pd_dog_sal_betas, ud_dog_sal_betas)
theme_set(theme_classic())
mcmc_areas(dog_changes, prob = 0.90) + geom_vline(xintercept = 0, color = "red", linetype = "dashed") + labs(x = "β (Dog Salivary Cortisol SD change from baseline)") + theme(axis.text.y = element_text(size = 14), axis.text.x = element_text(size = 12), axis.title = element_text(size = 16)) + coord_cartesian(ylim = c(1, 2), xlim = c(-1.1, 0.15)) + annotate("text", x = -1.08, y = 1.5, label = "Unfamiliar dog", size = 6, hjust = 0) + annotate("text", x = -1.08, y = 2.5, label = "Pet dog", size = 6, hjust = 0) + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

ggsave("figures/dog_salivary_f_change_across_time.png", width = 6, height = 4)

## behavioral predictors of dog cort change ----
# separate out the data by condition and scale first

# Make AUCi and AUCg based on log concentrations
d_saliva_wide <- mutate(d_saliva_wide, AUCi_log = (log(T3) - log(T1)) / 2, AUCg_log = (log(T1)+log(T3)) / 2)

pd_dsal_wide <- filter(d_saliva_wide, condition == "PD") %>% mutate(across(where(is.numeric), scale))
hist(pd_dsal_wide$AUCi_log)
ud_dsal_wide <- filter(d_saliva_wide, condition == "UD") %>% mutate(across(where(is.numeric), scale))
hist(ud_dsal_wide$AUCi_log)

# Dog salivary AUC behavioral predictors
pd_dsal_AUCi_behav_predictors <- brm(AUCi_log~ sex + age + weight + hai_pc + prop_cooriented + dog_locomoting, data = pd_dsal_wide, prior = b_prior, refresh = 0, seed = my_seed, control = list(adapt_delta = 0.99))
summary(pd_dsal_AUCi_behav_predictors, prob = 0.90)

pd_dsal_AUCg_behav_predictors <- brm(AUCg_log~ sex + age + weight + hai_pc + prop_cooriented + dog_locomoting, data = pd_dsal_wide, prior = b_prior, refresh = 0, seed = my_seed, control = list(adapt_delta = 0.99))
summary(pd_dsal_AUCi_behav_predictors, prob = 0.90)

ud_dsal_AUCi_behav_predictors <- brm(AUCi_log~ hai_pc + prop_cooriented + dog_locomoting, data = ud_dsal_wide, prior = b_prior, refresh = 0, seed = my_seed, control = list(adapt_delta = 0.99))
summary(ud_dsal_AUCi_behav_predictors, prob = 0.90)

ud_dsal_AUCg_behav_predictors <- brm(AUCg_log~ hai_pc + prop_cooriented + dog_locomoting, data = ud_dsal_wide, prior = b_prior, refresh = 0, control = list(adapt_delta = 0.99))
summary(ud_dsal_AUCi_behav_predictors, prob = 0.90, seed = my_seed)

# Make dog salivary AUC behavioral predictor tables - merge tables for AUCi and AUCg - later we will stack all the tables from dog and child.

## dog behavioral predictors of cort tables ----

pd_dog_auci_tab <- tbl_regression(pd_dsal_AUCi_behav_predictors, conf.level = 0.90, label = list(dog_locomoting ~ "locomotion", prop_cooriented ~ "time cooriented", hai_pc ~ "affectionate interaction"), include = c(dog_locomoting, prop_cooriented, hai_pc)) %>% modify_header(label = "**Predictor**") %>% modify_footnote(everything() ~ NA, abbreviation = T) %>% modify_column_alignment(columns = label, align = "right")

pd_dog_aucg_tab <- tbl_regression(pd_dsal_AUCg_behav_predictors, conf.level = 0.90, label = list(dog_locomoting ~ "locomotion", prop_cooriented ~ "time cooriented", hai_pc ~ "affectionate interaction"), include = c(dog_locomoting, prop_cooriented, hai_pc)) %>% modify_header(label = "**Predictor**") %>% modify_footnote(everything() ~ NA, abbreviation = T) %>% modify_column_alignment(columns = label, align = "right")

pd_behav_tbl_merged <- tbl_merge(list(pd_dog_auci_tab, pd_dog_aucg_tab), tab_spanner = c("AUCi","AUCg")) 

ud_dog_auci_tab <- tbl_regression(ud_dsal_AUCi_behav_predictors, conf.level = 0.90, label = list(dog_locomoting ~ "locomotion", prop_cooriented ~ "time cooriented", hai_pc ~ "affectionate interaction"), include = c(dog_locomoting, prop_cooriented, hai_pc)) %>% modify_header(label = "**Predictor**") %>% modify_footnote(everything() ~ NA, abbreviation = T) %>% modify_column_alignment(columns = label, align = "right")

ud_dog_aucg_tab <- tbl_regression(ud_dsal_AUCg_behav_predictors, conf.level = 0.90, label = list(dog_locomoting ~ "locomotion", prop_cooriented ~ "time cooriented", hai_pc ~ "affectionate interaction"), include = c(dog_locomoting, prop_cooriented, hai_pc)) %>% modify_header(label = "**Predictor**") %>% modify_footnote(everything() ~ NA, abbreviation = T) %>% modify_column_alignment(columns = label, align = "right")

ud_behav_tbl_merged <- tbl_merge(list(ud_dog_auci_tab, ud_dog_aucg_tab), tab_spanner = c("AUCi","AUCg")) 


# Child saliva----

# Log child F
h_saliva <- h_saliva %>% ungroup() %>% mutate(log_f = log(corr_conc), log_f_z = as.vector(scale(log_f)))
# 

# Note to self, I'm using a numeric encoding of session number below since if I do it as a factor, level 2 will be completely confounded with the CT condition and that's no good!

long_pd_first <- filter(h_saliva, visit_1 == "PD")
long_ud_first <- filter(h_saliva, visit_1 == "UD")

long_pd_first <- mutate(long_pd_first, session_num = case_when(
  condition == "PD" ~ 1,
  condition == "CT" ~ 2,
  condition == "UD" ~ 3
))

long_ud_first <- mutate(long_ud_first, session_num = case_when(
  condition == "PD" ~ 3,
  condition == "CT" ~ 2,
  condition == "UD" ~ 1
))

h_saliva <- bind_rows(long_pd_first, long_ud_first)

## child time x condition interaction model ----

child_sal_f_mod <- brm(log_f_z~timepoint*condition + sex + age + session_num + (1|id), data = h_saliva, prior = b_prior, sample_prior = TRUE, seed = my_seed, refresh = 0)
pp_check(child_sal_f_mod)
summary(child_sal_f_mod)
# 
means_frame <- emmeans(child_sal_f_mod, pairwise ~ timepoint | condition) 
summary(means_frame, point.est = "mean")

plot(means_frame)

# look at differences between conditions at each timepoint
cond_compare_by_time <- emmeans(child_sal_f_mod, pairwise ~ condition | timepoint)
summary(cond_compare_by_time, point.est = "mean")

## child AUC by condition models ----

# child AUCi by condition
h_saliva_wide <- mutate(h_saliva_wide, AUCi_log = (log(T3) - log(T1)) / 2, AUCg_log = (log(T1)+log(T3)) / 2)

cor(log(h_saliva_wide$T1), log(h_saliva_wide$T3), use = 'pairwise.complete.obs')
cor(log(h_saliva_wide$T1), h_saliva_wide$AUCg_log, use = 'pairwise.complete.obs')


# scale AUCi and AUCg
groups(h_saliva_wide)
h_saliva_wide <- mutate(h_saliva_wide, AUCi_norm = as.vector(scale((AUCi_log))), AUCg_norm = as.vector(scale(AUCg_raw)))
                                 
pd_first <- filter(h_saliva_wide, visit_1 == "PD")
ud_first <- filter(h_saliva_wide, visit_1 == "UD")

# Note to self, I'm using a numeric encoding of session number below since if I do it as a factor, level 2 will be completely confounded with the CT condition and that's no good!

pd_first <- mutate(pd_first, session_num = case_when(
  condition == "PD" ~ 1,
  condition == "CT" ~ 2,
  condition == "UD" ~ 3
))

ud_first <- mutate(ud_first, session_num = case_when(
  condition == "PD" ~ 3,
  condition == "CT" ~ 2,
  condition == "UD" ~ 1
))

h_saliva_wide <- bind_rows(pd_first, ud_first) %>% arrange(id, condition)
h_saliva_wide <- mutate(h_saliva_wide, session_num = as.vector(scale(session_num)))          

# AUCi mod
child_sal_AUCi_mod <- brm(AUCi_norm ~ condition + session_num + sex + age + (1|id), data = h_saliva_wide, prior = b_prior, sample_prior = TRUE, refresh = 0, seed = my_seed)
summary(child_sal_AUCi_mod, prob = 0.90)
pp_check(child_sal_AUCi_mod, ndraws = 50)

my_means <- emmeans(child_sal_AUCi_mod, specs = pairwise ~ condition)
summary(my_means, point.est = "mean")

my_means <- emmeans(child_sal_AUCi_mod, specs = trt.vs.ctrl ~ condition, ref = 2:3)
summary(my_means, point.est = "mean")

# Assemble three plots - I am intentionally making the two individual betas their own panels. 
theme_set(theme_classic(14))

p1 <- mcmc_areas(child_sal_AUCi_mod, pars = c('b_conditionPD'), prob = 0.90) + geom_vline(xintercept = 0, color = "red", linetype = "dashed") + coord_cartesian(ylim = c(0.95, 1), xlim = c(-1, 0.5)) + labs(x = expression(Delta ~ "Child Salivary Cortisol AUCi")) + annotate("text", x = -.75, y = 1.75, label = str_wrap("pet dog relative to control", width = 18), size = 4) + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

p2 <- mcmc_areas(child_sal_AUCi_mod, pars = c('b_conditionUD'), prob = 0.90) + geom_vline(xintercept = 0, color = "red", linetype = "dashed") + coord_cartesian(ylim = c(0.95, 1), xlim = c(-1, 0.5)) + labs(x = expression(Delta ~ "Child Salivary Cortisol AUCi")) + annotate("text", x = -.75, y = 1.75, label = str_wrap("unfamiliar dog relative to control", width = 20), size = 4) + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

# This includes the draws for the contrast between the dog and control conditions
p3 <- emmeans(child_sal_AUCi_mod, specs = trt.vs.ctrl ~ condition, ref = 2:3, epred = TRUE) %>% gather_emmeans_draws() %>% filter(contrast == "CT - avg(PD,UD)") %>% ungroup() %>% select(.value) %>% mutate(dog_v_no_dog = .value * -1) %>% select(-.value) %>% mcmc_areas(prob = 0.90) + geom_vline(xintercept = 0, color = "red", linetype = "dashed") + coord_cartesian(ylim = c(0.95, 1), xlim = c(-1, 0.5)) + labs(x = expression(Delta ~ "Child Salivary Cortisol AUCi")) + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) + annotate("text", x = -0.65, y = 1.75, label = str_wrap("treatment (dog conditions) relative to control", width = 26), size = 4)

p1 + p2 + p3 + plot_layout(ncol = 1) + plot_annotation(tag_levels = 'A')
ggsave("figures/child_salivary_f_AUCi_by_condition.png", width = 5, height = 7)


# AUCg mod
child_sal_AUCg_mod <- brm(AUCg_norm ~ condition + session_num + sex + age + (1|id), data = h_saliva_wide, prior = b_prior, sample_prior = TRUE, control = list(adapt_delta = 0.95), seed = my_seed, refresh = 0)
summary(child_sal_AUCg_mod, prob = 0.90)

my_means <- emmeans(child_sal_AUCg_mod, pairwise~condition)
summary(my_means, point.est = "mean")

my_means <- emmeans(child_sal_AUCg_mod, specs = trt.vs.ctrl ~ condition, ref = 2:3)
summary(my_means, point.est = "mean")

## child behavioral predictors of AUC models ----

# Effect of locomotion controlling for condition
three_cond_child_behav_mod <- brm(AUCi_norm ~ condition + session_num + sex + age + child_locomoting + (1|id), data = h_saliva_wide, prior = b_prior, refresh = 0, seed = my_seed)
summary(three_cond_child_behav_mod, prob = 0.90)

# HAI specific variables 
hsal_pd_cond <- filter(h_saliva_wide, condition == "PD") %>% mutate(across(where(is.numeric), scale))
hsal_ud_cond <- filter(h_saliva_wide, condition == "UD") %>% mutate(across(where(is.numeric), scale))

hsal_pd_AUCi_behav_mod <- brm(AUCi_norm ~ sex + age + session_num + hai_pc + prop_cooriented + child_locomoting, data = hsal_pd_cond, prior = b_prior, refresh = 0, seed = my_seed)
summary(hsal_pd_AUCi_behav_mod, prob = 0.90)

hsal_ud_AUCi_behav_mod <- brm(AUCi_norm ~ sex + age + session_num + hai_pc + prop_cooriented + child_locomoting, data = hsal_ud_cond, prior = b_prior, refresh = 0, seed = my_seed)
summary(hsal_ud_AUCi_behav_mod, prob = 0.90)

hsal_pd_AUCg_behav_mod <- brm(AUCg_norm ~ sex + age + session_num + hai_pc + prop_cooriented + child_locomoting, data = hsal_pd_cond, prior = b_prior, refresh = 0, seed = my_seed)
summary(hsal_pd_AUCg_behav_mod, prob = 0.90)

hsal_ud_AUCg_behav_mod <- brm(AUCg_norm ~ sex + age + session_num + hai_pc + prop_cooriented + child_locomoting, data = hsal_ud_cond, prior = b_prior, refresh = 0, seed = my_seed)
summary(hsal_ud_AUCg_behav_mod, prob = 0.90)

## child behavior predicting AUC tables ----

# Make dog salivary AUC behavioral predictor tables - merge tables for AUCi and AUCg - later we will stack all the tables from dog and child.

pd_child_auci_tab <- tbl_regression(hsal_pd_AUCi_behav_mod, conf.level = 0.90, label = list(child_locomoting ~ "locomotion", prop_cooriented ~ "time cooriented", hai_pc ~ "affectionate interaction"), include = c(child_locomoting, prop_cooriented, hai_pc)) %>% modify_header(label = "**Predictor**") %>% modify_footnote(everything() ~ NA, abbreviation = T) %>% modify_column_alignment(columns = label, align = "right")

pd_child_aucg_tab <- tbl_regression(hsal_pd_AUCg_behav_mod, conf.level = 0.90, label = list(child_locomoting ~ "locomotion", prop_cooriented ~ "time cooriented", hai_pc ~ "affectionate interaction"), include = c(child_locomoting, prop_cooriented, hai_pc)) %>% modify_header(label = "**Predictor**") %>% modify_footnote(everything() ~ NA, abbreviation = T) %>% modify_column_alignment(columns = label, align = "right")

pd_child_behav_tbl_merged <- tbl_merge(list(pd_child_auci_tab, pd_child_aucg_tab), tab_spanner = c("AUCi","AUCg")) 

ud_child_auci_tab <- tbl_regression(hsal_ud_AUCi_behav_mod, conf.level = 0.90, label = list(child_locomoting ~ "locomotion", prop_cooriented ~ "time cooriented", hai_pc ~ "affectionate interaction"), include = c(child_locomoting, prop_cooriented, hai_pc)) %>% modify_header(label = "**Predictor**") %>% modify_footnote(everything() ~ NA, abbreviation = T) %>% modify_column_alignment(columns = label, align = "right")

ud_child_aucg_tab <- tbl_regression(hsal_ud_AUCg_behav_mod, conf.level = 0.90, label = list(child_locomoting ~ "locomotion", prop_cooriented ~ "time cooriented", hai_pc ~ "affectionate interaction"), include = c(child_locomoting, prop_cooriented, hai_pc)) %>% modify_header(label = "**Predictor**") %>% modify_footnote(everything() ~ NA, abbreviation = T) %>% modify_column_alignment(columns = label, align = "right")

ud_child_behav_tbl_merged <- tbl_merge(list(ud_dog_auci_tab, ud_dog_aucg_tab), tab_spanner = c("AUCi","AUCg")) 

tbl_stack(list(pd_child_behav_tbl_merged, ud_child_behav_tbl_merged, pd_behav_tbl_merged, ud_behav_tbl_merged), group_header = c("child cortisol: pet dog condition", "child cortisol: unfamiliar dog condition", "dog cortisol: pet dog condition", "dog cortisol: unfamiliar dog condition")) %>% as_gt() %>% tab_options(row_group.font.weight = "bold") %>% gtsave("tables/dog_and_child_salivary_f_behavioral_predictors.html")


## child survey predicting cort change models ---- 

surveys <- read_rds('data/surveys.rds')
#h_saliva_wide <- left_join(h_saliva_wide, surveys)
hsal_pd_cond <- left_join(hsal_pd_cond, surveys)
hsal_ud_cond <- left_join(hsal_ud_cond, surveys)

# In models below with only two conditions I convert session_num to factor for simplicity

# PD AUCi with survey predictor
pd_child_sal_AUCi_survey <- brm(AUCi_norm ~ sex + age + hab_composite + mp_composite + djls + factor(session_num), data = hsal_pd_cond, prior = b_prior, sample_prior = TRUE, seed = my_seed, refresh = 0)
pp_check(pd_child_sal_AUCi_survey, ndraws = 50)
summary(pd_child_sal_AUCi_survey, prob = 0.90)

get_variables(pd_child_sal_AUCi_survey)
pd_hab_child_sal_betas <- gather_draws(pd_child_sal_AUCi_survey, b_hab_composite) %>% ungroup() %>% select(pd = .value)

# UD F AUCi with survey predictor
ud_child_sal_AUCi_survey <- brm(AUCi_norm ~ sex + age + hab_composite + mp_composite + djls + factor(session_num), data = hsal_ud_cond, prior = b_prior, sample_prior = TRUE, seed = my_seed, refresh = 0)
summary(ud_child_sal_AUCi_survey, prob = 0.90)

ud_hab_child_sal_betas <- gather_draws(ud_child_sal_AUCi_survey, b_hab_composite) %>% ungroup() %>% select(ud = .value)

# Save figure
child_sal_hab_effects <- bind_cols(pd_hab_child_sal_betas, ud_hab_child_sal_betas)
mcmc_areas(child_sal_hab_effects, prob = 0.90) + geom_vline(xintercept = 0, color = "red", linetype = "dashed") + labs(x = str_wrap("β (human-animal bond predicting cortisol response [AUCi])", width = 35)) + annotate("text", x = -0.65, y = 1.7, label = "interaction with \n unfamiliar dog") + annotate("text", x = -0.65, y = 2.7, label = "interaction with \n pet dog") + theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())
ggsave('figures/hab_predicting_child_cort_AUCi.png', width = 5, height = 4)

#PD AUCg with survey predictor
pd_child_sal_AUCg_survey <- brm(AUCg_norm ~ sex + age + hab_composite + mp_composite + djls + factor(session_num), data = hsal_pd_cond, prior = b_prior, sample_prior = TRUE, seed = my_seed, refresh = 0)
summary(pd_child_sal_AUCg_survey, prob = 0.90)

#UD AUCg with survey predictor
ud_child_sal_AUCg_survey <- brm(AUCg_norm ~ sex + age + hab_composite + mp_composite + djls + factor(session_num), data = hsal_ud_cond, prior = b_prior, sample_prior = TRUE, seed = my_seed, refresh = 0)
summary(ud_child_sal_AUCg_survey, prob = 0.90)

## child survey predictor of AUC tables ----

# Make hsal survey predictor tables
pd_child_auci_surv_tbl <- tbl_regression(pd_child_sal_AUCi_survey, conf.level = 0.90, include = c(hab_composite, mp_composite, djls), label = list(hab_composite ~ "human-animal bond", mp_composite ~ "meaning and purpose", djls ~ "loneliness")) %>% modify_header(label = "**Predictor**") %>% modify_footnote(everything() ~ NA, abbreviation = T) %>% modify_column_alignment(columns = label, align = "right")

# Make hsal survey predictor tables
pd_child_aucg_surv_tbl <- tbl_regression(pd_child_sal_AUCg_survey, conf.level = 0.90, include = c(hab_composite, mp_composite, djls), label = list(hab_composite ~ "human-animal bond", mp_composite ~ "meaning and purpose", djls ~ "loneliness")) %>% modify_header(label = "**Predictor**") %>% modify_footnote(everything() ~ NA, abbreviation = T) %>% modify_column_alignment(columns = label, align = "right")

pd_survey_joined <- tbl_merge(list(pd_child_auci_surv_tbl, pd_child_aucg_surv_tbl), tab_spanner = c("AUCi","AUCg"))

# Make hsal survey predictor tables
ud_child_auci_surv_tbl <- tbl_regression(ud_child_sal_AUCi_survey, conf.level = 0.90, include = c(hab_composite, mp_composite, djls), label = list(hab_composite ~ "human-animal bond", mp_composite ~ "meaning and purpose", djls ~ "loneliness")) %>% modify_header(label = "**Predictor**") %>% modify_footnote(everything() ~ NA, abbreviation = T) %>% modify_column_alignment(columns = label, align = "right")

ud_child_aucg_surv_tbl <- tbl_regression(ud_child_sal_AUCg_survey, conf.level = 0.90, include = c(hab_composite, mp_composite, djls), label = list(hab_composite ~ "human-animal bond", mp_composite ~ "meaning and purpose", djls ~ "loneliness")) %>% modify_header(label = "**Predictor**") %>% modify_footnote(everything() ~ NA, abbreviation = T) %>% modify_column_alignment(columns = label, align = "right")

ud_survey_joined <- tbl_merge(list(ud_child_auci_surv_tbl, ud_child_aucg_surv_tbl), tab_spanner = c("AUCi","AUCg"))

tbl_stack(list(pd_survey_joined, ud_survey_joined), group_header = c("child cortisol: pet dog condition", "child cortisol: unfamiliar dog condition")) %>% as_gt() %>% tab_options(row_group.font.weight = "bold") %>% gtsave("tables/child_salivary_f_survey_predictor.html")

# Save results ----
save.image('F_analysis_results.RData')

# if clearing is needed
rm(list=ls())
graphics.off()


# load necessary packages
library(brms)
library(tidyverse)
library(fs)
library(ggplot2)
library(dplyr)
library(broom.mixed)


# set path to dataframe folder
folder_path <- "Updated_Dataframes_CLEAN"
csv_files <- dir(folder_path, full.names = TRUE, pattern = "\\.csv$")

# create reading function
read_with_metadata <- function(file_path) {
  df <- read_csv(file_path, show_col_types = FALSE)
  
  # extract metadata from file names
  file_name <- basename(file_path)
  
  # elements of file names
  lineage <- str_extract(file_name, "Lineage_[AB]") %>% str_remove("Lineage_")
  window <- str_extract(file_name, "slidingWindow_[15]") %>% str_remove("slidingWindow_")
  tree_id <- str_extract(file_name, "tree_\\d{4}") %>% str_remove("tree_") %>% as.integer()
  
  df <- df %>%
    mutate(
      lineage = lineage,
      sliding_window = as.integer(window),
      replicate_id = tree_id
    )
  
  return(df)
}

# read data into single data frame
all_data <- map_dfr(csv_files, read_with_metadata)

# filter to a and b data frames
data_A <- all_data %>% filter(lineage == "A")
data_B <- all_data %>% filter(lineage == "B")

wdc_dist <- bind_rows(
  data_A %>% mutate(lineage = "A"),
  data_B %>% mutate(lineage = "B")
)

# check distribution
ggplot(wdc_dist, aes(x = weighted_diffusion_coefficient, fill = lineage)) +
  geom_density(alpha = 0.4) +
  scale_x_log10() + 
  labs(
    title = "Distribution of Weighted Diffusion Coefficient (WDC)",
    x = "Weighted Diffusion Coefficient (log scale)",
    y = "Density"
  ) +
  scale_fill_manual(values = c("A" = "#9d3d3d", "B" = "#00808a")) +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold"))



# LINEAGE A
# fix column issue
names(data_A)[names(data_A) == "hunter success"] <- "hunter_success"

# reduced formula check
gam_formula <- bf(
  weighted_diffusion_coefficient ~ 
    avg_temp +
    precipitation +
    Uni_Aedes +
    Multi_Aedes +
    Multi_Anopheles +
    Uni_Coquillettidia +
    Multi_Other +
    A_Uni_Aedes +
    A_Multi_Aedes +
    A_Multi_Anopheles +
    A_Uni_Coquillettidia +
    A_Multi_Other +
    hunter_success +
    s(number_of_branches, k = 5, bs = "tp") +
    (1 | sliding_window)
)

# set priors
priors <- c(
  set_prior("normal(0, 0.25)", class = "b"),      # stronger shrinkage for slopes (predictors scaled)
  set_prior("normal(0, 1.5)", class = "Intercept"), # intercept on link scale
  set_prior("exponential(1)", class = "shape"),  # Gamma shape > 0
  set_prior("student_t(3, 0, 1)", class = "sd")  # group-level SDs
)

# run quicker convergence check model
fit_gamma_A <- brm(
  formula = gam_formula,
  data = data_A,
  family = Gamma(link = "log"),
  prior = priors,
  chains = 2, cores = 2, iter = 1200, warmup = 600,
  control = list(adapt_delta = 0.98, max_treedepth = 15),
  backend = "cmdstanr", 
  seed = 1234
)

# summary(fit_gamma_A)
# Smoothing Spline Hyperparameters:
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sds(snumber_of_branches_1)     2.11      1.02     0.94     4.68 1.01      211      305

# Multilevel Hyperparameters:
#   ~sliding_window (Number of levels: 2) 
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)     0.55      0.68     0.03     2.44 1.02      157      282

# Regression Coefficients:
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept                29.06      0.72    27.04    29.98 1.01      175      352
# avg_temp                 -0.34      0.00    -0.35    -0.33 1.00     1199      796
# precipitation            -0.81      0.01    -0.84    -0.79 1.00     1169     1025
# Uni_Aedes                 0.09      0.00     0.09     0.10 1.00     1235     1136
# Multi_Aedes              -0.15      0.00    -0.16    -0.15 1.00     1208     1085
# Multi_Anopheles          -0.11      0.00    -0.11    -0.10 1.00     1244     1127
# Uni_Coquillettidia       -0.19      0.00    -0.20    -0.18 1.00     1050     1086
# Multi_Other               0.03      0.00     0.02     0.03 1.00     1382     1082
# A_Uni_Aedes               0.05      0.00     0.05     0.05 1.00      944      705
# A_Multi_Aedes             0.07      0.00     0.07     0.08 1.00     1148     1128
# A_Multi_Anopheles        -0.41      0.01    -0.42    -0.40 1.01     1141     1045
# A_Uni_Coquillettidia     -0.34      0.01    -0.35    -0.33 1.00     1092      916
# A_Multi_Other            -0.36      0.02    -0.39    -0.33 1.00     1054     1081
# hunter_success           -0.13      0.00    -0.14    -0.13 1.00     1198     1025
# snumber_of_branches_1     4.07      0.20     3.70     4.44 1.00     1137      951

# Further Distributional Parameters:
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# shape     1.14      0.01     1.13     1.16 1.00     1115      583

# Draws were sampled using sample(hmc). For each parameter, Bulk_ESS
# and Tail_ESS are effective sample size measures, and Rhat is the potential
# scale reduction factor on split chains (at convergence, Rhat = 1).


# saveRDS(fit_gamma_A, file = "fit_gamma_A.RDS")

# if you need to re-load
# fit_reduced_gamma_A <- readRDS("fit_reduced_gamma_A.RDS")

# now we're going to try reducing this
# reduced formula check
gam_formula_reduced <- bf(
  weighted_diffusion_coefficient ~ 
    avg_temp +
    precipitation +
    A_Uni_Aedes +
    A_Multi_Aedes +
    A_Multi_Anopheles +
    A_Uni_Coquillettidia +
    A_Multi_Other +
    hunter_success +
    s(number_of_branches, k = 5, bs = "tp") +   # moderate complexity
    (1 | sliding_window)
)

# set priors
priors <- c(
  set_prior("normal(0, 0.25)", class = "b"),      # stronger shrinkage for slopes (predictors scaled)
  set_prior("normal(0, 1.5)", class = "Intercept"), # intercept on link scale
  set_prior("exponential(1)", class = "shape"),  # Gamma shape > 0
  set_prior("student_t(3, 0, 1)", class = "sd")  # group-level SDs
)

# run quicker convergence check model
fit_reduced_gamma_A <- brm(
  formula = gam_formula_reduced,
  data = data_A,
  family = Gamma(link = "log"),
  prior = priors,
  chains = 2, cores = 2, iter = 1200, warmup = 600,
  control = list(adapt_delta = 0.98, max_treedepth = 15),
  backend = "cmdstanr", 
  seed = 1234
)

summary(fit_reduced_gamma_A)
# Smoothing Spline Hyperparameters:
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sds(snumber_of_branches_1)     2.00      0.95     0.93     4.50 1.00      335      623

# Multilevel Hyperparameters:
#   ~sliding_window (Number of levels: 2) 
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)     1.01      0.76     0.18     2.94 1.00      387      456

# Regression Coefficients:
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept                13.03      0.88    10.88    14.25 1.00      363      614
# avg_temp                 -0.11      0.00    -0.11    -0.10 1.00     1095      818
# precipitation            -0.48      0.01    -0.51    -0.46 1.00      999      917
# A_Uni_Aedes               0.07      0.00     0.07     0.07 1.00     1259      978
# A_Multi_Aedes             0.02      0.00     0.02     0.03 1.00      987     1091
# A_Multi_Anopheles        -0.24      0.01    -0.26    -0.23 1.01     1038      729
# A_Uni_Coquillettidia     -0.30      0.01    -0.31    -0.29 1.00      866      677
# A_Multi_Other            -0.45      0.02    -0.49    -0.42 1.00      899      812
# hunter_success           -0.10      0.00    -0.10    -0.10 1.00     1209     1064
# snumber_of_branches_1     5.62      0.19     5.25     5.99 1.01      995      771

# Further Distributional Parameters:
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# shape     0.94      0.01     0.93     0.95 1.00     1378      651

# Draws were sampled using sample(hmc). For each parameter, Bulk_ESS
# and Tail_ESS are effective sample size measures, and Rhat is the potential
# scale reduction factor on split chains (at convergence, Rhat = 1).

# try a longer version of the reduced formula
fit_gamma_A_full <- brm(
  formula = gam_formula_reduced,
  data = data_A,
  family = Gamma(link = "log"),
  prior = priors,
  chains = 4,
  cores = 4,
  iter = 4000,
  warmup = 2000,
  control = list(adapt_delta = 0.995, max_treedepth = 15),
  backend = "cmdstanr",
  seed = 1234
)

# save
# saveRDS(fit_gamma_A_full, "fit_gamma_A_reg.rds")

# combine uni and multi groups into single variables
# write formula
gam_formula_um_A <- bf(
  weighted_diffusion_coefficient ~ 
    avg_temp +
    precipitation +
    I(A_Uni_Aedes + A_Uni_Coquillettidia) +     # combined univoltine group
    I(A_Multi_Aedes + A_Multi_Anopheles + A_Multi_Other) +  # combined multivoltine group
    s(number_of_branches, k = 5, bs = "tp") +
    (1 | sliding_window)
)

# set priors
priors <- c(
  set_prior("normal(0, 0.25)", class = "b"),      # stronger shrinkage for slopes (predictors scaled)
  set_prior("normal(0, 1.5)", class = "Intercept"), # intercept on link scale
  set_prior("exponential(1)", class = "shape"),  # Gamma shape > 0
  set_prior("student_t(3, 0, 1)", class = "sd")  # group-level SDs
)

# fit gamma A
fit_gamma_A_um_full <- brm(
  formula = gam_formula_um_A,
  data = data_A,
  family = Gamma(link = "log"),
  prior = priors,
  chains = 4,
  cores = 4,
  iter = 2400,
  warmup = 1200,
  control = list(adapt_delta = 0.995, max_treedepth = 15),
  backend = "cmdstanr",
  seed = 1234
)
summary(fit_gamma_A_um_full)
# Smoothing Spline Hyperparameters:
#  Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sds(snumber_of_branches_1)     8.66      3.17     4.77    16.77 1.00     1232     1733

# Multilevel Hyperparameters:
#  ~sliding_window (Number of levels: 2) 
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)     1.71      1.02     0.56     4.32 1.00     2155     2228

# Regression Coefficients:
#  Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept                                          8.08      1.16     5.30     9.85 1.00     2110     2674
# avg_temp                                          -0.07      0.00    -0.07    -0.06 1.00     5179     3998
# precipitation                                     -0.73      0.01    -0.75    -0.70 1.00     4241     3034
# IA_Uni_AedesPA_Uni_Coquillettidia                  0.05      0.00     0.05     0.05 1.00     6046     4006
# IA_Multi_AedesPA_Multi_AnophelesPA_Multi_Other     0.01      0.00     0.01     0.01 1.00     6487     4228
# snumber_of_branches_1                              3.42      0.20     3.04     3.83 1.00     4418     3531

# Further Distributional Parameters:
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# shape     0.83      0.00     0.82     0.84 1.00     4857     3245

# Draws were sampled using sample(hmc). For each parameter, Bulk_ESS
# and Tail_ESS are effective sample size measures, and Rhat is the potential
# scale reduction factor on split chains (at convergence, Rhat = 1).


# effect summaries
coef_A <- tidy(fit_gamma_A_full, effects = "fixed", conf.int = TRUE, conf.level = 0.95) %>%
  # remove smoothing and intercept
  filter(!term %in% c("Intercept", "snumber_of_branches_1")) %>%
  # order terms by estimate effect size
  mutate(term = factor(term, levels = term[order(estimate)]))

# format labels
coef_A <- coef_A %>%
  mutate(term_label = case_when(
    term == "avg_temp" ~ "Average temperature",
    term == "precipitation" ~ "Precipitation",
    term == "A_Uni_Aedes" ~ "Positive Uni. Aedes",
    term == "A_Multi_Aedes" ~ "Positive Multi. Aedes",
    term == "A_Multi_Anopheles" ~ "Positive Multi. Anopheles",
    term == "A_Uni_Coquillettidia" ~ "Positive Uni. Coquillettidia",
    term == "A_Multi_Other" ~ "Positive Multi. Other",
    term == "hunter_success" ~ "Deer Activity",
    TRUE ~ term
  ))

# plot effect sizes for a
effects_A <- ggplot(
  subset(coef_A, term_label != "(Intercept)"), 
  aes(x = estimate, y = reorder(term_label, estimate))
) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(size = 3, color = "#9d3d3d") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2, color = "#9d3d3d") +
  labs(
    x = "Posterior estimate (log WDC, Gamma model)",
    y = NULL,
    title = "A: Effects of Ecological Variables \non Weighted Diffusion Coefficient",
    subtitle = "Points = posterior means, bars = 95% credible intervals"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )
effects_A

# combined uni and multi effects
# effect summaries
coef_A_um <- tidy(fit_gamma_A_um, effects = "fixed", conf.int = TRUE, conf.level = 0.95) %>%
  # remove smoothing and intercept
  filter(!term %in% c("Intercept")) %>%
  # order terms by estimate effect size
  mutate(term = factor(term, levels = term[order(estimate)]))

# format labels
coef_A_um <- coef_A_um %>%
  mutate(term_label = case_when(
    term == "avg_temp" ~ "Average temperature",
    term == "precipitation" ~ "Precipitation",
    term == "IA_Uni_AedesPA_Uni_Coquillettidia" ~ "Positive Univoltine",
    term == "IA_Multi_AedesPA_Multi_AnophelesPA_Multi_Other" ~ "Positive Multivoltine",
    term == "snumber_of_branches_1 " ~ "Number of Branches",
    TRUE ~ term
  ))

# plot effect sizes for a
effects_A_um <- ggplot(
  subset(coef_A_um, term_label != "(Intercept)"), 
  aes(x = estimate, y = reorder(term_label, estimate))
) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(size = 3, color = "#9d3d3d") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2, color = "#9d3d3d") +
  labs(
    x = "Posterior estimate (log WDC, Gamma model)",
    y = NULL,
    title = "A: Effects of Ecological Variables \non Weighted Diffusion Coefficient",
    subtitle = "Points = posterior means, bars = 95% credible intervals"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )
effects_A_um



# LINEAGE B
# fix column issue
names(data_B)[names(data_B) == "hunter success"] <- "hunter_success"

# reduced formula check
gam_formula_reduced_B <- bf(
  weighted_diffusion_coefficient ~ 
    avg_temp +
    precipitation +
    B_Uni_Aedes +
    B_Multi_Aedes +
    B_Multi_Anopheles +
    B_Uni_Coquillettidia +
    B_Multi_Other +
    hunter_success +
    s(number_of_branches, k = 5, bs = "tp") +   # moderate complexity
    (1 | sliding_window)
)

# set priors
priors <- c(
  set_prior("normal(0, 0.25)", class = "b"),      # stronger shrinkage for slopes (predictors scaled)
  set_prior("normal(0, 1.5)", class = "Intercept"), # intercept on link scale
  set_prior("exponential(1)", class = "shape"),  # Gamma shape > 0
  set_prior("student_t(3, 0, 1)", class = "sd")  # group-level SDs
)

# fit gamma B
fit_gamma_B <- brm(
  formula = gam_formula_reduced_B,
  data = data_B,
  family = Gamma(link = "log"),
  prior = priors,
  chains = 2, cores = 2, iter = 1200, warmup = 600,
  control = list(adapt_delta = 0.98, max_treedepth = 15),
  backend = "cmdstanr", 
  seed = 1234
)
summary(fit_gamma_B)
# Smoothing Spline Hyperparameters:
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sds(snumber_of_branches_1)     3.89      1.63     1.86     8.25 1.00      318      371

# Multilevel Hyperparameters:
#   ~sliding_window (Number of levels: 2) 
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)     0.95      0.74     0.21     2.73 1.00      562      713

# Regression Coefficients:
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept                 2.96      0.88     0.99     4.51 1.01      532      613
# avg_temp                 -0.23      0.01    -0.25    -0.20 1.00      907      861
# precipitation             3.94      0.05     3.83     4.04 1.00      736      897
# B_Uni_Aedes               0.04      0.01     0.02     0.05 1.00      808      742
# B_Multi_Aedes             0.38      0.02     0.34     0.42 1.00      662      720
# B_Multi_Anopheles         1.10      0.02     1.05     1.14 1.01      693      913
# B_Uni_Coquillettidia     -2.26      0.05    -2.35    -2.16 1.00      903      862
# B_Multi_Other            -0.00      0.25    -0.49     0.51 1.00     1648      766
# hunter_success           -0.22      0.00    -0.23    -0.22 1.00     1093      973
# snumber_of_branches_1    -0.29      0.23    -0.73     0.18 1.00     1287      997

# Further Distributional Parameters:
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# shape     0.75      0.01     0.74     0.76 1.01     1800      842

# Draws were sampled using sample(hmc). For each parameter, Bulk_ESS
# and Tail_ESS are effective sample size measures, and Rhat is the potential
# scale reduction factor on split chains (at convergence, Rhat = 1).

# saveRDS(fit_gamma_B, file = "fit_gamma_B.RDS")

# reduced formula
fit_gamma_B_full <- brm(
  formula = gam_formula_reduced_B,
  data = data_B,
  family = Gamma(link = "log"),
  prior = priors,
  chains = 4,
  cores = 4,
  iter = 4000,
  warmup = 2000,
  control = list(adapt_delta = 0.995, max_treedepth = 15),
  backend = "cmdstanr",
  seed = 1234
)
summary(fit_gamma_B_full)
# Smoothing Spline Hyperparameters:
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sds(snumber_of_branches_1)     3.92      1.54     2.01     7.73 1.00     2663     4133

# Multilevel Hyperparameters:
#   ~sliding_window (Number of levels: 2) 
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)     0.95      0.76     0.22     2.92 1.00     3498     4188

# Regression Coefficients:
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept                 2.96      0.91     0.83     4.47 1.00     3333     3293
# avg_temp                 -0.23      0.01    -0.25    -0.20 1.00     5051     5105
# precipitation             3.94      0.05     3.83     4.05 1.00     5721     5838
# B_Uni_Aedes               0.04      0.01     0.02     0.05 1.00     5622     5644
# B_Multi_Aedes             0.38      0.02     0.34     0.42 1.00     4436     4923
# B_Multi_Anopheles         1.10      0.02     1.05     1.14 1.00     4538     5130
# B_Uni_Coquillettidia     -2.26      0.05    -2.36    -2.16 1.00     5763     5958
# B_Multi_Other            -0.00      0.25    -0.50     0.49 1.00     9991     5417
# hunter_success           -0.22      0.00    -0.23    -0.22 1.00     6853     5601
# snumber_of_branches_1    -0.28      0.24    -0.75     0.19 1.00     9095     5719

# Further Distributional Parameters:
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# shape     0.75      0.01     0.74     0.76 1.00     9953     6047

# Draws were sampled using sample(hmc). For each parameter, Bulk_ESS
# and Tail_ESS are effective sample size measures, and Rhat is the potential
# scale reduction factor on split chains (at convergence, Rhat = 1).

gam_formula_um_B <- bf(
  weighted_diffusion_coefficient ~ 
    avg_temp +
    precipitation +
    I(B_Uni_Aedes + B_Uni_Coquillettidia) +     # combined univoltine group
    I(B_Multi_Aedes + B_Multi_Anopheles + B_Multi_Other) +  # combined multivoltine group
    s(number_of_branches, k = 5, bs = "tp") +
    (1 | sliding_window)
)

# set priors
priors <- c(
  set_prior("normal(0, 0.25)", class = "b"),      # stronger shrinkage for slopes (predictors scaled)
  set_prior("normal(0, 1.5)", class = "Intercept"), # intercept on link scale
  set_prior("exponential(1)", class = "shape"),  # Gamma shape > 0
  set_prior("student_t(3, 0, 1)", class = "sd")  # group-level SDs
)

# combine uni and multi groups into single variables
# fit gamma B
fit_gamma_B_um_full <- brm(
  formula = gam_formula_um_B,
  data = data_B,
  family = Gamma(link = "log"),
  prior = priors,
  chains = 4,
  cores = 4,
  iter = 2400,
  warmup = 1200,
  control = list(adapt_delta = 0.995, max_treedepth = 15),
  backend = "cmdstanr",
  seed = 1234
)
summary(fit_gamma_B_um_full)
# Smoothing Spline Hyperparameters:
#  Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sds(snumber_of_branches_1)     5.60      2.21     2.94    11.31 1.00     1648     1996

# Multilevel Hyperparameters:
#   ~sliding_window (Number of levels: 2) 
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)     1.56      0.94     0.53     4.00 1.00     2585     3115

# Regression Coefficients:
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept                                         -5.42      1.25    -8.08    -3.14 1.00     2300     2889
# avg_temp                                          -0.05      0.02    -0.08    -0.02 1.00     4186     3356
# precipitation                                      2.43      0.05     2.32     2.53 1.00     3962     3675
# IB_Uni_AedesPB_Uni_Coquillettidia                 -0.20      0.01    -0.23    -0.17 1.00     3175     3272
# IB_Multi_AedesPB_Multi_AnophelesPB_Multi_Other     0.48      0.03     0.42     0.53 1.00     3061     3357
# snumber_of_branches_1                              0.72      0.24     0.23     1.19 1.00     4449     3382

# Further Distributional Parameters:
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# shape     0.51      0.00     0.50     0.52 1.00     5733     3478

# Draws were sampled using sample(hmc). For each parameter, Bulk_ESS
# and Tail_ESS are effective sample size measures, and Rhat is the potential
# scale reduction factor on split chains (at convergence, Rhat = 1).

# saveRDS(fit_gamma_B_um_full, file = "fit_gamma_B_um.RDS")

# visualize coefficient results
# effect summaries
coef_B <- tidy(fit_gamma_B_full, effects = "fixed", conf.int = TRUE, conf.level = 0.95) %>%
  # remove smoothing and intercept
  filter(!term %in% c("Intercept", "snumber_of_branches_1")) %>%
  # order terms by estimate effect size
  mutate(term = factor(term, levels = term[order(estimate)]))

# format labels
coef_B <- coef_B %>%
  mutate(term_label = case_when(
    term == "avg_temp" ~ "Average temperature",
    term == "precipitation" ~ "Precipitation",
    term == "B_Uni_Aedes" ~ "Positive Uni. Aedes",
    term == "B_Multi_Aedes" ~ "Positive Multi. Aedes",
    term == "B_Multi_Anopheles" ~ "Positive Multi. Anopheles",
    term == "B_Uni_Coquillettidia" ~ "Positive Uni. Coquillettidia",
    term == "B_Multi_Other" ~ "Positive Multi. Other",
    term == "hunter_success" ~ "Deer Activity",
    TRUE ~ term
  ))

# plot effect sizes for b
effects_B <- ggplot(
  subset(coef_B, term_label != "(Intercept)"), 
  aes(x = estimate, y = reorder(term_label, estimate))
) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(size = 3, color = "#00808a") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2, color = "#00808a") +
  labs(
    x = "Posterior estimate (log WDC, Gamma model)",
    y = NULL,
    title = "B: Effects of Ecological Variables \non Weighted Diffusion Coefficient",
    subtitle = "Points = posterior means, bars = 95% credible intervals"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )
effects_B

# combined uni and multi effects
# effect summaries
coef_B_um <- tidy(fit_gamma_B_um, effects = "fixed", conf.int = TRUE, conf.level = 0.95) %>%
  # remove smoothing and intercept
  filter(!term %in% c("Intercept")) %>%
  # order terms by estimate effect size
  mutate(term = factor(term, levels = term[order(estimate)]))

# format labels
coef_B_um <- coef_B_um %>%
  mutate(term_label = case_when(
    term == "avg_temp" ~ "Average temperature",
    term == "precipitation" ~ "Precipitation",
    term == "IB_Uni_AedesPB_Uni_Coquillettidia" ~ "Positive Univoltine",
    term == "IB_Multi_AedesPB_Multi_AnophelesPB_Multi_Other" ~ "Positive Multivoltine",
    term == "number_of_branches_1   " ~ "Number of Branches",
    TRUE ~ term
  ))

# plot effect sizes for b
effects_B_um <- ggplot(
  subset(coef_B_um, term_label != "(Intercept)"), 
  aes(x = estimate, y = reorder(term_label, estimate))
) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(size = 3, color = "#00808a") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2, color = "#00808a") +
  labs(
    x = "Posterior estimate (log WDC, Gamma model)",
    y = NULL,
    title = "B: Effects of Ecological Variables \non Weighted Diffusion Coefficient",
    subtitle = "Points = posterior means, bars = 95% credible intervals"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )
effects_B_um

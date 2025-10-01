#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# MOSQUITO DATA
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# clear the current R environment of objects and variables
rm(list=ls())
graphics.off()

# check that the pacman package is already installed (and install if not)
if (!require("pacman")) install.packages("pacman")

# unload all currently loaded R packages
pacman::p_unload()

# load the following R packages
pacman::p_load(readr, readxl, tidyverse, janitor, lubridate, PooledInfRate, ggpubr, dplyr, patchwork, cowplot)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# create a function to clean the column names
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
tk_clean_names = function(x) {
  if (!require("pacman")) install.packages("pacman")
  pacman::p_load(stringr)
    
  if(!is.character(x)) {
    stop("x must be a character. Try using colnames(dataframe)")
  } else {
    x = str_to_lower(x)
    x = str_remove(x, "[[:punct:]]")
    x = str_trim(x)
    x = str_replace(x, " ", "_")
    x = str_replace(x, " ", "_")
    return(x)
  }

}



#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#> -------------------------- CONNECTICUT -------------------------------------
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Data Cleaning
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# set working directory and read in the file
setwd("/Users/elliebourgikos/Desktop")
fn = list.files("data_input/mosquito_data/", pattern = "Cumulative", full.names = T)
m = map_df(fn, read_xlsx) 

# clean column names using the tk_clean_names function
colnames(m) = tk_clean_names(colnames(m))

# remove all rows for which nothing was found in the trap
m_n <- m %>% filter(mosquitoes != 0)


# remove all rows that were not tested
# these are rows with the comment "not tested" in the comments
# great than 50 mosquitoes per pool
m_t <- m_n %>% filter(mosquitoes <= 50)

# transform all species names to lower case
m_t$species <- tolower(m_t$species)

# change all ochlerotatus species names to aedes
m_t$species <- gsub("ochlerotatus", "aedes", m_t$species, ignore.case = TRUE)

# change the genus name to upper case
m_t$species <- str_replace(m_t$species, "^(\\w+)", function(match) { str_to_title(match) })


# reorganize the data
m_m = m_t %>%
  mutate(number_of_traps = 1) %>% #to match the MA data
  mutate(st_grp = "CT") %>%
  dplyr::select(st_grp, # state
                species, # mosquito species
                site, # trapping sites
                town, # town
                county, # county
                trap_type, # type of trap used
                date,# date
                cdc_week, # cdc week
                accession, # if present
                mosquitoes, # number of mosquitoes
                number_of_traps, # set equal to 1
                virus # JC or no JC
  ) 


# check if there are any NA values in the date column
any(is.na(m_m$date))
# FALSE


# check the class of the date variable
class(date)
# function

# change the class to date
m_m$date <- as.character(m_m$date)


# add grouping variables to a new dataset, m1
# filter to mosquitoes (26 species)
m1 <- m_m %>%
  filter(species %in% c("Aedes abserratus", "Aedes aurifer", "Aedes canadensis", "Aedes cinereus", 
                        "Aedes communis", "Aedes excrucians", "Aedes provocans", "Aedes sticticus", 
                        "Aedes stimulans", "Aedes thibaulti", "Aedes cantator", "Aedes sollicitans", 
                        "Aedes taeniorhynchus", "Aedes triseriatus", "Aedes trivittatus", "Aedes vexans", 
                        "Anopheles punctipennis", "Anopheles quadrimaculatus", "Anopheles walkeri", 
                        "Psorophora ferox", "Coquillettidia perturbans", "Culex erraticus", "Culex salinarius", 
                        "Culex restuans", "Culiseta melanura", "Culiseta morsitans")) %>%
  mutate(
    year = year(date), # extract the year from date
    month = month(date, label = TRUE, abbr = TRUE), # extract month from date
    week = as.factor(week(date)),  # extract week from date
    virus = if_else(is.na(virus), "none", virus), # replace missing values with "none"
    pos = as.integer(str_detect(virus, "JC")), # new column, 1 for JC, 0 otherwise
    number_of_pools = 1 # creates new column for number of pools, fills with 1
  )

# create a data frame of unique sites to be used
sites = m1 %>%
   distinct(site, trap_type) %>%
   count(site)
# 372

# calculate the total number of traps
traps = sum(sites$n)

# find the total number of years
years = max(year(m1$date)) - min(year(m1$date))
# 25

# count the total number of distinct CDC weeks
weeks = m1 %>% distinct(cdc_week) %>% count()

# rough calculation of max total observations for abundance
traps*years*weeks
# 383125

# calculate the number of pools per trap per day
pools = m1 %>% get_dupes(date, site, trap_type, species)


# create a df to describe the number of traps per CDC week
traps_wk= m1 %>%
  group_by(year, week) %>%
  summarize(no_traps = sum(number_of_traps)) # summarize the group data
  # creates new data frame named traps_wk with columns "year", "week," 
  # and "no_traps" (total number of traps)


# plot the number of traps per CDC week to check code success
ggplot(traps_wk, aes(week, no_traps)) +
  geom_col() +
  facet_wrap(~year)
## looks good



# calculating abundance
# create new data frame to find total number of mosquitoes per trap per night
m2 = m1 %>%
  group_by(st_grp, date, site, trap_type, species) %>% # specify grouping variables
  summarize(mosq = sum(mosquitoes), # summarize grouped variables
            n_traps = sum(number_of_traps),
            n_pools = sum(number_of_pools)
  ) %>%
  ungroup() %>% # remove the grouping after summarizing
  mutate((n_traps = 1),
         mosq_per_trap = mosq/n_traps, # ratio of mosquitoes to traps
         week = as.factor(week(date)),
         year = year(date), # extract the year from date
         month = month(date, label = TRUE, abbr = TRUE))

# plot number of mosquitoes per trap
ggplot(m2, aes(mosq_per_trap, fill = st_grp)) +
  geom_density(alpha = 0.8)



#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#> ANALYSIS BY WEEK
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# abundance drop site info
m_abund_wk = m2 %>%
  group_by(st_grp,
           year, week, species) %>% # group by specified variables
  summarize(abundance = mean(mosq_per_trap, na.rm = T)) %>% # calculates the mean of mosq_per_trap for each group
  ungroup() # removes grouping

# mean and variance
m_abund_wk %>% 
  group_by(species) %>% # group by species
  summarize(mean(abundance), # calculate mean and standard deviation per species group
            sd(abundance))

# group the species into generation-genus groups
m_abund_gen = m2 %>%
  mutate(gen_genus = case_when(
    species %in% c("Aedes abserratus", "Aedes aurifer", "Aedes canadensis", "Aedes cinereus", 
                   "Aedes communis", "Aedes excrucians", "Aedes provocans", "Aedes sticticus", 
                   "Aedes stimulans", "Aedes taeniorhynchus", "Aedes thibaulti") ~ "Univoltine Aedes",
    species %in% c("Aedes cantator", "Aedes sollicitans", "Aedes triseriatus", "Aedes trivittatus", 
                   "Aedes vexans") ~ "Multivoltine Aedes",
    species %in% c("Anopheles punctipennis", "Anopheles quadrimaculatus", 
                   "Anopheles walkeri") ~ "Multivoltine Anopheles",
    species %in% c("Coquillettidia perturbans") ~ "Univoltine Coquillettidia", 
    species %in% c("Culex erraticus", "Culex restuans", "Culex salinarius") ~ "Multivoltine Culex",
    species %in% c("Culiseta melanura", "Culiseta morsitans") ~ "Multivoltine Culiseta", 
    species %in% c("Psorophora ferox") ~ "Multivoltine Psorophora",
    TRUE ~ "Other"
  ))
  
# abundance drop site info
m_abund_gen = m_abund_gen %>%
  group_by(st_grp,
           year, week, gen_genus) %>% # group by specified variables
  summarize(abundance = mean(mosq_per_trap, na.rm = T)) %>% # calculates the mean of mosq_per_trap for each group
  ungroup() # removes grouping

# mean and variance
m_abund_gen %>% 
  group_by(gen_genus) %>% # group by generation genus
  summarize(mean(abundance), # calculate mean and standard deviation per species group
            sd(abundance))

# create new data frame for positive samples
m_abund_pos = m1 %>%
  filter(pos == 1)

# calculate abundance for the new data frame
m_abund_pos = m_abund_pos %>%
  group_by(st_grp, date, site, trap_type, species) %>% # specify grouping variables
  summarize(mosq = sum(mosquitoes), # summarize grouped variables
            n_traps = sum(number_of_traps),
            n_pools = sum(number_of_pools)
  ) %>%
  ungroup() %>% # remove the grouping after summarizing
  mutate((n_traps = 1),
         mosq_per_trap = mosq/n_traps, # ratio of mosquitoes to traps
         week = as.factor(week(date)),
         year = year(date), # extract the year from date
         month = month(date, label = TRUE, abbr = TRUE))

# check the number of mosquitoes per trap
ggplot(m_abund_pos, aes(mosq_per_trap, fill = st_grp)) +
  geom_density(alpha = 0.8)

# separate the positive data frame into univoltine and multivoltine
m_pos_uni = m_abund_pos %>%
  filter(species %in% c("Aedes abserratus", "Aedes aurifer", "Aedes canadensis", 
                  "Aedes cinereus", "Aedes communis", "Aedes excrucians", 
                  "Aedes provocans", "Aedes sticticus", "Aedes stimulans", 
                  "Aedes taeniorhynchus", "Aedes thibaulti", "Coquillettidia perturbans"))

m_pos_multi = m_abund_pos %>%
  filter(species %in% c("Aedes cantator", "Aedes sollicitans", "Aedes triseriatus", "Aedes trivittatus", 
                        "Aedes vexans", "Anopheles punctipennis", "Anopheles quadrimaculatus", 
                        "Anopheles walkeri", "Culex erraticus", "Culex restuans", "Culex salinarius", 
                        "Culiseta melanura", "Culiseta morsitans", "Psorophora ferox"))

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#> DESCRIPTIVE PLOTS
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# create a boxplot of mosquitoes per week from 1997-2023
# png("mosq_per_trap_species.png")
ggplot(m_abund_wk, aes(week,abundance, fill = species, color = species)) +
  geom_boxplot(alpha = 0.5) +
  ggtitle("Average mosquitoes per trap per week 1997-2022") +
  theme_classic() +
  facet_wrap(~st_grp, ncol =1 )
# dev.off()

# stack barplot of total number of mosquitoes sampled per year, per species
ggplot(m_abund_wk, aes(x=year, y=abundance, fill= species)) +
  geom_bar(stat = "identity") +
  ggtitle("Mosquitoes Sampled, 1997-2022") +
  theme_classic() + 
  facet_wrap(~st_grp, ncol =1 )



#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#CALCULATE POOLED INFECTIVITY RATE BY WEEK
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# split to map pIR over month pools
m1_list = m1 %>%
  group_by(st_grp, year, week, species) %>% # group by specified variables
  group_split() # split the grouped data frame into a list

# iterate each element of m1 over the pIR function
mle = purrr::map(m1_list, ~ pIR(pos ~ mosquitoes, data = ., pt.method = "mle"))
# probability of infection as a function of mosquitoes (using maximum likelihood estimation)

# extract values from pIR 
pir = sapply(mle,"[[",1)
lci = sapply(mle,"[[",2)
uci = sapply(mle,"[[",3)

# create a new data frame incorporating the calculated values
m_abund_wk2 = m_abund_wk %>%
  mutate(pir = round(pir,4),
         pir_lci = round(lci,4),
         pir_uci = round(uci,4),
         vector_index = round(abundance * pir,4) # find vector index
  )

# add a column to describe the voltine group
m_abund_wk2 = m_abund_wk2 %>%
  mutate(gen_genus = case_when(
    species %in% c("Aedes abserratus", "Aedes aurifer", "Aedes canadensis", "Aedes cinereus", 
                   "Aedes communis", "Aedes excrucians", "Aedes provocans", "Aedes sticticus", 
                   "Aedes stimulans", "Aedes taeniorhynchus", "Aedes thibaulti") ~ "Univoltine Aedes",
    species %in% c("Aedes cantator", "Aedes sollicitans", "Aedes triseriatus", "Aedes trivittatus", 
                   "Aedes vexans") ~ "Multivoltine Aedes",
    species %in% c("Anopheles punctipennis", "Anopheles quadrimaculatus", 
                   "Anopheles walkeri") ~ "Multivoltine Anopheles",
    species %in% c("Coquillettidia perturbans") ~ "Univoltine Coquillettidia", 
    species %in% c("Culex erraticus", "Culex restuans", "Culex salinarius") ~ "Multivoltine Culex",
    species %in% c("Culiseta melanura", "Culiseta morsitans") ~ "Multivoltine Culiseta", 
    species %in% c("Psorophora ferox") ~ "Multivoltine Psorophora",
    TRUE ~ "Other"
  ))



#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#> PLOT PIR 
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# plot all years by week
ggplot(m_abund_wk2, aes(x = week, y = pir, fill = species, color = species)) +
  geom_boxplot(alpha = 0.5) +
  theme_classic() +
  facet_wrap(~st_grp, nrow = 1) +
  theme(axis.text.x = element_text(angle = 90))


# focusing on specific species
# aedes canadensis
m_abund_wk2_cana = m_abund_wk2 %>%
  filter(species == "Aedes canadensis") %>%
  mutate(year = as.factor(year)) # convert year to a factor variable

# create a line plot by year for aedes canadensis
p_cana_jcv_wk = ggarrange(
  ggplot(m_abund_wk2_cana, aes(x = week)) +
    geom_line(aes(y = abundance, group = year, color = year)) + 
    theme_classic() +
    facet_wrap(~st_grp, ncol = 1) +
    theme(axis.text.x = element_text(angle = 90)),
  
  ggplot(m_abund_wk2_cana, aes(x = week, y = pir, group = year, color = year)) +
    geom_line() +
    theme_classic() +
    facet_wrap(~st_grp, ncol = 1) +
    theme(axis.text.x = element_text(angle = 90)),
  
  ggplot(m_abund_wk2_cana, aes(x = week, y = vector_index, group = year,color = year)) +
    geom_line() +
    theme_classic() +
    facet_wrap(~st_grp, ncol = 1) +
    theme(axis.text.x = element_text(angle = 90)),
  nrow = 1
  
)
p_cana_jcv_wk

# save all of the information from canadensis analysis
# ggsave("data_output/p_cana_jcv_wk.png")
# write.csv(m_abund_wk2_cana, "data_mid/cana_abundance_by_week_CT.csv")
# write.csv(m_abund_wk2, "data_mid/mosquito_abundance_by_week_CT_all_mosq.csv")


# aedes cantator
m_abund_wk2_cant = m_abund_wk2 %>%
  filter(species == "Aedes cantator") %>%
  mutate(year = as.factor(year))

# create a line plot by year for aedes cantator
p_cant_jcv_wk = ggarrange(
  ggplot(m_abund_wk2_cant, aes(x = week)) +
    geom_line(aes(y = abundance, group = year, color = year)) + 
    theme_classic() +
    facet_wrap(~st_grp, ncol = 1) +
    theme(axis.text.x = element_text(angle = 90)),
  
  ggplot(m_abund_wk2_cant, aes(x = week, y = pir, group = year, color = year)) +
    geom_line() +
    theme_classic() +
    facet_wrap(~st_grp, ncol = 1) +
    theme(axis.text.x = element_text(angle = 90)),
  
  ggplot(m_abund_wk2_cant, aes(x = week, y = vector_index, group = year,color = year)) +
    geom_line() +
    theme_classic() +
    facet_wrap(~st_grp, ncol = 1) +
    theme(axis.text.x = element_text(angle = 90)),
  nrow = 1
  
)
p_cant_jcv_wk

# save all of the information from cantator analysis
# ggsave("data_output/p_cant_jcv_wk.png")
# write.csv(m_abund_wk2_cant, "data_mid/cant_abundance_by_week_CT.csv")


# aedes aurifer
m_abund_wk2_aur = m_abund_wk2 %>%
  filter(species == "Aedes aurifer") %>%
  mutate(year = as.factor(year))

# create a line plot by year for aedes aurifer
p_aur_jcv_wk = ggarrange(
  ggplot(m_abund_wk2_aur, aes(x = week)) +
    geom_line(aes(y = abundance, group = year, color = year)) + 
    theme_classic() +
    facet_wrap(~st_grp, ncol = 1) +
    theme(axis.text.x = element_text(angle = 90)),
  
  ggplot(m_abund_wk2_aur, aes(x = week, y = pir, group = year, color = year)) +
    geom_line() +
    theme_classic() +
    facet_wrap(~st_grp, ncol = 1) +
    theme(axis.text.x = element_text(angle = 90)),
  
  ggplot(m_abund_wk2_aur, aes(x = week, y = vector_index, group = year,color = year)) +
    geom_line() +
    theme_classic() +
    facet_wrap(~st_grp, ncol = 1) +
    theme(axis.text.x = element_text(angle = 90)),
  nrow = 1
  
)
p_aur_jcv_wk

# save all of the information from aurifer analysis
# ggsave("data_output/p_aur_jcv_wk.png")
# write.csv(m_abund_wk2_aur, "data_mid/aur_abundance_by_week_CT.csv")


# aedes abserratus
m_abund_wk2_abs = m_abund_wk2 %>%
  filter(species == "Aedes abserratus") %>%
  mutate(year = as.factor(year))

# create a line plot by year for aedes abserratus
p_abs_jcv_wk = ggarrange(
  ggplot(m_abund_wk2_abs, aes(x = week)) +
    geom_line(aes(y = abundance, group = year, color = year)) + 
    theme_classic() +
    facet_wrap(~st_grp, ncol = 1) +
    theme(axis.text.x = element_text(angle = 90)),
  
  ggplot(m_abund_wk2_abs, aes(x = week, y = pir, group = year, color = year)) +
    geom_line() +
    theme_classic() +
    facet_wrap(~st_grp, ncol = 1) +
    theme(axis.text.x = element_text(angle = 90)),
  
  ggplot(m_abund_wk2_abs, aes(x = week, y = vector_index, group = year,color = year)) +
    geom_line() +
    theme_classic() +
    facet_wrap(~st_grp, ncol = 1) +
    theme(axis.text.x = element_text(angle = 90)),
  nrow = 1
  
)
p_abs_jcv_wk

# save all of the information from abserratus analysis
# ggsave("data_output/p_abs_jcv_wk.png")
# write.csv(m_abund_wk2_abs, "data_mid/abs_abundance_by_week_CT.csv")


# anopheles punctipennis
m_abund_wk2_punct = m_abund_wk2 %>%
  filter(species == "Anopheles punctipennis") %>%
  mutate(year = as.factor(year))

# create a line plot by year for anopheles punctipennis
p_punct_jcv_wk = ggarrange(
  ggplot(m_abund_wk2_punct, aes(x = week)) +
    geom_line(aes(y = abundance, group = year, color = year)) + 
    theme_classic() +
    facet_wrap(~st_grp, ncol = 1) +
    theme(axis.text.x = element_text(angle = 90)),
  
  ggplot(m_abund_wk2_punct, aes(x = week, y = pir, group = year, color = year)) +
    geom_line() +
    theme_classic() +
    facet_wrap(~st_grp, ncol = 1) +
    theme(axis.text.x = element_text(angle = 90)),
  
  ggplot(m_abund_wk2_punct, aes(x = week, y = vector_index, group = year,color = year)) +
    geom_line() +
    theme_classic() +
    facet_wrap(~st_grp, ncol = 1) +
    theme(axis.text.x = element_text(angle = 90)),
  nrow = 1
  
)
p_punct_jcv_wk

# save all of the information from punctipennis analysis
# ggsave("data_output/p_punct_jcv_wk.png")
# write.csv(m_abund_wk2_punct, "data_mid/punct_abundance_by_week_CT.csv")


# coquillettidia perturbans
m_abund_wk2_coqp = m_abund_wk2 %>%
  filter(species == "Coquillettidia perturbans") %>%
  mutate(year = as.factor(year))

# create a line plot by year for coquillettidia perturbans
p_coqp_jcv_wk = ggarrange(
  ggplot(m_abund_wk2_coqp, aes(x = week)) +
    geom_line(aes(y = abundance, group = year, color = year)) + 
    theme_classic() +
    facet_wrap(~st_grp, ncol = 1) +
    theme(axis.text.x = element_text(angle = 90)),
  
  ggplot(m_abund_wk2_coqp, aes(x = week, y = pir, group = year, color = year)) +
    geom_line() +
    theme_classic() +
    facet_wrap(~st_grp, ncol = 1) +
    theme(axis.text.x = element_text(angle = 90)),
  
  ggplot(m_abund_wk2_coqp, aes(x = week, y = vector_index, group = year,color = year)) +
    geom_line() +
    theme_classic() +
    facet_wrap(~st_grp, ncol = 1) +
    theme(axis.text.x = element_text(angle = 90)),
  nrow = 1
  
)
p_coqp_jcv_wk

# save all of the information from coquillettidia analysis
# ggsave("data_output/p_coqp_jcv_wk.png")
# write.csv(m_abund_wk2_coqp, "data_mid/coqp_abundance_by_week_CT.csv")



#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#> ANALYSIS BY MONTH
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# find the abundance drop site info, grouping by month
m_abund_mo = m2 %>%
  group_by(st_grp,
           year, month, species) %>% #true grouping variable 24968 just accession)
  summarize(abundance = mean(mosq_per_trap)) %>%
  ungroup()

# find mean and variance values
m_abund_mo %>% 
  group_by(species) %>%
  summarize(mean(abundance),
            sd(abundance))



#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#> DESCRIPTIVE PLOTS
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# create a boxplot mosquitoes per trap per week 1997-2022
ggplot(m_abund_mo, aes(month,abundance, fill = species, color = species)) +
  geom_boxplot(alpha = 0.5) +
  ggtitle("Average mosquitoes per trap per Month 1997-2022") +
  facet_wrap(~st_grp, ncol = 1) +
  theme_classic()



#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# CALCULATE POOLED INFECTIVITY RATE BY MONTH
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# group by month
m1_list = m1 %>%
  group_by(st_grp, year, month, species) %>%
  group_split()

# iterate each element of m1 over the pIR function
mle = map(m1_list, ~ pIR(pos ~ mosquitoes, data = ., pt.method = "mle"))

# extract values from pIR 
pir = sapply(mle,"[[",1)
lci = sapply(mle,"[[",2)
uci = sapply(mle,"[[",3)

# create a new data frame incorporating the calculated values
m_abund_mo2 = m_abund_mo %>%
  ungroup() %>%
  mutate(pir = round(pir,4),
         pir_lci = round(lci,4),
         pir_uci = round(uci,4),
         vector_index = round(abundance * pir,4)
  )

# separate by species
# aedes canadensis
m_abund_mo2_cana = m_abund_mo2 %>%
  filter(species == "Aedes canadensis")

# aedes cantator
m_abund_mo2_cant = m_abund_mo2 %>%
  filter(species == "Aedes cantator")

# aedes aurifer
m_abund_mo2_aur = m_abund_mo2 %>%
  filter(species == "Aedes aurifer")

# aedes abserratus
m_abund_mo2_abs = m_abund_mo2 %>%
  filter(species == "Aedes abserratus")

# anopheles punctipennis
m_abund_mo2_punct = m_abund_mo2 %>%
  filter(species == "Anopheles punctipennis")

# coquillettidia perturbans
m_abund_mo2_coqp = m_abund_mo2 %>%
  filter(species == "Coquillettidia perturbans")


# save all of the information from analysis
# write.csv(m_abund_mo2_cana, "data_mid/cana_abundance_by_month_CT.csv")
# write.csv(m_abund_mo2_cant, "data_mid/cant_abundance_by_month_CT.csv")
# write.csv(m_abund_mo2_aur, "data_mid/aur_abundance_by_month_CT.csv")
# write.csv(m_abund_mo2_abs, "data_mid/abs_abundance_by_month_CT.csv")
# write.csv(m_abund_mo2_punct, "data_mid/punct_abundance_by_month_CT.csv")
# write.csv(m_abund_mo2_coqp, "data_mid/coqp_abundance_by_month_CT.csv")

# write.csv(m_abund_mo2, "data_mid/mosquito_abundance_by_month_CT_all_mosq.csv")


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#> DESCRIPTIVE PLOTS (OVERVIEW)
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# add a gen variable to the data plots
m_abund_wk2 = m_abund_wk2 %>%
  mutate(gen = case_when(
    gen_genus %in% c("Univoltine Aedes", "Univoltine Coquillettidia") ~ "Univoltine",
    gen_genus %in% c("Multivoltine Aedes", "Multivoltine Anopheles", "Multivoltine Culex", 
                     "Multivoltine Culiseta", "Multivoltine Culex", 
                     "Multivoltine Psorophora") ~ "Multivoltine",
    TRUE ~ "Other"
  ))


# stack barplot of total number of mosquitoes sampled per year, per generation genus
# order the mosquito species as desired
gen_genus_order <- c(
  "Multivoltine Aedes",
  "Multivoltine Anopheles",
  "Multivoltine Culex",
  "Multivoltine Culiseta",
  "Multivoltine Psorophora",
  "Univoltine Coquillettidia",
  "Univoltine Aedes"
)

# plot the values as a barplot
abun_plot <- ggplot(m_abund_gen, aes(x = year, y = abundance, 
                        fill = factor(gen_genus, levels = gen_genus_order))) +
  geom_bar(stat = "identity") +
  labs(title = "A", x = "Year", 
       y = "Abundance", 
       fill = "Generation Number and Genus") +
  theme_classic() +
  facet_wrap(~st_grp, ncol = 1) +
  scale_fill_manual(values = c(
    "Multivoltine Aedes" = "#000000",
    "Multivoltine Anopheles" = "#41242b",
    "Multivoltine Culex" = "#583b3e",
    "Multivoltine Culiseta" = "#775c5a",
    "Multivoltine Psorophora" = "#947e7c",
    "Univoltine Coquillettidia" = "#65A2A7",
    "Univoltine Aedes" = "#7BC3C5"
  ))
abun_plot
# abundance is mean number of mosquitoes per trapping location


# uni-voltine species seasonality
# uni voltine order
uni_genus_order <- c(
  "Aedes cinereus",
  "Aedes communis",
  "Aedes excrucians",
  "Aedes provocans",
  "Aedes sticticus",
  "Aedes stimulans",
  "Aedes taeniorhynchus",
  "Aedes thibaulti",
  "Coquillettidia perturbans",
  "Aedes abserratus",
  "Aedes aurifer",
  "Aedes canadensis"
)

# create a bar plot of univoltine species abundance per week over the past years
uni_bar_plot <- ggplot(m_abund_wk %>%
                     filter(species %in% c("Aedes abserratus", "Aedes aurifer", "Aedes canadensis", 
                                           "Aedes cinereus", "Aedes communis", "Aedes excrucians", 
                                           "Aedes provocans", "Aedes sticticus", "Aedes stimulans", 
                                           "Aedes taeniorhynchus", "Aedes thibaulti", "Coquillettidia perturbans")),
                   aes(week, abundance, fill = factor(species, levels = uni_genus_order))) +
  geom_bar(stat = "identity") +
  labs(title = "Univoltine Mosquito Seasonality", x = "Week", 
       y = "Abundance", 
       fill = "Species") +
  theme_classic() + 
  scale_y_continuous(limits = c(0, 3500, na.rm = TRUE)) +
  scale_fill_manual(values = c(
    "Aedes cinereus" = "#171e1f",
    "Aedes communis" = "#181f20",
    "Aedes excrucians" = "#112b32",
    "Aedes provocans" = "#002F3D",
    "Aedes sticticus" = "#00394d",
    "Aedes stimulans" = "#004E6E",
    "Aedes taeniorhynchus" = "#006080",
    "Aedes thibaulti" = "#3b7c8f",
    "Coquillettidia perturbans" = "#65A2A7",
    "Aedes abserratus" = "#71b2b7",
    "Aedes aurifer" = "#7cbfc5",
    "Aedes canadensis" = "#87ccd3"
  )) +
  scale_color_manual(values = c(
    "Aedes cinereus" = "#171e1f",
    "Aedes communis" = "#181f20",
    "Aedes excrucians" = "#112b32",
    "Aedes provocans" = "#002F3D",
    "Aedes sticticus" = "#00394d",
    "Aedes stimulans" = "#004E6E",
    "Aedes taeniorhynchus" = "#006080",
    "Aedes thibaulti" = "#3b7c8f",
    "Coquillettidia perturbans" = "#65A2A7",
    "Aedes abserratus" = "#71b2b7",
    "Aedes aurifer" = "#7cbfc5",
    "Aedes canadensis" = "#87ccd3"
  ))


# create a plot summarizing the average pir per week over the past years
uni_line_plot <- m_abund_wk2 %>%
  filter(species %in% c("Aedes abserratus", "Aedes aurifer", "Aedes canadensis", 
                        "Aedes cinereus", "Aedes communis", "Aedes excrucians", 
                        "Aedes provocans", "Aedes sticticus", "Aedes stimulans", 
                        "Aedes taeniorhynchus", "Aedes thibaulti", "Coquillettidia perturbans")) %>%
  group_by(year, week) %>%
  summarize(mean_pir = mean(pir, na.rm = TRUE),
            ci_low = max(0, mean_pir - qt(0.975, length(pir) - 1) * sd(pir, na.rm = TRUE) / sqrt(length(pir))),
            ci_high = mean_pir + qt(0.975, length(pir) - 1) * sd(pir, na.rm = TRUE) / sqrt(length(pir))) %>%
  ggplot(aes(x = week, y = mean_pir, group = year, color = as.factor(year))) +
  geom_line() +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high, fill = as.factor(year)), alpha = 0.3) +
  labs(title = "Univoltine PIR Seasonality",
       x = "Week",
       y = "Mean JCV PIR",
       color = "Year") +
  theme_classic() +
  scale_y_continuous(limits = c(0, 0.12, na.rm = TRUE)) + 
  scale_color_manual(values = c(
    "#91bfdb", "#a6cee3", "#cab2d6", "#6a51a3", "#6a3d9a", "#483897", 
    "#313696", "#1f78b4", "#4575b4", "#008080", "#009369", "#33a02c", 
    "#b2df8a", "#ffff99", "#fee08b", "#fdbf6f", "#fdae61", "#fe993b", 
    "#ff7f00", "#b15928", "#d6604d", "#d73027", "#e31a1c", "#b2182b", 
    "#67001f", "#560b0e")) +
  scale_fill_manual(values = c(
    "#91bfdb", "#a6cee3", "#cab2d6", "#6a51a3", "#6a3d9a", "#483897", 
    "#313696", "#1f78b4", "#4575b4", "#008080", "#009369", "#33a02c", 
    "#b2df8a", "#ffff99", "#fee08b", "#fdbf6f", "#fdae61", "#fe993b", 
    "#ff7f00", "#b15928", "#d6604d", "#d73027", "#e31a1c", "#b2182b", 
    "#67001f", "#560b0e")) +
  guides(color = guide_legend(title = "Year"), fill = guide_legend(title = "Year"))

# combine the univoltine plots
uni_combined_plot <- ggarrange(uni_bar_plot, uni_line_plot, nrow = 2, align = "v")

# print the combined univoltine plot
print(uni_combined_plot)



# multi voltine order
multi_genus_order <- c(
  "Culiseta morsitans",
  "Culiseta melanura",
  "Culex salinarius",
  "Culex restuans",
  "Culex erraticus",
  "Psorophora ferox",
  "Anopheles walkeri",
  "Anopheles quadrimaculatus",
  "Aedes vexans",
  "Aedes trivittatus",
  "Aedes triseriatus",
  "Aedes sollicitans",
  "Anopheles punctipennis",
  "Aedes cantator"
)


# create a bar plot of multivoltine species abundance per week over the past years
multi_bar_plot <- ggplot(m_abund_wk %>%
                         filter(species %in% c("Aedes cantator", "Aedes sollicitans", "Aedes triseriatus", "Aedes trivittatus", 
                                               "Aedes vexans", "Anopheles punctipennis", "Anopheles quadrimaculatus", 
                                               "Anopheles walkeri", "Culex erraticus", "Culex restuans", "Culex salinarius", 
                                               "Culiseta melanura", "Culiseta morsitans", "Psorophora ferox")),
                       aes(week, abundance, fill = factor(species, levels = multi_genus_order))) +
  geom_bar(stat = "identity") +
  labs(title = "Multivoltine Mosquito Seasonality", x = "Week", 
       y = "Abundance", 
       fill = "Species") +
  theme_classic() + 
  scale_y_continuous(limits = c(0, 3500, na.rm = TRUE)) +
  scale_fill_manual(values = c(
    "Culiseta morsitans" = "#000000",
    "Culiseta melanura" = "#32161f",
    "Culex salinarius" = "#41242b",
    "Culex restuans" = "#583b3e",
    "Culex erraticus" = "#684b4c",
    "Psorophora ferox" = "#775c5a",
    "Anopheles walkeri" = "#8d6866",
    "Anopheles quadrimaculatus" = "#a47574",
    "Aedes vexans" = "#bb8180",
    "Aedes trivittatus" = "#de9494",
    "Aedes triseriatus" = "#eea6a5",
    "Aedes sollicitans" = "#fdb5b4",
    "Anopheles punctipennis" = "#fac1c0",
    "Aedes cantator" = "#f6d6d5"
  )) + 
  scale_color_manual(values = c(
    "Culiseta morsitans" = "#000000",
    "Culiseta melanura" = "#32161f",
    "Culex salinarius" = "#41242b",
    "Culex restuans" = "#583b3e",
    "Culex erraticus" = "#684b4c",
    "Psorophora ferox" = "#775c5a",
    "Anopheles walkeri" = "#8d6866",
    "Anopheles quadrimaculatus" = "#a47574",
    "Aedes vexans" = "#bb8180",
    "Aedes trivittatus" = "#de9494",
    "Aedes triseriatus" = "#eea6a5",
    "Aedes sollicitans" = "#fdb5b4",
    "Anopheles punctipennis" = "#fac1c0",
    "Aedes cantator" = "#f6d6d5"
  ))

# create a plot summarizing the average pir per week over the past years
multi_line_plot <- m_abund_wk2 %>%
  filter(species %in% c("Aedes cantator", "Aedes sollicitans", "Aedes triseriatus", "Aedes trivittatus", 
                        "Aedes vexans", "Anopheles punctipennis", "Anopheles quadrimaculatus", 
                        "Anopheles walkeri", "Culex erraticus", "Culex restuans", "Culex salinarius", 
                        "Culiseta melanura", "Culiseta morsitans", "Psorophora ferox")) %>%
  group_by(year, week) %>%
  summarize(mean_pir = mean(pir, na.rm = TRUE),
            ci_low = max(0, mean_pir - qt(0.975, length(pir) - 1) * sd(pir, na.rm = TRUE) / sqrt(length(pir))),
            ci_high = mean_pir + qt(0.975, length(pir) - 1) * sd(pir, na.rm = TRUE) / sqrt(length(pir))) %>%
  ggplot(aes(x = week, y = mean_pir, group = year, color = as.factor(year))) +
  geom_line() +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high, fill = as.factor(year)), alpha = 0.3) +
  labs(title = "Multivoltine PIR Seasonality",
       x = "Week",
       y = "Mean JCV PIR",
       color = "Year") +
  theme_classic() +
  scale_y_continuous(limits = c(0, 0.12, na.rm = TRUE)) + 
  scale_color_manual(values = c(
    "#91bfdb", "#a6cee3", "#cab2d6", "#6a51a3", "#6a3d9a", "#483897", 
    "#313696", "#1f78b4", "#4575b4", "#008080", "#009369", "#33a02c", 
    "#b2df8a", "#ffff99", "#fee08b", "#fdbf6f", "#fdae61", "#fe993b", 
    "#ff7f00", "#b15928", "#d6604d", "#d73027", "#e31a1c", "#b2182b", 
    "#67001f", "#560b0e")) +
  scale_fill_manual(values = c(
    "#91bfdb", "#a6cee3", "#cab2d6", "#6a51a3", "#6a3d9a", "#483897", 
    "#313696", "#1f78b4", "#4575b4", "#008080", "#009369", "#33a02c", 
    "#b2df8a", "#ffff99", "#fee08b", "#fdbf6f", "#fdae61", "#fe993b", 
    "#ff7f00", "#b15928", "#d6604d", "#d73027", "#e31a1c", "#b2182b", 
    "#67001f", "#560b0e")) +
  guides(color = guide_legend(title = "Year"), fill = guide_legend(title = "Year"))

# combine the multivoltine plots
multi_combined_plot <- ggarrange(multi_bar_plot, multi_line_plot, nrow = 2, align = "v")

# print the combined multivoltine plot
print(multi_combined_plot)



#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#> FIGURE 3 VISUALIZATION
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# create a new data frame to include only univoltine species
m_abund_wk2_u = m_abund_wk2 %>%
  filter(gen == "Univoltine") %>%
  mutate(year = as.factor(year)) %>%
  group_by(gen, week, year) %>%
  summarize(total_abundance = sum(abundance, na.rm = TRUE)) %>% # sum total abundance per week per year
  summarize(mean_abundance = mean(total_abundance, na.rm = TRUE), # average abundance 
            ci_low = max(0, mean_abundance - qt(0.975, length(total_abundance) - 1) * sd(total_abundance, na.rm = TRUE) / sqrt(length(total_abundance))),
            ci_high = mean_abundance + qt(0.975, length(total_abundance) - 1) * sd(total_abundance, na.rm = TRUE) / sqrt(length(total_abundance)))

# create a new data frame to include only multivoltine species
m_abund_wk2_m = m_abund_wk2 %>%
  filter(gen == "Multivoltine") %>%
  mutate(year = as.factor(year)) %>%
  group_by(gen, week, year) %>%
  summarize(total_abundance = sum(abundance, na.rm = TRUE)) %>% # sum total abundance per week per year
  summarize(mean_abundance = mean(total_abundance, na.rm = TRUE), # average abundance 
            ci_low = max(0, mean_abundance - qt(0.975, length(total_abundance) - 1) * sd(total_abundance, na.rm = TRUE) / sqrt(length(total_abundance))),
            ci_high = mean_abundance + qt(0.975, length(total_abundance) - 1) * sd(total_abundance, na.rm = TRUE) / sqrt(length(total_abundance))) 

# combine the plots 
abundance_plot_combined <- ggplot() + 
  geom_line(data = m_abund_wk2_u, aes(x = 1:length(week), y = mean_abundance, color = gen)) + 
  geom_line(data = m_abund_wk2_m, aes(x = 1:length(week), y = mean_abundance, color = gen)) +
  geom_ribbon(data = m_abund_wk2_u, aes(x = 1:length(week), ymin = ci_low, ymax = ci_high, fill = gen), alpha = 0.4, colour = NA) + 
  geom_ribbon(data = m_abund_wk2_m, aes(x = 1:length(week), ymin = ci_low, ymax = ci_high, fill = gen), alpha = 0.4, colour = NA) + 
  labs(x = "Week", 
       y = "Mean Abundance", 
       color = "Generation",
       fill = "Generation") +
  theme_classic(base_size = 17, base_family = "Helvetica") +
  theme(axis.title.x = element_text(face = "bold"), 
        axis.title.y = element_text(face = "bold")) +
  scale_color_manual(values = c("#775c5a", "#7BC3C5")) +
  scale_fill_manual(values = c("#775c5a", "#7BC3C5")) +
  scale_x_continuous(name = "Month", 
                     limits = c(1, 25), breaks = c(1:25), 
                     labels = c("", "Jun", "", "", "", "", 
                                "Jul", "", "", "", "Aug", "", 
                                "", "", "Sept", "", "", "", 
                                "", "Oct", "", "", "", "", "Nov")) +
  theme(
    legend.position      = c(1, 1),
    legend.justification = c(1, 1) 
  )
abundance_plot_combined


# now we do the same thing for the pir values
# only univoltine species
m_pir_u = m_abund_wk2 %>%
  filter(gen == "Univoltine") %>%
  group_by(gen, week) %>%
  summarize(mean_pir = mean(pir, na.rm = TRUE),
            ci_low = mean_pir - qt(0.975, n() - 1) * sd(pir, na.rm = TRUE) / sqrt(n()),
            ci_high = mean_pir + qt(0.975, n() - 1) * sd(pir, na.rm = TRUE) / sqrt(n()))
m_pir_u$ci_low <- ifelse(m_pir_u$ci_low < 0, 0, m_pir_u$ci_low)

# only multivoltine species
m_pir_m = m_abund_wk2 %>%
  filter(gen == "Multivoltine") %>%
  group_by(gen, week) %>%
  summarize(mean_pir = mean(pir, na.rm = TRUE),
            ci_low = mean_pir - qt(0.975, n() - 1) * sd(pir, na.rm = TRUE) / sqrt(n()),
            ci_high = mean_pir + qt(0.975, n() - 1) * sd(pir, na.rm = TRUE) / sqrt(n()))
m_pir_m$ci_low <- ifelse(m_pir_m$ci_low < 0, 0, m_pir_m$ci_low)

# combine the plots 
pir_plot_combined <- ggplot() + 
  geom_line(data = m_pir_u, aes(x = 1:length(week), y = mean_pir, color = gen)) + 
  geom_line(data = m_pir_m, aes(x = 1:length(week), y = mean_pir, color = gen)) +
  geom_ribbon(data = m_pir_u, aes(x = 1:length(week), ymin = ci_low, ymax = ci_high, fill = gen), alpha = 0.4, colour = NA) + 
  geom_ribbon(data = m_pir_m, aes(x = 1:length(week), ymin = ci_low, ymax = ci_high, fill = gen), alpha = 0.4, colour = NA) + 
  labs(x = "Week", 
       y = "Mean Infection Rate", 
       color = "Generation",
       fill = "Generation") +
  theme_classic(base_size = 17, base_family = "Helvetica") +
  theme(axis.title.x = element_text(face = "bold"), 
        axis.title.y = element_text(face = "bold")) +
  theme(legend.position = "none") + 
  coord_cartesian(ylim = c(0, 0.004)) +
  scale_color_manual(values = c("#775c5a", "#7BC3C5")) +
  scale_fill_manual(values = c("#775c5a", "#7BC3C5")) + 
  scale_x_continuous(name = "Month", 
                     limits = c(1, 25), breaks = c(1:25), 
                     labels = c("", "Jun", "", "", "", "", 
                                "Jul", "", "", "", "Aug", "", 
                                "", "", "Sept", "", "", "", 
                                "", "Oct", "", "", "", "", "Nov"))
pir_plot_combined



# now we do the same thing for the pir values
# only univoltine species
m_vi_u = m_abund_wk2 %>%
  filter(gen == "Univoltine") %>%
  group_by(gen, week) %>%
  summarize(mean_vector_index = mean(vector_index, na.rm = TRUE),
            ci_low = mean_vector_index - qt(0.975, n() - 1) * sd(vector_index, na.rm = TRUE) / sqrt(n()),
            ci_high = mean_vector_index + qt(0.975, n() - 1) * sd(vector_index, na.rm = TRUE) / sqrt(n()))
m_vi_u$ci_low <- ifelse(m_vi_u$ci_low < 0, 0, m_vi_u$ci_low)

# only multivoltine species
m_vi_m = m_abund_wk2 %>%
  filter(gen == "Multivoltine") %>%
  group_by(gen, week) %>%
  summarize(mean_vector_index = mean(vector_index, na.rm = TRUE),
            ci_low = mean_vector_index - qt(0.975, n() - 1) * sd(vector_index, na.rm = TRUE) / sqrt(n()),
            ci_high = mean_vector_index + qt(0.975, n() - 1) * sd(vector_index, na.rm = TRUE) / sqrt(n()))
m_vi_m$ci_low <- ifelse(m_vi_m$ci_low < 0, 0, m_vi_m$ci_low)


# combine the plots 
vi_plot_combined <- ggplot() + 
  geom_line(data = m_vi_u, aes(x = 1:length(week), y = mean_vector_index, color = gen)) + 
  geom_line(data = m_vi_m, aes(x = 1:length(week), y = mean_vector_index, color = gen)) +
  geom_ribbon(data = m_vi_u, aes(x = 1:length(week), ymin = ci_low, ymax = ci_high, fill = gen), alpha = 0.4, colour = NA) + 
  geom_ribbon(data = m_vi_m, aes(x = 1:length(week), ymin = ci_low, ymax = ci_high, fill = gen), alpha = 0.4, colour = NA) + 
  labs(x = "Month", 
       y = "Mean Vector Index", 
       color = "Generation",
       fill = "Generation") +
  theme_classic(base_size = 17, base_family = "Helvetica") +
  theme(axis.title.x = element_text(face = "bold"), 
        axis.title.y = element_text(face = "bold")) +
  theme(legend.position = "none") +
  coord_cartesian(ylim = c(0, 0.04)) +
  scale_color_manual(values = c("#775c5a", "#7BC3C5")) +
  scale_fill_manual(values = c("#775c5a", "#7BC3C5")) + 
  scale_x_continuous(name = "Month", 
                     limits = c(1, 25), breaks = c(1:25), 
                     labels = c("", "Jun", "", "", "", "", 
                                "Jul", "", "", "", "Aug", "", 
                                "", "", "Sept", "", "", "", 
                                "", "Oct", "", "", "", "", "Nov"))
vi_plot_combined

# stack vertically
fig1_plot <- ggarrange(abundance_plot_combined, pir_plot_combined, vi_plot_combined, nrow = 3, align = "v")
fig1_plot

fig3a_plot <- ggarrange(abundance_plot_combined, pir_plot_combined, vi_plot_combined, ncol = 3, align = "h")
fig3a_plot

# save the plot
ggsave("fig3_surv_h.png", plot = fig3a_plot, width = 20, height = 5, dpi = 600, device = "png")




# remaking positivity abundance plot
# collapse other multivoltine groups into the other category
jcv_pos$voltine_genus_grouped <- recode(jcv_pos$voltine_genus,
                                        "Multi_Culex" = "Multi_Other",
                                        "Multi_Culiseta" = "Multi_Other",
                                        "Multi_Psorophora" = "Multi_Other",
                                        .default = jcv_pos$voltine_genus
)

# define plotting order again
volt_order_grouped <- c(
  "Multi_Aedes", "Multi_Anopheles", "Multi_Other", 
  "Uni_Coquillettidia", "Uni_Aedes"
)

# plot
abun_pos_plot <- ggplot(jcv_pos, aes(x = year, y = abundance,
                                     fill = factor(voltine_genus_grouped, levels = volt_order_grouped))) +
  geom_bar(stat = "identity") +
  labs(x = "Year", y = "Abundance", fill = "Generation-Genus") +
  theme_classic(base_family = "Helvetica", base_size = 15) +
  theme(
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 13),
    legend.position = c(0.02, 0.90),     # position in upper left
    legend.justification = c(0, 1)       # align legend box from top left corner
  ) +
  scale_fill_manual(
    values = c(
      "Multi_Aedes" = "#000000",
      "Multi_Anopheles" = "#593b40",
      "Multi_Other" = "#b9a3a3",
      "Uni_Coquillettidia" = "#4f7f84",
      "Uni_Aedes" = "#86cacc"
    ),
    labels = c(
      "Multivoltine Aedes", 
      "Multivoltine Anopheles", 
      "Multivoltine Other",
      "Univoltine Coquillettidia", 
      "Univoltine Aedes"
    )
  ) +
  xlim(1997, NA)

abun_pos_plot

# save the plot
ggsave("suppl_abund.png", plot = abun_pos_plot, width = 20, height = 6, dpi = 600, device = "png")


# additional, new figure 3a
# making a new panel a for fig3
summary_pos <- jcv_pos %>%
  count(species, voltine_genus_grouped)

# re order the bars to be in descending order based on isolation number
summary_pos <- summary_pos %>%
  group_by(species) %>%
  summarise(n = sum(n), voltine_genus_grouped = first(voltine_genus_grouped)) %>%
  ungroup() %>%
  mutate(species = reorder(species, -n))

# plot bars by species
fig3a <- ggplot(summary_pos, aes(x = species, y = n, fill = voltine_genus_grouped)) +
  geom_bar(stat = "identity") +
  labs(
    x = "Species",
    y = "Cumulative Number of JCV Isolations",
    fill = "Mosquito Type"
  ) +
  theme_classic(base_family = "Helvetica", base_size = 15) +
  theme(
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.text = element_text(size = 13),
    legend.position = c(0.98, 0.98),
    legend.justification = c(1, 1)
  ) +
  scale_x_discrete(labels = c(
    "Aedes_abserratus" = "Ae. abserratus",
    "Aedes_aurifer" = "Ae. aurifer", 
    "Aedes_canadensis" = "Ae. canadensis",
    "Aedes_cantator" = "Ae. cantator", 
    "Aedes_cinereus" = "Ae. cinereus", 
    "Aedes_communis" = "Ae. communis", 
    "Aedes_excrucians" = "Ae. excrucians", 
    "Aedes_provocans" = "Ae. provocans", 
    "Aedes_sollicitans" = "Ae. sollicitans", 
    "Aedes_sticticus" = "Ae. sticticus", 
    "Aedes_stimulans" = "Ae. stimulans", 
    "Aedes_taeniorhynchus" = "Ae. taeniorhynchus", 
    "Aedes_thibaulti" = "Ae. thibaulti",
    "Aedes_trivittatus" = "Ae. trivittatus",
    "Aedes_vexans" = "Ae. vexans",
    "Anopheles_punctipennis" = "An. punctipennis", 
    "Anopheles_quadrimaculatus" = "An. quadrimaculatus", 
    "Coquillettidia_perturbans" = "Cq. perturbans", 
    "Culex_pipiens" = "Cx. pipiens", 
    "Culex_restuans" = "Cx. restuans", 
    "Culex_salinarius" = "Cx. salinarius", 
    "Culiseta_melanura" = "Cs. melanura", 
    "Psorophora_ferox" = "Ps. ferox"
  )) +
  scale_fill_manual(
    values = c(
      "Multi_Aedes" = "#000000",
      "Multi_Anopheles" = "#593b40",
      "Multi_Other" = "#b9a3a3",
      "Uni_Coquillettidia" = "#4f7f84",
      "Uni_Aedes" = "#86cacc"
    ),
    labels = c(
      "Multivoltine Aedes", 
      "Multivoltine Anopheles", 
      "Multivoltine Other",
      "Univoltine Aedes",
      "Univoltine Coquillettidia"
    )) 


# save the plot
ggsave("fig3a.png", plot = fig3a, width = 18, height = 6, dpi = 600, device = "png")


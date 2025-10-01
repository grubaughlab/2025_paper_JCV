## making alluvial plots for mosquito DTA

# if clearing is needed
rm(list=ls())
graphics.off()

# load necessary libraries
library(ggalluvial)
library(ggplot2)
library(dplyr)
library(tidyr)

# read in lineage a data
a_genus <- read.table("a_genus_jumps.txt", header = TRUE, sep = "\t")

# set color palette for mosquito colors (and layering)
mosquito_categories <- c("Uni_Aedes", "Uni_Coquillettidia", "Multi_Aedes", "Multi_Anopheles", "Multi_Other")
cols <- c("#86cacc", "#4f7f84", "#000000", "#593b40", "#b9a3a3")
names(cols) <- mosquito_categories

# count cumulative jumps
a_summary <- a_genus %>%
  group_by(from, to) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(
    from = factor(from, levels = mosquito_categories),
    to = factor(to, levels = mosquito_categories)
  )

# reshape to match axes
a_long <- a_summary %>%
  pivot_longer(cols = c(from, to), names_to = "axis", values_to = "category")


# plot alluvial plot
a_alluvial <- ggplot(a_summary,
            aes(axis1 = from, axis2 = to, y = count)) +
  geom_alluvium(aes(fill = from), width = 0.25, alpha = 0.70) +
  geom_stratum(aes(fill = after_stat(stratum)), width = 0.25, color = NA) +
  scale_fill_manual(values = cols) +
  scale_x_discrete(limits = c("", ""), expand = c(0.05, 0.05)) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_family = "Helvetica") +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    # plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    legend.position = "none", 
    panel.grid = element_blank()
  )
a_alluvial

# save image
ggsave("a_jumps_alluvial.png", plot = a_alluvial, width = 6, height = 6, dpi = 600, bg = "white")




# move on to lineage b, repeating the process
b_genus <- read.table("b_genus_jumps.txt", header = TRUE, sep = "\t")

# count cumulative jumps
b_summary <- b_genus %>%
  group_by(from, to) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(
    from = factor(from, levels = mosquito_categories),
    to = factor(to, levels = mosquito_categories)
  )

# reshape to match axes
b_long <- b_summary %>%
  pivot_longer(cols = c(from, to), names_to = "axis", values_to = "category")

# plot alluvial plot
b_alluvial <- ggplot(b_summary,
                     aes(axis1 = from, axis2 = to, y = count)) +
  geom_alluvium(aes(fill = from), width = 0.25, alpha = 0.70) +
  geom_stratum(aes(fill = after_stat(stratum)), width = 0.25, color = NA) +
  scale_fill_manual(values = cols) +
  scale_x_discrete(limits = c("", ""), expand = c(0.05, 0.05)) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_family = "Helvetica") +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    # plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    legend.position = "none", 
    panel.grid = element_blank()
  )
b_alluvial

# save image
ggsave("b_jumps_alluvial.png", plot = b_alluvial, width = 6, height = 6, dpi = 600, bg = "white")



# calculate proportions for main text
# calculate % of jumps from uni_aedes in lineage a
a_uni_aedes_pct <- a_genus %>%
  summarise(
    from_Uni_Aedes = sum(from == "Uni_Aedes"), 
    total = n(), 
    percent_from_Uni_Aedes = 100 * from_Uni_Aedes / total
  )
print(a_uni_aedes_pct)

# repeat for lineage b
b_uni_aedes_pct <- b_genus %>%
  summarise(
    from_Uni_Aedes = sum(from == "Uni_Aedes"), 
    total = n(), 
    percent_from_Uni_Aedes = 100 * from_Uni_Aedes / total
  )
print(b_uni_aedes_pct)



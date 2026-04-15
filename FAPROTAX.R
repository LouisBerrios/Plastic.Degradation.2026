##### FAPROTAX Analysis #####

library(ggplot2);library(ggpubr);library(dplyr);library(stringr)

fun <- read.csv("func_table_clean.csv")

# add Compartment and Treatment.Fin from Treatment
fun2 <- fun %>%
  mutate(
    Compartment = case_when(
      str_detect(Treatment, "Control") ~ "Control",
      Treatment == "FNoP"              ~ "Bulk Soil",
      Treatment == "NoPNoF"            ~ "Bulk Soil",
      Treatment == "PF"                ~ "Bulk Soil",
      Treatment == "PF_PL"             ~ "Plastisphere",
      Treatment == "PnoF"              ~ "Bulk Soil",
      Treatment == "PnoF_PL"           ~ "Plastisphere",
      TRUE                             ~ NA_character_
    ),
    Treatment.Fin = case_when(
      str_detect(Treatment, "Control") ~ "Control",
      Treatment == "FNoP"              ~ "Fungus / No Plastic",
      Treatment == "NoPNoF"            ~ "No Plastic / No Fungus",
      Treatment == "PF"                ~ "Plastic+Fungus",
      Treatment == "PF_PL"             ~ "Plastic+Fungus",
      Treatment == "PnoF"              ~ "Plastic / No Fungus",
      Treatment == "PnoF_PL"           ~ "Plastic / No Fungus",
      TRUE                             ~ NA_character_
    )
  )


##### PLOTTING #####

### HEATMAP ###

library(tidyverse)

fun_long <- fun2 %>%
  pivot_longer(
    cols = -c(Sample, Treatment.Fin, Treatment, Compartment),  # keep these, melt everything else
    names_to = "Variable",
    values_to = "Value"
  )

fun_long_NC <- subset(fun_long, Compartment !="Control")
# Grouped by Treatment only
ggplot(fun_long_NC, aes(x = Variable, y = Treatment.Fin, fill = log(Value))) +
  geom_tile(color = "white") +
  facet_grid(rows = vars(Compartment), scales = "free_y", space = "free_y") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text.y = element_text(face = "bold")) +
  scale_fill_gradient2(
    low = "#313695", mid = "white", high = "#A50026",
    midpoint = 0, name = "log(Value)",
    limits = c(0, 8)  # optional; adjust to clip extremes
  )

library(tidyverse)

# ---------------------------
# 0) Base data (rename if needed)
# ---------------------------
dat <- fun_long_NC %>%
  mutate(
    Value_log = log1p(Value)  # safe log (handles zeros)
  )

# ---------------------------
# 1) Cluster VARIABLES (global, across all compartments)
# ---------------------------
# Build a matrix with rows = (Compartment,Treatment.Fin) combos, cols = Variables
mat_all <- dat %>%
  select(Compartment, Treatment.Fin, Variable, Value_log) %>%
  pivot_wider(
    id_cols   = c(Compartment, Treatment.Fin),
    names_from = Variable,
    values_from = Value_log,
    values_fn   = mean,
    values_fill = NA_real_
  ) %>%
  select(-Compartment, -Treatment.Fin) %>%
  as.matrix()

# Impute NAs with column means (so dist/hclust won’t fail)
if (ncol(mat_all) > 0) {
  cm <- colMeans(mat_all, na.rm = TRUE)
  for (j in seq_len(ncol(mat_all))) {
    idx <- is.na(mat_all[, j])
    mat_all[idx, j] <- ifelse(is.finite(cm[j]), cm[j], 0)
  }
}

# Z-score per variable (column)
if (ncol(mat_all) > 0) mat_all <- scale(mat_all)

# Hierarchical clustering on variables (columns)
if (ncol(mat_all) >= 2) {
  hc_cols   <- hclust(dist(t(mat_all)), method = "ward.D2")
  var_order <- colnames(mat_all)[hc_cols$order]
} else {
  var_order <- colnames(mat_all)
}

# ---------------------------
# 2) Cluster TREATMENTS within each Compartment (row order per facet)
# ---------------------------
row_orders <- dat %>%
  group_by(Compartment) %>%
  group_modify(~{
    # wide: rows = Treatment.Fin, cols = Variables
    wide <- .x %>%
      select(Treatment.Fin, Variable, Value_log) %>%
      pivot_wider(
        id_cols    = Treatment.Fin,
        names_from = Variable,
        values_from = Value_log,
        values_fn   = mean,
        values_fill = NA_real_
      )
    
    rn  <- wide$Treatment.Fin
    mat <- wide %>% select(-Treatment.Fin) %>% as.matrix()
    
    # Align columns to global var_order
    missing <- setdiff(var_order, colnames(mat))
    if (length(missing) > 0) {
      mat <- cbind(
        mat,
        matrix(NA_real_, nrow(mat), length(missing),
               dimnames = list(NULL, missing))
      )
    }
    mat <- mat[, var_order, drop = FALSE]
    
    # Impute NAs with column means (within this compartment), fallback to 0
    if (ncol(mat) > 0) {
      cm <- colMeans(mat, na.rm = TRUE)
      for (j in seq_len(ncol(mat))) {
        idx <- is.na(mat[, j])
        mat[idx, j] <- ifelse(is.finite(cm[j]), cm[j], 0)
      }
      # Z-score per variable (column)
      mat <- scale(mat)
    }
    
    # Cluster rows if possible; else keep original order
    ord <- if (nrow(mat) >= 2) {
      hc <- hclust(dist(mat), method = "ward.D2")
      rn[hc$order]
    } else rn
    
    tibble(Treatment.Fin = ord,
           order_y      = seq_along(ord))
  }) %>%
  ungroup()

# ---------------------------
# 3) Apply the orders to your plotting data
#    - Variables: global clustered order
#    - Treatment.Fin: clustered within each Compartment
#    - Use a panel-specific y factor so each facet can have its own order
# ---------------------------
plot_dat <- dat %>%
  left_join(row_orders, by = c("Compartment", "Treatment.Fin")) %>%
  mutate(
    Variable = factor(Variable, levels = var_order),
    # Panel-specific y (so each facet can have its own order)
    TF_panel = paste(Compartment, Treatment.Fin, sep = " | "),
    TF_panel = forcats::fct_reorder(TF_panel, order_y, .desc = FALSE)
  )

# ---------------------------
# 4) Plot (labels hide the compartment prefix on the y-axis)
# ---------------------------
library(ggplot2)

ggplot(plot_dat, aes(x = Variable, y = TF_panel, fill = Value_log)) +
  geom_tile(color = "black") +
  facet_grid(rows = vars(Compartment), scales = "free_y", space = "free_y") +
  scale_y_discrete(labels = function(x) sub("^.*\\|\\s*", "", x)) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x   = element_text(angle = 45, hjust = 1),
    strip.text.y  = element_text(face = "bold"),
    panel.spacing = unit(0.2, "lines")
  ) + scale_fill_gradient2(
    low = "white", high = "red3",
    name = "log(Value)") + xlab("") + ylab("") +
  scale_x_discrete(labels=pretty_map)

# save plot
ggsave(
  filename = "FAPROTAX.HeatMap.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 12.5,
  height = 6,
  units = c("in"),
  dpi = 300)

### ANOVA MODEL FOR ALL 47 FAPROTAX VARIABLES

# Assuming your data frame is called df
# and has columns: Treatment, Compartment, and 47 numeric variables

fun3 <- subset(fun2, Compartment !="Control")
### Run compartments separately to account for differences in treatments (i.e., plastisphere has two treatments, bulk soi = 4)
fun3.p <- subset(fun3, Compartment == "Plastisphere")
fun3.bs <- subset(fun3, Compartment == "Bulk Soil")

# List of quantitative variables
num_vars.p <- fun3.p %>%
  dplyr::select(where(is.numeric)) %>%
  colnames()

# Function to run ANOVA and extract p-values
get_pvals.p <- function(var) {
  model.p <- aov(reformulate(c("Treatment"), response = var), data = fun3.p)
  anova_df.p <- as.data.frame(summary(model.p)[[1]])
  data.frame(
    Variable = var,
    Treatment_p = anova_df.p["Treatment", "Pr(>F)"]
  )
}

# Apply across all 47 numeric variables
results.p <- do.call(rbind, lapply(num_vars.p, get_pvals.p))

# Optional: adjust p-values (e.g., Benjamini-Hochberg FDR)
results.p <- results.p %>%
  mutate(
    Treatment_p_adj = p.adjust(Treatment_p, method = "BH"))

# View the results
head(results.p)

#### BUlk soils ####
# List of quantitative variables
num_vars.bs <- fun3.bs %>%
  dplyr::select(where(is.numeric)) %>%
  colnames()

# Function to run ANOVA and extract p-values
get_pvals.bs <- function(var) {
  model.bs <- aov(reformulate(c("Treatment"), response = var), data = fun3.bs)
  anova_df.bs <- as.data.frame(summary(model.bs)[[1]])
  data.frame(
    Variable = var,
    Treatment_p = anova_df.bs["Treatment", "Pr(>F)"]
  )
}

# Apply across all 47 numeric variables
results.bs <- do.call(rbind, lapply(num_vars.bs, get_pvals.bs))

# Optional: adjust p-values (e.g., Benjamini-Hochberg FDR)
results.bs <- results.bs %>%
  mutate(
    Treatment_p_adj = p.adjust(Treatment_p, method = "BH"))

# View the results
head(results.bs)

### PLOTTING GLOBAL COMPs
library(ggplot2)

results_long.bs <- results.bs %>%
  tidyr::pivot_longer(cols = ends_with("_p_adj"),
                      names_to = "Effect",
                      values_to = "p_value")

results_long.bs$title <- "Bulk Soil FAPROTAX Functional Estimates"


pretty_map <- c(
  methanotrophy = "Methanotrophy",
  methylotrophy = "Methylotrophy",
  aerobic_ammonia_oxidation = "Aerobic Ammonia Oxidation",
  aerobic_nitrite_oxidation = "Aerobic Nitrite Oxidation",
  nitrification = "Nitrification",
  nitrate_denitrification = "Nitrate Denitrification",
  nitrite_denitrification = "Nitrite Denitrification",
  denitrification = "Denitrification",
  chitinolysis = "Chitinolysis",
  nitrogen_fixation = "Nitrogen Fixation",
  nitrite_respiration = "Nitrite Respiration",
  cellulolysis = "Cellulolysis",
  xylanolysis = "Xylanolysis",
  fermentation = "Fermentation",
  aerobic_chemoheterotrophy = "Aerobic Chemoheterotrophy",
  human_pathogens_all = "Human Pathogens (All)",
  human_associated = "Human Associated",
  animal_parasites_or_symbionts = "Animal Parasites Or Symbionts",
  aromatic_compound_degradation = "Aromatic Compound Degradation",
  hydrocarbon_degradation = "Hydrocarbon Degradation",
  nitrate_reduction = "Nitrate Reduction",
  predatory_or_exoparasitic = "Predatory Or Exoparasitic",
  nonphotosynthetic_cyanobacteria = "Nonphotosynthetic Cyanobacteria",
  'anoxygenic_photoautotrophy_s-oxidizing' = "Anoxygenic Photoautotrophy S-Oxidizing",
  anoxygenic_photoautotrophy = "Anoxygenic Photoautotrophy",
  photoautotrophy = "Photoautotrophy",
  aerobic_anoxygenic_phototrophy = "Aerobic Anoxygenic Phototrophy",
  ureolysis = "Ureolysis",
  chemoheterotrophy = "Chemoheterotrophy",
  other = "Other",
  anoxygenic_photoautotrophy_S_oxidizing = "Anoxygenic Photoautotrophy (S Oxidizing)",
  nitrogen_respiration = "Nitrogen Respiration",
  nitrate_respiration = "Nitrate Respiration",
  intracellular_parasites = "Intracellular Parasites",
  dark_oxidation_of_sulfur_compounds = "Dark Oxidation of Sulfur Compounds",
  mammal_gut = "Mammal Gut",
  human_gut = "Human Gut",
  manganese_oxidation = "Manganese Oxidation",
  aromatic_hydrocarbon_degradation = "Aromatic Hydrocarbon Degradation",
  aliphatic_non_methane_hydrocarbon_degradation = "Aliphatic Non-Methane Hydrocarbon Degradation",
  plant_pathogen = "Plant Pathogen",
  dark_thiosulfate_oxidation = "Dark Thiosulfate Oxidation",
  dark_iron_oxidation = "Dark Iron Oxidation",
  anoxygenic_photoautotrophy_H2_oxidizing = "Anoxygenic Phototrophy H2-Oxidizing",
  anoxygenic_photoautotrophy_Fe_oxidizing = "Anoxygenic Photoautotrophy Fe-Oxidizing",
  iron_respiration = "Iron Respiration",
  phototrophy = "Phototrophy",
  photoheterotrophy = "Photoheterotrophy"
)


ggplot(results_long.bs, aes(x = reorder(Variable,-log10(p_value)), y =-log10(p_value), fill = Effect)) +
  geom_col(position = "dodge", fill="burlywood4") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  theme_bw(base_size = 12) + xlab("") + ylab("-log10[p value]") +
  theme(axis.text.x = element_text(hjust = 1)) + facet_grid(~title) +
  theme(strip.background = element_rect(fill="burlywood4"),
        strip.text = element_text(color="white", size=10)) + coord_flip() +
  scale_x_discrete(labels=pretty_map)


# save plot
ggsave(
  filename = "ANOVA.FAPROTAX.BulkSoil.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 6.2,
  height = 7,
  units = c("in"),
  dpi = 300)

sig_vars <- results_long.bs %>%
  filter(p_value < 0.05) %>%
  pull(Variable)
df_bulk <- fun3.bs %>%
  filter(Compartment == "Bulk Soil")

treatment_means <- df_bulk %>%
  select(Treatment, all_of(sig_vars)) %>%
  pivot_longer(-Treatment, names_to = "Variable", values_to = "Value") %>%
  group_by(Variable, Treatment) %>%
  summarise(
    Mean = mean(Value, na.rm = TRUE),
    SE = sd(Value, na.rm = TRUE) / sqrt(n())
  ) %>%
  ungroup()
top_treatments <- treatment_means %>%
  group_by(Variable) %>%
  filter(Mean == max(Mean, na.rm = TRUE)) %>%
  arrange(Variable)
library(ggplot2)

ggplot(treatment_means, aes(x = Treatment, y = Variable, fill = log(Mean))) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(name = "Mean") +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text.y = element_text(face = "bold")
  ) + scale_y_discrete(labels=pretty_map)


### not p-adjusting since effectively a t.test
results_long.p <- results.p %>%
  tidyr::pivot_longer(cols = ends_with("_p"),
                      names_to = "Effect",
                      values_to = "p_value")

results_long.p$title <- "Plastisphere FAPROTAX Functional Estimates"
results_long.p <- subset(results_long.p, Variable !="plant_pathogen")

ggplot(results_long.p, aes(x = reorder(Variable,-log10(p_value)), y =-log10(p_value), fill = Effect)) +
  geom_col(position = "dodge", fill="azure3") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  theme_bw(base_size = 12) + xlab("") + ylab("-log10[p value]") +
  theme(axis.text.x = element_text(hjust = 1)) + facet_grid(~title) +
  theme(strip.background = element_rect(fill="azure3"),
        strip.text = element_text(color="black", size=9.5)) + coord_flip() +
  scale_x_discrete(labels=pretty_map)


# save plot
ggsave(
  filename = "ANOVA.FAPROTAX.Plastisphere.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 6.2,
  height = 7,
  units = c("in"),
  dpi = 300)

### BOX PLOTS ###

# Assign Tukey test letter codes to illustrate significant differences
library(multcompView)
library(tidyverse)

### Remove controls
fun2 <- subset(fun2, Treatment.Fin !="Control")
### Subset plastisphere and bulk soil communities
fun2.p <- subset(fun2, Compartment == "Plastisphere")
fun2.bs <- subset(fun2, Compartment == "Bulk Soil")

#### nitrification ####
nitri.letter <- data.frame(multcompLetters(TukeyHSD(aov(nitrification ~ Treatment.Fin, data = fun2.bs))$Treatment.Fin[,4])$Letters)
colnames(nitri.letter)[1] <- "Letter" #Reassign column name
nitri.letter$Treatment.Fin <- rownames(nitri.letter) #Create column based on rownames
fun2$nitri.facet <- "Nitrification"
# plot
ggplot(fun2.bs, aes(x=Treatment.Fin, y=nitrification, color=Treatment.Fin)) + 
  geom_boxplot() + theme_bw() + xlab("") + geom_text(data=nitri.letter, 
                                                     aes(label=Letter, x=Treatment.Fin, y = 40, size=12)) +
  ylab("Nitrification") + geom_jitter(alpha=0.4) + 
  scale_color_manual(values = c("slategray4", "red3", "goldenrod4", "#003C30")) + 
  theme(legend.position="none") + theme(axis.text.y = element_text(size=12, face="bold")) + 
  theme(axis.text.x = element_text(size = 10, angle=45, hjust=1)) + 
  theme(axis.title.x = element_text(size = 12)) +
  facet_grid(~nitri.facet) +
  theme(strip.text = element_text(size=14, face="bold", color="white")) + 
  theme(strip.background = element_rect(fill="burlywood4"))

# save plot
ggsave(
  filename = "Nitrification.BulkSoil.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 4,
  height = 5,
  units = c("in"),
  dpi = 300)

comps <- c("Plastic / No Fungus", "Plastic+Fungus")
t.test(nitrification ~ Treatment.Fin, data = fun2.p) ### t.test instead of AOV since only two groups
### Nitrification Plastisphere
ggplot(fun2.p, aes(x=Treatment.Fin, y=nitrification, color=Treatment.Fin)) + 
  geom_boxplot() + theme_bw() + xlab("") +
  ylab("Nitrification") + geom_jitter(alpha=0.4) + 
  scale_color_manual(values = c("goldenrod4", "#003C30")) + 
  theme(legend.position="none") + theme(axis.text.y = element_text(size=12, face="bold")) + 
  theme(axis.text.x = element_text(size = 10, angle=45, hjust=1)) + 
  theme(axis.title.x = element_text(size = 12)) +
  facet_grid(~nitri.facet) +
  theme(strip.text = element_text(size=14, face="bold", color="black")) + 
  theme(strip.background = element_rect(fill="azure3")) + 
  annotate("text", x = 1.5, y = 8, label = "p = 0.3621", size = 5)

# save plot
ggsave(
  filename = "Nitrification.Plastisphere.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 4,
  height = 5,
  units = c("in"),
  dpi = 300)

#### N fixation ####
nfix.letter <- data.frame(multcompLetters(TukeyHSD(aov(nitrogen_fixation ~ Treatment.Fin, data = fun2.bs))$Treatment.Fin[,4])$Letters)
colnames(nfix.letter)[1] <- "Letter" #Reassign column name
nfix.letter$Treatment.Fin <- rownames(nfix.letter) #Create column based on rownames
fun2.bs$nfix.facet <- "Nitrogen Fixation"
# plot
ggplot(fun2.bs, aes(x=Treatment.Fin, y=nitrogen_fixation, color=Treatment.Fin)) + 
  geom_boxplot() + theme_bw() + xlab("") + geom_text(data=nfix.letter, 
                                                     aes(label=Letter, x=Treatment.Fin, y = 250, size=12)) +
  ylab("Nitrogen Fixation") + geom_jitter(alpha=0.4) + 
  scale_color_manual(values = c("slategray4", "red3", "goldenrod4", "#003C30")) + 
  theme(legend.position="none") + theme(axis.text.y = element_text(size=12, face="bold")) + 
  theme(axis.text.x = element_text(size = 10, angle=45, hjust=1)) + 
  theme(axis.title.x = element_text(size = 12)) +
  facet_grid(~nfix.facet) +
  theme(strip.text = element_text(size=14, face="bold", color="white")) + 
  theme(strip.background = element_rect(fill="burlywood4"))

# save plot
ggsave(
  filename = "N.Fixation.bulksoil.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 4,
  height = 5,
  units = c("in"),
  dpi = 300)

### nfixation Plastisphere
t.test(nitrogen_fixation ~ Treatment.Fin, data=fun2.p)
ggplot(fun2.p, aes(x=Treatment.Fin, y=nitrogen_fixation, color=Treatment.Fin)) + 
  geom_boxplot() + theme_bw() + xlab("") +
  ylab("Nitrogen Fixation") + geom_jitter(alpha=0.4) + 
  scale_color_manual(values = c("goldenrod4", "#003C30")) + 
  theme(legend.position="none") + theme(axis.text.y = element_text(size=12, face="bold")) + 
  theme(axis.text.x = element_text(size = 10, angle=45, hjust=1)) + 
  theme(axis.title.x = element_text(size = 12)) +
  facet_grid(~nfix.facet) +
  theme(strip.text = element_text(size=14, face="bold", color="black")) + 
  theme(strip.background = element_rect(fill="azure3")) + 
  annotate("text", x = 1.5, y = 300, label = "p = 0.5805", size = 5)

# save plot
ggsave(
  filename = "NFixation.Plastisphere.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 4,
  height = 5,
  units = c("in"),
  dpi = 300)

##### Human pathogens all #####

HPA.letter <- data.frame(multcompLetters(TukeyHSD(aov(human_pathogens_all ~ Treatment.Fin, data = fun2.bs))$Treatment.Fin[,4])$Letters)
colnames(HPA.letter)[1] <- "Letter" #Reassign column name
HPA.letter$Treatment.Fin <- rownames(HPA.letter) #Create column based on rownames
fun2.bs$HPA.facet <- "Human Pathogens (Bulk Soil)"
# plot
ggplot(fun2.bs, aes(x=Treatment.Fin, y=human_pathogens_all, color=Treatment.Fin)) + 
  geom_boxplot() + theme_bw() + xlab("") + geom_text(data=HPA.letter, 
                                                     aes(label=Letter, x=Treatment.Fin, y = 35, size=12)) +
  ylab("Human Pathogens") + geom_jitter(alpha=0.4) + 
  scale_color_manual(values = c("slategray4", "red3", "goldenrod4", "#003C30")) + 
  theme(legend.position="none") + theme(axis.text.y = element_text(size=12, face="bold")) + 
  theme(axis.text.x = element_text(size = 10, angle=45, hjust=1)) + 
  theme(axis.title.x = element_text(size = 12)) +
  facet_grid(~HPA.facet) +
  theme(strip.text = element_text(size=14, face="bold", color="white")) + 
  theme(strip.background = element_rect(fill="burlywood4"))

# save plot
ggsave(
  filename = "HumanPATHOGENS.bulksoil.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 4,
  height = 5,
  units = c("in"),
  dpi = 300)

### Human Pathogens Plastisphere
t.test(human_pathogens_all ~ Treatment.Fin, data=fun2.p)
fun2.p$HPA.facet <- "Human Pathogens (Plastisphere)"
ggplot(fun2.p, aes(x=Treatment.Fin, y=human_pathogens_all, color=Treatment.Fin)) + 
  geom_boxplot() + theme_bw() + xlab("") +
  ylab("Human Pathogens") + geom_jitter(alpha=0.4) + 
  scale_color_manual(values = c("goldenrod4", "#003C30")) + 
  theme(legend.position="none") + theme(axis.text.y = element_text(size=12, face="bold")) + 
  theme(axis.text.x = element_text(size = 10, angle=45, hjust=1)) + 
  theme(axis.title.x = element_text(size = 12)) +
  facet_grid(~HPA.facet) +
  theme(strip.text = element_text(size=14, face="bold", color="black")) + 
  theme(strip.background = element_rect(fill="azure3")) + 
  annotate("text", x = 1.5, y = 35, label = "p = 0.0154", size = 5)

# save plot
ggsave(
  filename = "HumanPATHOS.Plastisphere.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 4,
  height = 5,
  units = c("in"),
  dpi = 300)


#### METHANOTROPHY ####

methano.letter <- data.frame(multcompLetters(TukeyHSD(aov(methanotrophy ~ Treatment.Fin, data = fun2.bs))$Treatment.Fin[,4])$Letters)
colnames(methano.letter)[1] <- "Letter" #Reassign column name
methano.letter$Treatment.Fin <- rownames(methano.letter) #Create column based on rownames
fun2.bs$methano.letter <- "Methanotrophy"
# plot
ggplot(fun2.bs, aes(x=Treatment.Fin, y=methanotrophy, color=Treatment.Fin)) + 
  geom_boxplot() + theme_bw() + xlab("") + geom_text(data=methano.letter, 
                                                     aes(label=Letter, x=Treatment.Fin, y = 4, size=12)) +
  ylab("Methanotrophy") + geom_jitter(alpha=0.4) + 
  scale_color_manual(values = c("slategray4", "red3", "goldenrod4", "#003C30")) + 
  theme(legend.position="none") + theme(axis.text.y = element_text(size=12, face="bold")) + 
  theme(axis.text.x = element_text(size = 10, angle=45, hjust=1)) + 
  theme(axis.title.x = element_text(size = 12)) +
  facet_grid(~HPA.facet) +
  theme(strip.text = element_text(size=14, face="bold", color="white")) + 
  theme(strip.background = element_rect(fill="burlywood4"))

# save plot
ggsave(
  filename = "Methanotrophy.bulksoil.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 4,
  height = 5,
  units = c("in"),
  dpi = 300)

### Human Pathogens Plastisphere
t.test(methanotrophy ~ Treatment.Fin, data=fun2.p)
fun2.p$methano.facet <- "Methanotrophy"
ggplot(fun2.p, aes(x=Treatment.Fin, y=methanotrophy, color=Treatment.Fin)) + 
  geom_boxplot() + theme_bw() + xlab("") +
  ylab("Methanotrophy") + geom_jitter(alpha=0.4) + 
  scale_color_manual(values = c("goldenrod4", "#003C30")) + 
  theme(legend.position="none") + theme(axis.text.y = element_text(size=12, face="bold")) + 
  theme(axis.text.x = element_text(size = 10, angle=45, hjust=1)) + 
  theme(axis.title.x = element_text(size = 12)) +
  facet_grid(~methano.facet) +
  theme(strip.text = element_text(size=14, face="bold", color="black")) + 
  theme(strip.background = element_rect(fill="azure3")) + 
  annotate("text", x = 1.5, y = 4, label = "p = 0.0324", size = 5)

# save plot
ggsave(
  filename = "methanotrophy.Plastisphere.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 4,
  height = 5,
  units = c("in"),
  dpi = 300)

unique(fun_long_NC$Variable)

##### predatory_or_exoparasitic #####


pred.letter <- data.frame(multcompLetters(TukeyHSD(aov(predatory_or_exoparasitic ~ Treatment.Fin, data = fun2.bs))$Treatment.Fin[,4])$Letters)
colnames(pred.letter)[1] <- "Letter" #Reassign column name
pred.letter$Treatment.Fin <- rownames(pred.letter) #Create column based on rownames
fun2.bs$pred.facet <- "Predatory / Exoparasitic"
# plot
ggplot(fun2.bs, aes(x=Treatment.Fin, y=predatory_or_exoparasitic, color=Treatment.Fin)) + 
  geom_boxplot() + theme_bw() + xlab("") + geom_text(data=pred.letter, 
                                                     aes(label=Letter, x=Treatment.Fin, y = 200, size=12)) +
  ylab("Predatory / Exoparasitic") + geom_jitter(alpha=0.4) + 
  scale_color_manual(values = c("slategray4", "red3", "goldenrod4", "#003C30")) + 
  theme(legend.position="none") + theme(axis.text.y = element_text(size=12, face="bold")) + 
  theme(axis.text.x = element_text(size = 10, angle=45, hjust=1)) + 
  theme(axis.title.x = element_text(size = 12)) +
  facet_grid(~pred.facet) +
  theme(strip.text = element_text(size=14, face="bold", color="white")) + 
  theme(strip.background = element_rect(fill="burlywood4"))

# save plot
ggsave(
  filename = "pred.bulksoil.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 4,
  height = 5,
  units = c("in"),
  dpi = 300)

### Human Pathogens Plastisphere
t.test(predatory_or_exoparasitic ~ Treatment.Fin, data=fun2.p)
fun2.p$pred.facet <- "Predatory / Exoparasitic"
ggplot(fun2.p, aes(x=Treatment.Fin, y=predatory_or_exoparasitic, color=Treatment.Fin)) + 
  geom_boxplot() + theme_bw() + xlab("") +
  ylab("Predatory / Exoparasitic") + geom_jitter(alpha=0.4) + 
  scale_color_manual(values = c("goldenrod4", "#003C30")) + 
  theme(legend.position="none") + theme(axis.text.y = element_text(size=12, face="bold")) + 
  theme(axis.text.x = element_text(size = 10, angle=45, hjust=1)) + 
  theme(axis.title.x = element_text(size = 12)) +
  facet_grid(~pred.facet) +
  theme(strip.text = element_text(size=14, face="bold", color="black")) + 
  theme(strip.background = element_rect(fill="azure3")) + 
  annotate("text", x = 1.5, y = 200, label = "p = 0.9374", size = 5)

# save plot
ggsave(
  filename = "Predatory.Plastisphere.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 4,
  height = 5,
  units = c("in"),
  dpi = 300)


#### DESeq2 to compare FAPROTAX outputs ####

## Do bulk vs. plastisphere with and without P. microspora

# ===== PACKAGES =====
library(dplyr)
library(tidyr)
library(tibble)
library(DESeq2)
library(ggplot2)

# ===== 0) START FROM YOUR DATA FRAME: fun3 =====
# Optional: drop any “controls”; keep this or remove it if not needed
# fun3 <- fun3 %>% filter(Treatment != "NA")

w.fungus <- subset(fun3, Treatment.Fin == "Plastic+Fungus" | Treatment.Fin == "Fungus / No Plastic")
wo.fungus <- subset(fun3, Treatment.Fin == "Plastic / No Fungus" | Treatment.Fin == "No Plastic / No Fungus")

# Identify the 47 function columns (numeric)
feature_cols <- w.fungus %>% select(where(is.numeric)) %>% names()

# Basic sanity (unique, non-missing sample IDs)
stopifnot("Sample" %in% names(w.fungus))
stopifnot(!any(is.na(w.fungus$Sample)), !any(duplicated(w.fungus$Sample)))

# ===== 1) Build counts matrix (features x samples) and sample table =====
# (Make sure there are no rownames conflicts)
w.fungus2 <- tibble::remove_rownames(w.fungus)

counts_mat <- w.fungus2 %>%
  select(Sample, all_of(feature_cols)) %>%
  column_to_rownames("Sample") %>%
  t() %>%
  as.matrix()

coldata <- w.fungus2 %>%
  distinct(Sample, .keep_all = TRUE) %>%
  select(Sample, Compartment) %>%
  as.data.frame() %>%
  `rownames<-`(.$Sample)

# Clean compartment factor (e.g., "Bulk Soil" -> "bulk"; "Plastisphere" -> "plastisphere")
coldata$Compartment <- tolower(gsub("[[:space:]]+", "_", coldata$Compartment))
coldata$Compartment <- factor(coldata$Compartment, levels = c("bulk_soil","plastisphere","bulk"))

# If your levels are "bulk_soil" and "plastisphere", set a simple two-level factor:
if (all(c("bulk_soil","plastisphere") %in% levels(coldata$Compartment))) {
  coldata$Compartment <- factor(coldata$Compartment, levels = c("bulk_soil","plastisphere"))
} else if (all(c("bulk","plastisphere") %in% levels(coldata$Compartment))) {
  coldata$Compartment <- factor(coldata$Compartment, levels = c("bulk","plastisphere"))
}

# Align order just in case
counts_mat <- counts_mat[, rownames(coldata), drop = FALSE]

# Quick checks
stopifnot(all(counts_mat >= 0))
stopifnot(all(abs(counts_mat - round(counts_mat)) < .Machine$double.eps^0.5))

# ===== 2) DESeq2 (~ Compartment) =====
dds <- DESeqDataSetFromMatrix(countData = counts_mat,
                              colData   = coldata,
                              design    = ~ Compartment)

# Run DESeq
dds <- DESeq(dds, test = "Wald", fitType = "parametric")

# ===== 3) Results: Bulk vs Plastisphere =====
# Contrast order: c("factor","numerator","denominator")
# Use what's in your factor levels above (either "bulk_soil" or "bulk").
levs <- levels(coldata$Compartment)
bulk_name <- if ("bulk_soil" %in% levs) "bulk_soil" else "bulk"

res_bulk_vs_plast <- results(dds,
                             contrast = c("Compartment", bulk_name, "plastisphere"),
                             cooksCutoff = FALSE)
alpha <- 0.1
sigtab <- as.data.frame(res_bulk_vs_plast) %>%
  rownames_to_column("Function") %>%
  filter(!is.na(padj), padj < alpha)

# Order by effect size for plotting
sigtab <- sigtab %>%
  arrange(desc(log2FoldChange)) %>%
  mutate(Function = factor(Function, levels = rev(Function)),
         label = "With Fungus")

# ===== 4) Plot (dot plot like your example) =====
ggplot(sigtab, aes(x = Function, y = log2FoldChange)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = 2, color = "red") +
  coord_flip() +
  theme_bw() +
  labs(x = NULL, y = "[Log2] Fold-Change") +
  facet_grid(~ label) +
  theme(
    axis.text.y = element_text(size = 10),
    strip.background = element_rect(fill = "gray85"),
    strip.text = element_text(size = 10, face = "bold")
  ) + scale_x_discrete(labels=pretty_map)

# ===== 5) Save plot =====
ggsave("FAPROTAX_WithFungus.DESeq2.png", width = 4, height = 3.5, dpi = 300)

# ===== (Optional) Volcano quickie =====
# ggplot(as.data.frame(res_bulk_vs_plast) %>% rownames_to_column("Function"),
#        aes(x = log2FoldChange, y = -log10(padj))) +
#   geom_point(alpha = 0.6) +
#   geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
#   theme_bw() +
#   labs(title = "Volcano: Bulk vs Plastisphere")

##### Without Fungus added #####

# Identify the 47 function columns (numeric)
feature_cols2 <- wo.fungus %>% select(where(is.numeric)) %>% names()

# Basic sanity (unique, non-missing sample IDs)
stopifnot("Sample" %in% names(wo.fungus))
stopifnot(!any(is.na(wo.fungus$Sample)), !any(duplicated(wo.fungus$Sample)))

# ===== 1) Build counts matrix (features x samples) and sample table =====
# (Make sure there are no rownames conflicts)
wo.fungus2 <- tibble::remove_rownames(wo.fungus)

counts_mat2 <- wo.fungus2 %>%
  select(Sample, all_of(feature_cols2)) %>%
  column_to_rownames("Sample") %>%
  t() %>%
  as.matrix()

coldata2 <- wo.fungus2 %>%
  distinct(Sample, .keep_all = TRUE) %>%
  select(Sample, Compartment) %>%
  as.data.frame() %>%
  `rownames<-`(.$Sample)

# Clean compartment factor (e.g., "Bulk Soil" -> "bulk"; "Plastisphere" -> "plastisphere")
coldata2$Compartment <- tolower(gsub("[[:space:]]+", "_", coldata2$Compartment))
coldata2$Compartment <- factor(coldata2$Compartment, levels = c("bulk_soil","plastisphere","bulk"))

# If your levels are "bulk_soil" and "plastisphere", set a simple two-level factor:
if (all(c("bulk_soil","plastisphere") %in% levels(coldata2$Compartment))) {
  coldata2$Compartment <- factor(coldata2$Compartment, levels = c("bulk_soil","plastisphere"))
} else if (all(c("bulk","plastisphere") %in% levels(coldata2$Compartment))) {
  coldata2$Compartment <- factor(coldata2$Compartment, levels = c("bulk","plastisphere"))
}

# Align order just in case
counts_mat2 <- counts_mat2[, rownames(coldata2), drop = FALSE]

# Quick checks
stopifnot(all(counts_mat2 >= 0))
stopifnot(all(abs(counts_mat2 - round(counts_mat2)) < .Machine$double.eps^0.5))

# ===== 2) DESeq2 (~ Compartment) =====
dds2 <- DESeqDataSetFromMatrix(countData = counts_mat2,
                              colData   = coldata2,
                              design    = ~ Compartment)

# Run DESeq
dds2 <- DESeq(dds2, test = "Wald", fitType = "parametric")

# ===== 3) Results: Bulk vs Plastisphere =====
# Contrast order: c("factor","numerator","denominator")
# Use what's in your factor levels above (either "bulk_soil" or "bulk").
levs2 <- levels(coldata2$Compartment)
bulk_name2 <- if ("bulk_soil" %in% levs2) "bulk_soil" else "bulk"

res_bulk_vs_plast2 <- results(dds2,
                             contrast = c("Compartment", bulk_name2, "plastisphere"),
                             cooksCutoff = FALSE)
alpha <- 0.1
sigtab2 <- as.data.frame(res_bulk_vs_plast2) %>%
  rownames_to_column("Function") %>%
  filter(!is.na(padj), padj < alpha)

# Order by effect size for plotting
sigtab2 <- sigtab2 %>%
  arrange(desc(log2FoldChange)) %>%
  mutate(Function = factor(Function, levels = rev(Function)),
         label = "Without Fungus")

# ===== 4) Plot (dot plot like your example) =====
ggplot(sigtab2, aes(x = Function, y = log2FoldChange)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = 2, color = "red") +
  coord_flip() +
  theme_bw() +
  labs(x = NULL, y = "[Log2] Fold-Change") +
  facet_grid(~ label) +
  theme(
    axis.text.y = element_text(size = 10),
    strip.background = element_rect(fill = "gray85"),
    strip.text = element_text(size = 10, face = "bold")
  ) + scale_x_discrete(labels=pretty_map)

# ===== 5) Save plot =====
ggsave("FAPROTAX_WithoutFungus.DESeq2.png", width = 4, height = 3.5, dpi = 300)

# ===== (Optional) Volcano quickie =====
# ggplot(as.data.frame(res_bulk_vs_plast) %>% rownames_to_column("Function"),
#        aes(x = log2FoldChange, y = -log10(padj))) +
#   geom_point(alpha = 0.6) +
#   geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
#   theme_bw() +
#   labs(title = "Volcano: Bulk vs Plastisphere")

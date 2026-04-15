### Plastic degradation ###

library(ggplot2)
library(ggpubr)
df <- read.csv("Plastic_Harvest.csv")
df2 <- subset(df, Analytical.Weight !="NA")
df2$Title <- "Soil Microcosm"

my_comps <- list(c("Plastic", "Plastic+Fungus"))
ggplot(df2, aes(x=Facet_Treatment, y=Analytical.Weight, color=Facet_Treatment)) +
  geom_boxplot(color=c("Darkgoldenrod4", "Cornsilk4")) +
  geom_jitter(alpha=0.4) + 
  scale_color_manual(values = c("Darkgoldenrod4", "Cornsilk4")) +
  theme_bw() + facet_grid(~Title) + 
  stat_compare_means(comparisons=my_comps) +
  xlab("") + ylab("Mass (g)") +
  theme(axis.text.x = element_text(face="bold", size=12)) +
  theme(strip.background = element_rect(fill = "burlywood4")) +
  theme(strip.text = element_text(size=12, face="bold", color="white")) +
  theme(axis.title.y = element_text(vjust = 1.5)) +
  theme(legend.position = 'none')

# Save plot
ggsave(
  filename = "Soil.PlasticDegradation.BoxPlotComp.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 4,
  height = 4.5,
  units = c("in"),
  dpi = 300)

#### Test Tube Experiment ####
tt.data <- read.csv("InVitro.data.csv")

tt.data$Title <- "In Vitro Microcosm"
my_comps2 <- list(c("Control", "P. microspora"))
ggplot(tt.data, aes(x=Treatment, y=End.Mass, color=Treatment)) +
  geom_boxplot(color=c("Darkgoldenrod4", "Cornsilk4")) +
  geom_jitter(alpha=0.4) + 
  scale_color_manual(values = c("Darkgoldenrod4", "Cornsilk4")) +
  theme_bw() + facet_grid(~Title) + 
  stat_compare_means(comparisons=my_comps2) +
  xlab("") + ylab("Mass (g)") +
  theme(axis.text.x = element_text(face="bold", size=12)) +
  theme(strip.background = element_rect(fill = "gray")) +
  theme(strip.text = element_text(size=12, face="bold")) +
  theme(axis.title.y = element_text(vjust = 1.5)) +
  theme(legend.position = "none") + 
  scale_x_discrete(labels=c("Plastic", "Plastic+Fungus"))

# Save plot
ggsave(
  filename = "TT.PlasticDegradation.BoxPlotComp.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 4,
  height = 4.5,
  units = c("in"),
  dpi = 300)

library(multcompView)
vit.Letter <- data.frame(multcompLetters(TukeyHSD(aov(mass.loss ~ Treatment, data = tt.data))$Treatment[,2])$Letters)
colnames(vit.Letter)[1] <- "Letter" #Reassign column name
vit.Letter$Treatment <- rownames(vit.Letter) #Create column based on rownames
vit.placement <- tt.data %>% #We want to create a dataframe to assign the letter position.
  group_by(Treatment) %>%
  summarise(mean=mean(mass.loss), sd=sd(mass.loss)) %>%
  arrange(desc(mean))
# merge the letters and placement dataframe for plotting
vit.letter2 <- left_join(vit.Letter, vit.placement)
# plot
ggplot(tt.data, aes(x=Treatment, y=mass.loss, color=Treatment)) + 
  geom_boxplot() + theme_bw() + xlab("") + geom_text(data=vit.letter2, 
                                                     aes(label=Letter, x=Treatment, y = 0.2, size=12)) +
  ylab("Mass (g)") + geom_jitter(alpha=0.4) + 
  scale_color_manual(values = c("slategray4", "#003C30")) + 
  theme(legend.position="none") + theme(axis.text.y = element_text(size=12, face="bold")) + 
  theme(axis.text.x = element_text(size = 10, angle=45, hjust=1)) + theme(axis.title.x = element_text(size = 12)) +
  facet_grid(~Title) +
  theme(strip.text = element_text(size=14, face="bold", color="white")) + 
  theme(strip.background = element_rect(fill="gray"))
# save plot
ggsave(
  filename = "Soil.TEC.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 4,
  height = 5,
  units = c("in"),
  dpi = 300)

#### Soil Chemistry Analyses

soil.chem <- read.csv("SoilChemFinal.csv")
soil.chem <- subset(soil.chem, pH != "NA")

# Assign Tukey test letter codes to illustrate significant differences
### install.packages("multcompView") ###
library(multcompView)
library(tidyverse)
library(ggplot2)

#### pH ####
pH.letter <- data.frame(multcompLetters(TukeyHSD(aov(pH ~ Treatment, data = soil.chem))$Treatment[,4])$Letters)
colnames(pH.letter)[1] <- "Letter" #Reassign column name
pH.letter$Treatment <- rownames(pH.letter) #Create column based on rownames
pH.placement <- soil.chem %>% #We want to create a dataframe to assign the letter position.
  group_by(Treatment) %>%
  summarise(mean=mean(pH), sd=sd(pH)) %>%
  arrange(desc(mean))
# merge the letters and placement dataframe for plotting
pH.letter2 <- left_join(pH.letter, pH.placement)
soil.chem$pH.facet <- "pH"
# plot
ggplot(soil.chem, aes(x=Treatment, y=pH, color=Treatment)) + 
  geom_boxplot() + theme_bw() + xlab("") + geom_text(data=pH.letter2, 
                                                     aes(label=Letter, x=Treatment, y = 7, size=12)) +
  ylab("pH") + geom_jitter(alpha=0.4) + 
  scale_color_manual(values = c("slategray4", "red3", "goldenrod4", "#003C30")) + 
  theme(legend.position="none") + theme(axis.text.y = element_text(size=12, face="bold")) + 
  theme(axis.text.x = element_text(size = 10, angle=45, hjust=1)) + theme(axis.title.x = element_text(size = 12)) +
  facet_grid(~pH.facet) +
  theme(strip.text = element_text(size=14, face="bold", color="white")) + 
  theme(strip.background = element_rect(fill="burlywood4"))
# save plot
ggsave(
  filename = "Soil.pH.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 4,
  height = 5,
  units = c("in"),
  dpi = 300)

### Total Exchange Capacity

TEC.letter <- data.frame(multcompLetters(TukeyHSD(aov(TEC ~ Treatment, data = soil.chem))$Treatment[,4])$Letters)
colnames(TEC.letter)[1] <- "Letter" #Reassign column name
TEC.letter$Treatment <- rownames(TEC.letter) #Create column based on rownames
TEC.placement <- soil.chem %>% #We want to create a dataframe to assign the letter position.
  group_by(Treatment) %>%
  summarise(mean=mean(TEC), sd=sd(TEC)) %>%
  arrange(desc(mean))
# merge the letters and placement dataframe for plotting
TEC.letter2 <- left_join(TEC.letter, TEC.placement)
soil.chem$TEC.facet <- "Total Exchange Capacity"
# plot
ggplot(soil.chem, aes(x=Treatment, y=TEC, color=Treatment)) + 
  geom_boxplot() + theme_bw() + xlab("") + geom_text(data=TEC.letter2, 
                                                     aes(label=Letter, x=Treatment, y = 10.5, size=12)) +
  ylab("Total Exchange Capacity (meq/100 g)") + geom_jitter(alpha=0.4) + 
  scale_color_manual(values = c("slategray4", "red3", "goldenrod4", "#003C30")) + 
  theme(legend.position="none") + theme(axis.text.y = element_text(size=12, face="bold")) + 
  theme(axis.text.x = element_text(size = 10, angle=45, hjust=1)) + theme(axis.title.x = element_text(size = 12)) +
  facet_grid(~TEC.facet) +
  theme(strip.text = element_text(size=14, face="bold", color="white")) + 
  theme(strip.background = element_rect(fill="burlywood4"))
# save plot
ggsave(
  filename = "Soil.TEC.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 4,
  height = 5,
  units = c("in"),
  dpi = 300)

#### Percent Organic Matter ####

OM.Letter <- data.frame(multcompLetters(TukeyHSD(aov(Perc.OM ~ Treatment, data = soil.chem))$Treatment[,4])$Letters)
colnames(OM.Letter)[1] <- "Letter" #Reassign column name
OM.Letter$Treatment <- rownames(OM.Letter) #Create column based on rownames
OM.placement <- soil.chem %>% #We want to create a dataframe to assign the letter position.
  group_by(Treatment) %>%
  summarise(mean=mean(Perc.OM), sd=sd(Perc.OM)) %>%
  arrange(desc(mean))
# merge the letters and placement dataframe for plotting
OM.letter2 <- left_join(OM.Letter, OM.placement)
soil.chem$OM.facet <- "Organic Matter (%)"
# plot
ggplot(soil.chem, aes(x=Treatment, y=Perc.OM, color=Treatment)) + 
  geom_boxplot() + theme_bw() + xlab("") + geom_text(data=OM.letter2, 
                                                     aes(label=Letter, x=Treatment, y = 3.5, size=12)) +
  ylab("Organic Matter (%)") + geom_jitter(alpha=0.4) + 
  scale_color_manual(values = c("slategray4", "red3", "goldenrod4", "#003C30")) + 
  theme(legend.position="none") + theme(axis.text.y = element_text(size=12, face="bold")) + 
  theme(axis.text.x = element_text(size = 10, angle=45, hjust=1)) + theme(axis.title.x = element_text(size = 12)) +
  facet_grid(~OM.facet) +
  theme(strip.text = element_text(size=14, face="bold", color="white")) + 
  theme(strip.background = element_rect(fill="burlywood4"))
# save plot
ggsave(
  filename = "Soil.OM.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 4,
  height = 5,
  units = c("in"),
  dpi = 300)

#### Estimated Nitrogen Release ####

Nre.Letter <- data.frame(multcompLetters(TukeyHSD(aov(Est.N.Release ~ Treatment, data = soil.chem))$Treatment[,4])$Letters)
colnames(Nre.Letter)[1] <- "Letter" #Reassign column name
Nre.Letter$Treatment <- rownames(Nre.Letter) #Create column based on rownames
Nre.placement <- soil.chem %>% #We want to create a dataframe to assign the letter position.
  group_by(Treatment) %>%
  summarise(mean=mean(Est.N.Release), sd=sd(Est.N.Release)) %>%
  arrange(desc(mean))
# merge the letters and placement dataframe for plotting
Nre.Letter2 <- left_join(Nre.Letter, Nre.placement)
soil.chem$Nre.facet <- "Estimated Nitrogen Release"
# plot
ggplot(soil.chem, aes(x=Treatment, y=Est.N.Release, color=Treatment)) + 
  geom_boxplot() + theme_bw() + xlab("") + geom_text(data=Nre.Letter2, 
                                                     aes(label=Letter, x=Treatment, y = 90, size=12)) +
  ylab("Estimated N Release (N/acre)") + geom_jitter(alpha=0.4) + 
  scale_color_manual(values = c("slategray4", "red3", "goldenrod4", "#003C30")) + 
  theme(legend.position="none") + theme(axis.text.y = element_text(size=12, face="bold")) + 
  theme(axis.text.x = element_text(size = 10, angle=45, hjust=1)) + theme(axis.title.x = element_text(size = 12)) +
  facet_grid(~Nre.facet) +
  theme(strip.text = element_text(size=14, face="bold", color="white")) + 
  theme(strip.background = element_rect(fill="burlywood4"))
# save plot
ggsave(
  filename = "Soil.N.Release.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 4,
  height = 5,
  units = c("in"),
  dpi = 300)


#### Sulfur PPM ####

S.Letter <- data.frame(multcompLetters(TukeyHSD(aov(S.PPM ~ Treatment, data = soil.chem))$Treatment[,4])$Letters)
colnames(S.Letter)[1] <- "Letter" #Reassign column name
S.Letter$Treatment <- rownames(S.Letter) #Create column based on rownames
S.placement <- soil.chem %>% #We want to create a dataframe to assign the letter position.
  group_by(Treatment) %>%
  summarise(mean=mean(S.PPM), sd=sd(S.PPM)) %>%
  arrange(desc(mean))
# merge the letters and placement dataframe for plotting
S.Letter2 <- left_join(S.Letter, S.placement)
soil.chem$S.facet <- "Sulfur PPM"
# plot
ggplot(soil.chem, aes(x=Treatment, y=S.PPM, color=Treatment)) + 
  geom_boxplot() + theme_bw() + xlab("") + geom_text(data=S.Letter2, 
                                                     aes(label=Letter, x=Treatment, y = 30, size=12)) +
  ylab("Sulfur (PPM)") + geom_jitter(alpha=0.4) + 
  scale_color_manual(values = c("slategray4", "red3", "goldenrod4", "#003C30")) + 
  theme(legend.position="none") + theme(axis.text.y = element_text(size=12, face="bold")) + 
  theme(axis.text.x = element_text(size = 10, angle=45, hjust=1)) + theme(axis.title.x = element_text(size = 12)) +
  facet_grid(~S.facet) +
  theme(strip.text = element_text(size=14, face="bold", color="white")) + 
  theme(strip.background = element_rect(fill="burlywood4"))
# save plot
ggsave(
  filename = "Soil.S.PPM.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 4,
  height = 5,
  units = c("in"),
  dpi = 300)

#### Phosphorus mg/kg or PPM ####

P.Letter <- data.frame(multcompLetters(TukeyHSD(aov(P.mg.kg ~ Treatment, data = soil.chem))$Treatment[,4])$Letters)
colnames(P.Letter)[1] <- "Letter" #Reassign column name
P.Letter$Treatment <- rownames(P.Letter) #Create column based on rownames
P.placement <- soil.chem %>% #We want to create a dataframe to assign the letter position.
  group_by(Treatment) %>%
  summarise(mean=mean(P.mg.kg), sd=sd(P.mg.kg)) %>%
  arrange(desc(mean))
# merge the letters and placement dataframe for plotting
P.Letter2 <- left_join(P.Letter, P.placement)
soil.chem$P.facet <- "Phosphorus PPM"
# plot
ggplot(soil.chem, aes(x=Treatment, y=P.mg.kg, color=Treatment)) + 
  geom_boxplot() + theme_bw() + xlab("") + geom_text(data=P.Letter2, 
                                                     aes(label=Letter, x=Treatment, y = 15, size=12)) +
  ylab("Phosphorus (PPM)") + geom_jitter(alpha=0.4) + 
  scale_color_manual(values = c("slategray4", "red3", "goldenrod4", "#003C30")) + 
  theme(legend.position="none") + theme(axis.text.y = element_text(size=12, face="bold")) + 
  theme(axis.text.x = element_text(size = 10, angle=45, hjust=1)) + theme(axis.title.x = element_text(size = 12)) +
  facet_grid(~P.facet) +
  theme(strip.text = element_text(size=14, face="bold", color="white")) + 
  theme(strip.background = element_rect(fill="burlywood4"))
# save plot
ggsave(
  filename = "Soil.P.PPM.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 4,
  height = 5,
  units = c("in"),
  dpi = 300)

#### Ca.mg.kg / PPM ####

Ca.Letter <- data.frame(multcompLetters(TukeyHSD(aov(Ca.mg.kg ~ Treatment, data = soil.chem))$Treatment[,4])$Letters)
colnames(Ca.Letter)[1] <- "Letter" #Reassign column name
Ca.Letter$Treatment <- rownames(Ca.Letter) #Create column based on rownames
Ca.placement <- soil.chem %>% #We want to create a dataframe to assign the letter position.
  group_by(Treatment) %>%
  summarise(mean=mean(Ca.mg.kg), sd=sd(Ca.mg.kg)) %>%
  arrange(desc(mean))
# merge the letters and placement dataframe for plotting
Ca.Letter2 <- left_join(Ca.Letter, Ca.placement)
soil.chem$Ca.facet <- "Calcium PPM"
# plot
ggplot(soil.chem, aes(x=Treatment, y=Ca.mg.kg, color=Treatment)) + 
  geom_boxplot() + theme_bw() + xlab("") + geom_text(data=Ca.Letter2, 
                                                     aes(label=Letter, x=Treatment, y = 800, size=12)) +
  ylab("Calcium (PPM)") + geom_jitter(alpha=0.4) + 
  scale_color_manual(values = c("slategray4", "red3", "goldenrod4", "#003C30")) + 
  theme(legend.position="none") + theme(axis.text.y = element_text(size=12, face="bold")) + 
  theme(axis.text.x = element_text(size = 10, angle=45, hjust=1)) + theme(axis.title.x = element_text(size = 12)) +
  facet_grid(~Ca.facet) +
  theme(strip.text = element_text(size=14, face="bold", color="white")) + 
  theme(strip.background = element_rect(fill="burlywood4"))
# save plot
ggsave(
  filename = "Soil.Ca.PPM.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 4,
  height = 5,
  units = c("in"),
  dpi = 300)

#### Mg.mg.kg / PPM ####

Mg.Letter <- data.frame(multcompLetters(TukeyHSD(aov(Mg.mg.kg ~ Treatment, data = soil.chem))$Treatment[,4])$Letters)
colnames(Mg.Letter)[1] <- "Letter" #Reassign column name
Mg.Letter$Treatment <- rownames(Mg.Letter) #Create column based on rownames
Mg.placement <- soil.chem %>% #We want to create a dataframe to assign the letter position.
  group_by(Treatment) %>%
  summarise(mean=mean(Mg.mg.kg), sd=sd(Mg.mg.kg)) %>%
  arrange(desc(mean))
# merge the letters and placement dataframe for plotting
Mg.Letter2 <- left_join(Mg.Letter, Mg.placement)
soil.chem$Mg.facet <- "Magnesium PPM"
# plot
ggplot(soil.chem, aes(x=Treatment, y=Mg.mg.kg, color=Treatment)) + 
  geom_boxplot() + theme_bw() + xlab("") + geom_text(data=Mg.Letter2, 
                                                     aes(label=Letter, x=Treatment, y = 600, size=12)) +
  ylab("Magnesium (PPM)") + geom_jitter(alpha=0.4) + 
  scale_color_manual(values = c("slategray4", "red3", "goldenrod4", "#003C30")) + 
  theme(legend.position="none") + theme(axis.text.y = element_text(size=12, face="bold")) + 
  theme(axis.text.x = element_text(size = 10, angle=45, hjust=1)) + theme(axis.title.x = element_text(size = 12)) +
  facet_grid(~Mg.facet) +
  theme(strip.text = element_text(size=14, face="bold", color="white")) + 
  theme(strip.background = element_rect(fill="burlywood4"))
# save plot
ggsave(
  filename = "Soil.Mg.PPM.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 4,
  height = 5,
  units = c("in"),
  dpi = 300)

#### K.mg.kg / PPM ####

K.Letter <- data.frame(multcompLetters(TukeyHSD(aov(K.mg.kg ~ Treatment, data = soil.chem))$Treatment[,4])$Letters)
colnames(K.Letter)[1] <- "Letter" #Reassign column name
K.Letter$Treatment <- rownames(K.Letter) #Create column based on rownames
K.placement <- soil.chem %>% #We want to create a dataframe to assign the letter position.
  group_by(Treatment) %>%
  summarise(mean=mean(K.mg.kg), sd=sd(K.mg.kg)) %>%
  arrange(desc(mean))
# merge the letters and placement dataframe for plotting
K.Letter2 <- left_join(K.Letter, K.placement)
soil.chem$K.facet <- "Potassium PPM"
# plot
ggplot(soil.chem, aes(x=Treatment, y=K.mg.kg, color=Treatment)) + 
  geom_boxplot() + theme_bw() + xlab("") + geom_text(data=K.Letter2, 
                                                     aes(label=Letter, x=Treatment, y = 270, size=12)) +
  ylab("Potassium (PPM)") + geom_jitter(alpha=0.4) + 
  scale_color_manual(values = c("slategray4", "red3", "goldenrod4", "#003C30")) + 
  theme(legend.position="none") + theme(axis.text.y = element_text(size=12, face="bold")) + 
  theme(axis.text.x = element_text(size = 10, angle=45, hjust=1)) + theme(axis.title.x = element_text(size = 12)) +
  facet_grid(~K.facet) +
  theme(strip.text = element_text(size=14, face="bold", color="white")) + 
  theme(strip.background = element_rect(fill="burlywood4"))
# save plot
ggsave(
  filename = "Soil.K.PPM.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 4,
  height = 5,
  units = c("in"),
  dpi = 300)

#### Na.mg.kg / PPM ####

Na.Letter <- data.frame(multcompLetters(TukeyHSD(aov(Na.mg.kg ~ Treatment, data = soil.chem))$Treatment[,4])$Letters)
colnames(Na.Letter)[1] <- "Letter" #Reassign column name
Na.Letter$Treatment <- rownames(Na.Letter) #Create column based on rownames
Na.placement <- soil.chem %>% #We want to create a dataframe to assign the letter position.
  group_by(Treatment) %>%
  summarise(mean=mean(Na.mg.kg), sd=sd(Na.mg.kg)) %>%
  arrange(desc(mean))
# merge the letters and placement dataframe for plotting
Na.Letter2 <- left_join(Na.Letter, Na.placement)
soil.chem$Na.facet <- "Sodium PPM"
# plot
ggplot(soil.chem, aes(x=Treatment, y=Na.mg.kg, color=Treatment)) + 
  geom_boxplot() + theme_bw() + xlab("") + geom_text(data=Na.Letter2, 
                                                     aes(label=Letter, x=Treatment, y = 150, size=12)) +
  ylab("Sodium (PPM)") + geom_jitter(alpha=0.4) + 
  scale_color_manual(values = c("slategray4", "red3", "goldenrod4", "#003C30")) + 
  theme(legend.position="none") + theme(axis.text.y = element_text(size=12, face="bold")) + 
  theme(axis.text.x = element_text(size = 10, angle=45, hjust=1)) + theme(axis.title.x = element_text(size = 12)) +
  facet_grid(~Na.facet) +
  theme(strip.text = element_text(size=14, face="bold", color="white")) + 
  theme(strip.background = element_rect(fill="burlywood4"))
# save plot
ggsave(
  filename = "Soil.Na.PPM.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 4,
  height = 5,
  units = c("in"),
  dpi = 300)

#### B.mg.kg / PPM ####

B.Letter <- data.frame(multcompLetters(TukeyHSD(aov(B.mg.kg ~ Treatment, data = soil.chem))$Treatment[,4])$Letters)
colnames(B.Letter)[1] <- "Letter" #Reassign column name
B.Letter$Treatment <- rownames(B.Letter) #Create column based on rownames
B.placement <- soil.chem %>% #We want to create a dataframe to assign the letter position.
  group_by(Treatment) %>%
  summarise(mean=mean(B.mg.kg), sd=sd(B.mg.kg)) %>%
  arrange(desc(mean))
# merge the letters and placement dataframe for plotting
B.Letter2 <- left_join(B.Letter, B.placement)
soil.chem$B.facet <- "Boron PPM"
# plot
ggplot(soil.chem, aes(x=Treatment, y=B.mg.kg, color=Treatment)) + 
  geom_boxplot() + theme_bw() + xlab("") + geom_text(data=B.Letter2, 
                                                     aes(label=Letter, x=Treatment, y = 1, size=12)) +
  ylab("Boron (PPM)") + geom_jitter(alpha=0.4) + 
  scale_color_manual(values = c("slategray4", "red3", "goldenrod4", "#003C30")) + 
  theme(legend.position="none") + theme(axis.text.y = element_text(size=12, face="bold")) + 
  theme(axis.text.x = element_text(size = 10, angle=45, hjust=1)) + theme(axis.title.x = element_text(size = 12)) +
  facet_grid(~B.facet) +
  theme(strip.text = element_text(size=14, face="bold", color="white")) + 
  theme(strip.background = element_rect(fill="burlywood4"))
# save plot
ggsave(
  filename = "Soil.B.PPM.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 4,
  height = 5,
  units = c("in"),
  dpi = 300)

#### Fe.mg.kg / PPM ####

Fe.Letter <- data.frame(multcompLetters(TukeyHSD(aov(Fe.mg.kg ~ Treatment, data = soil.chem))$Treatment[,4])$Letters)
colnames(Fe.Letter)[1] <- "Letter" #Reassign column name
Fe.Letter$Treatment <- rownames(Fe.Letter) #Create column based on rownames
Fe.placement <- soil.chem %>% #We want to create a dataframe to assign the letter position.
  group_by(Treatment) %>%
  summarise(mean=mean(Ca.mg.kg), sd=sd(Fe.mg.kg)) %>%
  arrange(desc(mean))
# merge the letters and placement dataframe for plotting
Fe.Letter2 <- left_join(Fe.Letter, Fe.placement)
soil.chem$Fe.facet <- "Iron PPM"
# plot
ggplot(soil.chem, aes(x=Treatment, y=Fe.mg.kg, color=Treatment)) + 
  geom_boxplot() + theme_bw() + xlab("") + geom_text(data=Fe.Letter2, 
                                                     aes(label=Letter, x=Treatment, y = 270, size=12)) +
  ylab("Iron (PPM)") + geom_jitter(alpha=0.4) + 
  scale_color_manual(values = c("slategray4", "red3", "goldenrod4", "#003C30")) + 
  theme(legend.position="none") + theme(axis.text.y = element_text(size=12, face="bold")) + 
  theme(axis.text.x = element_text(size = 10, angle=45, hjust=1)) + theme(axis.title.x = element_text(size = 12)) +
  facet_grid(~Fe.facet) +
  theme(strip.text = element_text(size=14, face="bold", color="white")) + 
  theme(strip.background = element_rect(fill="burlywood4"))
# save plot
ggsave(
  filename = "Soil.Fe.PPM.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 4,
  height = 5,
  units = c("in"),
  dpi = 300)

#### Mn.mg.kg / PPM ####

Mn.Letter <- data.frame(multcompLetters(TukeyHSD(aov(Mn.mg.kg ~ Treatment, data = soil.chem))$Treatment[,4])$Letters)
colnames(Mn.Letter)[1] <- "Letter" #Reassign column name
Mn.Letter$Treatment <- rownames(Mn.Letter) #Create column based on rownames
Mn.placement <- soil.chem %>% #We want to create a dataframe to assign the letter position.
  group_by(Treatment) %>%
  summarise(mean=mean(Mn.mg.kg), sd=sd(Mn.mg.kg)) %>%
  arrange(desc(mean))
# merge the letters and placement dataframe for plotting
Mn.Letter2 <- left_join(Mn.Letter, Mn.placement)
soil.chem$Mn.facet <- "Manganese PPM"
# plot
ggplot(soil.chem, aes(x=Treatment, y=Mn.mg.kg, color=Treatment)) + 
  geom_boxplot() + theme_bw() + xlab("") + geom_text(data=Mn.Letter2, 
                                                     aes(label=Letter, x=Treatment, y = 40, size=12)) +
  ylab("Manganese (PPM)") + geom_jitter(alpha=0.4) + 
  scale_color_manual(values = c("slategray4", "red3", "goldenrod4", "#003C30")) + 
  theme(legend.position="none") + theme(axis.text.y = element_text(size=12, face="bold")) + 
  theme(axis.text.x = element_text(size = 10, angle=45, hjust=1)) + theme(axis.title.x = element_text(size = 12)) +
  facet_grid(~Mn.facet) +
  theme(strip.text = element_text(size=14, face="bold", color="white")) + 
  theme(strip.background = element_rect(fill="burlywood4"))
# save plot
ggsave(
  filename = "Soil.Mn.PPM.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 4,
  height = 5,
  units = c("in"),
  dpi = 300)

#### Cu.mg.kg / PPM ####

Cu.Letter <- data.frame(multcompLetters(TukeyHSD(aov(Cu.mg.kg ~ Treatment, data = soil.chem))$Treatment[,4])$Letters)
colnames(Cu.Letter)[1] <- "Letter" #Reassign column name
Cu.Letter$Treatment <- rownames(Cu.Letter) #Create column based on rownames
Cu.placement <- soil.chem %>% #We want to create a dataframe to assign the letter position.
  group_by(Treatment) %>%
  summarise(mean=mean(Cu.mg.kg), sd=sd(Cu.mg.kg)) %>%
  arrange(desc(mean))
# merge the letters and placement dataframe for plotting
Cu.Letter2 <- left_join(Cu.Letter, Cu.placement)
soil.chem$Cu.facet <- "Copper PPM"
# plot
ggplot(soil.chem, aes(x=Treatment, y=Cu.mg.kg, color=Treatment)) + 
  geom_boxplot() + theme_bw() + xlab("") + geom_text(data=Cu.Letter2, 
                                                     aes(label=Letter, x=Treatment, y = 10, size=12)) +
  ylab("Copper (PPM)") + geom_jitter(alpha=0.4) + 
  scale_color_manual(values = c("slategray4", "red3", "goldenrod4", "#003C30")) + 
  theme(legend.position="none") + theme(axis.text.y = element_text(size=12, face="bold")) + 
  theme(axis.text.x = element_text(size = 10, angle=45, hjust=1)) + theme(axis.title.x = element_text(size = 12)) +
  facet_grid(~Cu.facet) +
  theme(strip.text = element_text(size=14, face="bold", color="white")) + 
  theme(strip.background = element_rect(fill="burlywood4"))
# save plot
ggsave(
  filename = "Soil.Cu.PPM.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 4,
  height = 5,
  units = c("in"),
  dpi = 300)

#### Zn.mg.kg / PPM #### (There's a < value, which throws the analysis.)

Zn.Letter <- data.frame(multcompLetters(TukeyHSD(aov(Zn.mg.kg ~ Treatment, data = soil.chem))$Treatment[,4])$Letters)
colnames(Zn.Letter)[1] <- "Letter" #Reassign column name
Zn.Letter$Treatment <- rownames(Zn.Letter) #Create column based on rownames
Zn.placement <- soil.chem %>% #We want to create a dataframe to assign the letter position.
  group_by(Treatment) %>%
  summarise(mean=mean(Zn.mg.kg), sd=sd(Zn.mg.kg)) %>%
  arrange(desc(mean))
# merge the letters and placement dataframe for plotting
Zn.Letter2 <- left_join(Zn.Letter, Zn.placement)
soil.chem$Zn.facet <- "Zinc PPM"
# plot
ggplot(soil.chem, aes(x=Treatment, y=Zn.mg.kg, color=Treatment)) + 
  geom_boxplot() + theme_bw() + xlab("") + geom_text(data=Zn.Letter2, 
                                                     aes(label=Letter, x=Treatment, y = 800, size=12)) +
  ylab("Zinc (PPM)") + geom_jitter(alpha=0.4) + 
  scale_color_manual(values = c("slategray4", "red3", "goldenrod4", "#003C30")) + 
  theme(legend.position="none") + theme(axis.text.y = element_text(size=12, face="bold")) + 
  theme(axis.text.x = element_text(size = 10, angle=45, hjust=1)) + theme(axis.title.x = element_text(size = 12)) +
  facet_grid(~Zn.facet) +
  theme(strip.text = element_text(size=14, face="bold", color="white")) + 
  theme(strip.background = element_rect(fill="burlywood4"))
# save plot
ggsave(
  filename = "Soil.Zn.PPM.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 4,
  height = 5,
  units = c("in"),
  dpi = 300)

#### Al.mg.kg / PPM ####

Al.Letter <- data.frame(multcompLetters(TukeyHSD(aov(Al.mg.kg ~ Treatment, data = soil.chem))$Treatment[,4])$Letters)
colnames(Al.Letter)[1] <- "Letter" #Reassign column name
Al.Letter$Treatment <- rownames(Al.Letter) #Create column based on rownames
Al.placement <- soil.chem %>% #We want to create a dataframe to assign the letter position.
  group_by(Treatment) %>%
  summarise(mean=mean(Al.mg.kg), sd=sd(Al.mg.kg)) %>%
  arrange(desc(mean))
# merge the letters and placement dataframe for plotting
Al.Letter2 <- left_join(Al.Letter, Al.placement)
soil.chem$Al.facet <- "Aluminum PPM"
# plot
ggplot(soil.chem, aes(x=Treatment, y=Al.mg.kg, color=Treatment)) + 
  geom_boxplot() + theme_bw() + xlab("") + geom_text(data=Al.Letter2, 
                                                     aes(label=Letter, x=Treatment, y = 600, size=12)) +
  ylab("Aluminum (PPM)") + geom_jitter(alpha=0.4) + 
  scale_color_manual(values = c("slategray4", "red3", "goldenrod4", "#003C30")) + 
  theme(legend.position="none") + theme(axis.text.y = element_text(size=12, face="bold")) + 
  theme(axis.text.x = element_text(size = 10, angle=45, hjust=1)) + theme(axis.title.x = element_text(size = 12)) +
  facet_grid(~Al.facet) +
  theme(strip.text = element_text(size=14, face="bold", color="white")) + 
  theme(strip.background = element_rect(fill="burlywood4"))
# save plot
ggsave(
  filename = "Soil.Al.PPM.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 4,
  height = 5,
  units = c("in"),
  dpi = 300)

#### NH4-N PPM #### (There's a < that throws the analysis.)

NH4.Letter <- data.frame(multcompLetters(TukeyHSD(aov(NH4.ppm ~ Treatment, data = soil.chem))$Treatment[,4])$Letters)
colnames(NH4.Letter)[1] <- "Letter" #Reassign column name
NH4.Letter$Treatment <- rownames(NH4.Letter) #Create column based on rownames
NH4.placement <- soil.chem %>% #We want to create a dataframe to assign the letter position.
  group_by(Treatment) %>%
  summarise(mean=mean(NH4.ppm), sd=sd(NH4.ppm)) %>%
  arrange(desc(mean))
# merge the letters and placement dataframe for plotting
NH4.Letter2 <- left_join(NH4.Letter, NH4.placement)
soil.chem$NH4.ppm.facet <- "Ammonium PPM"
# plot
ggplot(soil.chem, aes(x=Treatment, y=NH4.ppm, color=Treatment)) + 
  geom_boxplot() + theme_bw() + xlab("") + geom_text(data=NH4.Letter2, 
                                                     aes(label=Letter, x=Treatment, y = 800, size=12)) +
  ylab("Amonium (PPM)") + geom_jitter(alpha=0.4) + 
  scale_color_manual(values = c("slategray4", "red3", "goldenrod4", "#003C30")) + 
  theme(legend.position="none") + theme(axis.text.y = element_text(size=12, face="bold")) + 
  theme(axis.text.x = element_text(size = 10, angle=45, hjust=1)) + theme(axis.title.x = element_text(size = 12)) +
  facet_grid(~NH4.ppm.facet) +
  theme(strip.text = element_text(size=14, face="bold", color="white")) + 
  theme(strip.background = element_rect(fill="burlywood4"))
# save plot
ggsave(
  filename = "Soil.NH4.PPM.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 4,
  height = 5,
  units = c("in"),
  dpi = 300)

#### Bray.IP.ppm ####

Bray.Letter <- data.frame(multcompLetters(TukeyHSD(aov(Bray.IP.ppm ~ Treatment, data = soil.chem))$Treatment[,4])$Letters)
colnames(Bray.Letter)[1] <- "Letter" #Reassign column name
Bray.Letter$Treatment <- rownames(Bray.Letter) #Create column based on rownames
Bray.placement <- soil.chem %>% #We want to create a dataframe to assign the letter position.
  group_by(Treatment) %>%
  summarise(mean=mean(Bray.IP.ppm), sd=sd(Bray.IP.ppm)) %>%
  arrange(desc(mean))
# merge the letters and placement dataframe for plotting
Bray.Letter2 <- left_join(Bray.Letter, Bray.placement)
soil.chem$Bray.facet <- "Bray IP PPM"
# plot
ggplot(soil.chem, aes(x=Treatment, y=Bray.IP.ppm, color=Treatment)) + 
  geom_boxplot() + theme_bw() + xlab("") + geom_text(data=Bray.Letter2, 
                                                     aes(label=Letter, x=Treatment, y = 10, size=12)) +
  ylab("Bray IP (PPM)") + geom_jitter(alpha=0.4) + 
  scale_color_manual(values = c("slategray4", "red3", "goldenrod4", "#003C30")) + 
  theme(legend.position="none") + theme(axis.text.y = element_text(size=12, face="bold")) + 
  theme(axis.text.x = element_text(size = 10, angle=45, hjust=1)) + theme(axis.title.x = element_text(size = 12)) +
  facet_grid(~Bray.facet) +
  theme(strip.text = element_text(size=14, face="bold", color="white")) + 
  theme(strip.background = element_rect(fill="burlywood4"))
# save plot
ggsave(
  filename = "Soil.BrayIP.PPM.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 4,
  height = 5,
  units = c("in"),
  dpi = 300)


#### SPIDER PLOT ####

# install.packages("fmsb")  # once
library(fmsb)
library(dplyr)
library(tidyr)
library(stringr)

# --- 0) Read (tab-delimited) ---
# df_wide <- read.delim("soil_chemistry.tsv", check.names = FALSE)  # keep original column labels

# If your data frame is already in memory, ensure the column names are not altered:
# (We will reference them exactly as they appear in your snippet.)

# --- 1) Choose the 16 variables for the radar axes (edit if you prefer others) ---
vars16 <- c(
  "pH", "Perc.OM", "Est.N.Release",
  "S.PPM", "P.mg.kg", "Ca.mg.kg", "Mg.mg.kg", "K.mg.kg", "Na.mg.kg",
  "Fe.mg.kg", "Cu.mg.kg", "B.mg.kg", "TEC",
  "Al.mg.kg", "NH4.ppm", "Mn.mg.kg", "Bray.IP.ppm"
)

# --- 2) Clean numeric columns (handle values like '< 0.5') ---
to_num <- function(x) {
  if (is.numeric(x)) return(x)
  x <- trimws(as.character(x))
  is_lt <- grepl("^<\\s*([0-9.]+)$", x)
  # rule: treat below-LOD as half the LOD (common practice)
  out <- suppressWarnings(as.numeric(x))
  if (any(is_lt)) {
    lod <- as.numeric(sub("^<\\s*([0-9.]+)$", "\\1", x[is_lt]))
    out[is_lt] <- lod / 2
  }
  out
}

# Apply conversion to the 16 variables (keep all other columns intact)
df_wide <- soil.chem %>%
  mutate(across(all_of(vars16), to_num))


#### SPIDER PLOT 3.0 ####

# ====== PACKAGES ======
library(dplyr)
library(fmsb)
library(stats)

# ====== 0) CLEAN NAMES (only once per session) ======
# If your df was read with default check.names=TRUE you can skip this.
# If not, this makes names syntactic: "P.mg/kg" -> "P.mg.kg"
df_wide <- df_wide %>% rename_with(make.names)

# ====== 1) CHOOSE EXACTLY 16 VARIABLES (RAW COLUMN NAMES IN df_wide) ======
vars_raw <- c(
  "pH", "Perc.OM", "Est.N.Release",
  "S.PPM", "Bray.IP.ppm",
  "Mn.mg.kg", "NH4.ppm",
  "TEC", "P.mg.kg", "Ca.mg.kg", "Mg.mg.kg",
  "K.mg.kg", "Na.mg.kg", "Fe.mg.kg", "Cu.mg.kg", "B.mg.kg", "Al.mg.kg"
)
stopifnot(length(vars_raw) == 17, all(vars_raw %in% names(df_wide)))

# Pretty labels in the SAME order (raw -> pretty)
pretty_map <- setNames(c(
  "Soil pH", "% Organic Matter", "Est. N Release",
  "Sulfur", "Bray P",
  "Manganese", "Ammonium",
  "T.E.C", "Phosphorus", "Calcium", "Magnesium",
  "Potassium", "Sodium", "Iron", "Copper", "Boron", "Aluminum"
), vars_raw)

# ====== 2) SUMMARIZE -> ONE ROW PER TREATMENT, THEN BUILD fmsb TABLE ======
wide_med <- df_wide %>%
  group_by(Treatment) %>%
  summarize(across(all_of(vars_raw), ~ median(.x, na.rm = TRUE)), .groups = "drop")

# Axis ranges from all samples (so rings reflect observed spread)
rng <- df_wide %>% summarize(across(all_of(vars_raw), ~ range(.x, na.rm = TRUE)))
axis_max <- as.data.frame(rng[2, vars_raw, drop = FALSE])
axis_min <- as.data.frame(rng[1, vars_raw, drop = FALSE])

# fmsb input: max, min, then one row per treatment (keep order!)
radar_dat <- bind_rows(axis_max, axis_min, wide_med %>% select(all_of(vars_raw)))
rownames(radar_dat) <- c("max", "min", wide_med$Treatment)
# rename columns to pretty labels AFTER ordering
colnames(radar_dat) <- pretty_map[vars_raw]

# ====== 3) STATS: ANOVA -> BH/FDR -> STARS (same order as plotted labels) ======
pvals <- sapply(vars_raw, function(v) {
  summary(aov(reformulate("Treatment", response = v), data = df_wide))[[1]][["Pr(>F)"]][1]
})
padj <- setNames(p.adjust(pvals, method = "BH"), vars_raw)

# Map to pretty names in the plotted order
padj_plot <- setNames(padj[vars_raw], pretty_map[vars_raw])
stars_plot <- cut(padj_plot,
                  breaks = c(-Inf, .001, .01, .05, .10, Inf),
                  labels = c("***","**","*","·",""))

# Append stars to labels (so fmsb draws them)
labs_with_stars <- ifelse(stars_plot == "", colnames(radar_dat),
                          paste0(colnames(radar_dat), " ", stars_plot))

radar_dat_lab <- radar_dat
colnames(radar_dat_lab) <- labs_with_stars

# ====== 4) PLOT + LEGEND (far right) ======
earthy_cols <- c("slategray4", "red3", "goldenrod4", "#003C30")  # edit if >4 treatments
n_tr <- nrow(radar_dat_lab) - 2

png(filename="radar_plot.png", width=4300, height=2900, res=300, bg="white")


op <- par(xpd = NA, mar = c(1, 2, 2, 12))  # extra right margin so legend never overlaps
radarchart(
  radar_dat_lab,
  axistype = 1,
  caxislabels = rep("", 5),
  pcol  = earthy_cols[seq_len(n_tr)],
  pfcol = sapply(earthy_cols[seq_len(n_tr)], \(z) adjustcolor(z, 0.22)),
  plwd  = 2, plty = 1,
  cglcol = "grey80", cglty = 1, cglwd = 0.8,
  seg = 4,
  vlcex = 1.0  # fmsb draws the (starred) labels
)

legend(
  x = grconvertX(1.03, from = "ndc", to = "user"),   # flush to the right
  y = grconvertY(0.5,  from = "ndc", to = "user"),   # halfway up vertically
  legend = rownames(radar_dat_lab)[-c(1, 2)],
  col = earthy_cols[seq_len(nrow(radar_dat_lab)-2)],
  lwd = 2, bty = "n", cex = 0.9,
  yjust = 0.5   # center legend vertically on that y
)

mtext("* FDR<0.05; ** FDR<0.01; *** FDR<0.001; · FDR<0.10", side = 1, line = 0, adj = 0)
par(op)
# after radarchart(...)
k   <- ncol(radar_dat_lab)
ang <- seq(0, 2*pi, length.out = k+1)[1:k]

# extract min/max per variable
vmin <- as.numeric(radar_dat_lab[2, ])
vmax <- as.numeric(radar_dat_lab[1, ])

# Example: just show max value near outer edge
x1 <- 1.05 * sin(ang);  y1 <- 1.05 * cos(ang)
text(x1, y1, labels = format(vmax, digits = 3, trim = TRUE),
     cex = 0.7, col = "grey30")
dev.off()




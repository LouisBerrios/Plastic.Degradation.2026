### 16S Amplicon pipeline for plastic degradation experiment ###

# DADA2 workflow for processing 16S raw sequences
# Follows https://benjjneb.github.io/dada2/tutorial.html

library(dada2)
library(ShortRead)
library(Biostrings)


#############################

path <- "~/Documents/Plastic_Expt/16S_Amps"
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq.gz and SAMPLENAME_R2_001.fastq.gz
fnFs <- sort(list.files(path, pattern="_R1.fastq.gz", full.names = TRUE)) # If set to be FALSE, then working directory must contain the files
fnRs <- sort(list.files(path, pattern="_R2.fastq.gz", full.names = TRUE))

# Remove any forward files that don't have reverse counterparts, and vise versa
# (filterAndTrim will throw an error if fnFs and fnRs have any mismatches)
basefilenames_Fs <- sub("_R1.fastq.gz","",basename(fnFs))
basefilenames_Rs <- sub("_R2.fastq.gz","",basename(fnRs))
rm_from_fnFs <- basefilenames_Fs[which(!(basefilenames_Fs %in% basefilenames_Rs))]
rm_from_fnRs <- basefilenames_Rs[which(!(basefilenames_Rs %in% basefilenames_Fs))]

for(name in rm_from_fnFs) {
  print(paste(name, "does not have a reverse-reads counterpart. Omitting from this analysis."))
  fnFs <- fnFs[-which(fnFs == paste0(path, "/", name, "_R1.fastq.gz"))]
}
for(name in rm_from_fnRs) {
  print(paste(name, "does not have a forward-reads counterpart. Omitting from this analysis."))
  fnRs <- fnRs[-which(fnRs == paste0(path, "/", name, "_R2.fastq.gz"))]
}


# Identify primers - used 515F & 806R from Hiro's spreadsheet
FWD <- "GTGYCAGCMGCCGCGGTAA"  ## CHANGE ME to your forward primer seq
REV <- "GGACTACNVGGGTWTCTAAT"  ## CHANGE ME...


# Get all orientations of primers, just to be safe
allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients

# “pre-filter” the sequences just to remove those with Ns, but perform no other filtering
fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 1, multithread = TRUE)

# (From tutorial) We are now ready to count the number of times the primers appear in the 
# forward and reverse read, while considering all possible primer orientations. 
# Identifying and counting the primers on one set of paired end FASTQ files is
# sufficient, assuming all the files were created using the same library preparation,
# so we’ll just process the first sample.

primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))

# If you see the reverse-complement of the forward primer in the reverse reads (cells [2,4] and [3,4]),
# it's because the ITS region is short and it is reading part of the forward primer.

# Remove primers using cutadapt

cutadapt <- "/Users/louisberrios/Documents/Cutadapt/cutadapt" # CHANGE ME to the cutadapt path on your machine
system2(cutadapt, args = "--version") # Run shell commands from R

path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs.filtN))
fnRs.cut <- file.path(path.cut, basename(fnRs.filtN))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 
# Run Cutadapt
for(i in seq_along(fnFs.filtN)) {
  # for(i in 1:10) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i], # input files; fnFs.filtN replaced by fnFs.filtN, etc.
                             "--minimum-length", "1")) # min length of cutadapted reads: >0 
}

# Count primers in first post-cutadapt sample (should all be 0):
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[20]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[20]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[20]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[20]]))

# Since they are zero, skip step to remove other orientations of primers
## These were 1,2,2 but not far from zero so will move on


# Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path.cut, pattern = "_R1.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "_R2.fastq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have format:
get.sample.name <- function(fname) {
  paste(strsplit(basename(fname), split="_")[[1]][1:3], collapse="_")
}
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)

# Inspect read quality profiles of forward reads #1-2
plotQualityProfile(cutFs[1:12])

# Inspect read quality profiles of reverse reads #1-2
plotQualityProfile(cutRs[1:12])

# Filter and trim

# Assigning the filenames for the output of the filtered reads 
# to be stored as fastq.gz files.
filtFs <- file.path(path, "filtered", basename(fnFs.filtN))
filtRs <- file.path(path, "filtered", basename(fnRs.filtN))

out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, maxEE = c(2, 2), truncQ = 2, minLen = 100, rm.phix = TRUE, compress = TRUE, multithread = TRUE)  # on windows, set multithread = FALSE
head(out)

filtFs.out <- list.files(paste(path, "filtered", sep="/"), pattern="_R1.fastq.gz", full.names=TRUE)
filtRs.out <- list.files(paste(path, "filtered", sep="/"), pattern="_R2.fastq.gz", full.names=TRUE)

# Learn the error rates
errF <- learnErrors(filtFs.out, multithread = TRUE)
errR <- learnErrors(filtRs.out, multithread = TRUE)

# Visualize estimated error rates
plotErrors(errF, nominalQ = TRUE)

# NEED TO INCREASE MEMORY ALLOCATED TO R <-  do this in the terminal
#Step 1: Open terminal, #Step 2:
# cd ~
# touch .Renviron
# open .Renviron
#Step 3: Save the following as the first line of .Renviron:
#  R_MAX_VSIZE=200Gb 
#Step 4: Close RStudio and reopen

# Dereplicate identical reads
derepFs <- derepFastq(filtFs.out, verbose = TRUE)
derepRs <- derepFastq(filtRs.out, verbose = TRUE)
# Name the derep-class objects by the sample names
get.sample.name <- function(fname) {
  paste(strsplit(basename(fname), "_")[[1]][1:3], collapse="_")
}
sample.names <- unname(sapply(filtFs.out, get.sample.name))

# DADA2's core sample inference algorithm
dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
dadaRs <- dada(derepRs, err = errR, multithread = TRUE)


# Merge pairs
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE,trimOverhang = TRUE)

# Construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

saveRDS(seqtab,"Plastic.Deg_16S_seqtab.bac.rds")
saveRDS(seqtab.nochim,"Plasttic.Deg_16S_seqtab.bac.nochim.rds")
seqtab.nochim <- readRDS("Plasttic.Deg_16S_seqtab.bac.nochim.rds")
# Inspect distribution of sequence lengths
hist(nchar(getSequences(seqtab.nochim)))

# Track reads through pipeline
getN <- function(x) sum(getUniques(x))

#format out to accommodate dropped samples
raw.sample.names <- unname(sapply(row.names(out), get.sample.name))

track <- cbind(sapply(dadaFs, getN), sapply(dadaRs, getN), 
               sapply(mergers, getN), rowSums(seqtab.nochim))

track2<-cbind(out,track[match(row.names(out),row.names(track)),])

# If processing a single sample, remove the sapply calls: e.g. replace
# sapply(dadaFs, getN) with getN(dadaFs)
colnames(track2) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", 
                      "nonchim")
write.csv(track2,"Plastic.Deg_16S_summary.csv")

rownames(track) <- sample.names
head(track2)


# Assign taxonomy using the UNITE database
silva.ref<-"silva_nr99_v138.1_train_set.fa.gz"
silva.species<-"silva_species_assignment_v138.1.fa.gz"

taxa <- assignTaxonomy(seqtab.nochim, silva.ref, multithread = TRUE, tryRC = TRUE)

taxa <- addSpecies(taxa, silva.species)

taxa.print <- taxa  # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)


# Save OTU table and taxonomic table as RDS files
# to hand off to dada2_to_phyloseq.R
saveRDS(seqtab.nochim, "Plastic.Deg_16S_seqtab_nochim_taxa.rds")
saveRDS(taxa, "Plastic.Deg_16S_taxa.rds")
seqtab.nochim <- readRDS("Plastic.Deg_16S_seqtab_nochim_taxa.rds")
########Construct phyloseq object with sample data variables#########

# make sample table
library(phyloseq)
sample <- as.data.frame(rownames(seqtab.nochim))
colnames(sample)[1] <- "Sample"
library(dplyr)
library(tibble) 
sample2 <- sample %>%
  mutate(Compartment = ifelse(grepl("PL", Sample), "Plastisphere", "Bulk"))
sample2$Treatment <- ifelse(grepl("FNoP", sample2$Sample), "Fungus_NoPlastic",
                       ifelse(grepl("PF", sample2$Sample), "Plastic+Fungus",
                              ifelse(grepl("NoPNoF", sample2$Sample), "NoPlastic_NoFungus",
                                     ifelse(grepl("PNoF", sample2$Sample), "Plastic_NoFungus", NA))))
# Load and merge soil chemistry data
soil.chem <- read.csv("16S.SoilChem.csv")
SD <- right_join(sample2,soil.chem, by="Sample") # this is only bulk soil since they are the only ones with soil chem
write.csv(SD, "16S.Sample.Data.All.csv")

#read in environmental sample table (this is used when building sample data outside of R)

# PD.16S <- readRDS("Plastic.Deg_16S_seqtab_nochim_taxa.rds")
colnames(sample2)[1] <- "Sample"
rownames(sample2) <- sample2$Sample # this code is important for joining sample data with taxa and seq data
SAM <- sample_data(sample2, errorIfNULL = T) # use this for the alpha diversity sample data matching
# PD.16S.taxa <- readRDS("Plastic.Deg_16S_taxa.rds")

#create a phyloseq object
bac.ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), sample_data(SAM), tax_table(taxa))

#filter out unwanted taxa (e.g., mitochondira, chloroplast, and archaea sequences)
bac.ps.filt<-subset_taxa(bac.ps,Family!="Mitochondria")
bac.ps.filt<-subset_taxa(bac.ps.filt,Genus!="Chloroplast")
bac.ps.filt <- subset_taxa(bac.ps.filt,Kingdom!="Archaea")

# Removing sequence rownames for display only
taxa.print <- tax_table(bac.ps.filt)
rownames(taxa.print) <- NULL
head(taxa.print)

#save the filtered dataset 
saveRDS(bac.ps.filt,"PlasticDeg.16S.filtered.rds")

#filter out low abundant sequences
bac.ps.filt2 = prune_taxa(taxa_sums(bac.ps.filt) > 10, bac.ps.filt) 
bac.ps.filt2 = prune_samples(sample_sums(bac.ps.filt2)>1000, bac.ps.filt2)

#save the filtered+pruned dataset
saveRDS(bac.ps.filt2, "PlasticDeg.16S.filt-prune.rds")

bac.ps.filt2 <- readRDS("PlasticDeg.16S.filt-prune.rds")

#rarefy the dataset
bac.ps.rare <- rarefy_even_depth(bac.ps.filt2, sample.size = min(sample_sums(bac.ps.filt2)),
                                 rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE) #removed 1,617 OTUs
# remove controls
#rare2 <- subset_samples(bac.ps.rare, Condition != "NA") # use this file for the DESeq2 analysis
#raredf <- psmelt(rare2)
#save rarefied dataset
saveRDS(bac.ps.rare, "PlasticDeg.RARE.Bacteria.rds")
bac.ps.rare <- readRDS("PlasticDeg.RARE.Bacteria.rds")
#create ASV-TAX table
taxa_names(bac.ps.rare) <- paste0("Seq", seq(ntaxa(bac.ps.rare)))
ASV_BAC <- data.frame(otu_table(bac.ps.rare))
TAX_ASV <- data.frame(tax_table(bac.ps.rare))
ASV_BAC_T <- t(ASV_BAC)
MERGER2 <- merge(ASV_BAC_T, TAX_ASV, by = "row.names")
write.csv(MERGER2, "ASV_Bacteria_PD_Rare_Final.csv")


#transform data [compositional]
pseq <- microbiome::transform(bac.ps.rare, "compositional")

#save as RDS
saveRDS(pseq, "PlasticDeg-Bacteria-PS_filter+rare+trans.rds")

# create data frame 
pseq.df <- psmelt(pseq)
earthy_colors <- c(
  "#A6611A", # Brownish orange
  "#DFC27D", # Tan
  "slategray4", # Muted teal
  "#018571", # Deep green
  "#D73027", # Earthy red
  "#FC8D59", # Salmon
  "#91BFDB", # Sky blue
  "#4575B4", # Deep blue
  "#543005", # Dark brown
  "#8C510A", # Rust brown
  "#01665E", # Dark teal
  "#003C30", # Forest green
  "#F6E8C3", # Light sand
  "#C7EAE5", # Light teal
  "#F46D43", # Orange-red
  "gray", # Crimson
  "#762A83", # Purple
  "#998EC3"  # Lavender
)

# remove controls from data frame
pseq.df2 <- subset(pseq.df, Treatment !="NA")
# subset bulk soil and plastisphere
bulk.df <- subset(pseq.df2, Compartment == "Bulk")
plast.df <- subset(pseq.df2, Compartment == "Plastisphere")
# add facet label for aesthetic 
bulk.df$label <- "Bulk soil"
# change Treatment names
custom_labels <- c(
  "Fungus_NoPlastic" = "Fungus (No Plastic)",
  "NoPlastic_NoFungus" = "No Fungus / No Plastic",
  "Plastic_NoFungus" = "Plastic (No Fungus)",
  "Plastic+Fungus" = "Plastic + Fungus"
)
custom_labels2 <- c("Bulk" = "Bulk soil", "Plastisphere" = "Plastisphere")
#plot the relative abundance by Phylum
library(ggplot2)
ggplot(bulk.df, aes(fill=Phylum, y=Abundance, x=Treatment)) + 
  geom_bar(position="fill", stat="identity") + theme_bw() + 
  ylab("Relative Abundance") +  
  scale_fill_manual(values = earthy_colors) +
  scale_y_continuous(labels = scales::percent) + xlab("") + 
  theme(axis.text.x = element_text(size=12, face="bold", angle=45, hjust = 1)) + 
  theme(axis.text.y = element_text(face = "bold", size = 12)) + 
  theme(axis.title.y = element_text(size = 16, face="bold")) + facet_grid(~label) +
  theme(strip.text = element_text(size=14, face="bold")) + 
  theme(strip.background = element_rect(fill = "burlywood4")) +
  scale_x_discrete(labels=custom_labels)

# save
ggsave(
  filename = "PD.BulkSoil.Bacteria.CompBarPlot.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 5,
  height = 6.5,
  units = c("in"),
  dpi = 300)

# plastisphere comp plot
ggplot(plast.df, aes(fill=Phylum, y=Abundance, x=Treatment)) + 
  geom_bar(position="fill", stat="identity") + theme_bw() + 
  ylab("Relative Abundance") +  
  scale_fill_manual(values = earthy_colors) +
  scale_y_continuous(labels = scales::percent) + xlab("") + 
  theme(axis.text.x = element_text(size=12, face="bold", angle=45, hjust = 1)) + 
  theme(axis.text.y = element_text(face = "bold", size = 12)) + 
  theme(axis.title.y = element_text(size = 16, face="bold")) + facet_grid(~Compartment) +
  theme(strip.text = element_text(size=14, face="bold")) + 
  theme(strip.background = element_rect(fill = "azure3")) +
  scale_x_discrete(labels=custom_labels)

# save
ggsave(
  filename = "PD.Plastisphere.Bacteria.CompBarPlot.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 5,
  height = 6.5,
  units = c("in"),
  dpi = 300)


# create dataframe by Phylum
Agg.phy <- aggregate(Abundance ~ Phylum + Condition, data=df_subset, FUN=sum)
# create aggregated by genus dataframe
AGG <- aggregate(Abundance ~ Genus + Condition, data = raredf, FUN=sum)
# subset top 20 most abundant genera for each condition (i.e., EcM + No EcM)
EcM <- subset(AGG, Condition == "EcM")
noEcM <- subset(AGG, Condition == "No-EcM")
top20EcM <- subset(EcM, Abundance > 450) # top 20 
topnoECM <- subset(noEcM, Abundance > 480) # top 20
# merge dataframes 
top20 <- merge(top20EcM, topnoECM, all=TRUE)
# save dataframe 
write.csv(top20, "Merged_Top20_Genus-Condition.csv")
# plot barplot
ggplot(top20, aes(x=Genus, y=Abundance)) + xlab("") + 
  ylab("No. of Sequences") + geom_bar(stat="identity", fill = c("black")) + 
  coord_flip() + theme_bw() + theme(strip.text = element_text(size = 8, face = "bold")) + 
  theme(axis.text.x = element_text(size = 6, face = "bold")) + 
  theme(axis.text.y = element_text(size = 6, face = "bold.italic")) + 
  theme(axis.title.x = element_text(size = 7, face = "bold", vjust = 0.5)) + 
  facet_grid(~Condition)
ggsave(
  filename = "EMSL_16S_Top20_by_ConditionBarPlot.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 5,
  height = 5.5,
  units = c("in"),
  dpi = 300)

####### Alpha Diversity ###########
# rarefied data use non-rareified data
# remove controls from data set
rare.no.con <- subset_samples(bac.ps.rare, Treatment !="NA")
ad <- estimate_richness(rare.no.con)
write.csv(ad, "Plastic.Deg.Bacteria_Alpha_Diversity.FINAL.csv")
ad <- read.csv("Plastic.Deg.Bacteria_Alpha_Diversity.FINAL.csv")

# match AD to sample data metrics
colnames(ad)[1] <- "Sample"
ad.sam <- right_join(ad,SAM, by="Sample")
ad.sam <- subset(ad.sam, Treatment !="NA")# remove controls from the dataset
# SAVE final sample data
write.csv(ad.sam, "PD.16S.FinalSampleData.csv")
ad.sam.bulk <- subset(ad.sam, Compartment == "Bulk")
ad.sam.plast <- subset(ad.sam, Compartment == "Plastisphere")
# model to determine effect
# BULK soil
ad.bulk.model <- aov(Observed ~ Treatment, data = ad.sam.bulk)
summary(ad.bulk.model)
# Plastisphere
ad.plast.model <- aov(Observed ~ Treatment, data = ad.sam.plast)
summary(ad.plast.model)
TukeyHSD(ad.plast.model) # not significant (have to change code from bulk soil comparison)
 # Bulk vs. Plastisphere
ad.sam.all <- aov(Observed ~ Compartment, data = ad.sam)
summary(ad.sam.all) # not significant (repeat steps taken for plastisphere)

# Assign Tukey test letter codes to illustrate significant differences
### install.packages("multcompView") ###
library(multcompView)
library(tidyverse)

### BULK SOIL ONLY ###
PC.letters.df.bulk <- data.frame(multcompLetters(TukeyHSD(aov(Observed ~ Treatment, data = ad.sam.bulk))$Treatment[,4])$Letters)
colnames(PC.letters.df.bulk)[1] <- "Letter" #Reassign column name
PC.letters.df.bulk$Treatment <- rownames(PC.letters.df.bulk) #Create column based on rownames
PC.placement.bulk <- ad.sam.bulk %>% #We want to create a dataframe to assign the letter position.
  group_by(Treatment) %>%
  summarise(mean=mean(Observed), sd=sd(Observed)) %>%
  arrange(desc(mean))
# merge the letters and placement dataframe for plotting
PC.letter.df2.bulk <- left_join(PC.letters.df.bulk, PC.placement.bulk)
ad.sam.bulk$label <- "Bacterial Richness (Bulk Soil)"
# plot for Bulk Soil only
ggplot(ad.sam.bulk, aes(x=Treatment, y=Observed, color=Treatment)) + 
  geom_boxplot() + theme_bw() + xlab("") + geom_text(data=PC.letter.df2.bulk, 
                                                     aes(label=Letter, x=Treatment, y = mean, vjust =-8.5, size=12)) +
  ylab("Observed Richness") + geom_jitter(alpha=0.4) + 
  scale_color_manual(values = c("slategray4", "red3", "goldenrod4", "#003C30")) + 
  theme(legend.position="none") + theme(axis.text.y = element_text(size=12, face="bold")) + 
  theme(axis.text.x = element_text(size = 10, angle=45, hjust=1)) + theme(axis.title.x = element_text(size = 12)) +
  facet_grid(~label) +   scale_x_discrete(labels=custom_labels) +
  theme(strip.text = element_text(size=14, face="bold")) + theme(strip.background = element_rect(fill="burlywood4"))
# save plot
ggsave(
  filename = "PlasticDeg.16S.Richness.BulkSoil.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 4,
  height = 5.5,
  units = c("in"),
  dpi = 300)

### PLASTISPHERE ONLY ### wasn't working with standard TukeyHSD: just made new DF with appropriate letter codes
PC.letters.df.plast <- data.frame(Letter = c("A","A"))
colnames(PC.letters.df.plast)[1] <- "Letter" #Reassign column name
PC.letters.df.plast$Treatment <- c("Plastic+Fungus", "Plastic_NoFungus") #Create column based on rownames
PC.placement.plast <- ad.sam.plast %>% #We want to create a dataframe to assign the letter position.
  group_by(Treatment) %>%
  summarise(mean=mean(Observed), sd=sd(Observed)) %>%
  arrange(desc(mean))
# merge the letters and placement dataframe for plotting
PC.letter.df2.plast <- left_join(PC.letters.df.plast, PC.placement.plast)
ad.sam.plast$label <- "Bacterial Richness (Plastisphere)"
#### Plastisphere Only ####
ggplot(ad.sam.plast, aes(x=Treatment, y=Observed, color=Treatment)) + 
  geom_boxplot() + theme_bw() + xlab("") + geom_text(data=PC.letter.df2.plast, 
                                                     aes(label=Letter, x=Treatment, y = mean,vjust=-12, size=12)) +
  ylab("Observed Richness") + geom_jitter(alpha=0.4) + 
  scale_color_manual(values = c("goldenrod4", "#003C30")) + 
  theme(legend.position="none") + theme(axis.text.y = element_text(size=12, face="bold")) + 
  theme(axis.text.x = element_text(size = 10, angle=45, hjust=1)) + theme(axis.title.x = element_text(size = 12)) +
  facet_grid(~label) + scale_x_discrete(labels=custom_labels) +
  theme(strip.text = element_text(size=14, face="bold")) + 
  theme(strip.background = element_rect(fill="azure3"))

# save plot
ggsave(
  filename = "PlasticDeg.16S.Richness.Plastisphere.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 4,
  height = 5.5,
  units = c("in"),
  dpi = 300)

#### Comparison between bulk soil and plastisphere ####
## Since there was no significant difference between plastisphere communities, 
# let's compare bulk soil and plastisphere richness

PC.letters.df.all <- data.frame(Letter = c("A","A"))
colnames(PC.letters.df.all)[1] <- "Letter" #Reassign column name
PC.letters.df.all$Compartment <- c("Bulk", "Plastisphere") #Create column based on rownames
PC.placement.all <- ad.sam %>% #We want to create a dataframe to assign the letter position.
  group_by(Compartment) %>%
  summarise(mean=mean(Observed), sd=sd(Observed)) %>%
  arrange(desc(mean))
# merge the letters and placement dataframe for plotting
PC.letter.df2.all <- left_join(PC.letters.df.all, PC.placement.all)
ad.sam$label <- "Bacterial Richness (Bulk soil versus Plastisphere)"
#### Bulk v Plastisphere ####
ggplot(ad.sam, aes(x=Compartment, y=Observed, color=Compartment)) + 
  geom_boxplot() + theme_bw() + xlab("") + geom_text(data=PC.letter.df2.all, 
                                                     aes(label=Letter, x=Compartment, y = mean,vjust=-12, size=12)) +
  ylab("Observed Richness") + geom_jitter(alpha=0.4) + 
  scale_color_manual(values = c("burlywood4", "azure4")) + 
  theme(legend.position="none") + theme(axis.text.y = element_text(size=12, face="bold")) + 
  theme(axis.text.x = element_text(size = 10, angle=45, hjust=1)) + theme(axis.title.x = element_text(size = 12)) +
  facet_grid(~label) + scale_x_discrete(labels=custom_labels) +
  theme(strip.text = element_text(size=9.5, face="bold")) + 
  theme(strip.background = element_rect(fill="gray")) + scale_x_discrete(labels=custom_labels2)

# save plot
ggsave(
  filename = "PlasticDeg.16S.Richness.BulkSoilvPlastisphere.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 4,
  height = 5.5,
  units = c("in"),
  dpi = 300)

# Beta Diversity (i.e., Ordinations)

########PERMANOVA###########

#calculate bray-curtis distance matrix
library(vegan)
pseq.no.con <- subset_samples(pseq, Treatment !="NA")
# remove FNoP_3_16S_R1.fastq.gz and FNoP_4_16S_R1.fastq.gz because they are NMDS outliers
pseq.no.con1 <- subset_samples(pseq.no.con, Sample != "FNoP_3_16S_R1.fastq.gz")
pseq.no.con1 <- subset_samples(pseq.no.con1, Sample != "FNoP_4_16S_R1.fastq.gz")

# separate bulk soil and plastisphere samples
pseq.plast <- subset_samples(pseq, Compartment == "Plastisphere")
pseq.plast <- subset_samples(pseq.plast, Treatment !="NA")
pseq.bulk <- subset_samples(pseq.no.con1, Compartment == "Bulk")
pseq.bulk <- subset_samples(pseq.bulk, Treatment !="NA") # remove controls


ord.bac.bray <- phyloseq::distance(pseq.no.con1, method = "bray", k=3)
ord.bac.bray.plast <- phyloseq::distance(pseq.plast, method = "bray", k=3)
ord.bac.bray.bulk <- phyloseq::distance(pseq.bulk, method = "bray", k=3)

#make a dataframe from the sample_data
ord.bac.bray.sampleDF <- data.frame(sample_data(pseq.no.con1)) # no controls
ord.bac.bray.sampleDF.plast <- data.frame(sample_data(pseq.plast)) # only plastisphere
ord.bac.bray.sampleDF.bulk <- data.frame(sample_data(pseq.bulk)) # only bulk soils

#adonis test
PERMANOVA.all <- adonis2(ord.bac.bray ~ Treatment + Compartment, data = ord.bac.bray.sampleDF)
PERMANOVA.plast <- adonis2(ord.bac.bray.plast ~ Treatment, data = ord.bac.bray.sampleDF.plast)
PERMANOVA.bulk <- adonis2(ord.bac.bray.bulk ~ Treatment, data = ord.bac.bray.sampleDF.bulk)

#write PERMANOVA results to table CSV
bac.per.res.all <- as.data.frame(PERMANOVA.all)
write.csv(bac.per.res.all, "PD.AllSamps.Bac.PERMANOVA.Results.csv")
bac.perm.res.plast <- as.data.frame(PERMANOVA.plast)
write.csv(bac.perm.res.plast, "PD.PLASTISPHERE.Bac.PERMANOVA.Results.csv")
bac.perm.res.bulk <- as.data.frame(PERMANOVA.bulk)
write.csv(bac.perm.res.bulk, "PD.BULK.Bac.PERMANOVA.Results.csv")


ord.bac.all <- ordinate(pseq.no.con1, "NMDS", "bray", k=4)
ord.bac.plast <- ordinate(pseq.plast, "NMDS", "bray", k=3)
ord.bac.bulk <- ordinate(pseq.bulk, "NMDS", "bray", k=3)
scores.bulk <- as.data.frame(scores(ord.bac.bulk)$sites)


#plot NMDS ordination [axes 1:2] w/ PERMANOVA results 
# plast and bulk
sample_data(pseq.no.con1)$label1 <- "Bulk soil + Plastisphere"
plot_ordination(pseq.no.con1, ord.bac.all, type="sample", shape="Compartment", color = "Compartment", axes=1:2) + 
  theme_bw()  + scale_color_manual(values=c("burlywood4", "azure3")) +
  geom_point(size = 4, alpha = 1) + stat_ellipse(type ="norm", linetype = 2) +
  facet_grid(~label1) + 
  theme(strip.background = element_rect(fill = "gray")) +
  theme(strip.text = element_text(size=14, face = "bold")) +
annotate(geom = 'text',
         x = -0.6,
         y = 0.6, 
         label = paste("F[\"1,117\"] ==", 27.57), 
         parse = TRUE,
         hjust = 0) +
  annotate(geom = 'text',
           x = -0.6,
           y = 0.55, 
           label = paste("p < 0.001"), 
           parse = TRUE,
           hjust = 0) +
  annotate(geom = 'text',
           x = -0.6,
           y = 0.5, 
           label = paste(
             "R^2",
             "==",
             0.145), 
           parse = TRUE,
           hjust = 0)
  
#save
ggsave(
  filename = "NMDS_Plastisphere.BULK.16S.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 5,
  height = 5,
  units = c("in"),
  dpi = 300)

#plastisphere
plot_ordination(pseq.plast, ord.bac.plast, type="sample", color="Treatment", axes=1:2) + 
  theme_bw() + scale_color_manual(values = c("goldenrod4", "#003C30"), labels=custom_labels) + 
  geom_point(size = 4, alpha = 1) + stat_ellipse(type ="norm", linetype = 2) + 
  facet_grid(~Compartment) + theme(strip.background = element_rect(fill = "azure3")) +
  theme(strip.text = element_text(size=14, face = "bold")) +
  annotate(geom = 'text',
           x = -0.5,
           y = 0.6, 
           label = paste("F[\"1,39\"] ==", 7.52), 
           parse = TRUE,
           hjust = 0) +
  annotate(geom = 'text',
           x = -0.5,
           y = 0.55, 
           label = paste("p < 0.001"), 
           parse = TRUE,
           hjust = 0) +
  annotate(geom = 'text',
           x = -0.5,
           y = 0.5, 
           label = paste(
             "R^2",
             "==",
             0.165), 
           parse = TRUE,
           hjust = 0)
#save
ggsave(
  filename = "NMDS_Plastisphere.16S.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 5,
  height = 5,
  units = c("in"),
  dpi = 300)

#bulk soil
bulk.label <- c("Bulk soil")
names(bulk.label) <- "Bulk"
plot_ordination(pseq.bulk, ord.bac.bulk, type="sample", color="Treatment", axes=1:2) + 
  theme_bw() + scale_color_manual(values = c("slategray4", "red3", "goldenrod4", "#003C30"), labels=custom_labels) + 
  geom_point(size = 4, alpha = 1) + stat_ellipse(type ="norm", linetype = 2) + 
  facet_grid(~Compartment, labeller = labeller(Compartment=bulk.label)) + 
  theme(strip.background = element_rect(fill = "burlywood4")) +
  theme(strip.text = element_text(size=14, face = "bold")) +
  annotate(geom = 'text',
           x = -0.6,
           y = 0.5, 
           label = paste("F[\"1,77\"] ==", 9.39), 
           parse = TRUE,
           hjust = 0) +
  annotate(geom = 'text',
           x = -0.6,
           y = 0.45, 
           label = paste("p < 0.001"), 
           parse = TRUE,
           hjust = 0) +
  annotate(geom = 'text',
           x = -0.6,
           y = 0.4, 
           label = paste(
             "R^2",
             "==",
             0.276), 
           parse = TRUE,
           hjust = 0)
#save
ggsave(
  filename = "NMDS_BulkSoil.16S.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 5,
  height = 5,
  units = c("in"),
  dpi = 300)

##### DESeq2 Analysis #####

library(DESeq2)

# use rarefied (non-transformed dataset) and remove controls
rare.dseq <- subset_samples(bac.ps.rare, Treatment !="NA")
# compare fungal plastisphere to bulk soil first
#convert phyloseq to DESeq2 format
dds <- phyloseq_to_deseq2(rare.dseq, ~ Compartment)

#calculate size factors using edgeR
library(edgeR)
sizeFactors(dds) <- calcNormFactors(counts(dds))

#run DESeq function
dds = DESeq(dds, test="Wald", fitType="parametric")

res <- results(dds, cooksCutoff = FALSE)
alpha = 0.1
sigtab <- res[which(res$padj < alpha), ]
sigtab <- cbind(as(sigtab, "data.frame"), as(tax_table(rare.dseq)[rownames(sigtab), ], "matrix"))

#prepare data to plot DESeq2 output
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))

library(RColorBrewer)
getPalette <- colorRampPalette(brewer.pal(9, "Paired"))
colourCount = length(unique(sigtab$Phylum))

library(ggplot2)
sigtab2 <- subset(sigtab, Genus !="NA")
sigtab2$label <- "Bulk soil versus Plastisphere"
ggplot(sigtab2, aes(x=Genus, y=log2FoldChange, color=Phylum)) + 
  geom_point(size=4, alpha = 0.6) + theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + 
  theme_bw() + coord_flip() + scale_color_manual(values = earthy_colors) + xlab("") + 
  ylab("[Log2] Fold-Change") + theme(axis.text.y = element_text(size = 7)) + 
  theme(plot.title = element_text(size = 7)) + geom_hline(yintercept=0, col="red", linetype = 2) +
  theme(axis.text.y = element_text(face = "bold.italic")) +
  facet_grid(~label) +
  theme(strip.background = element_rect(fill="gray")) +
  theme(strip.text = element_text(size=6, face="bold"))

# save plot 
ggsave(
  filename = "P.Deg_Bacteria_DESEq2.BulkvPlast.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 7.5,
  height = 13,
  units = c("in"),
  dpi = 300)

#### Repeat DESeq2 to test enriched fungi in plastisphere (w/ and w/o fungi added)

rare.dseq.plast <- subset_samples(rare.dseq, Compartment == "Plastisphere")
dds.plast <- phyloseq_to_deseq2(rare.dseq.plast, ~ Treatment)

#calculate size factors using edgeR
library(edgeR)
sizeFactors(dds.plast) <- calcNormFactors(counts(dds.plast))

#run DESeq function
dds.plast = DESeq(dds.plast, test="Wald", fitType="parametric")

res.plast <- results(dds.plast, cooksCutoff = FALSE)
alpha = 0.1
sigtab.plast <- res[which(res.plast$padj < alpha), ]
sigtab.plast <- cbind(as(sigtab.plast, "data.frame"), as(tax_table(rare.dseq.plast)[rownames(sigtab.plast), ], "matrix"))

#prepare data to plot DESeq2 output
x = tapply(sigtab.plast$log2FoldChange, sigtab.plast$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab.plast$Phylum = factor(as.character(sigtab.plast$Phylum), levels=names(x))
x = tapply(sigtab.plast$log2FoldChange, sigtab.plast$Genus, function(x) max(x))
sigtab.plast$Genus = factor(as.character(sigtab.plast$Genus), levels=names(x))


library(ggplot2)
sigtab2.plast <- subset(sigtab.plast, Genus !="NA")
sigtab2.plast$label <- "Plastisphere"
ggplot(sigtab2.plast, aes(x=Genus, y=log2FoldChange, color=Phylum)) + 
  geom_point(size=4, alpha = 0.6) + theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + 
  theme_bw() + coord_flip() + scale_color_manual(values = earthy_colors) + xlab("") + 
  ylab("[Log2] Fold-Change") + theme(axis.text.y = element_text(size = 7)) + 
  theme(plot.title = element_text(size = 7)) + geom_hline(yintercept=0, col="red", linetype = 2) +
  theme(axis.text.y = element_text(face = "bold.italic")) +
  facet_grid(~label) +
  theme(strip.background = element_rect(fill="azure3")) +
  theme(strip.text = element_text(size=10, face="bold"))

# save plot 
ggsave(
  filename = "P.Deg_Bacteria_DESEq2.Plast.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 6.5,
  height = 13,
  units = c("in"),
  dpi = 300)

#### Compare BULK SOILS with plastic WITH and WITHOUT FUNGI ####

rare.dseq.bulk.FnoF <- subset_samples(rare.dseq, Compartment == "Bulk")
rare.dseq.bulk.FnoF <- subset_samples(rare.dseq.bulk.FnoF, Treatment == "Plastic+Fungus" |
                                        Treatment == "Plastic_NoFungus")

dds.bulk.FnoF <- phyloseq_to_deseq2(rare.dseq.bulk.FnoF, ~ Treatment)

#calculate size factors using edgeR
library(edgeR)
sizeFactors(dds.bulk.FnoF) <- calcNormFactors(counts(dds.bulk.FnoF))

#run DESeq function
dds.bulk.FnoF = DESeq(dds.bulk.FnoF, test="Wald", fitType="parametric")

res.bulk.FnoF <- results(dds.bulk.FnoF, cooksCutoff = FALSE)
alpha = 0.1
sigtab.bulk.FnoF <- res[which(res.bulk.FnoF$padj < alpha), ]
sigtab.bulk.FnoF <- cbind(as(sigtab.bulk.FnoF, "data.frame"), 
                          as(tax_table(rare.dseq.bulk.FnoF)[rownames(sigtab.bulk.FnoF), ], "matrix"))

#prepare data to plot DESeq2 output
x = tapply(sigtab.bulk.FnoF$log2FoldChange, sigtab.bulk.FnoF$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab.bulk.FnoF$Phylum = factor(as.character(sigtab.bulk.FnoF$Phylum), levels=names(x))
x = tapply(sigtab.bulk.FnoF$log2FoldChange, sigtab.bulk.FnoF$Genus, function(x) max(x))
sigtab.bulk.FnoF$Genus = factor(as.character(sigtab.bulk.FnoF$Genus), levels=names(x))


library(ggplot2)
sigtab.bulk.FnoF <- subset(sigtab.bulk.FnoF, Genus !="NA")
sigtab.bulk.FnoF$label <- "Bulk Soil (Plastic+Fungus vs. Plastic Only)"
ggplot(sigtab.bulk.FnoF, aes(x=Genus, y=log2FoldChange, color=Phylum)) + 
  geom_point(size=4, alpha = 0.6) + theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + 
  theme_bw() + coord_flip() + scale_color_manual(values = earthy_colors) + xlab("") + 
  ylab("[Log2] Fold-Change") + theme(axis.text.y = element_text(size = 7)) + 
  theme(plot.title = element_text(size = 7)) + geom_hline(yintercept=0, col="red", linetype = 2) +
  theme(axis.text.y = element_text(face = "bold.italic")) +
  facet_grid(~label) +
  theme(strip.background = element_rect(fill="burlywood4")) +
  theme(strip.text = element_text(size=4, face="bold"))

# save plot 
ggsave(
  filename = "P.Deg_Bacteria_DESEq2.Bulk.PF.PnoF.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 8.5,
  height = 15,
  units = c("in"),
  dpi = 300)


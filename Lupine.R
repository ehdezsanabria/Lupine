# source("https://bioconductor.org/biocLite.R")
# 
# biocLite("phyloseq")
# biocLite("Biobase")

devtools::install_github("karthik/wesanderson")

library("phyloseq")
library("ggplot2")
library("RColorBrewer")
library("SpiecEasi")
library("dplyr")
library("Phenoflow")
library("grid")
library("vegan")
library("igraph")
library("DESeq2")
source("functions.R")

#### Convert the OTU table from csv to phyloseq object, to make the abundance plots ####

# Add interaction factor levels to pool replicates
metadata <- read.csv('Metadata Rhizosheath_lupine.csv')

metadata <- data.frame(metadata)
metadata$merge_factor <- interaction(metadata$Timepoint, metadata$Fertilizer)

metadata$Timepoint <- plyr::revalue(metadata$Timepoint, replace = 
                                c("T2" = "First harvest",
                                  "T3" = "Second harvest"))
metadata$Timepoint <- factor(as.character(metadata$Timepoint), levels = c("T2" = "First harvest",
                                                              "T3" = "Second harvest"))

# Create phyloseq object from csv file

emma_df <- read.csv("./Counts_rhizosheath.csv")

# Taxonomy table
tax_df <-emma_df %>%  dplyr::select(OTU:Genus)
rownames(tax_df) <- tax_df$OTU


# Create OTU table with all samples
otu_df <- emma_df %>% dplyr::select(Oliver35:Oliver133)

# Name the rows with the names of the OTUs
rownames(otu_df) <- tax_df$OTU

tax_df <- tax_df %>% dplyr::select(-OTU)

# Make phyloseq object
physeq <- phyloseq(otu_table(otu_df, taxa_are_rows = TRUE),
                   tax_table(as.matrix(tax_df)))

colnames(tax_table(physeq)) <- c("Kingdom", "Phylum" ,"Class", "Order", "Family", "Genus")
rownames(metadata) <- metadata$Sample # make sure rownames are matching with physeq
sample_data(physeq) <- sample_data(metadata) # add metadata

# Merge replicate samples for visualization purposes
phy_genus <- merge_samples(physeq, "merge_factor")

# Check library sizes
hist(sample_sums(phy_genus), breaks = 100)
print(sample_sums(phy_genus))

# rescale to minimum sampling depth
phy_genus_scale <- scale_reads(phy_genus)
print(sample_sums(phy_genus_scale))

# Pool at genus level
phy_genus_scale = tax_glom(phy_genus_scale, "Genus")

# Calculate relative abundances
phy_genus_scale <- transform_sample_counts(phy_genus_scale, function(x) 100*x/sum(x))

# Select top 12 Genera
TopNGenera <- names(sort(taxa_sums(phy_genus_scale), TRUE)[1:12])
tax_table(phy_genus_scale)[!taxa_names(phy_genus_scale) %in% TopNGenera, "Genus"] <- "Other"

# psmelt into dataframe
df_phy_genus_scale_pruned <- psmelt(phy_genus_scale)

# Add the following part if the version of phyloseq is giving issues and does not add metadata
# correctly

# The following instractions are to replace metadata again 
# since the merging of samples removes this information

df_phy_genus_scale_pruned <- Filter(function(x)!all(is.na(x)), df_phy_genus_scale_pruned)
df_phy_genus_scale_pruned <- df_phy_genus_scale_pruned[, -c(4:11)]
df_phy_genus_scale_pruned <- dplyr::left_join(df_phy_genus_scale_pruned, metadata[-c(1)],
                                              by = c("Sample" = "merge_factor")) %>% distinct()

# adjust order Genus label so that "other" group is last and the other
# Genera are ordered according to their total abundance across data set
sum_table_genera <- df_phy_genus_scale_pruned %>% group_by(Genus) %>% summarize(sum_abund = sum(Abundance))

names_sorted_genera <- as.character(sum_table_genera$Genus)[order(sum_table_genera$sum_abund,
                                                                  decreasing = TRUE)]
names_sorted_genera <- c(names_sorted_genera[names_sorted_genera != "Other"],
                         names_sorted_genera[names_sorted_genera == "Other"]) # Put "other" group last

df_phy_genus_scale_pruned$Genus <- factor(as.character(df_phy_genus_scale_pruned$Genus),
                                          levels = names_sorted_genera)

# Finally we merge the Genera present in "Others"
df_phy_genus_scale_pruned_remerged <- df_phy_genus_scale_pruned %>% 
  group_by(Genus, Sample) %>% summarize(sum_abund = sum(Abundance))

# Add metadata AGAIN
df_phy_genus_scale_pruned_remerged <- left_join(df_phy_genus_scale_pruned_remerged,
                                                metadata[-c(1)], 
                                                by = c("Sample" = "merge_factor")) %>% distinct()

# Make barplot #
library (RColorBrewer)
#getPalette <- colorRampPalette(brewer.pal(9, "Set1"))
#getPalette <- colorRampPalette(brewer.pal(12, "Set3"))
getPalette <- colorRampPalette(brewer.pal(12, "Paired"))
#getPalette <- colorRampPalette(wes_palette(n=5, "Zissou1"))

# Generate a plot with the relative abundances pooled by location
library(ggplot2)
library(wesanderson)


plot.taxa <- ggplot(aes(Fertilizer, y = sum_abund, 
                        fill = Genus, color = Genus), data = df_phy_genus_scale_pruned_remerged)+
  geom_bar(stat="identity", alpha = 1)+
  scale_color_manual(values = rep("gray50", 13))+
  facet_grid(Timepoint~.)+
  ylab("Relative Abundance (%)")+ xlab("")+
  scale_fill_manual(values = (getPalette(13)))+
  theme_bw()+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=16),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        strip.text = element_text(size = 14))

print(plot.taxa)

################################################################################
### Alpha diversity analysis
################################################################################

# Species richness and eveness indices
library(vegan)
library(ggplot2)
counts <- t(otu_df)
shannon <- diversity(counts)
totalspecies<-specnumber(counts)
Pielou <- diversity(counts)/log(totalspecies)

# Calculate diversity from rescaled OTU table
diversity_results <- Diversity_16S(scale_reads(physeq), brea = TRUE, thresh = 1000, R = 100)
diversity_results <- data.frame(Sample = rownames(diversity_results), diversity_results)

# Merge diversity results with metadata
diversity_results <- left_join(diversity_results, metadata, by = c("Sample")) %>% 
  distinct()

# Plot diversity results
## Overall difference between locations - D2 (Inverse Simpson)
p.div1_D2 <- ggplot(aes(x = Timepoint, y = D2, fill = Fertilizer), data = diversity_results)+
  geom_jitter(shape = 21, size = 4, position=position_dodge(0.8),
              aes(shape=Fertilizer))+
  geom_boxplot(alpha = 0.7, outlier.shape = NA)+
  ylab("Inverse Simpson")+ xlab("Timepoint")+
  scale_fill_manual(values=brewer.pal(n=8,"Accent"))+
  theme_bw()+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=16),
        axis.text.x = element_text(size = 14),
        strip.text = element_text(size = 14),
        legend.position = "right")


print(p.div1_D2)

## Observed diversity

p.div1_D0 <- ggplot(aes(x = Timepoint, y = D0, fill = Fertilizer), data = diversity_results)+
  geom_jitter(shape = 21, size = 4, position=position_dodge(0.8),
              aes(shape=Fertilizer))+
  geom_boxplot(alpha = 0.7, outlier.shape = NA)+
  ylab("Richness")+ xlab("Timepoint")+
  scale_fill_manual(values=brewer.pal(n=8,"Accent"))+
  theme_bw()+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=16),
        axis.text.x = element_text(size = 14),
        strip.text = element_text(size = 14),
        legend.position = "right")
print(p.div1_D0)

# Export diversity values with metadata
indices <- cbind(Pielou, shannon, totalspecies)

# You can use the values in these files for statistical analysis
# Try to run a mixed model in SAS or GraphPad or any other software
write.csv(file = "evenness_roots lupine.csv", indices)
write.csv(file = "diversity_roots lupine.csv", diversity_results)

## Plot Evenness

# Change factor names

evenness_results <- read.csv('evenness_roots lupine.csv', stringsAsFactors = TRUE)

evenness_results$Timepoint <- plyr::revalue(evenness_results$Timepoint, replace = 
                                     c("T2" = "First harvest",
                                       "T3" = "Second harvest"))

p.div1_ev <- ggplot(aes(x = Timepoint, y = Pielou, fill = Fertilizer), data = evenness_results)+
  geom_jitter(shape = 21, size = 4, position=position_dodge(0.8),
              aes(shape=Fertilizer))+
  geom_boxplot(alpha = 0.7, outlier.shape = NA)+
  ylab("Evenness")+ xlab("Timepoint")+
  scale_fill_manual(values=brewer.pal(n=8,"Accent"))+
  theme_bw()+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=16),
        axis.text.x = element_text(size = 14),
        strip.text = element_text(size = 14),
        legend.position = "right")

print(p.div1_ev)




################################################################################
### Beta diversity analysis
################################################################################

# Add metadata to phyloseq object

sample_data(physeq) <- metadata

# Rescale OTU table to account for library size differences
df_phy_scaled <- scale_reads(physeq)

# Start making the PCoA
pcoa <- ordinate(
  physeq = df_phy_scaled, 
  method = "PCoA", 
  distance = "bray",
  correction = "lingoes",
  k=2
)

# Now lets transform this to dataframe
pcoa.df <- data.frame(pcoa$vectors, sample_data(df_phy_scaled))

# And calculate the variance explained
var <- round(pcoa$values$Eigenvalues/sum(pcoa$values$Eigenvalues)*100,1)

# Lets run an exploratory permanova to see suggested effect sizes.
# Similar variances across treatments: 
dist.seq <- vegan::vegdist(t(otu_table(df_phy_scaled)))
disper.seq_tpt <- vegan::betadisper(dist.seq, group = sample_data(df_phy_scaled)$Timepoint)

anova(disper.seq_tpt)
print(disper.seq_tpt)

disper.seq_fer <- vegan::betadisper(dist.seq, group = sample_data(df_phy_scaled)$Fertilizer)

anova(disper.seq_fer)
print(disper.seq_fer)

# Permutations are constrained within each Location 
perm_results <- vegan::adonis(dist.seq ~ Timepoint + Fertilizer, 
                              data = data.frame(sample_data(df_phy_scaled)))

# Add this information on plots
my_grob = grid::grobTree(textGrob(bquote(paste(r[Timepoint]^2 == 
                                                 .(round(100 * perm_results$aov.tab[1, 5], 1)), 
                                               "%")), x = 0.7, y = 0.87, 
                                  hjust = 0, gp = gpar(col = "black", 
                                                       fontsize = 14, fontface = "italic")))
my_grob2 = grid::grobTree(textGrob(bquote(paste(r[Fertilizer]^2 == 
                                                  .(format(round(100 * perm_results$aov.tab[2, 5], 
                                                                 1), nsmall = 1)), "%")), x = 0.7, y = 0.80, 
                                   hjust = 0, gp = gpar(col = "black", fontsize = 14, 
                                                        fontface = "italic")))
# Now we can plot the beta diversity plot
library(wesanderson)
p_beta_roots <- ggplot(aes(x=Axis.1, y=Axis.2, color=Fertilizer), data = pcoa.df)+
  geom_point(alpha=0.7, size=6)+
  scale_shape_manual(values = c(23))+
  theme_bw()+
  scale_color_manual(values=wes_palette(n=3, name="FantasticFox1"))+
  #scale_fill_brewer(palette = "Dark2")+
  labs(x = paste0("PCoA axis 1 (",var[1], "%)"), 
       y = paste0("PCoA axis 2 (",var[2], "%)"), 
       colour="")+  
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=16))+
  guides(fill = guide_legend(override.aes = list(shape = 22)))+
  annotation_custom(my_grob)+
  annotation_custom(my_grob2)+
  facet_grid(.~Timepoint)+
    theme(strip.text.x = element_text(size=16))

print(p_beta_roots)

tiff(file = "Betadispersion rhizosheath.tiff", width = 12, height = 9, units = "in", res = 300)
print(p_beta_roots)
dev.off()

png(p_beta_bulk, file = "Betadispersion rhizosheath.png", width = 12, height = 9, units = "in", res = 300)
print(p_beta_roots)
dev.off()

###############################################################
##################### BULK SOIL ###############################
###############################################################

#### Convert the OTU table from csv to phyloseq object, to make the abundance plots ####

# Add interaction factor levels to pool replicates
metadata <- read.csv('Metadata Rhizosphere_lupine2.csv')

metadata <- data.frame(metadata)
metadata$merge_factor <- interaction(metadata$Condition, metadata$Plant,
                                     metadata$Timepoint, metadata$Fertilizer)

metadata$Timepoint <- plyr::revalue(metadata$Timepoint, replace = 
                                      c("T1" = "Seeding",
                                        "T2" = "First harvest",
                                        "T3" = "Second harvest"))
metadata$Timepoint <- factor(as.character(metadata$Timepoint), levels = c("T1" = "Seeding",
                                                                          "T2" = "First harvest",
                                                                          "T3" = "Second harvest"))

# Create phyloseq object from csv file

emma_df <- read.csv("./Counts_rhizosphere.csv")

# Taxonomy table
tax_df <-emma_df %>%  dplyr::select(OTU:Genus)
rownames(tax_df) <- tax_df$OTU


# Create OTU table with all samples
otu_df <- emma_df %>% dplyr::select(Oliver1:Oliver207)

# Name the rows with the names of the OTUs
rownames(otu_df) <- tax_df$OTU

tax_df <- tax_df %>% dplyr::select(-OTU)

# Make phyloseq object
physeq <- phyloseq(otu_table(otu_df, taxa_are_rows = TRUE),
                   tax_table(as.matrix(tax_df)))

colnames(tax_table(physeq)) <- c("Kingdom", "Phylum" ,"Class", "Order", "Family", "Genus")
rownames(metadata) <- metadata$Sample # make sure rownames are matching with physeq
sample_data(physeq) <- sample_data(metadata) # add metadata

# Merge replicate samples for visualization purposes
phy_genus <- merge_samples(physeq, "merge_factor")

# Check library sizes
hist(sample_sums(phy_genus), breaks = 100)
print(sample_sums(phy_genus))

# rescale to minimum sampling depth
phy_genus_scale <- scale_reads(phy_genus)
print(sample_sums(phy_genus_scale))

# Pool at genus level
phy_genus_scale = tax_glom(phy_genus_scale, "Genus")

# Calculate relative abundances
phy_genus_scale <- transform_sample_counts(phy_genus_scale, function(x) 100*x/sum(x))

# Select top 12 Genera
TopNGenera <- names(sort(taxa_sums(phy_genus_scale), TRUE)[1:12])
tax_table(phy_genus_scale)[!taxa_names(phy_genus_scale) %in% TopNGenera, "Genus"] <- "Other"

# psmelt into dataframe
df_phy_genus_scale_pruned <- psmelt(phy_genus_scale)

# Add the following part if the version of phyloseq is giving issues and does not add metadata
# correctly

# The following instractions are to replace metadata again 
# since the merging of samples removes this information

#df_phy_genus_scale_pruned <- Filter(function(x)!all(is.na(x)), df_phy_genus_scale_pruned)
#df_phy_genus_scale_pruned <- df_phy_genus_scale_pruned[, -c(4:11)]
#df_phy_genus_scale_pruned <- dplyr::left_join(df_phy_genus_scale_pruned, metadata[-c(1)],
#                                             by = c("Sample" = "merge_factor")) %>% distinct()

# adjust order Genus label so that "other" group is last and the other
# Genera are ordered according to their total abundance across data set
sum_table_genera <- df_phy_genus_scale_pruned %>% group_by(Genus) %>% summarize(sum_abund = sum(Abundance))

names_sorted_genera <- as.character(sum_table_genera$Genus)[order(sum_table_genera$sum_abund,
                                                                  decreasing = TRUE)]
names_sorted_genera <- c(names_sorted_genera[names_sorted_genera != "Other"],
                         names_sorted_genera[names_sorted_genera == "Other"]) # Put "other" group last

df_phy_genus_scale_pruned$Genus <- factor(as.character(df_phy_genus_scale_pruned$Genus),
                                          levels = names_sorted_genera)

# Finally we merge the Genera present in "Others"
df_phy_genus_scale_pruned_remerged <- df_phy_genus_scale_pruned %>% 
  group_by(Genus, Sample) %>% summarize(sum_abund = sum(Abundance))

# Add metadata AGAIN
df_phy_genus_scale_pruned_remerged <- left_join(df_phy_genus_scale_pruned_remerged,
                                                metadata[-c(1)], 
                                                by = c("Sample" = "merge_factor")) %>% distinct()

# Make barplot #
library (RColorBrewer)
#getPalette <- colorRampPalette(brewer.pal(9, "Set1"))
#getPalette <- colorRampPalette(brewer.pal(12, "Set3"))
getPalette <- colorRampPalette(brewer.pal(12, "Paired"))
#getPalette <- colorRampPalette(wes_palette(n=5, "Zissou1"))

# Generate a plot with the relative abundances pooled by location
library(ggplot2)
library(wesanderson)


plot.taxa1 <- df_phy_genus_scale_pruned_remerged %>% 
  dplyr::filter(Plant == "No_plant") %>% 
  ggplot(aes(Fertilizer, y = sum_abund, 
                        fill = Genus, color = Genus), data = .)+
  geom_bar(stat="identity", alpha = 1)+
  scale_color_manual(values = rep("gray50", 13))+
  facet_grid(Condition~Timepoint)+
  ylab("Relative Abundance (%)")+ xlab("")+
  scale_fill_manual(values = rev(getPalette(13)))+
  theme_bw()+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=16),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 14))

print(plot.taxa1)

plot.taxa2 <- df_phy_genus_scale_pruned_remerged %>% 
  dplyr::filter(Plant == "Lupine") %>% 
  ggplot(aes(Fertilizer, y = sum_abund, 
             fill = Genus, color = Genus), data = .)+
  geom_bar(stat="identity", alpha = 1)+
  scale_color_manual(values = rep("gray50", 13))+
  facet_grid(Plant~Timepoint)+
  ylab("Relative Abundance (%)")+ xlab("")+
  scale_fill_manual(values = rev(getPalette(13)))+
  theme_bw()+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=16),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 14))

print(plot.taxa2)


tiff(file = "Relative abundance rhizosphere.tiff", width = 9, height = 12, units = "in", res = 300)
cowplot::plot_grid(plot.taxa1, plot.taxa2, ncol = 1)
dev.off()

png(file = "Relative abundance rhizosphere.tiff", width = 9, height = 12, units = "in", res = 300)
cowplot::plot_grid(plot.taxa1, plot.taxa2, ncol = 1)
dev.off()

################################################################################
### Alpha diversity analysis
################################################################################

# Species richness and eveness indices
library(vegan)
library(ggplot2)
counts <- t(otu_df)
shannon <- diversity(counts)
totalspecies<-specnumber(counts)
Pielou <- diversity(counts)/log(totalspecies)

# Calculate diversity from rescaled OTU table
library(Phenoflow)
library(dplyr)
library(RColorBrewer)
diversity_results <- Diversity_16S(scale_reads(physeq), brea = TRUE, thresh = 1000, R = 100)
diversity_results <- data.frame(Sample = rownames(diversity_results), diversity_results)

# Merge diversity results with metadata
diversity_results <- left_join(diversity_results, metadata, by = c("Sample")) %>% 
  distinct()

# Plot diversity results
## Overall difference between locations - D2 (Inverse Simpson)
plot.D2_1 <- diversity_results %>% 
  dplyr::filter(Plant == "No_plant") %>% 
  ggplot(aes(Timepoint, y = D2, 
             fill = Fertilizer), data = .)+
  geom_jitter(shape = 21, size = 4, position=position_dodge(0.8),
              aes(shape=Fertilizer))+
  geom_boxplot(alpha = 0.7, outlier.shape = NA)+
  ylab("Inverse Simpson")+ xlab(".")+
  scale_fill_manual(values=brewer.pal(n=8,"Accent"))+
  facet_grid(Condition~Plant)+
  theme_bw()+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=16),
        axis.text.x = element_text(size = 14, angle = 30, hjust = 1),
        strip.text = element_text(size = 14),
        legend.position = "right")


print(plot.D2_1)

plot.D2_2 <- diversity_results %>% 
  dplyr::filter(Plant == "Lupine") %>% 
  ggplot(aes(Timepoint, y = D2, 
             fill = Fertilizer), data = .)+
  geom_jitter(shape = 21, size = 4, position=position_dodge(0.8),
              aes(shape=Fertilizer))+
  geom_boxplot(alpha = 0.7, outlier.shape = NA)+
  ylab("Inverse Simpson")+ xlab(".")+
  scale_fill_manual(values=brewer.pal(n=8,"Accent"))+
  facet_grid(Condition~Plant)+
  theme_bw()+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=16),
        axis.text.x = element_text(size = 14, angle = 30, hjust = 1),
        strip.text = element_text(size = 14),
        legend.position = "right")


print(plot.D2_2)

tiff(file = "Alpha diversity rhizosphere.tiff", width = 9, height = 12, units = "in", res = 300)
cowplot::plot_grid(plot.D2_1, plot.D2_2, ncol = 1)
dev.off()

png(file = "Alpha diversity rhizosphere.png", width = 9, height = 12, units = "in", res = 300)
cowplot::plot_grid(plot.D2_1, plot.D2_2, ncol = 1)
dev.off()

## Observed diversity

plot.D0_1 <- diversity_results %>% 
  dplyr::filter(Plant == "No_plant") %>% 
  ggplot(aes(Timepoint, y = D0, 
             fill = Fertilizer), data = .)+
  geom_jitter(shape = 21, size = 4, position=position_dodge(0.8),
              aes(shape=Fertilizer))+
  geom_boxplot(alpha = 0.7, outlier.shape = NA)+
  ylab("Richness")+ xlab(".")+
  scale_fill_manual(values=brewer.pal(n=8,"Accent"))+
  facet_grid(Condition~Plant)+
  theme_bw()+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=16),
        axis.text.x = element_text(size = 14,angle = 30, hjust = 1),
        strip.text = element_text(size = 14),
        legend.position = "right")
print(plot.D0_1)

plot.D0_2 <- diversity_results %>% 
  dplyr::filter(Plant == "Lupine") %>% 
  ggplot(aes(Timepoint, y = D0, 
             fill = Fertilizer), data = .)+
  geom_jitter(shape = 21, size = 4, position=position_dodge(0.8),
              aes(shape=Fertilizer))+
  geom_boxplot(alpha = 0.7, outlier.shape = NA)+
  ylab("Richness")+ xlab(".")+
  scale_fill_manual(values=brewer.pal(n=8,"Accent"))+
  facet_grid(Condition~Plant)+
  theme_bw()+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=16),
        axis.text.x = element_text(size = 14,angle = 30, hjust = 1),
        strip.text = element_text(size = 14),
        legend.position = "right")
print(plot.D0_2)

tiff(file = "Richness.tiff", width = 9, height = 12, units = "in", res = 300)
cowplot::plot_grid(plot.D0_1, plot.D0_2, ncol = 1)
dev.off

png(file = "Richness.png", width = 9, height = 12, units = "in", res = 300)
cowplot::plot_grid(plot.D0_1, plot.D0_2, ncol = 1)
dev.off()

# Export diversity values with metadata
indices <- cbind(Pielou, shannon, totalspecies)

# You can use the values in these files for statistical analysis
# Try to run a mixed model in SAS or GraphPad or any other software
write.csv(file = "evenness_rhizosphere lupine.csv", indices)
write.csv(file = "diversity_rhizosphere lupine.csv", diversity_results)

## Plot Evenness

# Change factor names

evenness_results <- read.csv('evenness_rhizosphere lupine.csv', stringsAsFactors = TRUE)

evenness_results$Timepoint <- plyr::revalue(metadata$Timepoint, replace = 
                                      c("T1" = "Seeding",
                                        "T2" = "First harvest",
                                        "T3" = "Second harvest"))
evenness_results$Timepoint <- factor(as.character(metadata$Timepoint), levels = c("T1" = "Seeding",
                                                                          "T2" = "First harvest",
                                                                          "T3" = "Second harvest"))

plot.ev1 <- evenness_results %>% 
  dplyr::filter(Plant == "No_plant") %>% 
  ggplot(aes(Timepoint, y = Pielou,fill = Fertilizer), data = .)+
  geom_jitter(shape = 21, size = 4, position=position_dodge(0.8),
              aes(shape=Fertilizer))+
  geom_boxplot(alpha = 0.7, outlier.shape = NA)+
  ylab("Evenness")+ xlab(".")+
  scale_fill_manual(values=brewer.pal(n=8,"Accent"))+
  facet_grid(Condition~Plant)+
  theme_bw()+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=16),
        axis.text.x = element_text(size = 14,angle = 30, hjust = 1),
        strip.text = element_text(size = 14),
        legend.position = "right")

print(plot.ev1)

plot.ev2 <- evenness_results %>% 
  dplyr::filter(Plant == "Lupine") %>% 
  ggplot(aes(Timepoint, y = Pielou,fill = Fertilizer), data = .)+
  geom_jitter(shape = 21, size = 4, position=position_dodge(0.8),
              aes(shape=Fertilizer))+
 geom_boxplot(alpha = 0.7, outlier.shape = NA)+
  ylab("Evenness")+ xlab(".")+
  scale_fill_manual(values=brewer.pal(n=8,"Accent"))+
  facet_grid(Condition~Plant)+
  theme_bw()+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=16),
        axis.text.x = element_text(size = 14,angle = 30, hjust = 1),
        strip.text = element_text(size = 14),
        legend.position = "right")

print(plot.ev2)


tiff(file = "Evenness rhizosphere.tiff", width = 9, height = 12, units = "in", res = 300)
cowplot::plot_grid(plot.ev1, plot.ev2, ncol = 1)
dev.off()

png(file = "Evenness rhizosphere.png", width = 9, height = 12, units = "in", res = 300)
cowplot::plot_grid(plot.ev1, plot.ev2, ncol = 1)
dev.off()

################################################################################
### Beta diversity analysis
################################################################################

# Add metadata to phyloseq object

sample_data(physeq) <- metadata

# Rescale OTU table to account for library size differences
df_phy_scaled <- scale_reads(physeq)

# Start making the PCoA
pcoa <- ordinate(
  physeq = df_phy_scaled, 
  method = "PCoA", 
  distance = "bray",
  correction = "lingoes",
  k=2
)

# Now lets transform this to dataframe
pcoa.df <- data.frame(pcoa$vectors, sample_data(df_phy_scaled))

# And calculate the variance explained
var <- round(pcoa$values$Eigenvalues/sum(pcoa$values$Eigenvalues)*100,1)

# Lets run an exploratory permanova to see suggested effect sizes.
# Similar variances across treatments: 
dist.seq <- vegan::vegdist(t(otu_table(df_phy_scaled)))
disper.seq_tpt <- vegan::betadisper(dist.seq, group = sample_data(df_phy_scaled)$Timepoint)

anova(disper.seq_tpt)
print(disper.seq_tpt)

disper.seq_fer <- vegan::betadisper(dist.seq, group = sample_data(df_phy_scaled)$Fertilizer)

anova(disper.seq_fer)
print(disper.seq_fer)

disper.seq_con <- vegan::betadisper(dist.seq, group = sample_data(df_phy_scaled)$Condition)

anova(disper.seq_con)
print(disper.seq_con)

disper.seq_pla <- vegan::betadisper(dist.seq, group = sample_data(df_phy_scaled)$Plant)

anova(disper.seq_pla)
print(disper.seq_pla)

# Permutations are constrained within each Location 
perm_results <- vegan::adonis(dist.seq ~ Timepoint + Fertilizer + Condition + Plant, 
                              data = data.frame(sample_data(df_phy_scaled)))

# Add this information on plots
my_grob = grid::grobTree(textGrob(bquote(paste(r[Timepoint]^2 == 
                                                 .(round(100 * perm_results$aov.tab[1, 5], 1)), 
                                               "%")), x = 0.55, y = 0.40, 
                                  hjust = 0, gp = gpar(col = "black", 
                                                       fontsize = 12, fontface = "italic")))
my_grob2 = grid::grobTree(textGrob(bquote(paste(r[Fertilizer]^2 == 
                                                  .(format(round(100 * perm_results$aov.tab[2, 5],1), 
                                                           nsmall = 1)), "%")), x = 0.55, y = 0.35, 
                                   hjust = 0, gp = gpar(col = "black", fontsize = 12, 
                                                        fontface = "italic")))
my_grob3 = grobTree(textGrob(bquote(paste(r[Condition]^2 == 
                                            .(round(100 * perm_results$aov.tab[3, 5], 1)), 
                                          "%")), x = 0.55, y = 0.30, hjust = 0, gp = gpar(col = "black", 
                                                                                         fontsize = 12, fontface = "italic")))
my_grob4 = grobTree(textGrob(bquote(paste(r[Plant]^2 == 
                                            .(round(100 * perm_results$aov.tab[4, 5]),1), 
                                          "%")), x = 0.55, y = 0.25, hjust = 0, gp = gpar(col = "black", 
                                                                                         fontsize = 12, fontface = "italic")))
# Now we can plot the beta diversity plot
library(wesanderson)
p_beta_bulk <- ggplot(data = pcoa.df, aes(x=Axis.1, 
                                          y=Axis.2,shape=Plant))+
  geom_point(alpha=0.7, size=4, aes(fill=Fertilizer))+
  scale_shape_manual(values = c(21,23))+
  theme_bw()+
  scale_fill_manual(values=wes_palette(n=3, name="FantasticFox1"))+
  #scale_fill_brewer(palette = "Accent")+
  labs(x = paste0("PCoA axis 1 (",var[1], "%)"), 
       y = paste0("PCoA axis 2 (",var[2], "%)"), 
       colour="")+  
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=16))+
  guides(fill = guide_legend(override.aes = list(shape = 22)))+
  annotation_custom(my_grob)+
  annotation_custom(my_grob2)+
  annotation_custom(my_grob3)+
  annotation_custom(my_grob4)+
  facet_grid(.~Timepoint)+
    theme(strip.text.x = element_text(size=16))



tiff(file = "Betadispersion rhizosphere.tiff", width = 12, height = 9, units = "in", res = 300)
print(p_beta_bulk)
dev.off()

png(p_beta_bulk, file = "Betadispersion rhizosphere.png", width = 12, height = 9, units = "in", res = 300)
print(p_beta_bulk)
dev.off()
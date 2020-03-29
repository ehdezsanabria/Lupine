## Run DEseq2 between Donors (default pipeline)
################################################################################

### Analysis between donors for TPT = 16h controlling for treatment ###
df_phy_T3 <- prune_samples(sample_data(df_phy)$TPT == "16h", df_phy)
df_phy_T16CX <- prune_samples(sample_data(df_phy_T16)$TRT == "CX", df_phy_T16)

diagdds <- phyloseq_to_deseq2(df_phy_T16CX, ~ Donor + TRT)
diagdds <- DESeq2::DESeq(diagdds, test = "Wald", fitType = "parametric")

# Summarize results
res = DESeq2::results(dds, cooksCutoff = FALSE)
alpha = 0.01
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(df_phy_T16)[rownames(sigtab), ], "matrix"))
head(sigtab)

# Make default result plot
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}

# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))

p_deseq_donor_T16 <- ggplot(sigtab, aes(x=Genus, y=log2FoldChange, fill=Genus)) + 
  geom_point(shape = 21, size=6, alpha = 0.7, color = "black") + 
  theme(axis.text.x = element_text(angle = -60, hjust = 0, vjust=1))+
  scale_color_brewer(palette = "Paired")+
  ylab("Log2 fold change")

png(file = "./Figures_DESEQ/DESEQ-DONOR-P0.01-T16.png", width = 8, height = 6, res = 500, units = "in")
print(p_deseq_donor_T16)
dev.off()



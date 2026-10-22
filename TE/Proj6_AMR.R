#load packages
library(phyloseq);library(metagenomeSeq);library(dplyr);library(scales);
library(pairwiseAdonis); library(vegan); library(metagMisc); library(stringr)
library(ggplot2); library(btools); library(randomcoloR); library(cowplot)
library(pairwiseAdonis); library(picante); library(gridExtra); library(plyr)
library(ggalt); library(ggforce); library(concaveman); library(ggdendro);
library(microbiome)

# setwd
setwd("/Users/ljpinnell/Documents/VERO/Project6/TE/phyloseq/")

# import data
# elk
elk <- import_biom("elk_AMR_counts.biom")
bison_other <- import_biom("Bison_AMR_counts.biom")
ybf2 <- import_biom("YBF2_AMR_counts.biom")
sample_names(ybf2) <- "YBF2"
cow <- import_biom("cow_AMR_counts.biom")
sample_names(cow) <- c("BFe10","BFe3","BFe33","BFe17","BFe28","BFe16","BFe23","BFe8")
map_file <- import_qiime_sample_data("../../metadata/Project6_metadata.txt")

### want the number of classes etc. for just bison for an unassocaited comparison
colnames(tax_table(bison_other)) <- c("Type","Class","Mechanism","Group","Gene","SNP")
bison_other_class <- tax_glom(bison_other, taxrank = "Class")
bison_other_class
bison_other_group <- tax_glom(bison_other, taxrank = "Group")
bison_other_group
colnames(tax_table(ybf2)) <- c("Type","Class","Mechanism","Group","Gene","SNP")
ybf2_class <- tax_glom(ybf2, taxrank = "Class")
ybf2_class
ybf2_group <- tax_glom(ybf2, taxrank = "Group")
ybf2_group

data <- merge_phyloseq(elk, bison_other, ybf2, cow)
data # 1775 AMR genes, 56 samples
data <- merge_phyloseq(data, map_file)
data

# check the names of our ranks
rank_names(data) # "Rank1" - "Rank7" not ideal, lets change em
colnames(tax_table(data)) <- c("Type","Class","Mechanism","Group","Gene","SNP")
rank_names(data) # beauty, now they are named properly

### need to split up the genes needing SNP confirmation:
data_SNPconfirm <- subset_taxa(data, SNP=="Yes", T) # 205 of the 1775 genes
data_noSNP <- subset_taxa(data, SNP=="No") # 1570 of the 1775 genes

# some QC checks
min(sample_sums(data_noSNP)) # 8166
max(sample_sums(data_noSNP)) # 771256
sort(sample_sums(data_noSNP))

## exporting SNP ratios
soil_noSNPs <- subset_samples(data_noSNP, matrix=="Soil")
soil_SNPSs <- subset_samples(data_SNPconfirm, matrix=="Soil")
feces_noSNPs <- subset_samples(data_noSNP, matrix=="Feces")
feces_SNPs <- subset_samples(data_SNPconfirm, matrix=="Feces")

write.csv(sample_sums(soil_noSNPs),"soil_noSNP_counts.csv")
write.csv(sample_sums(soil_SNPSs),"soil_SNPs_counts.csv")
write.csv(sample_sums(feces_noSNPs),"feces_noSNP_counts.csv")
write.csv(sample_sums(feces_SNPs),"feces_SNPs_counts.csv")

# any taxa with no counted reads
any(taxa_sums(data_noSNP)==0) # nope
any(taxa_sums(data_SNPconfirm)==0) # nope

##### Sequencing depth (based on reads passing QC that were put into ARM++)
data.df <- as(sample_data(data), "data.frame")

# just compare the "groups"
ggplot(data.df, aes(x= group, y= NumberOfInputReads, fill = group)) + theme_bw() +
  geom_bar(stat="summary") +
  geom_errorbar(stat = "summary")

pairwise.wilcox.test(data.df$NumberOfInputReads, data.df$group)
# one difference: YNP Bison soil vs. YNP elk feces

four_colour_palette <- c("#638652","#4d5660","#8d391e","#ab6116")
matrix.labs <- c(Feces = "FECES", Soil = "SOIL")

#############################################################################################
##############################         ALPHA DIVERSITY         ##############################
#############################################################################################
#############################################################################################

alpha_div <- estimate_richness(data_noSNP, measures = c("Observed","Shannon","Simpson","InvSimpson"))
alpha_div
alpha_div.df <- as(sample_data(data_noSNP), "data.frame")
alpha_div_meta <- cbind(alpha_div, alpha_div.df)
alpha_div_meta

ggplot(alpha_div_meta, aes(x= location_species, y = Observed, fill = location_species, colour = location_species)) + 
  theme_bw() + facet_wrap(~matrix, ncol=2, labeller = labeller(matrix = matrix.labs)) +
  labs(title= "", y= "No. unique ARGs", x= "") +
  geom_boxplot(alpha = 0.5, size = 1) +
  geom_point(size = 2.5) +
  scale_fill_manual(values = four_colour_palette) +
  scale_colour_manual(values = four_colour_palette) +
  scale_y_continuous(expand = c(0.0271,0,0.2,0), breaks = c(100,400,700,1000)) +
  scale_x_discrete(limits = c("Cow; Farm", "Bison; Yellowstone National Park ","Elk; Yellowstone National Park ","Elk; Rocky Mountain National Park")) +
  theme(legend.position = "none",
        plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
        strip.background = element_blank(),
        strip.text = element_text(size =24, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 28, vjust = 1.75),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(),
        panel.border = element_rect(colour = "black", size = 1.0),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())

ggplot(alpha_div_meta, aes(x= location_species, y = Shannon, fill = location_species, colour = location_species)) + 
  theme_bw() + facet_wrap(~matrix, ncol=2, labeller = labeller(matrix = matrix.labs)) +
  labs(title= "", y= "Shannon", x= "") +
  geom_boxplot(alpha = 0.5, size = 1) +
  geom_point(size = 2.5) +
  scale_fill_manual(values = four_colour_palette) +
  scale_colour_manual(values = four_colour_palette) +
  scale_y_continuous(expand = c(0.1,0,0.16,0), labels = scales::number_format(accuracy = 1.0)) +
  scale_x_discrete(limits = c("Cow; Farm", "Bison; Yellowstone National Park ","Elk; Yellowstone National Park ","Elk; Rocky Mountain National Park")) +
  theme(legend.position = "none",
        plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
        strip.background = element_blank(),
        strip.text = element_text(size =24, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 28, vjust = 1.75),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(),
        panel.border = element_rect(colour = "black", size = 1.0),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())

richness_stats <- pairwise.wilcox.test(alpha_div_meta$Observed, alpha_div_meta$group, p.adjust.method = "BH")
shannon_stats <- pairwise.wilcox.test(alpha_div_meta$Shannon, alpha_div_meta$group, p.adjust.method = "BH")

write.csv(richness_stats[["p.value"]],"richness_stats.csv")
write.csv(shannon_stats[["p.value"]],"shannon_stats.csv")

#############################################################################################
##############################         BETA DIVERSITY         ###############################
#############################################################################################
#############################################################################################

# feces alone to check it
sum(taxa_sums(feces_noSNPs)==0) # 298 taxa not present in feces by itself
feces_noSNPs <- prune_taxa(taxa_sums(feces_noSNPs) > 0, feces_noSNPs)
feces_noSNPS.css <- phyloseq_transform_css(feces_noSNPs, log = F)
# now soil
soil_noSNPs <- prune_taxa(taxa_sums(soil_noSNPs) > 0, soil_noSNPs)
soil_noSNPs.css <- phyloseq_transform_css(soil_noSNPs, log = F)

feces_noSNP.dist <- vegdist(t(otu_table(feces_noSNPS.css)), method = "bray")

feces_noSNP.ord <- vegan::metaMDS(comm = t(otu_table(feces_noSNPS.css)), distance = "bray", try = 20, trymax = 50, autotransform = F)
soil_noSNP.ord <- vegan::metaMDS(comm = t(otu_table(soil_noSNPs.css)), distance = "bray", try = 20, trymax = 50, autotransform = F)


plot_ordination(feces_noSNPS.css,feces_noSNP.ord, type = "samples", color = "location_species") +
  stat_ellipse()

# definitely cow is unique, but not quite as much as 16S

## CSS
data_noSNP.css <- phyloseq_transform_css(data_noSNP, log = F)
data_SNPconfirm.css <- phyloseq_transform_css(data_SNPconfirm, log = F)
data.css <- phyloseq_transform_css(data, log = F)

data_noSNP.df <- as(sample_data(data_noSNP.css), "data.frame")
data_SNPconfirm.df <- as(sample_data(data_SNPconfirm.css), "data.frame")
data.df <- as(sample_data(data.css), "data.frame")

# distance
noSNP.dist <- vegdist(t(otu_table(data_noSNP.css)), method = "bray")

# ORDINATE
noSNP.ord <- vegan::metaMDS(comm = t(otu_table(data_noSNP.css)), distance = "bray", try = 20, trymax = 50, autotransform = F)
SNPconfirm.ord <- ordinate(data_SNPconfirm.css, method = "NMDS", distance = "bray")
data.ord <- ordinate(data.css, method = "NMDS", distance = "bray")

plot_ordination(data_noSNP.css, noSNP.ord, type = "samples", 
                color = "location_species") +
  theme_bw() + facet_wrap(~matrix, ncol=1, labeller = labeller(matrix = matrix.labs)) +
  geom_point(size = 6) +
  #geom_mark_hull(aes(fill=location_species), concavity = 5, expand = 0, radius = 0, size =0.2, alpha = 0.35) +
  stat_ellipse(geom = "polygon", aes(fill=location_species), level = 0.9, lty =4, lwd = 2, alpha= 0.5) +
  scale_fill_manual(values = four_colour_palette) +
  scale_colour_manual(values = four_colour_palette) +
  theme(legend.position = "none",
        plot.margin = unit(c(1,1,1,1), "cm"),
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.border = element_rect(colour = "black", size = 1.5),
        axis.text = element_text(size = 28, colour = "black"),
        axis.title.y = element_text(size = 42, vjust = 1.75),
        axis.title.x = element_text(size = 42, vjust = -1.5))

# stats
# no SNPs
noSNP.adonis <- pairwise.adonis(noSNP.dist, data_noSNP.df$group, perm = 9999, p.adjust.m = "BH")
noSNP.adonis
write.csv(noSNP.adonis, "TE_noSNP_adonis.csv")

noSNP.disper <- betadisper(noSNP.dist, data_noSNP.df$group)
plot(noSNP.disper)
noSNP.permdisp <- permutest(noSNP.disper, permutations = 9999, pairwise = T)
write.csv(noSNP.permdisp[["pairwise"]][["permuted"]],"TE_bray_permdisp.csv")


#### CLUSTERING ON BRAY-CURTIS
noSNP.hclust <- hclust(noSNP.dist, method = "ward.D2")
plot(noSNP.hclust)
noSNP.dendro <- as.dendrogram(noSNP.hclust)
noSNP.dendro.data <- dendro_data(noSNP.dendro, type = "rectangle")
noSNP.dendro.data

### need to add some metadata columns for plotting, easier to export and do in excel
## then copy into a the dataframe
#write.csv(noSNP.dendro.data[["labels"]][["label"]], "TE_noSNP_dendro_sample_order.csv")

species_col <- data.frame(species = c("elk",	"bison",	"elk",	"elk",
                                      "bison",	"elk",	"bison",	"elk",
                                      "elk",	"elk",	"bison",	"elk",
                                      "bison",	"bison",	"elk",	"bison",
                                      "elk",	"elk",	"elk",	"elk",
                                      "elk",	"elk",	"elk", "bison",
                                      "elk",	"elk",	"elk",	"bison",
                                      "elk",	"elk",	"bison",	"bison",
                                      "bison",	"bison",	"elk",	"elk",
                                      "elk",	"bison",	"elk",	"bison",
                                      "cow",	"cow",	"cow",	"cow",
                                      "cow",	"cow",	"cow",	"cow",
                                      "elk",	"bison",	"elk",	"elk",
                                      "elk",	"elk",	"elk",	"elk"))

matrix_col <- data.frame(matrix = c("soil",	"soil",	"soil",	"soil",
                                    "soil",	"soil",	"soil",	"soil",
                                    "soil",	"soil",	"soil",	"soil",
                                    "soil",	"soil",	"soil",	"soil",
                                    "soil",	"soil",	"soil",	"soil",
                                    "soil",	"soil",	"soil",	"soil",
                                    "feces",	"feces",	"feces",	"feces",
                                    "feces",	"feces",	"feces",	"feces",
                                    "feces",	"feces",	"feces",	"feces",
                                    "feces",	"feces",	"feces",	"feces",
                                    "feces",	"feces",	"feces",	"feces",
                                    "feces",	"feces",	"feces",	"feces",
                                    "feces",	"feces",	"feces",	"feces",
                                    "feces",	"feces",	"feces",	"feces"))

location_species_col <- data.frame(location_species = c("elk_RMNP",	"bison_YNP",	"elk_RMNP",	"elk_YNP",
                                                        "bison_YNP",	"elk_YNP",	"bison_YNP",	"elk_RMNP",
                                                        "elk_RMNP",	"elk_YNP",	"bison_YNP",	"elk_RMNP",
                                                        "bison_YNP",	"bison_YNP",	"elk_RMNP",	"bison_YNP",
                                                        "elk_YNP",	"elk_RMNP",	"elk_YNP",	"elk_YNP",
                                                        "elk_RMNP",	"elk_YNP",	"elk_YNP",	"bison_YNP",
                                                        "elk_YNP",	"elk_RMNP",	"elk_RMNP",	"bison_YNP",
                                                        "elk_RMNP",	"elk_YNP",	"bison_YNP",	"bison_YNP",
                                                        "bison_YNP",	"bison_YNP",	"elk_RMNP",	"elk_RMNP",
                                                        "elk_YNP",	"bison_YNP",	"elk_YNP",	"bison_YNP",
                                                        "cow_farm",	"cow_farm",	"cow_farm",	"cow_farm",
                                                        "cow_farm",	"cow_farm",	"cow_farm",	"cow_farm",
                                                        "elk_YNP",	"bison_YNP",	"elk_YNP",	"elk_YNP",
                                                        "elk_YNP",	"elk_RMNP",	"elk_RMNP",	"elk_RMNP"))

location_col <- data.frame(location = c("RMNP",	"YNP",	"RMNP",	"YNP",	"YNP",	"YNP",	"YNP",	"RMNP",	"RMNP",	"YNP",	"YNP",	"RMNP",	"YNP",	"YNP",	"RMNP",	"YNP",	"YNP",	"RMNP",	"YNP",	"YNP",	"RMNP",	"YNP",	"YNP",	"YNP",	"YNP",	"RMNP",	"RMNP",	"YNP",	"RMNP",	"YNP",	"YNP",	"YNP",	"YNP",	"YNP",	"RMNP",	"RMNP",	"YNP",	"YNP",	"YNP",	"YNP",	"farm",	"farm",	"farm",	"farm",	"farm",	"farm",	"farm",	"farm",	"YNP",	"YNP",	"YNP",	"YNP",	"YNP",	"RMNP",	"RMNP",	"RMNP"))

dendro_sample_order <- c("RES3",	"YBS4",	"RES7",	"YES5",	"YBS6",	"YES6",	"YBS8",	"RES1",	"RES5",	"YES1",	"YBS2",	"RES2",	"YBS5",	"YBS3",	"RES6",	"YBS1",	"YES2",	"RES8",	"YES3",	"YES4",	"RES4",	"YES8",	"YES7",	"YBS7",	"YEF2",	"REF6",	"REF4",	"YBF8",	"REF2",	"YEF7",	"YBF4",	"YBF6",	"YBF5",	"YBF1",	"REF3",	"REF5",	"YEF5",	"YBF2",	"YEF1",	"YBF7",	"BFe10",	"BFe16",	"BFe8",	"BFe3",	"BFe33",	"BFe28",	"BFe17",	"BFe23",	"YEF4",	"YBF3",	"YEF8",	"YEF6",	"YEF3",	"REF1",	"REF7",	"REF8")

noSNP.dendro.data$labels <- cbind(noSNP.dendro.data$labels, species_col)
noSNP.dendro.data$labels <- cbind(noSNP.dendro.data$labels, matrix_col)
noSNP.dendro.data$labels <- cbind(noSNP.dendro.data$labels, location_species_col)
noSNP.dendro.data$labels <- cbind(noSNP.dendro.data$labels, location_col)

ggplot(noSNP.dendro.data$segments) + theme_bw() +
  labs(y= "Ward's Distance") +
  geom_segment(aes(x=x,y=y, xend=xend,yend=yend),
               size = 1, lineend = "round", linejoin = "round") +
  geom_point(data = noSNP.dendro.data$labels, aes(x,y, shape = location, colour= location_species, fill = location_species), 
             size = c(8,8,8,8,8,8,8,8,
                      8,8,8,8,8,8,8,8,
                      8,8,8,8,8,8,8,8,
                      8,8,8,8,8,8,8,8,
                      8,8,8,8,8,8,8,8,
                      6,6,6,6,6,6,6,6,
                      8,8,8,8,8,8,8,8), position = position_nudge(y= -0.16), stroke = 1.3) +
  scale_y_continuous(expand = c(0.35,0,0.2,0)) +
  scale_colour_manual(values = four_colour_palette) +
  scale_fill_manual(values = ggplot2::alpha(c(four_colour_palette), 0.5)) +
  scale_shape_manual(values = c(23,22,21)) +
  geom_text(data= noSNP.dendro.data$labels, aes(x,y, label = matrix),
            position = position_nudge(y=-0.4, x=0), size = 8, angle = 90, hjust = 1) +
  theme(legend.position = "none",
        plot.margin = unit(c(1,1,1,1),"cm"),
        panel.border = element_blank(),
        axis.line.y = element_line(size = 1.5, lineend = "square", colour = "black"), 
        axis.ticks.y = element_line(size = 1.5, lineend = "square", colour = "black"),
        axis.ticks.length.y = unit(0.25,"lines"), panel.grid = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(), axis.text.x = element_blank(), 
        axis.title.y = element_text(size = 42, vjust = 3, hjust = 0.59),
        axis.text.y = element_text(size = 28, colour = "black"))

## RA plot for under dendro
ggplot(ra_class_melt, aes(x= Sample, y= Abundance, fill= Class)) +
  theme_bw() +
  labs(y= "Relative Abundance (%)") +
  geom_bar(stat = "summary", colour = "black") +
  scale_x_discrete(limits = dendro_sample_order) +
  scale_fill_manual(values = ra_class_palette) +
  theme(legend.position = "none",
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.y = element_line(size = 1.25, lineend = "square", colour = "black"), 
    axis.ticks.y = element_line(size = 1.5, lineend = "square", colour = "black"),
    axis.ticks.length.y = unit(0.25,"lines"),
    axis.title.y = element_text(size = 42),
    axis.text.y = element_text(size = 28, colour = "black"))

#############################################################################################
############################         RELATIVE ABUNDANCE         #############################
#############################################################################################
#############################################################################################

rel_abund <- transform_sample_counts(data_noSNP.css, function(x) {x/sum(x)}*100)
ra_200 <- prune_taxa(names(sort(taxa_sums(rel_abund),T)[1:200]), rel_abund)
ra_500 <- prune_taxa(names(sort(taxa_sums(rel_abund),T)[1:500]), rel_abund)
ra_75 <- prune_taxa(names(sort(taxa_sums(rel_abund),T)[1:75]), rel_abund)

# agglomerate at different levels
rel_abund_group <- tax_glom(rel_abund, taxrank = "Group") # 497 groups
ra_group_melt <- psmelt(rel_abund_group)
ra_group_palette <- distinctColorPalette(497)
ra_group_50 <- prune_taxa(names(sort(taxa_sums(rel_abund_group),T)[1:50]), rel_abund_group)
ra_group_100 <- prune_taxa(names(sort(taxa_sums(rel_abund_group),T)[1:100]), rel_abund_group)


rel_abund_mechanism <- tax_glom(rel_abund, taxrank = "Mechanism") #115 mechs
ra_mechanism_melt <- psmelt(rel_abund_mechanism)
ra_mechanism_palette <- distinctColorPalette(115)

rel_abund_class <- tax_glom(rel_abund, taxrank = "Class") # 39 classes
ra_class_melt <- psmelt(rel_abund_class)
ra_class_palette <- distinctColorPalette(26)
ra_class_melt_feces <- ra_class_melt[which(ra_class_melt$matrix=="Feces"),]

write.csv(otu_table(rel_abund_class),"class_otus.csv")
write.csv(tax_table(rel_abund_class),"class_taxa.csv")

ggplot(ra_class_melt, aes(x= location_species, y= Abundance, fill= Class)) +
  theme_bw() + 
  facet_wrap(~matrix, scales = "free_x", labeller = labeller(matrix = matrix.labs)) +
  labs(y= "Relative abundance (%)") +
  geom_bar(stat = "summary", colour = "black") +
  scale_fill_manual(values = ra_class_palette) +
  scale_y_continuous(expand = c(0.0005,0,0.0005,0)) +
  scale_x_discrete(limits = c("Cow; Farm", "Bison; Yellowstone National Park ","Elk; Yellowstone National Park ","Elk; Rocky Mountain National Park"),
                   labels = c("Cow: SW USA","Bison; YNP","Elk; YNP","Elk; RMNP")) +
  theme(#legend.position = "none",
    panel.border = element_rect(colour = "black", size = 1),
    plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
    strip.background = element_blank(),
    strip.text = element_text(size =24, colour = "black"),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.text.x = element_text(size =32, colour = "black",angle = 45, hjust = 0.95, vjust = 0.95),
    axis.title.y = element_text(size = 36, vjust = 1.75),
    axis.title.x = element_blank(),
    plot.title = element_text(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank())

### ra plot of just feces

## heatmap at the group level
plot_heatmap(ra_75, "NMDS","bray","group","Class")
# heatmap not great, not going to use


#### Medically important genes

important_AMR_ra <- subset_taxa(rel_abund, Group=="CTX" |
                               Group=="CPH" |
                               Group=="GES" |
                               Group=="IMI" |
                               Group=="KPC" |
                               Group=="SHV" |
                               Group=="TEM" |
                               Group=="IMP" |
                               Group=="NDM" |
                               Group=="CMY" |
                               Group=="OXA" |
                               Group=="MEC" |
                               Group=="MCR" |
                               Group=="VAT" |
                               Group=="VGB" |
                               Group=="CFR" |
                               Group=="VGA" |
                               Group=="SME")


## look at OXA SME specifically at the gene level
oxa_sme <- subset_taxa(important_AMR_ra, Group=="SME" | Group=="OXA")
oxa_sme_melt <- psmelt(oxa_sme)
oxa_sme_palette <- distinctColorPalette(41)

write.csv(otu_table(oxa_sme),"oxa_sem_RA.csv")

ggplot(oxa_sme_melt, aes(x= species, y= Abundance, fill = Gene)) +
  geom_bar(stat = "summary", colour = "black") + facet_grid(~matrix) +
  scale_fill_manual(values = oxa_sme_palette)

important_AMR_ra_group <- tax_glom(important_AMR_ra, taxrank = "Group")
write.csv(otu_table(important_AMR_ra_group),"important_AMR_GroupRA_otus.csv")
write.csv(tax_table(important_AMR_ra_group),"important_AMR_GroupRA_taxa.csv")

# want to add an overall 'taxonomy' so I can get total counts for this
# going to call the 'taxonomy' Total and they will all be 'AMR'
Total <- c("AMR","AMR","AMR","AMR","AMR","AMR","AMR","AMR",
           "AMR","AMR","AMR","AMR","AMR","AMR","AMR","AMR",
           "AMR","AMR","AMR","AMR","AMR","AMR","AMR","AMR",
           "AMR","AMR","AMR","AMR","AMR","AMR","AMR","AMR",
           "AMR","AMR","AMR","AMR","AMR","AMR","AMR","AMR",
           "AMR","AMR","AMR","AMR","AMR","AMR","AMR","AMR",
           "AMR","AMR","AMR","AMR","AMR","AMR","AMR","AMR")

tax_table(important_AMR_ra) <- cbind(Total=Total, tax_table(important_AMR_ra))

important_AMR_ra_melt_group <- important_AMR_ra %>%
  tax_glom(taxrank = "Group") %>%
  psmelt()

important_AMR_ra_melt_total <- important_AMR_ra %>%
  tax_glom(taxrank = "Total") %>%
  psmelt()

## want to get the mean values for each of the feces
importantTotal_cow <- important_AMR_ra_melt_total[which(important_AMR_ra_melt_total$species_matrix=="Cow Feces"),]
mean(importantTotal_cow$Abundance) # 0.78
sd(importantTotal_cow$Abundance)/sqrt(8) # SEM 0.11

importantTotal_bison <- important_AMR_ra_melt_total[which(important_AMR_ra_melt_total$species_matrix=="Bison Feces"),]
mean(importantTotal_bison$Abundance) # 0.78
sd(importantTotal_bison$Abundance)/sqrt(8) # SEM 0.11

importantTotal_elk <- important_AMR_ra_melt_total[which(important_AMR_ra_melt_total$species_matrix=="Elk Feces"),]
mean(importantTotal_elk$Abundance) # 0.78
sd(importantTotal_elk$Abundance)/sqrt(16) # SEM 0.11

importantTotal_bison_soil <- important_AMR_ra_melt_total[which(important_AMR_ra_melt_total$species_matrix=="Bison Soil"),]
mean(importantTotal_bison_soil$Abundance) # 0.78
sd(importantTotal_bison_soil$Abundance)/sqrt(8) # SEM 0.11

importantTotal_elk_soil <- important_AMR_ra_melt_total[which(important_AMR_ra_melt_total$species_matrix=="Elk Soil"),]
mean(importantTotal_elk_soil$Abundance) # 0.78
sd(importantTotal_elk_soil$Abundance)/sqrt(16) # SEM 0.11

#important genes palette
importantAMR_palette <- c("dodgerblue2","firebrick4","indianred3","chocolate4","darkorange3","#7DE56E","brown3","goldenrod2")

ggplot(important_AMR_ra_melt_group, aes(x= species, y= Abundance)) + 
  theme_bw() +
  labs(y="Relative Abundance (%)") +
  facet_wrap(~matrix, ncol=2, labeller = labeller(matrix = matrix.labs)) +
  geom_bar(aes(fill = Group), stat = "summary", colour = "black", size = 0.75) +
  geom_errorbar(data = important_AMR_ra_melt_total, aes(x= species, y= Abundance), stat = "summary", width = 0.45, size = 0.75) +
  scale_x_discrete(limits = c("Cow","Bison","Elk")) +
  scale_fill_manual(values = importantAMR_palette) +
  scale_y_continuous(expand = c(0.002,0,0.13,0)) +
  theme(#legend.position = "none",
    plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
    panel.border = element_rect(size = 1.5, colour = "black"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size =28, colour = "black"),
    axis.text = element_text(size = 14, colour = "black"),
    axis.title.y = element_text(size = 40, vjust = 1.75),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 28, colour = "black", angle = 45, hjust = 0.99, vjust = 0.99),
    axis.ticks = element_line(size = 0.9, colour = "black"))

## need a new species_matrix column for stats
important_AMR_ra_melt_total$species_matrix <- paste(important_AMR_ra_melt_total$species, important_AMR_ra_melt_total$matrix)

pairwise.wilcox.test(important_AMR_ra_melt_total$Abundance, important_AMR_ra_melt_total$species_matrix, p.adjust.method = "BH")
  
## diversity and richness of medically important AMR genes
importantAMR_alpha <- subset_taxa(data_noSNP, Group=="CTX" |
                                    Group=="GES" |
                                    Group=="IMI" |
                                    Group=="KPC" |
                                    Group=="SHV" |
                                    Group=="TEM" |
                                    Group=="IMP" |
                                    Group=="NDM" |
                                    Group=="CMY" |
                                    Group=="OXA" |
                                    Group=="MEC" |
                                    Group=="MCR" |
                                    Group=="VAT" |
                                    Group=="VGB" |
                                    Group=="CFR" |
                                    Group=="VGA" |
                                    Group=="SME")

## to make a table of counts
write.csv(otu_table(importantAMR_alpha), "importantAMR_otus.csv")
write.csv(tax_table(importantAMR_alpha), "importantAMR_taxa.csv")

importantAMR_alpha1 <- estimate_richness(importantAMR_alpha, measures = c("Observed","Shannon"))
importantAMR_alpha.df <- as(sample_data(importantAMR_alpha), "data.frame")
important_alpha_meta <- cbind(importantAMR_alpha1, importantAMR_alpha.df)
important_alpha_meta

#richness
ggplot(important_alpha_meta, aes(x= species, y = Observed, fill = species, colour = species)) + 
  theme_bw() + facet_wrap(~matrix, ncol=1, labeller = labeller(matrix = matrix.labs)) +
  coord_flip() +
  labs(title= "", y= "No. unique ARGs", x= "") +
  #geom_bar(stat = "summary", alpha = 0.5, size = 0.75) +
  #geom_errorbar(stat = "summary", width = 0.45, size = 0.75) +
  geom_boxplot(alpha = 0.5, size = 1) +
  geom_point(size = 2) +
  scale_fill_manual(values = four_colour_palette) +
  scale_colour_manual(values = four_colour_palette) +
  scale_y_continuous(expand = c(0,0,0.2,0)) +
  scale_x_discrete(limits = c("Elk","Bison","Cow")) +
  theme(legend.position = "none",
        plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
        strip.background = element_blank(),
        strip.text = element_text(size =28, colour = "black", hjust = -0.005),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 16, colour = "black"),
        #axis.text.x = element_blank(),
        axis.title.x = element_text(size = 40, vjust = -0.5),
        axis.title.y = element_blank(),
        axis.ticks.x = element_line(colour = "black"),
        axis.ticks.y = element_blank(),
        plot.title = element_text(),
        panel.border = element_rect(colour = "black", size = 1.5),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank())

# need a new column for species_matrix
important_alpha_meta$species_matrix <- paste(important_alpha_meta$species, important_alpha_meta$matrix)

pairwise.wilcox.test(important_alpha_meta$Observed, important_alpha_meta$species_matrix, p.adjust.method = "BH")
pairwise.wilcox.test(important_alpha_meta$Shannon, important_alpha_meta$species_matrix, p.adjust.method = "BH")

#Shannon
ggplot(important_alpha_meta, aes(x= species, y = Shannon, fill = species, colour = species)) + 
  theme_bw() + facet_wrap(~matrix, ncol=1, labeller = labeller(matrix = matrix.labs)) +
  coord_flip() +
  labs(title= "", y= "Shannon", x= "") +
  #geom_bar(stat = "summary", alpha = 0.5, size = 0.75) +
  #geom_errorbar(stat = "summary", width = 0.45, size = 0.75) +
  geom_boxplot(alpha = 0.5, size = 0.75) +
  geom_point(size = 2) +
  scale_fill_manual(values = four_colour_palette) +
  scale_colour_manual(values = four_colour_palette) +
  scale_y_continuous(expand = c(0,0,0.2,0)) +
  scale_x_discrete(limits = c("Elk","Bison","Cow")) +
  theme(#legend.position = "none",
        plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
        strip.background = element_blank(),
        strip.text = element_text(size =28, colour = "black", hjust = -0.005),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 16, colour = "black"),
        #axis.text.x = element_blank(),
        axis.title.x = element_text(size = 40, vjust = -0.5),
        axis.title.y = element_blank(),
        axis.ticks.x = element_line(colour = "black"),
        axis.ticks.y = element_blank(),
        plot.title = element_text(),
        panel.border = element_rect(colour = "black", size = 1.5),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank())

## ordination of medically important genes
importantAMR.css <- subset_taxa(data_noSNP.css, Group=="CTX" |
                                    Group=="GES" |
                                    Group=="IMI" |
                                    Group=="KPC" |
                                    Group=="SHV" |
                                    Group=="TEM" |
                                    Group=="IMP" |
                                    Group=="NDM" |
                                    Group=="CMY" |
                                    Group=="OXA" |
                                    Group=="MEC" |
                                    Group=="MCR" |
                                    Group=="VAT" |
                                    Group=="VGB" |
                                    Group=="CFR" |
                                    Group=="VGA" |
                                    Group=="SME")

CPH <- subset_ta

important_ARGs_group <- tax_glom(important)

write.csv(otu_table(importantAMR.css),"importantAMR_otus.csv")
write.csv(tax_table(importantAMR.css),"importantAMR_taxa.csv")

sum(taxa_sums(importantAMR.css)==0) # 0
sum(sample_sums(importantAMR.css)==0) # 19
importantAMR.css <- prune_samples(sample_sums(importantAMR.css) > 0, importantAMR.css)
sum(sample_sums(importantAMR.css)==0) # 0

## need to add the species_matrix category to sample_data
important_alpha_meta$species_matrix <- paste(important_alpha_meta$species, important_alpha_meta$matrix)

sample_data(importantAMR.css)$species_matrix <- paste(sample_data(importantAMR.css)$species, sample_data(importantAMR.css)$matrix)

importantAMR.css.df <- as(sample_data(importantAMR.css),"data.frame")

# distance
importantAMR.dist <- vegdist(t(otu_table(importantAMR.css)), method = "bray")

# ORDINATE
importantAMR.ord <- vegan::metaMDS(comm = t(otu_table(importantAMR.css)), distance = "bray", try = 20, trymax = 500, autotransform = F)

plot_ordination(importantAMR.css, importantAMR.ord, type = "samples", color = "species_matrix") +
  theme_bw() +
  geom_point(size = 4) +
  facet_wrap(~matrix, ncol=2, labeller = labeller(matrix = matrix.labs)) +
  #geom_mark_hull(aes(fill=species_matrix), concavity = 5, expand = 0, radius = 0, size =0.5, alpha = 0.35) +
  stat_ellipse(geom = "polygon", alpha = 0.5, aes(fill = species_matrix), lty = 4, lwd = 1) +
  scale_fill_manual(values = c("#638652","#638652","#4d5660","#8d391e","#8d391e")) +
  scale_colour_manual(values = c("#638652","#638652","#4d5660","#8d391e","#8d391e")) +
  theme(legend.position = "none",
        plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
        panel.border = element_rect(colour = "black", size =1.7),
        strip.background = element_blank(),
        strip.text = element_text(size =30, colour = "black"),
        axis.text = element_text(size = 18, colour = "black"),
        axis.title.y = element_text(size = 30, vjust = 1.75),
        axis.title.x = element_text(size = 30, vjust = -1.5))

# permanova
pairwise.adonis(importantAMR.dist, importantAMR.css.df$species_matrix, perm = 9999, p.adjust.m = "BH")
importantAMR.disper <- betadisper(importantAMR.dist, importantAMR.css.df$species_matrix)
plot(importantAMR.disper)
importantAMR.permdisp <- permutest(importantAMR.disper, permutations = 9999, pairwise = T)
importantAMR.permdisp

### rel abund stuff
rel_abund <- transform_sample_counts(data_noSNP.css, function(x) {x/sum(x)} * 100)

rel_abund_class <- tax_glom(rel_abund, taxrank = "Class", NArm = F)
rel_abund_class

### we need some kinda core/pan resistome thing for feces
core_ra <- transform(feces_noSNPS.css,"compositional")
core_ra_group <- aggregate_taxa(core_ra,"Group")
core_ra_group_filt <- core(core_ra_group, detection= 0.01, prevalence =0)
core_ra_mechanism <- aggregate_taxa(core_ra,"Mechanism")
core_ra_mechanism_filt <- core(core_ra_mechanism, detection= 0.01, prevalence =0)

# core with compositionals
prevalences <- seq(0.1,1,0.1)
detections <- 10^seq(log10(2e-3), log10(0.1), length = 30)

p1 <- plot_core(core_ra_group,
                plot.type = "heatmap",
                colours = rev(brewer.pal(10, "Spectral")),
                prevalences = prevalences,
                detections = detections, min.prevalence = .90)
p1 + theme_bw() +
  labs(y= "", x= "Detection Threshold \n (Relative Abundance (%))", title = "") +
  theme(legend.key.size = unit(4, "lines"),
        plot.title = element_text(size = 24),
        plot.margin = unit(c(1,1,1,1), "cm"),
        panel.border = element_rect(colour = "black", size = 1),
        axis.ticks = element_line(size = 0.75, colour = "black"),
        axis.text.y = element_text(colour = "black", size = 26),
        axis.text.x = element_text(colour = "black", size = 24),
        axis.title.x = element_text(size = 40, vjust = -0.75))

### RA plot for feces only
ra_feces <- transform_sample_counts(feces_noSNPS.css, function(x) {x/sum(x)} * 100)
sum(taxa_sums(ra_feces)==0) # none, good

ra_feces_class <- tax_glom(ra_feces, taxrank = "Class", NArm = F) %>%
  psmelt() # 36 classes
ra_feces_class_palette <- distinctColorPalette(36)
ra_feces_mech <- tax_glom(ra_feces, taxrank = "Mechanism", NArm = F) %>%
  psmelt() # 97 mech

ra_feces_group <- tax_glom(ra_feces, taxrank = "Group", NArm = F) %>%
  psmelt() #393 groups, toooooo many

ggplot(ra_feces_class, aes(x= group, y = Abundance, fill = Class)) +
  theme_bw() +
  geom_bar(stat = "summary", colour = "black") +
  scale_fill_manual(values = ra_feces_class_palette)

### want to look at tetracyclines, drug_biocide, aminoglycosides, MLS

## tetracyclines
feces_tetracyclines <- subset_taxa(ra_feces, Class=="Tetracyclines")

feces_tetra_class <- tax_glom(feces_tetracyclines, taxrank = "Class", NArm = F) %>%
  psmelt()
feces_tetra_group <- tax_glom(feces_tetracyclines, taxrank = "Group", NArm = F) %>%
  psmelt() # 25 groups

feces_tetra_group_export <- tax_glom(feces_tetracyclines, taxrank = "Group", NArm = F)
write.csv(otu_table(feces_tetra_group_export),"feces_tetra_groupRA_otus.csv")
write.csv(tax_table(feces_tetra_group_export),"feces_tetra_groupRA_taxa.csv")

feces_tetra_group_palette <- c("#F0FFFF","#89CFF0","#0000FF","#7393B3","#088F8F",
                               "#0096FF","#0047AB","#6495ED","#00FFFF","#00008B",
                               "#1434A4","#7DF9FF","#6082B6","#5D3FD3","#191970",
                               "#1F51FF","#CCCCFF","#B6D0E2","#96DED1","#87CEEB",
                               "#4682B4","#0F52BA","#008080","#40B5AD","#0818A8")
feces_tetra_mech <- tax_glom(feces_tetracyclines, taxrank = "Mechanism", NArm = F) %>%
  psmelt() # 3

ggplot(feces_tetra_group, aes(x= group, y= Abundance)) + 
  theme_bw() +
  labs(y= "Relative Abundance (%)") +
  geom_bar(aes(fill= Group), stat = "summary", colour = "black") +
  geom_errorbar(data = feces_tetra_class, aes(x= group, y= Abundance), stat = "summary", width = 0.55, size = 0.6) +
  scale_fill_manual(values = feces_tetra_group_palette) +
  scale_x_discrete(limits = c("Farm Cow feces", "YNP Bison feces","YNP Elk feces","RMNP Elk feces"),
                   labels = c("Cow; SW USA","Bison; YNP","Elk; YNP","Elk; RMNP")) +
  scale_y_continuous(expand = c(0.002,0,0.1,0)) +
  theme(#legend.position = "none",
    plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
    panel.border = element_rect(size = 1.5, colour = "black"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text = element_text(size = 14, colour = "black"),
    axis.title.y = element_text(size = 28, vjust = 1.5),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 28, colour = "black", angle = 45, hjust= 0.95, vjust = 0.95),
    axis.ticks = element_line(size = 0.9, colour = "black"))

pairwise.wilcox.test(feces_tetra_class$Abundance, feces_tetra_class$group, p.adjust.method = "BH")
# NS

## b-lactams
feces_betalactams <- subset_taxa(ra_feces, Class=="betalactams")

feces_betalactams_class <- tax_glom(feces_betalactams, taxrank = "Class", NArm = F) %>%
  psmelt()

feces_betalactams_group <- tax_glom(feces_betalactams, taxrank = "Group", NArm = F) %>%
  psmelt() # 22 groups

feces_betalactams_group_palette <- c("#F8DE7E","#FADA5E","#F9A602","#FFD300","#D2B55B","#C3B091",
                                     "#DAA520","#FCF4A3","#FCD12A","#FFC30B","#C49102","#FCE205",
                                     "#FDA50F","#CC7722","#FFBF00","#EEDC82","#FEDC56","#FFFDD0","#F5F5DC",
                                     "#EFFDF5","#F8E473","#CED180")

ggplot(feces_betalactams_group, aes(x= group, y= Abundance)) +
  theme_bw() +
  labs(y= "Relative Abundance (%)") +
  geom_bar(aes(fill= Group), stat = "summary", colour = "black") +
  geom_errorbar(data = feces_betalactams_class, aes(x= group, y= Abundance), stat = "summary", width = 0.55, size = 0.6) +
  scale_fill_manual(values = feces_betalactams_group_palette) +
  scale_y_continuous(expand = c(0.002,0,0.1,0)) +
  scale_x_discrete(limits = c("Farm Cow feces", "YNP Bison feces","YNP Elk feces","RMNP Elk feces"),
                   labels = c("Cow; SW USA","Bison; YNP","Elk; YNP","Elk; RMNP")) +
  theme(#legend.position = "none",
    plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
    panel.border = element_rect(size = 1.5, colour = "black"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text = element_text(size = 14, colour = "black"),
    axis.title.y = element_text(size = 28, vjust = 1.5),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 28, colour = "black", angle = 45, hjust= 0.95, vjust = 0.95),
    axis.ticks = element_line(size = 0.9, colour = "black"))

pairwise.wilcox.test(feces_betalactams_class$Abundance, feces_betalactams_class$group, p.adjust.method = "BH")

## aminoglycosides
feces_aminoglycosides <- subset_taxa(ra_feces, Class=="Aminoglycosides")

feces_amino_class <- tax_glom(feces_aminoglycosides, taxrank = "Class", NArm = F) %>%
  psmelt()
feces_amino_group <- tax_glom(feces_aminoglycosides, taxrank = "Group", NArm = F) %>%
  psmelt() # 14 groups
feces_amino_group_export <- tax_glom(feces_aminoglycosides, taxrank = "Group", NArm = F)
write.csv(otu_table(feces_amino_group_export),"feces_amino_group_otus.csv")
write.csv(tax_table(feces_amino_group_export),"feces_amino_group_taxa.csv")

feces_amino_group_palette <- c("#EADDCA","#E1C16E","#CD7F32","#A52A2A","#4A0404",
                               "#800020","#E97451","#6E260E","#C19A6B","#954535",
                               "#D27D2D","#6F4E37","#8B0000","#E5AA70")

ggplot(feces_amino_group, aes(x= group, y= Abundance)) + 
  theme_bw() +
  labs(y= "Relative Abundance (%)") +
  geom_bar(aes(fill= Group), stat = "summary", colour = "black") +
  geom_errorbar(data = feces_amino_class, aes(x= group, y= Abundance), stat = "summary", width = 0.55, size = 0.6) +
  scale_fill_manual(values = feces_amino_group_palette) +
  scale_x_discrete(limits = c("Farm Cow feces", "YNP Bison feces","YNP Elk feces","RMNP Elk feces"),
                   labels = c("Cow; SW USA","Bison; YNP","Elk; YNP","Elk; RMNP")) +
  scale_y_continuous(expand = c(0.002,0,0.1,0)) +
  theme(#legend.position = "none",
        plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
        panel.border = element_rect(size = 1.5, colour = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text = element_text(size = 14, colour = "black"),
        axis.title.y = element_text(size = 28, vjust = 1.5),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 28, colour = "black", angle = 45, hjust= 0.95, vjust = 0.95),
        axis.ticks = element_line(size = 0.9, colour = "black"))

pairwise.wilcox.test(feces_amino_class$Abundance, feces_amino_class$group, p.adjust.method = "BH")
# everything but bison and elk

## drug_biocide
feces_drugbiocide <- subset_taxa(ra_feces, Class=="Drug_and_biocide_resistance")

feces_drugbio_class <- tax_glom(feces_drugbiocide, taxrank = "Class", NArm = F) %>%
  psmelt()
feces_drugbio_group <- tax_glom(feces_drugbiocide, taxrank = "Group", NArm = F) %>%
  psmelt() # 65 groups
feces_drugbio_group_export <- tax_glom(feces_drugbiocide, taxrank = "Group", NArm = F)
write.csv(otu_table(feces_drugbio_group_export),"drugbio_groupRA_otus.csv")
write.csv(tax_table(feces_drugbio_group_export),"drugbio_groupRA_taxa.csv")
feces_drugbio_mech_export <- tax_glom(feces_drugbiocide, taxrank = "Mechanism", NArm = F)
write.csv(otu_table(feces_drugbio_mech_export),"drugbio_mechRA_otus.csv")
write.csv(tax_table(feces_drugbio_mech_export),"drugbio_mechRA_taxa.csv")


feces_drugbio_mech <- tax_glom(feces_drugbiocide, taxrank = "Mechanism", NArm = F) %>%
  psmelt() # 8 mechanisms
feces_drugbio_mech_palette <- c("#9F2B68","#DA70D6","#F88379","#F89880",
                                "#AA336A","#FF69B4","#FFB6C1","#FAA0A0")

ggplot(feces_drugbio_mech, aes(x= group, y= Abundance)) + 
  theme_bw() +
  labs(y= "Relative Abundance (%)") +
  geom_bar(aes(fill= Mechanism), stat = "summary", colour = "black") +
  geom_errorbar(data = feces_drugbio_class, aes(x= group, y= Abundance), stat = "summary", width = 0.55, size = 0.6) +
  scale_fill_manual(values = feces_drugbio_mech_palette) +
  scale_x_discrete(limits = c("Farm Cow feces", "YNP Bison feces","YNP Elk feces","RMNP Elk feces"),
                              labels = c("Cow; SW USA","Bison; YNP","Elk; YNP","Elk; RMNP")) +
  scale_y_continuous(expand = c(0.002,0,0.1,0)) +
  theme(#legend.position = "none",
        plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
        panel.border = element_rect(size = 1.5, colour = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text = element_text(size = 14, colour = "black"),
        axis.title.y = element_text(size = 28, vjust = 1.5),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 28, colour = "black", angle = 45, hjust= 0.95, vjust = 0.95),
        axis.ticks = element_line(size = 0.9, colour = "black"))

pairwise.wilcox.test(feces_drugbio_class$Abundance, feces_drugbio_class$group, p.adjust.method = "BH")
# bison and cow Sig.

## MLS
feces_mls <- subset_taxa(ra_feces, Class=="MLS")

feces_mls_class <- tax_glom(feces_mls, taxrank = "Class", NArm = F) %>%
  psmelt()
feces_mls_class_export <- tax_glom(feces_mls, taxrank = "Class", NArm = F) 
write.csv(otu_table(feces_mls_class_export),"feces_mls_class_otus.csv")

feces_mls_group <- tax_glom(feces_mls, taxrank = "Group", NArm = F) %>%
  psmelt() # 29 groups

feces_mls_group_export <- tax_glom(feces_mls, taxrank = "Group", NArm = F)
write.csv(otu_table(feces_mls_group_export),"feces_MLS_groupRA_otus.csv")
write.csv(tax_table(feces_mls_group_export),"feces_MLS_groupRA_taxa.csv")

feces_drugbio_mech <- tax_glom(feces_mls, taxrank = "Mechanism", NArm = F) %>%
  psmelt() # 8 mechanisms
feces_mls_group_palette <- c("#7FFFD4","#454B1B","#088F8F","#AAFF00","#097969",
                             "#AFE1AF","#32CD32","#00FF7F","#023020","#50C878",
                             "#5F8575","#4F7942","#228B22","#7CFC00","#008000",
                             "#355E3B","#00A36C","#4CBB17","#2AAA8A","#90EE90",
                             "#478778","#32CD32","#0BDA51","#98FB98","#8A9A5B",
                             "#ECFFDC","#C1E1C1","#C9CC3F","#93C572")

ggplot(feces_mls_group, aes(x= group, y= Abundance)) + 
  theme_bw() +
  labs(y= "Relative Abundance (%)") +
  geom_bar(aes(fill= Group), stat = "summary", colour = "black") +
  geom_errorbar(data = feces_mls_class, aes(x= group, y= Abundance), stat = "summary", width = 0.55, size = 0.6) +
  scale_fill_manual(values = feces_mls_group_palette) +
  scale_x_discrete(limits = c("Farm Cow feces", "YNP Bison feces","YNP Elk feces","RMNP Elk feces"),
                   labels = c("Cow; SW USA","Bison; YNP","Elk; YNP","Elk; RMNP")) +
  scale_y_continuous(expand = c(0.002,0,0.1,0)) +
  theme(#legend.position = "none",
        plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
        panel.border = element_rect(size = 1.5, colour = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text = element_text(size = 14, colour = "black"),
        axis.title.y = element_text(size = 28, vjust = 1.5),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 28, colour = "black", angle = 45, hjust= 0.95, vjust = 0.95),
        axis.ticks = element_line(size = 0.9, colour = "black"))

pairwise.wilcox.test(feces_mls_class$Abundance, feces_mls_class$group, p.adjust.method = "BH")
# all sig. diff.

## need filtered to see what is in legend
feces_mls_group_filter <- tax_glom(feces_mls, taxrank = "Group", NArm = F)
feces_mls_group_filter <- filter_taxa(feces_mls_group_filter, function(x) mean(x) > 0.1, T) %>%
  psmelt() # 29 groups

ggplot(feces_mls_group_filter, aes(x= group, y= Abundance)) + 
    theme_bw() +
    labs(y= "Relative Abundance (%)") +
    geom_bar(aes(fill= Group), stat = "summary", colour = "black") +
    #geom_errorbar(data = feces_mls_class, aes(x= group, y= Abundance), stat = "summary", width = 0.55, size = 0.6) +
    scale_fill_manual(values = feces_mls_group_palette) +
    scale_x_discrete(limits = c("Farm Cow feces", "YNP Bison feces","YNP Elk feces","RMNP Elk feces"),
                     labels = c("Cow; SW USA","Bison; YNP","Elk; YNP","Elk; RMNP")) +
    scale_y_continuous(expand = c(0.002,0,0.1,0)) +
    theme(#legend.position = "none",
      plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
      panel.border = element_rect(size = 1.5, colour = "black"),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.text = element_text(size = 14, colour = "black"),
      axis.title.y = element_text(size = 28, vjust = 1.5),
      axis.title.x = element_blank(),
      axis.text.x = element_text(size = 28, colour = "black", angle = 45, hjust= 0.95, vjust = 0.95),
      axis.ticks = element_line(size = 0.9, colour = "black"))

feces_tetracyclines_class <- tax_glom(feces_tetracyclines, taxrank = "Class") %>%
  psmelt()

write.csv(feces_tetracyclines_class, "feces_tetracyclineRA.csv")  

### correlation between S24-7 and MLS resistance test
S247MLS.df <- read.table("S247-MLS_RAs.txt", header = T, row.names = 1)
S247MLS.df

S247MLS.cor <- cor.test(S247MLS.df$S24.7_RA, S247MLS.df$MLS_RA, method = "pearson")
S247MLS.cor  

ggplot(S247MLS.df, aes(x= ERM_all, y= S24.7_RA)) +
  geom_point()

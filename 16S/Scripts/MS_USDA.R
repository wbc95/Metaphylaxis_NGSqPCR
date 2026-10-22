#### setwd ####
setwd("~/Documents/Research/Projects/USDA_NIFA_AMR/Sequencing/16S")

#### load libraries ####
library(phyloseq); library(ggplot2); library(btools); library(dplyr);
library(tidyr); library(stringr); library(randomcoloR); library(metagMisc);
library(metagenomeSeq); library(GUniFrac); library(randomForest); library(pairwiseAdonis)
library(knitr); library(readr); library(kableExtra); library(scales); library(vegan);
library(ggdendro); library(forcats); library(ggbreak); library(BiodiversityR); library(ggrepel); library(ggforce);
library(ggprism); library(rstatix); library(ggpubr)

#### source some functions ####
source("~/Documents/Research/Projects/USDA_NIFA_AMR/Pooling Pilot/16S Pooling/changeGGTaxaNames.R")
source("~/Documents/Research/Projects/USDA_NIFA_AMR/Pooling Pilot/16S Pooling/w_unifrac.R")
source("~/Documents/Research/Projects/USDA_NIFA_AMR/Pooling Pilot/16S Pooling/g_unifrac.R")
source("~/Documents/Research/Projects/USDA_NIFA_AMR/Pooling Pilot/16S Pooling/uw_unifrac.R")
#source("/Users/ljpinnell/Documents/R_Functions_Scripts/MergeLowAbund.R")

#### IMPORTING DATA ####
qiimedata <- import_biom("table-with-taxonomy.biom","tree.nwk","dna-sequences.fasta")
qiimedata # 181978 taxa and 313 samples
write.csv(sample_names(qiimedata),"samplenamesformetadata.csv")

map_file <- import_qiime_sample_data("metadata_16S.txt")

View(sample_data(data))

# combining data with metadata
data <- merge_phyloseq(qiimedata, map_file)
data # 181978 taxa and 99 samples

sample_names(data) <- sample_data(data)$ID
sample_names(data)


#### PRE-PROCESSING ####
## renaming ranks
rank_names(data) # "Rank1" - "Rank7" not ideal, lets change em
colnames(tax_table(data)) <- c("Domain","Kingdom","Phylum","Class","Order","Family","Genus")
rank_names(data) # beauty, now they are named properly

# changing the SILVA style naming (i.e., getting rid of the 'k' in k__Bacteria, etc.)
tax.data <- data.frame(tax_table(data)) # extract the taxonomy table as a data frame
tax.data.names <- changeGGtaxa(tax.data) # this gets rid of the GG format

# now to change the NAs to a better naming scheme by putting in the highest assigned taxonomy
for (i in 1:7){ tax.data.names[,i] <- as.character(tax.data.names[,i])} # converting all columns to characters
tax.data.names[is.na(tax.data.names)] <- "" # replacing the NAs with an empty string

# now filling in the empty slots with the highest assigned taxonomy
for (i in 1:nrow(tax.data.names)){
  if (tax.data.names[i,2] == ""){
    domain <- paste("unclassified ", tax.data.names[i,1], sep = "")
    tax.data.names[i, 2:7] <- domain
  } else if (tax.data.names[i,3] == ""){
    kingdom <- paste("unclassified ", tax.data.names[i,2], sep = "")
    tax.data.names[i, 3:7] <- kingdom
  } else if (tax.data.names[i,4] == ""){
    phylum <- paste("unclassified ", tax.data.names[i,3], sep = "")
    tax.data.names[i, 4:7] <- phylum
  } else if (tax.data.names[i,5] == ""){
    class <- paste("unclassified ", tax.data.names[i,4], sep = "")
    tax.data.names[i, 5:7] <- class
  } else if (tax.data.names[i,6] == ""){
    order <- paste("unclassified ", tax.data.names[i,5], sep = "")
    tax.data.names[i, 6:7] <- order
  } else if (tax.data.names[i,7] == ""){
    tax.data.names$Genus[i] <- paste("unclassified ",tax.data.names$Family[i], sep = "")
  }
}

head(tax.data.names) # great, no more NAs and no more k__
tax_table(data) <- as.matrix(tax.data.names) # re-insert the taxonomy table into the phyloseq object
head(tax_table(data), 20) # sweet, lookin good!

# removing Eukaryota
#data <- subset_taxa(data, Domain!="Eukaryota")
#data # lost two taxa

# some QC checks
min(sample_sums(data)) #286
max(sample_sums(data)) # 3858775
mean(sample_sums(data)) # 838797
median(sample_sums(data)) # 890796
sort(sample_sums(data)) # ...2872,2874,6546,13022, 69352 .... going to cutoff at 15K

data <- prune_samples(sample_sums(data) > 15000, data)
data # left with 83 samples
data <- prune_taxa(taxa_sums(data) > 0, data)
any(taxa_sums(data)==0)
data # 64293 taxa and 83 samples remain

sample_data(data)

# splitting samples and controls
samples <- subset_samples(data, Type=="Sample")
samples # 64293 and 82 samples, so one of the control samples made it through pre-processing...problem?

min(sample_sums(samples)) #69352
max(sample_sums(samples)) #3858775
mean(sample_sums(samples)) #1000944

any(taxa_sums(samples)==0)
samples <- prune_taxa(taxa_sums(samples)>0, samples)


########################################### COMPARING SEQ DEPTH #####
sample_sum_df_new <- data.frame(ASV_count = sample_sums(samples))
metadata.df <- as(sample_data(samples),"data.frame")
SeqDepth_metadata <- cbind(metadata.df, sample_sum_df_new)

# comparing treatment groups
wilcox_test(SeqDepth_metadata, ASV_count~TxGroup) #p=0.01075, difference
ggplot(SeqDepth_metadata, aes(x=TxGroup,y=ASV_count))+
  geom_boxplot()+
  facet_wrap(~BRD, ncol = 3)+
  stat_compare_means()



########################################### ALPHA DIVERSITY #####
alpha_div1 <- estimate_richness(samples, measures = c("Observed", "Shannon", "Simpson","InvSimpson")) # richness, div
alpha_div2 <- estimate_pd(samples) # faith's pd

# combine metrics with metadata
alpha_div <- cbind(alpha_div1, alpha_div2)
alpha_div
alpha_div <- alpha_div[,c(1:5)]
alpha_div
alpha_div_meta <- cbind(metadata.df, alpha_div)
alpha_div_meta # metadata and div metrics

#comparisons for brd
brdcomparisons <- list(c("LATE","NO"),c("LATE","YES"),c("NO","YES"))

#### RICHNESS FIGURES ####
#Treatment group
ggplot(alpha_div_meta, aes(x=TxGroup,y=Observed, fill = TxGroup, colour = TxGroup)) +
  theme_bw() +
  labs(y= "Observed ASVs") +
  facet_wrap(~Time, ncol = 2, scales = "free") +
  #scale_y_continuous(limits = c(250,800), breaks = c(250,500,750)) +
  geom_boxplot(alpha = 0.5, linewidth = 0.75) +
  geom_point(shape = 18, size = 4) +
  #eom_bar(stat = "summary") +
  #geom_errorbar(stat = "summary") +
  #scale_fill_manual(values = tylosin_palette) +
  #scale_colour_manual(values = tylosin_palette) +
  theme(#legend.position = "none",
        panel.border = element_rect(colour = "black", linewidth = 1),
        panel.grid.major.x = element_blank(),
        strip.background = element_blank(),
        #strip.text = element_blank(),
        axis.title.y = element_text(size = 24),
        axis.text.y = element_text(colour = "black", size = 14),
        axis.ticks.y = element_line(linewidth = 0.75, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())+
  stat_compare_means(method = "wilcox")
wilcox_test(alpha_div_meta, Observed~TxGroup) #NS

#BRD
ggplot(alpha_div_meta, aes(x=BRD,y=Observed, fill = BRD, colour = BRD)) +
  theme_bw() +
  labs(y= "Observed ASVs") +
  facet_wrap(~TxGroup+Time, scales = "free") +
  #scale_y_continuous(limits = c(250,800), breaks = c(250,500,750)) +
  geom_boxplot(alpha = 0.5, linewidth = 0.75) +
  geom_point(shape = 18, size = 4) +
  #eom_bar(stat = "summary") +
  #geom_errorbar(stat = "summary") +
  #scale_fill_manual(values = tylosin_palette) +
  #scale_colour_manual(values = tylosin_palette) +
  theme(#legend.position = "none",
    panel.border = element_rect(colour = "black", linewidth = 1),
    panel.grid.major.x = element_blank(),
    strip.background = element_blank(),
    #strip.text = element_blank(),
    axis.title.y = element_text(size = 24),
    axis.text.y = element_text(colour = "black", size = 14),
    axis.ticks.y = element_line(linewidth = 0.75, colour = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank())+
  stat_compare_means()
kruskal_test(alpha_div_meta, Observed~BRD) #NS



#### FPD FIGURES ####
ggplot(alpha_div_meta, aes(x=TxGroup,y=PD, fill = TxGroup, colour = TxGroup)) +
  theme_bw() +
  labs(y= "Faith's PD") +
  facet_wrap(~Time, ncol = 2, scales = "free") +
  #scale_y_continuous(limits = c(250,800), breaks = c(250,500,750)) +
  geom_boxplot(alpha = 0.5, linewidth = 0.75) +
  geom_point(shape = 18, size = 4) +
  #eom_bar(stat = "summary") +
  #geom_errorbar(stat = "summary") +
  #scale_fill_manual(values = tylosin_palette) +
  #scale_colour_manual(values = tylosin_palette) +
  theme(#legend.position = "none",
    panel.border = element_rect(colour = "black", linewidth = 1),
    panel.grid.major.x = element_blank(),
    strip.background = element_blank(),
    #strip.text = element_blank(),
    axis.title.y = element_text(size = 24),
    axis.text.y = element_text(colour = "black", size = 14),
    axis.ticks.y = element_line(linewidth = 0.75, colour = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank())+
  stat_compare_means(method = "wilcox")
wilcox_test(alpha_div_meta, PD~TxGroup) #NS

#BRD
ggplot(alpha_div_meta, aes(x=BRD,y=PD, fill = BRD, colour = BRD)) +
  theme_bw() +
  labs(y= "Faith's PD") +
  facet_wrap(~TxGroup+Time, scales = "free") +
  #scale_y_continuous(limits = c(250,800), breaks = c(250,500,750)) +
  geom_boxplot(alpha = 0.5, linewidth = 0.75) +
  geom_point(shape = 18, size = 4) +
  #eom_bar(stat = "summary") +
  #geom_errorbar(stat = "summary") +
  #scale_fill_manual(values = tylosin_palette) +
  #scale_colour_manual(values = tylosin_palette) +
  theme(#legend.position = "none",
    panel.border = element_rect(colour = "black", linewidth = 1),
    panel.grid.major.x = element_blank(),
    strip.background = element_blank(),
    #strip.text = element_blank(),
    axis.title.y = element_text(size = 24),
    axis.text.y = element_text(colour = "black", size = 14),
    axis.ticks.y = element_line(linewidth = 0.75, colour = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank())+
  stat_compare_means()+
  stat_compare_means(comparisons=brdcomparisons)
kruskal_test(alpha_div_meta, PD~BRD) #NS

########################################### BETA DIVERSITY #####

#### CSS TRANSFORM ################
any(taxa_sums(samples)==0) # triple checking, none - good!
samples.css <- phyloseq_transform_css(samples, log = F)

#keeping for code skeleton if need to split out later
#### SPLITTING CSS COUNTS BY BODY SITE #######
### All LIVER ABSCESSES
#liv_abs.css <- subset_samples(samples.css, body_location=="liver abscess")
#liv_abs.css <- prune_taxa(taxa_sums(liv_abs.css) > 0, liv_abs.css)
#liv_abs.css # 131 taxa and 53 samples
#liv_abs.css.df <- as(sample_data(liv_abs.css),"data.frame")



#################################### REL ABUND CALCULATIONS ######
rel_abund <- transform_sample_counts(samples.css, function(x) {x/sum(x)} * 100)

ra_phylum <- tax_glom(rel_abund, taxrank = "Phylum", NArm = F)
ra_phylum #60 taxa
write.csv(tax_table(ra_phylum),"phylum_taxa.csv"); write.csv(otu_table(ra_phylum),"pylum_otus.csv")

#### ALL LIVER ABSCESSES
#ra_LAs <- subset_samples(rel_abund, body_location=="liver abscess")
#ra_LAs <- prune_taxa(taxa_sums(ra_LAs) > 0, ra_LAs)
#ra_LAs #131 taxa, 53 samples
## AGGLOM.
LAs_phylum <- tax_glom(ra_LAs, taxrank = "Phylum", NArm = F)
#write.csv(tax_table(LAs_phylum),"LA_Phylum_taxa.csv"); write.csv(otu_table(LAs_phylum),"LA_Phylum_otus.csv")
LAs_class <- tax_glom(ra_LAs, taxrank = "Class", NArm = F)
#write.csv(tax_table(LAs_class),"LA_class_taxa.csv"); write.csv(otu_table(LAs_class),"LA_class_otus.csv")
LAs_order <- tax_glom(ra_LAs, taxrank = "Order", NArm = F)
#write.csv(tax_table(LAs_order),"LA_order_taxa.csv"); write.csv(otu_table(LAs_order),"LA_order_otus.csv")
LAs_family <- tax_glom(ra_LAs, taxrank = "Family", NArm = F)
#write.csv(tax_table(LAs_family),"LA_family_taxa.csv"); write.csv(otu_table(LAs_family),"LA_family_otus.csv")
LAs_genus <- tax_glom(ra_LAs, taxrank = "Genus", NArm = F)
#write.csv(tax_table(LAs_genus),"LA_genus_taxa.csv"); write.csv(otu_table(LAs_genus),"LA_genus_otus.csv")

#### MERGED LAs
mergedLAs_ra <- transform_sample_counts(liver_abscesses_merged.css, function(x) {x/sum(x)} * 100)
mergedLAs_genus <- tax_glom(mergedLAs_ra, taxrank = "Genus", NArm = F)
#write.csv(tax_table(mergedLAs_genus), "LA_genus_taxa.csv")
#write.csv(t(otu_table(mergedLAs_genus)), "LA_genus_otus.csv")

#### RUMEN OVERALL
ra_rumen <- transform_sample_counts(rumen.css,function(x) {x/sum(x)} * 100)
ra_rumen # 5212 taxa, 71 samples
rumen_phylum <- tax_glom(ra_rumen, taxrank = "Phylum", NArm = F)
#write.csv(tax_table(rumen_phylum),"rumen_Phylum_taxa.csv"); write.csv(otu_table(rumen_phylum),"rumen_Phylum_otus.csv")
rumen_class <- tax_glom(ra_rumen, taxrank = "Class", NArm = F)
#write.csv(tax_table(rumen_class),"rumen_class_taxa.csv"); write.csv(otu_table(rumen_class),"rumen_class_otus.csv")
rumen_order <- tax_glom(ra_rumen, taxrank = "Order", NArm = F)
#write.csv(tax_table(rumen_order),"rumen_order_taxa.csv"); write.csv(otu_table(rumen_order),"rumen_order_otus.csv")
rumen_family <- tax_glom(ra_rumen, taxrank = "Family", NArm = F)
#write.csv(tax_table(rumen_family),"rumen_family_taxa.csv"); write.csv(otu_table(rumen_family),"rumen_family_otus.csv")
rumen_genus <- tax_glom(ra_rumen, taxrank = "Genus", NArm = F)
#write.csv(tax_table(rumen_genus),"rumen_genus_taxa.csv"); write.csv(otu_table(rumen_genus),"rumen_genus_otus.csv")

#### RUMEN LUMEN
ra_rumen_lumen <- transform_sample_counts(rumen_lumen.css, function(x) {x/sum(x)} * 100)
rumen_lumen_phylum <- tax_glom(ra_rumen_lumen, taxrank = "Phylum", NArm = F) #21 phyla
rumen_lumen_class <- tax_glom(ra_rumen_lumen, taxrank = "Class", NArm = F) #35 classes
rumen_lumen_order <- tax_glom(ra_rumen_lumen, taxrank = "Order", NArm = F) # 74 orders
rumen_lumen_family <- tax_glom(ra_rumen_lumen, taxrank = "Family", NArm = F) #120 families
#write.csv(otu_table(rumen_lumen_family),"rumen_lumen_otus.csv")
#write.csv(tax_table(rumen_lumen_family),"rumen_lumen_taxa.csv")
rumen_lumen_family_melt <- psmelt(rumen_lumen_family)
rumen_lumen_genus <- tax_glom(ra_rumen_lumen, taxrank = "Genus", NArm = F) # 257 genera
rumen_lumen_genus_melt <- psmelt(rumen_lumen_genus)

#### RUMEN EPITHELIUM
ra_rumen_epithelium <- transform_sample_counts(rumen_epithelium.css, function(x) {x/sum(x)} * 100)
rumen_epithelium_phylum <- tax_glom(ra_rumen_epithelium, taxrank = "Phylum", NArm = F) #41 phyla
rumen_epithelium_class <- tax_glom(ra_rumen_epithelium, taxrank = "Class", NArm = F) #40 classes
rumen_epithelium_order <- tax_glom(ra_rumen_epithelium, taxrank = "Order", NArm = F) # 74 orders
rumen_epithelium_family <- tax_glom(ra_rumen_epithelium, taxrank = "Family", NArm = F) #128 families
#write.csv(otu_table(rumen_epithelium_family),"rumen_epithelium_otus.csv")
#write.csv(tax_table(rumen_epithelium_family),"rumen_epithelium_taxa.csv")
rumen_epithelium_family_melt <- psmelt(rumen_epithelium_family)
rumen_epithelium_genus <- tax_glom(ra_rumen_epithelium, taxrank = "Genus", NArm = F) # 286 genera
rumen_epithelium_genus_melt <- psmelt(rumen_epithelium_genus)

######### SM INTESTINE LUMEN
ra_sm_intestine_lumen <- transform_sample_counts(sm_intestine_lumen.css, function(x) {x/sum(x)} * 100)
sm_intestine_lumen_phylum <- tax_glom(ra_sm_intestine_lumen, taxrank = "Phylum", NArm = F) #21 phyla
sm_intestine_lumen_class <- tax_glom(ra_sm_intestine_lumen, taxrank = "Class", NArm = F) #35 classes
sm_intestine_lumen_order <- tax_glom(ra_sm_intestine_lumen, taxrank = "Order", NArm = F) # 74 orders
sm_intestine_lumen_family <- tax_glom(ra_sm_intestine_lumen, taxrank = "Family", NArm = F) #120 families
#write.csv(otu_table(sm_intestine_lumen_family),"sm_intestine_lumen_otus.csv")
#write.csv(tax_table(sm_intestine_lumen_family),"sm_intestine_lumen_taxa.csv")
sm_intestine_lumen_family_melt <- psmelt(sm_intestine_lumen_family)
sm_intestine_lumen_genus <- tax_glom(ra_sm_intestine_lumen, taxrank = "Genus", NArm = F) # 257 genera
sm_intestine_lumen_genus_melt <- psmelt(sm_intestine_lumen_genus)

#### SM INTESTINE EPITHELIUM
ra_sm_intestine_epithelium <- transform_sample_counts(sm_intestine_epithelium.css, function(x) {x/sum(x)} * 100)
sm_intestine_epithelium_phylum <- tax_glom(ra_sm_intestine_epithelium, taxrank = "Phylum", NArm = F) #41 phyla
sm_intestine_epithelium_class <- tax_glom(ra_sm_intestine_epithelium, taxrank = "Class", NArm = F) #40 classes
sm_intestine_epithelium_order <- tax_glom(ra_sm_intestine_epithelium, taxrank = "Order", NArm = F) # 74 orders
sm_intestine_epithelium_family <- tax_glom(ra_sm_intestine_epithelium, taxrank = "Family", NArm = F) #128 families
#write.csv(otu_table(sm_intestine_epithelium_family),"sm_intestine_epithelium_otus.csv")
#write.csv(tax_table(sm_intestine_epithelium_family),"sm_intestine_epithelium_taxa.csv")
sm_intestine_epithelium_family_melt <- psmelt(sm_intestine_epithelium_family)
sm_intestine_epithelium_genus <- tax_glom(ra_sm_intestine_epithelium, taxrank = "Genus", NArm = F) # 286 genera
sm_intestine_epithelium_genus_melt <- psmelt(sm_intestine_epithelium_genus)

#### LG INTESTINE LUMEN
ra_lg_intestine_lumen <- transform_sample_counts(lg_intestine_lumen.css, function(x) {x/sum(x)} * 100)
lg_intestine_lumen_phylum <- tax_glom(ra_lg_intestine_lumen, taxrank = "Phylum", NArm = F) #21 phyla
lg_intestine_lumen_class <- tax_glom(ra_lg_intestine_lumen, taxrank = "Class", NArm = F) #35 classes
lg_intestine_lumen_order <- tax_glom(ra_lg_intestine_lumen, taxrank = "Order", NArm = F) # 74 orders
lg_intestine_lumen_family <- tax_glom(ra_lg_intestine_lumen, taxrank = "Family", NArm = F) #120 families
#write.csv(otu_table(lg_intestine_lumen_family),"lg_intestine_lumen_otus.csv")
#write.csv(tax_table(lg_intestine_lumen_family),"lg_intestine_lumen_taxa.csv")
lg_intestine_lumen_family_melt <- psmelt(lg_intestine_lumen_family)
lg_intestine_lumen_genus <- tax_glom(ra_lg_intestine_lumen, taxrank = "Genus", NArm = F) # 257 genera
lg_intestine_lumen_genus_melt <- psmelt(lg_intestine_lumen_genus)

#### LG INTESTINE EPITHELIUM
ra_lg_intestine_epithelium <- transform_sample_counts(lg_intestine_epithelium.css, function(x) {x/sum(x)} * 100)
lg_intestine_epithelium_phylum <- tax_glom(ra_lg_intestine_epithelium, taxrank = "Phylum", NArm = F) #41 phyla
lg_intestine_epithelium_class <- tax_glom(ra_lg_intestine_epithelium, taxrank = "Class", NArm = F) #40 classes
lg_intestine_epithelium_order <- tax_glom(ra_lg_intestine_epithelium, taxrank = "Order", NArm = F) # 74 orders
lg_intestine_epithelium_family <- tax_glom(ra_lg_intestine_epithelium, taxrank = "Family", NArm = F) #128 families
write.csv(otu_table(lg_intestine_epithelium_family),"lg_intestine_epithelium_otus.csv")
write.csv(tax_table(lg_intestine_epithelium_family),"lg_intestine_epithelium_taxa.csv")
lg_intestine_epithelium_family_melt <- psmelt(lg_intestine_epithelium_family)
lg_intestine_epithelium_genus <- tax_glom(ra_lg_intestine_epithelium, taxrank = "Genus", NArm = F) # 286 genera
lg_intestine_epithelium_genus_melt <- psmelt(lg_intestine_epithelium_genus)

#### ALL GUT SAMPLES
gut_samples_ra <- transform_sample_counts(gut_samples.css, function(x) {x/sum(x)} * 100)

######## UNIFRAC AND ORDINATIONS SPLIT BY BODY SITE #######

#### LIV ABS
liv_abs.gunifrac <- gunifrac(liv_abs.css)
liv_abs.ord <- ordinate(liv_abs.css, method = "NMDS", distance = liv_abs.gunifrac)

#### LA AVERAGE PER ANIMAL
# for some reason my function/wrapper for GUniFrac spazzes out with the merged table so I'm using the original version from 'GUniFrac' here
liver_abscesses_merged.dist <- GUniFrac(otu_table(liver_abscesses_merged.css), tree = phy_tree(liver_abscesses_merged.css), alpha = 0.5)$unifracs
liver_abscesses_merged.gunifrac.matrix <- liver_abscesses_merged.dist[,,"d_0.5"]
liver_abscesses_merged.gunifrac <-as.dist(liver_abscesses_merged.gunifrac.matrix)
liver_abscesses_merged.ord <- ordinate(liver_abscesses_merged.css, method = "NMDS", distance = liver_abscesses_merged.gunifrac)

#### MULTIPLE LAs
multiple_LAs.gunifrac <- gunifrac(multiple_LAs.css)
multiple_LAs.ord <- ordinate(multiple_LAs.css, method = "NMDS", distance = multiple_LAs.gunifrac)

#### RUMEN
rumen.gunifrac <- gunifrac(rumen.css)
rumen.gunifrac.ord <- ordinate(rumen.css, method = "NMDS", dist = rumen.gunifrac)

#### SM. INTESTINE
sm_intestine.gunifrac <- gunifrac(sm_intestine.css)
sm_intestine.gunifrac.ord <- ordinate(sm_intestine.css, method = "NMDS", dist = sm_intestine.gunifrac)

#### LG. INTESTINE
lg_intestine.gunifrac <- gunifrac(lg_intestine.css)
lg_intestine.gunifrac.ord <- ordinate(lg_intestine.css, method = "NMDS", dist = lg_intestine.gunifrac)


#### CHECKING IF WE SHOULD CONSIDER LUMEN & EPI SEPARATELY #######
#### RUMEN
plot_ordination(rumen.css, rumen.gunifrac.ord, type = "samples", color = "specific_body_loc") +
  theme_bw() +
  labs(title = "RUMEN", x= "Axis 1 (30.2% variation explained)", y= "Axis 2 (21.1% variation explained)") +
  geom_point(size = 5, shape = 18) +
  stat_ellipse(geom = "polygon", aes(fill = specific_body_loc), alpha = 0.5, lty = 2, size = 1) +
  scale_colour_manual(values = lumen_epi_palette) +
  scale_fill_manual(values = lumen_epi_palette) +
  theme(legend.position = "none",
        panel.border = element_rect(colour= "black", size = 1),
        title = element_text(size = 28),
        axis.ticks = element_line(colour = "black", size = 0.75),
        axis.text = element_text(colour = "black", size = 12),
        axis.title = element_text(size = 24))

adonis2(rumen.gunifrac ~ specific_body_loc, rumen.css.df) # sig.

#### SM. INTESTINE
plot_ordination(sm_intestine.css, sm_intestine.gunifrac.ord, type = "samples", color = "specific_body_loc") +
  theme_bw() +
  labs(title = "SM. INTESTINE", x= "Axis 1 (39.5% variation explained", y= "Axis 2 (12.6% variation explained)") +
  geom_point(size = 5, shape = 18) +
  stat_ellipse(geom = "polygon", aes(fill = specific_body_loc), alpha = 0.5, lty = 2, size = 1) +
  scale_colour_manual(values = lumen_epi_palette) +
  scale_fill_manual(values = lumen_epi_palette) +
  theme(legend.position = "none",
        panel.border = element_rect(colour= "black", size = 1),
        title = element_text(size = 28),
        axis.ticks = element_line(colour = "black", size = 0.75),
        axis.text = element_text(colour = "black", size = 12),
        axis.title = element_text(size = 24))

adonis2(sm_intestine.gunifrac ~ specific_body_loc, sm_intestine.css.df) # NS

#### LG. INTESTINE
plot_ordination(lg_intestine.css, lg_intestine.gunifrac.ord, type = "samples", color = "specific_body_loc") +
  theme_bw() +
  labs(title = "LG. INTESTINE", x= "Axis 1 (22.1% variation explained", y= "Axis 2 (15.5% variation explained") +
  geom_point(size = 5, shape = 18) +
  stat_ellipse(geom = "polygon", aes(fill = specific_body_loc), alpha = 0.5, lty = 2, size = 1) +
  scale_colour_manual(values = lumen_epi_palette) +
  scale_fill_manual(values = lumen_epi_palette) +
  theme(legend.position = "none",
        panel.border = element_rect(colour= "black", size = 1),
        title = element_text(size = 28),
        axis.ticks = element_line(colour = "black", size = 0.75),
        axis.text = element_text(colour = "black", size = 12),
        axis.title = element_text(size = 24))

adonis2(lg_intestine.gunifrac ~ specific_body_loc, lg_intestine.css.df) # NS

#### going to consider them different even though they are similar in intestine
#### since they are very different in rumen and it might be important for diff. abundance of taxa


######## SPLITTING GIT CSS LOCATIONS BY LUMEN & EPI ########

#### ALL GIT LUMEN
gut_lumen.css <- subset_samples(gut_samples.css, matrix=="lumen")
gut_lumen.css <- prune_taxa(taxa_sums(gut_lumen.css) > 0, gut_lumen.css)
gut_lumen.css # 101 samples, 7465 samples

#### ALL GIT EPITHELIUM
gut_epithelium.css <- subset_samples(gut_samples.css, matrix=="epithelium")
gut_epithelium.css <- prune_taxa(taxa_sums(gut_epithelium.css) > 0, gut_epithelium.css)
gut_epithelium.css # 102 samples, 9270 samples

#### RUMEN EPITHELIUM
rumen_epithelium.css <- subset_samples(samples.css, specific_body_loc=="rumen epithelium")
rumen_epithelium.css <- prune_taxa(taxa_sums(rumen_epithelium.css) > 0, rumen_epithelium.css)
rumen_epithelium.css # 4,329 taxa and 37 samples
rumen_epithelium.css.df <- as(sample_data(rumen_epithelium.css),"data.frame")

#### RUMEN LUMEN
rumen_lumen.css <- subset_samples(samples.css, specific_body_loc=="rumen lumen")
rumen_lumen.css <- prune_taxa(taxa_sums(rumen_lumen.css) > 0, rumen_lumen.css)
rumen_lumen.css # 3,114 taxa and 34 samples
rumen_lumen.css.df <- as(sample_data(rumen_lumen.css),"data.frame")

#### SM. INTESTINE EPITHELIUM
sm_intestine_epithelium.css <- subset_samples(samples.css, specific_body_loc=="small intestine epithelium")
sm_intestine_epithelium.css <- prune_taxa(taxa_sums(sm_intestine_epithelium.css) > 0, sm_intestine_epithelium.css)
sm_intestine_epithelium.css # 1612 taxa and 30 samples
sm_intestine_epithelium.css.df <- as(sample_data(sm_intestine_epithelium.css),"data.frame")

#### SM. INTESTINE LUMEN
sm_intestine_lumen.css <- subset_samples(samples.css, specific_body_loc=="small intestine lumen")
sm_intestine_lumen.css <- prune_taxa(taxa_sums(sm_intestine_lumen.css) > 0, sm_intestine_lumen.css)
sm_intestine_lumen.css # 1231 taxa and 31 samples
sm_intestine_lumen.css.df <- as(sample_data(sm_intestine_lumen.css),"data.frame")

#### LG. INTESTINE EPITHELIUM
lg_intestine_epithelium.css <- subset_samples(samples.css, specific_body_loc=="large intestine epithelium")
lg_intestine_epithelium.css <- prune_taxa(taxa_sums(lg_intestine_epithelium.css) > 0, lg_intestine_epithelium.css)
lg_intestine_epithelium.css # 5284 taxa and 35 samples
lg_intestine_epithelium.css.df <- as(sample_data(lg_intestine_epithelium.css),"data.frame")

#### LG. INTESTINE LUMEN
lg_intestine_lumen.css <- subset_samples(samples.css, specific_body_loc=="large intestine lumen")
lg_intestine_lumen.css <- prune_taxa(taxa_sums(lg_intestine_lumen.css) > 0, lg_intestine_lumen.css)
lg_intestine_lumen.css # 4768 taxa and 36 samples
lg_intestine_lumen.css.df <- as(sample_data(lg_intestine_lumen.css),"data.frame")

#### UNIFRAC AND ORDINATIONS ON BODY SITE LUMEN & EPITH ########

### ALL GIT
gut_lumen.dist <- gunifrac(gut_lumen.css)
gut_lumen.ord <- ordinate(gut_lumen.css, method = "NMDS", distance = gut_lumen.dist)

gut_epithelium.dist <- gunifrac(gut_epithelium.css)
gut_epithelium.ord <- ordinate(gut_epithelium.css, method = "NMDS", distance = gut_epithelium.dist)

#### RUMEN
rumen_epithelium.dist <- gunifrac(rumen_epithelium.css)
rumen_epithelium.uwunifrac <- uwunifrac(rumen_epithelium.css)
rumen_epithelium.wunifrac <- wunifrac(rumen_epithelium.css)
rumen_epithelium.bray <- distance(rumen_epithelium.css, method = "bray")
rumen_epithelium.ord <- ordinate(rumen_epithelium.css, method = "NMDS", distance = rumen_epithelium.dist)
rumen_epithelium_uw.ord <- ordinate(rumen_epithelium.css, method = "NMDS", distance = rumen_epithelium.uwunifrac)
rumen_epithelium_w.ord <- ordinate(rumen_epithelium.css, method = "NMDS", distance = rumen_epithelium.wunifrac)

rumen_lumen.dist <- gunifrac(rumen_lumen.css)
rumen_lumen.uwunifrac <- uwunifrac(rumen_lumen.css)
rumen_lumen.wunifrac <- wunifrac(rumen_lumen.css)
rumen_lumen.ord <- ordinate(rumen_lumen.css, method = "NMDS", distance = rumen_lumen.dist)
rumen_lumen_uw.ord <- ordinate(rumen_lumen.css, method = "NMDS", distance = rumen_lumen.uwunifrac)
rumen_lumen_w.ord <- ordinate(rumen_lumen.css, method = "NMDS", distance = rumen_lumen.wunifrac)

#### SM. INTESTINE
sm_intestine_epithelium.dist <- gunifrac(sm_intestine_epithelium.css)
sm_intestine_epithelium.ord <- ordinate(sm_intestine_epithelium.css, method = "NMDS", distance = sm_intestine_epithelium.dist)

sm_intestine_lumen.dist <- gunifrac(sm_intestine_lumen.css)
sm_intestine_lumen.ord <- ordinate(sm_intestine_lumen.css, method = "NMDS", distance = sm_intestine_lumen.dist)

#### LG. INTESTINE
lg_intestine_epithelium.dist <- gunifrac(lg_intestine_epithelium.css)
lg_intestine_epithelium.ord <- ordinate(lg_intestine_epithelium.css, method = "NMDS", distance = lg_intestine_epithelium.dist)

lg_intestine_lumen.dist <- gunifrac(lg_intestine_lumen.css)
lg_intestine_lumen.ord <- ordinate(lg_intestine_lumen.css, method = "NMDS", distance = lg_intestine_lumen.dist)



#### ORDINATIONS OF GIT LOCATIONS VS EACH OTHER ####
plot_ordination(gut_lumen.css, gut_lumen.ord, color = "specific_body_loc") +
  stat_ellipse() +
  geom_text(label = gut_lumen.css@sam_data$Animal)

plot_ordination(gut_epithelium.css, gut_epithelium.ord, color = "specific_body_loc") +
  stat_ellipse() +
  geom_text(label = gut_epithelium.css@sam_data$Animal)

################## Tylosin - NMDS PLOTS#######

#### RUMEN EPITH.
# find centroids
rumen_epithelium_tylosin_plot <- ordiplot(rumen_epithelium.ord$points)
rumen_epithelium_tylosin_siteslong <- sites.long(rumen_epithelium_tylosin_plot, rumen_epithelium.css.df)
rumen_epithelium_tylosin_centroids <- envfit(rumen_epithelium.ord ~ rumen_epithelium.css.df$treatment)
rumen_epithelium_tylosin_centroids
# make df for centroids
rumen_epithelium_tylosin_col <- c("control","tylosin")
rumen_epithelium_tylosin_NMDS1_col1 <- c(-0.0649,0.0551)
rumen_epithelium_tylosin_NMDS2_col2 <- c(-0.0123,0.0104)

rumen_epithelium_tylosin_centroids.df <- data.frame(rumen_epithelium_tylosin_col,rumen_epithelium_tylosin_NMDS1_col1, rumen_epithelium_tylosin_NMDS2_col2)
rumen_epithelium_tylosin_centroids.df

ggplot() + theme_bw() +
  labs(x= "NMDS1", y= "NMDS2") +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL), expand = c(0.08,0,0.027,0)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL), expand = c(0.05,0,0.15,0)) +
  scale_shape_manual(values = c(18,18)) +
  geom_point(data = rumen_epithelium_tylosin_siteslong, aes(x=axis1,y=axis2, colour= treatment, shape = treatment, alpha = treatment, size = treatment)) +
  scale_alpha_manual(values = c(0.7,0.7), guide = F) +
  scale_size_manual(values = c(5,5), guide = "none") +
  stat_ellipse(data = rumen_epithelium_tylosin_siteslong, geom = "polygon", aes(x=axis1,y=axis2, colour =treatment, fill = treatment), alpha = c(0.3), lty=2, level = 0.90, linewidth = 1) +
  geom_point(data = rumen_epithelium_tylosin_centroids.df, aes(x=rumen_epithelium_tylosin_NMDS1_col1, y=rumen_epithelium_tylosin_NMDS2_col2), fill = tylosin_palette, colour = tylosin_palette, size = 18, shape = c(18,18)) +
  geom_text(data = rumen_epithelium_tylosin_centroids.df, aes(x=rumen_epithelium_tylosin_NMDS1_col1, y=rumen_epithelium_tylosin_NMDS2_col2, label = c("C","T")), colour = "white", size = 9, fontface = "bold") +
  scale_colour_manual(values = tylosin_palette) +
  scale_fill_manual(values = tylosin_palette) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        axis.ticks = element_line(colour = "black", size = 0.75),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 12, colour = "black"))

adonis2(rumen_epithelium.dist ~ treatment, rumen_epithelium.css.df) # significant 
rumen_epithelium.disper <- betadisper(rumen_epithelium.dist, rumen_epithelium.css.df$treatment)
plot(rumen_epithelium.disper)
rumen_epithelium.permdisp <- permutest(rumen_epithelium.disper, permutations = 9999)
rumen_epithelium.permdisp #NS

#### RUMEN LUMEN
# find centroids
rumen_lumen_tylosin_plot <- ordiplot(rumen_lumen.ord$points)
rumen_lumen_tylosin_siteslong <- sites.long(rumen_lumen_tylosin_plot, rumen_lumen.css.df)
rumen_lumen_tylosin_centroids <- envfit(rumen_lumen.ord ~ rumen_lumen.css.df$treatment)
rumen_lumen_tylosin_centroids
# make df for centroids
rumen_lumen_tylosin_col <- c("control","tylosin")
rumen_lumen_tylosin_NMDS1_col1 <- c(0.0751,-0.0593)
rumen_lumen_tylosin_NMDS2_col2 <- c(-0.0245,0.0193)

rumen_lumen_tylosin_centroids.df <- data.frame(rumen_lumen_tylosin_col,rumen_lumen_tylosin_NMDS1_col1, rumen_lumen_tylosin_NMDS2_col2)
rumen_lumen_tylosin_centroids.df

ggplot() + theme_bw() +
  labs(x= "NMDS1", y= "NMDS2") +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL), limits = c(-0.17,0.2)) +
  scale_shape_manual(values = c(18,18)) +
  geom_point(data = rumen_lumen_tylosin_siteslong, aes(x=axis1,y=axis2, colour= treatment, shape = treatment, alpha = treatment, size = treatment)) +
  scale_alpha_manual(values = c(0.7,0.7), guide = F) +
  scale_size_manual(values = c(5,5), guide = "none") +
  stat_ellipse(data = rumen_lumen_tylosin_siteslong, geom = "polygon", aes(x=axis1,y=axis2, colour =treatment, fill = treatment), alpha = c(0.3), lty=2, level = 0.90, linewidth = 1) +
  geom_point(data = rumen_lumen_tylosin_centroids.df, aes(x=rumen_lumen_tylosin_NMDS1_col1, y=rumen_lumen_tylosin_NMDS2_col2), fill = tylosin_palette, colour = tylosin_palette, size = 18, shape = c(18,18)) +
  geom_text(data = rumen_lumen_tylosin_centroids.df, aes(x=rumen_lumen_tylosin_NMDS1_col1, y=rumen_lumen_tylosin_NMDS2_col2, label = c("C","T")), colour = "white", size = 9, fontface = "bold") +
  scale_colour_manual(values = tylosin_palette) +
  scale_fill_manual(values = tylosin_palette) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        axis.ticks = element_line(colour = "black", size = 0.75),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 12, colour = "black"))


plot_ordination(rumen_lumen.css, rumen_lumen.ord, type = "samples", color = "treatment") +
  theme_bw() +
  labs(x= "NMDS1", y= "NMDS2") +
  geom_point(size = 5, shape = 18) + 
  stat_ellipse(geom = "polygon", aes(fill= treatment), alpha = 0.5, lty = 2, size = 1, level =0.90) +
  scale_colour_manual(values = tylosin_palette) +
  scale_fill_manual(values = tylosin_palette) +
  scale_y_continuous(limits = c(-0.17,0.2)) +
  theme(legend.position = "none",
        panel.border = element_rect(colour= "black", size = 1),
        title = element_text(size = 28),
        axis.ticks = element_line(colour = "black", size = 0.75),
        axis.text = element_text(colour = "black", size = 12),
        axis.title = element_text(size = 24))

adonis2(rumen_lumen.dist ~ treatment, rumen_lumen.css.df) # significant 
rumen_lumen.disper <- betadisper(rumen_lumen.dist, rumen_lumen.css.df$treatment)
plot(rumen_lumen.disper)
rumen_lumen.permdisp <- permutest(rumen_lumen.disper, permutations = 9999)
rumen_lumen.permdisp #NS

#### SM. INTESTINE EPITH.
# find centroids
sm_intestine_epithelium_tylosin_plot <- ordiplot(sm_intestine_epithelium.ord$points)
sm_intestine_epithelium_tylosin_siteslong <- sites.long(sm_intestine_epithelium_tylosin_plot, sm_intestine_epithelium.css.df)
sm_intestine_epithelium_tylosin_centroids <- envfit(sm_intestine_epithelium.ord ~ sm_intestine_epithelium.css.df$treatment)
sm_intestine_epithelium_tylosin_centroids
# make df for centroids
sm_intestine_epithelium_tylosin_col <- c("control","tylosin")
sm_intestine_epithelium_tylosin_NMDS1_col1 <- c(-0.0041,0.0036)
sm_intestine_epithelium_tylosin_NMDS2_col2 <- c(0.0795,-0.0695)

sm_intestine_epithelium_tylosin_centroids.df <- data.frame(sm_intestine_epithelium_tylosin_col,sm_intestine_epithelium_tylosin_NMDS1_col1, sm_intestine_epithelium_tylosin_NMDS2_col2)
sm_intestine_epithelium_tylosin_centroids.df

ggplot() + theme_bw() +
  labs(x= "NMDS1", y= "NMDS2") +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL), expand = c(0.005,0.05,0.07,0.05)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL), expand = c(0.005,0.05,0.18,0.05)) +
  scale_shape_manual(values = c(18,18)) +
  geom_point(data = sm_intestine_epithelium_tylosin_siteslong, aes(x=axis1,y=axis2, colour= treatment, shape = treatment, alpha = treatment, size = treatment)) +
  scale_alpha_manual(values = c(0.7,0.7), guide = F) +
  scale_size_manual(values = c(5,5), guide = "none") +
  stat_ellipse(data = sm_intestine_epithelium_tylosin_siteslong, geom = "polygon", aes(x=axis1,y=axis2, colour =treatment, fill = treatment), alpha = c(0.3), lty=2, level = 0.90, linewidth = 1) +
  geom_point(data = sm_intestine_epithelium_tylosin_centroids.df, aes(x=sm_intestine_epithelium_tylosin_NMDS1_col1, y=sm_intestine_epithelium_tylosin_NMDS2_col2), fill = tylosin_palette, colour = tylosin_palette, size = 18, shape = c(18,18)) +
  geom_text(data = sm_intestine_epithelium_tylosin_centroids.df, aes(x=sm_intestine_epithelium_tylosin_NMDS1_col1, y=sm_intestine_epithelium_tylosin_NMDS2_col2, label = c("C","T")), colour = "white", size = 9, fontface = "bold") +
  scale_colour_manual(values = tylosin_palette) +
  scale_fill_manual(values = tylosin_palette) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        axis.ticks = element_line(colour = "black", size = 0.75),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 12, colour = "black"))

adonis2(sm_intestine_epithelium.dist ~ treatment, sm_intestine_epithelium.css.df) # significant 
sm_intestine_epithelium.disper <- betadisper(sm_intestine_epithelium.dist, sm_intestine_epithelium.css.df$treatment)
plot(sm_intestine_epithelium.disper)
sm_intestine_epithelium.permdisp <- permutest(sm_intestine_epithelium.disper, permutations = 9999)
sm_intestine_epithelium.permdisp #NS

#### SM. INTESTINE LUMEN
# find centroids
sm_intestine_lumen_tylosin_plot <- ordiplot(sm_intestine_lumen.ord$points)
sm_intestine_lumen_tylosin_siteslong <- sites.long(sm_intestine_lumen_tylosin_plot, sm_intestine_lumen.css.df)
sm_intestine_lumen_tylosin_centroids <- envfit(sm_intestine_lumen.ord ~ sm_intestine_lumen.css.df$treatment)
sm_intestine_lumen_tylosin_centroids
# make df for centroids
sm_intestine_lumen_tylosin_col <- c("control","tylosin")
sm_intestine_lumen_tylosin_NMDS1_col1 <- c(-0.0212,0.0119)
sm_intestine_lumen_tylosin_NMDS2_col2 <- c(0.0199,-0.0111)

sm_intestine_lumen_tylosin_centroids.df <- data.frame(sm_intestine_lumen_tylosin_col,sm_intestine_lumen_tylosin_NMDS1_col1, sm_intestine_lumen_tylosin_NMDS2_col2)
sm_intestine_lumen_tylosin_centroids.df

ggplot() + theme_bw() +
  labs(x= "NMDS1", y= "NMDS2") +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL), expand = c(0.005,0.01,0.01,0.01)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL), expand = c(0.005,0.01,0.02,0.01)) +
  scale_shape_manual(values = c(18,18)) +
  geom_point(data = sm_intestine_lumen_tylosin_siteslong, aes(x=axis1,y=axis2, colour= treatment, shape = treatment, alpha = treatment, size = treatment)) +
  scale_alpha_manual(values = c(0.7,0.7), guide = F) +
  scale_size_manual(values = c(5,5), guide = "none") +
  stat_ellipse(data = sm_intestine_lumen_tylosin_siteslong, geom = "polygon", aes(x=axis1,y=axis2, colour =treatment, fill = treatment), alpha = c(0.3), lty=2, level = 0.90, linewidth = 1) +
  geom_point(data = sm_intestine_lumen_tylosin_centroids.df, aes(x=sm_intestine_lumen_tylosin_NMDS1_col1, y=sm_intestine_lumen_tylosin_NMDS2_col2), fill = tylosin_palette, colour = tylosin_palette, size = 18, shape = c(18,18)) +
  geom_text(data = sm_intestine_lumen_tylosin_centroids.df, aes(x=sm_intestine_lumen_tylosin_NMDS1_col1, y=sm_intestine_lumen_tylosin_NMDS2_col2, label = c("C","T")), colour = "white", size = 9, fontface = "bold") +
  scale_colour_manual(values = tylosin_palette) +
  scale_fill_manual(values = tylosin_palette) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        axis.ticks = element_line(colour = "black", size = 0.75),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 12, colour = "black"))

# stats
adonis2(sm_intestine_lumen.dist ~ treatment, sm_intestine_lumen.css.df) # sig.
sm_intestine_lumen.disper <- betadisper(sm_intestine_lumen.dist, sm_intestine_lumen.css.df$treatment)
plot(sm_intestine_lumen.disper)
sm_intestine_lumen.permdisp <- permutest(sm_intestine_lumen.disper, permutations = 9999)
sm_intestine_lumen.permdisp #NS

#### LG. INTESTINE EPITH.
# find centroids
lg_intestine_epithelium_tylosin_plot <- ordiplot(lg_intestine_epithelium.ord$points)
lg_intestine_epithelium_tylosin_siteslong <- sites.long(lg_intestine_epithelium_tylosin_plot, lg_intestine_epithelium.css.df)
lg_intestine_epithelium_tylosin_centroids <- envfit(lg_intestine_epithelium.ord ~ lg_intestine_epithelium.css.df$treatment)
lg_intestine_epithelium_tylosin_centroids
# make df for centroids
lg_intestine_epithelium_tylosin_col <- c("control","tylosin")
lg_intestine_epithelium_tylosin_NMDS1_col1 <- c(-0.0125,0.0118)
lg_intestine_epithelium_tylosin_NMDS2_col2 <- c(0.0059,-0.0056)

lg_intestine_epithelium_tylosin_centroids.df <- data.frame(lg_intestine_epithelium_tylosin_col,lg_intestine_epithelium_tylosin_NMDS1_col1, lg_intestine_epithelium_tylosin_NMDS2_col2)
lg_intestine_epithelium_tylosin_centroids.df

ggplot() + theme_bw() +
  labs(x= "NMDS1", y= "NMDS2") +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL), expand = c(0.02,0.02,0.02,0.02)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL), expand = c(0.02,0.02,0.02,0.02)) +
  scale_shape_manual(values = c(18,18)) +
  geom_point(data = lg_intestine_epithelium_tylosin_siteslong, aes(x=axis1,y=axis2, colour= treatment, shape = treatment, alpha = treatment, size = treatment)) +
  scale_alpha_manual(values = c(0.7,0.7), guide = F) +
  scale_size_manual(values = c(5,5), guide = "none") +
  stat_ellipse(data = lg_intestine_epithelium_tylosin_siteslong, geom = "polygon", aes(x=axis1,y=axis2, colour =treatment, fill = treatment), alpha = c(0.3), lty=2, level = 0.90, linewidth = 1) +
  geom_point(data = lg_intestine_epithelium_tylosin_centroids.df, aes(x=lg_intestine_epithelium_tylosin_NMDS1_col1, y=lg_intestine_epithelium_tylosin_NMDS2_col2), fill = tylosin_palette, colour = tylosin_palette, size = 18, shape = c(18,18)) +
  geom_text(data = lg_intestine_epithelium_tylosin_centroids.df, aes(x=lg_intestine_epithelium_tylosin_NMDS1_col1, y=lg_intestine_epithelium_tylosin_NMDS2_col2, label = c("C","T")), colour = "white", size = 9, fontface = "bold") +
  scale_colour_manual(values = tylosin_palette) +
  scale_fill_manual(values = tylosin_palette) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        axis.ticks = element_line(colour = "black", size = 0.75),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 12, colour = "black"))

adonis2(lg_intestine_epithelium.dist ~ treatment, lg_intestine_epithelium.css.df, permutations = 9999) # NS
lg_intestine_epithelium.disper <- betadisper(lg_intestine_epithelium.dist, lg_intestine_epithelium.css.df$treatment)
plot(lg_intestine_epithelium.disper)
lg_intestine_epithelium.permdisp <- permutest(lg_intestine_epithelium.disper, permutations = 9999)
lg_intestine_epithelium.permdisp #sig

#### SM. INTESTINE LUMEN
# find centroids
lg_intestine_lumen_tylosin_plot <- ordiplot(lg_intestine_lumen.ord$points)
lg_intestine_lumen_tylosin_siteslong <- sites.long(lg_intestine_lumen_tylosin_plot, lg_intestine_lumen.css.df)
lg_intestine_lumen_tylosin_centroids <- envfit(lg_intestine_lumen.ord ~ lg_intestine_lumen.css.df$treatment)
lg_intestine_lumen_tylosin_centroids
# make df for centroids
lg_intestine_lumen_tylosin_col <- c("control","tylosin")
lg_intestine_lumen_tylosin_NMDS1_col1 <- c(0.0167,0.0098)
lg_intestine_lumen_tylosin_NMDS2_col2 <- c(-0.0150,-0.0088)

lg_intestine_lumen_tylosin_centroids.df <- data.frame(lg_intestine_lumen_tylosin_col,lg_intestine_lumen_tylosin_NMDS1_col1, lg_intestine_lumen_tylosin_NMDS2_col2)
lg_intestine_lumen_tylosin_centroids.df

ggplot() + theme_bw() +
  labs(x= "NMDS1", y= "NMDS2") +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL), expand = c(0.01,0.01,0.02,0.01)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL), expand = c(0.01,0.02,0.02,0.02)) +
  scale_shape_manual(values = c(18,18)) +
  geom_point(data = lg_intestine_lumen_tylosin_siteslong, aes(x=axis1,y=axis2, colour= treatment, shape = treatment, alpha = treatment, size = treatment)) +
  scale_alpha_manual(values = c(0.7,0.7), guide = F) +
  scale_size_manual(values = c(5,5), guide = "none") +
  stat_ellipse(data = lg_intestine_lumen_tylosin_siteslong, geom = "polygon", aes(x=axis1,y=axis2, colour =treatment, fill = treatment), alpha = c(0.3), lty=2, level = 0.90, linewidth = 1) +
  geom_point(data = lg_intestine_lumen_tylosin_centroids.df, aes(x=lg_intestine_lumen_tylosin_NMDS1_col1, y=lg_intestine_lumen_tylosin_NMDS2_col2), fill = tylosin_palette, colour = tylosin_palette, size = 18, shape = c(18,18)) +
  geom_text(data = lg_intestine_lumen_tylosin_centroids.df, aes(x=lg_intestine_lumen_tylosin_NMDS1_col1, y=lg_intestine_lumen_tylosin_NMDS2_col2, label = c("C","T")), colour = "white", size = 9, fontface = "bold") +
  scale_colour_manual(values = tylosin_palette) +
  scale_fill_manual(values = tylosin_palette) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        axis.ticks = element_line(colour = "black", size = 0.75),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 12, colour = "black"))

# stats
adonis2(lg_intestine_lumen.dist ~ treatment, lg_intestine_lumen.css.df, permutations = 9999) # sig.
lg_intestine_lumen.disper <- betadisper(lg_intestine_lumen.dist, lg_intestine_lumen.css.df$treatment)
plot(lg_intestine_lumen.disper)
lg_intestine_lumen.permdisp <- permutest(lg_intestine_lumen.disper, permutations = 9999)
lg_intestine_lumen.permdisp #sig, believe this is cause for sig PERMANOVA



#### LIVER ABSCESSES
## Effect of tylosin ##
# find centroids
liver_abscesses_merged_tylosin_plot <- ordiplot(liver_abscesses_merged.ord$points)
liver_abscesses_merged_tylosin_siteslong <- sites.long(liver_abscesses_merged_tylosin_plot, liver_abscesses_merged.css.df)
liver_abscesses_merged_tylosin_centroids <- envfit(liver_abscesses_merged.ord ~ liver_abscesses_merged.css.df$treatment)
liver_abscesses_merged_tylosin_centroids
# make df for centroids
liver_abscesses_merged_tylosin_col <- c("control","tylosin")
liver_abscesses_merged_tylosin_NMDS1_col1 <- c(-0.0091,0.0091)
liver_abscesses_merged_tylosin_NMDS2_col2 <- c(-0.0065,0.0065)

liver_abscesses_merged_tylosin_centroids.df <- data.frame(liver_abscesses_merged_tylosin_col,liver_abscesses_merged_tylosin_NMDS1_col1, liver_abscesses_merged_tylosin_NMDS2_col2)
liver_abscesses_merged_tylosin_centroids.df

ggplot() + theme_bw() +
  labs(x= "NMDS1", y= "NMDS2") +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL), expand = c(0.02,0.02,0.02,0.02)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL), expand = c(0.02,0.02,0.02,0.02)) +
  scale_shape_manual(values = c(18,18)) +
  geom_point(data = liver_abscesses_merged_tylosin_siteslong, aes(x=axis1,y=axis2, colour= treatment, shape = treatment, alpha = treatment, size = treatment)) +
  scale_alpha_manual(values = c(0.7,0.7), guide = F) +
  scale_size_manual(values = c(5,5), guide = "none") +
  stat_ellipse(data = liver_abscesses_merged_tylosin_siteslong, geom = "polygon", aes(x=axis1,y=axis2, colour =treatment, fill = treatment), alpha = c(0.3), lty=2, level = 0.90, linewidth = 1) +
  geom_point(data = liver_abscesses_merged_tylosin_centroids.df, aes(x=liver_abscesses_merged_tylosin_NMDS1_col1, y=liver_abscesses_merged_tylosin_NMDS2_col2), fill = tylosin_palette, colour = tylosin_palette, size = 18, shape = c(18,18)) +
  geom_text(data = liver_abscesses_merged_tylosin_centroids.df, aes(x=liver_abscesses_merged_tylosin_NMDS1_col1, y=liver_abscesses_merged_tylosin_NMDS2_col2, label = c("C","T")), colour = "white", size = 9, fontface = "bold") +
  scale_colour_manual(values = tylosin_palette) +
  scale_fill_manual(values = tylosin_palette) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        axis.ticks = element_line(colour = "black", size = 0.75),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 12, colour = "black"))

adonis2(liver_abscesses_merged.gunifrac ~ treatment, liver_abscesses_merged.css.df, permutations = 9999) #NS

################## Tylosin - DENDROGRAM + RA PLOTS ######
### RUMEN####
## LUMEN
rumen_lumen.hclust <- hclust(rumen_lumen.dist, method = "ward.D2")
rumen_lumen.dendro <- as.dendrogram(rumen_lumen.hclust)
rumen_lumen.dendro.data <- dendro_data(rumen_lumen.dendro, type = "rectangle")
rumen_lumen.metadata_for_dendro <- as_tibble(rumen_lumen.css@sam_data)
rumen_lumen.dendro.data$labels <- rumen_lumen.dendro.data$labels %>%
  left_join(rumen_lumen.metadata_for_dendro, by = c("label" = "Sample"))

ggplot(rumen_lumen.dendro.data$segments) +
  theme_minimal() +
  labs(y= "Ward's Distance") +
  geom_segment(aes(x=x,y=y,xend=xend,yend=yend)) +
  geom_point(data = rumen_lumen.dendro.data$labels, 
             aes(x=x,y=y, colour = treatment, fill= treatment),
             size = 11, shape=22, stroke =1.5, position = position_nudge(y=-0.18)) +
  scale_y_continuous(limits = c(-.25,1.27)) +
  scale_x_discrete(expand = c(0.03,0,0.03,0)) +
  scale_fill_manual(values = alpha(tylosin_palette,0.5)) +
  scale_colour_manual(values = tylosin_palette) +
  theme(legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.line.y = element_line(size = 0.7, colour = "black"),
        axis.ticks.y = element_line(size = 0.75, colour = "black"),
        axis.title.y = element_text(size = 24),
        axis.text.y = element_text(size = 14, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())

# family RA plot for under dendro
rumen_lumen_family_filt_melt <- rumen_lumen_family_melt %>%
  group_by(Family) %>%
  mutate(mean_family_ra = mean(Abundance))

rumen_lumen_family_filt_melt$Family[rumen_lumen_family_filt_melt$mean_family_ra < 0.01] <- "zzzzOther"
length(unique(rumen_lumen_family_filt_melt$Family)) #64 families > 0.01%
write.csv(rumen_lumen_family_filt_melt$Family,"rumen_lumen_families.csv")
write.csv(rumen_lumen_family_filt_melt$mean_family_ra,"rumen_lumen_meanRA.csv")
rumen_lumen_family_palette <-c("#A4AE33",	"#9FDC5B",	"#C664A8",	"#F0E9DE",	"#8AD8F3",	"#48CFCA",	"#E7C57D",	"#9CF5E2",	"#B39AD0",	"palegreen2",	"#C9B2BF",	"#466D92",	"#C1EF9B",	"#90B4F2",	"#A3748D",	"#A1D676",	"#E9906D",	"#B6F56A",	"#593050",	"#DCE080",	"#7C7CC3",	"#4DA176",	"#F0D3B0",	"#F0EDA4",	"#5BA845",	"#7E589A",	"darkgoldenrod3",	"#6F97E9",	"yellow4",	"#ACD0EF",	"#5CF4F9",	"#CDC5B5",	"#D0E850",	"springgreen3",	"#CD90DF",	"bisque3",	"#CFF685",	"#9ED6DB",	"#9E82EA",	"darkorange1",	"#D5B0D6",	"#D08992",	"#F3B2E3",	"#8FF4B6",	"dodgerblue3",	"#8DE9EA",	"mediumorchid",	"firebrick2",	"#DEB2F3",	"#E6F132",	"chocolate3",	"#BCC0EA",	"#EA3C5D",	"hotpink",	"#B8EFC1",	"#EE3DD4",	"honeydew3",	"#C76ED7",	"#EACC62",	"#5FB0B1",	"chartreuse2",	"#71C9B2",	"#EF9A95",	"grey87")
rumen_lumen_dendro_sample_order <- rumen_lumen.dendro.data$labels$label

ggplot(rumen_lumen_family_filt_melt, aes(x= sample_Sample, y= Abundance, fill= Family)) +
  theme_minimal() +
  labs(y= "Relative Abundance (%)") +
  geom_bar(stat = "identity", colour = "black") +
  scale_x_discrete(limits = rumen_lumen_dendro_sample_order, expand = c(0.03,0,0.015,0)) +
  scale_fill_manual(values = rumen_lumen_family_palette) +
  theme(legend.position = "none",
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.line.y = element_line(size = 0.7, colour = "black"),
        axis.ticks.y = element_line(colour = "black", size = 0.75),
        axis.title.y = element_text(size = 24),
        axis.text.y = element_text(size = 14, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())

## stats
#Lachnospiraceae
rumen_lumen_Lachnospiraceae <- rumen_lumen_family_melt[which(rumen_lumen_family_melt$Family=="Lachnospiraceae"),]
kruskal.test(rumen_lumen_Lachnospiraceae$Abundance, rumen_lumen_Lachnospiraceae$treatment) # sig
#Prevotellaceae
rumen_lumen_Prevotellaceae <- rumen_lumen_family_melt[which(rumen_lumen_family_melt$Family=="Prevotellaceae"),]
kruskal.test(rumen_lumen_Prevotellaceae$Abundance, rumen_lumen_Prevotellaceae$treatment) # NS
#Atopobiaceae
rumen_lumen_Atopobiaceae <- rumen_lumen_family_melt[which(rumen_lumen_family_melt$Family=="Atopobiaceae"),]
kruskal.test(rumen_lumen_Atopobiaceae$Abundance, rumen_lumen_Atopobiaceae$treatment) #NS
#Oscillospiraceae
rumen_lumen_Oscillospiraceae <- rumen_lumen_family_melt[which(rumen_lumen_family_melt$Family=="Oscillospiraceae"),]
kruskal.test(rumen_lumen_Oscillospiraceae$Abundance, rumen_lumen_Oscillospiraceae$treatment) # NS 0.074
#Methanobacteriaceae
rumen_lumen_Methanobacteriaceae <- rumen_lumen_family_melt[which(rumen_lumen_family_melt$Family=="Methanobacteriaceae"),]
kruskal.test(rumen_lumen_Methanobacteriaceae$Abundance, rumen_lumen_Methanobacteriaceae$treatment) # NS 0.218
#Anaerovoracaceae
rumen_lumen_Anaerovoracaceae <- rumen_lumen_family_melt[which(rumen_lumen_family_melt$Family=="Anaerovoracaceae"),]
kruskal.test(rumen_lumen_Anaerovoracaceae$Abundance, rumen_lumen_Anaerovoracaceae$treatment) # NS (0.069)
#Ruminococcaceae
rumen_lumen_Ruminococcaceae <- rumen_lumen_family_melt[which(rumen_lumen_family_melt$Family=="Ruminococcaceae"),]
kruskal.test(rumen_lumen_Ruminococcaceae$Abundance, rumen_lumen_Ruminococcaceae$treatment) # NS (0.245)
#Muribaculaceae
rumen_lumen_Muribaculaceae <- rumen_lumen_family_melt[which(rumen_lumen_family_melt$Family=="Muribaculaceae"),]
kruskal.test(rumen_lumen_Muribaculaceae$Abundance, rumen_lumen_Muribaculaceae$treatment) # NS (0.107)
#Clostridia UCG-014
rumen_lumen_Clostridia_UCG014 <- rumen_lumen_family_melt[which(rumen_lumen_family_melt$Family=="Clostridia_UCG-014"),]
kruskal.test(rumen_lumen_Clostridia_UCG014$Abundance, rumen_lumen_Clostridia_UCG014$treatment) # sig (0.003)
#Erysipelotrichaceae
rumen_lumen_Erysipelotrichaceae <- rumen_lumen_family_melt[which(rumen_lumen_family_melt$Family=="Erysipelotrichaceae"),]
kruskal.test(rumen_lumen_Erysipelotrichaceae$Abundance, rumen_lumen_Erysipelotrichaceae$treatment) # NS (0.245)
#Acidaminococcaceae
rumen_lumen_Acidaminococcaceae <- rumen_lumen_family_melt[which(rumen_lumen_family_melt$Family=="Acidaminococcaceae"),]
kruskal.test(rumen_lumen_Acidaminococcaceae$Abundance, rumen_lumen_Acidaminococcaceae$treatment) # NS (0.499)
#Rikenellaceae
rumen_lumen_Rikenellaceae<- rumen_lumen_family_melt[which(rumen_lumen_family_melt$Family=="Rikenellaceae"),]
kruskal.test(rumen_lumen_Rikenellaceae$Abundance, rumen_lumen_Rikenellaceae$treatment) # NS (0.252)
#Selenomonadaceae
rumen_lumen_Selenomonadaceae <- rumen_lumen_family_melt[which(rumen_lumen_family_melt$Family=="Selenomonadaceae"),]
kruskal.test(rumen_lumen_Selenomonadaceae$Abundance, rumen_lumen_Selenomonadaceae$treatment) # NS (0.396)
#Bacteroidales_RF16_group
rumen_lumen_Bacteroidales_RF16_group <- rumen_lumen_family_melt[which(rumen_lumen_family_melt$Family=="Bacteroidales_RF16_group"),]
kruskal.test(rumen_lumen_Bacteroidales_RF16_group$Abundance, rumen_lumen_Bacteroidales_RF16_group$treatment) # NS (0.127)
#Eubacterium_coprostanoligenes_group
rumen_lumen_Eubacterium_coprostanoligenes_group <- rumen_lumen_family_melt[which(rumen_lumen_family_melt$Family=="[Eubacterium]_coprostanoligenes_group"),]
kruskal.test(rumen_lumen_Eubacterium_coprostanoligenes_group$Abundance, rumen_lumen_Eubacterium_coprostanoligenes_group$treatment) # NS (0.377)

#### RUMEN EPITHELIUM
rumen_epithelium.hclust <- hclust(rumen_epithelium.dist, method = "ward.D2")
rumen_epithelium.dendro <- as.dendrogram(rumen_epithelium.hclust)
rumen_epithelium.dendro.data <- dendro_data(rumen_epithelium.dendro, type = "rectangle")
rumen_epithelium.metadata_for_dendro <- as_tibble(rumen_epithelium.css@sam_data)
rumen_epithelium.dendro.data$labels <- rumen_epithelium.dendro.data$labels %>%
  left_join(rumen_epithelium.metadata_for_dendro, by = c("label" = "Sample"))

ggplot(rumen_epithelium.dendro.data$segments) +
  theme_minimal() +
  labs(y= "Ward's Distance") +
  geom_segment(aes(x=x,y=y,xend=xend,yend=yend)) +
  geom_point(data = rumen_epithelium.dendro.data$labels, 
             aes(x=x,y=y, colour = treatment, fill= treatment),
             size = 10, shape=22, stroke =1.5, position = position_nudge(y=-0.18)) +
  scale_y_continuous(limits = c(-.25,1.27)) +
  scale_x_discrete(expand = c(0.03,0,0.03,0)) +
  scale_fill_manual(values = alpha(tylosin_palette,0.5)) +
  scale_colour_manual(values = tylosin_palette) +
  theme(legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.line.y = element_line(size = 0.7, colour = "black"),
        axis.ticks.y = element_line(size = 0.75, colour = "black"),
        axis.title.y = element_text(size = 24),
        axis.text.y = element_text(size = 14, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())

# family RA plot for under dendro
rumen_epithelium_family_filt_melt <- rumen_epithelium_family_melt %>%
  group_by(Family) %>%
  mutate(mean_family_ra = mean(Abundance))

rumen_epithelium_family_filt_melt$Family[rumen_epithelium_family_filt_melt$mean_family_ra < 0.01] <- "zzzzOther"
length(unique(rumen_epithelium_family_filt_melt$Family)) #67 families > 0.01%
write.csv(rumen_epithelium_family_filt_melt$Family,"rumen_epithelium_families.csv")
write.csv(rumen_epithelium_family_filt_melt$mean_family_ra,"rumen_epithelium_meanRA.csv")
rumen_epithelium_family_palette <- c("#A4AE33",	"#9FDC5B",	"#E33898",	"#C664A8",	"#42EDC1",	"#F0E9DE",	"#8AD8F3",	"#48CFCA",	"#E7C57D",	"springgreen4",	"#B39AD0",	"palegreen2",	"#C9B2BF",	"#466D92",	"#D94DEF",	"#50EDDA",	"#90B4F2",	"#A3748D",	"#A1D676",	"#E9906D",	"#55BDE5",	"#B6F56A",	"#593050",	"#DCE080",	"#7C7CC3",	"#4DA176",	"#F0EDA4",	"#5BA845",	"#7E589A",	"darkgoldenrod3",	"#6F97E9",	"yellow4",	"#ACD0EF",	"#5CF4F9",	"#CDC5B5",	"springgreen3",	"bisque3",	"#CFF685",	"#9ED6DB",	"#9E82EA",	"#9D26A4",	"darkorange1",	"#D5B0D6",	"#D08992",	"#F3B2E3",	"#8FF4B6",	"dodgerblue3",	"#8DE9EA",	"#DDF2BA",	"mediumorchid",	"firebrick2",	"#DEB2F3",	"#E6F132",	"chocolate3",	"#BCC0EA",	"#EA3C5D",	"hotpink",	"#B8EFC1",	"#90EB8C",	"#EE3DD4",	"honeydew3",	"#EACC62",	"#5FB0B1",	"chartreuse2",	"#71C9B2",	"#EF9A95",	"grey87")
rumen_epithelium_dendro_sample_order <- rumen_epithelium.dendro.data$labels$label

ggplot(rumen_epithelium_family_filt_melt, aes(x= sample_Sample, y= Abundance, fill= Family)) +
  theme_minimal() +
  labs(y= "Relative Abundance (%)") +
  geom_bar(stat = "identity", colour = "black") +
  scale_x_discrete(limits = rumen_epithelium_dendro_sample_order, expand = c(0.03,0,0.015,0)) +
  scale_fill_manual(values = rumen_epithelium_family_palette) +
  theme(legend.position = "none",
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.line.y = element_line(size = 0.7, colour = "black"),
        axis.ticks.y = element_line(colour = "black", size = 0.75),
        axis.title.y = element_text(size = 24),
        axis.text.y = element_text(size = 14, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())

## stats
#Prevotellaceae
rumen_epithelium_Prevotellaceae <- rumen_epithelium_family_melt[which(rumen_epithelium_family_melt$Family=="Prevotellaceae"),]
kruskal.test(rumen_epithelium_Prevotellaceae$Abundance, rumen_epithelium_Prevotellaceae$treatment) # NS (0.976)
#Lachnospiraceae
rumen_epithelium_Lachnospiraceae <- rumen_epithelium_family_melt[which(rumen_epithelium_family_melt$Family=="Lachnospiraceae"),]
kruskal.test(rumen_epithelium_Lachnospiraceae$Abundance, rumen_epithelium_Lachnospiraceae$treatment) # NS(0.784)
#Succinivibrionaceae
rumen_epithelium_Succinivibrionaceae <- rumen_epithelium_family_melt[which(rumen_epithelium_family_melt$Family=="Succinivibrionaceae"),]
kruskal.test(rumen_epithelium_Succinivibrionaceae$Abundance, rumen_epithelium_Succinivibrionaceae$treatment) #sig. (0.005)
#Selenomonadaceae
rumen_epithelium_Selenomonadaceae <- rumen_epithelium_family_melt[which(rumen_epithelium_family_melt$Family=="Selenomonadaceae"),]
kruskal.test(rumen_epithelium_Selenomonadaceae$Abundance, rumen_epithelium_Selenomonadaceae$treatment) #NS (0.161)
#Veillonellaceae
rumen_epithelium_Veillonellaceae <- rumen_epithelium_family_melt[which(rumen_epithelium_family_melt$Family=="Veillonellaceae"),]
kruskal.test(rumen_epithelium_Veillonellaceae$Abundance, rumen_epithelium_Veillonellaceae$treatment) # sig. 0.04
#Atopobiaceae
rumen_epithelium_Atopobiaceae <- rumen_epithelium_family_melt[which(rumen_epithelium_family_melt$Family=="Atopobiaceae"),]
kruskal.test(rumen_epithelium_Atopobiaceae$Abundance, rumen_epithelium_Atopobiaceae$treatment) #NS 0.855
#Oscillospiraceae
rumen_epithelium_Oscillospiraceae <- rumen_epithelium_family_melt[which(rumen_epithelium_family_melt$Family=="Oscillospiraceae"),]
kruskal.test(rumen_epithelium_Oscillospiraceae$Abundance, rumen_epithelium_Oscillospiraceae$treatment) # sig. 0.02
#Ruminococcaceae
rumen_epithelium_Ruminococcaceae <- rumen_epithelium_family_melt[which(rumen_epithelium_family_melt$Family=="Ruminococcaceae"),]
kruskal.test(rumen_epithelium_Ruminococcaceae$Abundance, rumen_epithelium_Ruminococcaceae$treatment) # sig. 0.048
#Muribaculaceae
rumen_epithelium_Muribaculaceae <- rumen_epithelium_family_melt[which(rumen_epithelium_family_melt$Family=="Muribaculaceae"),]
kruskal.test(rumen_epithelium_Muribaculaceae$Abundance, rumen_epithelium_Muribaculaceae$treatment) # NS (0.273)
#Anaerovoracaceae
rumen_epithelium_Anaerovoracaceae <- rumen_epithelium_family_melt[which(rumen_epithelium_family_melt$Family=="Anaerovoracaceae"),]
kruskal.test(rumen_epithelium_Anaerovoracaceae$Abundance, rumen_epithelium_Anaerovoracaceae$treatment) # NS (0.051)
#Acidaminococcaceae
rumen_epithelium_Acidaminococcaceae <- rumen_epithelium_family_melt[which(rumen_epithelium_family_melt$Family=="Acidaminococcaceae"),]
kruskal.test(rumen_epithelium_Acidaminococcaceae$Abundance, rumen_epithelium_Acidaminococcaceae$treatment) # NS (0.273)
#Rikenellaceae
rumen_epithelium_Rikenellaceae<- rumen_epithelium_family_melt[which(rumen_epithelium_family_melt$Family=="Rikenellaceae"),]
kruskal.test(rumen_epithelium_Rikenellaceae$Abundance, rumen_epithelium_Rikenellaceae$treatment) # sig (0.035)
#Methanobacteriaceae
rumen_epithelium_Methanobacteriaceae <- rumen_epithelium_family_melt[which(rumen_epithelium_family_melt$Family=="Methanobacteriaceae"),]
kruskal.test(rumen_epithelium_Methanobacteriaceae$Abundance, rumen_epithelium_Methanobacteriaceae$treatment) # NS 0.190
#Clostridia UCG-014
rumen_epithelium_Clostridia_UCG014 <- rumen_epithelium_family_melt[which(rumen_epithelium_family_melt$Family=="Clostridia_UCG-014"),]
kruskal.test(rumen_epithelium_Clostridia_UCG014$Abundance, rumen_epithelium_Clostridia_UCG014$treatment) # sig (0.015)



### SM. INTESTINE ####
### LUMEN 
sm_intestine_lumen.hclust <- hclust(sm_intestine_lumen.dist, method = "ward.D2")
sm_intestine_lumen.dendro <- as.dendrogram(sm_intestine_lumen.hclust)
sm_intestine_lumen.dendro.data <- dendro_data(sm_intestine_lumen.dendro, type = "rectangle")
sm_intestine_lumen.metadata_for_dendro <- as_tibble(sm_intestine_lumen.css@sam_data)
sm_intestine_lumen.dendro.data$labels <- sm_intestine_lumen.dendro.data$labels %>%
  left_join(sm_intestine_lumen.metadata_for_dendro, by = c("label" = "Sample"))

ggplot(sm_intestine_lumen.dendro.data$segments) +
  theme_minimal() +
  labs(y= "Ward's Distance") +
  geom_segment(aes(x=x,y=y,xend=xend,yend=yend)) +
  geom_point(data = sm_intestine_lumen.dendro.data$labels, 
             aes(x=x,y=y, colour = trt_LAtype_specifci_loc, fill= trt_LAtype_specifci_loc),
             size = 12.5, shape=22, stroke =1.5, position = position_nudge(y=-0.06)) +
  scale_y_continuous(limits = c(-.1,0.4)) +
  scale_x_discrete(expand = c(0.03,0,0.03,0)) +
  #scale_fill_manual(values = alpha(tylosin_palette,0.5)) +
  #scale_colour_manual(values = tylosin_palette) +
  scale_fill_manual(values = c("mediumorchid1","mediumorchid4","darkorchid4","goldenrod1","goldenrod3","darkgoldenrod4")) +
  scale_colour_manual(values = c("mediumorchid1","mediumorchid4","darkorchid4","goldenrod1","goldenrod3","darkgoldenrod4")) +
  theme(#legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.line.y = element_line(size = 0.7, colour = "black"),
        axis.ticks.y = element_line(size = 0.75, colour = "black"),
        axis.title.y = element_text(size = 24),
        axis.text.y = element_text(size = 14, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())

# family RA plot for under dendro
sm_intestine_lumen_family_filt_melt <- sm_intestine_lumen_family_melt %>%
  group_by(Family) %>%
  mutate(mean_family_ra = mean(Abundance))

sm_intestine_lumen_family_filt_melt$mean_family_ra <- as.numeric(sm_intestine_lumen_family_filt_melt$mean_family_ra)
sm_intestine_lumen_family_filt_melt$Family[sm_intestine_lumen_family_filt_melt$mean_family_ra < 0.01] <- "zzzzOther"
length(unique(sm_intestine_lumen_family_filt_melt$Family)) #46 families > 0.01%
write.csv(sm_intestine_lumen_family_filt_melt$Family,"sm_intestine_lumen_families.csv")
write.csv(sm_intestine_lumen_family_filt_melt$mean_family_ra,"mean_SI.csv")
sm_intestine_lumen_family_palette <-c("#A4AE33",	"#9FDC5B",	"#F0E9DE",	"goldenrod2",	"#48CFCA",	"#E7C57D",	"#9CF5E2",	"#466D92",	"#90B4F2",	"#A3748D",	"#E9906D",	"#B6F56A",	"#593050",	"#DCE080",	"#4DA176",	"#F0D3B0",	"#F0EDA4",	"#8FB77D",	"#DE9B37",	"#7E589A",	"darkgoldenrod3",	"#6F97E9",	"#5CF4F9",	"springgreen3",	"#CD90DF",	"#EE61C9",	"bisque3",	"#CFF685",	"#9E82EA",	"darkorange1",	"#D08992",	"#F3B2E3",	"#8FF4B6",	"dodgerblue3",	"#8DE9EA",	"firebrick2",	"#DEB2F3",	"#E6F132",	"chocolate3",	"#EA3C5D",	"hotpink",	"#D9DBEF",	"#5FB0B1",	"chartreuse2",	"#C5D356",	"grey87")
sm_intestine_lumen_dendro_sample_order <- sm_intestine_lumen.dendro.data$labels$label

ggplot(sm_intestine_lumen_family_filt_melt, aes(x= sample_Sample, y= Abundance, fill= Family)) +
  theme_minimal() +
  labs(y= "Relative Abundance (%)") +
  geom_bar(stat = "identity", colour = "black") +
  scale_x_discrete(limits = sm_intestine_lumen_dendro_sample_order, expand = c(0.03,0,0.015,0)) +
  scale_fill_manual(values = sm_intestine_lumen_family_palette) +
  theme(legend.position = "none",
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.line.y = element_line(size = 0.7, colour = "black"),
        axis.ticks.y = element_line(colour = "black", size = 0.75),
        axis.title.y = element_text(size = 24),
        axis.text.y = element_text(size = 14, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())

## stats
#Lachnospiraceae
sm_intestine_lumen_Lachnospiraceae <- sm_intestine_lumen_family_melt[which(sm_intestine_lumen_family_melt$Family=="Lachnospiraceae"),]
kruskal.test(sm_intestine_lumen_Lachnospiraceae$Abundance, sm_intestine_lumen_Lachnospiraceae$treatment) # NS (0.843)
#Peptostreptococcaceae
sm_intestine_lumen_Peptostreptococcaceae <- sm_intestine_lumen_family_melt[which(sm_intestine_lumen_family_melt$Family=="Peptostreptococcaceae"),]
kruskal.test(sm_intestine_lumen_Peptostreptococcaceae$Abundance, sm_intestine_lumen_Peptostreptococcaceae$treatment) # NS (0.058)
#Clostridiaceae
sm_intestine_lumen_Clostridiaceae <- sm_intestine_lumen_family_melt[which(sm_intestine_lumen_family_melt$Family=="Clostridiaceae"),]
kruskal.test(sm_intestine_lumen_Clostridiaceae$Abundance, sm_intestine_lumen_Clostridiaceae$treatment) # NS (0.477)
#Atopobiaceae
sm_intestine_lumen_Atopobiaceae <- sm_intestine_lumen_family_melt[which(sm_intestine_lumen_family_melt$Family=="Atopobiaceae"),]
kruskal.test(sm_intestine_lumen_Atopobiaceae$Abundance, sm_intestine_lumen_Atopobiaceae$treatment) #NS (0.286)
#Erysipelotrichaceae
sm_intestine_lumen_Erysipelotrichaceae <- sm_intestine_lumen_family_melt[which(sm_intestine_lumen_family_melt$Family=="Erysipelotrichaceae"),]
kruskal.test(sm_intestine_lumen_Erysipelotrichaceae$Abundance, sm_intestine_lumen_Erysipelotrichaceae$treatment) # NS (0.429)
#Methanobacteriaceae
sm_intestine_lumen_Methanobacteriaceae <- sm_intestine_lumen_family_melt[which(sm_intestine_lumen_family_melt$Family=="Methanobacteriaceae"),]
kruskal.test(sm_intestine_lumen_Methanobacteriaceae$Abundance, sm_intestine_lumen_Methanobacteriaceae$treatment) # sig. 0.005
#Bifidobacteriaceae
sm_intestine_lumen_Bifidobacteriaceae <- sm_intestine_lumen_family_melt[which(sm_intestine_lumen_family_melt$Family=="Bifidobacteriaceae"),]
kruskal.test(sm_intestine_lumen_Bifidobacteriaceae$Abundance, sm_intestine_lumen_Bifidobacteriaceae$treatment) # NS 0.843
#Enterobacteriaceae
sm_intestine_lumen_Enterobacteriaceae <- sm_intestine_lumen_family_melt[which(sm_intestine_lumen_family_melt$Family=="Enterobacteriaceae"),]
kruskal.test(sm_intestine_lumen_Enterobacteriaceae$Abundance, sm_intestine_lumen_Enterobacteriaceae$treatment) # NS (0.357)
#Ruminococcaceae
sm_intestine_lumen_Ruminococcaceae <- sm_intestine_lumen_family_melt[which(sm_intestine_lumen_family_melt$Family=="Ruminococcaceae"),]
kruskal.test(sm_intestine_lumen_Ruminococcaceae$Abundance, sm_intestine_lumen_Ruminococcaceae$treatment) # sig. 0.018
#Anaerovoracaceae
sm_intestine_lumen_Anaerovoracaceae <- sm_intestine_lumen_family_melt[which(sm_intestine_lumen_family_melt$Family=="Anaerovoracaceae"),]
kruskal.test(sm_intestine_lumen_Anaerovoracaceae$Abundance, sm_intestine_lumen_Anaerovoracaceae$treatment) # sig. 0.003


#### SM. INTESTINE EPITHELIUM
sm_intestine_epithelium.hclust <- hclust(sm_intestine_epithelium.dist, method = "ward.D2")
sm_intestine_epithelium.dendro <- as.dendrogram(sm_intestine_epithelium.hclust)
sm_intestine_epithelium.dendro.data <- dendro_data(sm_intestine_epithelium.dendro, type = "rectangle")
sm_intestine_epithelium.metadata_for_dendro <- as_tibble(sm_intestine_epithelium.css@sam_data)
sm_intestine_epithelium.dendro.data$labels <- sm_intestine_epithelium.dendro.data$labels %>%
  left_join(sm_intestine_epithelium.metadata_for_dendro, by = c("label" = "Sample"))

ggplot(sm_intestine_epithelium.dendro.data$segments) +
  theme_minimal() +
  labs(y= "Ward's Distance") +
  geom_segment(aes(x=x,y=y,xend=xend,yend=yend)) +
  geom_point(data = sm_intestine_epithelium.dendro.data$labels, 
             aes(x=x,y=y, colour = treatment, fill= treatment),
             size = 13, shape=22, stroke =1.5, position = position_nudge(y=-0.24)) +
  scale_y_continuous(limits = c(-.35,1.7)) +
  scale_x_discrete(expand = c(0.03,0,0.03,0)) +
  scale_fill_manual(values = alpha(tylosin_palette,0.5)) +
  scale_colour_manual(values = tylosin_palette) +
  theme(legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.line.y = element_line(size = 0.7, colour = "black"),
        axis.ticks.y = element_line(size = 0.75, colour = "black"),
        axis.title.y = element_text(size = 24),
        axis.text.y = element_text(size = 14, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())

# family RA plot for under dendro
sm_intestine_epithelium_family_filt_melt <- sm_intestine_epithelium_family_melt %>%
  group_by(Family) %>%
  mutate(mean_family_ra = mean(Abundance))

sm_intestine_epithelium_family_filt_melt$Family[sm_intestine_epithelium_family_filt_melt$mean_family_ra < 0.01] <- "zzzzOther"
length(unique(sm_intestine_epithelium_family_filt_melt$Family)) #52 families > 0.01%
write.csv(sm_intestine_epithelium_family_filt_melt$Family,"sm_intestine_epithelium_families.csv")
write.csv(sm_intestine_epithelium_family_filt_melt$mean_family_ra,"sm_intestine_epithelium_meanRA.csv")
sm_intestine_epithelium_family_palette <- c("#A4AE33",	"#9FDC5B",	"#F0E9DE",	"goldenrod2",	"#48CFCA",	"#E7C57D",	"#9CF5E2",	"springgreen4",	"#466D92",	"#90B4F2",	"#A3748D",	"#E9906D",	"#B6F56A",	"#593050",	"#DCE080",	"#4DA176",	"#F0D3B0",	"#F0EDA4",	"#8FB77D",	"#DE9B37",	"#7E589A",	"darkgoldenrod3",	"#6F97E9",	"#5CF4F9",	"springgreen3",	"#CD90DF",	"#EE61C9",	"bisque3",	"#AC8B80",	"#CFF685",	"#B1E8D4",	"#9ED6DB",	"#9E82EA",	"#3DE96D",	"darkorange1",	"#D08992",	"#F3B2E3",	"dodgerblue3",	"#8DE9EA",	"firebrick2",	"#DEB2F3",	"#E6F132",	"chocolate3",	"#5771DE",	"#EA3C5D",	"hotpink",	"#6BA7C8",	"#D9DBEF",	"#5FB0B1",	"chartreuse2",	"#A2B9B0",	"grey87")
sm_intestine_epithelium_dendro_sample_order <- sm_intestine_epithelium.dendro.data$labels$label

ggplot(sm_intestine_epithelium_family_filt_melt, aes(x= sample_Sample, y= Abundance, fill= Family)) +
  theme_minimal() +
  labs(y= "Relative Abundance (%)") +
  geom_bar(stat = "identity", colour = "black") +
  scale_x_discrete(limits = sm_intestine_epithelium_dendro_sample_order, expand = c(0.03,0,0.015,0)) +
  scale_fill_manual(values = sm_intestine_epithelium_family_palette) +
  theme(legend.position = "none",
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.line.y = element_line(size = 0.7, colour = "black"),
        axis.ticks.y = element_line(colour = "black", size = 0.75),
        axis.title.y = element_text(size = 24),
        axis.text.y = element_text(size = 14, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())

## stats
#Lachnospiraceae
sm_intestine_epithelium_Lachnospiraceae <- sm_intestine_epithelium_family_melt[which(sm_intestine_epithelium_family_melt$Family=="Lachnospiraceae"),]
kruskal.test(sm_intestine_epithelium_Lachnospiraceae$Abundance, sm_intestine_epithelium_Lachnospiraceae$treatment) # NS (0.618)
#Peptostreptococcaceae
sm_intestine_epithelium_Peptostreptococcaceae <- sm_intestine_epithelium_family_melt[which(sm_intestine_epithelium_family_melt$Family=="Peptostreptococcaceae"),]
kruskal.test(sm_intestine_epithelium_Peptostreptococcaceae$Abundance, sm_intestine_epithelium_Peptostreptococcaceae$treatment) # NS (0.868)
#Atopobiaceae
sm_intestine_epithelium_Atopobiaceae <- sm_intestine_epithelium_family_melt[which(sm_intestine_epithelium_family_melt$Family=="Atopobiaceae"),]
kruskal.test(sm_intestine_epithelium_Atopobiaceae$Abundance, sm_intestine_epithelium_Atopobiaceae$treatment) #NS (0.934)
#Clostridiaceae
sm_intestine_epithelium_Clostridiaceae <- sm_intestine_epithelium_family_melt[which(sm_intestine_epithelium_family_melt$Family=="Clostridiaceae"),]
kruskal.test(sm_intestine_epithelium_Clostridiaceae$Abundance, sm_intestine_epithelium_Clostridiaceae$treatment) # NS (0.803)
#Erysipelotrichaceae
sm_intestine_epithelium_Erysipelotrichaceae <- sm_intestine_epithelium_family_melt[which(sm_intestine_epithelium_family_melt$Family=="Erysipelotrichaceae"),]
kruskal.test(sm_intestine_epithelium_Erysipelotrichaceae$Abundance, sm_intestine_epithelium_Erysipelotrichaceae$treatment) # NS (0.708)
#Ruminococcaceae
sm_intestine_epithelium_Ruminococcaceae <- sm_intestine_epithelium_family_melt[which(sm_intestine_epithelium_family_melt$Family=="Ruminococcaceae"),]
kruskal.test(sm_intestine_epithelium_Ruminococcaceae$Abundance, sm_intestine_epithelium_Ruminococcaceae$treatment) # sig. 0.006
#Methanobacteriaceae
sm_intestine_epithelium_Methanobacteriaceae <- sm_intestine_epithelium_family_melt[which(sm_intestine_epithelium_family_melt$Family=="Methanobacteriaceae"),]
kruskal.test(sm_intestine_epithelium_Methanobacteriaceae$Abundance, sm_intestine_epithelium_Methanobacteriaceae$treatment) # NS 0.803
#Anaerovoracaceae
sm_intestine_epithelium_Anaerovoracaceae <- sm_intestine_epithelium_family_melt[which(sm_intestine_epithelium_family_melt$Family=="Anaerovoracaceae"),]
kruskal.test(sm_intestine_epithelium_Anaerovoracaceae$Abundance, sm_intestine_epithelium_Anaerovoracaceae$treatment) # NS (0.051)
#Bifidobacteriaceae
sm_intestine_epithelium_Bifidobacteriaceae <- sm_intestine_epithelium_family_melt[which(sm_intestine_epithelium_family_melt$Family=="Bifidobacteriaceae"),]
kruskal.test(sm_intestine_epithelium_Bifidobacteriaceae$Abundance, sm_intestine_epithelium_Bifidobacteriaceae$treatment) # NS 0.061


### LG INTESTINE ####
## LUMEN
lg_intestine_lumen.hclust <- hclust(lg_intestine_lumen.dist, method = "ward.D2")
lg_intestine_lumen.dendro <- as.dendrogram(lg_intestine_lumen.hclust)
lg_intestine_lumen.dendro.data <- dendro_data(lg_intestine_lumen.dendro, type = "rectangle")
lg_intestine_lumen.metadata_for_dendro <- as_tibble(lg_intestine_lumen.css@sam_data)
lg_intestine_lumen.dendro.data$labels <- lg_intestine_lumen.dendro.data$labels %>%
  left_join(lg_intestine_lumen.metadata_for_dendro, by = c("label" = "Sample"))

ggplot(lg_intestine_lumen.dendro.data$segments) +
  theme_minimal() +
  labs(y= "Ward's Distance") +
  geom_segment(aes(x=x,y=y,xend=xend,yend=yend)) +
  geom_point(data = lg_intestine_lumen.dendro.data$labels, 
             aes(x=x,y=y, colour = treatment, fill= treatment),
             size = 10.5, shape=22, stroke =1.5, position = position_nudge(y=-0.18)) +
  scale_y_continuous(limits = c(-.25,1.27)) +
  scale_x_discrete(expand = c(0.03,0,0.03,0)) +
  scale_fill_manual(values = alpha(tylosin_palette,0.5)) +
  scale_colour_manual(values = tylosin_palette) +
  theme(legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.line.y = element_line(size = 0.7, colour = "black"),
        axis.ticks.y = element_line(size = 0.75, colour = "black"),
        axis.title.y = element_text(size = 24),
        axis.text.y = element_text(size = 14, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())

# family RA plot for under dendro
lg_intestine_lumen_family_filt_melt <- lg_intestine_lumen_family_melt %>%
  group_by(Family) %>%
  mutate(mean_family_ra = mean(Abundance))

lg_intestine_lumen_family_filt_melt$Family[lg_intestine_lumen_family_filt_melt$mean_family_ra < 0.01] <- "zzzzOther"
length(unique(lg_intestine_lumen_family_filt_melt$Family)) #60 families > 0.01%
write.csv(lg_intestine_lumen_family_filt_melt$Family,"lg_intestine_lumen_families.csv")
write.csv(lg_intestine_lumen_family_filt_melt$mean_family_ra,"lg_intestine_lumen_meanRA.csv")
lg_intestine_lumen_family_palette <-c("#9FDC5B",	"#42EDC1",	"#F0E9DE",	"#78CE38",	"#AF4B71",	"#48CFCA",	"#E7C57D",	"#9CF5E2",	"springgreen4",	"palegreen2",	"#9CD62D",	"#466D92",	"#3084D5",	"#90B4F2",	"#A3748D",	"#A1D676",	"#E9906D",	"#593050",	"#DCE080",	"#4DA176",	"#F0D3B0",	"#F0EDA4",	"#8FB77D",	"#7E589A",	"darkgoldenrod3",	"#6F97E9",	"yellow4",	"#ACD0EF",	"darkorchid4",	"#5CF4F9",	"#EB6393",	"springgreen3",	"#CD90DF",	"#B19DEF",	"#CFF685",	"#9ED6DB",	"darkorange1",	"#D8EEEE",	"#3DE49F",	"#D08992",	"#F3B2E3",	"#8FF4B6",	"dodgerblue3",	"#8DE9EA",	"mediumorchid",	"firebrick2",	"#DEB2F3",	"#E6F132",	"chocolate3",	"#BCC0EA",	"#EA3C5D",	"hotpink",	"#8761E0",	"#90EB8C",	"#EE3DD4",	"honeydew3",	"#C76ED7",	"#EC592B",	"chartreuse2",	"grey87")
lg_intestine_lumen_dendro_sample_order <- lg_intestine_lumen.dendro.data$labels$label

ggplot(lg_intestine_lumen_family_filt_melt, aes(x= sample_Sample, y= Abundance, fill= Family)) +
  theme_minimal() +
  labs(y= "Relative Abundance (%)") +
  geom_bar(stat = "identity", colour = "black") +
  scale_x_discrete(limits = lg_intestine_lumen_dendro_sample_order, expand = c(0.03,0,0.015,0)) +
  scale_fill_manual(values = lg_intestine_lumen_family_palette) +
  theme(#legend.position = "none",
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.line.y = element_line(size = 0.7, colour = "black"),
        axis.ticks.y = element_line(colour = "black", size = 0.75),
        axis.title.y = element_text(size = 24),
        axis.text.y = element_text(size = 14, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())

#Lachnospiraceae
lg_intestine_lumen_Lachnospiraceae <- lg_intestine_lumen_family_melt[which(lg_intestine_lumen_family_melt$Family=="Lachnospiraceae"),]
kruskal.test(lg_intestine_lumen_Lachnospiraceae$Abundance, lg_intestine_lumen_Lachnospiraceae$treatment) # NS (0.862)
#Peptostreptococcaceae
lg_intestine_lumen_Peptostreptococcaceae <- lg_intestine_lumen_family_melt[which(lg_intestine_lumen_family_melt$Family=="Peptostreptococcaceae"),]
kruskal.test(lg_intestine_lumen_Peptostreptococcaceae$Abundance, lg_intestine_lumen_Peptostreptococcaceae$treatment) # NS (0.168)
#Oscillospiraceae
lg_intestine_lumen_Oscillospiraceae <- lg_intestine_lumen_family_melt[which(lg_intestine_lumen_family_melt$Family=="Oscillospiraceae"),]
kruskal.test(lg_intestine_lumen_Oscillospiraceae$Abundance, lg_intestine_lumen_Oscillospiraceae$treatment) # NS 0.159
#Prevotellaceae
lg_intestine_lumen_Prevotellaceae <- lg_intestine_lumen_family_melt[which(lg_intestine_lumen_family_melt$Family=="Prevotellaceae"),]
kruskal.test(lg_intestine_lumen_Prevotellaceae$Abundance, lg_intestine_lumen_Prevotellaceae$treatment) # NS 0.199
#Bacteroidaceae
lg_intestine_lumen_Bacteroidaceae <- lg_intestine_lumen_family_melt[which(lg_intestine_lumen_family_melt$Family=="Bacteroidaceae"),]
kruskal.test(lg_intestine_lumen_Bacteroidaceae$Abundance, lg_intestine_lumen_Bacteroidaceae$treatment) #NS 0.103
#Clostridiaceae
lg_intestine_lumen_Clostridiaceae <- lg_intestine_lumen_family_melt[which(lg_intestine_lumen_family_melt$Family=="Clostridiaceae"),]
kruskal.test(lg_intestine_lumen_Clostridiaceae$Abundance, lg_intestine_lumen_Clostridiaceae$treatment) #NS 0.103
#Erysipelotrichaceae
lg_intestine_lumen_Erysipelotrichaceae <- lg_intestine_lumen_family_melt[which(lg_intestine_lumen_family_melt$Family=="Erysipelotrichaceae"),]
kruskal.test(lg_intestine_lumen_Erysipelotrichaceae$Abundance, lg_intestine_lumen_Erysipelotrichaceae$treatment) # sig. 0.003
#Rikenellaceae
lg_intestine_lumen_Rikenellaceae<- lg_intestine_lumen_family_melt[which(lg_intestine_lumen_family_melt$Family=="Rikenellaceae"),]
kruskal.test(lg_intestine_lumen_Rikenellaceae$Abundance, lg_intestine_lumen_Rikenellaceae$treatment) # NS (0.962)
#Muribaculaceae
lg_intestine_lumen_Muribaculaceae <- lg_intestine_lumen_family_melt[which(lg_intestine_lumen_family_melt$Family=="Muribaculaceae"),]
kruskal.test(lg_intestine_lumen_Muribaculaceae$Abundance, lg_intestine_lumen_Muribaculaceae$treatment) # NS (0.887)
#Atopobiaceae
lg_intestine_lumen_Atopobiaceae <- lg_intestine_lumen_family_melt[which(lg_intestine_lumen_family_melt$Family=="Atopobiaceae"),]
kruskal.test(lg_intestine_lumen_Atopobiaceae$Abundance, lg_intestine_lumen_Atopobiaceae$treatment) #sig 0.006
#Ruminococcaceae
lg_intestine_lumen_Ruminococcaceae <- lg_intestine_lumen_family_melt[which(lg_intestine_lumen_family_melt$Family=="Ruminococcaceae"),]
kruskal.test(lg_intestine_lumen_Ruminococcaceae$Abundance, lg_intestine_lumen_Ruminococcaceae$treatment) # NS (0.261)
#Eubacterium_coprostanoligenes_group
lg_intestine_lumen_Eubacterium_coprostanoligenes_group <- lg_intestine_lumen_family_melt[which(lg_intestine_lumen_family_melt$Family=="[Eubacterium]_coprostanoligenes_group"),]
kruskal.test(lg_intestine_lumen_Eubacterium_coprostanoligenes_group$Abundance, lg_intestine_lumen_Eubacterium_coprostanoligenes_group$treatment) # NS (0.646)
#Spirochaetaceae
lg_intestine_lumen_Spirochaetaceae <- lg_intestine_lumen_family_melt[which(lg_intestine_lumen_family_melt$Family=="Spirochaetaceae"),]
kruskal.test(lg_intestine_lumen_Spirochaetaceae$Abundance, lg_intestine_lumen_Spirochaetaceae$treatment) # NS 0.962
#Methanobacteriaceae
lg_intestine_lumen_Methanobacteriaceae <- lg_intestine_lumen_family_melt[which(lg_intestine_lumen_family_melt$Family=="Methanobacteriaceae"),]
kruskal.test(lg_intestine_lumen_Methanobacteriaceae$Abundance, lg_intestine_lumen_Methanobacteriaceae$treatment) # sig. 0.044
#Anaerovoracaceae
lg_intestine_lumen_Anaerovoracaceae <- lg_intestine_lumen_family_melt[which(lg_intestine_lumen_family_melt$Family=="Anaerovoracaceae"),]
kruskal.test(lg_intestine_lumen_Anaerovoracaceae$Abundance, lg_intestine_lumen_Anaerovoracaceae$treatment) # sig. 0.030
#Oscillospirales UCG-010
lg_intestine_lumen_Oscillospirales_UCG010 <- lg_intestine_lumen_family_melt[which(lg_intestine_lumen_family_melt$Family=="UCG-010"),]
kruskal.test(lg_intestine_lumen_Oscillospirales_UCG010$Abundance, lg_intestine_lumen_Oscillospirales_UCG010$treatment) # NS (0.669)
#Christensenellaceae
lg_intestine_lumen_Christensenellaceae <- lg_intestine_lumen_family_melt[which(lg_intestine_lumen_family_melt$Family=="Christensenellaceae"),]
kruskal.test(lg_intestine_lumen_Christensenellaceae$Abundance, lg_intestine_lumen_Christensenellaceae$treatment) # NS (0.516)
#Monoglobaceae
lg_intestine_lumen_Monoglobaceae <- lg_intestine_lumen_family_melt[which(lg_intestine_lumen_family_melt$Family=="Monoglobaceae"),]
kruskal.test(lg_intestine_lumen_Monoglobaceae$Abundance, lg_intestine_lumen_Monoglobaceae$treatment) # NS (0.288)

#### LG INTESTINE EPITHELIUM
lg_intestine_epithelium.hclust <- hclust(lg_intestine_epithelium.dist, method = "ward.D2")
lg_intestine_epithelium.dendro <- as.dendrogram(lg_intestine_epithelium.hclust)
lg_intestine_epithelium.dendro.data <- dendro_data(lg_intestine_epithelium.dendro, type = "rectangle")
lg_intestine_epithelium.metadata_for_dendro <- as_tibble(lg_intestine_epithelium.css@sam_data)
lg_intestine_epithelium.dendro.data$labels <- lg_intestine_epithelium.dendro.data$labels %>%
  left_join(lg_intestine_epithelium.metadata_for_dendro, by = c("label" = "Sample"))

ggplot(lg_intestine_epithelium.dendro.data$segments) +
  theme_minimal() +
  labs(y= "Ward's Distance") +
  geom_segment(aes(x=x,y=y,xend=xend,yend=yend)) +
  geom_point(data = lg_intestine_epithelium.dendro.data$labels, 
             aes(x=x,y=y, colour = treatment, fill= treatment),
             size = 10, shape=22, stroke =1.5, position = position_nudge(y=-0.18)) +
  scale_y_continuous(limits = c(-.25,1.27)) +
  scale_x_discrete(expand = c(0.03,0,0.03,0)) +
  scale_fill_manual(values = alpha(tylosin_palette,0.5)) +
  scale_colour_manual(values = tylosin_palette) +
  theme(legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.line.y = element_line(size = 0.7, colour = "black"),
        axis.ticks.y = element_line(size = 0.75, colour = "black"),
        axis.title.y = element_text(size = 24),
        axis.text.y = element_text(size = 14, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())

# family RA plot for under dendro
lg_intestine_epithelium_family_filt_melt <- lg_intestine_epithelium_family_melt %>%
  group_by(Family) %>%
  mutate(mean_family_ra = mean(Abundance))

lg_intestine_epithelium_family_filt_melt$Family[lg_intestine_epithelium_family_filt_melt$mean_family_ra < 0.01] <- "zzzzOther"
length(unique(lg_intestine_epithelium_family_filt_melt$Family)) #64 families > 0.01%
write.csv(lg_intestine_epithelium_family_filt_melt$Family,"lg_intestine_epithelium_families.csv")
write.csv(lg_intestine_epithelium_family_filt_melt$mean_family_ra,"lg_intestine_epithelium_meanRA.csv")
lg_intestine_epithelium_family_palette <- c("#9FDC5B",	"#E7CB30",	"#42EDC1",	"#F0E9DE",	"#78CE38",	"#AF4B71",	"#48CFCA",	"#E7C57D",	"#9CF5E2",	"springgreen4",	"palegreen2",	"#9CD62D",	"#466D92",	"#3084D5",	"#90B4F2",	"#A3748D",	"#A1D676",	"#E9906D",	"#593050",	"#DCE080",	"#4DA176",	"#F0D3B0",	"#F0EDA4",	"#8FB77D",	"#7E589A",	"darkgoldenrod3",	"#6F97E9",	"yellow4",	"#ACD0EF",	"#5CF4F9",	"#4D4DE3",	"#EB6393",	"springgreen3",	"#CD90DF",	"bisque3",	"#CFF685",	"#9ED6DB",	"#9E82EA",	"darkorange1",	"#D8EEEE",	"#3DE49F",	"#D08992",	"#F3B2E3",	"#8FF4B6",	"dodgerblue3",	"#8DE9EA",	"mediumorchid",	"firebrick2",	"#DEB2F3",	"#E6F132",	"chocolate3",	"#BCC0EA",	"hotpink",	"#8761E0",	"#90EB8C",	"#EE3DD4",	"honeydew3",	"#E8CFF1",	"#C76ED7",	"#5FB0B1",	"#EC592B",	"chartreuse2",	"#71C9B2",	"grey87")
lg_intestine_epithelium_dendro_sample_order <- lg_intestine_epithelium.dendro.data$labels$label

ggplot(lg_intestine_epithelium_family_filt_melt, aes(x= sample_Sample, y= Abundance, fill= Family)) +
  theme_minimal() +
  labs(y= "Relative Abundance (%)") +
  geom_bar(stat = "identity", colour = "black") +
  scale_x_discrete(limits = lg_intestine_epithelium_dendro_sample_order, expand = c(0.03,0,0.015,0)) +
  scale_fill_manual(values = lg_intestine_epithelium_family_palette) +
  theme(legend.position = "none",
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.line.y = element_line(size = 0.7, colour = "black"),
    axis.ticks.y = element_line(colour = "black", size = 0.75),
    axis.title.y = element_text(size = 24),
    axis.text.y = element_text(size = 14, colour = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_blank())

# stats
#Lachnospiraceae
lg_intestine_epithelium_Lachnospiraceae <- lg_intestine_epithelium_family_melt[which(lg_intestine_epithelium_family_melt$Family=="Lachnospiraceae"),]
kruskal.test(lg_intestine_epithelium_Lachnospiraceae$Abundance, lg_intestine_epithelium_Lachnospiraceae$treatment) # NS (0.468)
#Prevotellaceae
lg_intestine_epithelium_Prevotellaceae <- lg_intestine_epithelium_family_melt[which(lg_intestine_epithelium_family_melt$Family=="Prevotellaceae"),]
kruskal.test(lg_intestine_epithelium_Prevotellaceae$Abundance, lg_intestine_epithelium_Prevotellaceae$treatment) # NS 0.355
#Oscillospiraceae
lg_intestine_epithelium_Oscillospiraceae <- lg_intestine_epithelium_family_melt[which(lg_intestine_epithelium_family_melt$Family=="Oscillospiraceae"),]
kruskal.test(lg_intestine_epithelium_Oscillospiraceae$Abundance, lg_intestine_epithelium_Oscillospiraceae$treatment) # NS 0.488
#Bacteroidaceae
lg_intestine_epithelium_Bacteroidaceae <- lg_intestine_epithelium_family_melt[which(lg_intestine_epithelium_family_melt$Family=="Bacteroidaceae"),]
kruskal.test(lg_intestine_epithelium_Bacteroidaceae$Abundance, lg_intestine_epithelium_Bacteroidaceae$treatment) #NS 0.391
#Peptostreptococcaceae
lg_intestine_epithelium_Peptostreptococcaceae <- lg_intestine_epithelium_family_melt[which(lg_intestine_epithelium_family_melt$Family=="Peptostreptococcaceae"),]
kruskal.test(lg_intestine_epithelium_Peptostreptococcaceae$Abundance, lg_intestine_epithelium_Peptostreptococcaceae$treatment) # NS (0.262)
#Clostridiaceae
lg_intestine_epithelium_Clostridiaceae <- lg_intestine_epithelium_family_melt[which(lg_intestine_epithelium_family_melt$Family=="Clostridiaceae"),]
kruskal.test(lg_intestine_epithelium_Clostridiaceae$Abundance, lg_intestine_epithelium_Clostridiaceae$treatment) #NS 0.080
#Erysipelotrichaceae
lg_intestine_epithelium_Erysipelotrichaceae <- lg_intestine_epithelium_family_melt[which(lg_intestine_epithelium_family_melt$Family=="Erysipelotrichaceae"),]
kruskal.test(lg_intestine_epithelium_Erysipelotrichaceae$Abundance, lg_intestine_epithelium_Erysipelotrichaceae$treatment) # sig. 0.008
#Rikenellaceae
lg_intestine_epithelium_Rikenellaceae<- lg_intestine_epithelium_family_melt[which(lg_intestine_epithelium_family_melt$Family=="Rikenellaceae"),]
kruskal.test(lg_intestine_epithelium_Rikenellaceae$Abundance, lg_intestine_epithelium_Rikenellaceae$treatment) # NS (0.947)
#Ruminococcaceae
lg_intestine_epithelium_Ruminococcaceae <- lg_intestine_epithelium_family_melt[which(lg_intestine_epithelium_family_melt$Family=="Ruminococcaceae"),]
kruskal.test(lg_intestine_epithelium_Ruminococcaceae$Abundance, lg_intestine_epithelium_Ruminococcaceae$treatment) # NS (0.553)
#Muribaculaceae
lg_intestine_epithelium_Muribaculaceae <- lg_intestine_epithelium_family_melt[which(lg_intestine_epithelium_family_melt$Family=="Muribaculaceae"),]
kruskal.test(lg_intestine_epithelium_Muribaculaceae$Abundance, lg_intestine_epithelium_Muribaculaceae$treatment) # NS (0.869)
#Atopobiaceae
lg_intestine_epithelium_Atopobiaceae <- lg_intestine_epithelium_family_melt[which(lg_intestine_epithelium_family_melt$Family=="Atopobiaceae"),]
kruskal.test(lg_intestine_epithelium_Atopobiaceae$Abundance, lg_intestine_epithelium_Atopobiaceae$treatment) #sig 0.011
#Eubacterium_coprostanoligenes_group
lg_intestine_epithelium_Eubacterium_coprostanoligenes_group <- lg_intestine_epithelium_family_melt[which(lg_intestine_epithelium_family_melt$Family=="[Eubacterium]_coprostanoligenes_group"),]
kruskal.test(lg_intestine_epithelium_Eubacterium_coprostanoligenes_group$Abundance, lg_intestine_epithelium_Eubacterium_coprostanoligenes_group$treatment) # NS (0.843)
#Oscillospirales UCG-010
lg_intestine_epithelium_Oscillospirales_UCG010 <- lg_intestine_epithelium_family_melt[which(lg_intestine_epithelium_family_melt$Family=="UCG-010"),]
kruskal.test(lg_intestine_epithelium_Oscillospirales_UCG010$Abundance, lg_intestine_epithelium_Oscillospirales_UCG010$treatment) # NS (0.921)
#Spirochaetaceae
lg_intestine_epithelium_Spirochaetaceae <- lg_intestine_epithelium_family_melt[which(lg_intestine_epithelium_family_melt$Family=="Spirochaetaceae"),]
kruskal.test(lg_intestine_epithelium_Spirochaetaceae$Abundance, lg_intestine_epithelium_Spirochaetaceae$treatment) # NS 0.817
#Anaerovoracaceae
lg_intestine_epithelium_Anaerovoracaceae <- lg_intestine_epithelium_family_melt[which(lg_intestine_epithelium_family_melt$Family=="Anaerovoracaceae"),]
kruskal.test(lg_intestine_epithelium_Anaerovoracaceae$Abundance, lg_intestine_epithelium_Anaerovoracaceae$treatment) # NS 0.766
#Methanobacteriaceae
lg_intestine_epithelium_Methanobacteriaceae <- lg_intestine_epithelium_family_melt[which(lg_intestine_epithelium_family_melt$Family=="Methanobacteriaceae"),]
kruskal.test(lg_intestine_epithelium_Methanobacteriaceae$Abundance, lg_intestine_epithelium_Methanobacteriaceae$treatment) # sig. 0.032
#Bifidobacteriaceae
lg_intestine_epithelium_Bifidobacteriaceae <- lg_intestine_epithelium_family_melt[which(lg_intestine_epithelium_family_melt$Family=="Bifidobacteriaceae"),]
kruskal.test(lg_intestine_epithelium_Bifidobacteriaceae$Abundance, lg_intestine_epithelium_Bifidobacteriaceae$treatment) # NS (0.478)
#Monoglobaceae
lg_intestine_epithelium_Monoglobaceae <- lg_intestine_epithelium_family_melt[which(lg_intestine_epithelium_family_melt$Family=="Monoglobaceae"),]
kruskal.test(lg_intestine_epithelium_Monoglobaceae$Abundance, lg_intestine_epithelium_Monoglobaceae$treatment) # NS (0.869)



################## Tylosin- DIFF ABUNDANT TAXA ####

### LACHNOSPIRACEAE ####
### RUMEN
## RUMEN LUMEN
rumen_lumen_lachnospiraceae_genus <- subset_taxa(rumen_lumen_genus, Family=="Lachnospiraceae")
rumen_lumen_lachnospiraceae_genus_filt <- merge_low_abundance(rumen_lumen_lachnospiraceae_genus, threshold = 0.1) # 13 taxa
write.csv(tax_table(rumen_lumen_lachnospiraceae_genus_filt), "rumen_lumen_lachno_taxa.csv")
write.csv(otu_table(rumen_lumen_lachnospiraceae_genus_filt), "rumen_lumen_lachno_otus.csv")
rumen_lumen_lachnospiraceae_genus_melt <- psmelt(rumen_lumen_lachnospiraceae_genus_filt)

rumen_lumen_lachno_pallete <- c("#A8D8A8",	"#00F0A8",	"#007848",	"#00C0A8",	"#78D878",	"#30C078",	"#30FF60",	"#186048",	"#C0C030",	"#309030",	"#00A890",	"#A8FFAB",	"grey88")
rumen_lumen_lachnospiraceae_family <- subset_taxa(rumen_lumen_family, Family=="Lachnospiraceae") %>%
  psmelt()

ggplot() + theme_bw() + labs(y= "Lachnospiraceae (RA%)", title = "LUMEN") +
  coord_flip() + theme_bw() +
  geom_bar(rumen_lumen_lachnospiraceae_genus_melt, colour = "black", linewidth = 0.75,
           mapping = aes(x= treatment, y= Abundance, fill=Genus), stat = "summary") +
  scale_fill_manual(values = rumen_lumen_lachno_pallete) +
  scale_x_discrete(limits = c("tylosin","control"), labels =c("TYLOSIN","NO TYLOSIN")) +
  scale_y_continuous(expand = c(0.001,0,0.1,0)) +
  geom_errorbar(rumen_lumen_lachnospiraceae_family, mapping= aes(x= treatment, y= Abundance),
                stat = "summary", width = 0.5, linewidth = 0.75) +
  theme(legend.position = "none",
        plot.title = element_text(size = 44),
        panel.border = element_rect(colour = "black", size = 1.25),
        panel.grid.major.y = element_blank(),
        axis.ticks = element_line(size = 0.75, colour= "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 34, colour = "black"),
        axis.title.x = element_text(size =44),
        axis.text.x = element_text(size = 16, colour = "black"))

kruskal.test(rumen_lumen_lachnospiraceae_family$Abundance, rumen_lumen_lachnospiraceae_family$treatment) #SIG

# now the individual genera

## RUMEN EPITH
rumen_epithelium_lachnospiraceae_genus <- subset_taxa(rumen_epithelium_genus, Family=="Lachnospiraceae")
rumen_epithelium_lachnospiraceae_genus_filt <- merge_low_abundance(rumen_epithelium_lachnospiraceae_genus, threshold = 0.1) # 15 taxa
write.csv(tax_table(rumen_epithelium_lachnospiraceae_genus_filt), "rumen_epith_lachno_taxa.csv")
write.csv(otu_table(rumen_epithelium_lachnospiraceae_genus_filt), "rumen_epith_lachno_otus.csv")
rumen_epithelium_lachnospiraceae_genus_melt <- psmelt(rumen_epithelium_lachnospiraceae_genus_filt)

rumen_epithelium_lachno_pallete <- c("#A8D8A8",	"#D8F0C0",	"#00F0A8",	"#007848",	"#00C0A8",	"#78D878",	"#007878",	"#30FF60",	"#183030",	"#186048",	"#C0C030",	"#309030",	"#00A890",	"#A8FFAB",	"grey88")
rumen_epithelium_lachnospiraceae_family <- subset_taxa(rumen_epithelium_family, Family=="Lachnospiraceae") %>%
  psmelt()

ggplot() + theme_bw() + labs(y= "Lachnospiraceae (RA%)", title = "EPITHELIUM") +
  coord_flip() + theme_bw() +
  geom_bar(rumen_epithelium_lachnospiraceae_genus_melt, colour = "black", linewidth = 0.75,
           mapping = aes(x= treatment, y= Abundance, fill=Genus), stat = "summary") +
  scale_fill_manual(values = rumen_epithelium_lachno_pallete) +
  scale_x_discrete(limits = c("tylosin","control"), labels =c("TYLOSIN","NO TYLOSIN")) +
  scale_y_continuous(expand = c(0.001,0,0.918,0)) +
  geom_errorbar(rumen_epithelium_lachnospiraceae_family, mapping= aes(x= treatment, y= Abundance),
                stat = "summary", width = 0.5, linewidth = 0.75) +
  theme(legend.position = "none",
        plot.title = element_text(size = 44),
        panel.border = element_rect(colour = "black", size = 1.25),
        panel.grid.major.y = element_blank(),
        axis.ticks = element_line(size = 0.75, colour= "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 34, colour = "black"),
        axis.title.x = element_text(size =44),
        axis.text.x = element_text(size = 16, colour = "black"))
kruskal.test(rumen_epithelium_lachnospiraceae_family$Abundance, rumen_epithelium_lachnospiraceae_family$treatment) #SIG

# now the individual genera
#ruminococcus_gauvreauii grp.
rumen_epithelium_ruminococcus_gauvreauii <- rumen_epithelium_lachno_genus_melt[which(rumen_epithelium_lachno_genus_melt$Genus=="[Ruminococcus]_gauvreauii_group"),]
kruskal.test(rumen_epithelium_ruminococcus_gauvreauii$Abundance, rumen_epithelium_ruminococcus_gauvreauii$treatment)

#Acetitomaculum
rumen_epithelium_Acetitomaculum<- rumen_epithelium_lachno_genus_melt[which(rumen_epithelium_lachno_genus_melt$Genus=="Acetitomaculum"),]
kruskal.test(rumen_epithelium_Acetitomaculum$Abundance, rumen_epithelium_Acetitomaculum$treatment)

#Eubacteirum hallii group
rumen_epithlium_Eubacterium_hallii_group <- rumen_epithelium_lachno_genus_melt[which(rumen_epithelium_lachno_genus_melt$Genus=="[Eubacterium]_hallii_group"),]
kruskal.test(rumen_epithlium_Eubacterium_hallii_group$Abundance, rumen_epithlium_Eubacterium_hallii_group$treatment) #NS


#Oribacterium
rumen_epithelium_Oribacterium <- rumen_epithelium_lachno_genus_melt[which(rumen_epithelium_lachno_genus_melt$Genus=="Oribacterium"),]
kruskal.test(rumen_epithelium_Oribacterium$Abundance, rumen_epithelium_Oribacterium$treatment)

#Lachnospiraceae NK3A20 group
rumen_epithelium_NK3A20_group <- rumen_epithelium_lachno_genus_melt[which(rumen_epithelium_lachno_genus_melt$Genus=="Lachnospiraceae_NK3A20_group"),]
kruskal.test(rumen_epithelium_NK3A20_group$Abundance, rumen_epithelium_NK3A20_group$treatment)

#Shuttleworthia
rumen_epithelium_Shuttleworthia <- rumen_epithelium_lachno_genus_melt[which(rumen_epithelium_lachno_genus_melt$Genus=="Shuttleworthia"),]
kruskal.test(rumen_epithelium_Shuttleworthia$Abundance, rumen_epithelium_Shuttleworthia$treatment)

#Syntrophococcus
rumen_epithelium_Syntrophococcus <- rumen_epithelium_lachno_genus_melt[which(rumen_epithelium_lachno_genus_melt$Genus=="Syntrophococcus"),]
kruskal.test(rumen_epithelium_Syntrophococcus$Abundance, rumen_epithelium_Syntrophococcus$treatment)

# unclassified_Lachnospiraceae
rumen_epithelium_unclassified_Lachnospiraceae <- rumen_epithelium_lachno_genus_melt[which(rumen_epithelium_lachno_genus_melt$Genus=="unclassified Lachnospiraceae"),]
kruskal.test(rumen_epithelium_unclassified_Lachnospiraceae$Abundance, rumen_epithelium_unclassified_Lachnospiraceae$treatment)

### SM INTESTINE
## SM LUMEN
sm_intestine_lumen_lachnospiraceae_genus <- subset_taxa(sm_intestine_lumen_genus, Family=="Lachnospiraceae")
sm_intestine_lumen_lachnospiraceae_genus_filt <- merge_low_abundance(sm_intestine_lumen_lachnospiraceae_genus, threshold = 0.1) # 13 taxa
write.csv(tax_table(sm_intestine_lumen_lachnospiraceae_genus_filt), "sm_intestine_lumen_lachno_taxa.csv")
write.csv(otu_table(sm_intestine_lumen_lachnospiraceae_genus_filt), "sm_intestine_lumen_lachno_otus.csv")
sm_intestine_lumen_lachnospiraceae_genus_melt <- psmelt(sm_intestine_lumen_lachnospiraceae_genus_filt)

sm_intestine_lumen_lachno_pallete <- c("#A8D8A8",	"#00F0A8",	"#007848",	"#78D878",	"#007878",	"#30FF60",	"#309030",	"#00A890",	"#A8FFAB",	"grey88")
sm_intestine_lumen_lachnospiraceae_family <- subset_taxa(sm_intestine_lumen_family, Family=="Lachnospiraceae") %>%
  psmelt()

ggplot() + theme_bw() + labs(y= "Lachnospiraceae (RA%)", title = "LUMEN") +
  coord_flip() + theme_bw() +
  geom_bar(sm_intestine_lumen_lachnospiraceae_genus_melt, colour = "black", linewidth = 0.75,
           mapping = aes(x= treatment, y= Abundance, fill=Genus), stat = "summary") +
  scale_fill_manual(values = sm_intestine_lumen_lachno_pallete) +
  scale_x_discrete(limits = c("tylosin","control"), labels =c("TYLOSIN","NO TYLOSIN")) +
  scale_y_continuous(expand = c(0.001,0,0.2355,0)) +
  geom_errorbar(sm_intestine_lumen_lachnospiraceae_family, mapping= aes(x= treatment, y= Abundance),
                stat = "summary", width = 0.5, linewidth = 0.75) +
  theme(#legend.position = "none",
    plot.title = element_text(size = 44),
    panel.border = element_rect(colour = "black", size = 1.25),
    panel.grid.major.y = element_blank(),
    axis.ticks = element_line(size = 0.75, colour= "black"),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 32, colour = "black"),
    axis.title.x = element_text(size =36),
    axis.text.x = element_text(size = 16, colour = "black"))

kruskal.test(sm_intestine_lumen_lachnospiraceae_family$Abundance, sm_intestine_lumen_lachnospiraceae_family$treatment) #SIG

# now the individual genera
#ruminococcus_gauvreauii grp.
sm_intestine_lumen_ruminococcus_gauvreauii <- sm_intestine_lumen_lachno_genus_melt[which(sm_intestine_lumen_lachno_genus_melt$Genus=="[Ruminococcus]_gauvreauii_group"),]
kruskal.test(sm_intestine_lumen_ruminococcus_gauvreauii$Abundance, sm_intestine_lumen_ruminococcus_gauvreauii$treatment) # sig.

#Acetitomaculum
sm_intestine_lumen_Acetitomaculum <- sm_intestine_lumen_lachno_genus_melt[which(sm_intestine_lumen_lachno_genus_melt$Genus=="Acetitomaculum"),]
kruskal.test(sm_intestine_lumen_Acetitomaculum$Abundance, sm_intestine_lumen_Acetitomaculum$treatment) #NS

#Eubacteirum hallii group
sm_intestine_lumen_Eubacterium_hallii_group <- sm_intestine_lumen_lachno_genus_melt[which(sm_intestine_lumen_lachno_genus_melt$Genus=="[Eubacterium]_hallii_group"),]
kruskal.test(sm_intestine_lumen_Eubacterium_hallii_group$Abundance, sm_intestine_lumen_Eubacterium_hallii_group$treatment) #NS

#Lachnospiraceae NK3A20 group
sm_intestine_lumen_NK3A20_group <- sm_intestine_lumen_lachno_genus_melt[which(sm_intestine_lumen_lachno_genus_melt$Genus=="Lachnospiraceae_NK3A20_group"),]
kruskal.test(sm_intestine_lumen_NK3A20_group$Abundance, sm_intestine_lumen_NK3A20_group$treatment) #NS

#Syntrophococcus
sm_intestine_lumen_Syntrophococcus <- sm_intestine_lumen_lachno_genus_melt[which(sm_intestine_lumen_lachno_genus_melt$Genus=="Syntrophococcus"),]
kruskal.test(sm_intestine_lumen_Syntrophococcus$Abundance, sm_intestine_lumen_Syntrophococcus$treatment) #NS

# unclassified_Lachnospiraceae
sm_intestine_lumen_unclassified_Lachnospiraceae <- sm_intestine_lumen_lachno_genus_melt[which(sm_intestine_lumen_lachno_genus_melt$Genus=="unclassified Lachnospiraceae"),]
kruskal.test(sm_intestine_lumen_unclassified_Lachnospiraceae$Abundance, sm_intestine_lumen_unclassified_Lachnospiraceae$treatment) #sig.

# oribacterium
sm_intestine_lumen_oribacterium <- sm_intestine_lumen_lachno_genus_melt[which(sm_intestine_lumen_lachno_genus_melt$Genus=="Oribacterium"),]
kruskal.test(sm_intestine_lumen_oribacterium$Abundance, sm_intestine_lumen_oribacterium$treatment) #NS



## SM EPITH
sm_intestine_epithelium_lachnospiraceae_genus <- subset_taxa(sm_intestine_epithelium_genus, Family=="Lachnospiraceae")
sm_intestine_epithelium_lachnospiraceae_genus_filt <- merge_low_abundance(sm_intestine_epithelium_lachnospiraceae_genus, threshold = 0.1) # 15 taxa
write.csv(tax_table(sm_intestine_epithelium_lachnospiraceae_genus_filt), "sm_intestine_epith_lachno_taxa.csv")
write.csv(otu_table(sm_intestine_epithelium_lachnospiraceae_genus_filt), "sm_intestine_epith_lachno_otus.csv")
sm_intestine_epithelium_lachnospiraceae_genus_melt <- psmelt(sm_intestine_epithelium_lachnospiraceae_genus_filt)

sm_intestine_epithelium_lachno_pallete <- c("#A8D8A8",	"#00F0A8",	"#007848",	"#90C0A8",	"#78D878",	"#007878",	"#30FF60",	"#309030",	"#00A890",	"#A8FFAB",	"grey88")
sm_intestine_epithelium_lachnospiraceae_family <- subset_taxa(sm_intestine_epithelium_family, Family=="Lachnospiraceae") %>%
  psmelt()

ggplot() + theme_bw() + labs(y= "Lachnospiraceae (RA%)", title = "EPITHELIUM") +
  coord_flip() + theme_bw() +
  geom_bar(sm_intestine_epithelium_lachnospiraceae_genus_melt, colour = "black", linewidth = 0.75,
           mapping = aes(x= treatment, y= Abundance, fill=Genus), stat = "summary") +
  scale_fill_manual(values = sm_intestine_epithelium_lachno_pallete) +
  scale_x_discrete(limits = c("tylosin","control"), labels =c("TYLOSIN","NO TYLOSIN")) +
  scale_y_continuous(expand = c(0.001,0,0.092,0)) +
  geom_errorbar(sm_intestine_epithelium_lachnospiraceae_family, mapping= aes(x= treatment, y= Abundance),
                stat = "summary", width = 0.5, linewidth = 0.75) +
  theme(legend.position = "none",
    plot.title = element_text(size = 44),
    panel.border = element_rect(colour = "black", size = 1.25),
    panel.grid.major.y = element_blank(),
    axis.ticks = element_line(size = 0.75, colour= "black"),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 32, colour = "black"),
    axis.title.x = element_text(size =36),
    axis.text.x = element_text(size = 16, colour = "black"))
kruskal.test(sm_intestine_epithelium_lachnospiraceae_family$Abundance, sm_intestine_epithelium_lachnospiraceae_family$treatment) #SIG

### LG INTESTINE
## LG LUMEN
lg_intestine_lumen_lachnospiraceae_genus <- subset_taxa(lg_intestine_lumen_genus, Family=="Lachnospiraceae")
lg_intestine_lumen_lachnospiraceae_genus_filt <- merge_low_abundance(lg_intestine_lumen_lachnospiraceae_genus, threshold = 0.1) # 13 taxa
write.csv(tax_table(lg_intestine_lumen_lachnospiraceae_genus_filt), "lg_intestine_lumen_lachno_taxa.csv")
write.csv(otu_table(lg_intestine_lumen_lachnospiraceae_genus_filt), "lg_intestine_lumen_lachno_otus.csv")
lg_intestine_lumen_lachnospiraceae_genus_melt <- psmelt(lg_intestine_lumen_lachnospiraceae_genus_filt)

lg_intestine_lumen_lachno_pallete <- c("#A8D8A8",	"#00F0A8",	"#909078",	"#007848",	"#609078",	"#609048",	"#90C0A8",	"#A8C090",	"#487830",	"#A8FF48",	"#78D878",	"#007878",	"#30FF60",	"#486060",	"#48A890",	"#487800",	"#309030",	"#00A890",	"#A8FFAB",	"grey88")
lg_intestine_lumen_lachnospiraceae_family <- subset_taxa(lg_intestine_lumen_family, Family=="Lachnospiraceae") %>%
  psmelt()

ggplot() + theme_bw() + labs(y= "Lachnospiraceae (RA%)", title = "LUMEN") +
  coord_flip() + theme_bw() +
  geom_bar(lg_intestine_lumen_lachnospiraceae_genus_melt, colour = "black", linewidth = 0.75,
           mapping = aes(x= treatment, y= Abundance, fill=Genus), stat = "summary") +
  scale_fill_manual(values = lg_intestine_lumen_lachno_pallete) +
  scale_x_discrete(limits = c("tylosin","control"), labels =c("TYLOSIN","NO TYLOSIN")) +
  scale_y_continuous(expand = c(0.001,0,1.215,0)) +
  geom_errorbar(lg_intestine_lumen_lachnospiraceae_family, mapping= aes(x= treatment, y= Abundance),
                stat = "summary", width = 0.5, linewidth = 0.75) +
  theme(legend.position = "none",
    plot.title = element_text(size = 44),
    panel.border = element_rect(colour = "black", size = 1.25),
    panel.grid.major.y = element_blank(),
    axis.ticks = element_line(size = 0.75, colour= "black"),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 32, colour = "black"),
    axis.title.x = element_text(size =36),
    axis.text.x = element_text(size = 16, colour = "black"))

kruskal.test(lg_intestine_lumen_lachnospiraceae_family$Abundance, lg_intestine_lumen_lachnospiraceae_family$treatment) #SIG


## LG EPITH
lg_intestine_epithelium_lachnospiraceae_genus <- subset_taxa(lg_intestine_epithelium_genus, Family=="Lachnospiraceae")
lg_intestine_epithelium_lachnospiraceae_genus_filt <- merge_low_abundance(lg_intestine_epithelium_lachnospiraceae_genus, threshold = 0.1) # 15 taxa
write.csv(tax_table(lg_intestine_epithelium_lachnospiraceae_genus_filt), "lg_intestine_epith_lachno_taxa.csv")
write.csv(otu_table(lg_intestine_epithelium_lachnospiraceae_genus_filt), "lg_intestine_epith_lachno_otus.csv")
lg_intestine_epithelium_lachnospiraceae_genus_melt <- psmelt(lg_intestine_epithelium_lachnospiraceae_genus_filt)

lg_intestine_epithelium_lachno_pallete <- c("#A8D8A8",	"#00F0A8",	"#909078",	"#007848",	"#609078",	"#609048",	"#90C0A8",	"#A8C090",	"#487830",	"#A8FF48",	"#78D878",	"#609060",	"#30FF60",	"#486060",	"#48A890",	"#487800",	"#309030",	"#00A890",	"#A8FFAB",	"grey88")
lg_intestine_epithelium_lachnospiraceae_family <- subset_taxa(lg_intestine_epithelium_family, Family=="Lachnospiraceae") %>%
  psmelt()

ggplot() + theme_bw() + labs(y= "Lachnospiraceae (RA%)", title = "EPITHELIUM") +
  coord_flip() + theme_bw() +
  geom_bar(lg_intestine_epithelium_lachnospiraceae_genus_melt, colour = "black", linewidth = 0.75,
           mapping = aes(x= treatment, y= Abundance, fill=Genus), stat = "summary") +
  scale_fill_manual(values = lg_intestine_epithelium_lachno_pallete) +
  scale_x_discrete(limits = c("tylosin","control"), labels =c("TYLOSIN","NO TYLOSIN")) +
  scale_y_continuous(expand = c(0.001,0,1.685,0)) +
  geom_errorbar(lg_intestine_epithelium_lachnospiraceae_family, mapping= aes(x= treatment, y= Abundance),
                stat = "summary", width = 0.5, linewidth = 0.75) +
  theme(legend.position = "none",
    plot.title = element_text(size = 44),
    panel.border = element_rect(colour = "black", size = 1.25),
    panel.grid.major.y = element_blank(),
    axis.ticks = element_line(size = 0.75, colour= "black"),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 32, colour = "black"),
    axis.title.x = element_text(size =36),
    axis.text.x = element_text(size = 16, colour = "black"))
kruskal.test(lg_intestine_epithelium_lachnospiraceae_family$Abundance, lg_intestine_epithelium_lachnospiraceae_family$treatment) #SIG




## genus plot themes
genus_plot_details = theme(legend.position = "none",
                         panel.border = element_rect(linewidth = 1, colour = "black"),
                         panel.grid.major.x = element_blank(),
                         axis.title.y = element_text(size = 36),
                         axis.text.y = element_text(size = 14, colour = "black"),
                         axis.ticks.y = element_line(linewidth = 0.75, colour= "black"),
                         axis.ticks.x = element_blank(),
                         axis.title.x = element_blank(),
                         axis.text.x = element_blank())

### Ruminococcus_gauvreauii_group
rumen_lumen_ruminococcus_gauvreauii <- rumen_lumen_lachno_genus_melt[which(rumen_lumen_lachno_genus_melt$Genus=="[Ruminococcus]_gauvreauii_group"),]
rumen_epithelium_ruminococcus_gauvreauii <- rumen_epithelium_lachno_genus_melt[which(rumen_epithelium_lachno_genus_melt$Genus=="[Ruminococcus]_gauvreauii_group"),]
sm_intestine_lumen_ruminococcus_gauvreauii <- sm_intestine_lumen_lachno_genus_melt[which(sm_intestine_lumen_lachno_genus_melt$Genus=="[Ruminococcus]_gauvreauii_group"),]
sm_intestine_epithelium_ruminococcus_gauvreauii <- sm_intestine_epithelium_lachno_genus_melt[which(sm_intestine_epithelium_lachno_genus_melt$Genus=="[Ruminococcus]_gauvreauii_group"),]
lg_intestine_lumen_ruminococcus_gauvreauii <- lg_intestine_lumen_lachno_genus_melt[which(lg_intestine_lumen_lachno_genus_melt$Genus=="[Ruminococcus]_gauvreauii_group"),]
lg_intestine_epithelium_ruminococcus_gauvreauii <- lg_intestine_epithelium_lachnospiraceae_genus_melt[which(lg_intestine_epithelium_lachnospiraceae_genus_melt$Genus=="[Ruminococcus]_gauvreauii_group"),]

ggplot(rumen_lumen_ruminococcus_gauvreauii, aes(x=treatment, y= Abundance, fill = treatment, colour = treatment)) +
  theme_bw() + labs(y= "RUMEN LUMEN (RA%)") +
  geom_bar(stat = "summary", alpha=0.5, linewidth = 1) + geom_errorbar(stat = "summary", linewidth = 1, width = .5) +
  scale_y_continuous(expand = c(0.001,0,0.8495,0)) + scale_fill_manual(values = tylosin_palette) + scale_colour_manual(values = tylosin_palette) +
  genus_plot_details

ggplot(rumen_epithelium_ruminococcus_gauvreauii, aes(x=treatment, y= Abundance, fill = treatment, colour = treatment)) +
  theme_bw() + labs(y= "RUMEN EPITH. (RA%)") +
  geom_bar(stat = "summary", alpha=0.5, linewidth = 1) + geom_errorbar(stat = "summary", linewidth = 1, width = .5) +
  scale_y_continuous(expand = c(0.001,0,3.095,0)) + scale_fill_manual(values = tylosin_palette) + scale_colour_manual(values = tylosin_palette) +
  genus_plot_details

ggplot(sm_intestine_lumen_ruminococcus_gauvreauii, aes(x=treatment, y= Abundance, fill = treatment, colour = treatment)) +
  theme_bw() + labs(y= "SM. INT. LUMEN (RA%)") +
  geom_bar(stat = "summary", alpha=0.5, linewidth = 1) + geom_errorbar(stat = "summary", linewidth = 1, width = .5) +
  scale_y_continuous(expand = c(0.001,0,0.235,0)) + scale_fill_manual(values = tylosin_palette) + scale_colour_manual(values = tylosin_palette) +
  genus_plot_details

ggplot(sm_intestine_epithelium_ruminococcus_gauvreauii, aes(x=treatment, y= Abundance, fill = treatment, colour = treatment)) +
  theme_bw() + labs(y= "SM. INT. EPITH. (RA%)") +
  geom_bar(stat = "summary", alpha=0.5, linewidth = 1) + geom_errorbar(stat = "summary", linewidth = 1, width = .5) +
  scale_y_continuous(expand = c(0.001,0,0.1,0)) + scale_fill_manual(values = tylosin_palette) + scale_colour_manual(values = tylosin_palette) +
  genus_plot_details

ggplot(lg_intestine_lumen_ruminococcus_gauvreauii, aes(x=treatment, y= Abundance, fill = treatment, colour = treatment)) +
  theme_bw() + labs(y= "LG. INT. LUMEN (RA%)") +
  geom_bar(stat = "summary", alpha=0.5, linewidth = 1) + geom_errorbar(stat = "summary", linewidth = 1, width = .5) +
  scale_y_continuous(expand = c(0.001,0,4.905,0)) + scale_fill_manual(values = tylosin_palette) + scale_colour_manual(values = tylosin_palette) +
  genus_plot_details

ggplot(lg_intestine_epithelium_ruminococcus_gauvreauii, aes(x=treatment, y= Abundance, fill = treatment, colour = treatment)) +
  theme_bw() + labs(y= "LG. INT. LUMEN (RA%)") +
  geom_bar(stat = "summary", alpha=0.5, linewidth = 1) + geom_errorbar(stat = "summary", linewidth = 1, width = .5) +
  scale_y_continuous(expand = c(0.001,0,5.77,0)) + scale_fill_manual(values = tylosin_palette) + scale_colour_manual(values = tylosin_palette) +
  genus_plot_details




### PREVOTELLACEAE ####
## RUMEN LUMEN
rumen_lumen_prevotellaceae_genus <- subset_taxa(rumen_lumen_genus, Family=="Prevotellaceae")
rumen_lumen_prevotellaceae_genus_filt <- merge_low_abundance(rumen_lumen_prevotellaceae_genus, threshold = 0.1) # 7 taxa
write.csv(tax_table(rumen_lumen_prevotellaceae_genus_filt), "rumen_lumen_prevotellaceae_taxa.csv")
write.csv(otu_table(rumen_lumen_prevotellaceae_genus_filt), "rumen_lumen_prevotellaceae_otus.csv")
rumen_lumen_prevotellaceae_genus_melt <- psmelt(rumen_lumen_prevotellaceae_genus_filt)

rumen_lumen_prevotellaceae_palette <- c("#006090",	"#7890D8",	"#486078",	"#003048",	"#78C0F0",	"#A8D8F0",	"grey88")
rumen_lumen_prevotellaceae_melt <-  psmelt(rumen_lumen_prevotellaceae_genus) # 37 genera
rumen_lumen_prevotellaceae_family <- subset_taxa(rumen_lumen_family, Family=="Prevotellaceae") %>%
  psmelt()

ggplot() + theme_bw() + labs(y= "Prevotellaceae (RA%)", title = "LUMEN") +
  coord_flip() + theme_bw() +
  geom_bar(rumen_lumen_prevotellaceae_genus_melt, colour = "black", linewidth = 0.75,
           mapping = aes(x= treatment, y= Abundance, fill=Genus), stat = "summary") +
  scale_fill_manual(values = rumen_lumen_prevotellaceae_palette) +
  scale_x_discrete(limits = c("tylosin","control"), labels =c("TYLOSIN","NO TYLOSIN")) +
  scale_y_continuous(expand = c(0.001,0,1.13,0)) +
  geom_errorbar(rumen_lumen_prevotellaceae_family, mapping= aes(x= treatment, y= Abundance),
                stat = "summary", width = 0.5, linewidth = 0.75) +
  theme(#legend.position = "none",
        plot.title = element_text(size = 44),
    panel.border = element_rect(colour = "black", size = 1.25),
    panel.grid.major.y = element_blank(),
    axis.ticks = element_line(size = 0.75, colour= "black"),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 34, colour = "black"),
    axis.title.x = element_text(size =44),
    axis.text.x = element_text(size = 16, colour = "black"))

kruskal.test(rumen_lumen_prevotellaceae_family$Abundance, rumen_lumen_prevotellaceae_family$treatment) #NS

## RUMEN EPITHELIUM
rumen_epithelium_prevotellaceae_genus <- subset_taxa(rumen_epithelium_genus, Family=="Prevotellaceae")
rumen_epithelium_prevotellaceae_genus_filt <- merge_low_abundance(rumen_epithelium_prevotellaceae_genus, threshold = 0.1) # 7 taxa
write.csv(tax_table(rumen_epithelium_prevotellaceae_genus_filt), "rumen_epithelium_prevotellaceae_taxa.csv")
write.csv(otu_table(rumen_epithelium_prevotellaceae_genus_filt), "rumen_epithelium_prevotellaceae_otus.csv")
rumen_epithelium_prevotellaceae_genus_melt <- psmelt(rumen_epithelium_prevotellaceae_genus_filt)

rumen_epithelium_prevotellaceae_palette <- c("#006090",	"#7890D8",	"#486078",	"#003048",	"#78C0F0",	"#A8D8F0",	"#D8F0FF",	"grey88")
rumen_epithelium_prevotellaceae_melt <-  psmelt(rumen_epithelium_prevotellaceae_genus) # 37 genera
rumen_epithelium_prevotellaceae_family <- subset_taxa(rumen_epithelium_family, Family=="Prevotellaceae") %>%
  psmelt()

ggplot() + theme_bw() + labs(y= "Prevotellaceae (RA%)", title = "EPITHELIUM") +
  coord_flip() + theme_bw() +
  geom_bar(rumen_epithelium_prevotellaceae_genus_melt , colour = "black", linewidth = 0.75,
           mapping = aes(x= treatment, y= Abundance, fill=Genus), stat = "summary") +
  scale_fill_manual(values = rumen_epithelium_prevotellaceae_palette) +
  scale_x_discrete(limits = c("tylosin","control"), labels =c("TYLOSIN","NO TYLOSIN")) +
  scale_y_continuous(expand = c(0.001,0,0.199,0)) +
  geom_errorbar(rumen_epithelium_prevotellaceae_family, mapping= aes(x= treatment, y= Abundance),
                stat = "summary", width = 0.5, linewidth = 0.75) +
  theme(legend.position = "none",
        plot.title = element_text(size = 44),
        panel.border = element_rect(colour = "black", size = 1.25),
        panel.grid.major.y = element_blank(),
        axis.ticks = element_line(size = 0.75, colour= "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 34, colour = "black"),
        axis.title.x = element_text(size =44),
        axis.text.x = element_text(size = 16, colour = "black"))

kruskal.test(rumen_epithelium_prevotellaceae_family$Abundance, rumen_epithelium_prevotellaceae_family$treatment) #NS

## SM INT LUMEN
sm_intestine_lumen_prevotellaceae_genus <- subset_taxa(sm_intestine_lumen_genus, Family=="Prevotellaceae")
sm_intestine_lumen_prevotellaceae_genus_filt <- merge_low_abundance(sm_intestine_lumen_prevotellaceae_genus, threshold = 0.1) # 7 taxa
write.csv(tax_table(sm_intestine_lumen_prevotellaceae_genus_filt), "sm_intestine_lumen_prevotellaceae_taxa.csv")
write.csv(otu_table(sm_intestine_lumen_prevotellaceae_genus_filt), "sm_intestine_lumen_prevotellaceae_otus.csv")
sm_intestine_lumen_prevotellaceae_genus_melt <- psmelt(sm_intestine_lumen_prevotellaceae_genus_filt)

sm_intestine_lumen_prevotellaceae_palette <- c("#CBD794",	"dodgerblue2",	"#77E37D",	"#D7ACD2",	"cyan2",	"cadetblue",	"steelblue2",	"blue4",	"cornflowerblue",	"lightskyblue1",	"#7FE1CD")
sm_intestine_lumen_prevotellaceae_melt <-  psmelt(sm_intestine_lumen_prevotellaceae_genus) # 37 genera
sm_intestine_lumen_prevotellaceae_family <- subset_taxa(sm_intestine_lumen_family, Family=="Prevotellaceae") %>%
  psmelt()

ggplot() + theme_bw() + labs(y= "Prevotellaceae (RA%)") +
  coord_flip() + theme_bw() +
  geom_bar(sm_intestine_lumen_prevotellaceae_melt, colour = "black", linewidth = 0.75,
           mapping = aes(x= treatment, y= Abundance, fill=Genus), stat = "summary") +
  scale_fill_manual(values = sm_intestine_lumen_prevotellaceae_palette) +
  scale_x_discrete(limits = c("tylosin","control"), labels =c("TYLOSIN","NO TYLOSIN")) +
  scale_y_continuous(expand = c(0.001,0,1.13,0)) +
  geom_errorbar(sm_intestine_lumen_prevotellaceae_family, mapping= aes(x= treatment, y= Abundance),
                stat = "summary", width = 0.5, linewidth = 0.75) +
  theme(#legend.position = "none",
    panel.border = element_rect(colour = "black", size = 1.25),
    panel.grid.major.y = element_blank(),
    axis.ticks = element_line(size = 0.75, colour= "black"),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 32, colour = "black"),
    axis.title.x = element_text(size =36),
    axis.text.x = element_text(size = 16, colour = "black"))

kruskal.test(sm_intestine_lumen_prevotellaceae_family$Abundance, sm_intestine_lumen_prevotellaceae_family$treatment) #NS

## SM INT EPITHELIUM
sm_intestine_epithelium_prevotellaceae_genus <- subset_taxa(sm_intestine_epithelium_genus, Family=="Prevotellaceae")
sm_intestine_epithelium_prevotellaceae_genus_filt <- merge_low_abundance(sm_intestine_epithelium_prevotellaceae_genus, threshold = 0.1) # 7 taxa
write.csv(tax_table(sm_intestine_epithelium_prevotellaceae_genus_filt), "sm_intestine_epithelium_prevotellaceae_taxa.csv")
write.csv(otu_table(sm_intestine_epithelium_prevotellaceae_genus_filt), "sm_intestine_epithelium_prevotellaceae_otus.csv")
sm_intestine_epithelium_prevotellaceae_genus_melt <- psmelt(sm_intestine_epithelium_prevotellaceae_genus_filt)

sm_intestine_epithelium_prevotellaceae_palette <- c("#CBD794",	"#80B4CD",	"dodgerblue2",	"#77E37D",	"#D7ACD2",	"#D65B69",	"cadetblue",	"steelblue3",	"blue4",	"cornflowerblue",	"lightskyblue1",	"#7FE1CD")
sm_intestine_epithelium_prevotellaceae_melt <-  psmelt(sm_intestine_epithelium_prevotellaceae_genus) # 37 genera
sm_intestine_epithelium_prevotellaceae_family <- subset_taxa(sm_intestine_epithelium_family, Family=="Prevotellaceae") %>%
  psmelt()

ggplot() + theme_bw() + labs(y= "Prevotellaceae (RA%)") +
  coord_flip() + theme_bw() +
  geom_bar(sm_intestine_epithelium_prevotellaceae_melt, colour = "black", linewidth = 0.75,
           mapping = aes(x= treatment, y= Abundance, fill=Genus), stat = "summary") +
  scale_fill_manual(values = sm_intestine_epithelium_prevotellaceae_palette) +
  scale_x_discrete(limits = c("tylosin","control"), labels =c("TYLOSIN","NO TYLOSIN")) +
  scale_y_continuous(expand = c(0.001,0,0.199,0)) +
  geom_errorbar(sm_intestine_epithelium_prevotellaceae_family, mapping= aes(x= treatment, y= Abundance),
                stat = "summary", width = 0.5, linewidth = 0.75) +
  theme(#legend.position = "none",
    panel.border = element_rect(colour = "black", size = 1.25),
    panel.grid.major.y = element_blank(),
    axis.ticks = element_line(size = 0.75, colour= "black"),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 32, colour = "black"),
    axis.title.x = element_text(size =36),
    axis.text.x = element_text(size = 16, colour = "black"))

kruskal.test(sm_intestine_epithelium_prevotellaceae_family$Abundance, sm_intestine_epithelium_prevotellaceae_family$treatment) #NS

## LG INT LUMEN
lg_intestine_lumen_prevotellaceae_genus <- subset_taxa(lg_intestine_lumen_genus, Family=="Prevotellaceae")
lg_intestine_lumen_prevotellaceae_genus_filt <- merge_low_abundance(lg_intestine_lumen_prevotellaceae_genus, threshold = 0.1) # 7 taxa
write.csv(tax_table(lg_intestine_lumen_prevotellaceae_genus_filt), "lg_intestine_lumen_prevotellaceae_taxa.csv")
write.csv(otu_table(lg_intestine_lumen_prevotellaceae_genus_filt), "lg_intestine_lumen_prevotellaceae_otus.csv")
lg_intestine_lumen_prevotellaceae_genus_melt <- psmelt(lg_intestine_lumen_prevotellaceae_genus_filt)

lg_intestine_lumen_prevotellaceae_palette <- c("#CBD794",	"dodgerblue2",	"#77E37D",	"#D7ACD2",	"cyan2",	"cadetblue",	"steelblue2",	"blue4",	"cornflowerblue",	"lightskyblue1",	"#7FE1CD")
lg_intestine_lumen_prevotellaceae_melt <-  psmelt(lg_intestine_lumen_prevotellaceae_genus) # 37 genera
lg_intestine_lumen_prevotellaceae_family <- subset_taxa(lg_intestine_lumen_family, Family=="Prevotellaceae") %>%
  psmelt()

ggplot() + theme_bw() + labs(y= "Prevotellaceae (RA%)") +
  coord_flip() + theme_bw() +
  geom_bar(lg_intestine_lumen_prevotellaceae_melt, colour = "black", linewidth = 0.75,
           mapping = aes(x= treatment, y= Abundance, fill=Genus), stat = "summary") +
  scale_fill_manual(values = lg_intestine_lumen_prevotellaceae_palette) +
  scale_x_discrete(limits = c("tylosin","control"), labels =c("TYLOSIN","NO TYLOSIN")) +
  scale_y_continuous(expand = c(0.001,0,1.13,0)) +
  geom_errorbar(lg_intestine_lumen_prevotellaceae_family, mapping= aes(x= treatment, y= Abundance),
                stat = "summary", width = 0.5, linewidth = 0.75) +
  theme(#legend.position = "none",
    panel.border = element_rect(colour = "black", size = 1.25),
    panel.grid.major.y = element_blank(),
    axis.ticks = element_line(size = 0.75, colour= "black"),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 32, colour = "black"),
    axis.title.x = element_text(size =36),
    axis.text.x = element_text(size = 16, colour = "black"))

kruskal.test(lg_intestine_lumen_prevotellaceae_family$Abundance, lg_intestine_lumen_prevotellaceae_family$treatment) #NS

## LG INTEPITHELIUM
lg_intestine_epithelium_prevotellaceae_genus <- subset_taxa(lg_intestine_epithelium_genus, Family=="Prevotellaceae")
lg_intestine_epithelium_prevotellaceae_genus_filt <- merge_low_abundance(lg_intestine_epithelium_prevotellaceae_genus, threshold = 0.1) # 7 taxa
write.csv(tax_table(lg_intestine_epithelium_prevotellaceae_genus_filt), "lg_intestine_epithelium_prevotellaceae_taxa.csv")
write.csv(otu_table(lg_intestine_epithelium_prevotellaceae_genus_filt), "lg_intestine_epithelium_prevotellaceae_otus.csv")
lg_intestine_epithelium_prevotellaceae_genus_melt <- psmelt(lg_intestine_epithelium_prevotellaceae_genus_filt)

lg_intestine_epithelium_prevotellaceae_palette <- c("#CBD794",	"#80B4CD",	"dodgerblue2",	"#77E37D",	"#D7ACD2",	"#D65B69",	"cadetblue",	"steelblue3",	"blue4",	"cornflowerblue",	"lightskyblue1",	"#7FE1CD")
lg_intestine_epithelium_prevotellaceae_melt <-  psmelt(lg_intestine_epithelium_prevotellaceae_genus) # 37 genera
lg_intestine_epithelium_prevotellaceae_family <- subset_taxa(lg_intestine_epithelium_family, Family=="Prevotellaceae") %>%
  psmelt()

ggplot() + theme_bw() + labs(y= "Prevotellaceae (RA%)") +
  coord_flip() + theme_bw() +
  geom_bar(lg_intestine_epithelium_prevotellaceae_melt, colour = "black", linewidth = 0.75,
           mapping = aes(x= treatment, y= Abundance, fill=Genus), stat = "summary") +
  scale_fill_manual(values = lg_intestine_epithelium_prevotellaceae_palette) +
  scale_x_discrete(limits = c("tylosin","control"), labels =c("TYLOSIN","NO TYLOSIN")) +
  scale_y_continuous(expand = c(0.001,0,0.199,0)) +
  geom_errorbar(lg_intestine_epithelium_prevotellaceae_family, mapping= aes(x= treatment, y= Abundance),
                stat = "summary", width = 0.5, linewidth = 0.75) +
  theme(#legend.position = "none",
    panel.border = element_rect(colour = "black", size = 1.25),
    panel.grid.major.y = element_blank(),
    axis.ticks = element_line(size = 0.75, colour= "black"),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 32, colour = "black"),
    axis.title.x = element_text(size =36),
    axis.text.x = element_text(size = 16, colour = "black"))

kruskal.test(lg_intestine_epithelium_prevotellaceae_family$Abundance, lg_intestine_epithelium_prevotellaceae_family$treatment) #NS


### ATOPOBIACEAE ####
## RUMEN LUMEN
rumen_lumen_Atopobiaceae_genus <- subset_taxa(rumen_lumen_genus, Family=="Atopobiaceae")
rumen_lumen_Atopobiaceae_genus_filt <- merge_low_abundance(rumen_lumen_Atopobiaceae_genus, threshold = 0.1) # 7 taxa
write.csv(tax_table(rumen_lumen_Atopobiaceae_genus_filt), "rumen_lumen_Atopobiaceae_taxa.csv")
write.csv(otu_table(rumen_lumen_Atopobiaceae_genus_filt), "rumen_lumen_Atopobiaceae_otus.csv")
rumen_lumen_Atopobiaceae_genus_melt <- psmelt(rumen_lumen_Atopobiaceae_genus_filt)

rumen_lumen_atopobiaceae_palette <- c("#E7C57D","ivory3","cornsilk1")
rumen_lumen_atopobiaceae_melt <-  psmelt(rumen_lumen_atopobiaceae_genus) # 37 genera
rumen_lumen_atopobiaceae_family <- subset_taxa(rumen_lumen_family, Family=="Atopobiaceae") %>%
  psmelt()

ggplot() + theme_bw() + labs(y= "Atopobiaceae (RA%)", title = "LUMEN") +
  coord_flip() + theme_bw() +
  geom_bar(rumen_lumen_Atopobiaceae_genus_melt, colour = "black", linewidth = 0.75,
           mapping = aes(x= treatment, y= Abundance, fill=Genus), stat = "summary") +
  scale_fill_manual(values = rumen_lumen_atopobiaceae_palette) +
  scale_x_discrete(limits = c("tylosin","control"), labels =c("TYLOSIN","NO TYLOSIN")) +
  scale_y_continuous(expand = c(0.001,0,2.89,0)) +
  geom_errorbar(rumen_lumen_atopobiaceae_family, mapping= aes(x= treatment, y= Abundance),
                stat = "summary", width = 0.5, linewidth = 0.75) +
  theme(legend.position = "none",
        plot.title = element_text(size = 44),
        panel.border = element_rect(colour = "black", size = 1.25),
        panel.grid.major.y = element_blank(),
        axis.ticks = element_line(size = 0.75, colour= "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 34, colour = "black"),
        axis.title.x = element_text(size = 44),
        axis.text.x = element_text(size = 16, colour = "black"))

kruskal.test(rumen_lumen_atopobiaceae_family$Abundance, rumen_lumen_atopobiaceae_family$treatment) #NS

## RUMEN EPITHELIUM
rumen_epithelium_Atopobiaceae_genus <- subset_taxa(rumen_epithelium_genus, Family=="Atopobiaceae")
rumen_epithelium_Atopobiaceae_genus_filt <- merge_low_abundance(rumen_epithelium_Atopobiaceae_genus, threshold = 0.1) # 7 taxa
write.csv(tax_table(rumen_epithelium_Atopobiaceae_genus_filt), "rumen_epithelium_Atopobiaceae_taxa.csv")
write.csv(otu_table(rumen_epithelium_Atopobiaceae_genus_filt), "rumen_epithelium_Atopobiaceae_otus.csv")
rumen_epithelium_Atopobiaceae_genus_melt <- psmelt(rumen_epithelium_Atopobiaceae_genus_filt)

rumen_epithelium_atopobiaceae_palette <- c("#E7C57D","cornsilk1")
rumen_epithelium_atopobiaceae_melt <-  psmelt(rumen_epithelium_atopobiaceae_genus) # 37 genera
rumen_epithelium_atopobiaceae_family <- subset_taxa(rumen_epithelium_family, Family=="Atopobiaceae") %>%
  psmelt()

ggplot() + theme_bw() + labs(y= "Atopobiaceae (RA%)", title = "EPITHELIUM") +
  coord_flip() + theme_bw() +
  geom_bar(rumen_epithelium_Atopobiaceae_genus_melt, colour = "black", linewidth = 0.75,
           mapping = aes(x= treatment, y= Abundance, fill=Genus), stat = "summary") +
  scale_fill_manual(values = rumen_epithelium_atopobiaceae_palette) +
  scale_x_discrete(limits = c("tylosin","control"), labels =c("TYLOSIN","NO TYLOSIN")) +
  scale_y_continuous(expand = c(0.001,0,11.819,0)) +
  geom_errorbar(rumen_epithelium_atopobiaceae_family, mapping= aes(x= treatment, y= Abundance),
                stat = "summary", width = 0.5, linewidth = 0.75) +
  theme(legend.position = "none",
        plot.title = element_text(size = 44),
        panel.border = element_rect(colour = "black", size = 1.25),
        panel.grid.major.y = element_blank(),
        axis.ticks = element_line(size = 0.75, colour= "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 34, colour = "black"),
        axis.title.x = element_text(size =44),
        axis.text.x = element_text(size = 16, colour = "black"))

kruskal.test(rumen_epithelium_atopobiaceae_family$Abundance, rumen_epithelium_atopobiaceae_family$treatment) #NS

### SM INTESTINE LUMEN
sm_intestine_lumen_Atopobiaceae_genus <- subset_taxa(sm_intestine_lumen_genus, Family=="Atopobiaceae")
sm_intestine_lumen_Atopobiaceae_genus_filt <- merge_low_abundance(sm_intestine_lumen_Atopobiaceae_genus, threshold = 0.1) # 7 taxa
write.csv(tax_table(sm_intestine_lumen_Atopobiaceae_genus_filt), "sm_intestine_lumen_Atopobiaceae_taxa.csv")
write.csv(otu_table(sm_intestine_lumen_Atopobiaceae_genus_filt), "sm_intestine_lumen_Atopobiaceae_otus.csv")
sm_intestine_lumen_Atopobiaceae_genus_melt <- psmelt(sm_intestine_lumen_Atopobiaceae_genus_filt)

sm_intestine_lumen_atopobiaceae_palette <- c("#E7C57D","cornsilk1","grey88")
sm_intestine_lumen_atopobiaceae_family <- subset_taxa(sm_intestine_lumen_family, Family=="Atopobiaceae") %>%
  psmelt()

ggplot() + theme_bw() + labs(y= "Atopobiaceae (RA%)") +
  coord_flip() + theme_bw() +
  geom_bar(sm_intestine_lumen_Atopobiaceae_genus_melt, colour = "black", linewidth = 0.75,
           mapping = aes(x= treatment, y= Abundance, fill=Genus), stat = "summary") +
  scale_fill_manual(values = sm_intestine_lumen_atopobiaceae_palette) +
  scale_x_discrete(limits = c("tylosin","control"), labels =c("TYLOSIN","NO TYLOSIN")) +
  scale_y_continuous(expand = c(0.001,0,0.17,0)) +
  geom_errorbar(sm_intestine_lumen_atopobiaceae_family, mapping= aes(x= treatment, y= Abundance),
                stat = "summary", width = 0.5, linewidth = 0.75) +
  theme(#legend.position = "none",
        panel.border = element_rect(colour = "black", size = 1.25),
        panel.grid.major.y = element_blank(),
        axis.ticks = element_line(size = 0.75, colour= "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 32, colour = "black"),
        axis.title.x = element_text(size =36),
        axis.text.x = element_text(size = 16, colour = "black"))

kruskal.test(sm_intestine_lumen_atopobiaceae_family$Abundance, sm_intestine_lumen_atopobiaceae_family$treatment) #NS

## sm_intestine EPITHELIUM
sm_intestine_epithelium_Atopobiaceae_genus <- subset_taxa(sm_intestine_epithelium_genus, Family=="Atopobiaceae")
sm_intestine_epithelium_Atopobiaceae_genus_filt <- merge_low_abundance(sm_intestine_epithelium_Atopobiaceae_genus, threshold = 0.1) # 7 taxa
write.csv(tax_table(sm_intestine_epithelium_Atopobiaceae_genus_filt), "sm_intestine_epithelium_Atopobiaceae_taxa.csv")
write.csv(otu_table(sm_intestine_epithelium_Atopobiaceae_genus_filt), "sm_intestine_epithelium_Atopobiaceae_otus.csv")
sm_intestine_epithelium_Atopobiaceae_genus_melt <- psmelt(sm_intestine_epithelium_Atopobiaceae_genus_filt)

sm_intestine_epithelium_atopobiaceae_palette <- c("#CBD794",	"#80B4CD",	"dodgerblue2",	"#77E37D",	"#D7ACD2",	"#D65B69",	"cadetblue",	"steelblue3",	"blue4",	"cornflowerblue",	"lightskyblue1",	"#7FE1CD")
sm_intestine_epithelium_atopobiaceae_family <- subset_taxa(sm_intestine_epithelium_family, Family=="Atopobiaceae") %>%
  psmelt()

ggplot() + theme_bw() + labs(y= "Atopobiaceae (RA%)") +
  coord_flip() + theme_bw() +
  geom_bar(sm_intestine_epithelium_Atopobiaceae_genus_melt, colour = "black", linewidth = 0.75,
           mapping = aes(x= treatment, y= Abundance, fill=Genus), stat = "summary") +
  scale_fill_manual(values = sm_intestine_epithelium_atopobiaceae_palette) +
  scale_x_discrete(limits = c("tylosin","control"), labels =c("TYLOSIN","NO TYLOSIN")) +
  scale_y_continuous(expand = c(0.001,0,1,0)) +
  geom_errorbar(sm_intestine_epithelium_atopobiaceae_family, mapping= aes(x= treatment, y= Abundance),
                stat = "summary", width = 0.5, linewidth = 0.75) +
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", size = 1.25),
        panel.grid.major.y = element_blank(),
        axis.ticks = element_line(size = 0.75, colour= "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 32, colour = "black"),
        axis.title.x = element_text(size =36),
        axis.text.x = element_text(size = 16, colour = "black"))

kruskal.test(sm_intestine_epithelium_atopobiaceae_family$Abundance, sm_intestine_epithelium_atopobiaceae_family$treatment) #NS

### LG INT LUMEN
lg_intestine_lumen_Atopobiaceae_genus <- subset_taxa(lg_intestine_lumen_genus, Family=="Atopobiaceae")
lg_intestine_lumen_Atopobiaceae_genus_filt <- merge_low_abundance(lg_intestine_lumen_Atopobiaceae_genus, threshold = 0.1) # 7 taxa
write.csv(tax_table(lg_intestine_lumen_Atopobiaceae_genus_filt), "lg_intestine_lumen_Atopobiaceae_taxa.csv")
write.csv(otu_table(lg_intestine_lumen_Atopobiaceae_genus_filt), "lg_intestine_lumen_Atopobiaceae_otus.csv")
lg_intestine_lumen_Atopobiaceae_genus_melt <- psmelt(lg_intestine_lumen_Atopobiaceae_genus_filt)

lg_intestine_lumen_atopobiaceae_palette <- c("#B5AAD8",	"bisque",	"#E7C57D",	"#89E379",	"palegoldenrod")
lg_intestine_lumen_atopobiaceae_melt <-  psmelt(lg_intestine_lumen_atopobiaceae_genus) # 37 genera
lg_intestine_lumen_atopobiaceae_family <- subset_taxa(lg_intestine_lumen_family, Family=="Atopobiaceae") %>%
  psmelt()

ggplot() + theme_bw() + labs(y= "Atopobiaceae (RA%)") +
  coord_flip() + theme_bw() +
  geom_bar(lg_intestine_lumen_Atopobiaceae_genus_melt, colour = "black", linewidth = 0.75,
           mapping = aes(x= treatment, y= Abundance, fill=Genus), stat = "summary") +
  scale_fill_manual(values = lg_intestine_lumen_atopobiaceae_palette) +
  scale_x_discrete(limits = c("tylosin","control"), labels =c("TYLOSIN","NO TYLOSIN")) +
  scale_y_continuous(expand = c(0.001,0,0.17,0)) +
  geom_errorbar(lg_intestine_lumen_atopobiaceae_family, mapping= aes(x= treatment, y= Abundance),
                stat = "summary", width = 0.5, linewidth = 0.75) +
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", size = 1.25),
        panel.grid.major.y = element_blank(),
        axis.ticks = element_line(size = 0.75, colour= "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 32, colour = "black"),
        axis.title.x = element_text(size =36),
        axis.text.x = element_text(size = 16, colour = "black"))

kruskal.test(lg_intestine_lumen_atopobiaceae_family$Abundance, lg_intestine_lumen_atopobiaceae_family$treatment) #NS

## lg_intestine EPITHELIUM
lg_intestine_epithelium_Atopobiaceae_genus <- subset_taxa(lg_intestine_epithelium_genus, Family=="Atopobiaceae")
lg_intestine_epithelium_Atopobiaceae_genus_filt <- merge_low_abundance(lg_intestine_epithelium_Atopobiaceae_genus, threshold = 0.1) # 7 taxa
write.csv(tax_table(lg_intestine_epithelium_Atopobiaceae_genus_filt), "lg_intestine_epithelium_Atopobiaceae_taxa.csv")
write.csv(otu_table(lg_intestine_epithelium_Atopobiaceae_genus_filt), "lg_intestine_epithelium_Atopobiaceae_otus.csv")
lg_intestine_epithelium_Atopobiaceae_genus_melt <- psmelt(lg_intestine_epithelium_Atopobiaceae_genus_filt)

lg_intestine_epithelium_atopobiaceae_palette <- c("#CBD794",	"#80B4CD",	"dodgerblue2",	"#77E37D",	"#D7ACD2",	"#D65B69",	"cadetblue",	"steelblue3",	"blue4",	"cornflowerblue",	"lightskyblue1",	"#7FE1CD")
lg_intestine_epithelium_atopobiaceae_melt <-  psmelt(lg_intestine_epithelium_atopobiaceae_genus) # 37 genera
lg_intestine_epithelium_atopobiaceae_family <- subset_taxa(lg_intestine_epithelium_family, Family=="Atopobiaceae") %>%
  psmelt()

ggplot() + theme_bw() + labs(y= "Atopobiaceae (RA%)") +
  coord_flip() + theme_bw() +
  geom_bar(lg_intestine_epithelium_Atopobiaceae_genus_melt, colour = "black", linewidth = 0.75,
           mapping = aes(x= treatment, y= Abundance, fill=Genus), stat = "summary") +
  scale_fill_manual(values = lg_intestine_epithelium_atopobiaceae_palette) +
  scale_x_discrete(limits = c("tylosin","control"), labels =c("TYLOSIN","NO TYLOSIN")) +
  scale_y_continuous(expand = c(0.001,0,11.85,0)) +
  geom_errorbar(lg_intestine_epithelium_atopobiaceae_family, mapping= aes(x= treatment, y= Abundance),
                stat = "summary", width = 0.5, linewidth = 0.75) +
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", size = 1.25),
        panel.grid.major.y = element_blank(),
        axis.ticks = element_line(size = 0.75, colour= "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 32, colour = "black"),
        axis.title.x = element_text(size =36),
        axis.text.x = element_text(size = 16, colour = "black"))

kruskal.test(lg_intestine_epithelium_atopobiaceae_family$Abundance, lg_intestine_epithelium_atopobiaceae_family$treatment) #NS


### OSCILLOSPIRACEAE ####
## RUMEN LUMEN
rumen_lumen_Oscillospiraceae_genus <- subset_taxa(rumen_lumen_genus, Family=="Oscillospiraceae")
rumen_lumen_Oscillospiraceae_genus_filt <- merge_low_abundance(rumen_lumen_Oscillospiraceae_genus, threshold = 0.1) # 7 taxa
write.csv(tax_table(rumen_lumen_Oscillospiraceae_genus_filt), "rumen_lumen_Oscillospiraceae_taxa.csv")
write.csv(otu_table(rumen_lumen_Oscillospiraceae_genus_filt), "rumen_lumen_Oscillospiraceae_otus.csv")
rumen_lumen_Oscillospiraceae_genus_melt <- psmelt(rumen_lumen_Oscillospiraceae_genus_filt)

rumen_lumen_oscillospiraceae_palette <- c("#FF7518",	"#9c4b08",	"#FF6030",	"grey88")
rumen_lumen_oscillospiraceae_family <- subset_taxa(rumen_lumen_family, Family=="Oscillospiraceae") %>%
  psmelt()

ggplot() + theme_bw() + labs(y= "Oscillospiraceae (RA%)", title = "LUMEN") +
  coord_flip() + theme_bw() +
  geom_bar(rumen_lumen_Oscillospiraceae_genus_melt, colour = "black", linewidth = 0.75,
           mapping = aes(x= treatment, y= Abundance, fill=Genus), stat = "summary") +
  scale_fill_manual(values = rumen_lumen_oscillospiraceae_palette) +
  scale_x_discrete(limits = c("tylosin","control"), labels =c("TYLOSIN","NO TYLOSIN")) +
  scale_y_continuous(expand = c(0.001,0,3.849,0)) +
  geom_errorbar(rumen_lumen_oscillospiraceae_family, mapping= aes(x= treatment, y= Abundance),
                stat = "summary", width = 0.5, linewidth = 0.75) +
  theme(#legend.position = "none",
        plot.title = element_text(size = 44),
        panel.border = element_rect(colour = "black", size = 1.25),
        panel.grid.major.y = element_blank(),
        axis.ticks = element_line(size = 0.75, colour= "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 34, colour = "black"),
        axis.title.x = element_text(size = 44),
        axis.text.x = element_text(size = 16, colour = "black"))

kruskal.test(rumen_lumen_oscillospiraceae_family$Abundance, rumen_lumen_oscillospiraceae_family$treatment) #NS

## RUMEN EPITHELIUM
rumen_epithelium_Oscillospiraceae_genus <- subset_taxa(rumen_epithelium_genus, Family=="Oscillospiraceae")
rumen_epithelium_Oscillospiraceae_genus_filt <- merge_low_abundance(rumen_epithelium_Oscillospiraceae_genus, threshold = 0.1) # 7 taxa
write.csv(tax_table(rumen_epithelium_Oscillospiraceae_genus_filt), "rumen_epithelium_Oscillospiraceae_taxa.csv")
write.csv(otu_table(rumen_epithelium_Oscillospiraceae_genus_filt), "rumen_epithelium_Oscillospiraceae_otus.csv")
rumen_epithelium_Oscillospiraceae_genus_melt <- psmelt(rumen_epithelium_Oscillospiraceae_genus_filt)

rumen_epithelium_oscillospiraceae_palette <- c("#FF7518",	"#9c4b08",	"#FF6030",	"grey88")
rumen_epithelium_oscillospiraceae_family <- subset_taxa(rumen_epithelium_family, Family=="Oscillospiraceae") %>%
  psmelt()

ggplot() + theme_bw() + labs(y= "Oscillospiraceae (RA%)", title = "EPITHELIUM") +
  coord_flip() + theme_bw() +
  geom_bar(rumen_epithelium_Oscillospiraceae_genus_melt, colour = "black", linewidth = 0.75,
           mapping = aes(x= treatment, y= Abundance, fill=Genus), stat = "summary") +
  scale_fill_manual(values = rumen_epithelium_oscillospiraceae_palette) +
  scale_x_discrete(limits = c("tylosin","control"), labels =c("TYLOSIN","NO TYLOSIN")) +
  scale_y_continuous(expand = c(0.001,0,9.67,0)) +
  geom_errorbar(rumen_epithelium_oscillospiraceae_family, mapping= aes(x= treatment, y= Abundance),
                stat = "summary", width = 0.5, linewidth = 0.75) +
  theme(legend.position = "none",
        plot.title = element_text(size = 44),
        panel.border = element_rect(colour = "black", size = 1.25),
        panel.grid.major.y = element_blank(),
        axis.ticks = element_line(size = 0.75, colour= "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 34, colour = "black"),
        axis.title.x = element_text(size =44),
        axis.text.x = element_text(size = 16, colour = "black"))

kruskal.test(rumen_epithelium_oscillospiraceae_family$Abundance, rumen_epithelium_oscillospiraceae_family$treatment) #NS

### LG INT LUMEN
lg_intestine_lumen_Oscillospiraceae_genus <- subset_taxa(lg_intestine_lumen_genus, Family=="Oscillospiraceae")
lg_intestine_lumen_Oscillospiraceae_genus_filt <- merge_low_abundance(lg_intestine_lumen_Oscillospiraceae_genus, threshold = 0.1) # 7 taxa
write.csv(tax_table(lg_intestine_lumen_Oscillospiraceae_genus_filt), "lg_intestine_lumen_Oscillospiraceae_taxa.csv")
write.csv(otu_table(lg_intestine_lumen_Oscillospiraceae_genus_filt), "lg_intestine_lumen_Oscillospiraceae_otus.csv")
lg_intestine_lumen_Oscillospiraceae_genus_melt <- psmelt(lg_intestine_lumen_Oscillospiraceae_genus_filt)

lg_intestine_lumen_oscillospiraceae_palette <- c("#FF7518",	"#E89362",	"#9c4b08",	"orangered3",	"#FBCEB1",	"#FF6030",	"grey88")
lg_intestine_lumen_oscillospiraceae_melt <-  psmelt(lg_intestine_lumen_oscillospiraceae_genus) # 37 genera
lg_intestine_lumen_oscillospiraceae_family <- subset_taxa(lg_intestine_lumen_family, Family=="Oscillospiraceae") %>%
  psmelt()

ggplot() + theme_bw() + labs(y= "Oscillospiraceae (RA%)") +
  coord_flip() + theme_bw() +
  geom_bar(lg_intestine_lumen_Oscillospiraceae_genus_melt, colour = "black", linewidth = 0.75,
           mapping = aes(x= treatment, y= Abundance, fill=Genus), stat = "summary") +
  scale_fill_manual(values = lg_intestine_lumen_oscillospiraceae_palette) +
  scale_x_discrete(limits = c("tylosin","control"), labels =c("TYLOSIN","NO TYLOSIN")) +
  scale_y_continuous(expand = c(0.001,0,2.5,0)) +
  geom_errorbar(lg_intestine_lumen_oscillospiraceae_family, mapping= aes(x= treatment, y= Abundance),
                stat = "summary", width = 0.5, linewidth = 0.75) +
  theme(#legend.position = "none",
        panel.border = element_rect(colour = "black", size = 1.25),
        panel.grid.major.y = element_blank(),
        axis.ticks = element_line(size = 0.75, colour= "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 32, colour = "black"),
        axis.title.x = element_text(size =36),
        axis.text.x = element_text(size = 16, colour = "black"))

kruskal.test(lg_intestine_lumen_oscillospiraceae_family$Abundance, lg_intestine_lumen_oscillospiraceae_family$treatment) #NS

## lg_intestine EPITHELIUM
lg_intestine_epithelium_Oscillospiraceae_genus <- subset_taxa(lg_intestine_epithelium_genus, Family=="Oscillospiraceae")
lg_intestine_epithelium_Oscillospiraceae_genus_filt <- merge_low_abundance(lg_intestine_epithelium_Oscillospiraceae_genus, threshold = 0.1) # 7 taxa
write.csv(tax_table(lg_intestine_epithelium_Oscillospiraceae_genus_filt), "lg_intestine_epithelium_Oscillospiraceae_taxa.csv")
write.csv(otu_table(lg_intestine_epithelium_Oscillospiraceae_genus_filt), "lg_intestine_epithelium_Oscillospiraceae_otus.csv")
lg_intestine_epithelium_Oscillospiraceae_genus_melt <- psmelt(lg_intestine_epithelium_Oscillospiraceae_genus_filt)

lg_intestine_epithelium_oscillospiraceae_palette <- c("#FF7518",	"#E89362",	"#9c4b08",	"orangered3",	"#FBCEB1",	"#FF6030",	"grey88")
lg_intestine_epithelium_oscillospiraceae_family <- subset_taxa(lg_intestine_epithelium_family, Family=="Oscillospiraceae") %>%
  psmelt()

ggplot() + theme_bw() + labs(y= "Oscillospiraceae (RA%)") +
  coord_flip() + theme_bw() +
  geom_bar(lg_intestine_epithelium_Oscillospiraceae_genus_melt, colour = "black", linewidth = 0.75,
           mapping = aes(x= treatment, y= Abundance, fill=Genus), stat = "summary") +
  scale_fill_manual(values = lg_intestine_epithelium_oscillospiraceae_palette) +
  scale_x_discrete(limits = c("tylosin","control"), labels =c("TYLOSIN","NO TYLOSIN")) +
  scale_y_continuous(expand = c(0.001,0,2.5,0)) +
  geom_errorbar(lg_intestine_epithelium_oscillospiraceae_family, mapping= aes(x= treatment, y= Abundance),
                stat = "summary", width = 0.5, linewidth = 0.75) +
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", size = 1.25),
        panel.grid.major.y = element_blank(),
        axis.ticks = element_line(size = 0.75, colour= "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 32, colour = "black"),
        axis.title.x = element_text(size =36),
        axis.text.x = element_text(size = 16, colour = "black"))

kruskal.test(lg_intestine_epithelium_oscillospiraceae_family$Abundance, lg_intestine_epithelium_oscillospiraceae_family$treatment) #NS

### Succinivibrionaceae ####
## RUMEN LUMEN
rumen_lumen_Succinivibrionaceae_genus <- subset_taxa(rumen_lumen_genus, Family=="Succinivibrionaceae")
rumen_lumen_Succinivibrionaceae_genus_filt <- merge_low_abundance(rumen_lumen_Succinivibrionaceae_genus, threshold = 0.1) # 7 taxa
write.csv(tax_table(rumen_lumen_Succinivibrionaceae_genus_filt), "rumen_lumen_Succinivibrionaceae_taxa.csv")
write.csv(otu_table(rumen_lumen_Succinivibrionaceae_genus_filt), "rumen_lumen_Succinivibrionaceae_otus.csv")
rumen_lumen_Succinivibrionaceae_genus_melt <- psmelt(rumen_lumen_Succinivibrionaceae_genus_filt)

rumen_lumen_Succinivibrionaceae_palette <- c("hotpink","grey88")
rumen_lumen_Succinivibrionaceae_melt <-  psmelt(rumen_lumen_Succinivibrionaceae_genus) # 37 genera
rumen_lumen_Succinivibrionaceae_family <- subset_taxa(rumen_lumen_family, Family=="Succinivibrionaceae") %>%
  psmelt()

ggplot() + theme_bw() + labs(y= "Succinivibrionaceae (RA%)", title = "LUMEN") +
  coord_flip() + theme_bw() +
  geom_bar(rumen_lumen_Succinivibrionaceae_genus_melt, colour = "black", linewidth = 0.75,
           mapping = aes(x= treatment, y= Abundance, fill=Genus), stat = "summary") +
  scale_fill_manual(values = rumen_lumen_Succinivibrionaceae_palette) +
  scale_x_discrete(limits = c("tylosin","control"), labels =c("TYLOSIN","NO TYLOSIN")) +
  scale_y_continuous(expand = c(0.001,0,24.8,0)) +
  geom_errorbar(rumen_lumen_Succinivibrionaceae_family, mapping= aes(x= treatment, y= Abundance),
                stat = "summary", width = 0.5, linewidth = 0.75) +
  theme(#legend.position = "none",
        plot.title = element_text(size = 44),
        panel.border = element_rect(colour = "black", size = 1.25),
        panel.grid.major.y = element_blank(),
        axis.ticks = element_line(size = 0.75, colour= "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 34, colour = "black"),
        axis.title.x = element_text(size = 44),
        axis.text.x = element_text(size = 16, colour = "black"))

kruskal.test(rumen_lumen_Succinivibrionaceae_family$Abundance, rumen_lumen_Succinivibrionaceae_family$treatment) #NS

## RUMEN EPITHELIUM
rumen_epithelium_Succinivibrionaceae_genus <- subset_taxa(rumen_epithelium_genus, Family=="Succinivibrionaceae")
rumen_epithelium_Succinivibrionaceae_genus_filt <- merge_low_abundance(rumen_epithelium_Succinivibrionaceae_genus, threshold = 0.1) # 7 taxa
write.csv(tax_table(rumen_epithelium_Succinivibrionaceae_genus_filt), "rumen_epithelium_Succinivibrionaceae_taxa.csv")
write.csv(otu_table(rumen_epithelium_Succinivibrionaceae_genus_filt), "rumen_epithelium_Succinivibrionaceae_otus.csv")
rumen_epithelium_Succinivibrionaceae_genus_melt <- psmelt(rumen_epithelium_Succinivibrionaceae_genus_filt)

rumen_epithelium_Succinivibrionaceae_palette <- c("deeppink4","hotpink","rosybrown1","grey88")
rumen_epithelium_Succinivibrionaceae_family <- subset_taxa(rumen_epithelium_family, Family=="Succinivibrionaceae") %>%
  psmelt()

ggplot() + theme_bw() + labs(y= "Succinivibrionaceae (RA%)", title = "EPITHELIUM") +
  coord_flip() + theme_bw() +
  geom_bar(rumen_epithelium_Succinivibrionaceae_genus_melt, colour = "black", linewidth = 0.75,
           mapping = aes(x= treatment, y= Abundance, fill=Genus), stat = "summary") +
  scale_fill_manual(values = rumen_epithelium_Succinivibrionaceae_palette) +
  scale_x_discrete(limits = c("tylosin","control"), labels =c("TYLOSIN","NO TYLOSIN")) +
  scale_y_continuous(expand = c(0.001,0,1.542,0)) +
  geom_errorbar(rumen_epithelium_Succinivibrionaceae_family, mapping= aes(x= treatment, y= Abundance),
                stat = "summary", width = 0.5, linewidth = 0.75) +
  theme(legend.position = "none",
        plot.title = element_text(size = 44),
        panel.border = element_rect(colour = "black", size = 1.25),
        panel.grid.major.y = element_blank(),
        axis.ticks = element_line(size = 0.75, colour= "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 34, colour = "black"),
        axis.title.x = element_text(size =44),
        axis.text.x = element_text(size = 16, colour = "black"))

kruskal.test(rumen_epithelium_Succinivibrionaceae_family$Abundance, rumen_epithelium_Succinivibrionaceae_family$treatment) #NS

### SM INTESTINE LUMEN
sm_intestine_lumen_Succinivibrionaceae_genus <- subset_taxa(sm_intestine_lumen_genus, Family=="Succinivibrionaceae")
sm_intestine_lumen_Succinivibrionaceae_genus_filt <- merge_low_abundance(sm_intestine_lumen_Succinivibrionaceae_genus, threshold = 0.1) # 7 taxa
write.csv(tax_table(sm_intestine_lumen_Succinivibrionaceae_genus_filt), "sm_intestine_lumen_Succinivibrionaceae_taxa.csv")
write.csv(otu_table(sm_intestine_lumen_Succinivibrionaceae_genus_filt), "sm_intestine_lumen_Succinivibrionaceae_otus.csv")
sm_intestine_lumen_Succinivibrionaceae_genus_melt <- psmelt(sm_intestine_lumen_Succinivibrionaceae_genus_filt)

sm_intestine_lumen_Succinivibrionaceae_palette <- c("#6030A8","grey88")
sm_intestine_lumen_Succinivibrionaceae_family <- subset_taxa(sm_intestine_lumen_family, Family=="Succinivibrionaceae") %>%
  psmelt()

ggplot() + theme_bw() + labs(y= "Succinivibrionaceae (RA%)", title = "LUMEN") +
  coord_flip() + theme_bw() +
  geom_bar(sm_intestine_lumen_Succinivibrionaceae_genus_melt, colour = "black", linewidth = 0.75,
           mapping = aes(x= treatment, y= Abundance, fill=Genus), stat = "summary") +
  scale_fill_manual(values = sm_intestine_lumen_Succinivibrionaceae_palette) +
  scale_x_discrete(limits = c("tylosin","control"), labels =c("TYLOSIN","NO TYLOSIN")) +
  scale_y_continuous(expand = c(0.001,0,93.6,0)) +
  geom_errorbar(sm_intestine_lumen_Succinivibrionaceae_family, mapping= aes(x= treatment, y= Abundance),
                stat = "summary", width = 0.5, linewidth = 0.75) +
  theme(legend.position = "none",
        plot.title = element_text(size = 44),
        panel.border = element_rect(colour = "black", size = 1.25),
        panel.grid.major.y = element_blank(),
        axis.ticks = element_line(size = 0.75, colour= "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 34, colour = "black"),
        axis.title.x = element_text(size =44),
        axis.text.x = element_text(size = 16, colour = "black"))

kruskal.test(sm_intestine_lumen_Succinivibrionaceae_family$Abundance, sm_intestine_lumen_Succinivibrionaceae_family$treatment) #NS

## sm_intestine EPITHELIUM
sm_intestine_epithelium_Succinivibrionaceae_genus <- subset_taxa(sm_intestine_epithelium_genus, Family=="Succinivibrionaceae")
sm_intestine_epithelium_Succinivibrionaceae_genus_filt <- merge_low_abundance(sm_intestine_epithelium_Succinivibrionaceae_genus, threshold = 0.1) # 7 taxa
write.csv(tax_table(sm_intestine_epithelium_Succinivibrionaceae_genus_filt), "sm_intestine_epithelium_Succinivibrionaceae_taxa.csv")
write.csv(otu_table(sm_intestine_epithelium_Succinivibrionaceae_genus_filt), "sm_intestine_epithelium_Succinivibrionaceae_otus.csv")
sm_intestine_epithelium_Succinivibrionaceae_genus_melt <- psmelt(sm_intestine_epithelium_Succinivibrionaceae_genus_filt)

sm_intestine_epithelium_Succinivibrionaceae_palette <- c("#6030A8","grey88")
sm_intestine_epithelium_Succinivibrionaceae_family <- subset_taxa(sm_intestine_epithelium_family, Family=="Succinivibrionaceae") %>%
  psmelt()

ggplot() + theme_bw() + labs(y= "Succinivibrionaceae (RA%)", title = "EPITHELIUM") +
  coord_flip() + theme_bw() +
  geom_bar(sm_intestine_epithelium_Succinivibrionaceae_genus_melt, colour = "black", linewidth = 0.75,
           mapping = aes(x= treatment, y= Abundance, fill=Genus), stat = "summary") +
  scale_fill_manual(values = sm_intestine_epithelium_Succinivibrionaceae_palette) +
  scale_x_discrete(limits = c("tylosin","control"), labels =c("TYLOSIN","NO TYLOSIN")) +
  scale_y_continuous(expand = c(0.001,0,112.5,0)) +
  geom_errorbar(sm_intestine_epithelium_Succinivibrionaceae_family, mapping= aes(x= treatment, y= Abundance),
                stat = "summary", width = 0.5, linewidth = 0.75) +
  theme(legend.position = "none",
        plot.title = element_text(size = 44),
        panel.border = element_rect(colour = "black", size = 1.25),
        panel.grid.major.y = element_blank(),
        axis.ticks = element_line(size = 0.75, colour= "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 34, colour = "black"),
        axis.title.x = element_text(size =44),
        axis.text.x = element_text(size = 16, colour = "black"))

kruskal.test(sm_intestine_epithelium_Succinivibrionaceae_family$Abundance, sm_intestine_epithelium_Succinivibrionaceae_family$treatment) #NS

### LG INT LUMEN
lg_intestine_lumen_Succinivibrionaceae_genus <- subset_taxa(lg_intestine_lumen_genus, Family=="Succinivibrionaceae")
lg_intestine_lumen_Succinivibrionaceae_genus_filt <- merge_low_abundance(lg_intestine_lumen_Succinivibrionaceae_genus, threshold = 0.1) # 7 taxa
write.csv(tax_table(lg_intestine_lumen_Succinivibrionaceae_genus_filt), "lg_intestine_lumen_Succinivibrionaceae_taxa.csv")
write.csv(otu_table(lg_intestine_lumen_Succinivibrionaceae_genus_filt), "lg_intestine_lumen_Succinivibrionaceae_otus.csv")
lg_intestine_lumen_Succinivibrionaceae_genus_melt <- psmelt(lg_intestine_lumen_Succinivibrionaceae_genus_filt)

lg_intestine_lumen_Succinivibrionaceae_palette <- c("#6030A8","grey88")
lg_intestine_lumen_Succinivibrionaceae_family <- subset_taxa(lg_intestine_lumen_family, Family=="Succinivibrionaceae") %>%
  psmelt()

ggplot() + theme_bw() + labs(y= "Succinivibrionaceae (RA%)", title = "LUMEN") +
  coord_flip() + theme_bw() +
  geom_bar(lg_intestine_lumen_Succinivibrionaceae_genus_melt, colour = "black", linewidth = 0.75,
           mapping = aes(x= treatment, y= Abundance, fill=Genus), stat = "summary") +
  scale_fill_manual(values = lg_intestine_lumen_Succinivibrionaceae_palette) +
  scale_x_discrete(limits = c("tylosin","control"), labels =c("TYLOSIN","NO TYLOSIN")) +
  scale_y_continuous(expand = c(0.001,0,23,0)) +
  geom_errorbar(lg_intestine_lumen_Succinivibrionaceae_family, mapping= aes(x= treatment, y= Abundance),
                stat = "summary", width = 0.5, linewidth = 0.75) +
  theme(legend.position = "none",
        plot.title = element_text(size = 44),
        panel.border = element_rect(colour = "black", size = 1.25),
        panel.grid.major.y = element_blank(),
        axis.ticks = element_line(size = 0.75, colour= "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 34, colour = "black"),
        axis.title.x = element_text(size =44),
        axis.text.x = element_text(size = 16, colour = "black"))

kruskal.test(lg_intestine_lumen_Succinivibrionaceae_family$Abundance, lg_intestine_lumen_Succinivibrionaceae_family$treatment) #NS

## lg_intestine EPITHELIUM
lg_intestine_epithelium_Succinivibrionaceae_genus <- subset_taxa(lg_intestine_epithelium_genus, Family=="Succinivibrionaceae")
lg_intestine_epithelium_Succinivibrionaceae_genus_filt <- merge_low_abundance(lg_intestine_epithelium_Succinivibrionaceae_genus, threshold = 0.1) # 7 taxa
write.csv(tax_table(lg_intestine_epithelium_Succinivibrionaceae_genus_filt), "lg_intestine_epithelium_Succinivibrionaceae_taxa.csv")
write.csv(otu_table(lg_intestine_epithelium_Succinivibrionaceae_genus_filt), "lg_intestine_epithelium_Succinivibrionaceae_otus.csv")
lg_intestine_epithelium_Succinivibrionaceae_genus_melt <- psmelt(lg_intestine_epithelium_Succinivibrionaceae_genus_filt)

lg_intestine_epithelium_Succinivibrionaceae_palette <- c("#4B0082","#6030A8","grey88")
lg_intestine_epithelium_Succinivibrionaceae_melt <-  psmelt(lg_intestine_epithelium_Succinivibrionaceae_genus) # 37 genera
lg_intestine_epithelium_Succinivibrionaceae_family <- subset_taxa(lg_intestine_epithelium_family, Family=="Succinivibrionaceae") %>%
  psmelt()

ggplot() + theme_bw() + labs(y= "Succinivibrionaceae (RA%)", title = "EPITHELIUM") +
  coord_flip() + theme_bw() +
  geom_bar(lg_intestine_epithelium_Succinivibrionaceae_genus_melt, colour = "black", linewidth = 0.75,
           mapping = aes(x= treatment, y= Abundance, fill=Genus), stat = "summary") +
  scale_fill_manual(values = lg_intestine_epithelium_Succinivibrionaceae_palette) +
  scale_x_discrete(limits = c("tylosin","control"), labels =c("TYLOSIN","NO TYLOSIN")) +
  scale_y_continuous(expand = c(0.001,0,45.3,0)) +
  geom_errorbar(lg_intestine_epithelium_Succinivibrionaceae_family, mapping= aes(x= treatment, y= Abundance),
                stat = "summary", width = 0.5, linewidth = 0.75) +
  theme(legend.position = "none",
        plot.title = element_text(size = 44),
        panel.border = element_rect(colour = "black", size = 1.25),
        panel.grid.major.y = element_blank(),
        axis.ticks = element_line(size = 0.75, colour= "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 34, colour = "black"),
        axis.title.x = element_text(size =44),
        axis.text.x = element_text(size = 16, colour = "black"))

kruskal.test(lg_intestine_epithelium_Succinivibrionaceae_family$Abundance, lg_intestine_epithelium_Succinivibrionaceae_family$treatment) #NS


### Selemonadaceae ####
## RUMEN LUMEN
rumen_lumen_Selenomonadaceae_genus <- subset_taxa(rumen_lumen_genus, Family=="Selenomonadaceae")
rumen_lumen_Selenomonadaceae_genus_filt <- merge_low_abundance(rumen_lumen_Selenomonadaceae_genus, threshold = 0.1) # 7 taxa
write.csv(tax_table(rumen_lumen_Selenomonadaceae_genus_filt), "rumen_lumen_Selenomonadaceae_taxa.csv")
write.csv(otu_table(rumen_lumen_Selenomonadaceae_genus_filt), "rumen_lumen_Selenomonadaceae_otus.csv")
rumen_lumen_Selenomonadaceae_genus_melt <- psmelt(rumen_lumen_Selenomonadaceae_genus_filt)

rumen_lumen_Selenomonadaceae_palette <- c("chocolate4","grey88")
rumen_lumen_Selenomonadaceae_family <- subset_taxa(rumen_lumen_family, Family=="Selenomonadaceae") %>%
  psmelt()

ggplot() + theme_bw() + labs(y= "Selenomonadaceae (RA%)", title = "LUMEN") +
  coord_flip() + theme_bw() +
  geom_bar(rumen_lumen_Selenomonadaceae_genus_melt, colour = "black", linewidth = 0.75,
           mapping = aes(x= treatment, y= Abundance, fill=Genus), stat = "summary") +
  scale_fill_manual(values = rumen_lumen_Selenomonadaceae_palette) +
  scale_x_discrete(limits = c("tylosin","control"), labels =c("TYLOSIN","NO TYLOSIN")) +
  scale_y_continuous(expand = c(0.001,0,22.96,0)) +
  geom_errorbar(rumen_lumen_Selenomonadaceae_family, mapping= aes(x= treatment, y= Abundance),
                stat = "summary", width = 0.5, linewidth = 0.75) +
  theme(legend.position = "none",
        plot.title = element_text(size = 44),
        panel.border = element_rect(colour = "black", size = 1.25),
        panel.grid.major.y = element_blank(),
        axis.ticks = element_line(size = 0.75, colour= "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 34, colour = "black"),
        axis.title.x = element_text(size = 44),
        axis.text.x = element_text(size = 16, colour = "black"))

kruskal.test(rumen_lumen_Selenomonadaceae_family$Abundance, rumen_lumen_Selenomonadaceae_family$treatment) #NS

## RUMEN EPITHELIUM
rumen_epithelium_Selenomonadaceae_genus <- subset_taxa(rumen_epithelium_genus, Family=="Selenomonadaceae")
rumen_epithelium_Selenomonadaceae_genus_filt <- merge_low_abundance(rumen_epithelium_Selenomonadaceae_genus, threshold = 0.1) # 7 taxa
write.csv(tax_table(rumen_epithelium_Selenomonadaceae_genus_filt), "rumen_epithelium_Selenomonadaceae_taxa.csv")
write.csv(otu_table(rumen_epithelium_Selenomonadaceae_genus_filt), "rumen_epithelium_Selenomonadaceae_otus.csv")
rumen_epithelium_Selenomonadaceae_genus_melt <- psmelt(rumen_epithelium_Selenomonadaceae_genus_filt)

rumen_epithelium_Selenomonadaceae_palette <- c("chocolate4","grey88")
rumen_epithelium_Selenomonadaceae_family <- subset_taxa(rumen_epithelium_family, Family=="Selenomonadaceae") %>%
  psmelt()

ggplot() + theme_bw() + labs(y= "Selenomonadaceae (RA%)", title = "EPITHELIUM") +
  coord_flip() + theme_bw() +
  geom_bar(rumen_epithelium_Selenomonadaceae_genus_melt, colour = "black", linewidth = 0.75,
           mapping = aes(x= treatment, y= Abundance, fill=Genus), stat = "summary") +
  scale_fill_manual(values = rumen_epithelium_Selenomonadaceae_palette) +
  scale_x_discrete(limits = c("tylosin","control"), labels =c("TYLOSIN","NO TYLOSIN")) +
  scale_y_continuous(expand = c(0.001,0,6.457,0)) +
  geom_errorbar(rumen_epithelium_Selenomonadaceae_family, mapping= aes(x= treatment, y= Abundance),
                stat = "summary", width = 0.5, linewidth = 0.75) +
  theme(legend.position = "none",
        plot.title = element_text(size = 44),
        panel.border = element_rect(colour = "black", size = 1.25),
        panel.grid.major.y = element_blank(),
        axis.ticks = element_line(size = 0.75, colour= "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 34, colour = "black"),
        axis.title.x = element_text(size =44),
        axis.text.x = element_text(size = 16, colour = "black"))

kruskal.test(rumen_epithelium_Selenomonadaceae_family$Abundance, rumen_epithelium_Selenomonadaceae_family$treatment) #NS

### SM INTESTINE LUMEN
sm_intestine_lumen_Selenomonadaceae_genus <- subset_taxa(sm_intestine_lumen_genus, Family=="Selenomonadaceae")
sm_intestine_lumen_Selenomonadaceae_genus_filt <- merge_low_abundance(sm_intestine_lumen_Selenomonadaceae_genus, threshold = 0.1) # 7 taxa
write.csv(tax_table(sm_intestine_lumen_Selenomonadaceae_genus_filt), "sm_intestine_lumen_Selenomonadaceae_taxa.csv")
write.csv(otu_table(sm_intestine_lumen_Selenomonadaceae_genus_filt), "sm_intestine_lumen_Selenomonadaceae_otus.csv")
sm_intestine_lumen_Selenomonadaceae_genus_melt <- psmelt(sm_intestine_lumen_Selenomonadaceae_genus_filt)

sm_intestine_lumen_Selenomonadaceae_palette <- c("#6030A8","grey88")
sm_intestine_lumen_Selenomonadaceae_family <- subset_taxa(sm_intestine_lumen_family, Family=="Selenomonadaceae") %>%
  psmelt()

ggplot() + theme_bw() + labs(y= "Selenomonadaceae (RA%)", title = "LUMEN") +
  coord_flip() + theme_bw() +
  geom_bar(sm_intestine_lumen_Selenomonadaceae_genus_melt, colour = "black", linewidth = 0.75,
           mapping = aes(x= treatment, y= Abundance, fill=Genus), stat = "summary") +
  scale_fill_manual(values = sm_intestine_lumen_Selenomonadaceae_palette) +
  scale_x_discrete(limits = c("tylosin","control"), labels =c("TYLOSIN","NO TYLOSIN")) +
  scale_y_continuous(expand = c(0.001,0,93.6,0)) +
  geom_errorbar(sm_intestine_lumen_Selenomonadaceae_family, mapping= aes(x= treatment, y= Abundance),
                stat = "summary", width = 0.5, linewidth = 0.75) +
  theme(legend.position = "none",
        plot.title = element_text(size = 44),
        panel.border = element_rect(colour = "black", size = 1.25),
        panel.grid.major.y = element_blank(),
        axis.ticks = element_line(size = 0.75, colour= "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 34, colour = "black"),
        axis.title.x = element_text(size =44),
        axis.text.x = element_text(size = 16, colour = "black"))

kruskal.test(sm_intestine_lumen_Selenomonadaceae_family$Abundance, sm_intestine_lumen_Selenomonadaceae_family$treatment) #NS

## sm_intestine EPITHELIUM
sm_intestine_epithelium_Selenomonadaceae_genus <- subset_taxa(sm_intestine_epithelium_genus, Family=="Selenomonadaceae")
sm_intestine_epithelium_Selenomonadaceae_genus_filt <- merge_low_abundance(sm_intestine_epithelium_Selenomonadaceae_genus, threshold = 0.1) # 7 taxa
write.csv(tax_table(sm_intestine_epithelium_Selenomonadaceae_genus_filt), "sm_intestine_epithelium_Selenomonadaceae_taxa.csv")
write.csv(otu_table(sm_intestine_epithelium_Selenomonadaceae_genus_filt), "sm_intestine_epithelium_Selenomonadaceae_otus.csv")
sm_intestine_epithelium_Selenomonadaceae_genus_melt <- psmelt(sm_intestine_epithelium_Selenomonadaceae_genus_filt)

sm_intestine_epithelium_Selenomonadaceae_palette <- c("#6030A8","grey88")
sm_intestine_epithelium_Selenomonadaceae_family <- subset_taxa(sm_intestine_epithelium_family, Family=="Selenomonadaceae") %>%
  psmelt()

ggplot() + theme_bw() + labs(y= "Selenomonadaceae (RA%)", title = "EPITHELIUM") +
  coord_flip() + theme_bw() +
  geom_bar(sm_intestine_epithelium_Selenomonadaceae_genus_melt, colour = "black", linewidth = 0.75,
           mapping = aes(x= treatment, y= Abundance, fill=Genus), stat = "summary") +
  scale_fill_manual(values = sm_intestine_epithelium_Selenomonadaceae_palette) +
  scale_x_discrete(limits = c("tylosin","control"), labels =c("TYLOSIN","NO TYLOSIN")) +
  scale_y_continuous(expand = c(0.001,0,112.5,0)) +
  geom_errorbar(sm_intestine_epithelium_Selenomonadaceae_family, mapping= aes(x= treatment, y= Abundance),
                stat = "summary", width = 0.5, linewidth = 0.75) +
  theme(legend.position = "none",
        plot.title = element_text(size = 44),
        panel.border = element_rect(colour = "black", size = 1.25),
        panel.grid.major.y = element_blank(),
        axis.ticks = element_line(size = 0.75, colour= "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 34, colour = "black"),
        axis.title.x = element_text(size =44),
        axis.text.x = element_text(size = 16, colour = "black"))

kruskal.test(sm_intestine_epithelium_Selenomonadaceae_family$Abundance, sm_intestine_epithelium_Selenomonadaceae_family$treatment) #NS

### LG INT LUMEN
lg_intestine_lumen_Selenomonadaceae_genus <- subset_taxa(lg_intestine_lumen_genus, Family=="Selenomonadaceae")
lg_intestine_lumen_Selenomonadaceae_genus_filt <- merge_low_abundance(lg_intestine_lumen_Selenomonadaceae_genus, threshold = 0.1) # 7 taxa
write.csv(tax_table(lg_intestine_lumen_Selenomonadaceae_genus_filt), "lg_intestine_lumen_Selenomonadaceae_taxa.csv")
write.csv(otu_table(lg_intestine_lumen_Selenomonadaceae_genus_filt), "lg_intestine_lumen_Selenomonadaceae_otus.csv")
lg_intestine_lumen_Selenomonadaceae_genus_melt <- psmelt(lg_intestine_lumen_Selenomonadaceae_genus_filt)

lg_intestine_lumen_Selenomonadaceae_palette <- c("#6030A8","grey88")
lg_intestine_lumen_Selenomonadaceae_family <- subset_taxa(lg_intestine_lumen_family, Family=="Selenomonadaceae") %>%
  psmelt()

ggplot() + theme_bw() + labs(y= "Selenomonadaceae (RA%)", title = "LUMEN") +
  coord_flip() + theme_bw() +
  geom_bar(lg_intestine_lumen_Selenomonadaceae_genus_melt, colour = "black", linewidth = 0.75,
           mapping = aes(x= treatment, y= Abundance, fill=Genus), stat = "summary") +
  scale_fill_manual(values = lg_intestine_lumen_Selenomonadaceae_palette) +
  scale_x_discrete(limits = c("tylosin","control"), labels =c("TYLOSIN","NO TYLOSIN")) +
  scale_y_continuous(expand = c(0.001,0,23,0)) +
  geom_errorbar(lg_intestine_lumen_Selenomonadaceae_family, mapping= aes(x= treatment, y= Abundance),
                stat = "summary", width = 0.5, linewidth = 0.75) +
  theme(legend.position = "none",
        plot.title = element_text(size = 44),
        panel.border = element_rect(colour = "black", size = 1.25),
        panel.grid.major.y = element_blank(),
        axis.ticks = element_line(size = 0.75, colour= "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 34, colour = "black"),
        axis.title.x = element_text(size =44),
        axis.text.x = element_text(size = 16, colour = "black"))

kruskal.test(lg_intestine_lumen_Selenomonadaceae_family$Abundance, lg_intestine_lumen_Selenomonadaceae_family$treatment) #NS

## lg_intestine EPITHELIUM
lg_intestine_epithelium_Selenomonadaceae_genus <- subset_taxa(lg_intestine_epithelium_genus, Family=="Selenomonadaceae")
lg_intestine_epithelium_Selenomonadaceae_genus_filt <- merge_low_abundance(lg_intestine_epithelium_Selenomonadaceae_genus, threshold = 0.1) # 7 taxa
write.csv(tax_table(lg_intestine_epithelium_Selenomonadaceae_genus_filt), "lg_intestine_epithelium_Selenomonadaceae_taxa.csv")
write.csv(otu_table(lg_intestine_epithelium_Selenomonadaceae_genus_filt), "lg_intestine_epithelium_Selenomonadaceae_otus.csv")
lg_intestine_epithelium_Selenomonadaceae_genus_melt <- psmelt(lg_intestine_epithelium_Selenomonadaceae_genus_filt)

lg_intestine_epithelium_Selenomonadaceae_palette <- c("#4B0082","#6030A8","grey88")
lg_intestine_epithelium_Selenomonadaceae_melt <-  psmelt(lg_intestine_epithelium_Selenomonadaceae_genus) # 37 genera
lg_intestine_epithelium_Selenomonadaceae_family <- subset_taxa(lg_intestine_epithelium_family, Family=="Selenomonadaceae") %>%
  psmelt()

ggplot() + theme_bw() + labs(y= "Selenomonadaceae (RA%)", title = "EPITHELIUM") +
  coord_flip() + theme_bw() +
  geom_bar(lg_intestine_epithelium_Selenomonadaceae_genus_melt, colour = "black", linewidth = 0.75,
           mapping = aes(x= treatment, y= Abundance, fill=Genus), stat = "summary") +
  scale_fill_manual(values = lg_intestine_epithelium_Selenomonadaceae_palette) +
  scale_x_discrete(limits = c("tylosin","control"), labels =c("TYLOSIN","NO TYLOSIN")) +
  scale_y_continuous(expand = c(0.001,0,45.3,0)) +
  geom_errorbar(lg_intestine_epithelium_Selenomonadaceae_family, mapping= aes(x= treatment, y= Abundance),
                stat = "summary", width = 0.5, linewidth = 0.75) +
  theme(legend.position = "none",
        plot.title = element_text(size = 44),
        panel.border = element_rect(colour = "black", size = 1.25),
        panel.grid.major.y = element_blank(),
        axis.ticks = element_line(size = 0.75, colour= "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 34, colour = "black"),
        axis.title.x = element_text(size =44),
        axis.text.x = element_text(size = 16, colour = "black"))

kruskal.test(lg_intestine_epithelium_Selenomonadaceae_family$Abundance, lg_intestine_epithelium_Selenomonadaceae_family$treatment) #NS

### METHANOBACTERIACEAE ####
## RUMEN LUMEN
rumen_lumen_Methanobacteriaceae_genus <- subset_taxa(rumen_lumen_genus, Family=="Methanobacteriaceae")
rumen_lumen_Methanobacteriaceae_genus_filt <- merge_low_abundance(rumen_lumen_Methanobacteriaceae_genus, threshold = 0.1) # 7 taxa
write.csv(tax_table(rumen_lumen_Methanobacteriaceae_genus_filt), "rumen_lumen_Methanobacteriaceae_taxa.csv")
write.csv(otu_table(rumen_lumen_Methanobacteriaceae_genus_filt), "rumen_lumen_Methanobacteriaceae_otus.csv")
rumen_lumen_Methanobacteriaceae_genus_melt <- psmelt(rumen_lumen_Methanobacteriaceae_genus_filt)

rumen_lumen_Methanobacteriaceae_palette <- c("firebrick3","red4","indianred1")
rumen_lumen_Methanobacteriaceae_family <- subset_taxa(rumen_lumen_family, Family=="Methanobacteriaceae") %>%
  psmelt()

ggplot() + theme_bw() + labs(y= "Methanobacteriaceae (RA%)", title = "LUMEN") +
  coord_flip() + theme_bw() +
  geom_bar(rumen_lumen_Methanobacteriaceae_genus_melt, colour = "black", linewidth = 0.75,
           mapping = aes(x= treatment, y= Abundance, fill=Genus), stat = "summary") +
  scale_fill_manual(values = rumen_lumen_Methanobacteriaceae_palette) +
  scale_x_discrete(limits = c("tylosin","control"), labels =c("TYLOSIN","NO TYLOSIN")) +
  scale_y_continuous(expand = c(0.001,0,6.41,0)) +
  geom_errorbar(rumen_lumen_Methanobacteriaceae_family, mapping= aes(x= treatment, y= Abundance),
                stat = "summary", width = 0.5, linewidth = 0.75) +
  theme(#legend.position = "none",
        plot.title = element_text(size = 44),
        panel.border = element_rect(colour = "black", size = 1.25),
        panel.grid.major.y = element_blank(),
        axis.ticks = element_line(size = 0.75, colour= "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 34, colour = "black"),
        axis.title.x = element_text(size = 44),
        axis.text.x = element_text(size = 16, colour = "black"))

kruskal.test(rumen_lumen_Methanobacteriaceae_family$Abundance, rumen_lumen_Methanobacteriaceae_family$treatment) #NS

## RUMEN EPITHELIUM
rumen_epithelium_Methanobacteriaceae_genus <- subset_taxa(rumen_epithelium_genus, Family=="Methanobacteriaceae")
rumen_epithelium_Methanobacteriaceae_genus_filt <- merge_low_abundance(rumen_epithelium_Methanobacteriaceae_genus, threshold = 0.1) # 7 taxa
write.csv(tax_table(rumen_epithelium_Methanobacteriaceae_genus_filt), "rumen_epithelium_Methanobacteriaceae_taxa.csv")
write.csv(otu_table(rumen_epithelium_Methanobacteriaceae_genus_filt), "rumen_epithelium_Methanobacteriaceae_otus.csv")
rumen_epithelium_Methanobacteriaceae_genus_melt <- psmelt(rumen_epithelium_Methanobacteriaceae_genus_filt)

rumen_epithelium_Methanobacteriaceae_palette <- c("red4","grey88")
rumen_epithelium_Methanobacteriaceae_family <- subset_taxa(rumen_epithelium_family, Family=="Methanobacteriaceae") %>%
  psmelt()

ggplot() + theme_bw() + labs(y= "Methanobacteriaceae (RA%)", title = "EPITHELIUM") +
  coord_flip() + theme_bw() +
  geom_bar(rumen_epithelium_Methanobacteriaceae_genus_melt, colour = "black", linewidth = 0.75,
           mapping = aes(x= treatment, y= Abundance, fill=Genus), stat = "summary") +
  scale_fill_manual(values = rumen_epithelium_Methanobacteriaceae_palette) +
  scale_x_discrete(limits = c("tylosin","control"), labels =c("TYLOSIN","NO TYLOSIN")) +
  scale_y_continuous(expand = c(0.001,0,20.353,0)) +
  geom_errorbar(rumen_epithelium_Methanobacteriaceae_family, mapping= aes(x= treatment, y= Abundance),
                stat = "summary", width = 0.5, linewidth = 0.75) +
  theme(#legend.position = "none",
        plot.title = element_text(size = 44),
        panel.border = element_rect(colour = "black", size = 1.25),
        panel.grid.major.y = element_blank(),
        axis.ticks = element_line(size = 0.75, colour= "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 34, colour = "black"),
        axis.title.x = element_text(size =44),
        axis.text.x = element_text(size = 16, colour = "black"))

kruskal.test(rumen_epithelium_Methanobacteriaceae_family$Abundance, rumen_epithelium_Methanobacteriaceae_family$treatment) #NS

### SM INTESTINE LUMEN
sm_intestine_lumen_Methanobacteriaceae_genus <- subset_taxa(sm_intestine_lumen_genus, Family=="Methanobacteriaceae")
sm_intestine_lumen_Methanobacteriaceae_genus_filt <- merge_low_abundance(sm_intestine_lumen_Methanobacteriaceae_genus, threshold = 0.1) # 7 taxa
write.csv(tax_table(sm_intestine_lumen_Methanobacteriaceae_genus_filt), "sm_intestine_lumen_Methanobacteriaceae_taxa.csv")
write.csv(otu_table(sm_intestine_lumen_Methanobacteriaceae_genus_filt), "sm_intestine_lumen_Methanobacteriaceae_otus.csv")
sm_intestine_lumen_Methanobacteriaceae_genus_melt <- psmelt(sm_intestine_lumen_Methanobacteriaceae_genus_filt)

sm_intestine_lumen_Methanobacteriaceae_palette <- c("red4","indianred1")
sm_intestine_lumen_Methanobacteriaceae_family <- subset_taxa(sm_intestine_lumen_family, Family=="Methanobacteriaceae") %>%
  psmelt()

ggplot() + theme_bw() + labs(y= "Methanobacteriaceae (RA%)", title = "LUMEN") +
  coord_flip() + theme_bw() +
  geom_bar(sm_intestine_lumen_Methanobacteriaceae_genus_melt, colour = "black", linewidth = 0.75,
           mapping = aes(x= treatment, y= Abundance, fill=Genus), stat = "summary") +
  scale_fill_manual(values = sm_intestine_lumen_Methanobacteriaceae_palette) +
  scale_x_discrete(limits = c("tylosin","control"), labels =c("TYLOSIN","NO TYLOSIN")) +
  scale_y_continuous(expand = c(0.001,0,4.6,0)) +
  geom_errorbar(sm_intestine_lumen_Methanobacteriaceae_family, mapping= aes(x= treatment, y= Abundance),
                stat = "summary", width = 0.5, linewidth = 0.75) +
  theme(legend.position = "none",
        plot.title = element_text(size = 44),
        panel.border = element_rect(colour = "black", size = 1.25),
        panel.grid.major.y = element_blank(),
        axis.ticks = element_line(size = 0.75, colour= "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 34, colour = "black"),
        axis.title.x = element_text(size =44),
        axis.text.x = element_text(size = 16, colour = "black"))

kruskal.test(sm_intestine_lumen_Methanobacteriaceae_family$Abundance, sm_intestine_lumen_Methanobacteriaceae_family$treatment) #NS

## sm_intestine EPITHELIUM
sm_intestine_epithelium_Methanobacteriaceae_genus <- subset_taxa(sm_intestine_epithelium_genus, Family=="Methanobacteriaceae")
sm_intestine_epithelium_Methanobacteriaceae_genus_filt <- merge_low_abundance(sm_intestine_epithelium_Methanobacteriaceae_genus, threshold = 0.1) # 7 taxa
write.csv(tax_table(sm_intestine_epithelium_Methanobacteriaceae_genus_filt), "sm_intestine_epithelium_Methanobacteriaceae_taxa.csv")
write.csv(otu_table(sm_intestine_epithelium_Methanobacteriaceae_genus_filt), "sm_intestine_epithelium_Methanobacteriaceae_otus.csv")
sm_intestine_epithelium_Methanobacteriaceae_genus_melt <- psmelt(sm_intestine_epithelium_Methanobacteriaceae_genus_filt)

sm_intestine_epithelium_Methanobacteriaceae_palette <- c("red4","indianred1")
sm_intestine_epithelium_Methanobacteriaceae_family <- subset_taxa(sm_intestine_epithelium_family, Family=="Methanobacteriaceae") %>%
  psmelt()

ggplot() + theme_bw() + labs(y= "Methanobacteriaceae (RA%)", title = "EPITHELIUM") +
  coord_flip() + theme_bw() +
  geom_bar(sm_intestine_epithelium_Methanobacteriaceae_genus_melt, colour = "black", linewidth = 0.75,
           mapping = aes(x= treatment, y= Abundance, fill=Genus), stat = "summary") +
  scale_fill_manual(values = sm_intestine_epithelium_Methanobacteriaceae_palette) +
  scale_x_discrete(limits = c("tylosin","control"), labels =c("TYLOSIN","NO TYLOSIN")) +
  scale_y_continuous(expand = c(0.001,0,4.5,0)) +
  geom_errorbar(sm_intestine_epithelium_Methanobacteriaceae_family, mapping= aes(x= treatment, y= Abundance),
                stat = "summary", width = 0.5, linewidth = 0.75) +
  theme(#legend.position = "none",
        plot.title = element_text(size = 44),
        panel.border = element_rect(colour = "black", size = 1.25),
        panel.grid.major.y = element_blank(),
        axis.ticks = element_line(size = 0.75, colour= "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 34, colour = "black"),
        axis.title.x = element_text(size =44),
        axis.text.x = element_text(size = 16, colour = "black"))

kruskal.test(sm_intestine_epithelium_Methanobacteriaceae_family$Abundance, sm_intestine_epithelium_Methanobacteriaceae_family$treatment) #NS

### LG INT LUMEN
lg_intestine_lumen_Methanobacteriaceae_genus <- subset_taxa(lg_intestine_lumen_genus, Family=="Methanobacteriaceae")
lg_intestine_lumen_Methanobacteriaceae_genus_filt <- merge_low_abundance(lg_intestine_lumen_Methanobacteriaceae_genus, threshold = 0.1) # 7 taxa
write.csv(tax_table(lg_intestine_lumen_Methanobacteriaceae_genus_filt), "lg_intestine_lumen_Methanobacteriaceae_taxa.csv")
write.csv(otu_table(lg_intestine_lumen_Methanobacteriaceae_genus_filt), "lg_intestine_lumen_Methanobacteriaceae_otus.csv")
lg_intestine_lumen_Methanobacteriaceae_genus_melt <- psmelt(lg_intestine_lumen_Methanobacteriaceae_genus_filt)

lg_intestine_lumen_Methanobacteriaceae_palette <- c("#6030A8","grey88")
lg_intestine_lumen_Methanobacteriaceae_family <- subset_taxa(lg_intestine_lumen_family, Family=="Methanobacteriaceae") %>%
  psmelt()

ggplot() + theme_bw() + labs(y= "Methanobacteriaceae (RA%)", title = "LUMEN") +
  coord_flip() + theme_bw() +
  geom_bar(lg_intestine_lumen_Methanobacteriaceae_genus_melt, colour = "black", linewidth = 0.75,
           mapping = aes(x= treatment, y= Abundance, fill=Genus), stat = "summary") +
  scale_fill_manual(values = lg_intestine_lumen_Methanobacteriaceae_palette) +
  scale_x_discrete(limits = c("tylosin","control"), labels =c("TYLOSIN","NO TYLOSIN")) +
  scale_y_continuous(expand = c(0.001,0,20,0)) +
  geom_errorbar(lg_intestine_lumen_Methanobacteriaceae_family, mapping= aes(x= treatment, y= Abundance),
                stat = "summary", width = 0.5, linewidth = 0.75) +
  theme(legend.position = "none",
        plot.title = element_text(size = 44),
        panel.border = element_rect(colour = "black", size = 1.25),
        panel.grid.major.y = element_blank(),
        axis.ticks = element_line(size = 0.75, colour= "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 34, colour = "black"),
        axis.title.x = element_text(size =44),
        axis.text.x = element_text(size = 16, colour = "black"))

kruskal.test(lg_intestine_lumen_Methanobacteriaceae_family$Abundance, lg_intestine_lumen_Methanobacteriaceae_family$treatment) #NS

## lg_intestine EPITHELIUM
lg_intestine_epithelium_Methanobacteriaceae_genus <- subset_taxa(lg_intestine_epithelium_genus, Family=="Methanobacteriaceae")
lg_intestine_epithelium_Methanobacteriaceae_genus_filt <- merge_low_abundance(lg_intestine_epithelium_Methanobacteriaceae_genus, threshold = 0.1) # 7 taxa
write.csv(tax_table(lg_intestine_epithelium_Methanobacteriaceae_genus_filt), "lg_intestine_epithelium_Methanobacteriaceae_taxa.csv")
write.csv(otu_table(lg_intestine_epithelium_Methanobacteriaceae_genus_filt), "lg_intestine_epithelium_Methanobacteriaceae_otus.csv")
lg_intestine_epithelium_Methanobacteriaceae_genus_melt <- psmelt(lg_intestine_epithelium_Methanobacteriaceae_genus_filt)

lg_intestine_epithelium_Methanobacteriaceae_palette <- c("#4B0082","#6030A8","grey88")
lg_intestine_epithelium_Methanobacteriaceae_family <- subset_taxa(lg_intestine_epithelium_family, Family=="Methanobacteriaceae") %>%
  psmelt()

ggplot() + theme_bw() + labs(y= "Methanobacteriaceae (RA%)", title = "EPITHELIUM") +
  coord_flip() + theme_bw() +
  geom_bar(lg_intestine_epithelium_Methanobacteriaceae_genus_melt, colour = "black", linewidth = 0.75,
           mapping = aes(x= treatment, y= Abundance, fill=Genus), stat = "summary") +
  scale_fill_manual(values = lg_intestine_epithelium_Methanobacteriaceae_palette) +
  scale_x_discrete(limits = c("tylosin","control"), labels =c("TYLOSIN","NO TYLOSIN")) +
  scale_y_continuous(expand = c(0.001,0,5.3,0)) +
  geom_errorbar(lg_intestine_epithelium_Methanobacteriaceae_family, mapping= aes(x= treatment, y= Abundance),
                stat = "summary", width = 0.5, linewidth = 0.75) +
  theme(legend.position = "none",
        plot.title = element_text(size = 44),
        panel.border = element_rect(colour = "black", size = 1.25),
        panel.grid.major.y = element_blank(),
        axis.ticks = element_line(size = 0.75, colour= "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 34, colour = "black"),
        axis.title.x = element_text(size =44),
        axis.text.x = element_text(size = 16, colour = "black"))

kruskal.test(lg_intestine_epithelium_Methanobacteriaceae_family$Abundance, lg_intestine_epithelium_Methanobacteriaceae_family$treatment) #NS

### Veillonellaceae ####
## RUMEN LUMEN
rumen_lumen_Veillonellaceae_genus <- subset_taxa(rumen_lumen_genus, Family=="Veillonellaceae")
rumen_lumen_Veillonellaceae_genus_filt <- merge_low_abundance(rumen_lumen_Veillonellaceae_genus, threshold = 0.1) # 7 taxa
write.csv(tax_table(rumen_lumen_Veillonellaceae_genus_filt), "rumen_lumen_Veillonellaceae_taxa.csv")
write.csv(otu_table(rumen_lumen_Veillonellaceae_genus_filt), "rumen_lumen_Veillonellaceae_otus.csv")
rumen_lumen_Veillonellaceae_genus_melt <- psmelt(rumen_lumen_Veillonellaceae_genus_filt)

rumen_lumen_Veillonellaceae_palette <- c("chartreuse2","grey88")
rumen_lumen_Veillonellaceae_family <- subset_taxa(rumen_lumen_family, Family=="Veillonellaceae") %>%
  psmelt()

ggplot() + theme_bw() + labs(y= "Veillonellaceae (RA%)", title = "LUMEN") +
  coord_flip() + theme_bw() +
  geom_bar(rumen_lumen_Veillonellaceae_genus_melt, colour = "black", linewidth = 0.75,
           mapping = aes(x= treatment, y= Abundance, fill=Genus), stat = "summary") +
  scale_fill_manual(values = rumen_lumen_Veillonellaceae_palette) +
  scale_x_discrete(limits = c("tylosin","control"), labels =c("TYLOSIN","NO TYLOSIN")) +
  scale_y_continuous(expand = c(0.001,0,34.4,0)) +
  geom_errorbar(rumen_lumen_Veillonellaceae_family, mapping= aes(x= treatment, y= Abundance),
                stat = "summary", width = 0.5, linewidth = 0.75) +
  theme(legend.position = "none",
        plot.title = element_text(size = 44),
        panel.border = element_rect(colour = "black", size = 1.25),
        panel.grid.major.y = element_blank(),
        axis.ticks = element_line(size = 0.75, colour= "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 34, colour = "black"),
        axis.title.x = element_text(size = 44),
        axis.text.x = element_text(size = 16, colour = "black"))

kruskal.test(rumen_lumen_Veillonellaceae_family$Abundance, rumen_lumen_Veillonellaceae_family$treatment) #NS

## RUMEN EPITHELIUM
rumen_epithelium_Veillonellaceae_genus <- subset_taxa(rumen_epithelium_genus, Family=="Veillonellaceae")
rumen_epithelium_Veillonellaceae_genus_filt <- merge_low_abundance(rumen_epithelium_Veillonellaceae_genus, threshold = 0.1) # 7 taxa
write.csv(tax_table(rumen_epithelium_Veillonellaceae_genus_filt), "rumen_epithelium_Veillonellaceae_taxa.csv")
write.csv(otu_table(rumen_epithelium_Veillonellaceae_genus_filt), "rumen_epithelium_Veillonellaceae_otus.csv")
rumen_epithelium_Veillonellaceae_genus_melt <- psmelt(rumen_epithelium_Veillonellaceae_genus_filt)

rumen_epithelium_Veillonellaceae_palette <- c("chartreuse2","chartreuse4","#b2fea5")
rumen_epithelium_Veillonellaceae_family <- subset_taxa(rumen_epithelium_family, Family=="Veillonellaceae") %>%
  psmelt()

ggplot() + theme_bw() + labs(y= "Veillonellaceae (RA%)", title = "EPITHELIUM") +
  coord_flip() + theme_bw() +
  geom_bar(rumen_epithelium_Veillonellaceae_genus_melt, colour = "black", linewidth = 0.75,
           mapping = aes(x= treatment, y= Abundance, fill=Genus), stat = "summary") +
  scale_fill_manual(values = rumen_epithelium_Veillonellaceae_palette) +
  scale_x_discrete(limits = c("tylosin","control"), labels =c("TYLOSIN","NO TYLOSIN")) +
  scale_y_continuous(expand = c(0.001,0,9.59,0)) +
  geom_errorbar(rumen_epithelium_Veillonellaceae_family, mapping= aes(x= treatment, y= Abundance),
                stat = "summary", width = 0.5, linewidth = 0.75) +
  theme(legend.position = "none",
        plot.title = element_text(size = 44),
        panel.border = element_rect(colour = "black", size = 1.25),
        panel.grid.major.y = element_blank(),
        axis.ticks = element_line(size = 0.75, colour= "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 34, colour = "black"),
        axis.title.x = element_text(size =44),
        axis.text.x = element_text(size = 16, colour = "black"))

kruskal.test(rumen_epithelium_Veillonellaceae_family$Abundance, rumen_epithelium_Veillonellaceae_family$treatment) #NS


### Ruminococcaceae ####
## RUMEN LUMEN
rumen_lumen_Ruminococcaceae_genus <- subset_taxa(rumen_lumen_genus, Family=="Ruminococcaceae")
rumen_lumen_Ruminococcaceae_genus_filt <- merge_low_abundance(rumen_lumen_Ruminococcaceae_genus, threshold = 0.1) # 7 taxa
write.csv(tax_table(rumen_lumen_Ruminococcaceae_genus_filt), "rumen_lumen_Ruminococcaceae_taxa.csv")
write.csv(otu_table(rumen_lumen_Ruminococcaceae_genus_filt), "rumen_lumen_Ruminococcaceae_otus.csv")
rumen_lumen_Ruminococcaceae_genus_melt <- psmelt(rumen_lumen_Ruminococcaceae_genus_filt)

rumen_lumen_Ruminococcaceae_palette <- c("firebrick2",	"#901818",	"#ffb7b7",	"#d6a49e",	"grey88")
rumen_lumen_Ruminococcaceae_family <- subset_taxa(rumen_lumen_family, Family=="Ruminococcaceae") %>%
  psmelt()

ggplot() + theme_bw() + labs(y= "Ruminococcaceae (RA%)", title = "LUMEN") +
  coord_flip() + theme_bw() +
  geom_bar(rumen_lumen_Ruminococcaceae_genus_melt, colour = "black", linewidth = 0.75,
           mapping = aes(x= treatment, y= Abundance, fill=Genus), stat = "summary") +
  scale_fill_manual(values = rumen_lumen_Ruminococcaceae_palette) +
  scale_x_discrete(limits = c("tylosin","control"), labels =c("TYLOSIN","NO TYLOSIN")) +
  scale_y_continuous(expand = c(0.001,0,8.74,0)) +
  geom_errorbar(rumen_lumen_Ruminococcaceae_family, mapping= aes(x= treatment, y= Abundance),
                stat = "summary", width = 0.5, linewidth = 0.75) +
  theme(#legend.position = "none",
        plot.title = element_text(size = 44),
        panel.border = element_rect(colour = "black", size = 1.25),
        panel.grid.major.y = element_blank(),
        axis.ticks = element_line(size = 0.75, colour= "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 34, colour = "black"),
        axis.title.x = element_text(size = 44),
        axis.text.x = element_text(size = 16, colour = "black"))

kruskal.test(rumen_lumen_Ruminococcaceae_family$Abundance, rumen_lumen_Ruminococcaceae_family$treatment) #NS

## RUMEN EPITHELIUM
rumen_epithelium_Ruminococcaceae_genus <- subset_taxa(rumen_epithelium_genus, Family=="Ruminococcaceae")
rumen_epithelium_Ruminococcaceae_genus_filt <- merge_low_abundance(rumen_epithelium_Ruminococcaceae_genus, threshold = 0.1) # 7 taxa
write.csv(tax_table(rumen_epithelium_Ruminococcaceae_genus_filt), "rumen_epithelium_Ruminococcaceae_taxa.csv")
write.csv(otu_table(rumen_epithelium_Ruminococcaceae_genus_filt), "rumen_epithelium_Ruminococcaceae_otus.csv")
rumen_epithelium_Ruminococcaceae_genus_melt <- psmelt(rumen_epithelium_Ruminococcaceae_genus_filt)

rumen_epithelium_Ruminococcaceae_palette <- c("firebrick2",	"#901818",	"#ffb7b7",	"#d6a49e",	"grey88")
rumen_epithelium_Ruminococcaceae_family <- subset_taxa(rumen_epithelium_family, Family=="Ruminococcaceae") %>%
  psmelt()

ggplot() + theme_bw() + labs(y= "Ruminococcaceae (RA%)", title = "EPITHELIUM") +
  coord_flip() + theme_bw() +
  geom_bar(rumen_epithelium_Ruminococcaceae_genus_melt, colour = "black", linewidth = 0.75,
           mapping = aes(x= treatment, y= Abundance, fill=Genus), stat = "summary") +
  scale_fill_manual(values = rumen_epithelium_Ruminococcaceae_palette) +
  scale_x_discrete(limits = c("tylosin","control"), labels =c("TYLOSIN","NO TYLOSIN")) +
  scale_y_continuous(expand = c(0.001,0,10.69,0)) +
  geom_errorbar(rumen_epithelium_Ruminococcaceae_family, mapping= aes(x= treatment, y= Abundance),
                stat = "summary", width = 0.5, linewidth = 0.75) +
  theme(legend.position = "none",
        plot.title = element_text(size = 44),
        panel.border = element_rect(colour = "black", size = 1.25),
        panel.grid.major.y = element_blank(),
        axis.ticks = element_line(size = 0.75, colour= "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 34, colour = "black"),
        axis.title.x = element_text(size =44),
        axis.text.x = element_text(size = 16, colour = "black"))

kruskal.test(rumen_epithelium_Ruminococcaceae_family$Abundance, rumen_epithelium_Ruminococcaceae_family$treatment) #NS

### SM INTESTINE LUMEN
sm_intestine_lumen_Ruminococcaceae_genus <- subset_taxa(sm_intestine_lumen_genus, Family=="Ruminococcaceae")
sm_intestine_lumen_Ruminococcaceae_genus_filt <- merge_low_abundance(sm_intestine_lumen_Ruminococcaceae_genus, threshold = 0.1) # 7 taxa
write.csv(tax_table(sm_intestine_lumen_Ruminococcaceae_genus_filt), "sm_intestine_lumen_Ruminococcaceae_taxa.csv")
write.csv(otu_table(sm_intestine_lumen_Ruminococcaceae_genus_filt), "sm_intestine_lumen_Ruminococcaceae_otus.csv")
sm_intestine_lumen_Ruminococcaceae_genus_melt <- psmelt(sm_intestine_lumen_Ruminococcaceae_genus_filt)

sm_intestine_lumen_Ruminococcaceae_palette <- c("red4","indianred1")
sm_intestine_lumen_Ruminococcaceae_family <- subset_taxa(sm_intestine_lumen_family, Family=="Ruminococcaceae") %>%
  psmelt()

ggplot() + theme_bw() + labs(y= "Ruminococcaceae (RA%)", title = "LUMEN") +
  coord_flip() + theme_bw() +
  geom_bar(sm_intestine_lumen_Ruminococcaceae_genus_melt, colour = "black", linewidth = 0.75,
           mapping = aes(x= treatment, y= Abundance, fill=Genus), stat = "summary") +
  scale_fill_manual(values = sm_intestine_lumen_Ruminococcaceae_palette) +
  scale_x_discrete(limits = c("tylosin","control"), labels =c("TYLOSIN","NO TYLOSIN")) +
  scale_y_continuous(expand = c(0.001,0,4.6,0)) +
  geom_errorbar(sm_intestine_lumen_Ruminococcaceae_family, mapping= aes(x= treatment, y= Abundance),
                stat = "summary", width = 0.5, linewidth = 0.75) +
  theme(legend.position = "none",
        plot.title = element_text(size = 44),
        panel.border = element_rect(colour = "black", size = 1.25),
        panel.grid.major.y = element_blank(),
        axis.ticks = element_line(size = 0.75, colour= "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 34, colour = "black"),
        axis.title.x = element_text(size =44),
        axis.text.x = element_text(size = 16, colour = "black"))

kruskal.test(sm_intestine_lumen_Ruminococcaceae_family$Abundance, sm_intestine_lumen_Ruminococcaceae_family$treatment) #NS

## sm_intestine EPITHELIUM
sm_intestine_epithelium_Ruminococcaceae_genus <- subset_taxa(sm_intestine_epithelium_genus, Family=="Ruminococcaceae")
sm_intestine_epithelium_Ruminococcaceae_genus_filt <- merge_low_abundance(sm_intestine_epithelium_Ruminococcaceae_genus, threshold = 0.1) # 7 taxa
write.csv(tax_table(sm_intestine_epithelium_Ruminococcaceae_genus_filt), "sm_intestine_epithelium_Ruminococcaceae_taxa.csv")
write.csv(otu_table(sm_intestine_epithelium_Ruminococcaceae_genus_filt), "sm_intestine_epithelium_Ruminococcaceae_otus.csv")
sm_intestine_epithelium_Ruminococcaceae_genus_melt <- psmelt(sm_intestine_epithelium_Ruminococcaceae_genus_filt)

sm_intestine_epithelium_Ruminococcaceae_palette <- c("red4","indianred1")
sm_intestine_epithelium_Ruminococcaceae_family <- subset_taxa(sm_intestine_epithelium_family, Family=="Ruminococcaceae") %>%
  psmelt()

ggplot() + theme_bw() + labs(y= "Ruminococcaceae (RA%)", title = "EPITHELIUM") +
  coord_flip() + theme_bw() +
  geom_bar(sm_intestine_epithelium_Ruminococcaceae_genus_melt, colour = "black", linewidth = 0.75,
           mapping = aes(x= treatment, y= Abundance, fill=Genus), stat = "summary") +
  #scale_fill_manual(values = sm_intestine_epithelium_Ruminococcaceae_palette) +
  scale_x_discrete(limits = c("tylosin","control"), labels =c("TYLOSIN","NO TYLOSIN")) +
  scale_y_continuous(expand = c(0.001,0,4.5,0)) +
  geom_errorbar(sm_intestine_epithelium_Ruminococcaceae_family, mapping= aes(x= treatment, y= Abundance),
                stat = "summary", width = 0.5, linewidth = 0.75) +
  theme(#legend.position = "none",
    plot.title = element_text(size = 44),
    panel.border = element_rect(colour = "black", size = 1.25),
    panel.grid.major.y = element_blank(),
    axis.ticks = element_line(size = 0.75, colour= "black"),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 34, colour = "black"),
    axis.title.x = element_text(size =44),
    axis.text.x = element_text(size = 16, colour = "black"))

kruskal.test(sm_intestine_epithelium_Ruminococcaceae_family$Abundance, sm_intestine_epithelium_Ruminococcaceae_family$treatment) #NS

### LG INT LUMEN
lg_intestine_lumen_Ruminococcaceae_genus <- subset_taxa(lg_intestine_lumen_genus, Family=="Ruminococcaceae")
lg_intestine_lumen_Ruminococcaceae_genus_filt <- merge_low_abundance(lg_intestine_lumen_Ruminococcaceae_genus, threshold = 0.1) # 7 taxa
write.csv(tax_table(lg_intestine_lumen_Ruminococcaceae_genus_filt), "lg_intestine_lumen_Ruminococcaceae_taxa.csv")
write.csv(otu_table(lg_intestine_lumen_Ruminococcaceae_genus_filt), "lg_intestine_lumen_Ruminococcaceae_otus.csv")
lg_intestine_lumen_Ruminococcaceae_genus_melt <- psmelt(lg_intestine_lumen_Ruminococcaceae_genus_filt)

lg_intestine_lumen_Ruminococcaceae_palette <- c("#6030A8","grey88")
lg_intestine_lumen_Ruminococcaceae_family <- subset_taxa(lg_intestine_lumen_family, Family=="Ruminococcaceae") %>%
  psmelt()

ggplot() + theme_bw() + labs(y= "Ruminococcaceae (RA%)", title = "LUMEN") +
  coord_flip() + theme_bw() +
  geom_bar(lg_intestine_lumen_Ruminococcaceae_genus_melt, colour = "black", linewidth = 0.75,
           mapping = aes(x= treatment, y= Abundance, fill=Genus), stat = "summary") +
  #scale_fill_manual(values = lg_intestine_lumen_Ruminococcaceae_palette) +
  scale_x_discrete(limits = c("tylosin","control"), labels =c("TYLOSIN","NO TYLOSIN")) +
  #scale_y_continuous(expand = c(0.001,0,20,0)) +
  geom_errorbar(lg_intestine_lumen_Ruminococcaceae_family, mapping= aes(x= treatment, y= Abundance),
                stat = "summary", width = 0.5, linewidth = 0.75) +
  theme(#legend.position = "none",
        plot.title = element_text(size = 44),
        panel.border = element_rect(colour = "black", size = 1.25),
        panel.grid.major.y = element_blank(),
        axis.ticks = element_line(size = 0.75, colour= "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 34, colour = "black"),
        axis.title.x = element_text(size =44),
        axis.text.x = element_text(size = 16, colour = "black"))

kruskal.test(lg_intestine_lumen_Ruminococcaceae_family$Abundance, lg_intestine_lumen_Ruminococcaceae_family$treatment) #NS

## lg_intestine EPITHELIUM
lg_intestine_epithelium_Ruminococcaceae_genus <- subset_taxa(lg_intestine_epithelium_genus, Family=="Ruminococcaceae")
lg_intestine_epithelium_Ruminococcaceae_genus_filt <- merge_low_abundance(lg_intestine_epithelium_Ruminococcaceae_genus, threshold = 0.1) # 7 taxa
write.csv(tax_table(lg_intestine_epithelium_Ruminococcaceae_genus_filt), "lg_intestine_epithelium_Ruminococcaceae_taxa.csv")
write.csv(otu_table(lg_intestine_epithelium_Ruminococcaceae_genus_filt), "lg_intestine_epithelium_Ruminococcaceae_otus.csv")
lg_intestine_epithelium_Ruminococcaceae_genus_melt <- psmelt(lg_intestine_epithelium_Ruminococcaceae_genus_filt)

lg_intestine_epithelium_Ruminococcaceae_palette <- c("#4B0082","#6030A8","grey88")
lg_intestine_epithelium_Ruminococcaceae_family <- subset_taxa(lg_intestine_epithelium_family, Family=="Ruminococcaceae") %>%
  psmelt()

ggplot() + theme_bw() + labs(y= "Ruminococcaceae (RA%)", title = "EPITHELIUM") +
  coord_flip() + theme_bw() +
  geom_bar(lg_intestine_epithelium_Ruminococcaceae_genus_melt, colour = "black", linewidth = 0.75,
           mapping = aes(x= treatment, y= Abundance, fill=Genus), stat = "summary") +
  #scale_fill_manual(values = lg_intestine_epithelium_Ruminococcaceae_palette) +
  scale_x_discrete(limits = c("tylosin","control"), labels =c("TYLOSIN","NO TYLOSIN")) +
  scale_y_continuous(expand = c(0.001,0,5.3,0)) +
  geom_errorbar(lg_intestine_epithelium_Ruminococcaceae_family, mapping= aes(x= treatment, y= Abundance),
                stat = "summary", width = 0.5, linewidth = 0.75) +
  theme(#legend.position = "none",
        plot.title = element_text(size = 44),
        panel.border = element_rect(colour = "black", size = 1.25),
        panel.grid.major.y = element_blank(),
        axis.ticks = element_line(size = 0.75, colour= "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 34, colour = "black"),
        axis.title.x = element_text(size =44),
        axis.text.x = element_text(size = 16, colour = "black"))

kruskal.test(lg_intestine_epithelium_Ruminococcaceae_family$Abundance, lg_intestine_epithelium_Ruminococcaceae_family$treatment) #NS

### Anaerovoracaceae ####
## RUMEN LUMEN
rumen_lumen_Anaerovoracaceae_genus <- subset_taxa(rumen_lumen_genus, Family=="Anaerovoracaceae")
rumen_lumen_Anaerovoracaceae_genus_filt <- merge_low_abundance(rumen_lumen_Anaerovoracaceae_genus, threshold = 0.1) # 7 taxa
write.csv(tax_table(rumen_lumen_Anaerovoracaceae_genus_filt), "rumen_lumen_Anaerovoracaceae_taxa.csv")
write.csv(otu_table(rumen_lumen_Anaerovoracaceae_genus_filt), "rumen_lumen_Anaerovoracaceae_otus.csv")
rumen_lumen_Anaerovoracaceae_genus_melt <- psmelt(rumen_lumen_Anaerovoracaceae_genus_filt)

rumen_lumen_Anaerovoracaceae_palette <- c("#7dffeb",	"#48CFCA",	"#01a6ac",	"#b4f8f5",	"grey88")
rumen_lumen_Anaerovoracaceae_family <- subset_taxa(rumen_lumen_family, Family=="Anaerovoracaceae") %>%
  psmelt()

ggplot() + theme_bw() + labs(y= "Anaerovoracaceae (RA%)", title = "LUMEN") +
  coord_flip() + theme_bw() +
  geom_bar(rumen_lumen_Anaerovoracaceae_genus_melt, colour = "black", linewidth = 0.75,
           mapping = aes(x= treatment, y= Abundance, fill=Genus), stat = "summary") +
  scale_fill_manual(values = rumen_lumen_Anaerovoracaceae_palette) +
  scale_x_discrete(limits = c("tylosin","control"), labels =c("TYLOSIN","NO TYLOSIN")) +
  scale_y_continuous(expand = c(0.001,0,8.6,0)) +
  geom_errorbar(rumen_lumen_Anaerovoracaceae_family, mapping= aes(x= treatment, y= Abundance),
                stat = "summary", width = 0.5, linewidth = 0.75) +
  theme(#legend.position = "none",
        plot.title = element_text(size = 44),
        panel.border = element_rect(colour = "black", size = 1.25),
        panel.grid.major.y = element_blank(),
        axis.ticks = element_line(size = 0.75, colour= "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 34, colour = "black"),
        axis.title.x = element_text(size = 44),
        axis.text.x = element_text(size = 16, colour = "black"))

kruskal.test(rumen_lumen_Anaerovoracaceae_family$Abundance, rumen_lumen_Anaerovoracaceae_family$treatment) #NS

## RUMEN EPITHELIUM
rumen_epithelium_Anaerovoracaceae_genus <- subset_taxa(rumen_epithelium_genus, Family=="Anaerovoracaceae")
rumen_epithelium_Anaerovoracaceae_genus_filt <- merge_low_abundance(rumen_epithelium_Anaerovoracaceae_genus, threshold = 0.1) # 7 taxa
write.csv(tax_table(rumen_epithelium_Anaerovoracaceae_genus_filt), "rumen_epithelium_Anaerovoracaceae_taxa.csv")
write.csv(otu_table(rumen_epithelium_Anaerovoracaceae_genus_filt), "rumen_epithelium_Anaerovoracaceae_otus.csv")
rumen_epithelium_Anaerovoracaceae_genus_melt <- psmelt(rumen_epithelium_Anaerovoracaceae_genus_filt)

rumen_epithelium_Anaerovoracaceae_palette <- c("#7dffeb",	"#48CFCA",	"#01a6ac",	"#b4f8f5",	"grey88")
rumen_epithelium_Anaerovoracaceae_family <- subset_taxa(rumen_epithelium_family, Family=="Anaerovoracaceae") %>%
  psmelt()

ggplot() + theme_bw() + labs(y= "Anaerovoracaceae (RA%)", title = "EPITHELIUM") +
  coord_flip() + theme_bw() +
  geom_bar(rumen_epithelium_Anaerovoracaceae_genus_melt, colour = "black", linewidth = 0.75,
           mapping = aes(x= treatment, y= Abundance, fill=Genus), stat = "summary") +
  scale_fill_manual(values = rumen_epithelium_Anaerovoracaceae_palette) +
  scale_x_discrete(limits = c("tylosin","control"), labels =c("TYLOSIN","NO TYLOSIN")) +
  scale_y_continuous(expand = c(0.001,0,13.4,0)) +
  geom_errorbar(rumen_epithelium_Anaerovoracaceae_family, mapping= aes(x= treatment, y= Abundance),
                stat = "summary", width = 0.5, linewidth = 0.75) +
  theme(legend.position = "none",
        plot.title = element_text(size = 44),
        panel.border = element_rect(colour = "black", size = 1.25),
        panel.grid.major.y = element_blank(),
        axis.ticks = element_line(size = 0.75, colour= "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 34, colour = "black"),
        axis.title.x = element_text(size =44),
        axis.text.x = element_text(size = 16, colour = "black"))

kruskal.test(rumen_epithelium_Anaerovoracaceae_family$Abundance, rumen_epithelium_Anaerovoracaceae_family$treatment) #NS

### SM INTESTINE LUMEN
sm_intestine_lumen_Anaerovoracaceae_genus <- subset_taxa(sm_intestine_lumen_genus, Family=="Anaerovoracaceae")
sm_intestine_lumen_Anaerovoracaceae_genus_filt <- merge_low_abundance(sm_intestine_lumen_Anaerovoracaceae_genus, threshold = 0.1) # 7 taxa
write.csv(tax_table(sm_intestine_lumen_Anaerovoracaceae_genus_filt), "sm_intestine_lumen_Anaerovoracaceae_taxa.csv")
write.csv(otu_table(sm_intestine_lumen_Anaerovoracaceae_genus_filt), "sm_intestine_lumen_Anaerovoracaceae_otus.csv")
sm_intestine_lumen_Anaerovoracaceae_genus_melt <- psmelt(sm_intestine_lumen_Anaerovoracaceae_genus_filt)

sm_intestine_lumen_Anaerovoracaceae_palette <- c("red4","indianred1")
sm_intestine_lumen_Anaerovoracaceae_family <- subset_taxa(sm_intestine_lumen_family, Family=="Anaerovoracaceae") %>%
  psmelt()

ggplot() + theme_bw() + labs(y= "Anaerovoracaceae (RA%)", title = "LUMEN") +
  coord_flip() + theme_bw() +
  geom_bar(sm_intestine_lumen_Anaerovoracaceae_genus_melt, colour = "black", linewidth = 0.75,
           mapping = aes(x= treatment, y= Abundance, fill=Genus), stat = "summary") +
  #scale_fill_manual(values = sm_intestine_lumen_Anaerovoracaceae_palette) +
  scale_x_discrete(limits = c("tylosin","control"), labels =c("TYLOSIN","NO TYLOSIN")) +
  scale_y_continuous(expand = c(0.001,0,4.6,0)) +
  geom_errorbar(sm_intestine_lumen_Anaerovoracaceae_family, mapping= aes(x= treatment, y= Abundance),
                stat = "summary", width = 0.5, linewidth = 0.75) +
  theme(#legend.position = "none",
        plot.title = element_text(size = 44),
        panel.border = element_rect(colour = "black", size = 1.25),
        panel.grid.major.y = element_blank(),
        axis.ticks = element_line(size = 0.75, colour= "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 34, colour = "black"),
        axis.title.x = element_text(size =44),
        axis.text.x = element_text(size = 16, colour = "black"))

kruskal.test(sm_intestine_lumen_Anaerovoracaceae_family$Abundance, sm_intestine_lumen_Anaerovoracaceae_family$treatment) #NS

## sm_intestine EPITHELIUM
sm_intestine_epithelium_Anaerovoracaceae_genus <- subset_taxa(sm_intestine_epithelium_genus, Family=="Anaerovoracaceae")
sm_intestine_epithelium_Anaerovoracaceae_genus_filt <- merge_low_abundance(sm_intestine_epithelium_Anaerovoracaceae_genus, threshold = 0.1) # 7 taxa
write.csv(tax_table(sm_intestine_epithelium_Anaerovoracaceae_genus_filt), "sm_intestine_epithelium_Anaerovoracaceae_taxa.csv")
write.csv(otu_table(sm_intestine_epithelium_Anaerovoracaceae_genus_filt), "sm_intestine_epithelium_Anaerovoracaceae_otus.csv")
sm_intestine_epithelium_Anaerovoracaceae_genus_melt <- psmelt(sm_intestine_epithelium_Anaerovoracaceae_genus_filt)

sm_intestine_epithelium_Anaerovoracaceae_palette <- c("red4","indianred1")
sm_intestine_epithelium_Anaerovoracaceae_family <- subset_taxa(sm_intestine_epithelium_family, Family=="Anaerovoracaceae") %>%
  psmelt()

ggplot() + theme_bw() + labs(y= "Anaerovoracaceae (RA%)", title = "EPITHELIUM") +
  coord_flip() + theme_bw() +
  geom_bar(sm_intestine_epithelium_Anaerovoracaceae_genus_melt, colour = "black", linewidth = 0.75,
           mapping = aes(x= treatment, y= Abundance, fill=Genus), stat = "summary") +
  #scale_fill_manual(values = sm_intestine_epithelium_Anaerovoracaceae_palette) +
  scale_x_discrete(limits = c("tylosin","control"), labels =c("TYLOSIN","NO TYLOSIN")) +
  scale_y_continuous(expand = c(0.001,0,4.5,0)) +
  geom_errorbar(sm_intestine_epithelium_Anaerovoracaceae_family, mapping= aes(x= treatment, y= Abundance),
                stat = "summary", width = 0.5, linewidth = 0.75) +
  theme(#legend.position = "none",
    plot.title = element_text(size = 44),
    panel.border = element_rect(colour = "black", size = 1.25),
    panel.grid.major.y = element_blank(),
    axis.ticks = element_line(size = 0.75, colour= "black"),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 34, colour = "black"),
    axis.title.x = element_text(size =44),
    axis.text.x = element_text(size = 16, colour = "black"))

kruskal.test(sm_intestine_epithelium_Anaerovoracaceae_family$Abundance, sm_intestine_epithelium_Anaerovoracaceae_family$treatment) #NS

### LG INT LUMEN
lg_intestine_lumen_Anaerovoracaceae_genus <- subset_taxa(lg_intestine_lumen_genus, Family=="Anaerovoracaceae")
lg_intestine_lumen_Anaerovoracaceae_genus_filt <- merge_low_abundance(lg_intestine_lumen_Anaerovoracaceae_genus, threshold = 0.1) # 7 taxa
write.csv(tax_table(lg_intestine_lumen_Anaerovoracaceae_genus_filt), "lg_intestine_lumen_Anaerovoracaceae_taxa.csv")
write.csv(otu_table(lg_intestine_lumen_Anaerovoracaceae_genus_filt), "lg_intestine_lumen_Anaerovoracaceae_otus.csv")
lg_intestine_lumen_Anaerovoracaceae_genus_melt <- psmelt(lg_intestine_lumen_Anaerovoracaceae_genus_filt)

lg_intestine_lumen_Anaerovoracaceae_palette <- c("#6030A8","grey88")
lg_intestine_lumen_Anaerovoracaceae_family <- subset_taxa(lg_intestine_lumen_family, Family=="Anaerovoracaceae") %>%
  psmelt()

ggplot() + theme_bw() + labs(y= "Anaerovoracaceae (RA%)", title = "LUMEN") +
  coord_flip() + theme_bw() +
  geom_bar(lg_intestine_lumen_Anaerovoracaceae_genus_melt, colour = "black", linewidth = 0.75,
           mapping = aes(x= treatment, y= Abundance, fill=Genus), stat = "summary") +
  scale_fill_manual(values = lg_intestine_lumen_Anaerovoracaceae_palette) +
  scale_x_discrete(limits = c("tylosin","control"), labels =c("TYLOSIN","NO TYLOSIN")) +
  scale_y_continuous(expand = c(0.001,0,20,0)) +
  geom_errorbar(lg_intestine_lumen_Anaerovoracaceae_family, mapping= aes(x= treatment, y= Abundance),
                stat = "summary", width = 0.5, linewidth = 0.75) +
  theme(legend.position = "none",
        plot.title = element_text(size = 44),
        panel.border = element_rect(colour = "black", size = 1.25),
        panel.grid.major.y = element_blank(),
        axis.ticks = element_line(size = 0.75, colour= "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 34, colour = "black"),
        axis.title.x = element_text(size =44),
        axis.text.x = element_text(size = 16, colour = "black"))

kruskal.test(lg_intestine_lumen_Anaerovoracaceae_family$Abundance, lg_intestine_lumen_Anaerovoracaceae_family$treatment) #NS

## lg_intestine EPITHELIUM
lg_intestine_epithelium_Anaerovoracaceae_genus <- subset_taxa(lg_intestine_epithelium_genus, Family=="Anaerovoracaceae")
lg_intestine_epithelium_Anaerovoracaceae_genus_filt <- merge_low_abundance(lg_intestine_epithelium_Anaerovoracaceae_genus, threshold = 0.1) # 7 taxa
write.csv(tax_table(lg_intestine_epithelium_Anaerovoracaceae_genus_filt), "lg_intestine_epithelium_Anaerovoracaceae_taxa.csv")
write.csv(otu_table(lg_intestine_epithelium_Anaerovoracaceae_genus_filt), "lg_intestine_epithelium_Anaerovoracaceae_otus.csv")
lg_intestine_epithelium_Anaerovoracaceae_genus_melt <- psmelt(lg_intestine_epithelium_Anaerovoracaceae_genus_filt)

lg_intestine_epithelium_Anaerovoracaceae_palette <- c("#4B0082","#6030A8","grey88")
lg_intestine_epithelium_Anaerovoracaceae_family <- subset_taxa(lg_intestine_epithelium_family, Family=="Anaerovoracaceae") %>%
  psmelt()

ggplot() + theme_bw() + labs(y= "Anaerovoracaceae (RA%)", title = "EPITHELIUM") +
  coord_flip() + theme_bw() +
  geom_bar(lg_intestine_epithelium_Anaerovoracaceae_genus_melt, colour = "black", linewidth = 0.75,
           mapping = aes(x= treatment, y= Abundance, fill=Genus), stat = "summary") +
  scale_fill_manual(values = lg_intestine_epithelium_Anaerovoracaceae_palette) +
  scale_x_discrete(limits = c("tylosin","control"), labels =c("TYLOSIN","NO TYLOSIN")) +
  scale_y_continuous(expand = c(0.001,0,5.3,0)) +
  geom_errorbar(lg_intestine_epithelium_Anaerovoracaceae_family, mapping= aes(x= treatment, y= Abundance),
                stat = "summary", width = 0.5, linewidth = 0.75) +
  theme(legend.position = "none",
        plot.title = element_text(size = 44),
        panel.border = element_rect(colour = "black", size = 1.25),
        panel.grid.major.y = element_blank(),
        axis.ticks = element_line(size = 0.75, colour= "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 34, colour = "black"),
        axis.title.x = element_text(size =44),
        axis.text.x = element_text(size = 16, colour = "black"))

kruskal.test(lg_intestine_epithelium_Anaerovoracaceae_family$Abundance, lg_intestine_epithelium_Anaerovoracaceae_family$treatment) #NS

### Muribaculaceae ####
## RUMEN LUMEN
rumen_lumen_Muribaculaceae_genus <- subset_taxa(rumen_lumen_genus, Family=="Muribaculaceae")
rumen_lumen_Muribaculaceae_genus_filt <- merge_low_abundance(rumen_lumen_Muribaculaceae_genus, threshold = 0.1) # 7 taxa
write.csv(tax_table(rumen_lumen_Muribaculaceae_genus_filt), "rumen_lumen_Muribaculaceae_taxa.csv")
write.csv(otu_table(rumen_lumen_Muribaculaceae_genus_filt), "rumen_lumen_Muribaculaceae_otus.csv")
rumen_lumen_Muribaculaceae_genus_melt <- psmelt(rumen_lumen_Muribaculaceae_genus_filt)

rumen_lumen_Muribaculaceae_palette <- c("#9ED6DB")
rumen_lumen_Muribaculaceae_family <- subset_taxa(rumen_lumen_family, Family=="Muribaculaceae") %>%
  psmelt()

ggplot() + theme_bw() + labs(y= "Muribaculaceae (RA%)", title = "LUMEN") +
  coord_flip() + theme_bw() +
  geom_bar(rumen_lumen_Muribaculaceae_genus_melt, colour = "black", linewidth = 0.75,
           mapping = aes(x= treatment, y= Abundance, fill=Genus), stat = "summary") +
  scale_fill_manual(values = rumen_lumen_Muribaculaceae_palette) +
  scale_x_discrete(limits = c("tylosin","control"), labels =c("TYLOSIN","NO TYLOSIN")) +
  scale_y_continuous(expand = c(0.001,0,9.33,0)) +
  geom_errorbar(rumen_lumen_Muribaculaceae_family, mapping= aes(x= treatment, y= Abundance),
                stat = "summary", width = 0.5, linewidth = 0.75) +
  theme(legend.position = "none",
        plot.title = element_text(size = 44),
        panel.border = element_rect(colour = "black", size = 1.25),
        panel.grid.major.y = element_blank(),
        axis.ticks = element_line(size = 0.75, colour= "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 34, colour = "black"),
        axis.title.x = element_text(size = 44),
        axis.text.x = element_text(size = 16, colour = "black"))

kruskal.test(rumen_lumen_Muribaculaceae_family$Abundance, rumen_lumen_Muribaculaceae_family$treatment) #NS

## RUMEN EPITHELIUM
rumen_epithelium_Muribaculaceae_genus <- subset_taxa(rumen_epithelium_genus, Family=="Muribaculaceae")
rumen_epithelium_Muribaculaceae_genus_filt <- merge_low_abundance(rumen_epithelium_Muribaculaceae_genus, threshold = 0.1) # 7 taxa
write.csv(tax_table(rumen_epithelium_Muribaculaceae_genus_filt), "rumen_epithelium_Muribaculaceae_taxa.csv")
write.csv(otu_table(rumen_epithelium_Muribaculaceae_genus_filt), "rumen_epithelium_Muribaculaceae_otus.csv")
rumen_epithelium_Muribaculaceae_genus_melt <- psmelt(rumen_epithelium_Muribaculaceae_genus_filt)

rumen_epithelium_Muribaculaceae_palette <- c("grey55","#9ED6DB")
rumen_epithelium_Muribaculaceae_family <- subset_taxa(rumen_epithelium_family, Family=="Muribaculaceae") %>%
  psmelt()

ggplot() + theme_bw() + labs(y= "Muribaculaceae (RA%)", title = "EPITHELIUM") +
  coord_flip() + theme_bw() +
  geom_bar(rumen_epithelium_Muribaculaceae_genus_melt, colour = "black", linewidth = 0.75,
           mapping = aes(x= treatment, y= Abundance, fill=Genus), stat = "summary") +
  scale_fill_manual(values = rumen_epithelium_Muribaculaceae_palette) +
  scale_x_discrete(limits = c("tylosin","control"), labels =c("TYLOSIN","NO TYLOSIN")) +
  scale_y_continuous(expand = c(0.001,0,11,0)) +
  geom_errorbar(rumen_epithelium_Muribaculaceae_family, mapping= aes(x= treatment, y= Abundance),
                stat = "summary", width = 0.5, linewidth = 0.75) +
  theme(legend.position = "none",
        plot.title = element_text(size = 44),
        panel.border = element_rect(colour = "black", size = 1.25),
        panel.grid.major.y = element_blank(),
        axis.ticks = element_line(size = 0.75, colour= "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 34, colour = "black"),
        axis.title.x = element_text(size =44),
        axis.text.x = element_text(size = 16, colour = "black"))

kruskal.test(rumen_epithelium_Muribaculaceae_family$Abundance, rumen_epithelium_Muribaculaceae_family$treatment) #NS

### SM INTESTINE LUMEN
sm_intestine_lumen_Muribaculaceae_genus <- subset_taxa(sm_intestine_lumen_genus, Family=="Muribaculaceae")
sm_intestine_lumen_Muribaculaceae_genus_filt <- merge_low_abundance(sm_intestine_lumen_Muribaculaceae_genus, threshold = 0.1) # 7 taxa
write.csv(tax_table(sm_intestine_lumen_Muribaculaceae_genus_filt), "sm_intestine_lumen_Muribaculaceae_taxa.csv")
write.csv(otu_table(sm_intestine_lumen_Muribaculaceae_genus_filt), "sm_intestine_lumen_Muribaculaceae_otus.csv")
sm_intestine_lumen_Muribaculaceae_genus_melt <- psmelt(sm_intestine_lumen_Muribaculaceae_genus_filt)

sm_intestine_lumen_Muribaculaceae_palette <- c("red4","indianred1")
sm_intestine_lumen_Muribaculaceae_family <- subset_taxa(sm_intestine_lumen_family, Family=="Muribaculaceae") %>%
  psmelt()

ggplot() + theme_bw() + labs(y= "Muribaculaceae (RA%)", title = "LUMEN") +
  coord_flip() + theme_bw() +
  geom_bar(sm_intestine_lumen_Muribaculaceae_genus_melt, colour = "black", linewidth = 0.75,
           mapping = aes(x= treatment, y= Abundance, fill=Genus), stat = "summary") +
  scale_fill_manual(values = sm_intestine_lumen_Muribaculaceae_palette) +
  scale_x_discrete(limits = c("tylosin","control"), labels =c("TYLOSIN","NO TYLOSIN")) +
  scale_y_continuous(expand = c(0.001,0,4.6,0)) +
  geom_errorbar(sm_intestine_lumen_Muribaculaceae_family, mapping= aes(x= treatment, y= Abundance),
                stat = "summary", width = 0.5, linewidth = 0.75) +
  theme(legend.position = "none",
        plot.title = element_text(size = 44),
        panel.border = element_rect(colour = "black", size = 1.25),
        panel.grid.major.y = element_blank(),
        axis.ticks = element_line(size = 0.75, colour= "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 34, colour = "black"),
        axis.title.x = element_text(size =44),
        axis.text.x = element_text(size = 16, colour = "black"))

kruskal.test(sm_intestine_lumen_Muribaculaceae_family$Abundance, sm_intestine_lumen_Muribaculaceae_family$treatment) #NS

## sm_intestine EPITHELIUM
sm_intestine_epithelium_Muribaculaceae_genus <- subset_taxa(sm_intestine_epithelium_genus, Family=="Muribaculaceae")
sm_intestine_epithelium_Muribaculaceae_genus_filt <- merge_low_abundance(sm_intestine_epithelium_Muribaculaceae_genus, threshold = 0.1) # 7 taxa
write.csv(tax_table(sm_intestine_epithelium_Muribaculaceae_genus_filt), "sm_intestine_epithelium_Muribaculaceae_taxa.csv")
write.csv(otu_table(sm_intestine_epithelium_Muribaculaceae_genus_filt), "sm_intestine_epithelium_Muribaculaceae_otus.csv")
sm_intestine_epithelium_Muribaculaceae_genus_melt <- psmelt(sm_intestine_epithelium_Muribaculaceae_genus_filt)

sm_intestine_epithelium_Muribaculaceae_palette <- c("red4","indianred1")
sm_intestine_epithelium_Muribaculaceae_family <- subset_taxa(sm_intestine_epithelium_family, Family=="Muribaculaceae") %>%
  psmelt()

ggplot() + theme_bw() + labs(y= "Muribaculaceae (RA%)", title = "EPITHELIUM") +
  coord_flip() + theme_bw() +
  geom_bar(sm_intestine_epithelium_Muribaculaceae_genus_melt, colour = "black", linewidth = 0.75,
           mapping = aes(x= treatment, y= Abundance, fill=Genus), stat = "summary") +
  scale_fill_manual(values = sm_intestine_epithelium_Muribaculaceae_palette) +
  scale_x_discrete(limits = c("tylosin","control"), labels =c("TYLOSIN","NO TYLOSIN")) +
  scale_y_continuous(expand = c(0.001,0,4.5,0)) +
  geom_errorbar(sm_intestine_epithelium_Muribaculaceae_family, mapping= aes(x= treatment, y= Abundance),
                stat = "summary", width = 0.5, linewidth = 0.75) +
  theme(#legend.position = "none",
    plot.title = element_text(size = 44),
    panel.border = element_rect(colour = "black", size = 1.25),
    panel.grid.major.y = element_blank(),
    axis.ticks = element_line(size = 0.75, colour= "black"),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 34, colour = "black"),
    axis.title.x = element_text(size =44),
    axis.text.x = element_text(size = 16, colour = "black"))

kruskal.test(sm_intestine_epithelium_Muribaculaceae_family$Abundance, sm_intestine_epithelium_Muribaculaceae_family$treatment) #NS

### LG INT LUMEN
lg_intestine_lumen_Muribaculaceae_genus <- subset_taxa(lg_intestine_lumen_genus, Family=="Muribaculaceae")
lg_intestine_lumen_Muribaculaceae_genus_filt <- merge_low_abundance(lg_intestine_lumen_Muribaculaceae_genus, threshold = 0.1) # 7 taxa
write.csv(tax_table(lg_intestine_lumen_Muribaculaceae_genus_filt), "lg_intestine_lumen_Muribaculaceae_taxa.csv")
write.csv(otu_table(lg_intestine_lumen_Muribaculaceae_genus_filt), "lg_intestine_lumen_Muribaculaceae_otus.csv")
lg_intestine_lumen_Muribaculaceae_genus_melt <- psmelt(lg_intestine_lumen_Muribaculaceae_genus_filt)

lg_intestine_lumen_Muribaculaceae_palette <- c("#6030A8","grey88")
lg_intestine_lumen_Muribaculaceae_family <- subset_taxa(lg_intestine_lumen_family, Family=="Muribaculaceae") %>%
  psmelt()

ggplot() + theme_bw() + labs(y= "Muribaculaceae (RA%)", title = "LUMEN") +
  coord_flip() + theme_bw() +
  geom_bar(lg_intestine_lumen_Muribaculaceae_genus_melt, colour = "black", linewidth = 0.75,
           mapping = aes(x= treatment, y= Abundance, fill=Genus), stat = "summary") +
  scale_fill_manual(values = lg_intestine_lumen_Muribaculaceae_palette) +
  scale_x_discrete(limits = c("tylosin","control"), labels =c("TYLOSIN","NO TYLOSIN")) +
  scale_y_continuous(expand = c(0.001,0,4,0)) +
  geom_errorbar(lg_intestine_lumen_Muribaculaceae_family, mapping= aes(x= treatment, y= Abundance),
                stat = "summary", width = 0.5, linewidth = 0.75) +
  theme(legend.position = "none",
        plot.title = element_text(size = 44),
        panel.border = element_rect(colour = "black", size = 1.25),
        panel.grid.major.y = element_blank(),
        axis.ticks = element_line(size = 0.75, colour= "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 34, colour = "black"),
        axis.title.x = element_text(size =44),
        axis.text.x = element_text(size = 16, colour = "black"))

kruskal.test(lg_intestine_lumen_Muribaculaceae_family$Abundance, lg_intestine_lumen_Muribaculaceae_family$treatment) #NS

## lg_intestine EPITHELIUM
lg_intestine_epithelium_Muribaculaceae_genus <- subset_taxa(lg_intestine_epithelium_genus, Family=="Muribaculaceae")
lg_intestine_epithelium_Muribaculaceae_genus_filt <- merge_low_abundance(lg_intestine_epithelium_Muribaculaceae_genus, threshold = 0.1) # 7 taxa
write.csv(tax_table(lg_intestine_epithelium_Muribaculaceae_genus_filt), "lg_intestine_epithelium_Muribaculaceae_taxa.csv")
write.csv(otu_table(lg_intestine_epithelium_Muribaculaceae_genus_filt), "lg_intestine_epithelium_Muribaculaceae_otus.csv")
lg_intestine_epithelium_Muribaculaceae_genus_melt <- psmelt(lg_intestine_epithelium_Muribaculaceae_genus_filt)

lg_intestine_epithelium_Muribaculaceae_palette <- c("#4B0082","#6030A8","grey88")
lg_intestine_epithelium_Muribaculaceae_family <- subset_taxa(lg_intestine_epithelium_family, Family=="Muribaculaceae") %>%
  psmelt()

ggplot() + theme_bw() + labs(y= "Muribaculaceae (RA%)", title = "EPITHELIUM") +
  coord_flip() + theme_bw() +
  geom_bar(lg_intestine_epithelium_Muribaculaceae_genus_melt, colour = "black", linewidth = 0.75,
           mapping = aes(x= treatment, y= Abundance, fill=Genus), stat = "summary") +
  scale_fill_manual(values = lg_intestine_epithelium_Muribaculaceae_palette) +
  scale_x_discrete(limits = c("tylosin","control"), labels =c("TYLOSIN","NO TYLOSIN")) +
  scale_y_continuous(expand = c(0.001,0,5.3,0)) +
  geom_errorbar(lg_intestine_epithelium_Muribaculaceae_family, mapping= aes(x= treatment, y= Abundance),
                stat = "summary", width = 0.5, linewidth = 0.75) +
  theme(legend.position = "none",
        plot.title = element_text(size = 44),
        panel.border = element_rect(colour = "black", size = 1.25),
        panel.grid.major.y = element_blank(),
        axis.ticks = element_line(size = 0.75, colour= "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 34, colour = "black"),
        axis.title.x = element_text(size =44),
        axis.text.x = element_text(size = 16, colour = "black"))

kruskal.test(lg_intestine_epithelium_Muribaculaceae_family$Abundance, lg_intestine_epithelium_Muribaculaceae_family$treatment) #NS

### Clostridia_UCG-014 ####
## RUMEN LUMEN
rumen_lumen_Clostridia_UCG014_genus <- subset_taxa(rumen_lumen_genus, Family=="Clostridia_UCG-014")
rumen_lumen_Clostridia_UCG014_genus_filt <- merge_low_abundance(rumen_lumen_Clostridia_UCG014_genus, threshold = 0.1) # 7 taxa
write.csv(tax_table(rumen_lumen_Clostridia_UCG014_genus_filt), "rumen_lumen_Clostridia_UCG-014_taxa.csv")
write.csv(otu_table(rumen_lumen_Clostridia_UCG014_genus_filt), "rumen_lumen_Clostridia_UCG-014_otus.csv")
rumen_lumen_Clostridia_UCG014_genus_melt <- psmelt(rumen_lumen_Clostridia_UCG014_genus_filt)

rumen_lumen_Clostridia_UCG014_palette <- c("#A3748D")
rumen_lumen_Clostridia_UCG014_family <- subset_taxa(rumen_lumen_family, Family=="Clostridia_UCG-014") %>%
  psmelt()

ggplot() + theme_bw() + labs(y= "Clostridia UCG-014 (RA%)", title = "LUMEN") +
  coord_flip() + theme_bw() +
  geom_bar(rumen_lumen_Clostridia_UCG014_genus_melt, colour = "black", linewidth = 0.75,
           mapping = aes(x= treatment, y= Abundance, fill=Genus), stat = "summary") +
  scale_fill_manual(values = rumen_lumen_Clostridia_UCG014_palette) +
  scale_x_discrete(limits = c("tylosin","control"), labels =c("TYLOSIN","NO TYLOSIN")) +
  scale_y_continuous(expand = c(0.001,0,11.9,0)) +
  geom_errorbar(rumen_lumen_Clostridia_UCG014_family, mapping= aes(x= treatment, y= Abundance),
                stat = "summary", width = 0.5, linewidth = 0.75) +
  theme(#legend.position = "none",
        plot.title = element_text(size = 44),
        panel.border = element_rect(colour = "black", size = 1.25),
        panel.grid.major.y = element_blank(),
        axis.ticks = element_line(size = 0.75, colour= "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 34, colour = "black"),
        axis.title.x = element_text(size = 44),
        axis.text.x = element_text(size = 16, colour = "black"))

kruskal.test(rumen_lumen_Clostridia_UCG014_family$Abundance, rumen_lumen_Clostridia_UCG014_family$treatment) #NS

## RUMEN EPITHELIUM
rumen_epithelium_Clostridia_UCG014_genus <- subset_taxa(rumen_epithelium_genus, Family=="Clostridia_UCG-014")
rumen_epithelium_Clostridia_UCG014_genus_filt <- merge_low_abundance(rumen_epithelium_Clostridia_UCG014_genus, threshold = 0.1) # 7 taxa
write.csv(tax_table(rumen_epithelium_Clostridia_UCG014_genus_filt), "rumen_epithelium_Clostridia_UCG-014_taxa.csv")
write.csv(otu_table(rumen_epithelium_Clostridia_UCG014_genus_filt), "rumen_epithelium_Clostridia_UCG-014_otus.csv")
rumen_epithelium_Clostridia_UCG014_genus_melt <- psmelt(rumen_epithelium_Clostridia_UCG014_genus_filt)

rumen_epithelium_Clostridia_UCG014_palette <- c("#A3748D")
rumen_epithelium_Clostridia_UCG014_family <- subset_taxa(rumen_epithelium_family, Family=="Clostridia_UCG-014") %>%
  psmelt()

ggplot() + theme_bw() + labs(y= "Clostridia UCG-014 (RA%)", title = "EPITHELIUM") +
  coord_flip() + theme_bw() +
  geom_bar(rumen_epithelium_Clostridia_UCG014_genus_melt, colour = "black", linewidth = 0.75,
           mapping = aes(x= treatment, y= Abundance, fill=Genus), stat = "summary") +
  scale_fill_manual(values = rumen_epithelium_Clostridia_UCG014_palette) +
  scale_x_discrete(limits = c("tylosin","control"), labels =c("TYLOSIN","NO TYLOSIN")) +
  scale_y_continuous(expand = c(0.001,0,29.29,0)) +
  geom_errorbar(rumen_epithelium_Clostridia_UCG014_family, mapping= aes(x= treatment, y= Abundance),
                stat = "summary", width = 0.5, linewidth = 0.75) +
  theme(legend.position = "none",
        plot.title = element_text(size = 44),
        panel.border = element_rect(colour = "black", size = 1.25),
        panel.grid.major.y = element_blank(),
        axis.ticks = element_line(size = 0.75, colour= "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 34, colour = "black"),
        axis.title.x = element_text(size =44),
        axis.text.x = element_text(size = 16, colour = "black"))

kruskal.test(rumen_epithelium_Clostridia_UCG014_family$Abundance, rumen_epithelium_Clostridia_UCG014_family$treatment) #NS

### SM INTESTINE LUMEN
sm_intestine_lumen_Clostridia_UCG014_genus <- subset_taxa(sm_intestine_lumen_genus, Family=="Clostridia_UCG-014")
sm_intestine_lumen_Clostridia_UCG014_genus_filt <- merge_low_abundance(sm_intestine_lumen_Clostridia_UCG014_genus, threshold = 0.1) # 7 taxa
write.csv(tax_table(sm_intestine_lumen_Clostridia_UCG014_genus_filt), "sm_intestine_lumen_Clostridia_UCG-014_taxa.csv")
write.csv(otu_table(sm_intestine_lumen_Clostridia_UCG014_genus_filt), "sm_intestine_lumen_Clostridia_UCG-014_otus.csv")
sm_intestine_lumen_Clostridia_UCG014_genus_melt <- psmelt(sm_intestine_lumen_Clostridia_UCG014_genus_filt)

sm_intestine_lumen_Clostridia_UCG014_palette <- c("red4","indianred1")
sm_intestine_lumen_Clostridia_UCG014_family <- subset_taxa(sm_intestine_lumen_family, Family=="Clostridia_UCG-014") %>%
  psmelt()

ggplot() + theme_bw() + labs(y= "Clostridia_UCG-014 (RA%)", title = "LUMEN") +
  coord_flip() + theme_bw() +
  geom_bar(sm_intestine_lumen_Clostridia_UCG014_genus_melt, colour = "black", linewidth = 0.75,
           mapping = aes(x= treatment, y= Abundance, fill=Genus), stat = "summary") +
  scale_fill_manual(values = sm_intestine_lumen_Clostridia_UCG014_palette) +
  scale_x_discrete(limits = c("tylosin","control"), labels =c("TYLOSIN","NO TYLOSIN")) +
  scale_y_continuous(expand = c(0.001,0,4.6,0)) +
  geom_errorbar(sm_intestine_lumen_Clostridia_UCG014_family, mapping= aes(x= treatment, y= Abundance),
                stat = "summary", width = 0.5, linewidth = 0.75) +
  theme(legend.position = "none",
        plot.title = element_text(size = 44),
        panel.border = element_rect(colour = "black", size = 1.25),
        panel.grid.major.y = element_blank(),
        axis.ticks = element_line(size = 0.75, colour= "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 34, colour = "black"),
        axis.title.x = element_text(size =44),
        axis.text.x = element_text(size = 16, colour = "black"))

kruskal.test(sm_intestine_lumen_Clostridia_UCG-014_family$Abundance, sm_intestine_lumen_Clostridia_UCG-014_family$treatment) #NS

## sm_intestine EPITHELIUM
sm_intestine_epithelium_Clostridia_UCG-014_genus <- subset_taxa(sm_intestine_epithelium_genus, Family=="Clostridia_UCG-014")
sm_intestine_epithelium_Clostridia_UCG-014_genus_filt <- merge_low_abundance(sm_intestine_epithelium_Clostridia_UCG-014_genus, threshold = 0.1) # 7 taxa
write.csv(tax_table(sm_intestine_epithelium_Clostridia_UCG-014_genus_filt), "sm_intestine_epithelium_Clostridia_UCG-014_taxa.csv")
write.csv(otu_table(sm_intestine_epithelium_Clostridia_UCG-014_genus_filt), "sm_intestine_epithelium_Clostridia_UCG-014_otus.csv")
sm_intestine_epithelium_Clostridia_UCG-014_genus_melt <- psmelt(sm_intestine_epithelium_Clostridia_UCG-014_genus_filt)

sm_intestine_epithelium_Clostridia_UCG-014_palette <- c("red4","indianred1")
sm_intestine_epithelium_Clostridia_UCG-014_family <- subset_taxa(sm_intestine_epithelium_family, Family=="Clostridia_UCG-014") %>%
  psmelt()

ggplot() + theme_bw() + labs(y= "Clostridia_UCG-014 (RA%)", title = "EPITHELIUM") +
  coord_flip() + theme_bw() +
  geom_bar(sm_intestine_epithelium_Clostridia_UCG-014_genus_melt, colour = "black", linewidth = 0.75,
           mapping = aes(x= treatment, y= Abundance, fill=Genus), stat = "summary") +
  scale_fill_manual(values = sm_intestine_epithelium_Clostridia_UCG-014_palette) +
  scale_x_discrete(limits = c("tylosin","control"), labels =c("TYLOSIN","NO TYLOSIN")) +
  scale_y_continuous(expand = c(0.001,0,4.5,0)) +
  geom_errorbar(sm_intestine_epithelium_Clostridia_UCG-014_family, mapping= aes(x= treatment, y= Abundance),
                stat = "summary", width = 0.5, linewidth = 0.75) +
  theme(#legend.position = "none",
    plot.title = element_text(size = 44),
    panel.border = element_rect(colour = "black", size = 1.25),
    panel.grid.major.y = element_blank(),
    axis.ticks = element_line(size = 0.75, colour= "black"),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 34, colour = "black"),
    axis.title.x = element_text(size =44),
    axis.text.x = element_text(size = 16, colour = "black"))

kruskal.test(sm_intestine_epithelium_Clostridia_UCG-014_family$Abundance, sm_intestine_epithelium_Clostridia_UCG-014_family$treatment) #NS

### LG INT LUMEN
lg_intestine_lumen_Clostridia_UCG-014_genus <- subset_taxa(lg_intestine_lumen_genus, Family=="Clostridia_UCG-014")
lg_intestine_lumen_Clostridia_UCG-014_genus_filt <- merge_low_abundance(lg_intestine_lumen_Clostridia_UCG-014_genus, threshold = 0.1) # 7 taxa
write.csv(tax_table(lg_intestine_lumen_Clostridia_UCG-014_genus_filt), "lg_intestine_lumen_Clostridia_UCG-014_taxa.csv")
write.csv(otu_table(lg_intestine_lumen_Clostridia_UCG-014_genus_filt), "lg_intestine_lumen_Clostridia_UCG-014_otus.csv")
lg_intestine_lumen_Clostridia_UCG-014_genus_melt <- psmelt(lg_intestine_lumen_Clostridia_UCG-014_genus_filt)

lg_intestine_lumen_Clostridia_UCG-014_palette <- c("#6030A8","grey88")
lg_intestine_lumen_Clostridia_UCG-014_family <- subset_taxa(lg_intestine_lumen_family, Family=="Clostridia_UCG-014") %>%
  psmelt()

ggplot() + theme_bw() + labs(y= "Clostridia_UCG-014 (RA%)", title = "LUMEN") +
  coord_flip() + theme_bw() +
  geom_bar(lg_intestine_lumen_Clostridia_UCG-014_genus_melt, colour = "black", linewidth = 0.75,
           mapping = aes(x= treatment, y= Abundance, fill=Genus), stat = "summary") +
  scale_fill_manual(values = lg_intestine_lumen_Clostridia_UCG-014_palette) +
  scale_x_discrete(limits = c("tylosin","control"), labels =c("TYLOSIN","NO TYLOSIN")) +
  scale_y_continuous(expand = c(0.001,0,20,0)) +
  geom_errorbar(lg_intestine_lumen_Clostridia_UCG-014_family, mapping= aes(x= treatment, y= Abundance),
                stat = "summary", width = 0.5, linewidth = 0.75) +
  theme(legend.position = "none",
        plot.title = element_text(size = 44),
        panel.border = element_rect(colour = "black", size = 1.25),
        panel.grid.major.y = element_blank(),
        axis.ticks = element_line(size = 0.75, colour= "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 34, colour = "black"),
        axis.title.x = element_text(size =44),
        axis.text.x = element_text(size = 16, colour = "black"))

kruskal.test(lg_intestine_lumen_Clostridia_UCG-014_family$Abundance, lg_intestine_lumen_Clostridia_UCG-014_family$treatment) #NS

## lg_intestine EPITHELIUM
lg_intestine_epithelium_Clostridia_UCG-014_genus <- subset_taxa(lg_intestine_epithelium_genus, Family=="Clostridia_UCG-014")
lg_intestine_epithelium_Clostridia_UCG-014_genus_filt <- merge_low_abundance(lg_intestine_epithelium_Clostridia_UCG-014_genus, threshold = 0.1) # 7 taxa
write.csv(tax_table(lg_intestine_epithelium_Clostridia_UCG-014_genus_filt), "lg_intestine_epithelium_Clostridia_UCG-014_taxa.csv")
write.csv(otu_table(lg_intestine_epithelium_Clostridia_UCG-014_genus_filt), "lg_intestine_epithelium_Clostridia_UCG-014_otus.csv")
lg_intestine_epithelium_Clostridia_UCG-014_genus_melt <- psmelt(lg_intestine_epithelium_Clostridia_UCG-014_genus_filt)

lg_intestine_epithelium_Clostridia_UCG-014_palette <- c("#4B0082","#6030A8","grey88")
lg_intestine_epithelium_Clostridia_UCG-014_family <- subset_taxa(lg_intestine_epithelium_family, Family=="Clostridia_UCG-014") %>%
  psmelt()

ggplot() + theme_bw() + labs(y= "Clostridia_UCG-014 (RA%)", title = "EPITHELIUM") +
  coord_flip() + theme_bw() +
  geom_bar(lg_intestine_epithelium_Clostridia_UCG-014_genus_melt, colour = "black", linewidth = 0.75,
           mapping = aes(x= treatment, y= Abundance, fill=Genus), stat = "summary") +
  scale_fill_manual(values = lg_intestine_epithelium_Clostridia_UCG-014_palette) +
  scale_x_discrete(limits = c("tylosin","control"), labels =c("TYLOSIN","NO TYLOSIN")) +
  scale_y_continuous(expand = c(0.001,0,5.3,0)) +
  geom_errorbar(lg_intestine_epithelium_Clostridia_UCG-014_family, mapping= aes(x= treatment, y= Abundance),
                stat = "summary", width = 0.5, linewidth = 0.75) +
  theme(legend.position = "none",
        plot.title = element_text(size = 44),
        panel.border = element_rect(colour = "black", size = 1.25),
        panel.grid.major.y = element_blank(),
        axis.ticks = element_line(size = 0.75, colour= "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 34, colour = "black"),
        axis.title.x = element_text(size =44),
        axis.text.x = element_text(size = 16, colour = "black"))

kruskal.test(lg_intestine_epithelium_Clostridia_UCG-014_family$Abundance, lg_intestine_epithelium_Clostridia_UCG-014_family$treatment) #NS

### PEPTOSTREPTOCOCCACEAE ####
## RUMEN LUMEN
rumen_lumen_Peptostreptococcaceae_genus <- subset_taxa(rumen_lumen_genus, Family=="Peptostreptococcaceae")
rumen_lumen_Peptostreptococcaceae_genus_filt <- merge_low_abundance(rumen_lumen_Peptostreptococcaceae_genus, threshold = 0.1) # 7 taxa
write.csv(tax_table(rumen_lumen_Peptostreptococcaceae_genus_filt), "rumen_lumen_Peptostreptococcaceae_taxa.csv")
write.csv(otu_table(rumen_lumen_Peptostreptococcaceae_genus_filt), "rumen_lumen_Peptostreptococcaceae_otus.csv")
rumen_lumen_Peptostreptococcaceae_genus_melt <- psmelt(rumen_lumen_Peptostreptococcaceae_genus_filt)

rumen_lumen_Peptostreptococcaceae_palette <- c("#9ED6DB")
rumen_lumen_Peptostreptococcaceae_family <- subset_taxa(rumen_lumen_family, Family=="Peptostreptococcaceae") %>%
  psmelt()

ggplot() + theme_bw() + labs(y= "Peptostreptococcaceae (RA%)", title = "LUMEN") +
  coord_flip() + theme_bw() +
  geom_bar(rumen_lumen_Peptostreptococcaceae_genus_melt, colour = "black", linewidth = 0.75,
           mapping = aes(x= treatment, y= Abundance, fill=Genus), stat = "summary") +
  #scale_fill_manual(values = rumen_lumen_Peptostreptococcaceae_palette) +
  scale_x_discrete(limits = c("tylosin","control"), labels =c("TYLOSIN","NO TYLOSIN")) +
  scale_y_continuous(expand = c(0.001,0,9.33,0)) +
  geom_errorbar(rumen_lumen_Peptostreptococcaceae_family, mapping= aes(x= treatment, y= Abundance),
                stat = "summary", width = 0.5, linewidth = 0.75) +
  theme(#legend.position = "none",
    plot.title = element_text(size = 44),
    panel.border = element_rect(colour = "black", size = 1.25),
    panel.grid.major.y = element_blank(),
    axis.ticks = element_line(size = 0.75, colour= "black"),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 34, colour = "black"),
    axis.title.x = element_text(size = 44),
    axis.text.x = element_text(size = 16, colour = "black"))

kruskal.test(rumen_lumen_Peptostreptococcaceae_family$Abundance, rumen_lumen_Peptostreptococcaceae_family$treatment) #NS

## RUMEN EPITHELIUM
rumen_epithelium_Peptostreptococcaceae_genus <- subset_taxa(rumen_epithelium_genus, Family=="Peptostreptococcaceae")
rumen_epithelium_Peptostreptococcaceae_genus_filt <- merge_low_abundance(rumen_epithelium_Peptostreptococcaceae_genus, threshold = 0.1) # 7 taxa
write.csv(tax_table(rumen_epithelium_Peptostreptococcaceae_genus_filt), "rumen_epithelium_Peptostreptococcaceae_taxa.csv")
write.csv(otu_table(rumen_epithelium_Peptostreptococcaceae_genus_filt), "rumen_epithelium_Peptostreptococcaceae_otus.csv")
rumen_epithelium_Peptostreptococcaceae_genus_melt <- psmelt(rumen_epithelium_Peptostreptococcaceae_genus_filt)

rumen_epithelium_Peptostreptococcaceae_palette <- c("grey55","#9ED6DB")
rumen_epithelium_Peptostreptococcaceae_family <- subset_taxa(rumen_epithelium_family, Family=="Peptostreptococcaceae") %>%
  psmelt()

ggplot() + theme_bw() + labs(y= "Peptostreptococcaceae (RA%)", title = "EPITHELIUM") +
  coord_flip() + theme_bw() +
  geom_bar(rumen_epithelium_Peptostreptococcaceae_genus_melt, colour = "black", linewidth = 0.75,
           mapping = aes(x= treatment, y= Abundance, fill=Genus), stat = "summary") +
  #scale_fill_manual(values = rumen_epithelium_Peptostreptococcaceae_palette) +
  scale_x_discrete(limits = c("tylosin","control"), labels =c("TYLOSIN","NO TYLOSIN")) +
  scale_y_continuous(expand = c(0.001,0,11,0)) +
  geom_errorbar(rumen_epithelium_Peptostreptococcaceae_family, mapping= aes(x= treatment, y= Abundance),
                stat = "summary", width = 0.5, linewidth = 0.75) +
  theme(legend.position = "none",
        plot.title = element_text(size = 44),
        panel.border = element_rect(colour = "black", size = 1.25),
        panel.grid.major.y = element_blank(),
        axis.ticks = element_line(size = 0.75, colour= "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 34, colour = "black"),
        axis.title.x = element_text(size =44),
        axis.text.x = element_text(size = 16, colour = "black"))

kruskal.test(rumen_epithelium_Peptostreptococcaceae_family$Abundance, rumen_epithelium_Peptostreptococcaceae_family$treatment) #NS

### SM INTESTINE LUMEN
sm_intestine_lumen_Peptostreptococcaceae_genus <- subset_taxa(sm_intestine_lumen_genus, Family=="Peptostreptococcaceae")
sm_intestine_lumen_Peptostreptococcaceae_genus_filt <- merge_low_abundance(sm_intestine_lumen_Peptostreptococcaceae_genus, threshold = 0.1) # 7 taxa
write.csv(tax_table(sm_intestine_lumen_Peptostreptococcaceae_genus_filt), "sm_intestine_lumen_Peptostreptococcaceae_taxa.csv")
write.csv(otu_table(sm_intestine_lumen_Peptostreptococcaceae_genus_filt), "sm_intestine_lumen_Peptostreptococcaceae_otus.csv")
sm_intestine_lumen_Peptostreptococcaceae_genus_melt <- psmelt(sm_intestine_lumen_Peptostreptococcaceae_genus_filt)

sm_intestine_lumen_Peptostreptococcaceae_palette <- c("#e754ae","#F3B2E3","#ffe8f9","grey88")
sm_intestine_lumen_Peptostreptococcaceae_family <- subset_taxa(sm_intestine_lumen_family, Family=="Peptostreptococcaceae") %>%
  psmelt()

ggplot() + theme_bw() + labs(y= "Peptostreptococcaceae (RA%)", title = "LUMEN") +
  coord_flip() + theme_bw() +
  geom_bar(sm_intestine_lumen_Peptostreptococcaceae_genus_melt, colour = "black", linewidth = 0.75,
           mapping = aes(x= treatment, y= Abundance, fill=Genus), stat = "summary") +
  scale_fill_manual(values = sm_intestine_lumen_Peptostreptococcaceae_palette) +
  scale_x_discrete(limits = c("tylosin","control"), labels =c("TYLOSIN","NO TYLOSIN")) +
  scale_y_continuous(expand = c(0.001,0,0.84,0)) +
  geom_errorbar(sm_intestine_lumen_Peptostreptococcaceae_family, mapping= aes(x= treatment, y= Abundance),
                stat = "summary", width = 0.5, linewidth = 0.75) +
  theme(legend.position = "none",
        plot.title = element_text(size = 44),
        panel.border = element_rect(colour = "black", size = 1.25),
        panel.grid.major.y = element_blank(),
        axis.ticks = element_line(size = 0.75, colour= "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 34, colour = "black"),
        axis.title.x = element_text(size =44),
        axis.text.x = element_text(size = 16, colour = "black"))

kruskal.test(sm_intestine_lumen_Peptostreptococcaceae_family$Abundance, sm_intestine_lumen_Peptostreptococcaceae_family$treatment) #NS

## sm_intestine EPITHELIUM
sm_intestine_epithelium_Peptostreptococcaceae_genus <- subset_taxa(sm_intestine_epithelium_genus, Family=="Peptostreptococcaceae")
sm_intestine_epithelium_Peptostreptococcaceae_genus_filt <- merge_low_abundance(sm_intestine_epithelium_Peptostreptococcaceae_genus, threshold = 0.1) # 7 taxa
write.csv(tax_table(sm_intestine_epithelium_Peptostreptococcaceae_genus_filt), "sm_intestine_epithelium_Peptostreptococcaceae_taxa.csv")
write.csv(otu_table(sm_intestine_epithelium_Peptostreptococcaceae_genus_filt), "sm_intestine_epithelium_Peptostreptococcaceae_otus.csv")
sm_intestine_epithelium_Peptostreptococcaceae_genus_melt <- psmelt(sm_intestine_epithelium_Peptostreptococcaceae_genus_filt)

sm_intestine_epithelium_Peptostreptococcaceae_palette <- c("#e754ae","#F3B2E3","#ffe8f9","grey88")
sm_intestine_epithelium_Peptostreptococcaceae_family <- subset_taxa(sm_intestine_epithelium_family, Family=="Peptostreptococcaceae") %>%
  psmelt()

ggplot() + theme_bw() + labs(y= "Peptostreptococcaceae (RA%)", title = "EPITHELIUM") +
  coord_flip() + theme_bw() +
  geom_bar(sm_intestine_epithelium_Peptostreptococcaceae_genus_melt, colour = "black", linewidth = 0.75,
           mapping = aes(x= treatment, y= Abundance, fill=Genus), stat = "summary") +
  scale_fill_manual(values = sm_intestine_epithelium_Peptostreptococcaceae_palette) +
  scale_x_discrete(limits = c("tylosin","control"), labels =c("TYLOSIN","NO TYLOSIN")) +
  scale_y_continuous(expand = c(0.001,0,0.905,0)) +
  geom_errorbar(sm_intestine_epithelium_Peptostreptococcaceae_family, mapping= aes(x= treatment, y= Abundance),
                stat = "summary", width = 0.5, linewidth = 0.75) +
  theme(legend.position = "none",
    plot.title = element_text(size = 44),
    panel.border = element_rect(colour = "black", size = 1.25),
    panel.grid.major.y = element_blank(),
    axis.ticks = element_line(size = 0.75, colour= "black"),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 34, colour = "black"),
    axis.title.x = element_text(size =44),
    axis.text.x = element_text(size = 16, colour = "black"))

kruskal.test(sm_intestine_epithelium_Peptostreptococcaceae_family$Abundance, sm_intestine_epithelium_Peptostreptococcaceae_family$treatment) #NS

### LG INT LUMEN
lg_intestine_lumen_Peptostreptococcaceae_genus <- subset_taxa(lg_intestine_lumen_genus, Family=="Peptostreptococcaceae")
lg_intestine_lumen_Peptostreptococcaceae_genus_filt <- merge_low_abundance(lg_intestine_lumen_Peptostreptococcaceae_genus, threshold = 0.1) # 7 taxa
write.csv(tax_table(lg_intestine_lumen_Peptostreptococcaceae_genus_filt), "lg_intestine_lumen_Peptostreptococcaceae_taxa.csv")
write.csv(otu_table(lg_intestine_lumen_Peptostreptococcaceae_genus_filt), "lg_intestine_lumen_Peptostreptococcaceae_otus.csv")
lg_intestine_lumen_Peptostreptococcaceae_genus_melt <- psmelt(lg_intestine_lumen_Peptostreptococcaceae_genus_filt)

lg_intestine_lumen_Peptostreptococcaceae_palette <- c("#e754ae","#F3B2E3","#ffe8f9","grey88")
lg_intestine_lumen_Peptostreptococcaceae_family <- subset_taxa(lg_intestine_lumen_family, Family=="Peptostreptococcaceae") %>%
  psmelt()

ggplot() + theme_bw() + labs(y= "Peptostreptococcaceae (RA%)", title = "LUMEN") +
  coord_flip() + theme_bw() +
  geom_bar(lg_intestine_lumen_Peptostreptococcaceae_genus_melt, colour = "black", linewidth = 0.75,
           mapping = aes(x= treatment, y= Abundance, fill=Genus), stat = "summary") +
  scale_fill_manual(values = lg_intestine_lumen_Peptostreptococcaceae_palette) +
  scale_x_discrete(limits = c("tylosin","control"), labels =c("TYLOSIN","NO TYLOSIN")) +
  scale_y_continuous(expand = c(0.001,0,4,0)) +
  geom_errorbar(lg_intestine_lumen_Peptostreptococcaceae_family, mapping= aes(x= treatment, y= Abundance),
                stat = "summary", width = 0.5, linewidth = 0.75) +
  theme(#legend.position = "none",
        plot.title = element_text(size = 44),
        panel.border = element_rect(colour = "black", size = 1.25),
        panel.grid.major.y = element_blank(),
        axis.ticks = element_line(size = 0.75, colour= "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 34, colour = "black"),
        axis.title.x = element_text(size =44),
        axis.text.x = element_text(size = 16, colour = "black"))

kruskal.test(lg_intestine_lumen_Peptostreptococcaceae_family$Abundance, lg_intestine_lumen_Peptostreptococcaceae_family$treatment) #NS

## lg_intestine EPITHELIUM
lg_intestine_epithelium_Peptostreptococcaceae_genus <- subset_taxa(lg_intestine_epithelium_genus, Family=="Peptostreptococcaceae")
lg_intestine_epithelium_Peptostreptococcaceae_genus_filt <- merge_low_abundance(lg_intestine_epithelium_Peptostreptococcaceae_genus, threshold = 0.1) # 7 taxa
write.csv(tax_table(lg_intestine_epithelium_Peptostreptococcaceae_genus_filt), "lg_intestine_epithelium_Peptostreptococcaceae_taxa.csv")
write.csv(otu_table(lg_intestine_epithelium_Peptostreptococcaceae_genus_filt), "lg_intestine_epithelium_Peptostreptococcaceae_otus.csv")
lg_intestine_epithelium_Peptostreptococcaceae_genus_melt <- psmelt(lg_intestine_epithelium_Peptostreptococcaceae_genus_filt)

lg_intestine_epithelium_Peptostreptococcaceae_palette <- c("#e754ae","#F3B2E3","#ffe8f9","grey88")
lg_intestine_epithelium_Peptostreptococcaceae_family <- subset_taxa(lg_intestine_epithelium_family, Family=="Peptostreptococcaceae") %>%
  psmelt()

ggplot() + theme_bw() + labs(y= "Peptostreptococcaceae (RA%)", title = "EPITHELIUM") +
  coord_flip() + theme_bw() +
  geom_bar(lg_intestine_epithelium_Peptostreptococcaceae_genus_melt, colour = "black", linewidth = 0.75,
           mapping = aes(x= treatment, y= Abundance, fill=Genus), stat = "summary") +
  scale_fill_manual(values = lg_intestine_epithelium_Peptostreptococcaceae_palette) +
  scale_x_discrete(limits = c("tylosin","control"), labels =c("TYLOSIN","NO TYLOSIN")) +
  scale_y_continuous(expand = c(0.001,0,5.3,0)) +
  geom_errorbar(lg_intestine_epithelium_Peptostreptococcaceae_family, mapping= aes(x= treatment, y= Abundance),
                stat = "summary", width = 0.5, linewidth = 0.75) +
  theme(#legend.position = "none",
        plot.title = element_text(size = 44),
        panel.border = element_rect(colour = "black", size = 1.25),
        panel.grid.major.y = element_blank(),
        axis.ticks = element_line(size = 0.75, colour= "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 34, colour = "black"),
        axis.title.x = element_text(size =44),
        axis.text.x = element_text(size = 16, colour = "black"))

kruskal.test(lg_intestine_epithelium_Peptostreptococcaceae_family$Abundance, lg_intestine_epithelium_Peptostreptococcaceae_family$treatment) #NS

### CLOSTRIDIACEAE ####
## RUMEN LUMEN
rumen_lumen_Clostridiaceae_genus <- subset_taxa(rumen_lumen_genus, Family=="Clostridiaceae")
rumen_lumen_Clostridiaceae_genus_filt <- merge_low_abundance(rumen_lumen_Clostridiaceae_genus, threshold = 0.1) # 7 taxa
write.csv(tax_table(rumen_lumen_Clostridiaceae_genus_filt), "rumen_lumen_Clostridiaceae_taxa.csv")
write.csv(otu_table(rumen_lumen_Clostridiaceae_genus_filt), "rumen_lumen_Clostridiaceae_otus.csv")
rumen_lumen_Clostridiaceae_genus_melt <- psmelt(rumen_lumen_Clostridiaceae_genus_filt)

rumen_lumen_Clostridiaceae_palette <- c("#9ED6DB")
rumen_lumen_Clostridiaceae_family <- subset_taxa(rumen_lumen_family, Family=="Clostridiaceae") %>%
  psmelt()

ggplot() + theme_bw() + labs(y= "Clostridiaceae (RA%)", title = "LUMEN") +
  coord_flip() + theme_bw() +
  geom_bar(rumen_lumen_Clostridiaceae_genus_melt, colour = "black", linewidth = 0.75,
           mapping = aes(x= treatment, y= Abundance, fill=Genus), stat = "summary") +
  scale_fill_manual(values = rumen_lumen_Clostridiaceae_palette) +
  scale_x_discrete(limits = c("tylosin","control"), labels =c("TYLOSIN","NO TYLOSIN")) +
  scale_y_continuous(expand = c(0.001,0,9.33,0)) +
  geom_errorbar(rumen_lumen_Clostridiaceae_family, mapping= aes(x= treatment, y= Abundance),
                stat = "summary", width = 0.5, linewidth = 0.75) +
  theme(legend.position = "none",
        plot.title = element_text(size = 44),
        panel.border = element_rect(colour = "black", size = 1.25),
        panel.grid.major.y = element_blank(),
        axis.ticks = element_line(size = 0.75, colour= "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 34, colour = "black"),
        axis.title.x = element_text(size = 44),
        axis.text.x = element_text(size = 16, colour = "black"))

kruskal.test(rumen_lumen_Clostridiaceae_family$Abundance, rumen_lumen_Clostridiaceae_family$treatment) #NS

## RUMEN EPITHELIUM
rumen_epithelium_Clostridiaceae_genus <- subset_taxa(rumen_epithelium_genus, Family=="Clostridiaceae")
rumen_epithelium_Clostridiaceae_genus_filt <- merge_low_abundance(rumen_epithelium_Clostridiaceae_genus, threshold = 0.1) # 7 taxa
write.csv(tax_table(rumen_epithelium_Clostridiaceae_genus_filt), "rumen_epithelium_Clostridiaceae_taxa.csv")
write.csv(otu_table(rumen_epithelium_Clostridiaceae_genus_filt), "rumen_epithelium_Clostridiaceae_otus.csv")
rumen_epithelium_Clostridiaceae_genus_melt <- psmelt(rumen_epithelium_Clostridiaceae_genus_filt)

rumen_epithelium_Clostridiaceae_palette <- c("grey55","#9ED6DB")
rumen_epithelium_Clostridiaceae_family <- subset_taxa(rumen_epithelium_family, Family=="Clostridiaceae") %>%
  psmelt()

ggplot() + theme_bw() + labs(y= "Clostridiaceae (RA%)", title = "EPITHELIUM") +
  coord_flip() + theme_bw() +
  geom_bar(rumen_epithelium_Clostridiaceae_genus_melt, colour = "black", linewidth = 0.75,
           mapping = aes(x= treatment, y= Abundance, fill=Genus), stat = "summary") +
  scale_fill_manual(values = rumen_epithelium_Clostridiaceae_palette) +
  scale_x_discrete(limits = c("tylosin","control"), labels =c("TYLOSIN","NO TYLOSIN")) +
  scale_y_continuous(expand = c(0.001,0,11,0)) +
  geom_errorbar(rumen_epithelium_Clostridiaceae_family, mapping= aes(x= treatment, y= Abundance),
                stat = "summary", width = 0.5, linewidth = 0.75) +
  theme(legend.position = "none",
        plot.title = element_text(size = 44),
        panel.border = element_rect(colour = "black", size = 1.25),
        panel.grid.major.y = element_blank(),
        axis.ticks = element_line(size = 0.75, colour= "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 34, colour = "black"),
        axis.title.x = element_text(size =44),
        axis.text.x = element_text(size = 16, colour = "black"))

kruskal.test(rumen_epithelium_Clostridiaceae_family$Abundance, rumen_epithelium_Clostridiaceae_family$treatment) #NS

### SM INTESTINE LUMEN
sm_intestine_lumen_Clostridiaceae_genus <- subset_taxa(sm_intestine_lumen_genus, Family=="Clostridiaceae")
sm_intestine_lumen_Clostridiaceae_genus_filt <- merge_low_abundance(sm_intestine_lumen_Clostridiaceae_genus, threshold = 0.1) # 7 taxa
write.csv(tax_table(sm_intestine_lumen_Clostridiaceae_genus_filt), "sm_intestine_lumen_Clostridiaceae_taxa.csv")
write.csv(otu_table(sm_intestine_lumen_Clostridiaceae_genus_filt), "sm_intestine_lumen_Clostridiaceae_otus.csv")
sm_intestine_lumen_Clostridiaceae_genus_melt <- psmelt(sm_intestine_lumen_Clostridiaceae_genus_filt)

sm_intestine_lumen_Clostridiaceae_palette <- c("red4","indianred1")
sm_intestine_lumen_Clostridiaceae_family <- subset_taxa(sm_intestine_lumen_family, Family=="Clostridiaceae") %>%
  psmelt()

ggplot() + theme_bw() + labs(y= "Clostridiaceae (RA%)", title = "LUMEN") +
  coord_flip() + theme_bw() +
  geom_bar(sm_intestine_lumen_Clostridiaceae_genus_melt, colour = "black", linewidth = 0.75,
           mapping = aes(x= treatment, y= Abundance, fill=Genus), stat = "summary") +
  scale_fill_manual(values = sm_intestine_lumen_Clostridiaceae_palette) +
  scale_x_discrete(limits = c("tylosin","control"), labels =c("TYLOSIN","NO TYLOSIN")) +
  scale_y_continuous(expand = c(0.001,0,4.6,0)) +
  geom_errorbar(sm_intestine_lumen_Clostridiaceae_family, mapping= aes(x= treatment, y= Abundance),
                stat = "summary", width = 0.5, linewidth = 0.75) +
  theme(#legend.position = "none",
        plot.title = element_text(size = 44),
        panel.border = element_rect(colour = "black", size = 1.25),
        panel.grid.major.y = element_blank(),
        axis.ticks = element_line(size = 0.75, colour= "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 34, colour = "black"),
        axis.title.x = element_text(size =44),
        axis.text.x = element_text(size = 16, colour = "black"))

kruskal.test(sm_intestine_lumen_Clostridiaceae_family$Abundance, sm_intestine_lumen_Clostridiaceae_family$treatment) #NS

## sm_intestine EPITHELIUM
sm_intestine_epithelium_Clostridiaceae_genus <- subset_taxa(sm_intestine_epithelium_genus, Family=="Clostridiaceae")
sm_intestine_epithelium_Clostridiaceae_genus_filt <- merge_low_abundance(sm_intestine_epithelium_Clostridiaceae_genus, threshold = 0.1) # 7 taxa
write.csv(tax_table(sm_intestine_epithelium_Clostridiaceae_genus_filt), "sm_intestine_epithelium_Clostridiaceae_taxa.csv")
write.csv(otu_table(sm_intestine_epithelium_Clostridiaceae_genus_filt), "sm_intestine_epithelium_Clostridiaceae_otus.csv")
sm_intestine_epithelium_Clostridiaceae_genus_melt <- psmelt(sm_intestine_epithelium_Clostridiaceae_genus_filt)

sm_intestine_epithelium_Clostridiaceae_palette <- c("red4","indianred1")
sm_intestine_epithelium_Clostridiaceae_family <- subset_taxa(sm_intestine_epithelium_family, Family=="Clostridiaceae") %>%
  psmelt()

ggplot() + theme_bw() + labs(y= "Clostridiaceae (RA%)", title = "EPITHELIUM") +
  coord_flip() + theme_bw() +
  geom_bar(sm_intestine_epithelium_Clostridiaceae_genus_melt, colour = "black", linewidth = 0.75,
           mapping = aes(x= treatment, y= Abundance, fill=Genus), stat = "summary") +
  scale_fill_manual(values = sm_intestine_epithelium_Clostridiaceae_palette) +
  scale_x_discrete(limits = c("tylosin","control"), labels =c("TYLOSIN","NO TYLOSIN")) +
  scale_y_continuous(expand = c(0.001,0,4.5,0)) +
  geom_errorbar(sm_intestine_epithelium_Clostridiaceae_family, mapping= aes(x= treatment, y= Abundance),
                stat = "summary", width = 0.5, linewidth = 0.75) +
  theme(#legend.position = "none",
    plot.title = element_text(size = 44),
    panel.border = element_rect(colour = "black", size = 1.25),
    panel.grid.major.y = element_blank(),
    axis.ticks = element_line(size = 0.75, colour= "black"),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 34, colour = "black"),
    axis.title.x = element_text(size =44),
    axis.text.x = element_text(size = 16, colour = "black"))

kruskal.test(sm_intestine_epithelium_Clostridiaceae_family$Abundance, sm_intestine_epithelium_Clostridiaceae_family$treatment) #NS

### LG INT LUMEN
lg_intestine_lumen_Clostridiaceae_genus <- subset_taxa(lg_intestine_lumen_genus, Family=="Clostridiaceae")
lg_intestine_lumen_Clostridiaceae_genus_filt <- merge_low_abundance(lg_intestine_lumen_Clostridiaceae_genus, threshold = 0.1) # 7 taxa
write.csv(tax_table(lg_intestine_lumen_Clostridiaceae_genus_filt), "lg_intestine_lumen_Clostridiaceae_taxa.csv")
write.csv(otu_table(lg_intestine_lumen_Clostridiaceae_genus_filt), "lg_intestine_lumen_Clostridiaceae_otus.csv")
lg_intestine_lumen_Clostridiaceae_genus_melt <- psmelt(lg_intestine_lumen_Clostridiaceae_genus_filt)

lg_intestine_lumen_Clostridiaceae_palette <- c("#6030A8","grey88")
lg_intestine_lumen_Clostridiaceae_family <- subset_taxa(lg_intestine_lumen_family, Family=="Clostridiaceae") %>%
  psmelt()

ggplot() + theme_bw() + labs(y= "Clostridiaceae (RA%)", title = "LUMEN") +
  coord_flip() + theme_bw() +
  geom_bar(lg_intestine_lumen_Clostridiaceae_genus_melt, colour = "black", linewidth = 0.75,
           mapping = aes(x= treatment, y= Abundance, fill=Genus), stat = "summary") +
  scale_fill_manual(values = lg_intestine_lumen_Clostridiaceae_palette) +
  scale_x_discrete(limits = c("tylosin","control"), labels =c("TYLOSIN","NO TYLOSIN")) +
  scale_y_continuous(expand = c(0.001,0,4,0)) +
  geom_errorbar(lg_intestine_lumen_Clostridiaceae_family, mapping= aes(x= treatment, y= Abundance),
                stat = "summary", width = 0.5, linewidth = 0.75) +
  theme(legend.position = "none",
        plot.title = element_text(size = 44),
        panel.border = element_rect(colour = "black", size = 1.25),
        panel.grid.major.y = element_blank(),
        axis.ticks = element_line(size = 0.75, colour= "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 34, colour = "black"),
        axis.title.x = element_text(size =44),
        axis.text.x = element_text(size = 16, colour = "black"))

kruskal.test(lg_intestine_lumen_Clostridiaceae_family$Abundance, lg_intestine_lumen_Clostridiaceae_family$treatment) #NS

## lg_intestine EPITHELIUM
lg_intestine_epithelium_Clostridiaceae_genus <- subset_taxa(lg_intestine_epithelium_genus, Family=="Clostridiaceae")
lg_intestine_epithelium_Clostridiaceae_genus_filt <- merge_low_abundance(lg_intestine_epithelium_Clostridiaceae_genus, threshold = 0.1) # 7 taxa
write.csv(tax_table(lg_intestine_epithelium_Clostridiaceae_genus_filt), "lg_intestine_epithelium_Clostridiaceae_taxa.csv")
write.csv(otu_table(lg_intestine_epithelium_Clostridiaceae_genus_filt), "lg_intestine_epithelium_Clostridiaceae_otus.csv")
lg_intestine_epithelium_Clostridiaceae_genus_melt <- psmelt(lg_intestine_epithelium_Clostridiaceae_genus_filt)

lg_intestine_epithelium_Clostridiaceae_palette <- c("#4B0082","#6030A8","grey88")
lg_intestine_epithelium_Clostridiaceae_family <- subset_taxa(lg_intestine_epithelium_family, Family=="Clostridiaceae") %>%
  psmelt()

ggplot() + theme_bw() + labs(y= "Clostridiaceae (RA%)", title = "EPITHELIUM") +
  coord_flip() + theme_bw() +
  geom_bar(lg_intestine_epithelium_Clostridiaceae_genus_melt, colour = "black", linewidth = 0.75,
           mapping = aes(x= treatment, y= Abundance, fill=Genus), stat = "summary") +
  scale_fill_manual(values = lg_intestine_epithelium_Clostridiaceae_palette) +
  scale_x_discrete(limits = c("tylosin","control"), labels =c("TYLOSIN","NO TYLOSIN")) +
  scale_y_continuous(expand = c(0.001,0,5.3,0)) +
  geom_errorbar(lg_intestine_epithelium_Clostridiaceae_family, mapping= aes(x= treatment, y= Abundance),
                stat = "summary", width = 0.5, linewidth = 0.75) +
  theme(legend.position = "none",
        plot.title = element_text(size = 44),
        panel.border = element_rect(colour = "black", size = 1.25),
        panel.grid.major.y = element_blank(),
        axis.ticks = element_line(size = 0.75, colour= "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 34, colour = "black"),
        axis.title.x = element_text(size =44),
        axis.text.x = element_text(size = 16, colour = "black"))

kruskal.test(lg_intestine_epithelium_Clostridiaceae_family$Abundance, lg_intestine_epithelium_Clostridiaceae_family$treatment) #NS

### Erysipelotrichaceae ####
## RUMEN LUMEN
rumen_lumen_Erysipelotrichaceae_genus <- subset_taxa(rumen_lumen_genus, Family=="Erysipelotrichaceae")
rumen_lumen_Erysipelotrichaceae_genus_filt <- merge_low_abundance(rumen_lumen_Erysipelotrichaceae_genus, threshold = 0.1) # 7 taxa
write.csv(tax_table(rumen_lumen_Erysipelotrichaceae_genus_filt), "rumen_lumen_Erysipelotrichaceae_taxa.csv")
write.csv(otu_table(rumen_lumen_Erysipelotrichaceae_genus_filt), "rumen_lumen_Erysipelotrichaceae_otus.csv")
rumen_lumen_Erysipelotrichaceae_genus_melt <- psmelt(rumen_lumen_Erysipelotrichaceae_genus_filt)

rumen_lumen_Erysipelotrichaceae_palette <- c("#9ED6DB")
rumen_lumen_Erysipelotrichaceae_family <- subset_taxa(rumen_lumen_family, Family=="Erysipelotrichaceae") %>%
  psmelt()

ggplot() + theme_bw() + labs(y= "Erysipelotrichaceae (RA%)", title = "LUMEN") +
  coord_flip() + theme_bw() +
  geom_bar(rumen_lumen_Erysipelotrichaceae_genus_melt, colour = "black", linewidth = 0.75,
           mapping = aes(x= treatment, y= Abundance, fill=Genus), stat = "summary") +
  #scale_fill_manual(values = rumen_lumen_Erysipelotrichaceae_palette) +
  scale_x_discrete(limits = c("tylosin","control"), labels =c("TYLOSIN","NO TYLOSIN")) +
  scale_y_continuous(expand = c(0.001,0,9.33,0)) +
  geom_errorbar(rumen_lumen_Erysipelotrichaceae_family, mapping= aes(x= treatment, y= Abundance),
                stat = "summary", width = 0.5, linewidth = 0.75) +
  theme(#legend.position = "none",
        plot.title = element_text(size = 44),
        panel.border = element_rect(colour = "black", size = 1.25),
        panel.grid.major.y = element_blank(),
        axis.ticks = element_line(size = 0.75, colour= "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 34, colour = "black"),
        axis.title.x = element_text(size = 44),
        axis.text.x = element_text(size = 16, colour = "black"))

kruskal.test(rumen_lumen_Erysipelotrichaceae_family$Abundance, rumen_lumen_Erysipelotrichaceae_family$treatment) #NS

## RUMEN EPITHELIUM
rumen_epithelium_Erysipelotrichaceae_genus <- subset_taxa(rumen_epithelium_genus, Family=="Erysipelotrichaceae")
rumen_epithelium_Erysipelotrichaceae_genus_filt <- merge_low_abundance(rumen_epithelium_Erysipelotrichaceae_genus, threshold = 0.1) # 7 taxa
write.csv(tax_table(rumen_epithelium_Erysipelotrichaceae_genus_filt), "rumen_epithelium_Erysipelotrichaceae_taxa.csv")
write.csv(otu_table(rumen_epithelium_Erysipelotrichaceae_genus_filt), "rumen_epithelium_Erysipelotrichaceae_otus.csv")
rumen_epithelium_Erysipelotrichaceae_genus_melt <- psmelt(rumen_epithelium_Erysipelotrichaceae_genus_filt)

rumen_epithelium_Erysipelotrichaceae_palette <- c("grey55","#9ED6DB")
rumen_epithelium_Erysipelotrichaceae_family <- subset_taxa(rumen_epithelium_family, Family=="Erysipelotrichaceae") %>%
  psmelt()

ggplot() + theme_bw() + labs(y= "Erysipelotrichaceae (RA%)", title = "EPITHELIUM") +
  coord_flip() + theme_bw() +
  geom_bar(rumen_epithelium_Erysipelotrichaceae_genus_melt, colour = "black", linewidth = 0.75,
           mapping = aes(x= treatment, y= Abundance, fill=Genus), stat = "summary") +
  #scale_fill_manual(values = rumen_epithelium_Erysipelotrichaceae_palette) +
  scale_x_discrete(limits = c("tylosin","control"), labels =c("TYLOSIN","NO TYLOSIN")) +
  scale_y_continuous(expand = c(0.001,0,11,0)) +
  geom_errorbar(rumen_epithelium_Erysipelotrichaceae_family, mapping= aes(x= treatment, y= Abundance),
                stat = "summary", width = 0.5, linewidth = 0.75) +
  theme(#legend.position = "none",
        plot.title = element_text(size = 44),
        panel.border = element_rect(colour = "black", size = 1.25),
        panel.grid.major.y = element_blank(),
        axis.ticks = element_line(size = 0.75, colour= "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 34, colour = "black"),
        axis.title.x = element_text(size =44),
        axis.text.x = element_text(size = 16, colour = "black"))

kruskal.test(rumen_epithelium_Erysipelotrichaceae_family$Abundance, rumen_epithelium_Erysipelotrichaceae_family$treatment) #NS

### SM INTESTINE LUMEN
sm_intestine_lumen_Erysipelotrichaceae_genus <- subset_taxa(sm_intestine_lumen_genus, Family=="Erysipelotrichaceae")
sm_intestine_lumen_Erysipelotrichaceae_genus_filt <- merge_low_abundance(sm_intestine_lumen_Erysipelotrichaceae_genus, threshold = 0.1) # 7 taxa
write.csv(tax_table(sm_intestine_lumen_Erysipelotrichaceae_genus_filt), "sm_intestine_lumen_Erysipelotrichaceae_taxa.csv")
write.csv(otu_table(sm_intestine_lumen_Erysipelotrichaceae_genus_filt), "sm_intestine_lumen_Erysipelotrichaceae_otus.csv")
sm_intestine_lumen_Erysipelotrichaceae_genus_melt <- psmelt(sm_intestine_lumen_Erysipelotrichaceae_genus_filt)

sm_intestine_lumen_Erysipelotrichaceae_palette <- c("red4","indianred1")
sm_intestine_lumen_Erysipelotrichaceae_family <- subset_taxa(sm_intestine_lumen_family, Family=="Erysipelotrichaceae") %>%
  psmelt()

ggplot() + theme_bw() + labs(y= "Erysipelotrichaceae (RA%)", title = "LUMEN") +
  coord_flip() + theme_bw() +
  geom_bar(sm_intestine_lumen_Erysipelotrichaceae_genus_melt, colour = "black", linewidth = 0.75,
           mapping = aes(x= treatment, y= Abundance, fill=Genus), stat = "summary") +
  #scale_fill_manual(values = sm_intestine_lumen_Erysipelotrichaceae_palette) +
  scale_x_discrete(limits = c("tylosin","control"), labels =c("TYLOSIN","NO TYLOSIN")) +
  scale_y_continuous(expand = c(0.001,0,4.6,0)) +
  geom_errorbar(sm_intestine_lumen_Erysipelotrichaceae_family, mapping= aes(x= treatment, y= Abundance),
                stat = "summary", width = 0.5, linewidth = 0.75) +
  theme(#legend.position = "none",
        plot.title = element_text(size = 44),
        panel.border = element_rect(colour = "black", size = 1.25),
        panel.grid.major.y = element_blank(),
        axis.ticks = element_line(size = 0.75, colour= "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 34, colour = "black"),
        axis.title.x = element_text(size =44),
        axis.text.x = element_text(size = 16, colour = "black"))

kruskal.test(sm_intestine_lumen_Erysipelotrichaceae_family$Abundance, sm_intestine_lumen_Erysipelotrichaceae_family$treatment) #NS

## sm_intestine EPITHELIUM
sm_intestine_epithelium_Erysipelotrichaceae_genus <- subset_taxa(sm_intestine_epithelium_genus, Family=="Erysipelotrichaceae")
sm_intestine_epithelium_Erysipelotrichaceae_genus_filt <- merge_low_abundance(sm_intestine_epithelium_Erysipelotrichaceae_genus, threshold = 0.1) # 7 taxa
write.csv(tax_table(sm_intestine_epithelium_Erysipelotrichaceae_genus_filt), "sm_intestine_epithelium_Erysipelotrichaceae_taxa.csv")
write.csv(otu_table(sm_intestine_epithelium_Erysipelotrichaceae_genus_filt), "sm_intestine_epithelium_Erysipelotrichaceae_otus.csv")
sm_intestine_epithelium_Erysipelotrichaceae_genus_melt <- psmelt(sm_intestine_epithelium_Erysipelotrichaceae_genus_filt)

sm_intestine_epithelium_Erysipelotrichaceae_palette <- c("red4","indianred1")
sm_intestine_epithelium_Erysipelotrichaceae_family <- subset_taxa(sm_intestine_epithelium_family, Family=="Erysipelotrichaceae") %>%
  psmelt()

ggplot() + theme_bw() + labs(y= "Erysipelotrichaceae (RA%)", title = "EPITHELIUM") +
  coord_flip() + theme_bw() +
  geom_bar(sm_intestine_epithelium_Erysipelotrichaceae_genus_melt, colour = "black", linewidth = 0.75,
           mapping = aes(x= treatment, y= Abundance, fill=Genus), stat = "summary") +
  #scale_fill_manual(values = sm_intestine_epithelium_Erysipelotrichaceae_palette) +
  scale_x_discrete(limits = c("tylosin","control"), labels =c("TYLOSIN","NO TYLOSIN")) +
  scale_y_continuous(expand = c(0.001,0,4.5,0)) +
  geom_errorbar(sm_intestine_epithelium_Erysipelotrichaceae_family, mapping= aes(x= treatment, y= Abundance),
                stat = "summary", width = 0.5, linewidth = 0.75) +
  theme(#legend.position = "none",
    plot.title = element_text(size = 44),
    panel.border = element_rect(colour = "black", size = 1.25),
    panel.grid.major.y = element_blank(),
    axis.ticks = element_line(size = 0.75, colour= "black"),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 34, colour = "black"),
    axis.title.x = element_text(size =44),
    axis.text.x = element_text(size = 16, colour = "black"))

kruskal.test(sm_intestine_epithelium_Erysipelotrichaceae_family$Abundance, sm_intestine_epithelium_Erysipelotrichaceae_family$treatment) #NS

### LG INT LUMEN
lg_intestine_lumen_Erysipelotrichaceae_genus <- subset_taxa(lg_intestine_lumen_genus, Family=="Erysipelotrichaceae")
lg_intestine_lumen_Erysipelotrichaceae_genus_filt <- merge_low_abundance(lg_intestine_lumen_Erysipelotrichaceae_genus, threshold = 0.1) # 7 taxa
write.csv(tax_table(lg_intestine_lumen_Erysipelotrichaceae_genus_filt), "lg_intestine_lumen_Erysipelotrichaceae_taxa.csv")
write.csv(otu_table(lg_intestine_lumen_Erysipelotrichaceae_genus_filt), "lg_intestine_lumen_Erysipelotrichaceae_otus.csv")
lg_intestine_lumen_Erysipelotrichaceae_genus_melt <- psmelt(lg_intestine_lumen_Erysipelotrichaceae_genus_filt)

lg_intestine_lumen_Erysipelotrichaceae_palette <- c("#6030A8","grey88")
lg_intestine_lumen_Erysipelotrichaceae_family <- subset_taxa(lg_intestine_lumen_family, Family=="Erysipelotrichaceae") %>%
  psmelt()

ggplot() + theme_bw() + labs(y= "Erysipelotrichaceae (RA%)", title = "LUMEN") +
  coord_flip() + theme_bw() +
  geom_bar(lg_intestine_lumen_Erysipelotrichaceae_genus_melt, colour = "black", linewidth = 0.75,
           mapping = aes(x= treatment, y= Abundance, fill=Genus), stat = "summary") +
  #scale_fill_manual(values = lg_intestine_lumen_Erysipelotrichaceae_palette) +
  scale_x_discrete(limits = c("tylosin","control"), labels =c("TYLOSIN","NO TYLOSIN")) +
  scale_y_continuous(expand = c(0.001,0,4,0)) +
  geom_errorbar(lg_intestine_lumen_Erysipelotrichaceae_family, mapping= aes(x= treatment, y= Abundance),
                stat = "summary", width = 0.5, linewidth = 0.75) +
  theme(#legend.position = "none",
        plot.title = element_text(size = 44),
        panel.border = element_rect(colour = "black", size = 1.25),
        panel.grid.major.y = element_blank(),
        axis.ticks = element_line(size = 0.75, colour= "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 34, colour = "black"),
        axis.title.x = element_text(size =44),
        axis.text.x = element_text(size = 16, colour = "black"))

kruskal.test(lg_intestine_lumen_Erysipelotrichaceae_family$Abundance, lg_intestine_lumen_Erysipelotrichaceae_family$treatment) #NS

## lg_intestine EPITHELIUM
lg_intestine_epithelium_Erysipelotrichaceae_genus <- subset_taxa(lg_intestine_epithelium_genus, Family=="Erysipelotrichaceae")
lg_intestine_epithelium_Erysipelotrichaceae_genus_filt <- merge_low_abundance(lg_intestine_epithelium_Erysipelotrichaceae_genus, threshold = 0.1) # 7 taxa
write.csv(tax_table(lg_intestine_epithelium_Erysipelotrichaceae_genus_filt), "lg_intestine_epithelium_Erysipelotrichaceae_taxa.csv")
write.csv(otu_table(lg_intestine_epithelium_Erysipelotrichaceae_genus_filt), "lg_intestine_epithelium_Erysipelotrichaceae_otus.csv")
lg_intestine_epithelium_Erysipelotrichaceae_genus_melt <- psmelt(lg_intestine_epithelium_Erysipelotrichaceae_genus_filt)

lg_intestine_epithelium_Erysipelotrichaceae_palette <- c("#4B0082","#6030A8","grey88")
lg_intestine_epithelium_Erysipelotrichaceae_family <- subset_taxa(lg_intestine_epithelium_family, Family=="Erysipelotrichaceae") %>%
  psmelt()

ggplot() + theme_bw() + labs(y= "Erysipelotrichaceae (RA%)", title = "EPITHELIUM") +
  coord_flip() + theme_bw() +
  geom_bar(lg_intestine_epithelium_Erysipelotrichaceae_genus_melt, colour = "black", linewidth = 0.75,
           mapping = aes(x= treatment, y= Abundance, fill=Genus), stat = "summary") +
  #scale_fill_manual(values = lg_intestine_epithelium_Erysipelotrichaceae_palette) +
  scale_x_discrete(limits = c("tylosin","control"), labels =c("TYLOSIN","NO TYLOSIN")) +
  scale_y_continuous(expand = c(0.001,0,5.3,0)) +
  geom_errorbar(lg_intestine_epithelium_Erysipelotrichaceae_family, mapping= aes(x= treatment, y= Abundance),
                stat = "summary", width = 0.5, linewidth = 0.75) +
  theme(legend.position = "none",
        plot.title = element_text(size = 44),
        panel.border = element_rect(colour = "black", size = 1.25),
        panel.grid.major.y = element_blank(),
        axis.ticks = element_line(size = 0.75, colour= "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 34, colour = "black"),
        axis.title.x = element_text(size =44),
        axis.text.x = element_text(size = 16, colour = "black"))

kruskal.test(lg_intestine_epithelium_Erysipelotrichaceae_family$Abundance, lg_intestine_epithelium_Erysipelotrichaceae_family$treatment) #NS

################### MOST ABUNDANT GENERA IN RUMEN #####
#### LUMEN ####
rumen_lumen_genus_abund <- merge_low_abundance(rumen_lumen_genus, threshold = 0.5)
write.csv(otu_table(rumen_lumen_genus_abund),"top25genus_otus.csv")
write.csv(tax_table(rumen_lumen_genus_abund),"top25genus_taxa.csv")
rumen_lumen_genus_abund_melt <- psmelt(rumen_lumen_genus_abund)
rumen_lumen_genus_abund_order <- c("[Eubacterium]_coprostanoligenes_group",	"Succiniclasticum",	"[Eubacterium]_nodatum_group",	"Family_XIII_AD3011_group",	"Olsenella",	"Bacteroidales_RF16_group",	"Bifidobacterium",	"Clostridia_UCG-014",	"Desulfovibrio",	"UCG-004",	"Erysipelotrichaceae_UCG-009",	"Turicibacter",	"F082",	"[Ruminococcus]_gauvreauii_group",	"Acetitomaculum",	"unclassified Lachnospiraceae",	"Lachnospiraceae_NK3A20_group",	"Shuttleworthia",	"Syntrophococcus",	"[Eubacterium]_hallii_group",	"Oribacterium",	"Howardella",	"Methanobrevibacter",	"Muribaculaceae",	"NK4A214_group",	"UCG-002",	"p-2534-18B5_gut_group",	"Prevotella",	"unclassified Prevotellaceae",	"Prevotellaceae_UCG-001",	"Rikenellaceae_RC9_gut_group",	"Ruminococcus",	"unclassified Ruminococcaceae",	"uncultured",	"Succinivibrionaceae_UCG-001",	"Dialister",	"zzzOther ")

ggplot(rumen_lumen_genus_abund_melt, aes(x= Genus, y= Abundance)) +
  theme_bw() + labs(y= "Relative Abundance (%)", title = "LUMEN") +
  coord_flip() + 
  scale_x_discrete(limits = rev(rumen_lumen_genus_abund_order)) +
  scale_y_continuous(expand = c(0.002,0,0.987,0), labels = scales::number_format(accuracy = 0.1)) +
  geom_bar(aes(fill = treatment, colour = treatment), position = position_dodge2(padding = 0.28, reverse = T), stat = "summary", width = 0.72, alpha = 0.5, linewidth = 0.75) +
  geom_errorbar(aes(colour = treatment), position = position_dodge2(padding = 0.55, reverse = T), stat = "summary", width = 0.72, linewidth = 0.75) +
  scale_fill_manual(values = tylosin_palette) +
  scale_colour_manual(values = tylosin_palette) +
  theme(panel.border = element_rect(colour = "black", size = 1),
        plot.title = element_text(size = 34),
        legend.position = "none",
        panel.grid.major.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 16, colour = "black"),
        axis.title.x = element_text(size = 32),
        axis.text.x = element_text(size = 14, colour = "black"),
        axis.ticks = element_line(size = 0.75, colour = "black"))

### stats
#Eubacterium corpostanoligenes group
rumen_lumen_Eubacterium_corpostanoligenes_group <- rumen_lumen_genus_abund_melt[which(rumen_lumen_genus_abund_melt$Genus=="[Eubacterium]_coprostanoligenes_group"),]
kruskal.test(rumen_lumen_Eubacterium_corpostanoligenes_group$Abundance, rumen_lumen_Eubacterium_corpostanoligenes_group$treatment) # NS
#Succiniclasticum 
rumen_lumen_Succiniclasticum <- rumen_lumen_genus_abund_melt[which(rumen_lumen_genus_abund_melt$Genus=="Succiniclasticum"),]
kruskal.test(rumen_lumen_Succiniclasticum$Abundance, rumen_lumen_Succiniclasticum$treatment) # NS
#[Eubacterium]_nodatum_group
rumen_lumen_Eubacterium_nodatum_group <- rumen_lumen_genus_abund_melt[which(rumen_lumen_genus_abund_melt$Genus=="[Eubacterium]_nodatum_group"),]
kruskal.test(rumen_lumen_Eubacterium_nodatum_group$Abundance, rumen_lumen_Eubacterium_nodatum_group$treatment) # NS
#Family_XIII_AD3011_group
rumen_lumen_Family_XIII_AD3011_group <- rumen_lumen_genus_abund_melt[which(rumen_lumen_genus_abund_melt$Genus=="Family_XIII_AD3011_group"),]
kruskal.test(rumen_lumen_Family_XIII_AD3011_group$Abundance, rumen_lumen_Family_XIII_AD3011_group$treatment) # sig.
#Olsenella
rumen_lumen_Olsenella <- rumen_lumen_genus_abund_melt[which(rumen_lumen_genus_abund_melt$Genus=="Olsenella"),]
kruskal.test(rumen_lumen_Olsenella$Abundance, rumen_lumen_Olsenella$treatment) # NS
#Bacteroidales_RF16_group
rumen_lumen_Bacteroidales_RF16_group <- rumen_lumen_genus_abund_melt[which(rumen_lumen_genus_abund_melt$Genus=="Bacteroidales_RF16_group"),]
kruskal.test(rumen_lumen_Bacteroidales_RF16_group$Abundance, rumen_lumen_Bacteroidales_RF16_group$treatment) # NS
#Bifidobacterium
rumen_lumen_Bifidobacterium <- rumen_lumen_genus_abund_melt[which(rumen_lumen_genus_abund_melt$Genus=="Bifidobacterium"),]
kruskal.test(rumen_lumen_Bifidobacterium$Abundance, rumen_lumen_Bifidobacterium$treatment) # NS
#Clostridia_UCG-014
rumen_lumen_Clostridia_UCG014 <- rumen_lumen_genus_abund_melt[which(rumen_lumen_genus_abund_melt$Genus=="Clostridia_UCG-014"),]
kruskal.test(rumen_lumen_Clostridia_UCG014$Abundance, rumen_lumen_Clostridia_UCG014$treatment) # sig.
#Desulfovibrio
rumen_lumen_Desulfovibrio <- rumen_lumen_genus_abund_melt[which(rumen_lumen_genus_abund_melt$Genus=="Desulfovibrio"),]
kruskal.test(rumen_lumen_Desulfovibrio$Abundance, rumen_lumen_Desulfovibrio$treatment) # NS
#Erysipelatoclostridiaceae UCG004
rumen_lumen_UCG004 <- rumen_lumen_genus_abund_melt[which(rumen_lumen_genus_abund_melt$Genus=="UCG-004"),]
kruskal.test(rumen_lumen_UCG004$Abundance, rumen_lumen_UCG004$treatment) # sig.
#Erysipelotrichaceae_UCG-009
rumen_lumen_Erysipelotrichaceae_UCG009 <- rumen_lumen_genus_abund_melt[which(rumen_lumen_genus_abund_melt$Genus=="Erysipelotrichaceae_UCG-009"),]
kruskal.test(rumen_lumen_Erysipelotrichaceae_UCG009$Abundance, rumen_lumen_Erysipelotrichaceae_UCG009$treatment) # sig.
#Turicibacter
rumen_lumen_Turicibacter <- rumen_lumen_genus_abund_melt[which(rumen_lumen_genus_abund_melt$Genus=="Turicibacter"),]
kruskal.test(rumen_lumen_Turicibacter$Abundance, rumen_lumen_Turicibacter$treatment) # NS
#F082
rumen_lumen_F082 <- rumen_lumen_genus_abund_melt[which(rumen_lumen_genus_abund_melt$Genus=="F082"),]
kruskal.test(rumen_lumen_F082$Abundance, rumen_lumen_F082$treatment) # sig.
#ruminococcus_gauvreauii grp.
rumen_lumen_ruminococcus_gauvreauii <- rumen_lumen_genus_abund_melt[which(rumen_lumen_genus_abund_melt$Genus=="[Ruminococcus]_gauvreauii_group"),]
kruskal.test(rumen_lumen_ruminococcus_gauvreauii$Abundance, rumen_lumen_ruminococcus_gauvreauii$treatment) # sig.
#Acetitomaculum
rumen_lumen_Acetitomaculum <- rumen_lumen_genus_abund_melt[which(rumen_lumen_genus_abund_melt$Genus=="Acetitomaculum"),]
kruskal.test(rumen_lumen_Acetitomaculum$Abundance, rumen_lumen_Acetitomaculum$treatment) #NS
#Eubacteirum hallii group
rumen_lumen_Eubacterium_hallii_group <- rumen_lumen_genus_abund_melt[which(rumen_lumen_genus_abund_melt$Genus=="[Eubacterium]_hallii_group"),]
kruskal.test(rumen_lumen_Eubacterium_hallii_group$Abundance, rumen_lumen_Eubacterium_hallii_group$treatment) #NS
#Lachnospiraceae NK3A20 group
rumen_lumen_NK3A20_group <- rumen_lumen_genus_abund_melt[which(rumen_lumen_genus_abund_melt$Genus=="Lachnospiraceae_NK3A20_group"),]
kruskal.test(rumen_lumen_NK3A20_group$Abundance, rumen_lumen_NK3A20_group$treatment) #NS
#Shuttleworthia
rumen_lumen_Shuttleworthia <- rumen_lumen_genus_abund_melt[which(rumen_lumen_genus_abund_melt$Genus=="Shuttleworthia"),]
kruskal.test(rumen_lumen_Shuttleworthia$Abundance, rumen_lumen_Shuttleworthia$treatment) #sig
#Syntrophococcus
rumen_lumen_Syntrophococcus <- rumen_lumen_genus_abund_melt[which(rumen_lumen_genus_abund_melt$Genus=="Syntrophococcus"),]
kruskal.test(rumen_lumen_Syntrophococcus$Abundance, rumen_lumen_Syntrophococcus$treatment) #NS
# unclassified_Lachnospiraceae
rumen_lumen_unclassified_Lachnospiraceae <- rumen_lumen_genus_abund_melt[which(rumen_lumen_genus_abund_melt$Genus=="unclassified Lachnospiraceae"),]
kruskal.test(rumen_lumen_unclassified_Lachnospiraceae$Abundance, rumen_lumen_unclassified_Lachnospiraceae$treatment) #sig.
# oribacterium
rumen_lumen_oribacterium <- rumen_lumen_genus_abund_melt[which(rumen_lumen_genus_abund_melt$Genus=="Oribacterium"),]
kruskal.test(rumen_lumen_oribacterium$Abundance, rumen_lumen_oribacterium$treatment) #NS
# Howardella
rumen_lumen_Howardella <- rumen_lumen_genus_abund_melt[which(rumen_lumen_genus_abund_melt$Genus=="Howardella"),]
kruskal.test(rumen_lumen_Howardella$Abundance, rumen_lumen_Howardella$treatment) #NS
# Methanobrevibacter
rumen_lumen_Methanobrevibacter <- rumen_lumen_genus_abund_melt[which(rumen_lumen_genus_abund_melt$Genus=="Methanobrevibacter"),]
kruskal.test(rumen_lumen_Methanobrevibacter$Abundance, rumen_lumen_Methanobrevibacter$treatment) #NS
# Muribaculaceae
rumen_lumen_Muribaculaceae <- rumen_lumen_genus_abund_melt[which(rumen_lumen_genus_abund_melt$Genus=="Muribaculaceae"),]
kruskal.test(rumen_lumen_Muribaculaceae$Abundance, rumen_lumen_Muribaculaceae$treatment) #NS
# NK4A214_group
rumen_lumen_NK4A214_group <- rumen_lumen_genus_abund_melt[which(rumen_lumen_genus_abund_melt$Genus=="NK4A214_group"),]
kruskal.test(rumen_lumen_NK4A214_group$Abundance, rumen_lumen_NK4A214_group$treatment) #NS
# OscillospriraceaeUCG-002
rumen_lumen_UCG002 <- rumen_lumen_genus_abund_melt[which(rumen_lumen_genus_abund_melt$Genus=="UCG-002"),]
kruskal.test(rumen_lumen_UCG002$Abundance, rumen_lumen_UCG002$treatment) #NS
# p-2534-18B5_gut_group
rumen_lumen_p253418B5_gut_group <- rumen_lumen_genus_abund_melt[which(rumen_lumen_genus_abund_melt$Genus=="p-2534-18B5_gut_group"),]
kruskal.test(rumen_lumen_p253418B5_gut_group$Abundance, rumen_lumen_p253418B5_gut_group$treatment) #NS
# unclassified Prevotellaceae
rumen_lumen_unclassified_Prevotellaceae <- rumen_lumen_genus_abund_melt[which(rumen_lumen_genus_abund_melt$Genus=="unclassified Prevotellaceae"),]
kruskal.test(rumen_lumen_unclassified_Prevotellaceae$Abundance, rumen_lumen_unclassified_Prevotellaceae$treatment) #NS
# Prevotella
rumen_lumen_Prevotella <- rumen_lumen_genus_abund_melt[which(rumen_lumen_genus_abund_melt$Genus=="Prevotella"),]
kruskal.test(rumen_lumen_Prevotella$Abundance, rumen_lumen_Prevotella$treatment) #NS
# Prevotellaceae_UCG-001
rumen_lumen_Prevotellaceae_UCG001 <- rumen_lumen_genus_abund_melt[which(rumen_lumen_genus_abund_melt$Genus=="Prevotellaceae_UCG-001"),]
kruskal.test(rumen_lumen_Prevotellaceae_UCG001$Abundance, rumen_lumen_Prevotellaceae_UCG001$treatment) # sig
# Rikenellaceae_RC9_gut_group
rumen_lumen_Rikenellaceae_RC9_gut_group <- rumen_lumen_genus_abund_melt[which(rumen_lumen_genus_abund_melt$Genus=="Rikenellaceae_RC9_gut_group"),]
kruskal.test(rumen_lumen_Rikenellaceae_RC9_gut_group$Abundance, rumen_lumen_Rikenellaceae_RC9_gut_group$treatment) #NS
# Ruminococcus
rumen_lumen_Ruminococcus <- rumen_lumen_genus_abund_melt[which(rumen_lumen_genus_abund_melt$Genus=="Ruminococcus"),]
kruskal.test(rumen_lumen_Ruminococcus$Abundance, rumen_lumen_Ruminococcus$treatment) #NS
# unclassified Ruminococcaceae
rumen_lumen_unclassified_Ruminococcaceae <- rumen_lumen_genus_abund_melt[which(rumen_lumen_genus_abund_melt$Genus=="unclassified Ruminococcaceae"),]
kruskal.test(rumen_lumen_unclassified_Ruminococcaceae$Abundance, rumen_lumen_unclassified_Ruminococcaceae$treatment) #NS
# uncultured Selemonadaceae
rumen_lumen_uncultured_Selemonadaceae <- rumen_lumen_genus_abund_melt[which(rumen_lumen_genus_abund_melt$Genus=="uncultured"),]
kruskal.test(rumen_lumen_uncultured_Selemonadaceae$Abundance, rumen_lumen_uncultured_Selemonadaceae$treatment) #NS
# Succinivibrionaceae_UCG-001
rumen_lumen_uncultured_Succinivibrionaceae_UCG001 <- rumen_lumen_genus_abund_melt[which(rumen_lumen_genus_abund_melt$Genus=="Succinivibrionaceae_UCG-001"),]
kruskal.test(rumen_lumen_Succinivibrionaceae_UCG001$Abundance, rumen_lumen_Succinivibrionaceae_UCG001$treatment) #NS
# Dialister
rumen_lumen_Dialister <- rumen_lumen_genus_abund_melt[which(rumen_lumen_genus_abund_melt$Genus=="Dialister"),]
kruskal.test(rumen_lumen_Dialister$Abundance, rumen_lumen_Dialister$treatment) #NS
# low abundance genera (<0.5%)
rumen_lumen_zzzOther <- rumen_lumen_genus_abund_melt[which(rumen_lumen_genus_abund_melt$Genus=="zzzOther "),]
kruskal.test(rumen_lumen_zzzOther$Abundance, rumen_lumen_zzzOther$treatment) #NS

#### EPITHELIUM ####
rumen_epithelium_genus_abund <- merge_low_abundance(rumen_epithelium_genus, threshold = 0.5)
write.csv(otu_table(rumen_epithelium_genus_abund),"top25genus_otus.csv")
write.csv(tax_table(rumen_epithelium_genus_abund),"top25genus_taxa.csv")
rumen_epithelium_genus_abund_melt <- psmelt(rumen_epithelium_genus_abund)
rumen_epithelium_genus_abund_order <- c("[Eubacterium]_coprostanoligenes_group",	"Succiniclasticum",	"[Eubacterium]_nodatum_group",	"Family_XIII_AD3011_group",	"Olsenella",	"Bacteroidales_RF16_group",	"Clostridia_UCG-014",	"F082",	"[Ruminococcus]_gauvreauii_group",	"unclassified Lachnospiraceae",	"Shuttleworthia",	"Acetitomaculum",	"Lachnospiraceae_NK3A20_group",	"Oribacterium",	"Syntrophococcus",	"[Eubacterium]_hallii_group",	"Methanobrevibacter",	"Muribaculaceae",	"NK4A214_group",	"Prevotella",	"unclassified Prevotellaceae",	"Prevotellaceae_UCG-001",	"Rikenellaceae_RC9_gut_group",	"Ruminococcus",	"uncultured",	"Treponema",	"Succinivibrionaceae_UCG-001",	"Dialister",	"zzzOther ")

ggplot(rumen_epithelium_genus_abund_melt, aes(x= Genus, y= Abundance)) +
  theme_bw() + labs(y= "Relative Abundance (%)", title = "EPITHELIUM") +
  coord_flip() + 
  scale_x_discrete(limits = rev(rumen_epithelium_genus_abund_order)) +
  scale_y_continuous(expand = c(0.002,0,0.07,0), labels = scales::number_format(accuracy = 0.1)) +
  geom_bar(aes(fill = treatment, colour = treatment), position = position_dodge2(padding = 0.28, reverse = T), stat = "summary", width = 0.72, alpha = 0.5, linewidth = 0.75) +
  geom_errorbar(aes(colour = treatment), position = position_dodge2(padding = 0.55, reverse = T), stat = "summary", width = 0.72, linewidth = 0.75) +
  scale_fill_manual(values = tylosin_palette) +
  scale_colour_manual(values = tylosin_palette) +
  theme(panel.border = element_rect(colour = "black", size = 1),
        plot.title = element_text(size = 34),
        legend.position = "none",
        panel.grid.major.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 16, colour = "black"),
        axis.title.x = element_text(size = 32),
        axis.text.x = element_text(size = 14, colour = "black"),
        axis.ticks = element_line(size = 0.75, colour = "black"))

### stats
#Eubacterium corpostanoligenes group
rumen_epithelium_Eubacterium_corpostanoligenes_group <- rumen_epithelium_genus_abund_melt[which(rumen_epithelium_genus_abund_melt$Genus=="[Eubacterium]_coprostanoligenes_group"),]
kruskal.test(rumen_epithelium_Eubacterium_corpostanoligenes_group$Abundance, rumen_epithelium_Eubacterium_corpostanoligenes_group$treatment) # NS
#Succiniclasticum 
rumen_epithelium_Succiniclasticum <- rumen_epithelium_genus_abund_melt[which(rumen_epithelium_genus_abund_melt$Genus=="Succiniclasticum"),]
kruskal.test(rumen_epithelium_Succiniclasticum$Abundance, rumen_epithelium_Succiniclasticum$treatment) # NS
#[Eubacterium]_nodatum_group
rumen_epithelium_Eubacterium_nodatum_group <- rumen_epithelium_genus_abund_melt[which(rumen_epithelium_genus_abund_melt$Genus=="[Eubacterium]_nodatum_group"),]
kruskal.test(rumen_epithelium_Eubacterium_nodatum_group$Abundance, rumen_epithelium_Eubacterium_nodatum_group$treatment) # NS
#Family_XIII_AD3011_group
rumen_epithelium_Family_XIII_AD3011_group <- rumen_epithelium_genus_abund_melt[which(rumen_epithelium_genus_abund_melt$Genus=="Family_XIII_AD3011_group"),]
kruskal.test(rumen_epithelium_Family_XIII_AD3011_group$Abundance, rumen_epithelium_Family_XIII_AD3011_group$treatment) # sig.
#Olsenella
rumen_epithelium_Olsenella <- rumen_epithelium_genus_abund_melt[which(rumen_epithelium_genus_abund_melt$Genus=="Olsenella"),]
kruskal.test(rumen_epithelium_Olsenella$Abundance, rumen_epithelium_Olsenella$treatment) # NS
#Bacteroidales_RF16_group
rumen_epithelium_Bacteroidales_RF16_group <- rumen_epithelium_genus_abund_melt[which(rumen_epithelium_genus_abund_melt$Genus=="Bacteroidales_RF16_group"),]
kruskal.test(rumen_epithelium_Bacteroidales_RF16_group$Abundance, rumen_epithelium_Bacteroidales_RF16_group$treatment) # NS
#Bifidobacterium
rumen_epithelium_Bifidobacterium <- rumen_epithelium_genus_abund_melt[which(rumen_epithelium_genus_abund_melt$Genus=="Bifidobacterium"),]
kruskal.test(rumen_epithelium_Bifidobacterium$Abundance, rumen_epithelium_Bifidobacterium$treatment) # NS
#Clostridia_UCG-014
rumen_epithelium_Clostridia_UCG014 <- rumen_epithelium_genus_abund_melt[which(rumen_epithelium_genus_abund_melt$Genus=="Clostridia_UCG-014"),]
kruskal.test(rumen_epithelium_Clostridia_UCG014$Abundance, rumen_epithelium_Clostridia_UCG014$treatment) # sig.
#F082
rumen_epithelium_F082 <- rumen_epithelium_genus_abund_melt[which(rumen_epithelium_genus_abund_melt$Genus=="F082"),]
kruskal.test(rumen_epithelium_F082$Abundance, rumen_epithelium_F082$treatment) # sig.
#ruminococcus_gauvreauii grp.
rumen_epithelium_ruminococcus_gauvreauii <- rumen_epithelium_genus_abund_melt[which(rumen_epithelium_genus_abund_melt$Genus=="[Ruminococcus]_gauvreauii_group"),]
kruskal.test(rumen_epithelium_ruminococcus_gauvreauii$Abundance, rumen_epithelium_ruminococcus_gauvreauii$treatment) # sig.
#Acetitomaculum
rumen_epithelium_Acetitomaculum <- rumen_epithelium_genus_abund_melt[which(rumen_epithelium_genus_abund_melt$Genus=="Acetitomaculum"),]
kruskal.test(rumen_epithelium_Acetitomaculum$Abundance, rumen_epithelium_Acetitomaculum$treatment) #NS
#Eubacteirum hallii group
rumen_epithelium_Eubacterium_hallii_group <- rumen_epithelium_genus_abund_melt[which(rumen_epithelium_genus_abund_melt$Genus=="[Eubacterium]_hallii_group"),]
kruskal.test(rumen_epithelium_Eubacterium_hallii_group$Abundance, rumen_epithelium_Eubacterium_hallii_group$treatment) #NS
#Lachnospiraceae NK3A20 group
rumen_epithelium_NK3A20_group <- rumen_epithelium_genus_abund_melt[which(rumen_epithelium_genus_abund_melt$Genus=="Lachnospiraceae_NK3A20_group"),]
kruskal.test(rumen_epithelium_NK3A20_group$Abundance, rumen_epithelium_NK3A20_group$treatment) #NS
#Shuttleworthia
rumen_epithelium_Shuttleworthia <- rumen_epithelium_genus_abund_melt[which(rumen_epithelium_genus_abund_melt$Genus=="Shuttleworthia"),]
kruskal.test(rumen_epithelium_Shuttleworthia$Abundance, rumen_epithelium_Shuttleworthia$treatment) #sig
#Syntrophococcus
rumen_epithelium_Syntrophococcus <- rumen_epithelium_genus_abund_melt[which(rumen_epithelium_genus_abund_melt$Genus=="Syntrophococcus"),]
kruskal.test(rumen_epithelium_Syntrophococcus$Abundance, rumen_epithelium_Syntrophococcus$treatment) #NS
# unclassified_Lachnospiraceae
rumen_epithelium_unclassified_Lachnospiraceae <- rumen_epithelium_genus_abund_melt[which(rumen_epithelium_genus_abund_melt$Genus=="unclassified Lachnospiraceae"),]
kruskal.test(rumen_epithelium_unclassified_Lachnospiraceae$Abundance, rumen_epithelium_unclassified_Lachnospiraceae$treatment) #sig.
# oribacterium
rumen_epithelium_oribacterium <- rumen_epithelium_genus_abund_melt[which(rumen_epithelium_genus_abund_melt$Genus=="Oribacterium"),]
kruskal.test(rumen_epithelium_oribacterium$Abundance, rumen_epithelium_oribacterium$treatment) #NS
# Methanobrevibacter
rumen_epithelium_Methanobrevibacter <- rumen_epithelium_genus_abund_melt[which(rumen_epithelium_genus_abund_melt$Genus=="Methanobrevibacter"),]
kruskal.test(rumen_epithelium_Methanobrevibacter$Abundance, rumen_epithelium_Methanobrevibacter$treatment) #NS
# Muribaculaceae
rumen_epithelium_Muribaculaceae <- rumen_epithelium_genus_abund_melt[which(rumen_epithelium_genus_abund_melt$Genus=="Muribaculaceae"),]
kruskal.test(rumen_epithelium_Muribaculaceae$Abundance, rumen_epithelium_Muribaculaceae$treatment) #NS
# NK4A214_group
rumen_epithelium_NK4A214_group <- rumen_epithelium_genus_abund_melt[which(rumen_epithelium_genus_abund_melt$Genus=="NK4A214_group"),]
kruskal.test(rumen_epithelium_NK4A214_group$Abundance, rumen_epithelium_NK4A214_group$treatment) #NS
# unclassified Prevotellaceae
rumen_epithelium_unclassified_Prevotellaceae <- rumen_epithelium_genus_abund_melt[which(rumen_epithelium_genus_abund_melt$Genus=="unclassified Prevotellaceae"),]
kruskal.test(rumen_epithelium_unclassified_Prevotellaceae$Abundance, rumen_epithelium_unclassified_Prevotellaceae$treatment) #NS
# Prevotella
rumen_epithelium_Prevotella <- rumen_epithelium_genus_abund_melt[which(rumen_epithelium_genus_abund_melt$Genus=="Prevotella"),]
kruskal.test(rumen_epithelium_Prevotella$Abundance, rumen_epithelium_Prevotella$treatment) #NS
# Prevotellaceae_UCG-001
rumen_epithelium_Prevotellaceae_UCG001 <- rumen_epithelium_genus_abund_melt[which(rumen_epithelium_genus_abund_melt$Genus=="Prevotellaceae_UCG-001"),]
kruskal.test(rumen_epithelium_Prevotellaceae_UCG001$Abundance, rumen_epithelium_Prevotellaceae_UCG001$treatment) # sig
# Rikenellaceae_RC9_gut_group
rumen_epithelium_Rikenellaceae_RC9_gut_group <- rumen_epithelium_genus_abund_melt[which(rumen_epithelium_genus_abund_melt$Genus=="Rikenellaceae_RC9_gut_group"),]
kruskal.test(rumen_epithelium_Rikenellaceae_RC9_gut_group$Abundance, rumen_epithelium_Rikenellaceae_RC9_gut_group$treatment) #NS
# Ruminococcus
rumen_epithelium_Ruminococcus <- rumen_epithelium_genus_abund_melt[which(rumen_epithelium_genus_abund_melt$Genus=="Ruminococcus"),]
kruskal.test(rumen_epithelium_Ruminococcus$Abundance, rumen_epithelium_Ruminococcus$treatment) #NS
# uncultured Selemonadaceae
rumen_epithelium_uncultured_Selemonadaceae <- rumen_epithelium_genus_abund_melt[which(rumen_epithelium_genus_abund_melt$Genus=="uncultured"),]
kruskal.test(rumen_epithelium_uncultured_Selemonadaceae$Abundance, rumen_epithelium_uncultured_Selemonadaceae$treatment) #NS
# Treponema
rumen_epithelium_Treponema <- rumen_epithelium_genus_abund_melt[which(rumen_epithelium_genus_abund_melt$Genus=="Treponema"),]
kruskal.test(rumen_epithelium_Treponema$Abundance, rumen_epithelium_Treponema$treatment) #sig
# Succinivibrionaceae_UCG-001
rumen_epithelium_Succinivibrionaceae_UCG001 <- rumen_epithelium_genus_abund_melt[which(rumen_epithelium_genus_abund_melt$Genus=="Succinivibrionaceae_UCG-001"),]
kruskal.test(rumen_epithelium_Succinivibrionaceae_UCG001$Abundance, rumen_epithelium_Succinivibrionaceae_UCG001$treatment) #NS
# Dialister
rumen_epithelium_Dialister <- rumen_epithelium_genus_abund_melt[which(rumen_epithelium_genus_abund_melt$Genus=="Dialister"),]
kruskal.test(rumen_epithelium_Dialister$Abundance, rumen_epithelium_Dialister$treatment) #NS
# low abundance genera (<0.5%)
rumen_epithelium_zzzOther <- rumen_epithelium_genus_abund_melt[which(rumen_epithelium_genus_abund_melt$Genus=="zzzOther "),]
kruskal.test(rumen_epithelium_zzzOther$Abundance, rumen_epithelium_zzzOther$treatment) #NS

################### MOST ABUNDANT GENERA IN SM INTESTINE ####
#### LUMEN ####
sm_intestine_lumen_genus_abund <- merge_low_abundance(sm_intestine_lumen_genus, threshold = 0.5)
write.csv(otu_table(sm_intestine_lumen_genus_abund),"top25genus_otus.csv")
write.csv(tax_table(sm_intestine_lumen_genus_abund),"top25genus_taxa.csv")
sm_intestine_lumen_genus_abund_melt <- psmelt(sm_intestine_lumen_genus_abund)
sm_intestine_lumen_genus_abund_order <- c("Family_XIII_AD3011_group",	"Olsenella",	"Aeriscardovia",	"Bifidobacteriaceae",	"Clostridium_sensu_stricto_1",	"Escherichia-Shigella",	"Turicibacter",	"[Ruminococcus]_gauvreauii_group",	"unclassified Lachnospiraceae",	"Lachnospiraceae_NK3A20_group",	"Acetitomaculum",	"[Eubacterium]_hallii_group",	"Syntrophococcus",	"Howardella",	"Lachnobacterium",	"Methanobrevibacter",	"Romboutsia",	"Paeniclostridium",	"unclassified Peptostreptococcaceae",	"Ruminococcus",	"unclassified Ruminococcaceae",	"zzzOther ")

ggplot(sm_intestine_lumen_genus_abund_melt, aes(x= Genus, y= Abundance)) +
  theme_bw() + labs(y= "Relative Abundance (%)", title = "LUMEN") +
  coord_flip() + 
  scale_x_discrete(limits = rev(sm_intestine_lumen_genus_abund_order)) +
  scale_y_continuous(expand = c(0.002,0,0.258,0), labels = scales::number_format(accuracy = 0.1)) +
  geom_bar(aes(fill = treatment, colour = treatment), position = position_dodge2(padding = 0.28, reverse = T), stat = "summary", width = 0.72, alpha = 0.5, linewidth = 0.75) +
  geom_errorbar(aes(colour = treatment), position = position_dodge2(padding = 0.55, reverse = T), stat = "summary", width = 0.72, linewidth = 0.75) +
  scale_fill_manual(values = tylosin_palette) +
  scale_colour_manual(values = tylosin_palette) +
  theme(panel.border = element_rect(colour = "black", size = 1),
        plot.title = element_text(size = 34),
        legend.position = "none",
        panel.grid.major.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 16, colour = "black"),
        axis.title.x = element_text(size = 32),
        axis.text.x = element_text(size = 14, colour = "black"),
        axis.ticks = element_line(size = 0.75, colour = "black"))

### stats
#Family_XIII_AD3011_group
sm_intestine_lumen_Family_XIII_AD3011_group <- sm_intestine_lumen_genus_abund_melt[which(sm_intestine_lumen_genus_abund_melt$Genus=="Family_XIII_AD3011_group"),]
kruskal.test(sm_intestine_lumen_Family_XIII_AD3011_group$Abundance, sm_intestine_lumen_Family_XIII_AD3011_group$treatment) # sig.
#Olsenella
sm_intestine_lumen_Olsenella <- sm_intestine_lumen_genus_abund_melt[which(sm_intestine_lumen_genus_abund_melt$Genus=="Olsenella"),]
kruskal.test(sm_intestine_lumen_Olsenella$Abundance, sm_intestine_lumen_Olsenella$treatment) # NS
# Aeriscardovia
sm_intestine_lumen_Aeriscardovia <- sm_intestine_lumen_genus_abund_melt[which(sm_intestine_lumen_genus_abund_melt$Genus=="Aeriscardovia"),]
kruskal.test(sm_intestine_lumen_Aeriscardovia$Abundance, sm_intestine_lumen_Aeriscardovia$treatment) # NS
# Bifidobacterium
sm_intestine_lumen_Bifidobacterium <- sm_intestine_lumen_genus_abund_melt[which(sm_intestine_lumen_genus_abund_melt$Genus=="Bifidobacteriaceae"),]
kruskal.test(sm_intestine_lumen_Bifidobacterium$Abundance, sm_intestine_lumen_Bifidobacterium$treatment) # NS
# Clostridium_sensu_stricto_1
sm_intestine_lumen_Clostridium_sensu_stricto_1 <- sm_intestine_lumen_genus_abund_melt[which(sm_intestine_lumen_genus_abund_melt$Genus=="Clostridium_sensu_stricto_1"),]
kruskal.test(sm_intestine_lumen_Clostridium_sensu_stricto_1$Abundance, sm_intestine_lumen_Clostridium_sensu_stricto_1$treatment) # NS
# Escherichia-Shigella
sm_intestine_lumen_EscherichiaShigella <- sm_intestine_lumen_genus_abund_melt[which(sm_intestine_lumen_genus_abund_melt$Genus=="Escherichia-Shigella"),]
kruskal.test(sm_intestine_lumen_EscherichiaShigella$Abundance, sm_intestine_lumen_EscherichiaShigella$treatment) # NS
# Turicibacter
sm_intestine_lumen_Turicibacter <- sm_intestine_lumen_genus_abund_melt[which(sm_intestine_lumen_genus_abund_melt$Genus=="Turicibacter"),]
kruskal.test(sm_intestine_lumen_Turicibacter$Abundance, sm_intestine_lumen_Turicibacter$treatment) # NS
#ruminococcus_gauvreauii grp.
sm_intestine_lumen_ruminococcus_gauvreauii <- sm_intestine_lumen_genus_abund_melt[which(sm_intestine_lumen_genus_abund_melt$Genus=="[Ruminococcus]_gauvreauii_group"),]
kruskal.test(sm_intestine_lumen_ruminococcus_gauvreauii$Abundance, sm_intestine_lumen_ruminococcus_gauvreauii$treatment) # sig.
#Acetitomaculum
sm_intestine_lumen_Acetitomaculum <- sm_intestine_lumen_genus_abund_melt[which(sm_intestine_lumen_genus_abund_melt$Genus=="Acetitomaculum"),]
kruskal.test(sm_intestine_lumen_Acetitomaculum$Abundance, sm_intestine_lumen_Acetitomaculum$treatment) #NS
#Eubacteirum hallii group
sm_intestine_lumen_Eubacterium_hallii_group <- sm_intestine_lumen_genus_abund_melt[which(sm_intestine_lumen_genus_abund_melt$Genus=="[Eubacterium]_hallii_group"),]
kruskal.test(sm_intestine_lumen_Eubacterium_hallii_group$Abundance, sm_intestine_lumen_Eubacterium_hallii_group$treatment) #NS
#Lachnospiraceae NK3A20 group
sm_intestine_lumen_NK3A20_group <- sm_intestine_lumen_genus_abund_melt[which(sm_intestine_lumen_genus_abund_melt$Genus=="Lachnospiraceae_NK3A20_group"),]
kruskal.test(sm_intestine_lumen_NK3A20_group$Abundance, sm_intestine_lumen_NK3A20_group$treatment) #NS
#Syntrophococcus
sm_intestine_lumen_Syntrophococcus <- sm_intestine_lumen_genus_abund_melt[which(sm_intestine_lumen_genus_abund_melt$Genus=="Syntrophococcus"),]
kruskal.test(sm_intestine_lumen_Syntrophococcus$Abundance, sm_intestine_lumen_Syntrophococcus$treatment) #NS
# unclassified_Lachnospiraceae
sm_intestine_lumen_unclassified_Lachnospiraceae <- sm_intestine_lumen_genus_abund_melt[which(sm_intestine_lumen_genus_abund_melt$Genus=="unclassified Lachnospiraceae"),]
kruskal.test(sm_intestine_lumen_unclassified_Lachnospiraceae$Abundance, sm_intestine_lumen_unclassified_Lachnospiraceae$treatment) #sig.
# Lachnobacteriuma
sm_intestine_lumen_Lachnobacterium <- sm_intestine_lumen_genus_abund_melt[which(sm_intestine_lumen_genus_abund_melt$Genus=="Lachnobacterium"),]
kruskal.test(sm_intestine_lumen_Lachnobacterium$Abundance, sm_intestine_lumen_Lachnobacterium$treatment) #NS
# Howardella
sm_intestine_lumen_Howardella <- sm_intestine_lumen_genus_abund_melt[which(sm_intestine_lumen_genus_abund_melt$Genus=="Howardella"),]
kruskal.test(sm_intestine_lumen_Howardella$Abundance, sm_intestine_lumen_Howardella$treatment) #NS
# Methanobrevibacter
sm_intestine_lumen_Methanobrevibacter <- sm_intestine_lumen_genus_abund_melt[which(sm_intestine_lumen_genus_abund_melt$Genus=="Methanobrevibacter"),]
kruskal.test(sm_intestine_lumen_Methanobrevibacter$Abundance, sm_intestine_lumen_Methanobrevibacter$treatment) #NS
# Romboutsia
sm_intestine_lumen_Romboutsia <- sm_intestine_lumen_genus_abund_melt[which(sm_intestine_lumen_genus_abund_melt$Genus=="Romboutsia"),]
kruskal.test(sm_intestine_lumen_Romboutsia$Abundance, sm_intestine_lumen_Romboutsia$treatment) #NS
# Paeniclostridium
sm_intestine_lumen_Paeniclostridium <- sm_intestine_lumen_genus_abund_melt[which(sm_intestine_lumen_genus_abund_melt$Genus=="Paeniclostridium"),]
kruskal.test(sm_intestine_lumen_Paeniclostridium$Abundance, sm_intestine_lumen_Paeniclostridium$treatment) #NS
# unclassified Peptostreptococcaceae
sm_intestine_lumen_unclassified_Peptostreptococcaceae <- sm_intestine_lumen_genus_abund_melt[which(sm_intestine_lumen_genus_abund_melt$Genus=="unclassified Peptostreptococcaceae"),]
kruskal.test(sm_intestine_lumen_unclassified_Peptostreptococcaceae$Abundance, sm_intestine_lumen_unclassified_Peptostreptococcaceae$treatment) #NS
# Ruminococcus
sm_intestine_lumen_Ruminococcus <- sm_intestine_lumen_genus_abund_melt[which(sm_intestine_lumen_genus_abund_melt$Genus=="Ruminococcus"),]
kruskal.test(sm_intestine_lumen_Ruminococcus$Abundance, sm_intestine_lumen_Ruminococcus$treatment) #sig
# unclassified Ruminococcaceae
sm_intestine_lumen_unclassified_Ruminococcaceae <- sm_intestine_lumen_genus_abund_melt[which(sm_intestine_lumen_genus_abund_melt$Genus=="unclassified Ruminococcaceae"),]
kruskal.test(sm_intestine_lumen_unclassified_Ruminococcaceae$Abundance, sm_intestine_lumen_unclassified_Ruminococcaceae$treatment) #NS
# low abundance genera (<0.5%)
sm_intestine_lumen_zzzOther <- sm_intestine_lumen_genus_abund_melt[which(sm_intestine_lumen_genus_abund_melt$Genus=="zzzOther "),]
kruskal.test(sm_intestine_lumen_zzzOther$Abundance, sm_intestine_lumen_zzzOther$treatment) #NS

#### EPITHELIUM ####
sm_intestine_epithelium_genus_abund <- merge_low_abundance(sm_intestine_epithelium_genus, threshold = 0.5)
write.csv(otu_table(sm_intestine_epithelium_genus_abund),"top25genus_otus.csv")
write.csv(tax_table(sm_intestine_epithelium_genus_abund),"top25genus_taxa.csv")
sm_intestine_epithelium_genus_abund_melt <- psmelt(sm_intestine_epithelium_genus_abund)
sm_intestine_epithelium_genus_abund_order <- c("Family_XIII_AD3011_group",	"Mogibacterium",	"[Eubacterium]_nodatum_group",	"Olsenella",	"Clostridium_sensu_stricto_1",	"Escherichia-Shigella",	"Turicibacter",	"[Ruminococcus]_gauvreauii_group",	"unclassified Lachnospiraceae",	"Lachnospiraceae_NK3A20_group",	"Acetitomaculum",	"[Eubacterium]_hallii_group",	"Syntrophococcus",	"Howardella",	"Methanobrevibacter",	"Romboutsia",	"Paeniclostridium",	"unclassified Peptostreptococcaceae",	"Ruminococcus",	"unclassified Ruminococcaceae",	"zzzOther ")

ggplot(sm_intestine_epithelium_genus_abund_melt, aes(x= Genus, y= Abundance)) +
  theme_bw() + labs(y= "Relative Abundance (%)", title = "EPITHELIUM") +
  coord_flip() + 
  scale_x_discrete(limits = rev(sm_intestine_epithelium_genus_abund_order)) +
  scale_y_continuous(expand = c(0.002,0,0.12,0), labels = scales::number_format(accuracy = 0.1)) +
  geom_bar(aes(fill = treatment, colour = treatment), position = position_dodge2(padding = 0.28, reverse = T), stat = "summary", width = 0.72, alpha = 0.5, linewidth = 0.75) +
  geom_errorbar(aes(colour = treatment), position = position_dodge2(padding = 0.55, reverse = T), stat = "summary", width = 0.72, linewidth = 0.75) +
  scale_fill_manual(values = tylosin_palette) +
  scale_colour_manual(values = tylosin_palette) +
  theme(panel.border = element_rect(colour = "black", size = 1),
        plot.title = element_text(size = 34),
        legend.position = "none",
        panel.grid.major.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 16, colour = "black"),
        axis.title.x = element_text(size = 32),
        axis.text.x = element_text(size = 14, colour = "black"),
        axis.ticks = element_line(size = 0.75, colour = "black"))

### stats
#Family_XIII_AD3011_group
sm_intestine_epithelium_Family_XIII_AD3011_group <- sm_intestine_epithelium_genus_abund_melt[which(sm_intestine_epithelium_genus_abund_melt$Genus=="Family_XIII_AD3011_group"),]
kruskal.test(sm_intestine_epithelium_Family_XIII_AD3011_group$Abundance, sm_intestine_epithelium_Family_XIII_AD3011_group$treatment) # sig.
#Mogibacterium
sm_intestine_epithelium_Mogibacterium <- sm_intestine_epithelium_genus_abund_melt[which(sm_intestine_epithelium_genus_abund_melt$Genus=="Mogibacterium"),]
kruskal.test(sm_intestine_epithelium_Mogibacterium$Abundance, sm_intestine_epithelium_Mogibacterium$treatment) # NS
#[Eubacterium]_nodatum_group
sm_intestine_epithelium_Eubacterium_nodatum_group <- sm_intestine_epithelium_genus_abund_melt[which(sm_intestine_epithelium_genus_abund_melt$Genus=="[Eubacterium]_nodatum_group"),]
kruskal.test(sm_intestine_epithelium_Eubacterium_nodatum_group$Abundance, sm_intestine_epithelium_Eubacterium_nodatum_group$treatment) # NS
#Olsenella
sm_intestine_epithelium_Olsenella <- sm_intestine_epithelium_genus_abund_melt[which(sm_intestine_epithelium_genus_abund_melt$Genus=="Olsenella"),]
kruskal.test(sm_intestine_epithelium_Olsenella$Abundance, sm_intestine_epithelium_Olsenella$treatment) # NS
# Clostridium_sensu_stricto_1
sm_intestine_epithelium_Clostridium_sensu_stricto_1 <- sm_intestine_epithelium_genus_abund_melt[which(sm_intestine_epithelium_genus_abund_melt$Genus=="Clostridium_sensu_stricto_1"),]
kruskal.test(sm_intestine_epithelium_Clostridium_sensu_stricto_1$Abundance, sm_intestine_epithelium_Clostridium_sensu_stricto_1$treatment) # NS
# Escherichia-Shigella
sm_intestine_epithelium_EscherichiaShigella <- sm_intestine_epithelium_genus_abund_melt[which(sm_intestine_epithelium_genus_abund_melt$Genus=="Escherichia-Shigella"),]
kruskal.test(sm_intestine_epithelium_EscherichiaShigella$Abundance, sm_intestine_epithelium_EscherichiaShigella$treatment) # NS
# Turicibacter
sm_intestine_epithelium_Turicibacter <- sm_intestine_epithelium_genus_abund_melt[which(sm_intestine_epithelium_genus_abund_melt$Genus=="Turicibacter"),]
kruskal.test(sm_intestine_epithelium_Turicibacter$Abundance, sm_intestine_epithelium_Turicibacter$treatment) # NS
#ruminococcus_gauvreauii grp.
sm_intestine_epithelium_ruminococcus_gauvreauii <- sm_intestine_epithelium_genus_abund_melt[which(sm_intestine_epithelium_genus_abund_melt$Genus=="[Ruminococcus]_gauvreauii_group"),]
kruskal.test(sm_intestine_epithelium_ruminococcus_gauvreauii$Abundance, sm_intestine_epithelium_ruminococcus_gauvreauii$treatment) # sig.
#Acetitomaculum
sm_intestine_epithelium_Acetitomaculum <- sm_intestine_epithelium_genus_abund_melt[which(sm_intestine_epithelium_genus_abund_melt$Genus=="Acetitomaculum"),]
kruskal.test(sm_intestine_epithelium_Acetitomaculum$Abundance, sm_intestine_epithelium_Acetitomaculum$treatment) #NS
#Eubacteirum hallii group
sm_intestine_epithelium_Eubacterium_hallii_group <- sm_intestine_epithelium_genus_abund_melt[which(sm_intestine_epithelium_genus_abund_melt$Genus=="[Eubacterium]_hallii_group"),]
kruskal.test(sm_intestine_epithelium_Eubacterium_hallii_group$Abundance, sm_intestine_epithelium_Eubacterium_hallii_group$treatment) #NS
#Lachnospiraceae NK3A20 group
sm_intestine_epithelium_NK3A20_group <- sm_intestine_epithelium_genus_abund_melt[which(sm_intestine_epithelium_genus_abund_melt$Genus=="Lachnospiraceae_NK3A20_group"),]
kruskal.test(sm_intestine_epithelium_NK3A20_group$Abundance, sm_intestine_epithelium_NK3A20_group$treatment) #NS
#Syntrophococcus
sm_intestine_epithelium_Syntrophococcus <- sm_intestine_epithelium_genus_abund_melt[which(sm_intestine_epithelium_genus_abund_melt$Genus=="Syntrophococcus"),]
kruskal.test(sm_intestine_epithelium_Syntrophococcus$Abundance, sm_intestine_epithelium_Syntrophococcus$treatment) #NS
# unclassified_Lachnospiraceae
sm_intestine_epithelium_unclassified_Lachnospiraceae <- sm_intestine_epithelium_genus_abund_melt[which(sm_intestine_epithelium_genus_abund_melt$Genus=="unclassified Lachnospiraceae"),]
kruskal.test(sm_intestine_epithelium_unclassified_Lachnospiraceae$Abundance, sm_intestine_epithelium_unclassified_Lachnospiraceae$treatment) #sig.
# Howardella
sm_intestine_epithelium_Howardella <- sm_intestine_epithelium_genus_abund_melt[which(sm_intestine_epithelium_genus_abund_melt$Genus=="Howardella"),]
kruskal.test(sm_intestine_epithelium_Howardella$Abundance, sm_intestine_epithelium_Howardella$treatment) #NS
# Methanobrevibacter
sm_intestine_epithelium_Methanobrevibacter <- sm_intestine_epithelium_genus_abund_melt[which(sm_intestine_epithelium_genus_abund_melt$Genus=="Methanobrevibacter"),]
kruskal.test(sm_intestine_epithelium_Methanobrevibacter$Abundance, sm_intestine_epithelium_Methanobrevibacter$treatment) #NS
# Romboutsia
sm_intestine_epithelium_Romboutsia <- sm_intestine_epithelium_genus_abund_melt[which(sm_intestine_epithelium_genus_abund_melt$Genus=="Romboutsia"),]
kruskal.test(sm_intestine_epithelium_Romboutsia$Abundance, sm_intestine_epithelium_Romboutsia$treatment) #NS
# Paeniclostridium
sm_intestine_epithelium_Paeniclostridium <- sm_intestine_epithelium_genus_abund_melt[which(sm_intestine_epithelium_genus_abund_melt$Genus=="Paeniclostridium"),]
kruskal.test(sm_intestine_epithelium_Paeniclostridium$Abundance, sm_intestine_epithelium_Paeniclostridium$treatment) #NS
# unclassified Peptostreptococcaceae
sm_intestine_epithelium_unclassified_Peptostreptococcaceae <- sm_intestine_epithelium_genus_abund_melt[which(sm_intestine_epithelium_genus_abund_melt$Genus=="unclassified Peptostreptococcaceae"),]
kruskal.test(sm_intestine_epithelium_unclassified_Peptostreptococcaceae$Abundance, sm_intestine_epithelium_unclassified_Peptostreptococcaceae$treatment) #NS
# Ruminococcus
sm_intestine_epithelium_Ruminococcus <- sm_intestine_epithelium_genus_abund_melt[which(sm_intestine_epithelium_genus_abund_melt$Genus=="Ruminococcus"),]
kruskal.test(sm_intestine_epithelium_Ruminococcus$Abundance, sm_intestine_epithelium_Ruminococcus$treatment) #sig
# unclassified Ruminococcaceae
sm_intestine_epithelium_unclassified_Ruminococcaceae <- sm_intestine_epithelium_genus_abund_melt[which(sm_intestine_epithelium_genus_abund_melt$Genus=="unclassified Ruminococcaceae"),]
kruskal.test(sm_intestine_epithelium_unclassified_Ruminococcaceae$Abundance, sm_intestine_epithelium_unclassified_Ruminococcaceae$treatment) #NS
# low abundance genera (<0.5%)
sm_intestine_epithelium_zzzOther <- sm_intestine_epithelium_genus_abund_melt[which(sm_intestine_epithelium_genus_abund_melt$Genus=="zzzOther "),]
kruskal.test(sm_intestine_epithelium_zzzOther$Abundance, sm_intestine_epithelium_zzzOther$treatment) #NS

################### MOST ABUNDANT GENERA IN LG INTESTINE ####
#### LUMEN ####
lg_intestine_lumen_genus_abund <- merge_low_abundance(lg_intestine_lumen_genus, threshold = 0.5)
write.csv(otu_table(lg_intestine_lumen_genus_abund),"top25genus_otus.csv")
write.csv(tax_table(lg_intestine_lumen_genus_abund),"top25genus_taxa.csv")
lg_intestine_lumen_genus_abund_melt <- psmelt(lg_intestine_lumen_genus_abund)
lg_intestine_lumen_genus_abund_order <- c("[Eubacterium]_coprostanoligenes_group",	"Phascolarctobacterium",	"Family_XIII_AD3011_group",	"Olsenella",	"Bacteroides",	"Bifidobacterium",	"Christensenellaceae_R-7_group",	"Clostridium_sensu_stricto_1",	"Turicibacter",	"[Ruminococcus]_gauvreauii_group",	"unclassified Lachnospiraceae",	"Lachnospiraceae_NK3A20_group",	"Blautia",	"Acetitomaculum",	"[Eubacterium]_hallii_group",	"Methanobrevibacter",	"Monoglobus",	"Muribaculaceae",	"UCG-005",	"UCG-010",	"Romboutsia",	"Paeniclostridium",	"unclassified Peptostreptococcaceae",	"Prevotellaceae_UCG-003",	"Prevotella",	"Alloprevotella",	"Prevotellaceae_NK3B31_group",	"uncultured",	"Rikenellaceae_RC9_gut_group",	"Alistipes",	"Ruminococcus",	"Treponema",	"Parabacteroides",	"zzzOther ")

ggplot(lg_intestine_lumen_genus_abund_melt, aes(x= Genus, y= Abundance)) +
  theme_bw() + labs(y= "Relative Abundance (%)", title = "LUMEN") +
  coord_flip() + 
  scale_x_discrete(limits = rev(lg_intestine_lumen_genus_abund_order)) +
  scale_y_continuous(expand = c(0.002,0,0.113,0), labels = scales::number_format(accuracy = 0.1)) +
  geom_bar(aes(fill = treatment, colour = treatment), position = position_dodge2(padding = 0.28, reverse = T), stat = "summary", width = 0.72, alpha = 0.5, linewidth = 0.75) +
  geom_errorbar(aes(colour = treatment), position = position_dodge2(padding = 0.55, reverse = T), stat = "summary", width = 0.72, linewidth = 0.75) +
  scale_fill_manual(values = tylosin_palette) +
  scale_colour_manual(values = tylosin_palette) +
  theme(panel.border = element_rect(colour = "black", size = 1),
        plot.title = element_text(size = 34),
        legend.position = "none",
        panel.grid.major.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 16, colour = "black"),
        axis.title.x = element_text(size = 32),
        axis.text.x = element_text(size = 14, colour = "black"),
        axis.ticks = element_line(size = 0.75, colour = "black"))

### stats
#[Eubacterium]_coprostanoligenes_group
lg_intestine_lumen_Eubacterium_coprostanoligenes_group <- lg_intestine_lumen_genus_abund_melt[which(lg_intestine_lumen_genus_abund_melt$Genus=="[Eubacterium]_coprostanoligenes_group"),]
kruskal.test(lg_intestine_lumen_Eubacterium_coprostanoligenes_group$Abundance, lg_intestine_lumen_Eubacterium_coprostanoligenes_group$treatment) #NS
#Phascolarctobacterium
lg_intestine_lumen_Phascolarctobacterium <- lg_intestine_lumen_genus_abund_melt[which(lg_intestine_lumen_genus_abund_melt$Genus=="Phascolarctobacterium"),]
kruskal.test(lg_intestine_lumen_Phascolarctobacterium$Abundance, lg_intestine_lumen_Phascolarctobacterium$treatment) # sig.
#Family_XIII_AD3011_group
lg_intestine_lumen_Family_XIII_AD3011_group <- lg_intestine_lumen_genus_abund_melt[which(lg_intestine_lumen_genus_abund_melt$Genus=="Family_XIII_AD3011_group"),]
kruskal.test(lg_intestine_lumen_Family_XIII_AD3011_group$Abundance, lg_intestine_lumen_Family_XIII_AD3011_group$treatment) # sig.
#Olsenella
lg_intestine_lumen_Olsenella <- lg_intestine_lumen_genus_abund_melt[which(lg_intestine_lumen_genus_abund_melt$Genus=="Olsenella"),]
kruskal.test(lg_intestine_lumen_Olsenella$Abundance, lg_intestine_lumen_Olsenella$treatment) # sig
# Bacteroides
lg_intestine_lumen_Bacteroides <- lg_intestine_lumen_genus_abund_melt[which(lg_intestine_lumen_genus_abund_melt$Genus=="Bacteroides"),]
kruskal.test(lg_intestine_lumen_Bacteroides$Abundance, lg_intestine_lumen_Bacteroides$treatment) # NS
# Bifidobacterium
lg_intestine_lumen_Bifidobacterium <- lg_intestine_lumen_genus_abund_melt[which(lg_intestine_lumen_genus_abund_melt$Genus=="Bifidobacterium"),]
kruskal.test(lg_intestine_lumen_Bifidobacterium$Abundance, lg_intestine_lumen_Bifidobacterium$treatment) # NS
# Christensenellaceae_R-7_group
lg_intestine_lumen_Christensenellaceae_R7_group <- lg_intestine_lumen_genus_abund_melt[which(lg_intestine_lumen_genus_abund_melt$Genus=="Christensenellaceae_R-7_group"),]
kruskal.test(lg_intestine_lumen_Christensenellaceae_R7_group$Abundance, lg_intestine_lumen_Christensenellaceae_R7_group$treatment) # NS
# Turicibacter
lg_intestine_lumen_Turicibacter <- lg_intestine_lumen_genus_abund_melt[which(lg_intestine_lumen_genus_abund_melt$Genus=="Turicibacter"),]
kruskal.test(lg_intestine_lumen_Turicibacter$Abundance, lg_intestine_lumen_Turicibacter$treatment) # NS
#ruminococcus_gauvreauii grp.
lg_intestine_lumen_ruminococcus_gauvreauii <- lg_intestine_lumen_genus_abund_melt[which(lg_intestine_lumen_genus_abund_melt$Genus=="[Ruminococcus]_gauvreauii_group"),]
kruskal.test(lg_intestine_lumen_ruminococcus_gauvreauii$Abundance, lg_intestine_lumen_ruminococcus_gauvreauii$treatment) # sig.
# unclassified_Lachnospiraceae
lg_intestine_lumen_unclassified_Lachnospiraceae <- lg_intestine_lumen_genus_abund_melt[which(lg_intestine_lumen_genus_abund_melt$Genus=="unclassified Lachnospiraceae"),]
kruskal.test(lg_intestine_lumen_unclassified_Lachnospiraceae$Abundance, lg_intestine_lumen_unclassified_Lachnospiraceae$treatment) #sig.
#Lachnospiraceae NK3A20 group
lg_intestine_lumen_NK3A20_group <- lg_intestine_lumen_genus_abund_melt[which(lg_intestine_lumen_genus_abund_melt$Genus=="Lachnospiraceae_NK3A20_group"),]
kruskal.test(lg_intestine_lumen_NK3A20_group$Abundance, lg_intestine_lumen_NK3A20_group$treatment) #NS
#Blautia
lg_intestine_lumen_Blautia <- lg_intestine_lumen_genus_abund_melt[which(lg_intestine_lumen_genus_abund_melt$Genus=="Blautia"),]
kruskal.test(lg_intestine_lumen_Blautia$Abundance, lg_intestine_lumen_Blautia$treatment) #NS
#Acetitomaculum
lg_intestine_lumen_Acetitomaculum <- lg_intestine_lumen_genus_abund_melt[which(lg_intestine_lumen_genus_abund_melt$Genus=="Acetitomaculum"),]
kruskal.test(lg_intestine_lumen_Acetitomaculum$Abundance, lg_intestine_lumen_Acetitomaculum$treatment) #NS
#Eubacteirum hallii group
lg_intestine_lumen_Eubacterium_hallii_group <- lg_intestine_lumen_genus_abund_melt[which(lg_intestine_lumen_genus_abund_melt$Genus=="[Eubacterium]_hallii_group"),]
kruskal.test(lg_intestine_lumen_Eubacterium_hallii_group$Abundance, lg_intestine_lumen_Eubacterium_hallii_group$treatment) #NS
# Methanobrevibacter
lg_intestine_lumen_Methanobrevibacter <- lg_intestine_lumen_genus_abund_melt[which(lg_intestine_lumen_genus_abund_melt$Genus=="Methanobrevibacter"),]
kruskal.test(lg_intestine_lumen_Methanobrevibacter$Abundance, lg_intestine_lumen_Methanobrevibacter$treatment) #sig
# Monoglobus
lg_intestine_lumen_Monoglobus <- lg_intestine_lumen_genus_abund_melt[which(lg_intestine_lumen_genus_abund_melt$Genus=="Monoglobus"),]
kruskal.test(lg_intestine_lumen_Monoglobus$Abundance, lg_intestine_lumen_Monoglobus$treatment) #sig
# Muribaculaceae
lg_intestine_lumen_Muribaculaceae <- lg_intestine_lumen_genus_abund_melt[which(lg_intestine_lumen_genus_abund_melt$Genus=="Muribaculaceae"),]
kruskal.test(lg_intestine_lumen_Muribaculaceae$Abundance, lg_intestine_lumen_Muribaculaceae$treatment) #sig
# Oscillo. UCG005
lg_intestine_lumen_UCG005 <- lg_intestine_lumen_genus_abund_melt[which(lg_intestine_lumen_genus_abund_melt$Genus=="UCG-005"),]
kruskal.test(lg_intestine_lumen_UCG005$Abundance, lg_intestine_lumen_UCG005$treatment) #NS
# Oscillo. UCG010
lg_intestine_lumen_UCG010 <- lg_intestine_lumen_genus_abund_melt[which(lg_intestine_lumen_genus_abund_melt$Genus=="UCG-010"),]
kruskal.test(lg_intestine_lumen_UCG010$Abundance, lg_intestine_lumen_UCG010$treatment) #NS
# Romboutsia
lg_intestine_lumen_Romboutsia <- lg_intestine_lumen_genus_abund_melt[which(lg_intestine_lumen_genus_abund_melt$Genus=="Romboutsia"),]
kruskal.test(lg_intestine_lumen_Romboutsia$Abundance, lg_intestine_lumen_Romboutsia$treatment) #NS
# Paeniclostridium
lg_intestine_lumen_Paeniclostridium <- lg_intestine_lumen_genus_abund_melt[which(lg_intestine_lumen_genus_abund_melt$Genus=="Paeniclostridium"),]
kruskal.test(lg_intestine_lumen_Paeniclostridium$Abundance, lg_intestine_lumen_Paeniclostridium$treatment) #NS
# unclassified Peptostreptococcaceae
lg_intestine_lumen_unclassified_Peptostreptococcaceae <- lg_intestine_lumen_genus_abund_melt[which(lg_intestine_lumen_genus_abund_melt$Genus=="unclassified Peptostreptococcaceae"),]
kruskal.test(lg_intestine_lumen_unclassified_Peptostreptococcaceae$Abundance, lg_intestine_lumen_unclassified_Peptostreptococcaceae$treatment) #NS
# Prevotellaceae_UCG-003
lg_intestine_lumen_Prevotellaceae_UCG003 <- lg_intestine_lumen_genus_abund_melt[which(lg_intestine_lumen_genus_abund_melt$Genus=="Prevotellaceae_UCG-003"),]
kruskal.test(lg_intestine_lumen_Prevotellaceae_UCG003$Abundance, lg_intestine_lumen_Prevotellaceae_UCG003$treatment) #sig
# Prevotella
lg_intestine_lumen_Prevotella <- lg_intestine_lumen_genus_abund_melt[which(lg_intestine_lumen_genus_abund_melt$Genus=="Prevotella"),]
kruskal.test(lg_intestine_lumen_Prevotella$Abundance, lg_intestine_lumen_Prevotella$treatment) #sig
# Alloprevotella
lg_intestine_lumen_Alloprevotella <- lg_intestine_lumen_genus_abund_melt[which(lg_intestine_lumen_genus_abund_melt$Genus=="Alloprevotella"),]
kruskal.test(lg_intestine_lumen_Alloprevotella$Abundance, lg_intestine_lumen_Alloprevotella$treatment) #sig
# Prevotellaceae_NK3B31_group
lg_intestine_lumen_Prevotellaceae_NK3B31_group <- lg_intestine_lumen_genus_abund_melt[which(lg_intestine_lumen_genus_abund_melt$Genus=="Prevotellaceae_NK3B31_group"),]
kruskal.test(lg_intestine_lumen_Prevotellaceae_NK3B31_group$Abundance, lg_intestine_lumen_Prevotellaceae_NK3B31_group$treatment) #sig
# uncultured_Prevotellaceae
lg_intestine_lumen_uncultured_Prevotellaceae <- lg_intestine_lumen_genus_abund_melt[which(lg_intestine_lumen_genus_abund_melt$Genus=="uncultured"),]
kruskal.test(lg_intestine_lumen_uncultured_Prevotellaceae$Abundance, lg_intestine_lumen_uncultured_Prevotellaceae$treatment) #sig
# Rikenellaceae_RC9_gut_group
lg_intestine_lumen_Rikenellaceae_RC9_gut_group <- lg_intestine_lumen_genus_abund_melt[which(lg_intestine_lumen_genus_abund_melt$Genus=="Rikenellaceae_RC9_gut_group"),]
kruskal.test(lg_intestine_lumen_Rikenellaceae_RC9_gut_group$Abundance, lg_intestine_lumen_Rikenellaceae_RC9_gut_group$treatment) #sig
# Alistipes
lg_intestine_lumen_Alistipes <- lg_intestine_lumen_genus_abund_melt[which(lg_intestine_lumen_genus_abund_melt$Genus=="Alistipes"),]
kruskal.test(lg_intestine_lumen_Alistipes$Abundance, lg_intestine_lumen_Alistipes$treatment) #sig
# Ruminococcus
lg_intestine_lumen_Ruminococcus <- lg_intestine_lumen_genus_abund_melt[which(lg_intestine_lumen_genus_abund_melt$Genus=="Ruminococcus"),]
kruskal.test(lg_intestine_lumen_Ruminococcus$Abundance, lg_intestine_lumen_Ruminococcus$treatment) #sig
# Treponema
lg_intestine_lumen_Treponema <- lg_intestine_lumen_genus_abund_melt[which(lg_intestine_lumen_genus_abund_melt$Genus=="Treponema"),]
kruskal.test(lg_intestine_lumen_Treponema$Abundance, lg_intestine_lumen_Treponema$treatment) #NS
# Parabacteroides
lg_intestine_lumen_Parabacteroides <- lg_intestine_lumen_genus_abund_melt[which(lg_intestine_lumen_genus_abund_melt$Genus=="Parabacteroides"),]
kruskal.test(lg_intestine_lumen_Parabacteroides$Abundance, lg_intestine_lumen_Parabacteroides$treatment) #NS
# low abundance genera (<0.5%)
lg_intestine_lumen_zzzOther <- lg_intestine_lumen_genus_abund_melt[which(lg_intestine_lumen_genus_abund_melt$Genus=="zzzOther "),]
kruskal.test(lg_intestine_lumen_zzzOther$Abundance, lg_intestine_lumen_zzzOther$treatment) #NS


#### EPITHELIUM ####
lg_intestine_epithelium_genus_abund <- merge_low_abundance(lg_intestine_epithelium_genus, threshold = 0.5)
write.csv(otu_table(lg_intestine_epithelium_genus_abund),"top25genus_otus.csv")
write.csv(tax_table(lg_intestine_epithelium_genus_abund),"top25genus_taxa.csv")
lg_intestine_epithelium_genus_abund_melt <- psmelt(lg_intestine_epithelium_genus_abund)
lg_intestine_epithelium_genus_abund_order <- c("[Eubacterium]_coprostanoligenes_group",	"Phascolarctobacterium",	"Akkermansia",	"Family_XIII_AD3011_group",	"Olsenella",	"Bacteroides",	"Bacteroidales_RF16_group",	"Christensenellaceae_R-7_group",	"Clostridium_sensu_stricto_1",	"Turicibacter",	"unclassified Lachnospiraceae",	"[Ruminococcus]_gauvreauii_group",	"Blautia",	"Lachnospiraceae_NK3A20_group",	"Frisingicoccus",	"Marvinbryantia",	"Lachnospiraceae_UCG-010",	"Methanobrevibacter",	"Monoglobus",	"Muribaculaceae",	"UCG-005",	"UCG-010",	"Romboutsia",	"Paeniclostridium",	"unclassified Peptostreptococcaceae",	"Prevotella",	"Alloprevotella",	"Prevotellaceae_UCG-003",	"uncultured",	"Prevotellaceae_NK3B31_group",	"Rikenellaceae_RC9_gut_group",	"Alistipes",	"Ruminococcus",	"unclassified Ruminococcaceae",	"Treponema",	"Sutterella",	"Parabacteroides",	"zzzOther ")

ggplot(lg_intestine_epithelium_genus_abund_melt, aes(x= Genus, y= Abundance)) +
  theme_bw() + labs(y= "Relative Abundance (%)", title = "EPITHELIUM") +
  coord_flip() + 
  scale_x_discrete(limits = rev(lg_intestine_epithelium_genus_abund_order)) +
  scale_y_continuous(expand = c(0.002,0,0.12,0), labels = scales::number_format(accuracy = 0.1)) +
  geom_bar(aes(fill = treatment, colour = treatment), position = position_dodge2(padding = 0.28, reverse = T), stat = "summary", width = 0.72, alpha = 0.5, linewidth = 0.75) +
  geom_errorbar(aes(colour = treatment), position = position_dodge2(padding = 0.55, reverse = T), stat = "summary", width = 0.72, linewidth = 0.75) +
  scale_fill_manual(values = tylosin_palette) +
  scale_colour_manual(values = tylosin_palette) +
  theme(panel.border = element_rect(colour = "black", size = 1),
        plot.title = element_text(size = 34),
        legend.position = "none",
        panel.grid.major.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 16, colour = "black"),
        axis.title.x = element_text(size = 32),
        axis.text.x = element_text(size = 14, colour = "black"),
        axis.ticks = element_line(size = 0.75, colour = "black"))

### stats
#[Eubacterium]_coprostanoligenes_group
lg_intestine_epithelium_Eubacterium_coprostanoligenes_group <- lg_intestine_epithelium_genus_abund_melt[which(lg_intestine_epithelium_genus_abund_melt$Genus=="[Eubacterium]_coprostanoligenes_group"),]
kruskal.test(lg_intestine_epithelium_Eubacterium_coprostanoligenes_group$Abundance, lg_intestine_epithelium_Eubacterium_coprostanoligenes_group$treatment) #NS
#Phascolarctobacterium
lg_intestine_epithelium_Phascolarctobacterium <- lg_intestine_epithelium_genus_abund_melt[which(lg_intestine_epithelium_genus_abund_melt$Genus=="Phascolarctobacterium"),]
kruskal.test(lg_intestine_epithelium_Phascolarctobacterium$Abundance, lg_intestine_epithelium_Phascolarctobacterium$treatment) # sig.
#Akkermansia
lg_intestine_epithelium_Akkermansia <- lg_intestine_epithelium_genus_abund_melt[which(lg_intestine_epithelium_genus_abund_melt$Genus=="Akkermansia"),]
kruskal.test(lg_intestine_epithelium_Akkermansia$Abundance, lg_intestine_epithelium_Akkermansia$treatment) # sig.
#Family_XIII_AD3011_group
lg_intestine_epithelium_Family_XIII_AD3011_group <- lg_intestine_epithelium_genus_abund_melt[which(lg_intestine_epithelium_genus_abund_melt$Genus=="Family_XIII_AD3011_group"),]
kruskal.test(lg_intestine_epithelium_Family_XIII_AD3011_group$Abundance, lg_intestine_epithelium_Family_XIII_AD3011_group$treatment) # sig.
#Olsenella
lg_intestine_epithelium_Olsenella <- lg_intestine_epithelium_genus_abund_melt[which(lg_intestine_epithelium_genus_abund_melt$Genus=="Olsenella"),]
kruskal.test(lg_intestine_epithelium_Olsenella$Abundance, lg_intestine_epithelium_Olsenella$treatment) # sig
# Bacteroides
lg_intestine_epithelium_Bacteroides <- lg_intestine_epithelium_genus_abund_melt[which(lg_intestine_epithelium_genus_abund_melt$Genus=="Bacteroides"),]
kruskal.test(lg_intestine_epithelium_Bacteroides$Abundance, lg_intestine_epithelium_Bacteroides$treatment) # NS
# Bacteroidales_RF16_group
lg_intestine_epithelium_Bacteroidales_RF16_group <- lg_intestine_epithelium_genus_abund_melt[which(lg_intestine_epithelium_genus_abund_melt$Genus=="Bacteroidales_RF16_group"),]
kruskal.test(lg_intestine_epithelium_Bacteroidales_RF16_group$Abundance, lg_intestine_epithelium_Bacteroidales_RF16_group$treatment) # NS
# Christensenellaceae_R-7_group
lg_intestine_epithelium_Christensenellaceae_R7_group <- lg_intestine_epithelium_genus_abund_melt[which(lg_intestine_epithelium_genus_abund_melt$Genus=="Christensenellaceae_R-7_group"),]
kruskal.test(lg_intestine_epithelium_Christensenellaceae_R7_group$Abundance, lg_intestine_epithelium_Christensenellaceae_R7_group$treatment) # NS
# Clostridium_sensu_stricto_1
lg_intestine_epithelium_Clostridium_sensu_stricto_1 <- lg_intestine_epithelium_genus_abund_melt[which(lg_intestine_epithelium_genus_abund_melt$Genus=="Clostridium_sensu_stricto_1"),]
kruskal.test(lg_intestine_epithelium_Clostridium_sensu_stricto_1$Abundance, lg_intestine_epithelium_Clostridium_sensu_stricto_1$treatment) # NS
# Turicibacter
lg_intestine_epithelium_Turicibacter <- lg_intestine_epithelium_genus_abund_melt[which(lg_intestine_epithelium_genus_abund_melt$Genus=="Turicibacter"),]
kruskal.test(lg_intestine_epithelium_Turicibacter$Abundance, lg_intestine_epithelium_Turicibacter$treatment) # NS
# unclassified_Lachnospiraceae
lg_intestine_epithelium_unclassified_Lachnospiraceae <- lg_intestine_epithelium_genus_abund_melt[which(lg_intestine_epithelium_genus_abund_melt$Genus=="unclassified Lachnospiraceae"),]
kruskal.test(lg_intestine_epithelium_unclassified_Lachnospiraceae$Abundance, lg_intestine_epithelium_unclassified_Lachnospiraceae$treatment) #sig.
#ruminococcus_gauvreauii grp.
lg_intestine_epithelium_ruminococcus_gauvreauii <- lg_intestine_epithelium_genus_abund_melt[which(lg_intestine_epithelium_genus_abund_melt$Genus=="[Ruminococcus]_gauvreauii_group"),]
kruskal.test(lg_intestine_epithelium_ruminococcus_gauvreauii$Abundance, lg_intestine_epithelium_ruminococcus_gauvreauii$treatment) # sig.
#Lachnospiraceae NK3A20 group
lg_intestine_epithelium_NK3A20_group <- lg_intestine_epithelium_genus_abund_melt[which(lg_intestine_epithelium_genus_abund_melt$Genus=="Lachnospiraceae_NK3A20_group"),]
kruskal.test(lg_intestine_epithelium_NK3A20_group$Abundance, lg_intestine_epithelium_NK3A20_group$treatment) #NS
#Blautia
lg_intestine_epithelium_Blautia <- lg_intestine_epithelium_genus_abund_melt[which(lg_intestine_epithelium_genus_abund_melt$Genus=="Blautia"),]
kruskal.test(lg_intestine_epithelium_Blautia$Abundance, lg_intestine_epithelium_Blautia$treatment) #NS
#Frisingicoccus
lg_intestine_epithelium_Frisingicoccus <- lg_intestine_epithelium_genus_abund_melt[which(lg_intestine_epithelium_genus_abund_melt$Genus=="Frisingicoccus"),]
kruskal.test(lg_intestine_epithelium_Frisingicoccus$Abundance, lg_intestine_epithelium_Frisingicoccus$treatment) #NS
#Marvinbryantia
lg_intestine_epithelium_Marvinbryantia <- lg_intestine_epithelium_genus_abund_melt[which(lg_intestine_epithelium_genus_abund_melt$Genus=="Marvinbryantia"),]
kruskal.test(lg_intestine_epithelium_Marvinbryantia$Abundance, lg_intestine_epithelium_Marvinbryantia$treatment) #NS
#Lachnospiraceae_UCG-010
lg_intestine_epithelium_Lachnospiraceae_UCG010 <- lg_intestine_epithelium_genus_abund_melt[which(lg_intestine_epithelium_genus_abund_melt$Genus=="Lachnospiraceae_UCG-010"),]
kruskal.test(lg_intestine_epithelium_Lachnospiraceae_UCG010$Abundance, lg_intestine_epithelium_Lachnospiraceae_UCG010$treatment) #NS
# Methanobrevibacter
lg_intestine_epithelium_Methanobrevibacter <- lg_intestine_epithelium_genus_abund_melt[which(lg_intestine_epithelium_genus_abund_melt$Genus=="Methanobrevibacter"),]
kruskal.test(lg_intestine_epithelium_Methanobrevibacter$Abundance, lg_intestine_epithelium_Methanobrevibacter$treatment) #sig
# Monoglobus
lg_intestine_epithelium_Monoglobus <- lg_intestine_epithelium_genus_abund_melt[which(lg_intestine_epithelium_genus_abund_melt$Genus=="Monoglobus"),]
kruskal.test(lg_intestine_epithelium_Monoglobus$Abundance, lg_intestine_epithelium_Monoglobus$treatment) #sig
# Muribaculaceae
lg_intestine_epithelium_Muribaculaceae <- lg_intestine_epithelium_genus_abund_melt[which(lg_intestine_epithelium_genus_abund_melt$Genus=="Muribaculaceae"),]
kruskal.test(lg_intestine_epithelium_Muribaculaceae$Abundance, lg_intestine_epithelium_Muribaculaceae$treatment) #sig
# Oscillo. UCG005
lg_intestine_epithelium_UCG005 <- lg_intestine_epithelium_genus_abund_melt[which(lg_intestine_epithelium_genus_abund_melt$Genus=="UCG-005"),]
kruskal.test(lg_intestine_epithelium_UCG005$Abundance, lg_intestine_epithelium_UCG005$treatment) #NS
# Oscillo. UCG010
lg_intestine_epithelium_UCG010 <- lg_intestine_epithelium_genus_abund_melt[which(lg_intestine_epithelium_genus_abund_melt$Genus=="UCG-010"),]
kruskal.test(lg_intestine_epithelium_UCG010$Abundance, lg_intestine_epithelium_UCG010$treatment) #NS
# Romboutsia
lg_intestine_epithelium_Romboutsia <- lg_intestine_epithelium_genus_abund_melt[which(lg_intestine_epithelium_genus_abund_melt$Genus=="Romboutsia"),]
kruskal.test(lg_intestine_epithelium_Romboutsia$Abundance, lg_intestine_epithelium_Romboutsia$treatment) #NS
# Paeniclostridium
lg_intestine_epithelium_Paeniclostridium <- lg_intestine_epithelium_genus_abund_melt[which(lg_intestine_epithelium_genus_abund_melt$Genus=="Paeniclostridium"),]
kruskal.test(lg_intestine_epithelium_Paeniclostridium$Abundance, lg_intestine_epithelium_Paeniclostridium$treatment) #NS
# unclassified Peptostreptococcaceae
lg_intestine_epithelium_unclassified_Peptostreptococcaceae <- lg_intestine_epithelium_genus_abund_melt[which(lg_intestine_epithelium_genus_abund_melt$Genus=="unclassified Peptostreptococcaceae"),]
kruskal.test(lg_intestine_epithelium_unclassified_Peptostreptococcaceae$Abundance, lg_intestine_epithelium_unclassified_Peptostreptococcaceae$treatment) #NS
# Prevotellaceae_UCG-003
lg_intestine_epithelium_Prevotellaceae_UCG003 <- lg_intestine_epithelium_genus_abund_melt[which(lg_intestine_epithelium_genus_abund_melt$Genus=="Prevotellaceae_UCG-003"),]
kruskal.test(lg_intestine_epithelium_Prevotellaceae_UCG003$Abundance, lg_intestine_epithelium_Prevotellaceae_UCG003$treatment) #sig
# Prevotella
lg_intestine_epithelium_Prevotella <- lg_intestine_epithelium_genus_abund_melt[which(lg_intestine_epithelium_genus_abund_melt$Genus=="Prevotella"),]
kruskal.test(lg_intestine_epithelium_Prevotella$Abundance, lg_intestine_epithelium_Prevotella$treatment) #sig
# Alloprevotella
lg_intestine_epithelium_Alloprevotella <- lg_intestine_epithelium_genus_abund_melt[which(lg_intestine_epithelium_genus_abund_melt$Genus=="Alloprevotella"),]
kruskal.test(lg_intestine_epithelium_Alloprevotella$Abundance, lg_intestine_epithelium_Alloprevotella$treatment) #sig
# Prevotellaceae_NK3B31_group
lg_intestine_epithelium_Prevotellaceae_NK3B31_group <- lg_intestine_epithelium_genus_abund_melt[which(lg_intestine_epithelium_genus_abund_melt$Genus=="Prevotellaceae_NK3B31_group"),]
kruskal.test(lg_intestine_epithelium_Prevotellaceae_NK3B31_group$Abundance, lg_intestine_epithelium_Prevotellaceae_NK3B31_group$treatment) #sig
# uncultured_Prevotellaceae
lg_intestine_epithelium_uncultured_Prevotellaceae <- lg_intestine_epithelium_genus_abund_melt[which(lg_intestine_epithelium_genus_abund_melt$Genus=="uncultured"),]
kruskal.test(lg_intestine_epithelium_uncultured_Prevotellaceae$Abundance, lg_intestine_epithelium_uncultured_Prevotellaceae$treatment) #sig
# Rikenellaceae_RC9_gut_group
lg_intestine_epithelium_Rikenellaceae_RC9_gut_group <- lg_intestine_epithelium_genus_abund_melt[which(lg_intestine_epithelium_genus_abund_melt$Genus=="Rikenellaceae_RC9_gut_group"),]
kruskal.test(lg_intestine_epithelium_Rikenellaceae_RC9_gut_group$Abundance, lg_intestine_epithelium_Rikenellaceae_RC9_gut_group$treatment) #sig
# Alistipes
lg_intestine_epithelium_Alistipes <- lg_intestine_epithelium_genus_abund_melt[which(lg_intestine_epithelium_genus_abund_melt$Genus=="Alistipes"),]
kruskal.test(lg_intestine_epithelium_Alistipes$Abundance, lg_intestine_epithelium_Alistipes$treatment) #sig
# Ruminococcus
lg_intestine_epithelium_Ruminococcus <- lg_intestine_epithelium_genus_abund_melt[which(lg_intestine_epithelium_genus_abund_melt$Genus=="Ruminococcus"),]
kruskal.test(lg_intestine_epithelium_Ruminococcus$Abundance, lg_intestine_epithelium_Ruminococcus$treatment) #sig
# unclassified Ruminococcaceae
lg_intestine_epithelium_unclassified_Ruminococcaceae <- lg_intestine_epithelium_genus_abund_melt[which(lg_intestine_epithelium_genus_abund_melt$Genus=="unclassified Ruminococcaceae"),]
kruskal.test(lg_intestine_epithelium_unclassified_Ruminococcaceae$Abundance, lg_intestine_epithelium_unclassified_Ruminococcaceae$treatment) #sig
# Treponema
lg_intestine_epithelium_Treponema <- lg_intestine_epithelium_genus_abund_melt[which(lg_intestine_epithelium_genus_abund_melt$Genus=="Treponema"),]
kruskal.test(lg_intestine_epithelium_Treponema$Abundance, lg_intestine_epithelium_Treponema$treatment) #NS
# Parabacteroides
lg_intestine_epithelium_Parabacteroides <- lg_intestine_epithelium_genus_abund_melt[which(lg_intestine_epithelium_genus_abund_melt$Genus=="Parabacteroides"),]
kruskal.test(lg_intestine_epithelium_Parabacteroides$Abundance, lg_intestine_epithelium_Parabacteroides$treatment) #NS
# Sutterella
lg_intestine_epithelium_Sutterella <- lg_intestine_epithelium_genus_abund_melt[which(lg_intestine_epithelium_genus_abund_melt$Genus=="Sutterella"),]
kruskal.test(lg_intestine_epithelium_Sutterella$Abundance, lg_intestine_epithelium_Sutterella$treatment) #NS
# low abundance genera (<0.5%)
lg_intestine_epithelium_zzzOther <- lg_intestine_epithelium_genus_abund_melt[which(lg_intestine_epithelium_genus_abund_melt$Genus=="zzzOther "),]
kruskal.test(lg_intestine_epithelium_zzzOther$Abundance, lg_intestine_epithelium_zzzOther$treatment) #NS


################### MOST ABUNDANT GENERA IN LIV ABS ####
#### LAs ####
write.csv(t(otu_table(mergedLAs_genus)),"top25genus_otus.csv")
write.csv(tax_table(mergedLAs_genus),"top25genus_taxa.csv")
mergedLAs_genus_abund_melt <- psmelt(mergedLAs_genus)
mergedLAs_genus_abund_order <- c("Fusobacterium",	"Bacteroides",	"Trueperella",	"Porphyromonas",	"Parvimonas",	"Helcococcus")

ggplot(mergedLAs_genus_abund_melt, aes(x= Genus, y= Abundance)) +
  theme_bw() + labs(y= "Relative Abundance (%)") +
  coord_flip() + 
  scale_x_discrete(limits = rev(mergedLAs_genus_abund_order)) +
  scale_y_continuous(expand = c(0.002,0,0.113,0), labels = scales::number_format(accuracy = 0.1)) +
  geom_bar(aes(fill = treatment, colour = treatment), position = position_dodge2(padding = 0.28, reverse = T), stat = "summary", width = 0.72, alpha = 0.5, linewidth = 0.75) +
  geom_errorbar(aes(colour = treatment), position = position_dodge2(padding = 0.55, reverse = T), stat = "summary", width = 0.72, linewidth = 0.75) +
  scale_fill_manual(values = tylosin_palette) +
  scale_colour_manual(values = tylosin_palette) +
  theme(panel.border = element_rect(colour = "black", size = 1),
        plot.title = element_text(size = 34),
        legend.position = "none",
        panel.grid.major.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 16, colour = "black"),
        axis.title.x = element_text(size = 32),
        axis.text.x = element_text(size = 14, colour = "black"),
        axis.ticks = element_line(size = 0.75, colour = "black"))

### stats
# Fusobacterium
mergedLAs_fusobacterium <- mergedLAs_genus_abund_melt[which(mergedLAs_genus_abund_melt$Genus=="Fusobacterium"),]
kruskal.test(mergedLAs_fusobacterium$Abundance, mergedLAs_fusobacterium$treatment) # NS
# Bacteroides
mergedLAs_Bacteroides <- mergedLAs_genus_abund_melt[which(mergedLAs_genus_abund_melt$Genus=="Bacteroides"),]
kruskal.test(mergedLAs_Bacteroides$Abundance, mergedLAs_Bacteroides$treatment) # NS
# Trueperella
mergedLAs_Trueperella <- mergedLAs_genus_abund_melt[which(mergedLAs_genus_abund_melt$Genus=="Trueperella"),]
kruskal.test(mergedLAs_Trueperella$Abundance, mergedLAs_Trueperella$treatment) # NS
# Porphyromonas
mergedLAs_Porphyromonas <- mergedLAs_genus_abund_melt[which(mergedLAs_genus_abund_melt$Genus=="Porphyromonas"),]
kruskal.test(mergedLAs_Porphyromonas$Abundance, mergedLAs_Porphyromonas$treatment) # NS
# Parvimonas
mergedLAs_Parvimonas <- mergedLAs_genus_abund_melt[which(mergedLAs_genus_abund_melt$Genus=="Parvimonas"),]
kruskal.test(mergedLAs_Parvimonas$Abundance, mergedLAs_Parvimonas$treatment) # NS
# Streptococcus
mergedLAs_Streptococcus <- mergedLAs_genus_abund_melt[which(mergedLAs_genus_abund_melt$Genus=="Streptococcus"),]
kruskal.test(mergedLAs_Streptococcus$Abundance, mergedLAs_Streptococcus$treatment) # NS
# Helcococcus
mergedLAs_Helcococcus <- mergedLAs_genus_abund_melt[which(mergedLAs_genus_abund_melt$Genus=="Helcococcus"),]
kruskal.test(mergedLAs_Helcococcus$Abundance, mergedLAs_Helcococcus$treatment) # NS

################## DIFFERENCES IN COMMUNITY STRUCTURES IN GIT BETWEEN LA+ and LA- ANIMALS #######

LA_major_clade_palette <- c("springgreen4","olivedrab4","darkorchid4","dodgerblue3")

#### RUMEN EPITH.
rumen_epithelium_plot1 <- ordiplot(rumen_epithelium.ord$points)
rumen_epithelium_siteslong <- sites.long(rumen_epithelium_plot1, rumen_epithelium.css.df)
rumen_epithelium_siteslong

rumen_epithelium_centroids <- envfit(rumen_epithelium.ord ~ rumen_epithelium.css.df$LA_Major_Clade2)
rumen_epithelium_centroids

rumen_epithelium_col <- c("2","3","1","N")
rumen_epithelium_NMDS_col1 <- c(-0.0046,0.0530, -0.0506,0.0074)
rumen_epithelium_NMDS_col2 <- c(-0.0127,-0.0013,0.0538,-0.0113)

rumen_epithelium_centroids.df <- data.frame(rumen_epithelium_col, rumen_epithelium_NMDS_col1, rumen_epithelium_NMDS_col2)
rumen_epithelium_centroids.df

ggplot() + theme_bw() +
  labs(x= "NMDS1", y= "NMDS2") +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_shape_manual(values = c(19,19,19,18)) +
  geom_point(data = rumen_epithelium_siteslong, aes(x=axis1,y=axis2, colour= LA_Major_Clade2, shape = LA_Major_Clade2, alpha = LA_Major_Clade2, size = LA_Major_Clade2)) +
  scale_alpha_manual(values = c(0.5,0.5,0.5,0.5), guide = F) +
  scale_size_manual(values = c(5,5,5,6), guide = "none") +
  stat_ellipse(data = rumen_epithelium_siteslong, geom = "polygon", aes(x=axis1,y=axis2, colour = LA_Major_Clade2, fill = LA_Major_Clade2), alpha = c(0.06), lty=2, level = 0.90, linewidth = 1) +
  geom_point(data = rumen_epithelium_centroids.df, aes(x=rumen_epithelium_NMDS_col1, y=rumen_epithelium_NMDS_col2), fill = LA_major_clade_palette, colour = LA_major_clade_palette, size = 14, shape = c(19,19,19,23)) +
  geom_text(data = rumen_epithelium_centroids.df, aes(x=rumen_epithelium_NMDS_col1, y=rumen_epithelium_NMDS_col2, label = rumen_epithelium_col), colour = "white", size = 9, fontface = "bold") +
  scale_colour_manual(values = LA_major_clade_palette) +
  scale_fill_manual(values = LA_major_clade_palette) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        axis.title = element_text(size = 36),
        axis.text = element_text(size = 14, colour = "black"))

pairwise.adonis2(rumen_epithelium.dist ~ LA_Major_Clade2, rumen_epithelium.css.df, strata = 'treatment', p.adjust.methods = "BH", nperm = 9999) # NS

#### RUMEN LUMEN
rumen_lumen_plot1 <- ordiplot(rumen_lumen.ord$points)
rumen_lumen_siteslong <- sites.long(rumen_lumen_plot1, rumen_lumen.css.df)
rumen_lumen_siteslong

rumen_lumen_centroids <- envfit(rumen_lumen.ord ~ rumen_lumen.css.df$LA_Major_Clade2)
rumen_lumen_centroids

rumen_lumen_col <- c("2","3","1","N")
rumen_lumen_NMDS_col1 <- c(-0.01100297,0.02800255,0.05166175,-0.01799732)
rumen_lumen_NMDS_col2 <- c(-0.033373875,-0.011756927,0.025470254,0.008790603)

rumen_lumen_centroids.df <- data.frame(rumen_lumen_col, rumen_lumen_NMDS_col1, rumen_lumen_NMDS_col2)
rumen_lumen_centroids.df

ggplot() + theme_bw() +
  labs(x= "NMDS1", y= "NMDS2") +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_shape_manual(values = c(19,19,19,18)) +
  geom_point(data = rumen_lumen_siteslong, aes(x=axis1,y=axis2, colour= LA_Major_Clade2, shape = LA_Major_Clade2, alpha = LA_Major_Clade2, size = LA_Major_Clade2)) +
  scale_alpha_manual(values = c(0.5,0.5,0.5,0.5), guide = F) +
  scale_size_manual(values = c(5,5,5,6), guide = "none") +
  stat_ellipse(data = rumen_lumen_siteslong, geom = "polygon", aes(x=axis1,y=axis2, colour = LA_Major_Clade2, fill = LA_Major_Clade2), alpha = c(0.06), lty=2, level = 0.90, linewidth = 1) +
  geom_point(data = rumen_lumen_centroids.df, aes(x=rumen_lumen_NMDS_col1, y=rumen_lumen_NMDS_col2), fill = LA_major_clade_palette, colour = LA_major_clade_palette, size = 14, shape = c(19,19,19,23)) +
  geom_text(data = rumen_lumen_centroids.df, aes(x=rumen_lumen_NMDS_col1, y=rumen_lumen_NMDS_col2, label = rumen_lumen_col), colour = "white", size = 9, fontface = "bold") +
  scale_colour_manual(values = LA_major_clade_palette) +
  scale_fill_manual(values = LA_major_clade_palette) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        axis.title = element_text(size = 36),
        axis.text = element_text(size = 14, colour = "black"))

pairwise.adonis2(rumen_lumen.dist ~ LA_Major_Clade2, rumen_lumen.css.df, strata = 'treatment', p.adjust.methods = "BH", nperm = 9999) # NS

#### SM. INTESTINE EPITHELIUM
SI_epithelium_plot1 <- ordiplot(sm_intestine_epithelium.ord$points)
SI_epithelium_siteslong <- sites.long(SI_epithelium_plot1, sm_intestine_epithelium.css.df)
SI_epithelium_siteslong

sm_intestine_epithelium_centroids <- envfit(sm_intestine_epithelium.ord ~ sm_intestine_epithelium.css.df$LA_Major_Clade2)
sm_intestine_epithelium_centroids

SI_epithelium_col <- c("2","3","1","N")
SI_epithelium_NMDS_col1 <- c(-0.0780,0.1257, 0.1500,-0.0387)
SI_epithelium_NMDS_col2 <- c(-0.0021,-0.0212,0.0657,-0.0167)

SI_epithelium_centroids.df <- data.frame(SI_epithelium_col, SI_epithelium_NMDS_col1, SI_epithelium_NMDS_col2)
SI_epithelium_centroids.df

ggplot() + theme_bw() +
  labs(x= "NMDS1", y= "NMDS2") +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_shape_manual(values = c(19,19,19,18)) +
  geom_point(data = SI_epithelium_siteslong, aes(x=axis1,y=axis2, colour= LA_Major_Clade2, shape = LA_Major_Clade2, alpha = LA_Major_Clade2, size = LA_Major_Clade2)) +
  scale_alpha_manual(values = c(0.5,0.5,0.5,0.5), guide = F) +
  scale_size_manual(values = c(5,5,5,6), guide = "none") +
  stat_ellipse(data = SI_epithelium_siteslong, geom = "polygon", aes(x=axis1,y=axis2, colour = LA_Major_Clade2, fill = LA_Major_Clade2), alpha = c(0.06), lty=2, level = 0.90, linewidth = 1) +
  geom_point(data = SI_epithelium_centroids.df, aes(x=SI_epithelium_NMDS_col1, y=SI_epithelium_NMDS_col2), fill = LA_major_clade_palette, colour = LA_major_clade_palette, size = 14, shape = c(19,19,19,23)) +
  geom_text(data = SI_epithelium_centroids.df, aes(x=SI_epithelium_NMDS_col1, y=SI_epithelium_NMDS_col2, label = SI_epithelium_col), colour = "white", size = 9, fontface = "bold") +
  scale_colour_manual(values = LA_major_clade_palette) +
  scale_fill_manual(values = LA_major_clade_palette) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        axis.title = element_text(size = 36),
        axis.text = element_text(size = 14, colour = "black"))

pairwise.adonis2(sm_intestine_epithelium.dist ~ LA_Major_Clade2, sm_intestine_epithelium.css.df, strata = 'treatment', p.adjust.methods = "BH", nperm = 9999) # NS


#### SM. INTESTINE LUMEN
SI_lumen_plot1 <- ordiplot(sm_intestine_lumen.ord$points)
SI_lumen_siteslong <- sites.long(SI_lumen_plot1, sm_intestine_lumen.css.df)
SI_lumen_siteslong

sm_intestine_lumen_centroids <- envfit(sm_intestine_lumen.ord ~ sm_intestine_lumen.css.df$LA_Major_Clade2)
sm_intestine_lumen_centroids

SI_lumen_col <- c("2","3","1","N")
SI_lumen_NMDS_col1 <- c(0.0014,0.0119, 0.0462,-0.0162)
SI_lumen_NMDS_col2 <- c(-0.0021,-0.0079,0.0167,-0.0012)

SI_lumen_centroids.df <- data.frame(SI_lumen_col, SI_lumen_NMDS_col1, SI_lumen_NMDS_col2)
SI_lumen_centroids.df

ggplot() + theme_bw() +
  labs(x= "NMDS1", y= "NMDS2") +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_shape_manual(values = c(19,19,19,18)) +
  geom_point(data = SI_lumen_siteslong, aes(x=axis1,y=axis2, colour= LA_Major_Clade2, shape = LA_Major_Clade2, alpha = LA_Major_Clade2, size = LA_Major_Clade2)) +
  scale_alpha_manual(values = c(0.5,0.5,0.5,0.5), guide = F) +
  scale_size_manual(values = c(5,5,5,6), guide = "none") +
  stat_ellipse(data = SI_lumen_siteslong, geom = "polygon", aes(x=axis1,y=axis2, colour = LA_Major_Clade2, fill = LA_Major_Clade2), alpha = c(0.06), lty=2, level = 0.90, linewidth = 1) +
  geom_point(data = SI_lumen_centroids.df, aes(x=SI_lumen_NMDS_col1, y=SI_lumen_NMDS_col2), fill = LA_major_clade_palette, colour = LA_major_clade_palette, size = 14, shape = c(19,19,19,23)) +
  geom_text(data = SI_lumen_centroids.df, aes(x=SI_lumen_NMDS_col1, y=SI_lumen_NMDS_col2, label = SI_lumen_col), colour = "white", size = 9, fontface = "bold") +
  scale_colour_manual(values = LA_major_clade_palette) +
  scale_fill_manual(values = LA_major_clade_palette) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        axis.title = element_text(size = 36),
        axis.text = element_text(size = 14, colour = "black"))

pairwise.adonis2(sm_intestine_lumen.dist ~ LA_Major_Clade2, sm_intestine_lumen.css.df, strata = 'treatment', p.adjust.methods = "BH", nperm = 9999) # sig.
sm_intestine_lumen_LA_disper <- betadisper(sm_intestine_lumen.dist, sm_intestine_lumen.css.df$LA_Major_Clade2)
plot(sm_intestine_lumen_LA_disper)
permutest(sm_intestine_lumen_LA_disper, permutations = 9999, pairwise = T)

### small intestine lumen abund genera
ggplot(sm_intestine_lumen_genus_abund_melt, aes(x= Genus, y= Abundance)) +
  theme_bw() + labs(y= "Relative Abundance (%)") + 
  scale_y_continuous(expand = c(0.002,0,0.15,0), labels = scales::number_format(accuracy = 0.1)) +
  scale_x_discrete(limits = sm_intestine_lumen_genus_abund_order) +
  geom_bar(aes(fill = LA_Major_Clade2, colour = LA_Major_Clade2), position = position_dodge2(padding = 0.28, reverse = T), stat = "summary", width = 0.72, alpha = 0.5, linewidth = 0.75) +
  geom_errorbar(aes(colour = LA_Major_Clade2), position = position_dodge2(padding = 0.55, reverse = T), stat = "summary", width = 0.72, linewidth = 0.75) +
  scale_fill_manual(values = LA_major_clade_palette) +
  scale_colour_manual(values = LA_major_clade_palette) +
  theme(panel.border = element_rect(colour = "black", size = 1),
        plot.title = element_text(size = 34),
        legend.position = "none",
        panel.grid.major.x = element_blank(),
        axis.title.x = element_blank(),
        #axis.text.x = element_blank(),
        axis.text.x = element_text(angle = 45),
        axis.title.y = element_text(size = 32),
        axis.text.y = element_text(size = 14, colour = "black"),
        axis.ticks = element_line(size = 0.75, colour = "black"))

### stats
#Family_XIII_AD3011_group
pairwise.wilcox.test(sm_intestine_lumen_Family_XIII_AD3011_group$Abundance, sm_intestine_lumen_Family_XIII_AD3011_group$LA_Major_Clade2) # NS
#Olsenella
pairwise.wilcox.test(sm_intestine_lumen_Olsenella$Abundance, sm_intestine_lumen_Olsenella$LA_Major_Clade2) # NS
# Aeriscardovia
pairwise.wilcox.test(sm_intestine_lumen_Aeriscardovia$Abundance, sm_intestine_lumen_Aeriscardovia$LA_Major_Clade2) # NS
# Bifidobacterium
pairwise.wilcox.test(sm_intestine_lumen_Bifidobacterium$Abundance, sm_intestine_lumen_Bifidobacterium$LA_Major_Clade2) # NS
# Clostridium_sensu_stricto_1
pairwise.wilcox.test(sm_intestine_lumen_Clostridium_sensu_stricto_1$Abundance, sm_intestine_lumen_Clostridium_sensu_stricto_1$LA_Major_Clade2) # sig. between no LA and Fuso
# Escherichia-Shigella
pairwise.wilcox.test(sm_intestine_lumen_EscherichiaShigella$Abundance, sm_intestine_lumen_EscherichiaShigella$LA_Major_Clade2) # NS
# Turicibacter
pairwise.wilcox.test(sm_intestine_lumen_Turicibacter$Abundance, sm_intestine_lumen_Turicibacter$LA_Major_Clade2) # NS
#ruminococcus_gauvreauii grp.
pairwise.wilcox.test(sm_intestine_lumen_ruminococcus_gauvreauii$Abundance, sm_intestine_lumen_ruminococcus_gauvreauii$LA_Major_Clade2, p.adjust.method = "BH") # NS
#Acetitomaculum
pairwise.wilcox.test(sm_intestine_lumen_Acetitomaculum$Abundance, sm_intestine_lumen_Acetitomaculum$LA_Major_Clade2) #NS
#Eubacteirum hallii group
pairwise.wilcox.test(sm_intestine_lumen_Eubacterium_hallii_group$Abundance, sm_intestine_lumen_Eubacterium_hallii_group$LA_Major_Clade2) #NS
#Lachnospiraceae NK3A20 group
pairwise.wilcox.test(sm_intestine_lumen_NK3A20_group$Abundance, sm_intestine_lumen_NK3A20_group$LA_Major_Clade2) #NS
#Syntrophococcus
pairwise.wilcox.test(sm_intestine_lumen_Syntrophococcus$Abundance, sm_intestine_lumen_Syntrophococcus$LA_Major_Clade2) #NS
# unclassified_Lachnospiraceae
pairwise.wilcox.test(sm_intestine_lumen_unclassified_Lachnospiraceae$Abundance, sm_intestine_lumen_unclassified_Lachnospiraceae$LA_Major_Clade2) #sig.
# Lachnobacteriuma
pairwise.wilcox.test(sm_intestine_lumen_Lachnobacterium$Abundance, sm_intestine_lumen_Lachnobacterium$LA_Major_Clade2) #NS
# Howardella
pairwise.wilcox.test(sm_intestine_lumen_Howardella$Abundance, sm_intestine_lumen_Howardella$LA_Major_Clade2) #NS
# Methanobrevibacter
pairwise.wilcox.test(sm_intestine_lumen_Methanobrevibacter$Abundance, sm_intestine_lumen_Methanobrevibacter$LA_Major_Clade2) #NS
# Romboutsia
pairwise.wilcox.test(sm_intestine_lumen_Romboutsia$Abundance, sm_intestine_lumen_Romboutsia$LA_Major_Clade2) #NS
# Paeniclostridium
pairwise.wilcox.test(sm_intestine_lumen_Paeniclostridium$Abundance, sm_intestine_lumen_Paeniclostridium$LA_Major_Clade2) #NS
# unclassified Peptostreptococcaceae
pairwise.wilcox.test(sm_intestine_lumen_unclassified_Peptostreptococcaceae$Abundance, sm_intestine_lumen_unclassified_Peptostreptococcaceae$LA_Major_Clade2) #NS
# Ruminococcus
pairwise.wilcox.test(sm_intestine_lumen_Ruminococcus$Abundance, sm_intestine_lumen_Ruminococcus$LA_Major_Clade2) #sig
# unclassified Ruminococcaceae
pairwise.wilcox.test(sm_intestine_lumen_unclassified_Ruminococcaceae$Abundance, sm_intestine_lumen_unclassified_Ruminococcaceae$LA_Major_Clade2) #NS
# low abundance genera (<0.5%)
pairwise.wilcox.test(sm_intestine_lumen_zzzOther$Abundance, sm_intestine_lumen_zzzOther$LA_Major_Clade2) #NS

#### LARGE INTESTINE EPITHELIUM
LI_epithelium_plot1 <- ordiplot(lg_intestine_epithelium.ord$points)
LI_epithelium_siteslong <- sites.long(LI_epithelium_plot1, lg_intestine_epithelium.css.df)
LI_epithelium_siteslong

lg_intestine_epithelium_centroids <- envfit(lg_intestine_epithelium.ord ~ lg_intestine_epithelium.css.df$LA_Major_Clade2)
lg_intestine_epithelium_centroids

LI_epithelium_col <- c("2","3","1","N")
LI_epithelium_NMDS_col1 <- c(-0.0181,0.0209, 0.0163,-0.0012)
LI_epithelium_NMDS_col2 <- c(0.0418,-0.0075,0.0008,-0.0220)

LI_epithelium_centroids.df <- data.frame(LI_epithelium_col, LI_epithelium_NMDS_col1, LI_epithelium_NMDS_col2)
LI_epithelium_centroids.df

ggplot() + theme_bw() +
  labs(x= "NMDS1", y= "NMDS2") +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_shape_manual(values = c(19,19,19,18)) +
  geom_point(data = LI_epithelium_siteslong, aes(x=axis1,y=axis2, colour= LA_Major_Clade2, shape = LA_Major_Clade2, alpha = LA_Major_Clade2, size = LA_Major_Clade2)) +
  scale_alpha_manual(values = c(0.5,0.5,0.5,0.5), guide = F) +
  scale_size_manual(values = c(5,5,5,6), guide = "none") +
  stat_ellipse(data = LI_epithelium_siteslong, geom = "polygon", aes(x=axis1,y=axis2, colour = LA_Major_Clade2, fill = LA_Major_Clade2), alpha = c(0.06), lty=2, level = 0.90, linewidth = 1) +
  geom_point(data = LI_epithelium_centroids.df, aes(x=LI_epithelium_NMDS_col1, y=LI_epithelium_NMDS_col2), fill = LA_major_clade_palette, colour = LA_major_clade_palette, size = 14, shape = c(19,19,19,23)) +
  geom_text(data = LI_epithelium_centroids.df, aes(x=LI_epithelium_NMDS_col1, y=LI_epithelium_NMDS_col2, label = LI_epithelium_col), colour = "white", size = 9, fontface = "bold") +
  scale_colour_manual(values = LA_major_clade_palette) +
  scale_fill_manual(values = LA_major_clade_palette) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        axis.title = element_text(size = 36),
        axis.text = element_text(size = 14, colour = "black"))

pairwise.adonis2(lg_intestine_epithelium.dist ~ LA_Major_Clade2, lg_intestine_epithelium.css.df, strata = 'treatment', p.adjust.methods = "BH", nperm = 9999) # NS


#### LARGE INTESTINE LUMEN
LI_lumen_plot1 <- ordiplot(lg_intestine_lumen.ord$points)
LI_lumen_siteslong <- sites.long(LI_lumen_plot1, lg_intestine_lumen.css.df)
LI_lumen_siteslong

lg_intestine_lumen_centroids <- envfit(lg_intestine_lumen.ord ~ lg_intestine_lumen.css.df$LA_Major_Clade2)
lg_intestine_lumen_centroids

LI_lumen_col <- c("2","3","1","N")
LI_lumen_NMDS_col1 <- c(-0.0254,-0.0269, 0.0622,-0.0022)
LI_lumen_NMDS_col2 <- c(0.0564,-0.0151,-0.0225,-0.0184)

LI_lumen_centroids.df <- data.frame(LI_lumen_col, LI_lumen_NMDS_col1, LI_lumen_NMDS_col2)
LI_lumen_centroids.df

ggplot() + theme_bw() +
  labs(x= "NMDS1", y= "NMDS2") +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_shape_manual(values = c(19,19,19,18)) +
  geom_point(data = LI_lumen_siteslong, aes(x=axis1,y=axis2, colour= LA_Major_Clade2, shape = LA_Major_Clade2, alpha = LA_Major_Clade2, size = LA_Major_Clade2)) +
  scale_alpha_manual(values = c(0.5,0.5,0.5,0.5), guide = F) +
  scale_size_manual(values = c(5,5,5,6), guide = "none") +
  stat_ellipse(data = LI_lumen_siteslong, geom = "polygon", aes(x=axis1,y=axis2, colour = LA_Major_Clade2, fill = LA_Major_Clade2), alpha = c(0.06), lty=2, level = 0.90, linewidth = 1) +
  geom_point(data = LI_lumen_centroids.df, aes(x=LI_lumen_NMDS_col1, y=LI_lumen_NMDS_col2), fill = LA_major_clade_palette, colour = LA_major_clade_palette, size = 14, shape = c(19,19,19,23)) +
  geom_text(data = LI_lumen_centroids.df, aes(x=LI_lumen_NMDS_col1, y=LI_lumen_NMDS_col2, label = LI_lumen_col), colour = "white", size = 9, fontface = "bold") +
  scale_colour_manual(values = LA_major_clade_palette) +
  scale_fill_manual(values = LA_major_clade_palette) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        axis.title = element_text(size = 36),
        axis.text = element_text(size = 14, colour = "black"))

pairwise.adonis2(lg_intestine_lumen.dist ~ LA_Major_Clade2, lg_intestine_lumen.css.df, strata = 'treatment', p.adjust.methods = "BH", nperm = 9999) # NS

#### HIERARCHAL CLUSTERING + RA FOR ALL LAs ####
liv_abs.hclust <- hclust(liv_abs.gunifrac, method = "ward.D2")
plot(liv_abs.hclust)

# convert to dendrogram
liv_abs.dendro <- as.dendrogram(liv_abs.hclust)

# add metadata to dendrogram labels
liv_abs.dendro.data <- dendro_data(liv_abs.dendro, type = "rectangle")
metadata_for_liv_abs_dendro <- as_tibble(liv_abs.css@sam_data)
liv_abs.dendro.data$labels <- liv_abs.dendro.data$labels %>%
  left_join(metadata_for_liv_abs_dendro, by = c("label" = "Sample"))
liv_abs.dendro.data$labels

ggplot(liv_abs.dendro.data$segments) +
  theme_minimal() +
  labs(y= "Ward's Distance") +
  geom_segment(aes(x=x,y=y,xend=xend,yend=yend)) +
  geom_point(data = liv_abs.dendro.data$labels, aes(x=x,y=y, colour = LA_type, fill = LA_type), 
             pch = 22, size = 12, stroke = 1, position = position_nudge(y=-0.03)) +
  geom_text(data = liv_abs.dendro.data$labels, aes(x=x,y=y, label = LA_AnimalID, colour = LA_type),
            size = 4.4, position = position_nudge(y=-0.03), fontface = "bold") +
  scale_colour_manual(values = LA_type_palette) +
  scale_fill_manual(values = alpha(LA_type_palette, 0.45)) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_line(size = 0.75, colour = "black"),
        axis.title.y = element_text(size = 24),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_line(size = 0.7, colour ="black"))

# relative abundance plot for under dendro
# make object for sample order
liv_abs_dendro_order <- liv_abs.dendro.data$labels$label

# let's look at family, and filter super low abundance
LAs_family #31 taxa
LAs_family_filt <- merge_low_abundance(LAs_family, threshold = 0.1)
LAs_family_filt #11 families + other
LAs_family_filt_melt <- psmelt(LAs_family_filt)
LA_family_filt_palette <- c("goldenrod2",	"#E7C57D",	"springgreen4",	"darkorange2",	"#4DA176",	"darkorchid4",	"springgreen3",	"#F3B2E3",	"coral2",	"olivedrab4",	"#EA3C5D",	"grey87")

ggplot(LAs_family_filt_melt, aes(x= sample_Sample, y= Abundance, fill = Family)) +
  theme_minimal() +
  labs(y= "Relative Abundance (%)") +
  geom_bar(stat = "identity", colour = "black") +
  scale_x_discrete(limits = liv_abs_dendro_order, expand = c(0.035,0,0,0)) +
  scale_fill_manual(values = LA_family_filt_palette) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_line(size = 0.75, colour = "black"),
        axis.title.y = element_text(size = 24),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_line(colour = "black", size = 0.7))

#### MULTIPLE LAs RA PLOT for EACH ANIMAL ####
multiple_LAs_family_filt_melt <- LAs_family_filt_melt[which(LAs_family_filt_melt$multiple_LAs=="yes"),]

# plot each on own then combine
GY01_multLAs <- multiple_LAs_family_filt_melt[which(multiple_LAs_family_filt_melt$Animal=="GY01"),]
GY07_multLAs <- multiple_LAs_family_filt_melt[which(multiple_LAs_family_filt_melt$Animal=="GY07"),]
GY08_multLAs <- multiple_LAs_family_filt_melt[which(multiple_LAs_family_filt_melt$Animal=="GY08"),]
GY09_multLAs <- multiple_LAs_family_filt_melt[which(multiple_LAs_family_filt_melt$Animal=="GY09"),]
GY10_multLAs <- multiple_LAs_family_filt_melt[which(multiple_LAs_family_filt_melt$Animal=="GY10"),]
RY01_multLAs <- multiple_LAs_family_filt_melt[which(multiple_LAs_family_filt_melt$Animal=="RY01"),]
RY02_multLAs <- multiple_LAs_family_filt_melt[which(multiple_LAs_family_filt_melt$Animal=="RY02"),]
RY03_multLAs <- multiple_LAs_family_filt_melt[which(multiple_LAs_family_filt_melt$Animal=="RY03"),]
RY04_multLAs <- multiple_LAs_family_filt_melt[which(multiple_LAs_family_filt_melt$Animal=="RY04"),]
RY05_multLAs <- multiple_LAs_family_filt_melt[which(multiple_LAs_family_filt_melt$Animal=="RY05"),]
RY06_multLAs <- multiple_LAs_family_filt_melt[which(multiple_LAs_family_filt_melt$Animal=="RY06"),]
RY07_multLAs <- multiple_LAs_family_filt_melt[which(multiple_LAs_family_filt_melt$Animal=="RY07"),]
RY08_multLAs <- multiple_LAs_family_filt_melt[which(multiple_LAs_family_filt_melt$Animal=="RY08"),]
RY09_multLAs <- multiple_LAs_family_filt_melt[which(multiple_LAs_family_filt_melt$Animal=="RY09"),]
RY10_multLAs <- multiple_LAs_family_filt_melt[which(multiple_LAs_family_filt_melt$Animal=="RY10"),]

# plots
ggplot(GY01_multLAs, aes(x = Sample, y= Abundance, fill = Family)) +
  theme_bw() + labs(y="Relative Abundance (%)") + facet_grid(~Animal) +
  geom_bar(stat = "summary", colour = "black") +
  scale_y_continuous(expand = c(0,0.1,0.001,0)) +
  scale_fill_manual(values = LA_family_filt_palette) +
  theme(legend.position = "none",
        strip.background = element_rect(fill="black", colour = "black", linewidth = 1),
        strip.text = element_text(size = 20, colour = "black"),
        panel.border = element_rect(colour = "black", linewidth = 1),
        axis.title.y = element_text(size = 32),
        axis.text.y = element_text(colour = "black", size = 16),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

ggplot(GY07_multLAs, aes(x = Sample, y= Abundance, fill = Family)) +
  theme_bw() + labs(y="Relative Abundance (%)") + facet_grid(~Animal) +
  geom_bar(stat = "summary", colour = "black") +
  scale_y_continuous(expand = c(0,0.1,0.001,0)) +
  scale_fill_manual(values = LA_family_filt_palette) +
  theme(legend.position = "none",
        strip.background = element_rect(fill="black", colour = "black", linewidth = 1),
        strip.text = element_text(size = 20, colour = "black"),
        panel.border = element_rect(colour = "black", linewidth = 1),
        axis.title.y = element_text(size = 32),
        axis.text.y = element_text(colour = "black", size = 16),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

ggplot(GY08_multLAs, aes(x = Sample, y= Abundance, fill = Family)) +
  theme_bw() + labs(y="Relative Abundance (%)") + facet_grid(~Animal) +
  geom_bar(stat = "summary", colour = "black") +
  scale_y_continuous(expand = c(0,0.1,0.001,0)) +
  scale_fill_manual(values = LA_family_filt_palette) +
  theme(legend.position = "none",
        strip.background = element_rect(fill="black", colour = "black", linewidth = 1),
        strip.text = element_text(size = 20, colour = "black"),
        panel.border = element_rect(colour = "black", linewidth = 1),
        axis.title.y = element_text(size = 32),
        axis.text.y = element_text(colour = "black", size = 16),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

ggplot(GY10_multLAs, aes(x = Sample, y= Abundance, fill = Family)) +
  theme_bw() + labs(y="Relative Abundance (%)") + facet_grid(~Animal) +
  geom_bar(stat = "summary", colour = "black") +
  scale_y_continuous(expand = c(0,0.1,0.001,0)) +
  scale_fill_manual(values = LA_family_filt_palette) +
  theme(legend.position = "none",
        strip.background = element_rect(fill="black", colour = "black", linewidth = 1),
        strip.text = element_text(size = 20, colour = "black"),
        panel.border = element_rect(colour = "black", linewidth = 1),
        axis.title.y = element_text(size = 32),
        axis.text.y = element_text(colour = "black", size = 16),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

ggplot(RY01_multLAs, aes(x = Sample, y= Abundance, fill = Family)) +
  theme_bw() + labs(y="Relative Abundance (%)") + facet_grid(~Animal) +
  geom_bar(stat = "summary", colour = "black") +
  scale_y_continuous(expand = c(0,0.1,0.001,0)) +
  scale_fill_manual(values = LA_family_filt_palette) +
  theme(legend.position = "none",
        strip.background = element_rect(fill="black", colour = "black", linewidth = 1),
        strip.text = element_text(size = 20, colour = "black"),
        panel.border = element_rect(colour = "black", linewidth = 1),
        axis.title.y = element_text(size = 32),
        axis.text.y = element_text(colour = "black", size = 16),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

ggplot(RY02_multLAs, aes(x = Sample, y= Abundance, fill = Family)) +
  theme_bw() + labs(y="Relative Abundance (%)") + facet_grid(~Animal) +
  geom_bar(stat = "summary", colour = "black") +
  scale_y_continuous(expand = c(0,0.1,0.001,0)) +
  scale_fill_manual(values = LA_family_filt_palette) +
  theme(legend.position = "none",
        strip.background = element_rect(fill="black", colour = "black", linewidth = 1),
        strip.text = element_text(size = 20, colour = "black"),
        panel.border = element_rect(colour = "black", linewidth = 1),
        axis.title.y = element_text(size = 32),
        axis.text.y = element_text(colour = "black", size = 16),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

ggplot(RY03_multLAs, aes(x = Sample, y= Abundance, fill = Family)) +
  theme_bw() + labs(y="Relative Abundance (%)") + facet_grid(~Animal) +
  geom_bar(stat = "summary", colour = "black") +
  scale_y_continuous(expand = c(0,0.1,0.001,0)) +
  scale_fill_manual(values = LA_family_filt_palette) +
  theme(legend.position = "none",
        strip.background = element_rect(fill="black", colour = "black", linewidth = 1),
        strip.text = element_text(size = 20, colour = "black"),
        panel.border = element_rect(colour = "black", linewidth = 1),
        axis.title.y = element_text(size = 32),
        axis.text.y = element_text(colour = "black", size = 16),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

ggplot(RY04_multLAs, aes(x = Sample, y= Abundance, fill = Family)) +
  theme_bw() + labs(y="Relative Abundance (%)") + facet_grid(~Animal) +
  geom_bar(stat = "summary", colour = "black") +
  scale_y_continuous(expand = c(0,0.1,0.001,0)) +
  scale_fill_manual(values = LA_family_filt_palette) +
  theme(legend.position = "none",
        strip.background = element_rect(fill="black", colour = "black", linewidth = 1),
        strip.text = element_text(size = 20, colour = "black"),
        panel.border = element_rect(colour = "black", linewidth = 1),
        axis.title.y = element_text(size = 32),
        axis.text.y = element_text(colour = "black", size = 16),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

ggplot(RY05_multLAs, aes(x = Sample, y= Abundance, fill = Family)) +
  theme_bw() + labs(y="Relative Abundance (%)") + facet_grid(~Animal) +
  geom_bar(stat = "summary", colour = "black") +
  scale_y_continuous(expand = c(0,0.1,0.001,0)) +
  scale_fill_manual(values = LA_family_filt_palette) +
  theme(legend.position = "none",
        strip.background = element_rect(fill="black", colour = "black", linewidth = 1),
        strip.text = element_text(size = 20, colour = "black"),
        panel.border = element_rect(colour = "black", linewidth = 1),
        axis.title.y = element_text(size = 32),
        axis.text.y = element_text(colour = "black", size = 16),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

ggplot(RY06_multLAs, aes(x = Sample, y= Abundance, fill = Family)) +
  theme_bw() + labs(y="Relative Abundance (%)") + facet_grid(~Animal) +
  geom_bar(stat = "summary", colour = "black") +
  scale_y_continuous(expand = c(0,0.1,0.001,0)) +
  scale_fill_manual(values = LA_family_filt_palette) +
  theme(legend.position = "none",
        strip.background = element_rect(fill="black", colour = "black", linewidth = 1),
        strip.text = element_text(size = 20, colour = "black"),
        panel.border = element_rect(colour = "black", linewidth = 1),
        axis.title.y = element_text(size = 32),
        axis.text.y = element_text(colour = "black", size = 16),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

ggplot(RY07_multLAs, aes(x = Sample, y= Abundance, fill = Family)) +
  theme_bw() + labs(y="Relative Abundance (%)") + facet_grid(~Animal) +
  geom_bar(stat = "summary", colour = "black") +
  scale_y_continuous(expand = c(0,0.1,0.001,0)) +
  scale_fill_manual(values = LA_family_filt_palette) +
  theme(legend.position = "none",
        strip.background = element_rect(fill="black", colour = "black", linewidth = 1),
        strip.text = element_text(size = 20, colour = "black"),
        panel.border = element_rect(colour = "black", linewidth = 1),
        axis.title.y = element_text(size = 32),
        axis.text.y = element_text(colour = "black", size = 16),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

ggplot(RY08_multLAs, aes(x = Sample, y= Abundance, fill = Family)) +
  theme_bw() + labs(y="Relative Abundance (%)") + facet_grid(~Animal) +
  geom_bar(stat = "summary", colour = "black") +
  scale_y_continuous(expand = c(0,0.1,0.001,0)) +
  scale_fill_manual(values = LA_family_filt_palette) +
  theme(legend.position = "none",
        strip.background = element_rect(fill="black", colour = "black", linewidth = 1),
        strip.text = element_text(size = 20, colour = "black"),
        panel.border = element_rect(colour = "black", linewidth = 1),
        axis.title.y = element_text(size = 32),
        axis.text.y = element_text(colour = "black", size = 16),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

ggplot(RY09_multLAs, aes(x = Sample, y= Abundance, fill = Family)) +
  theme_bw() + labs(y="Relative Abundance (%)") + facet_grid(~Animal) +
  geom_bar(stat = "summary", colour = "black") +
  scale_y_continuous(expand = c(0,0.1,0.001,0)) +
  scale_fill_manual(values = LA_family_filt_palette) +
  theme(legend.position = "none",
        strip.background = element_rect(fill="black", colour = "black", linewidth = 1),
        strip.text = element_text(size = 20, colour = "black"),
        panel.border = element_rect(colour = "black", linewidth = 1),
        axis.title.y = element_text(size = 32),
        axis.text.y = element_text(colour = "black", size = 16),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

ggplot(RY10_multLAs, aes(x = Sample, y= Abundance, fill = Family)) +
  theme_bw() + labs(y="Relative Abundance (%)") + facet_grid(~Animal) +
  geom_bar(stat = "summary", colour = "black") +
  scale_y_continuous(expand = c(0,0.1,0.001,0)) +
  scale_fill_manual(values = LA_family_filt_palette) +
  theme(legend.position = "none",
        strip.background = element_rect(fill="black", colour = "black", linewidth = 1),
        strip.text = element_text(size = 20, colour = "black"),
        panel.border = element_rect(colour = "black", linewidth = 1),
        axis.title.y = element_text(size = 32),
        axis.text.y = element_text(colour = "black", size = 16),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

#### PREVALENT TAXA WITHIN EACH OF THE THREE CLADES ####
LAs_genus_prev <- merge_low_abundance(LAs_genus, threshold = 0.0001)
write.csv(tax_table(LAs_genus_lefse), "LA_lefse_taxa.csv")
write.csv(otu_table(LAs_genus_lefse), "LA_lefse_otus.csv")

LAs_family_prev <- merge_low_abundance(LAs_family, threshold = 0.0001)
write.csv(tax_table(LAs_family_lefse),"LAs_lefse_fam_taxa.csv")
write.csv(otu_table(LAs_family_lefse),"LAs_lefse_fam_otus.csv")

#### 11 families in prevalence > 0.4
## Fusobacteriaceae
LA_fusobacteriaceae <- subset_taxa(LAs_family, Family=="Fusobacteriaceae") %>%
  psmelt()

ggplot(LA_fusobacteriaceae, aes(x= LA_Major_Clade, y= Abundance, fill = LA_Major_Clade, colour = LA_Major_Clade)) +
  theme_bw() +
  labs(y= "Fusobacteriaceae (RA%)") +
  geom_boxplot(alpha = 0.5) +
  geom_point(shape = 18, size = 3) +
  scale_y_continuous(limits = c(0,100), expand = c(0.001,0,.1,0)) +
  scale_x_discrete(limits = c("Fuso","Bact-1","Bact-2"), labels = c("CLADE 1","CLADE 2","CLADE 3")) +
  scale_fill_manual(values = LA_clade_palette) +
  scale_colour_manual(values = LA_clade_palette) +
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", linewidth = 1),
        panel.grid.major.x = element_blank(),
        axis.title.y = element_text(size = 38),
        axis.text.y = element_text(colour = "black", size = 16),
        axis.ticks.y = element_line(colour = "black", linewidth = 1),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

pairwise.wilcox.test(LA_fusobacteriaceae$Abundance, LA_fusobacteriaceae$LA_Major_Clade, p.adjust.method = "BH") #all from all

## Bacteroidaceae
LA_bacteroidaceae <- subset_taxa(LAs_family, Family=="Bacteroidaceae") %>%
  psmelt()

ggplot(LA_bacteroidaceae, aes(x= LA_Major_Clade, y= Abundance, fill = LA_Major_Clade, colour = LA_Major_Clade)) +
  theme_bw() +
  labs(y= "Bacteroidaceae (RA%)") +
  geom_boxplot(alpha = 0.5) +
  geom_point(shape = 18, size = 3) +
  scale_y_continuous(limits = c(0,100), expand = c(0.001,0,.1,0)) +
  scale_x_discrete(limits = c("Fuso","Bact-1","Bact-2"), labels = c("CLADE 1","CLADE 2","CLADE 3")) +
  scale_fill_manual(values = LA_clade_palette) +
  scale_colour_manual(values = LA_clade_palette) +
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", linewidth = 1),
        panel.grid.major.x = element_blank(),
        axis.title.y = element_text(size = 38),
        axis.text.y = element_text(colour = "black", size = 16),
        axis.ticks.y = element_line(colour = "black", linewidth = 1),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

pairwise.wilcox.test(LA_bacteroidaceae$Abundance, LA_bacteroidaceae$LA_Major_Clade, p.adjust.method = "BH") #all from all

## Porphyromonadaceae
LA_porphyromonadaceae <- subset_taxa(LAs_family, Family=="Porphyromonadaceae") %>%
  psmelt()

ggplot(LA_porphyromonadaceae, aes(x= LA_Major_Clade, y= Abundance, fill = LA_Major_Clade, colour = LA_Major_Clade)) +
  theme_bw() +
  labs(y= "Porphyromonadaceae (RA%)") +
  geom_boxplot(alpha = 0.5) +
  geom_point(shape = 18, size = 3) +
  scale_y_continuous(limits = c(0,100), expand = c(0.001,0,.1,0)) +
  scale_x_discrete(limits = c("Fuso","Bact-1","Bact-2"), labels = c("CLADE 1","CLADE 2","CLADE 3")) +
  scale_fill_manual(values = LA_clade_palette) +
  scale_colour_manual(values = LA_clade_palette) +
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", linewidth = 1),
        panel.grid.major.x = element_blank(),
        axis.title.y = element_text(size = 38),
        axis.text.y = element_text(colour = "black", size = 16),
        axis.ticks.y = element_line(colour = "black", linewidth = 1),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

pairwise.wilcox.test(LA_porphyromonadaceae$Abundance, LA_porphyromonadaceae$LA_Major_Clade, p.adjust.method = "BH") #all from all

## actinomycetaceae
LA_actinomycetaceae <- subset_taxa(LAs_family, Family=="Actinomycetaceae") %>%
  psmelt()

ggplot(LA_actinomycetaceae, aes(x= LA_Major_Clade, y= Abundance, fill = LA_Major_Clade, colour = LA_Major_Clade)) +
  theme_bw() +
  labs(y= "Actinomycetaceae (RA%)") +
  geom_boxplot(alpha = 0.5) +
  geom_point(shape = 18, size = 3) +
  scale_y_continuous(limits = c(0,100), expand = c(0.001,0,.1,0)) +
  scale_x_discrete(limits = c("Fuso","Bact-1","Bact-2"), labels = c("CLADE 1","CLADE 2","CLADE 3")) +
  scale_fill_manual(values = LA_clade_palette) +
  scale_colour_manual(values = LA_clade_palette) +
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", linewidth = 1),
        panel.grid.major.x = element_blank(),
        axis.title.y = element_text(size = 38),
        axis.text.y = element_text(colour = "black", size = 16),
        axis.ticks.y = element_line(colour = "black", linewidth = 1),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

pairwise.wilcox.test(LA_actinomycetaceae$Abundance, LA_actinomycetaceae$LA_Major_Clade, p.adjust.method = "BH") #all from all

## peptostreptococcales_tissierellales
LA_peptostreptococcales_tissierellales <- subset_taxa(LAs_family, Family=="Peptostreptococcales-Tissierellales") %>%
  psmelt()

ggplot(LA_peptostreptococcales_tissierellales, aes(x= LA_Major_Clade, y= Abundance, fill = LA_Major_Clade, colour = LA_Major_Clade)) +
  theme_bw() +
  labs(y= "Peptostreptococcales-Tissierellales (RA%)") +
  geom_boxplot(alpha = 0.5) +
  geom_point(shape = 18, size = 3) +
  scale_y_continuous(limits = c(0,100), expand = c(0.001,0,.1,0)) +
  scale_x_discrete(limits = c("Fuso","Bact-1","Bact-2"), labels = c("CLADE 1","CLADE 2","CLADE 3")) +
  scale_fill_manual(values = LA_clade_palette) +
  scale_colour_manual(values = LA_clade_palette) +
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", linewidth = 1),
        panel.grid.major.x = element_blank(),
        axis.title.y = element_text(size = 38),
        axis.text.y = element_text(colour = "black", size = 16),
        axis.ticks.y = element_line(colour = "black", linewidth = 1),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

pairwise.wilcox.test(LA_peptostreptococcales_tissierellales$Abundance, LA_peptostreptococcales_tissierellales$LA_Major_Clade, p.adjust.method = "BH") #all from all

## peptostreptococcaceae
LA_peptostreptococcaceae <- subset_taxa(LAs_family, Family=="Peptostreptococcaceae") %>%
  psmelt()

ggplot(LA_peptostreptococcaceae, aes(x= LA_Major_Clade, y= Abundance, fill = LA_Major_Clade, colour = LA_Major_Clade)) +
  theme_bw() +
  labs(y= "Peptostreptococcaceae (RA%)") +
  geom_boxplot(alpha = 0.5) +
  geom_point(shape = 18, size = 3) +
  scale_y_continuous(limits = c(0,10), expand = c(0.001,0,.1,0)) +
  scale_x_discrete(limits = c("Fuso","Bact-1","Bact-2"), labels = c("CLADE 1","CLADE 2","CLADE 3")) +
  scale_fill_manual(values = LA_clade_palette) +
  scale_colour_manual(values = LA_clade_palette) +
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", linewidth = 1),
        panel.grid.major.x = element_blank(),
        axis.title.y = element_text(size = 38),
        axis.text.y = element_text(colour = "black", size = 16),
        axis.ticks.y = element_line(colour = "black", linewidth = 1),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

pairwise.wilcox.test(LA_peptostreptococcaceae$Abundance, LA_peptostreptococcaceae$LA_Major_Clade, p.adjust.method = "BH") #all from all

## Campylobacteraceae
LA_Campylobacteraceae <- subset_taxa(LAs_family, Family=="Campylobacteraceae") %>%
  psmelt()

ggplot(LA_Campylobacteraceae, aes(x= LA_Major_Clade, y= Abundance, fill = LA_Major_Clade, colour = LA_Major_Clade)) +
  theme_bw() +
  labs(y= "Campylobacteraceae (RA%)") +
  geom_boxplot(alpha = 0.5) +
  geom_point(shape = 18, size = 3) +
  scale_y_continuous(limits = c(0,5), expand = c(0.001,0,.1,0)) +
  scale_x_discrete(limits = c("Fuso","Bact-1","Bact-2"), labels = c("CLADE 1","CLADE 2","CLADE 3")) +
  scale_fill_manual(values = LA_clade_palette) +
  scale_colour_manual(values = LA_clade_palette) +
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", linewidth = 1),
        panel.grid.major.x = element_blank(),
        axis.title.y = element_text(size = 38),
        axis.text.y = element_text(colour = "black", size = 16),
        axis.ticks.y = element_line(colour = "black", linewidth = 1),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

pairwise.wilcox.test(LA_Campylobacteraceae$Abundance, LA_Campylobacteraceae$LA_Major_Clade, p.adjust.method = "BH") #all from all

## Atopobiaceae
LA_Atopobiaceae <- subset_taxa(LAs_family, Family=="Atopobiaceae") %>%
  psmelt()

ggplot(LA_Atopobiaceae, aes(x= LA_Major_Clade, y= Abundance, fill = LA_Major_Clade, colour = LA_Major_Clade)) +
  theme_bw() +
  labs(y= "Atopobiaceae (RA%)") +
  geom_boxplot(alpha = 0.5) +
  geom_point(shape = 18, size = 3) +
  scale_y_continuous(limits = c(0,5), expand = c(0.001,0,.1,0)) +
  scale_x_discrete(limits = c("Fuso","Bact-1","Bact-2"), labels = c("CLADE 1","CLADE 2","CLADE 3")) +
  scale_fill_manual(values = LA_clade_palette) +
  scale_colour_manual(values = LA_clade_palette) +
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", linewidth = 1),
        panel.grid.major.x = element_blank(),
        axis.title.y = element_text(size = 38),
        axis.text.y = element_text(colour = "black", size = 16),
        axis.ticks.y = element_line(colour = "black", linewidth = 1),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

pairwise.wilcox.test(LA_Atopobiaceae$Abundance, LA_Atopobiaceae$LA_Major_Clade, p.adjust.method = "BH") #all from all

## Desulfovibrionaceae
LA_Desulfovibrionaceae <- subset_taxa(LAs_family, Family=="Desulfovibrionaceae") %>%
  psmelt()

ggplot(LA_Desulfovibrionaceae, aes(x= LA_Major_Clade, y= Abundance, fill = LA_Major_Clade, colour = LA_Major_Clade)) +
  theme_bw() +
  labs(y= "Desulfovibrionaceae (RA%)") +
  geom_boxplot(alpha = 0.5) +
  geom_point(shape = 18, size = 3) +
  scale_y_continuous(limits = c(0,5), expand = c(0.001,0,.1,0)) +
  scale_x_discrete(limits = c("Fuso","Bact-1","Bact-2"), labels = c("CLADE 1","CLADE 2","CLADE 3")) +
  scale_fill_manual(values = LA_clade_palette) +
  scale_colour_manual(values = LA_clade_palette) +
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", linewidth = 1),
        panel.grid.major.x = element_blank(),
        axis.title.y = element_text(size = 38),
        axis.text.y = element_text(colour = "black", size = 16),
        axis.ticks.y = element_line(colour = "black", linewidth = 1),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

pairwise.wilcox.test(LA_Desulfovibrionaceae$Abundance, LA_Desulfovibrionaceae$LA_Major_Clade, p.adjust.method = "BH") #all from all

## Lachnospiraceae
LA_Lachnospiraceae <- subset_taxa(LAs_family, Family=="Lachnospiraceae") %>%
  psmelt()

ggplot(LA_Lachnospiraceae, aes(x= LA_Major_Clade, y= Abundance, fill = LA_Major_Clade, colour = LA_Major_Clade)) +
  theme_bw() +
  labs(y= "Lachnospiraceae (RA%)") +
  geom_boxplot(alpha = 0.5) +
  geom_point(shape = 18, size = 3) +
  scale_y_continuous(limits = c(0,5), expand = c(0.001,0,.1,0)) +
  scale_x_discrete(limits = c("Fuso","Bact-1","Bact-2"), labels = c("CLADE 1","CLADE 2","CLADE 3")) +
  scale_fill_manual(values = LA_clade_palette) +
  scale_colour_manual(values = LA_clade_palette) +
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", linewidth = 1),
        panel.grid.major.x = element_blank(),
        axis.title.y = element_text(size = 38),
        axis.text.y = element_text(colour = "black", size = 16),
        axis.ticks.y = element_line(colour = "black", linewidth = 1),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

pairwise.wilcox.test(LA_Lachnospiraceae$Abundance, LA_Lachnospiraceae$LA_Major_Clade, p.adjust.method = "BH") #all from all

ggplot(gut_lachno, aes(x= body_location, y= Abundance, fill= body_location)) +
  facet_grid(~matrix) +
  geom_bar(stat = "summary") + geom_errorbar(stat = "summary")

pairwise.wilcox.test(gut_lachno$Abundance, gut_lachno$body_location, p.adjust.method = "BH")


### FAMILIES DIF ABUND IN:
## CLADE 1: Fusobacteriaceae
## CLADE 2: Bacteroidaceae
## CLADE 3: Actinomycetaceae, Atopobiaceae, Porphyromonadaceae, Peptostreptococcales-Tissierellales, Peptostreptococcaceae, Lachnospiraceae, Desulfovibrionaceae 

## look at these in the gut
gut_samples_family_LA <- tax_glom(gut_samples_ra, taxrank = "Family", NArm = F)
gut_samples_family_LA <- prune_taxa(taxa_sums(gut_samples_family_LA) > 0, gut_samples_family_LA)
gut_samples_family_LA

### gathering discriminant taxa
## CLADE 1
gut_clade1_disfamily <- subset_taxa(gut_samples_family_LA, Family=="Fusobacteriaceae")
gut_clade1_disfamily
gut_clade1_disfamily_melt <- psmelt(gut_clade1_disfamily)

## CLADE2
gut_clade2_disfamily  <- subset_taxa(gut_samples_family_LA, Family=="Bacteroidaceae")
gut_clade2_disfamily
gut_clade2_disfamily_melt <- psmelt(gut_clade2_disfamily)

## CLADE 3
gut_clade3_disfamily <- subset_taxa(gut_samples_family_LA, Family=="Actinomycetaceae" |
                                      Family=="Porphyromonadaceae" |
                                      Family=="Peptostreptococcales-Tissierellales" |
                                      Family=="Peptostreptococcaceae" |
                                      Family=="Desulfovibrionaceae" |
                                      Family=="Lachnospiraceae")
gut_clade3_disfamily
gut_clade3_disfamily_melt <- psmelt(gut_clade3_disfamily)

#### SPLITTING BY LUMEN/EPI AND PLOTTING/STATS

## CLADE 1
gut_clade1_fam_lumen <- gut_clade1_disfamily_melt[which(gut_clade1_disfamily_melt$matrix=="lumen"),]
gut_clade1_fam_epith <- gut_clade1_disfamily_melt[which(gut_clade1_disfamily_melt$matrix=="epithelium"),]

ggplot(gut_clade1_fam_epith, aes(x= body_location, y= Abundance, fill = Family, colour = Family)) +
  theme_bw() +
  labs(title = "CLADE 1", y= "Fusobacteriaceae (RA %)") +
  scale_y_continuous(expand = c(0.001,0,0.1,0)) +
  scale_x_discrete(limits = c("rumen","small intestine","large intestine"), labels = c("RUMEN","ILEUM","COLON")) +
  scale_fill_manual(values = "darkorchid4") +
  scale_colour_manual(values = "darkorchid4") +
  geom_boxplot(alpha=0.5, linewidth = 0.8) +
  geom_point() +
  theme(legend.position = "none",
        plot.title = element_text(size = 32),
        panel.border = element_rect(colour = "black", linewidth = 1),
        panel.grid.major.x = element_blank(),
        axis.title.y = element_text(size = 32),
        axis.text.y = element_text(colour = "black", size = 16),
        axis.ticks = element_line(colour = "black", linewidth = 1),
        axis.title.x = element_blank(),
        axis.text.x = element_text(colour = "black", size = 28, angle = 45, hjust=0.95, vjust=0.95))

ggplot(gut_clade1_fam_lumen, aes(x= body_location, y= Abundance, fill = Family, colour = Family)) +
  theme_bw() +
  labs(title = "CLADE 1", y= "Fusobacteriaceae (RA %)") +
  scale_y_continuous(expand = c(0.001,0,0.1,0), labels = scales::number_format(accuracy = 0.01)) +
  scale_x_discrete(limits = c("rumen","small intestine","large intestine"), labels = c("RUMEN","ILEUM","COLON")) +
  scale_fill_manual(values = "darkorchid4") +
  scale_colour_manual(values = "darkorchid4") +
  geom_boxplot(alpha=0.5, linewidth = 0.8) +
  geom_point() +
  theme(legend.position = "none",
        plot.title = element_text(size = 32),
        panel.border = element_rect(colour = "black", linewidth = 1),
        panel.grid.major.x = element_blank(),
        axis.title.y = element_text(size = 32),
        axis.text.y = element_text(colour = "black", size = 16),
        axis.ticks = element_line(colour = "black", linewidth = 1),
        axis.title.x = element_blank(),
        axis.text.x = element_text(colour = "black", size = 28, angle = 45, hjust=0.95, vjust=0.95))


pairwise.wilcox.test(gut_clade1_fam_lumen$Abundance, gut_clade1_fam_lumen$body_location, p.adjust.method = "BH") #NS
pairwise.wilcox.test(gut_clade1_fam_epith$Abundance, gut_clade1_fam_epith$body_location, p.adjust.method = "BH") #rumen higher than LI

#### CLADE 2
gut_clade2_fam_lumen <- gut_clade2_disfamily_melt[which(gut_clade2_disfamily_melt$matrix=="lumen"),]
gut_clade2_fam_epith <- gut_clade2_disfamily_melt[which(gut_clade2_disfamily_melt$matrix=="epithelium"),]

ggplot(gut_clade2_fam_lumen, aes(x= body_location, y= Abundance, fill = Family, colour = Family)) +
  theme_bw() +
  labs(title = "CLADE 2", y= "Bacteroidaceae (RA %)") +
  scale_y_continuous(expand = c(0.001,0,0.1,0), labels = scales::number_format(accuracy = 0.1)) +
  scale_x_discrete(limits = c("rumen","small intestine","large intestine"), labels = c("RUMEN","ILEUM","COLON")) +
  scale_fill_manual(values = "springgreen4") +
  scale_colour_manual(values = "springgreen4") +
  geom_boxplot(alpha=0.5, linewidth = 0.8) +
  geom_point() +
  theme(legend.position = "none",
        plot.title = element_text(size = 32),
        panel.border = element_rect(colour = "black", linewidth = 1),
        panel.grid.major.x = element_blank(),
        axis.title.y = element_text(size = 32),
        axis.text.y = element_text(colour = "black", size = 16),
        axis.ticks = element_line(colour = "black", linewidth = 1),
        axis.title.x = element_blank(),
        axis.text.x = element_text(colour = "black", size = 28, angle = 45, hjust=0.95, vjust=0.95))

ggplot(gut_clade2_fam_epith, aes(x= body_location, y= Abundance, fill = Family, colour = Family)) +
  theme_bw() +
  labs(title = "CLADE 2", y= "Bacteroidaceae (RA %)") +
  scale_y_continuous(expand = c(0.001,0,0.1,0), labels = scales::number_format(accuracy = 0.1)) +
  scale_x_discrete(limits = c("rumen","small intestine","large intestine"), labels = c("RUMEN","ILEUM","COLON")) +
  scale_fill_manual(values = "springgreen4") +
  scale_colour_manual(values = "springgreen4") +
  geom_boxplot(alpha=0.5, linewidth = 0.8) +
  geom_point() +
  theme(legend.position = "none",
        plot.title = element_text(size = 32),
        panel.border = element_rect(colour = "black", linewidth = 1),
        panel.grid.major.x = element_blank(),
        axis.title.y = element_text(size = 32),
        axis.text.y = element_text(colour = "black", size = 16),
        axis.ticks = element_line(colour = "black", linewidth = 1),
        axis.title.x = element_blank(),
        axis.text.x = element_text(colour = "black", size = 28, angle = 45, hjust=0.95, vjust=0.95))

pairwise.wilcox.test(gut_clade2_fam_lumen$Abundance, gut_clade2_fam_lumen$body_location, p.adjust.method = "BH") #LI higher than both
pairwise.wilcox.test(gut_clade2_fam_epith$Abundance, gut_clade2_fam_epith$body_location, p.adjust.method = "BH") # LI higher than both

## CLADE 3
## need to make a kingdom bar to include all discriminant families
gut_clade3_families_kingdom <- tax_glom(gut_clade3_disfamily, taxrank = "Kingdom", NArm = F)
gut_clade3_families_kingdom_melt <- psmelt(gut_clade3_families_kingdom)
gut_clade3_lumen_fam_kingdom_melt <- gut_clade3_families_kingdom_melt[which(gut_clade3_families_kingdom_melt$matrix=="lumen"),]
gut_clade3_epith_fam_kingdom_melt <- gut_clade3_families_kingdom_melt[which(gut_clade3_families_kingdom_melt$matrix=="epithelium"),]

ggplot(gut_clade3_lumen_fam_kingdom_melt, aes(x= body_location, y= Abundance, fill = Kingdom, colour = Kingdom)) +
  theme_bw() +
  labs(title = "CLADE 3", y= "Discriminant Families (RA %)") +
  scale_y_continuous(limits =c(0,72.5), expand = c(0.001,0,0.1,0), labels = scales::number_format(accuracy = 0.1)) +
  scale_x_discrete(limits = c("rumen","small intestine","large intestine"), labels = c("RUMEN","ILEUM","COLON")) +
  scale_fill_manual(values = "olivedrab4") +
  scale_colour_manual(values = "olivedrab4") +
  geom_boxplot(alpha=0.5, linewidth = 0.8) +
  geom_point() +
  theme(legend.position = "none",
        plot.title = element_text(size = 32),
        panel.border = element_rect(colour = "black", linewidth = 1),
        panel.grid.major.x = element_blank(),
        axis.title.y = element_text(size = 32),
        axis.text.y = element_text(colour = "black", size = 16),
        axis.ticks = element_line(colour = "black", linewidth = 1),
        axis.title.x = element_blank(),
        axis.text.x = element_text(colour = "black", size = 28, angle = 45, hjust=0.95, vjust=0.95))

ggplot(gut_clade3_epith_fam_kingdom_melt, aes(x= body_location, y= Abundance, fill = Kingdom, colour = Kingdom)) +
  theme_bw() +
  labs(title = "CLADE 3", y= "Discriminant Families (RA %)") +
  scale_y_continuous(limits =c(0,85), expand = c(0.001,0,0.1,0), labels = scales::number_format(accuracy = 0.1)) +
  scale_x_discrete(limits = c("rumen","small intestine","large intestine"), labels = c("RUMEN","ILEUM","COLON")) +
  scale_fill_manual(values = "olivedrab4") +
  scale_colour_manual(values = "olivedrab4") +
  geom_boxplot(alpha=0.5, linewidth = 0.8) +
  geom_point() +
  theme(legend.position = "none",
        plot.title = element_text(size = 32),
        panel.border = element_rect(colour = "black", linewidth = 1),
        panel.grid.major.x = element_blank(),
        axis.title.y = element_text(size = 32),
        axis.text.y = element_text(colour = "black", size = 16),
        axis.ticks = element_line(colour = "black", linewidth = 1),
        axis.title.x = element_blank(),
        axis.text.x = element_text(colour = "black", size = 28, angle = 45, hjust=0.95, vjust=0.95))

pairwise.wilcox.test(gut_clade3_lumen_fam_kingdom_melt$Abundance, gut_clade3_lumen_fam_kingdom_melt$body_location, p.adjust.method = "BH") # all different
pairwise.wilcox.test(gut_clade3_epith_fam_kingdom_melt$Abundance, gut_clade3_epith_fam_kingdom_melt$body_location, p.adjust.method = "BH") # SI higher than other two

### LOOK AT INDIVIDUAL FAMILIES
# splitting families by lumen/epi
gut_clade3_disfamily_melt_lumen <- gut_clade3_disfamily_melt[which(gut_clade3_disfamily_melt$matrix=="lumen"),]
gut_clade3_disfamily_melt_epith <- gut_clade3_disfamily_melt[which(gut_clade3_disfamily_melt$matrix=="epithelium"),]

clade3_fam_palette <- c("goldenrod2","#48CFCA","#E7C57D","#BEC77C","#F3B2E3","coral2","olivedrab4")

# plot by family
gut_clade3_lumen_Actinomycetaceae <- gut_clade3_disfamily_melt_lumen[which(gut_clade3_disfamily_melt_lumen$Family=="Actinomycetaceae"),]
gut_clade3_lumen_Anaerovoracaceae <- gut_clade3_disfamily_melt_lumen[which(gut_clade3_disfamily_melt_lumen$Family=="Anaerovoracaceae"),]
gut_clade3_lumen_Atopobiaceae <- gut_clade3_disfamily_melt_lumen[which(gut_clade3_disfamily_melt_lumen$Family=="Atopobiaceae"),]
gut_clade3_lumen_Erysipelotrichaceae <- gut_clade3_disfamily_melt_lumen[which(gut_clade3_disfamily_melt_lumen$Family=="Erysipelotrichaceae"),]
gut_clade3_lumen_Peptostreptococcaceae <- gut_clade3_disfamily_melt_lumen[which(gut_clade3_disfamily_melt_lumen$Family=="Peptostreptococcaceae"),]
gut_clade3_lumen_Peptostreptococcales_Tissierellales <- gut_clade3_disfamily_melt_lumen[which(gut_clade3_disfamily_melt_lumen$Family=="Peptostreptococcales-Tissierellales"),]
gut_clade3_lumen_Porphyromonadaceae <- gut_clade3_disfamily_melt_lumen[which(gut_clade3_disfamily_melt_lumen$Family=="Porphyromonadaceae"),]
gut_clade3_epith_Actinomycetaceae <- gut_clade3_disfamily_melt_epith[which(gut_clade3_disfamily_melt_epith$Family=="Actinomycetaceae"),]
gut_clade3_epith_Anaerovoraceae <- gut_clade3_disfamily_melt_epith[which(gut_clade3_disfamily_melt_epith$Family=="Anaerovoraceae"),]
gut_clade3_epith_Atopobiaceae <- gut_clade3_disfamily_melt_epith[which(gut_clade3_disfamily_melt_epith$Family=="Atopobiaceae"),]
gut_clade3_epith_Erysipelotrichaceae <- gut_clade3_disfamily_melt_epith[which(gut_clade3_disfamily_melt_epith$Family=="Erysipelotrichaceae"),]
gut_clade3_epith_Peptostreptococcaceae <- gut_clade3_disfamily_melt_epith[which(gut_clade3_disfamily_melt_epith$Family=="Peptostreptococcaceae"),]
gut_clade3_epith_Peptostreptococcaceae_Tissierellales <- gut_clade3_disfamily_melt_epith[which(gut_clade3_disfamily_melt_epith$Family=="Peptostreptococcaceae-Tissierellales"),]
gut_clade3_epith_Porphyromonadaceae <- gut_clade3_disfamily_melt_epith[which(gut_clade3_disfamily_melt_epith$Family=="Porphyromonadaceae"),]

ggplot(gut_clade3_lumen_Actinomycetaceae, aes(x=body_location, y= Abundance, fill = Family)) +
  theme_bw() +
  labs(y= "Actinomycetaceae (RA%)") +
  geom_bar(stat = "summary",colour = "black", size = 0.7) +
  geom_errorbar(stat = "summary", size = 0.7, width = 0.5) +
  scale_fill_manual(values = "goldenrod2") +
  #scale_colour_manual(values = "goldenrod2") +
  scale_y_continuous(expand = c(0.001,0,0.2,0)) +
  scale_x_discrete(limits = c("rumen","small intestine","large intestine"), labels = c("RUMEN","SM. INT","LG. INT")) +
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", size = 1),
        panel.grid.major.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
        axis.title.y = element_text(size = 30),
        axis.text.y = element_text(size = 16, colour = "black"),
        axis.ticks = element_line(colour = "black", size = 0.75))

pairwise.wilcox.test(gut_clade3_lumen_Actinomycetaceae$Abundance, gut_clade3_lumen_Actinomycetaceae$body_location, p.adjust.method = "BH")


pairwise.wilcox.test(gut_clade3_lumen_Anaerovoracaceae$Abundance,gut_clade3_lumen_Anaerovoracaceae$body_location, p.adjust.method = "BH")

## Atopobiaceae
ggplot(gut_clade3_lumen_Atopobiaceae, aes(x=body_location, y= Abundance, fill = Family)) +
  theme_bw() +
  labs(y= "Atopobiaceae (RA%)") +
  geom_bar(stat = "summary",colour = "black", size = 0.7) +
  geom_errorbar(stat = "summary", size = 0.7, width = 0.5) +
  scale_fill_manual(values = "#E7C57D") +
  #scale_colour_manual(values = "goldenrod2") +
  scale_y_continuous(expand = c(0.001,0,0.2,0), labels = scales::number_format(accuracy = 0.01)) +
  scale_x_discrete(limits = c("rumen","small intestine","large intestine"), labels = c("RUMEN","SM. INT","LG. INT")) +
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", size = 1),
        panel.grid.major.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
        axis.title.y = element_text(size = 30),
        axis.text.y = element_text(size = 16, colour = "black"),
        axis.ticks = element_line(colour = "black", size = 0.75))

## Erysipelotrichaceae
ggplot(gut_clade3_lumen_Erysipelotrichaceae, aes(x=body_location, y= Abundance, fill = Family)) +
  theme_bw() +
  labs(y= "Erysipelotrichaceae (RA%)") +
  geom_bar(stat = "summary",colour = "black", size = 0.7) +
  geom_errorbar(stat = "summary", size = 0.7, width = 0.5) +
  scale_fill_manual(values = "#BEC77C") +
  #scale_colour_manual(values = "goldenrod2") +
  scale_y_continuous(expand = c(0.001,0,0.2,0), labels = scales::number_format(accuracy = 0.01)) +
  scale_x_discrete(limits = c("rumen","small intestine","large intestine"), labels = c("RUMEN","SM. INT","LG. INT")) +
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", size = 1),
        panel.grid.major.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
        axis.title.y = element_text(size = 30),
        axis.text.y = element_text(size = 16, colour = "black"),
        axis.ticks = element_line(colour = "black", size = 0.75))

## Peptostreptococcaceae
ggplot(gut_clade3_lumen_Peptostreptococcaceae, aes(x=body_location, y= Abundance, fill = Family)) +
  theme_bw() +
  labs(y= "Peptostreptococcaceae (RA%)") +
  geom_bar(stat = "summary",colour = "black", size = 0.7) +
  geom_errorbar(stat = "summary", size = 0.7, width = 0.5) +
  scale_fill_manual(values = "#F3B2E3") +
  #scale_colour_manual(values = "goldenrod2") +
  scale_y_continuous(expand = c(0.001,0,0.2,0), labels = scales::number_format(accuracy = 0.01)) +
  scale_x_discrete(limits = c("rumen","small intestine","large intestine"), labels = c("RUMEN","SM. INT","LG. INT")) +
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", size = 1),
        panel.grid.major.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
        axis.title.y = element_text(size = 30),
        axis.text.y = element_text(size = 16, colour = "black"),
        axis.ticks = element_line(colour = "black", size = 0.75))

## Peptostreptococcales-Tissierellales
ggplot(gut_clade3_lumen_Peptostreptococcales_Tissierellales, aes(x=body_location, y= Abundance, fill = Family)) +
  theme_bw() +
  labs(y= "Peptostreptococcales-Tissierellales (RA%)") +
  geom_bar(stat = "summary",colour = "black", size = 0.7) +
  geom_errorbar(stat = "summary", size = 0.7, width = 0.5) +
  scale_fill_manual(values = "coral2") +
  #scale_colour_manual(values = "goldenrod2") +
  scale_y_continuous(expand = c(0.001,0,0.2,0), labels = scales::number_format(accuracy = 0.01)) +
  scale_x_discrete(limits = c("rumen","small intestine","large intestine"), labels = c("RUMEN","SM. INT","LG. INT")) +
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", size = 1),
        panel.grid.major.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
        axis.title.y = element_text(size = 30),
        axis.text.y = element_text(size = 16, colour = "black"),
        axis.ticks = element_line(colour = "black", size = 0.75))

## Porphyromonadaceae
ggplot(gut_clade3_lumen_Porphyromonadaceae, aes(x=body_location, y= Abundance, fill = Family)) +
  theme_bw() +
  labs(y= "Porphyromonadaceae (RA%)") +
  geom_bar(stat = "summary",colour = "black", size = 0.7) +
  geom_errorbar(stat = "summary", size = 0.7, width = 0.5) +
  scale_fill_manual(values = "#F3B2E3") +
  #scale_colour_manual(values = "goldenrod2") +
  scale_y_continuous(expand = c(0.001,0,0.2,0), labels = scales::number_format(accuracy = 0.001)) +
  scale_x_discrete(limits = c("rumen","small intestine","large intestine"), labels = c("RUMEN","SM. INT","LG. INT")) +
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", size = 1),
        panel.grid.major.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
        axis.title.y = element_text(size = 30),
        axis.text.y = element_text(size = 16, colour = "black"),
        axis.ticks = element_line(colour = "black", size = 0.75))

## all together now
ggplot(gut_clade3_disfamily_melt_epith, aes(x= body_location, y= Abundance)) +
  theme_bw() + coord_flip() + 
  labs(y= "Relative Abundance (%)") +
  geom_bar(stat = "summary", aes(fill=Family), colour = "black", size = 1) +
  #geom_errorbar(data = gut_clade3_lumen_fam_kingdom_melt, aes(x= body_location, y= Abundance), stat = "summary", size = 1, width = .5) +
  #scale_fill_manual(values = clade3_fam_palette) +
  scale_x_discrete(limits = c("large intestine","small intestine","rumen"), labels = c("LG. INT","SM. INT","RUMEN")) +
  scale_y_continuous(expand = c(0.001,0,.15,0)) +
  theme(#legend.position = "none",
        panel.border = element_rect(colour = "black", size = 1.25),
        panel.grid.major.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(colour = "black", size = 32),
        axis.title.x = element_text(size = 38),
        axis.text.x = element_text(size = 20, colour = "black"),
        axis.ticks = element_line(colour = "black", size = 1))

pairwise.wilcox.test(gut_clade3_lumen_fam_kingdom_melt$Abundance, gut_clade3_lumen_fam_kingdom_melt$body_location, p.adjust.method = "BH")

####################### DIFFERENTIALLY ABUNDANT TAXA BETWEEN LA + and - from the gut####
### INDENTIFY THE TAXA WITH LEFSE ####
## RUMEN LUMEN
rumen_lumen_genus.lefse <- merge_low_abundance(rumen_lumen_genus, threshold = 0.25)
rumen_lumen_genus.lefse #53 genera
write.csv(tax_table(rumen_lumen_genus.lefse),"rumen_lumen_lefse_taxa.csv")
write.csv(otu_table(rumen_lumen_genus.lefse),"rumen_lumen_lefse_otus.csv")
### LEFSE identified: Bifidobacterium and Syntrophococcus

## RUMEN EPITHELIUM
rumen_epithelium_genus.lefse <- merge_low_abundance(rumen_epithelium_genus, threshold = 0.25)
rumen_epithelium_genus.lefse #45 genera
write.csv(tax_table(rumen_epithelium_genus.lefse),"rumen_epithelium_lefse_taxa.csv")
write.csv(otu_table(rumen_epithelium_genus.lefse),"rumen_epithelium_lefse_otus.csv")
### LEFSE identified: unclassified members of Bacteroidales

## SM intestine LUMEN
sm_intestine_lumen_genus.lefse <- merge_low_abundance(sm_intestine_lumen_genus, threshold = 0.25)
sm_intestine_lumen_genus.lefse #30 genera
write.csv(tax_table(sm_intestine_lumen_genus.lefse),"sm_intestine_lumen_lefse_taxa.csv")
write.csv(otu_table(sm_intestine_lumen_genus.lefse),"sm_intestine_lumen_lefse_otus.csv")
### LEFSE identified: Bifidobacterium, Acetitomaculum, unclassified Lachnospiraceae 

## SM intestine EPITHELIUM
sm_intestine_epithelium_genus.lefse <- merge_low_abundance(sm_intestine_epithelium_genus, threshold = 0.25)
sm_intestine_epithelium_genus.lefse #29 genera
write.csv(tax_table(sm_intestine_epithelium_genus.lefse),"sm_intestine_epithelium_lefse_taxa.csv")
write.csv(otu_table(sm_intestine_epithelium_genus.lefse),"sm_intestine_epithelium_lefse_otus.csv")
### LEFSE identified: no taxa

## lg_intestine LUMEN
lg_intestine_lumen_genus.lefse <- merge_low_abundance(lg_intestine_lumen_genus, threshold = 0.25)
lg_intestine_lumen_genus.lefse #59 genera
write.csv(tax_table(lg_intestine_lumen_genus.lefse),"lg_intestine_lumen_lefse_taxa.csv")
write.csv(otu_table(lg_intestine_lumen_genus.lefse),"lg_intestine_lumen_lefse_otus.csv")
### LEFSE identified: NK4A214 group (Oscillospiraceae)

## lg_intestine EPITHELIUM
lg_intestine_epithelium_genus.lefse <- merge_low_abundance(lg_intestine_epithelium_genus, threshold = 0.25)
lg_intestine_epithelium_genus.lefse #60 genera
write.csv(tax_table(lg_intestine_epithelium_genus.lefse),"lg_intestine_epithelium_lefse_taxa.csv")
write.csv(otu_table(lg_intestine_epithelium_genus.lefse),"lg_intestine_epithelium_lefse_otus.csv")
### LEFSE identified: no taxa

### 6 genera across all GIT locations
liver_abscess_palette <- c("#001149","#c64f00")
### Pull them for plotting and stats ###
gut_ra <- transform_sample_counts(gut_samples.css, function(x) {x/sum(x)} * 100)
gut_ra_genus <- tax_glom(gut_ra, taxrank = "Genus", NArm = F)
gut_ra_genus_lumen <- subset_samples(gut_ra_genus, matrix=="lumen")
gut_ra_genus_epithelium <- subset_samples(gut_ra_genus, matrix=="epithelium")
gut_lumen_LAdiffAbund_genera <- subset_taxa(gut_ra_genus_lumen, Genus=="Bifidobacterium" | 
                                        Genus=="Syntrophococcus" | 
                                        Genus=="Acetitomaculum" |
                                        Genus=="unclassified Lachnospiraceae" |
                                        Genus=="NK4A214_group")
gut_lumen_LAdiffAbund_genera_melt <- psmelt(gut_lumen_LAdiffAbund_genera)

gut_epith_LAdiffAbund_genera <- subset_taxa(gut_ra_genus_epithelium, Genus=="unclassified Bacteroidales")
gut_epith_LAdiffAbund_genera_melt <- psmelt(gut_epith_LAdiffAbund_genera)

# LUMEN DIFFABUND TAXA
ggplot(gut_lumen_LAdiffAbund_genera_melt, aes(x= body_location, y= Abundance, fill= liver_abscesses, colour = liver_abscesses)) +
  theme_bw() +
  facet_wrap(~Genus, scales = "free_y", nrow = 1) +
  geom_bar(stat = "summary", position = position_dodge2(), alpha = 0.5) +
  geom_errorbar(stat = "summary", position = position_dodge2(padding = 0.4)) +
  scale_x_discrete(limits=c("rumen","small intestine","large intestine"), labels = c("RUMEN","ILEUM","COLON")) +
  scale_y_continuous(expand = c(0,0,0.2,0), labels = number_format(accuracy = 0.1)) +
  scale_colour_manual(values = liver_abscess_palette) +
  scale_fill_manual(values = liver_abscess_palette) +
  theme(legend.position = "none",
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 1),
        panel.grid.major.x = element_blank(),
        axis.title.y = element_text(size =  24),
        axis.text.y = element_text(colour = "black", size = 16),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 30, colour = "black", angle = 45, hjust = 0.95, vjust = 0.95))

# stats
# Acetitomaculum
rumen_lumen_LADiff_acetitomaculum <- gut_lumen_LAdiffAbund_genera_melt[which(gut_lumen_LAdiffAbund_genera_melt$Genus=="Acetitomaculum" & gut_lumen_LAdiffAbund_genera_melt$body_location=="rumen"),]
kruskal.test(rumen_lumen_LADiff_acetitomaculum$Abundance, rumen_lumen_LADiff_acetitomaculum$liver_abscesses) # NS

sm_intestine_lumen_LADiff_acetitomaculum <- gut_lumen_LAdiffAbund_genera_melt[which(gut_lumen_LAdiffAbund_genera_melt$Genus=="Acetitomaculum" & gut_lumen_LAdiffAbund_genera_melt$body_location=="small intestine"),]
kruskal.test(sm_intestine_lumen_LADiff_acetitomaculum$Abundance, sm_intestine_lumen_LADiff_acetitomaculum$liver_abscesses) # sig. p = 0.03

lg_intestine_lumen_LADiff_acetitomaculum <- gut_lumen_LAdiffAbund_genera_melt[which(gut_lumen_LAdiffAbund_genera_melt$Genus=="Acetitomaculum" & gut_lumen_LAdiffAbund_genera_melt$body_location=="large intestine"),]
kruskal.test(lg_intestine_lumen_LADiff_acetitomaculum$Abundance, lg_intestine_lumen_LADiff_acetitomaculum$liver_abscesses) # NS

# NK4A214_group
rumen_lumen_LADiff_NK4A214_group <- gut_lumen_LAdiffAbund_genera_melt[which(gut_lumen_LAdiffAbund_genera_melt$Genus=="NK4A214_group" & gut_lumen_LAdiffAbund_genera_melt$body_location=="rumen"),]
kruskal.test(rumen_lumen_LADiff_NK4A214_group$Abundance, rumen_lumen_LADiff_NK4A214_group$liver_abscesses) # NS

sm_intestine_lumen_LADiff_NK4A214_group <- gut_lumen_LAdiffAbund_genera_melt[which(gut_lumen_LAdiffAbund_genera_melt$Genus=="NK4A214_group" & gut_lumen_LAdiffAbund_genera_melt$body_location=="small intestine"),]
kruskal.test(sm_intestine_lumen_LADiff_NK4A214_group$Abundance, sm_intestine_lumen_LADiff_NK4A214_group$liver_abscesses) # sig. p = 0.03

lg_intestine_lumen_LADiff_NK4A214_group <- gut_lumen_LAdiffAbund_genera_melt[which(gut_lumen_LAdiffAbund_genera_melt$Genus=="NK4A214_group" & gut_lumen_LAdiffAbund_genera_melt$body_location=="large intestine"),]
kruskal.test(lg_intestine_lumen_LADiff_NK4A214_group$Abundance, lg_intestine_lumen_LADiff_NK4A214_group$liver_abscesses) # sig. p 0.04

# unclassified_Lachnospiraceae
rumen_lumen_LADiff_unclassified_Lachnospiraceae <- gut_lumen_LAdiffAbund_genera_melt[which(gut_lumen_LAdiffAbund_genera_melt$Genus=="unclassified Lachnospiraceae" & gut_lumen_LAdiffAbund_genera_melt$body_location=="rumen"),]
kruskal.test(rumen_lumen_LADiff_unclassified_Lachnospiraceae$Abundance, rumen_lumen_LADiff_unclassified_Lachnospiraceae$liver_abscesses) # NS

sm_intestine_lumen_LADiff_unclassified_Lachnospiraceae <- gut_lumen_LAdiffAbund_genera_melt[which(gut_lumen_LAdiffAbund_genera_melt$Genus=="unclassified Lachnospiraceae" & gut_lumen_LAdiffAbund_genera_melt$body_location=="small intestine"),]
kruskal.test(sm_intestine_lumen_LADiff_unclassified_Lachnospiraceae$Abundance, sm_intestine_lumen_LADiff_unclassified_Lachnospiraceae$liver_abscesses) # sig. p = 0.02

lg_intestine_lumen_LADiff_unclassified_Lachnospiraceae <- gut_lumen_LAdiffAbund_genera_melt[which(gut_lumen_LAdiffAbund_genera_melt$Genus=="unclassified Lachnospiraceae" & gut_lumen_LAdiffAbund_genera_melt$body_location=="large intestine"),]
kruskal.test(lg_intestine_lumen_LADiff_unclassified_Lachnospiraceae$Abundance, lg_intestine_lumen_LADiff_unclassified_Lachnospiraceae$liver_abscesses) # NS

# Syntrophococcus
rumen_lumen_LADiff_Syntrophococcus <- gut_lumen_LAdiffAbund_genera_melt[which(gut_lumen_LAdiffAbund_genera_melt$Genus=="Syntrophococcus" & gut_lumen_LAdiffAbund_genera_melt$body_location=="rumen"),]
kruskal.test(rumen_lumen_LADiff_Syntrophococcus$Abundance, rumen_lumen_LADiff_Syntrophococcus$liver_abscesses) # NS (p=0.06)

sm_intestine_lumen_LADiff_Syntrophococcus <- gut_lumen_LAdiffAbund_genera_melt[which(gut_lumen_LAdiffAbund_genera_melt$Genus=="Syntrophococcus" & gut_lumen_LAdiffAbund_genera_melt$body_location=="small intestine"),]
kruskal.test(sm_intestine_lumen_LADiff_Syntrophococcus$Abundance, sm_intestine_lumen_LADiff_Syntrophococcus$liver_abscesses) # NS

lg_intestine_lumen_LADiff_Syntrophococcus <- gut_lumen_LAdiffAbund_genera_melt[which(gut_lumen_LAdiffAbund_genera_melt$Genus=="Syntrophococcus" & gut_lumen_LAdiffAbund_genera_melt$body_location=="large intestine"),]
kruskal.test(lg_intestine_lumen_LADiff_Syntrophococcus$Abundance, lg_intestine_lumen_LADiff_Syntrophococcus$liver_abscesses) # NS

# Bifidobacterium
rumen_lumen_LADiff_Bifidobacterium <- gut_lumen_LAdiffAbund_genera_melt[which(gut_lumen_LAdiffAbund_genera_melt$Genus=="Bifidobacterium" & gut_lumen_LAdiffAbund_genera_melt$body_location=="rumen"),]
kruskal.test(rumen_lumen_LADiff_Bifidobacterium$Abundance, rumen_lumen_LADiff_Bifidobacterium$liver_abscesses) # sig. p=0.04

sm_intestine_lumen_LADiff_Bifidobacterium <- gut_lumen_LAdiffAbund_genera_melt[which(gut_lumen_LAdiffAbund_genera_melt$Genus=="Bifidobacterium" & gut_lumen_LAdiffAbund_genera_melt$body_location=="small intestine"),]
kruskal.test(sm_intestine_lumen_LADiff_Bifidobacterium$Abundance, sm_intestine_lumen_LADiff_Bifidobacterium$liver_abscesses) # sig. p = 0.02

lg_intestine_lumen_LADiff_Bifidobacterium <- gut_lumen_LAdiffAbund_genera_melt[which(gut_lumen_LAdiffAbund_genera_melt$Genus=="Bifidobacterium" & gut_lumen_LAdiffAbund_genera_melt$body_location=="large intestine"),]
kruskal.test(lg_intestine_lumen_LADiff_Bifidobacterium$Abundance, lg_intestine_lumen_LADiff_Bifidobacterium$liver_abscesses) # NS


# EPITH DIFFABUND TAXA
ggplot(gut_epith_LAdiffAbund_genera_melt, aes(x= body_location, y= Abundance, fill= liver_abscesses, colour = liver_abscesses)) +
  theme_bw() +
  facet_wrap(~Genus, scales = "free_y", nrow = 1) +
  geom_bar(stat = "summary", position = position_dodge2(), alpha = 0.5) +
  geom_errorbar(stat = "summary", position = position_dodge2(padding = 0.4)) +
  scale_x_discrete(limits=c("rumen","small intestine","large intestine"), labels = c("RUMEN","ILEUM","COLON")) +
  scale_y_continuous(expand = c(0,0,0.2,0), labels = number_format(accuracy = 0.1)) +
  scale_colour_manual(values = liver_abscess_palette) +
  scale_fill_manual(values = liver_abscess_palette) +
  theme(legend.position = "none",
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 1),
        panel.grid.major.x = element_blank(),
        axis.title.y = element_text(size =  24),
        axis.text.y = element_text(colour = "black", size = 16),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 30, colour = "black", angle = 45, hjust = 0.95, vjust = 0.95))

# stats
rumen_epith_unBacteroidales <- gut_epith_LAdiffAbund_genera_melt[which(gut_epith_LAdiffAbund_genera_melt$body_location=="rumen"),]
kruskal.test(rumen_epith_unBacteroidales$Abundance, rumen_epith_unBacteroidales$liver_abscesses) # sig. p=0.03

SI_epith_unBacteroidales <- gut_epith_LAdiffAbund_genera_melt[which(gut_epith_LAdiffAbund_genera_melt$body_location=="small intestine"),]
kruskal.test(SI_epith_unBacteroidales$Abundance, SI_epith_unBacteroidales$liver_abscesses) # NS (only present in 1 sample)

LI_epith_unBacteroidales <- gut_epith_LAdiffAbund_genera_melt[which(gut_epith_LAdiffAbund_genera_melt$body_location=="large intestine"),]
kruskal.test(LI_epith_unBacteroidales$Abundance, LI_epith_unBacteroidales$liver_abscesses) # NS (only present in 13 samples)


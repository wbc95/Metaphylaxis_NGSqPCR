#### set WD ####
setwd("/Users/ljpinnell/Documents/VERO/WPT/R_analysis/")

#### load libraries ####
library(phyloseq); library(ggplot2); library(btools); library(dplyr);
library(tidyr); library(stringr); library(randomcoloR); library(metagMisc);
library(metagenomeSeq); library(GUniFrac); library(randomForest); library(pairwiseAdonis)
library(knitr); library(readr); library(kableExtra); library(scales); library(vegan);
library(ggdendro); library(forcats); library(ggbreak); library(BiodiversityR); library(ggrepel); library(ggforce)

#### source some functions ####
source("/Users/ljpinnell/Documents/R_Functions_Scripts/changeSILVATaxaNames.R")
source("/Users/ljpinnell/Documents/R_Functions_Scripts/changeUNITETaxaNames.R")
source("/Users/ljpinnell/Documents/R_Functions_Scripts/w_unifrac.R")
source("/Users/ljpinnell/Documents/R_Functions_Scripts/g_unifrac.R")
source("/Users/ljpinnell/Documents/R_Functions_Scripts/uw_unifrac.R")
source("/Users/ljpinnell/Documents/R_Functions_Scripts/MergeLowAbund.R")

#################### IMPORTING DATA ####
microbiome_data <- import_biom("../16S/table-with-taxonomy.biom","../16S/tree.nwk","../16S/dna-sequences.fasta")
microbiome_data # 664 samples, 17,401 taxa
#write.csv(sample_names(microbiome_data),"16S_sample_names.csv") # sample names to build metadata file

ITS_data <- import_biom("../ITS/table-with-taxonomy.biom","../ITS/tree.nwk","../ITS/dna-sequences.fasta")
ITS_data # 672 samples, 24,390 taxa
#write.csv(sample_names(ITS_data),"ITS_sample_names.csv") # sample names to build metadata file

### import mapping file
map_file <- import_qiime_sample_data("WPT_metadata.txt")
map_file

### combining mapping and data
microbiome <- merge_phyloseq(microbiome_data, map_file)
microbiome # 663 samples, 17401 taxa

ITS <- merge_phyloseq(ITS_data, map_file)
ITS # 672 samples, 24390 taxa

#################### PRE-PROCESSING ####
########### 16S ####
## renaming ranks
rank_names(microbiome) # "Rank1" - "Rank7" not ideal, lets change em
colnames(tax_table(microbiome)) <- c("Domain","Kingdom","Phylum","Class","Order","Family","Genus")
rank_names(microbiome) # beauty, now they are named properly

# changing the SILVA style naming (k__Bacteria, etc.)
tax.data.microbiome <- data.frame(tax_table(microbiome)) # extract the taxonomy table as a data frame
tax.data.names.microbiome <- changeSILVAtaxa(tax.data.microbiome) # this gets rid of the GG format

# now to change the NAs to a better naming scheme
for (i in 1:7){ tax.data.names.microbiome[,i] <- as.character(tax.data.names.microbiome[,i])} # converting all columns to characters
tax.data.names.microbiome[is.na(tax.data.names.microbiome)] <- "" # replacing the NAs with an empty string

### adding information to the 'uncultured' taxa
for (i in 1:nrow(tax.data.names.microbiome)){
  if (tax.data.names.microbiome[i,2] == "uncultured"){
    domain <- paste("uncultured", tax.data.names.microbiome[i,1])
    tax.data.names.microbiome[i, 2:7] <- domain
  } else if (tax.data.names.microbiome[i,3] == "uncultured"){
    kingdom <- paste("uncultured", tax.data.names.microbiome[i,2])
    tax.data.names.microbiome[i, 3:7] <- kingdom
  } else if (tax.data.names.microbiome[i,4] == "uncultured"){
    phylum <- paste("uncultured", tax.data.names.microbiome[i,3])
    tax.data.names.microbiome[i, 4:7] <- phylum
  } else if (tax.data.names.microbiome[i,5] == "uncultured"){
    class <- paste("uncultured", tax.data.names.microbiome[i,4])
    tax.data.names.microbiome[i, 5:7] <- class
  } else if (tax.data.names.microbiome[i,6] == "uncultured"){
    order <- paste("uncultured", tax.data.names.microbiome[i,5])
    tax.data.names.microbiome[i, 6:7] <- order
  } else if (tax.data.names.microbiome[i,7] == "uncultured"){
    tax.data.names.microbiome$Genus[i] <- paste("uncultured",tax.data.names.microbiome$Family[i])
  }
}

# now filling in the empty slots with the highest assigned taxonomy
for (i in 1:nrow(tax.data.names.microbiome)){
  if (tax.data.names.microbiome[i,2] == ""){
    domain <- paste("unclassified", tax.data.names.microbiome[i,1])
    tax.data.names.microbiome[i, 2:7] <- domain
  } else if (tax.data.names.microbiome[i,3] == ""){
    kingdom <- paste("unclassified", tax.data.names.microbiome[i,2])
    tax.data.names.microbiome[i, 3:7] <- kingdom
  } else if (tax.data.names.microbiome[i,4] == ""){
    phylum <- paste("unclassified", tax.data.names.microbiome[i,3])
    tax.data.names.microbiome[i, 4:7] <- phylum
  } else if (tax.data.names.microbiome[i,5] == ""){
    class <- paste("unclassified", tax.data.names.microbiome[i,4])
    tax.data.names.microbiome[i, 5:7] <- class
  } else if (tax.data.names.microbiome[i,6] == ""){
    order <- paste("unclassified", tax.data.names.microbiome[i,5])
    tax.data.names.microbiome[i, 6:7] <- order
  } else if (tax.data.names.microbiome[i,7] == ""){
    tax.data.names.microbiome$Genus[i] <- paste("unclassified",tax.data.names.microbiome$Family[i])
  }
}

head(tax.data.names.microbiome, 50) # great, no more NAs and no more k__
tax_table(microbiome) <- as.matrix(tax.data.names.microbiome) # re-insert the taxonomy table into the phyloseq object
tail(tax_table(microbiome), 50) # sweet, lookin good!

#### remove eukaryota
microbiome <- subset_taxa(microbiome, Domain!="Eukaryota")
microbiome # 17,349 taxa; lost 52 taxa

#### split samples from controls
microbiome_controls <- subset_samples(microbiome, sample_type=="control")
microbiome_controls <- prune_taxa(taxa_sums(controls) > 0, controls)
microbiome_controls #21 samples, 134 ASVs across all controls
sort(sample_sums(microbiome_controls)) # all good, but BLKP7 has 3,079 ASVs. check later

microbiome <- subset_samples(microbiome, sample_type=="sample")
microbiome <- prune_taxa(taxa_sums(microbiome) > 0, microbiome)
microbiome # 624 samples, 17,327 taxa

#### some QC checks
min(sample_sums(microbiome)) #6
max(sample_sums(microbiome)) # 16,990
mean(sample_sums(microbiome)) # 9065
median(sample_sums(microbiome)) # 9701
sort(sample_sums(microbiome)) # ... going to cut off at 1000 ASVs

microbiome <- prune_samples(sample_sums(microbiome) > 1000, microbiome)
microbiome <- prune_taxa(taxa_sums(microbiome) > 0, microbiome)
microbiome # 610 samples, 17,094 taxa remain

########### ITS ####
rank_names(ITS)
colnames(tax_table(ITS)) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
rank_names(ITS)

# changing the SILVA style naming (k__Bacteria, etc.)
tax.data.ITS <- data.frame(tax_table(ITS)) # extract the taxonomy table as a data frame
tax.data.names.ITS <- changeUNITEtaxa(tax.data.ITS) # this gets rid of the GG format

# now to change the NAs to a better naming scheme
for (i in 1:7){ tax.data.names.ITS[,i] <- as.character(tax.data.names.ITS[,i])} # converting all columns to characters
tax.data.names.ITS[is.na(tax.data.names.ITS)] <- "" # replacing the NAs with an empty string

### adding information to the 'uncultured' taxa
for (i in 1:nrow(tax.data.names.ITS)){
  if (tax.data.names.ITS[i,2] == "uncultured"){
    kingdom <- paste("uncultured", tax.data.names.ITS[i,1])
    tax.data.names.ITS[i, 2:7] <- kingdom
  } else if (tax.data.names.ITS[i,3] == "uncultured"){
    phylum <- paste("uncultured", tax.data.names.ITS[i,2])
    tax.data.names.ITS[i, 3:7] <- phylum
  } else if (tax.data.names.ITS[i,4] == "uncultured"){
    class<- paste("uncultured", tax.data.names.ITS[i,3])
    tax.data.names.ITS[i, 4:7] <- class
  } else if (tax.data.names.ITS[i,5] == "uncultured"){
    order <- paste("uncultured", tax.data.names.ITS[i,4])
    tax.data.names.ITS[i, 5:7] <- order
  } else if (tax.data.names.ITS[i,6] == "uncultured"){
    family <- paste("uncultured", tax.data.names.ITS[i,5])
    tax.data.names.ITS[i, 6:7] <- genus
  } else if (tax.data.names.ITS[i,7] == "uncultured"){
    tax.data.names.ITS$Species[i] <- paste("uncultured",tax.data.names.ITS$Genus[i])
  }
}

# now filling in the empty slots with the highest assigned taxonomy
for (i in 1:nrow(tax.data.names.ITS)){
  if (tax.data.names.ITS[i,2] == ""){
    kingdom <- paste("unclassified", tax.data.names.ITS[i,1])
    tax.data.names.ITS[i, 2:7] <- kingdom
  } else if (tax.data.names.ITS[i,3] == ""){
    phylum <- paste("unclassified", tax.data.names.ITS[i,2])
    tax.data.names.ITS[i, 3:7] <- phylum
  } else if (tax.data.names.ITS[i,4] == ""){
    class <- paste("unclassified", tax.data.names.ITS[i,3])
    tax.data.names.ITS[i, 4:7] <- class
  } else if (tax.data.names.ITS[i,5] == ""){
    order <- paste("unclassified", tax.data.names.ITS[i,4])
    tax.data.names.ITS[i, 5:7] <- order
  } else if (tax.data.names.ITS[i,6] == ""){
    family <- paste("unclassified", tax.data.names.ITS[i,5])
    tax.data.names.ITS[i, 6:7] <- family
  } else if (tax.data.names.ITS[i,7] == ""){
    tax.data.names.ITS$Species[i] <- paste("unclassified",tax.data.names.ITS$Genus[i])
  }
}

tail(tax.data.names.ITS, 50) # great, no more NAs and no more k__
tax_table(ITS) <- as.matrix(tax.data.names.ITS) # re-insert the taxonomy table into the phyloseq object
tail(tax_table(ITS), 20) # sweet, lookin good!

#### split samples from controls
ITS_controls <- subset_samples(ITS, sample_type=="control")
ITS_controls <- prune_taxa(taxa_sums(ITS_controls) > 0, ITS_controls)
ITS_controls #30 samples, 278 ASVs across all controls
sort(sample_sums(ITS_controls)) # all good, but this includes the Mata-Mata samples I should look at later

ITS <- subset_samples(ITS, sample_type=="sample")
ITS <- prune_taxa(taxa_sums(ITS) > 0, ITS)
ITS # 642 samples, 24,199 taxa

#### some QC checks
min(sample_sums(ITS)) #1
max(sample_sums(ITS)) # 28,132
mean(sample_sums(ITS)) # 7071
median(sample_sums(ITS)) # 6615
sort(sample_sums(ITS)) # ... going to cut off at 1000 ASVs

ITS <- prune_samples(sample_sums(ITS) > 1000, ITS)
ITS <- prune_taxa(taxa_sums(ITS) > 0, ITS)
ITS # 614 samples, 23,993 taxa remain


########################################### ALPHA DIVERSITY #####
alpha_div1_16S <- estimate_richness(microbiome, measures = c("Observed", "Shannon", "Simpson","InvSimpson"))
alpha_div2_16S <- estimate_pd(microbiome)

alpha_div_16S <- cbind(alpha_div1_16S, alpha_div2_16S)
alpha_div_16S
alpha_div_16S <- alpha_div_16S[,c(1:5)]
alpha_div_16S
alpha_div_16S_meta <- cbind(as.data.frame(microbiome@sam_data), alpha_div_16S)
alpha_div_16S_meta # metadata and div metrics

alpha_div1_ITS <- estimate_richness(ITS, measures = c("Observed", "Shannon", "Simpson","InvSimpson"))
alpha_div2_ITS <- estimate_pd(ITS)

alpha_div_ITS <- cbind(alpha_div1_ITS, alpha_div2_ITS)
alpha_div_ITS
alpha_div_ITS <- alpha_div_ITS[,c(1:5)]
alpha_div_ITS
alpha_div_ITS_meta <- cbind(as.data.frame(ITS@sam_data), alpha_div_ITS)
alpha_div_ITS_meta # metadata and div metrics

######## RICHNESS ####
#### 16S
ggplot(alpha_div_16S_meta, aes(x= washed_unwashed, y= Observed)) +
  facet_grid(~matrix_specific, scales = "free_x") +
  geom_boxplot(alpha = 0.5) +
  geom_point()

########################################### CSS NORMALIZATION ####
microbiome.css <- phyloseq_transform_css(microbiome, log = F)
ITS.css <- phyloseq_transform_css(ITS, log = F)
########################################### REL ABUND AND TAX GLOMS #####
microbiome_ra <- transform_sample_counts(microbiome.css, function(x) {x/sum(x)} * 100)
microbiome_ra_phylum <- tax_glom(microbiome_ra, taxrank = "Phylum", NArm = F)
microbiome_ra_order <- tax_glom(microbiome_ra, taxrank = "Order", NArm = F)
microbiome_ra_family <- tax_glom(microbiome_ra, taxrank = "Family", NArm = F)
microbiome_ra_genus <- tax_glom(microbiome_ra, taxrank = "Genus", NArm = F)

## making palettes
microbiome_ra_family_filt <- merge_low_abundance(microbiome_ra_family, threshold = 0.05) #137
microbiome_ra_family_filt_palette <- distinctColorPalette(137)
microbiome_ra_family_filt_melt <- psmelt(microbiome_ra_family_filt)

ggplot(microbiome_ra_family_filt_melt, aes(x= sample_type, y= Abundance, fill = Family)) +
  geom_bar(stat = "summary") +
  scale_fill_manual(values = microbiome_ra_family_filt_palette)

write.table(microbiome_ra_family_filt_palette,"16S_family_colours.txt", sep = "\t", col.names = 1)
write.table(otu_table(microbiome_ra_family_filt),"16S_family_otus.txt", sep = "\t")
write.table(tax_table(microbiome_ra_family_filt),"16S_family_taxa.txt", sep = "\t", col.names = T)


ITS_ra <- transform_sample_counts(ITS.css, function(x) {x/sum(x)} * 100)
ITS_ra_phylum <- tax_glom(ITS_ra, taxrank = "Phylum", NArm = F)
ITS_ra_order <- tax_glom(ITS_ra, taxrank = "Order", NArm = F)
ITS_ra_family <- tax_glom(ITS_ra, taxrank = "Family", NArm = F) #441 families

## make ITS palette
ITS_ra_family_filt <- merge_low_abundance(ITS_ra_family, threshold = 0.05) #111 families
ITS_ra_family_filt_palette <- distinctColorPalette(111)
write.csv(ITS_ra_family_filt_palette, "ITS_family_colours.csv")
write.csv(tax_table(ITS_ra_family_filt), "ITS_family_taxa.csv")


ITS_ra_genus <- tax_glom(ITS_ra, taxrank = "Genus", NArm = F)
########################################### BETA DIVERSITY #####
#### EGGS + TURTLES ####
egg_turtle_16S.css <- subset_samples(microbiome.css, matrix=="egg" | matrix=="turtle")
egg_turtle_16S.css <- prune_taxa(taxa_sums(egg_turtle_16S.css) > 0, egg_turtle_16S.css)


egg_turtle_ITS.css <- subset_samples(ITS.css, matrix=="egg" | matrix=="turtle")
egg_turtle_ITS.css <- prune_taxa(taxa_sums(egg_turtle_ITS.css) > 0, egg_turtle_ITS.css)

###### ordinations
egg_turtle_16S.gunifrac <- gunifrac(egg_turtle_16S.css)
egg_turtle_ITS.gunifrac <- gunifrac(egg_turtle_ITS.css)

egg_turtle_16S.ord <- ordinate(egg_turtle_16S.css, method = "NMDS", distance = egg_turtle_16S.gunifrac)
egg_turtle_ITS.ord <- ordinate(egg_turtle_ITS.css, method = "NMDS", distance = egg_turtle_ITS.gunifrac)

plot_ordination(egg_turtle_16S.css, egg_turtle_16S.ord, color = "matrix_specific") +
  facet_grid(~nest_cat)

plot_ordination(egg_turtle_ITS.css, egg_turtle_ITS.ord, color = "matrix_specific") +
  facet_grid(~nest_cat)

############################################ EACH NEST INDIVIDUALLY ####
# NEST A ####
nestA_matrix_timepoint_palette <- c("#e14f00","#ff7800","#ffbb7f","#522700","#256602","#7bcf37","#44880d","#a9f46d","#183500","#008078","#44d6cd","#009caa","#8bf8f2","#003f3b","#8e4da3","#4e205d")
nestA_matrix_timepoint_ITS_palette <- c("#e14f00","#ff7800","#ffbb7f","#522700","#256602","#7bcf37","#44880d","#a9f46d","#183500","#008078","#44d6cd","#009caa","#8bf8f2","#003f3b","#3b0042","#8e4da3","#4e205d")

microbiome_nestA.css <- subset_samples(microbiome.css, nest_cat=="A")
microbiome_nestA.css <- prune_taxa(taxa_sums(microbiome_nestA.css) > 0, microbiome_nestA.css)
ITS_nestA.css <- subset_samples(ITS.css, nest_cat=="A")
ITS_nestA.css <- prune_taxa(taxa_sums(ITS.css) > 0, ITS.css)

microbiome_nestA.gunifrac <- gunifrac(microbiome_nestA.css)
ITS_nestA.gunifrac <- gunifrac(ITS_nestA.css)

microbiome_nestA.ord <- ordinate(microbiome_nestA.css, method = "PCoA", distance = microbiome_nestA.gunifrac)
ITS_nestA.ord <- ordinate(ITS_nestA.css, method = "PCoA", distance = ITS_nestA.gunifrac)

plot_ordination(microbiome_nestA.css, microbiome_nestA.ord, color = "matrix_timepoint_cat") +
  theme_bw() +
  labs(x= "Axis 1 (31.8% variation explained)", y= "Axis 2 (23.2% variation explained)") +
  scale_colour_manual(values = nestA_matrix_timepoint_palette) +
  scale_fill_manual(values = nestA_matrix_timepoint_palette) +
  geom_point(size = 7, shape = 18) +
  stat_ellipse(geom = "polygon", aes(fill= matrix_timepoint_cat), alpha = 0.25, lty = 2, linewidth = 1) +
  theme(legend.position = "none",
        panel.border = element_rect(linewidth = 1, colour = "black"),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 12, colour = "black"))

plot_ordination(ITS_nestA.css, ITS_nestA.ord, color = "matrix_timepoint_cat") +
  theme_bw() +
  labs(x= "Axis 1 (23.9% variation explained)", y= "Axis 2 (13.6% variation explained)") +
  scale_colour_manual(values = nestA_matrix_timepoint_ITS_palette) +
  scale_fill_manual(values = nestA_matrix_timepoint_ITS_palette) +
  geom_point(size = 7, shape = 18) +
  stat_ellipse(geom = "polygon", aes(fill= matrix_timepoint_cat), alpha = 0.25, lty = 2, linewidth = 1) +
  theme(legend.position = "none",
        panel.border = element_rect(linewidth = 1, colour = "black"),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 12, colour = "black"))

# NEST B ####
nestB_matrix_timepoint_palette <- c("#e14f00","#ff7800","#ffbb7f","#522700","#256602","#7bcf37","#44880d","#a9f46d","#183500","#008078","#44d6cd","#009caa","#8bf8f2","#003f3b","#8e4da3","#4e205d")
nestB_matrix_timepoint_ITS_palette <- c("#e14f00","#ff7800","#ffbb7f","#522700","#256602","#7bcf37","#44880d","#a9f46d","#183500","#008078","#44d6cd","#009caa","#8bf8f2","#003f3b","#3b0042","#8e4da3","#4e205d")

microbiome_nestB.css <- subset_samples(microbiome.css, nest_cat=="B")
microbiome_nestB.css <- prune_taxa(taxa_sums(microbiome_nestB.css) > 0, microbiome_nestB.css)
ITS_nestB.css <- subset_samples(ITS.css, nest_cat=="B")
ITS_nestB.css <- prune_taxa(taxa_sums(ITS_nestB.css) > 0, ITS_nestB.css)

microbiome_nestB.gunifrac <- gunifrac(microbiome_nestB.css)
ITS_nestB.gunifrac <- gunifrac(ITS_nestB.css)

microbiome_nestB.ord <- ordinate(microbiome_nestB.css, method = "PCoA", distance = microbiome_nestB.gunifrac)
ITS_nestB.ord <- ordinate(ITS_nestB.css, method = "PCoA", distance = ITS_nestB.gunifrac)

plot_ordination(microbiome_nestB.css, microbiome_nestB.ord, color = "matrix_timepoint_cat") +
  theme_bw() +
  labs(x= "Axis 1 (30.8% variation explained)", y= "Axis 2 (18.1% variation explained)") +
  scale_colour_manual(values = nestB_matrix_timepoint_ITS_palette) +
  scale_fill_manual(values = nestB_matrix_timepoint_ITS_palette) +
  geom_point(size = 7, shape = 18) +
  stat_ellipse(geom = "polygon", aes(fill= matrix_timepoint_cat), alpha = 0.25, lty = 2, linewidth = 1) +
  theme(legend.position = "none",
        panel.border = element_rect(linewidth = 1, colour = "black"),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 12, colour = "black"))

plot_ordination(ITS_nestB.css, ITS_nestB.ord, color = "matrix_timepoint_cat") +
  theme_bw() +
  labs(x= "Axis 1 (23.3% variation explained)", y= "Axis 2 (13.9% variation explained)") +
  scale_colour_manual(values = nestB_matrix_timepoint_ITS_palette) +
  scale_fill_manual(values = nestB_matrix_timepoint_ITS_palette) +
  geom_point(size = 7, shape = 18) +
  stat_ellipse(geom = "polygon", aes(fill= matrix_timepoint_cat), alpha = 0.25, lty = 2, linewidth = 1) +
  theme(legend.position = "none",
        panel.border = element_rect(linewidth = 1, colour = "black"),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 12, colour = "black"))

# NEST C ####
nestC_matrix_timepoint_palette <- c("#e14f00","#ff7800","#ffbb7f","#522700","#256602","#7bcf37","#44880d","#a9f46d","#183500","#008078","#44d6cd","#009caa","#8bf8f2","#003f3b","grey35","grey35","#3b0042","#8e4da3","#4e205d")
nestC_matrix_timepoint_ITS_palette <- c("#e14f00","#ff7800","#ffbb7f","#522700","#256602","#7bcf37","#44880d","#a9f46d","#183500","#008078","#44d6cd","#009caa","#8bf8f2","#003f3b","#3b0042","#8e4da3","#4e205d")

microbiome_nestC.css <- subset_samples(microbiome.css, nest_cat=="C")
microbiome_nestC.css <- prune_taxa(taxa_sums(microbiome_nestC.css) > 0, microbiome_nestC.css)
ITS_nestC.css <- subset_samples(ITS.css, nest_cat=="C")
ITS_nestC.css <- prune_taxa(taxa_sums(ITS_nestC.css) > 0, ITS_nestC.css)

microbiome_nestC.gunifrac <- gunifrac(microbiome_nestC.css)
ITS_nestC.gunifrac <- gunifrac(ITS_nestC.css)

microbiome_nestC.ord <- ordinate(microbiome_nestC.css, method = "PCoA", distance = microbiome_nestC.gunifrac)
ITS_nestC.ord <- ordinate(ITS_nestC.css, method = "PCoA", distance = ITS_nestC.gunifrac)

plot_ordination(microbiome_nestC.css, microbiome_nestC.ord, color = "matrix_timepoint_cat") +
  theme_bw() +
  labs(x= "Axis 1 (29.2% variation explained)", y= "Axis 2 (21.8% variation explained)") +
  scale_colour_manual(values = nestC_matrix_timepoint_palette) +
  scale_fill_manual(values = nestC_matrix_timepoint_palette) +
  geom_point(size = 7, shape = 18) +
  stat_ellipse(geom = "polygon", aes(fill= matrix_timepoint_cat), alpha = 0.25, lty = 2, linewidth = 1) +
  theme(legend.position = "none",
        panel.border = element_rect(linewidth = 1, colour = "black"),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 12, colour = "black"))

plot_ordination(ITS_nestC.css, ITS_nestC.ord, color = "matrix_timepoint_cat") +
  theme_bw() +
  labs(x= "Axis 1 (22.2% variation explained)", y= "Axis 2 (12.3% variation explained)") +
  scale_colour_manual(values = nestC_matrix_timepoint_palette) +
  scale_fill_manual(values = nestC_matrix_timepoint_palette) +
  geom_point(size = 7, shape = 18) +
  stat_ellipse(geom = "polygon", aes(fill= matrix_timepoint_cat), alpha = 0.25, lty = 2, linewidth = 1) +
  theme(legend.position = "none",
        panel.border = element_rect(linewidth = 1, colour = "black"),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 12, colour = "black"))

# NEST D ####
nestD_matrix_timepoint_palette <- c("#e14f00","#ff7800","#ffbb7f","#522700","#256602","#7bcf37","#44880d","#a9f46d","#183500","#008078","#44d6cd","#009caa","#8bf8f2","#003f3b","#8e4da3","#4e205d")
nestD_matrix_timepoint_ITS_palette <- c("#e14f00","#ff7800","#ffbb7f","#522700","#256602","#7bcf37","#44880d","#a9f46d","#183500","#008078","#44d6cd","#009caa","#8bf8f2","#003f3b","#3b0042","#8e4da3","#4e205d")

microbiome_nestD.css <- subset_samples(microbiome.css, nest_cat=="D")
microbiome_nestD.css <- prune_taxa(taxa_sums(microbiome_nestD.css) > 0, microbiome_nestD.css)
ITS_nestD.css <- subset_samples(ITS.css, nest_cat=="D")
ITS_nestD.css <- prune_taxa(taxa_sums(ITS_nestD.css) > 0, ITS_nestD.css)

microbiome_nestD.gunifrac <- gunifrac(microbiome_nestD.css)
ITS_nestD.gunifrac <- gunifrac(ITS_nestD.css)

microbiome_nestD.ord <- ordinate(microbiome_nestD.css, method = "PCoA", distance = microbiome_nestD.gunifrac)
ITS_nestD.ord <- ordinate(ITS_nestD.css, method = "PCoA", distance = ITS_nestD.gunifrac)

plot_ordination(microbiome_nestD.css, microbiome_nestD.ord, color = "matrix_timepoint_cat") +
  theme_bw() +
  labs(x= "Axis 1 (29.0% variation explained)", y= "Axis 2 (17.0% variation explained)") +
  scale_colour_manual(values = nestD_matrix_timepoint_ITS_palette) +
  scale_fill_manual(values = nestD_matrix_timepoint_ITS_palette) +
  geom_point(size = 7, shape = 18) +
  stat_ellipse(geom = "polygon", aes(fill= matrix_timepoint_cat), alpha = 0.25, lty = 2, linewidth = 1) +
  theme(legend.position = "none",
        panel.border = element_rect(linewidth = 1, colour = "black"),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 12, colour = "black"))

plot_ordination(ITS_nestD.css, ITS_nestD.ord, color = "matrix_timepoint_cat") +
  theme_bw() +
  labs(x= "Axis 1 (21.5% variation explained)", y= "Axis 2 (13.4% variation explained)") +
  scale_colour_manual(values = nestD_matrix_timepoint_ITS_palette) +
  scale_fill_manual(values = nestD_matrix_timepoint_ITS_palette) +
  geom_point(size = 7, shape = 18) +
  stat_ellipse(geom = "polygon", aes(fill= matrix_timepoint_cat), alpha = 0.25, lty = 2, linewidth = 1) +
  theme(legend.position = "none",
        panel.border = element_rect(linewidth = 1, colour = "black"),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 12, colour = "black"))

# NEST E ####
nestE_matrix_timepoint_palette <- c("#e14f00","#ff7800","#ffbb7f","#522700","#256602","#7bcf37","#44880d","#a9f46d","#183500","#008078","#44d6cd","#009caa","#8bf8f2","#003f3b","#8e4da3","#4e205d")
nestE_matrix_timepoint_ITS_palette <- c("#e14f00","#ff7800","#ffbb7f","#522700","#256602","#7bcf37","#44880d","#a9f46d","#183500","#008078","#44d6cd","#009caa","#8bf8f2","#003f3b","#3b0042","#8e4da3","#4e205d")

microbiome_nestE.css <- subset_samples(microbiome.css, nest_cat=="E")
microbiome_nestE.css <- prune_taxa(taxa_sums(microbiome_nestE.css) > 0, microbiome_nestE.css)
ITS_nestE.css <- subset_samples(ITS.css, nest_cat=="E")
ITS_nestE.css <- prune_taxa(taxa_sums(ITS_nestE.css) > 0, ITS_nestE.css)

microbiome_nestE.gunifrac <- gunifrac(microbiome_nestE.css)
ITS_nestE.gunifrac <- gunifrac(ITS_nestE.css)

microbiome_nestE.ord <- ordinate(microbiome_nestE.css, method = "PCoA", distance = microbiome_nestE.gunifrac)
ITS_nestE.ord <- ordinate(ITS_nestE.css, method = "PCoA", distance = ITS_nestE.gunifrac)

plot_ordination(microbiome_nestE.css, microbiome_nestE.ord, color = "matrix_timepoint_cat") +
  theme_bw() +
  labs(x= "Axis 1 (29.9% variation explained)", y= "Axis 2 (20.7% variation explained)") +
  scale_colour_manual(values = nestE_matrix_timepoint_ITS_palette) +
  scale_fill_manual(values = nestE_matrix_timepoint_ITS_palette) +
  geom_point(size = 7, shape = 18) +
  stat_ellipse(geom = "polygon", aes(fill= matrix_timepoint_cat), alpha = 0.25, lty = 2, linewidth = 1) +
  theme(legend.position = "none",
        panel.border = element_rect(linewidth = 1, colour = "black"),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 12, colour = "black"))

plot_ordination(ITS_nestE.css, ITS_nestE.ord, color = "matrix_timepoint_cat") +
  theme_bw() +
  labs(x= "Axis 1 (23.1% variation explained)", y= "Axis 2 (14.8% variation explained)") +
  scale_colour_manual(values = nestE_matrix_timepoint_ITS_palette) +
  scale_fill_manual(values = nestE_matrix_timepoint_ITS_palette) +
  geom_point(size = 7, shape = 18) +
  stat_ellipse(geom = "polygon", aes(fill= matrix_timepoint_cat), alpha = 0.25, lty = 2, linewidth = 1) +
  theme(legend.position = "none",
        panel.border = element_rect(linewidth = 1, colour = "black"),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 12, colour = "black"))

# NEST F ####
nestF_matrix_timepoint_palette <- c("#e14f00","#ff7800","#ffbb7f","#522700","#256602","#7bcf37","#44880d","#a9f46d","#183500","#008078","#44d6cd","#009caa","#8bf8f2","#003f3b","#8e4da3","#4e205d")
nestF_matrix_timepoint_ITS_palette <- c("#e14f00","#ff7800","#ffbb7f","#522700","#256602","#7bcf37","#44880d","#a9f46d","#183500","#008078","#44d6cd","#009caa","#8bf8f2","#003f3b","#3b0042","#8e4da3","#4e205d")

microbiome_nestF.css <- subset_samples(microbiome.css, nest_cat=="F")
microbiome_nestF.css <- prune_taxa(taxa_sums(microbiome_nestF.css) > 0, microbiome_nestF.css)
ITS_nestF.css <- subset_samples(ITS.css, nest_cat=="F")
ITS_nestF.css <- prune_taxa(taxa_sums(ITS_nestF.css) > 0, ITS_nestF.css)

microbiome_nestF.gunifrac <- gunifrac(microbiome_nestF.css)
ITS_nestF.gunifrac <- gunifrac(ITS_nestF.css)

microbiome_nestF.ord <- ordinate(microbiome_nestF.css, method = "PCoA", distance = microbiome_nestF.gunifrac)
ITS_nestF.ord <- ordinate(ITS_nestF.css, method = "PCoA", distance = ITS_nestF.gunifrac)

plot_ordination(microbiome_nestF.css, microbiome_nestF.ord, color = "matrix_timepoint_cat") +
  theme_bw() +
  labs(x= "Axis 1 (31.2% variation explained)", y= "Axis 2 (19.8% variation explained)") +
  scale_colour_manual(values = nestF_matrix_timepoint_ITS_palette) +
  scale_fill_manual(values = nestF_matrix_timepoint_ITS_palette) +
  geom_point(size = 7, shape = 18) +
  stat_ellipse(geom = "polygon", aes(fill= matrix_timepoint_cat), alpha = 0.25, lty = 2, linewidth = 1) +
  theme(legend.position = "none",
        panel.border = element_rect(linewidth = 1, colour = "black"),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 12, colour = "black"))

plot_ordination(ITS_nestF.css, ITS_nestF.ord, color = "matrix_timepoint_cat") +
  theme_bw() +
  labs(x= "Axis 1 (23.9% variation explained)", y= "Axis 2 (19.2% variation explained)") +
  scale_colour_manual(values = nestF_matrix_timepoint_ITS_palette) +
  scale_fill_manual(values = nestF_matrix_timepoint_ITS_palette) +
  geom_point(size = 7, shape = 18) +
  stat_ellipse(geom = "polygon", aes(fill= matrix_timepoint_cat), alpha = 0.25, lty = 2, linewidth = 1) +
  theme(legend.position = "none",
        panel.border = element_rect(linewidth = 1, colour = "black"),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 12, colour = "black"))

# NEST G ####
nestG_matrix_timepoint_palette <- c("#e14f00","#ff7800","#ffbb7f","#522700","#256602","#7bcf37","#44880d","#a9f46d","#183500","#008078","#44d6cd","#009caa","#8bf8f2","#003f3b","#8e4da3","#4e205d")
nestG_matrix_timepoint_ITS_palette <- c("#e14f00","#ff7800","#ffbb7f","#522700","#256602","#7bcf37","#44880d","#a9f46d","#183500","#008078","#44d6cd","#009caa","#8bf8f2","#003f3b","#3b0042","#8e4da3","#4e205d")

microbiome_nestG.css <- subset_samples(microbiome.css, nest_cat=="G")
microbiome_nestG.css <- prune_taxa(taxa_sums(microbiome_nestG.css) > 0, microbiome_nestG.css)
ITS_nestG.css <- subset_samples(ITS.css, nest_cat=="G")
ITS_nestG.css <- prune_taxa(taxa_sums(ITS_nestG.css) > 0, ITS_nestG.css)

microbiome_nestG.gunifrac <- gunifrac(microbiome_nestG.css)
ITS_nestG.gunifrac <- gunifrac(ITS_nestG.css)

microbiome_nestG.ord <- ordinate(microbiome_nestG.css, method = "PCoA", distance = microbiome_nestG.gunifrac)
ITS_nestG.ord <- ordinate(ITS_nestG.css, method = "PCoA", distance = ITS_nestG.gunifrac)

plot_ordination(microbiome_nestG.css, microbiome_nestG.ord, color = "matrix_timepoint_cat") +
  theme_bw() +
  labs(x= "Axis 1 (33.6% variation explained)", y= "Axis 2 (21.4% variation explained)") +
  scale_colour_manual(values = nestG_matrix_timepoint_ITS_palette) +
  scale_fill_manual(values = nestG_matrix_timepoint_ITS_palette) +
  geom_point(size = 7, shape = 18) +
  stat_ellipse(geom = "polygon", aes(fill= matrix_timepoint_cat), alpha = 0.25, lty = 2, linewidth = 1) +
  theme(legend.position = "none",
        panel.border = element_rect(linewidth = 1, colour = "black"),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 12, colour = "black"))

plot_ordination(ITS_nestG.css, ITS_nestG.ord, color = "matrix_timepoint_cat") +
  theme_bw() +
  labs(x= "Axis 1 (20.9% variation explained)", y= "Axis 2 (11.4% variation explained)") +
  scale_colour_manual(values = nestG_matrix_timepoint_ITS_palette) +
  scale_fill_manual(values = nestG_matrix_timepoint_ITS_palette) +
  geom_point(size = 7, shape = 18) +
  stat_ellipse(geom = "polygon", aes(fill= matrix_timepoint_cat), alpha = 0.25, lty = 2, linewidth = 1) +
  theme(#legend.position = "none",
        panel.border = element_rect(linewidth = 1, colour = "black"),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 12, colour = "black"))

#### let's take a look at with all nests combined into one plot
#### I need to do this to run stats between matrix_timepoint stratified on nest anyway
## 16S
all_nests_microbiome.css <- subset_samples(microbiome.css, nest_cat=="A" |
                                  nest_cat=="B" |
                                  nest_cat=="C" |
                                  nest_cat=="D" |
                                  nest_cat=="E" |
                                  nest_cat=="F" |
                                  nest_cat=="G")
all_nests_microbiome.css <- subset_samples(all_nests_microbiome.css, 
                                           matrix_timepoint_cat!="unhatched_egg-hatching of egg" &
                                           matrix_timepoint_cat!="unatched_egg-hatching of egg")
all_nests_microbiome.css <- prune_taxa(taxa_sums(all_nests_microbiome.css) > 0, all_nests_microbiome.css)
all_nests_microbiome.css.df <- as(all_nests_microbiome.css@sam_data,"data.frame")
all_nests_microbiome.gunifrac <- gunifrac(all_nests_microbiome.css)
all_nests_microbiome.ord <- ordinate(all_nests_microbiome.css, method = "PCoA", distance = all_nests_microbiome.gunifrac)
all_nests_palette <- c("#e14f00","#ff7800","#ffbb7f","#522700","#256602","#7bcf37","#44880d","#a9f46d","#183500","#008078","#44d6cd","#009caa","#8bf8f2","#003f3b","#3b0042","#8e4da3","#4e205d")

plot_ordination(all_nests_microbiome.css, all_nests_microbiome.ord, color = "matrix_timepoint_cat") +
  theme_bw() +
  labs(x= "Axis 1 (25.9% variation explained)", y= "Axis 2 (18.2% variation explained)") +
  scale_colour_manual(values = all_nests_palette) +
  scale_fill_manual(values = all_nests_palette) +
  geom_point(size = 7, shape = 18) +
  stat_ellipse(geom = "polygon", aes(fill= matrix_timepoint_cat), alpha = 0.25, lty = 2, linewidth = 1) +
  theme(legend.position = "none",
    panel.border = element_rect(linewidth = 1, colour = "black"),
    axis.title = element_text(size = 24),
    axis.text = element_text(size = 12, colour = "black"))
# ya, its chaotic. Better when split by nest, and then with the stats run with stratified on nest

# stats
all_nests_microbiome_adonis <- pairwise.adonis2(all_nests_microbiome.gunifrac ~ matrix_timepoint_cat, data = all_nests_microbiome.css.df, strata = 'nest_cat', nperm = 9999, p.adjust.methods = "BH")
all_nests_microbiome_adonis

### ITS
all_nests_ITS.css <- subset_samples(ITS.css, nest_cat=="A" |
                                             nest_cat=="B" |
                                             nest_cat=="C" |
                                             nest_cat=="D" |
                                             nest_cat=="E" |
                                             nest_cat=="F" |
                                             nest_cat=="G")
all_nests_ITS.css <- subset_samples(all_nests_ITS.css, 
                                           matrix_timepoint_cat!="unhatched_egg-hatching of egg" &
                                             matrix_timepoint_cat!="unatched_egg-hatching of egg")
all_nests_ITS.css <- prune_taxa(taxa_sums(all_nests_ITS.css) > 0, all_nests_ITS.css)
all_nests_ITS.css.df <- as(all_nests_ITS.css@sam_data,"data.frame")
all_nests_ITS.gunifrac <- gunifrac(all_nests_ITS.css)

#stats
all_nests_ITS_adonis <- pairwise.adonis2(all_nests_ITS.gunifrac ~ matrix_timepoint_cat, data = all_nests_ITS.css.df, strata = 'nest_cat', nperm = 9999, p.adjust.methods = "BH")
all_nests_ITS_adonis

############################################ TURTLE SAMPLES ####
washed_palette <- c("darkgoldenrod3","dodgerblue4")
### 16S
turtles_microbiome.css <- subset_samples(microbiome.css, matrix=="turtle")
turtles_microbiome.css <- prune_taxa(taxa_sums(turtles_microbiome.css) > 0, turtles_microbiome.css)
turtles_microbiome.gunifrac <- gunifrac(turtles_microbiome.css)
turtles_microbiome.ord <- ordinate(turtles_microbiome.css, method = "PCoA", distance = turtles_microbiome.gunifrac)

plot_ordination(turtles_microbiome.css, turtles_microbiome.ord, color = "washed_unwashed") +
  labs(x= "Axis 1 (29.6% variation explained)", y= "Axis 2 (11.8% variation explained)") +
  facet_wrap(~factor(matrix_timepoint_cat, levels = c("turtle_carapace-hatching of egg","turtle_carapace-3 months post hatch","turtle_carapace-6 months post hatch","turtle_carapace-1 month post move outside","turtle_carapace-release to wild",
                                                      "turtle_plastron-hatching of egg","turtle_plastron-3 months post hatch","turtle_plastron-6 months post hatch","turtle_plastron-1 month post move outside","turtle_plastron-release to wild")), nrow = 2) + theme_bw() +
  stat_ellipse(geom = "polygon", alpha = 0.3, aes(fill= washed_unwashed), lty = 2, linewidth = 1) +
  geom_point(size = 6.5) +
  geom_text(aes(label = nest_cat), colour = "white", fontface = "bold") +
  scale_colour_manual(values = washed_palette) +
  scale_fill_manual(values = washed_palette) +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.border = element_rect(linewidth = 1, colour = "black"),
        axis.title = element_text(size = 32),
        axis.text = element_text(size = 12, colour = "black"))

# plastron
plastron_microbiome.css <- subset_samples(turtles_microbiome.css, matrix_specific=="turtle_plastron")
plastron_microbiome.css <- prune_taxa(taxa_sums(plastron_microbiome.css) > 0, plastron_microbiome.css)
plastron_microbiome.gunifrac <- gunifrac(plastron_microbiome.css)
plastron_microbiome.df <- as(plastron_microbiome.css@sam_data, "data.frame")

adonis2(plastron_microbiome.gunifrac ~ washed_unwashed, plastron_microbiome.df, strata = 'nest')

# carapace
plastron_microbiome.css <- subset_samples(turtles_microbiome.css, matrix_specific=="turtle_plastron")
plastron_microbiome.css <- prune_taxa(taxa_sums(plastron_microbiome.css) > 0, plastron_microbiome.css)

### ITS
turtles_ITS.css <- subset_samples(ITS.css, matrix=="turtle")
turtles_ITS.css <- prune_taxa(taxa_sums(turtles_ITS.css) > 0, turtles_ITS.css)
turtles_ITS.gunifrac <- gunifrac(turtles_ITS.css)
turtles_ITS.ord <- ordinate(turtles_ITS.css, method = "PCoA", distance = turtles_ITS.gunifrac)

plot_ordination(turtles_ITS.css, turtles_ITS.ord, color = "washed_unwashed") +
  labs(x= "Axis 1 (20.3% variation explained)", y= "Axis 2 (9.6% variation explained)") +
  facet_wrap(~factor(matrix_timepoint_cat, levels = c("turtle_carapace-hatching of egg","turtle_carapace-3 months post hatch","turtle_carapace-6 months post hatch","turtle_carapace-1 month post move outside","turtle_carapace-release to wild",
                                                      "turtle_plastron-hatching of egg","turtle_plastron-3 months post hatch","turtle_plastron-6 months post hatch","turtle_plastron-1 month post move outside","turtle_plastron-release to wild")), nrow = 2) + theme_bw() +
  stat_ellipse(geom = "polygon", alpha = 0.3, aes(fill= washed_unwashed), lty = 2, linewidth = 1) +
  geom_point(size = 6.5) +
  geom_text(aes(label = nest_cat), colour = "white", fontface = "bold") +
  scale_colour_manual(values = washed_palette) +
  scale_fill_manual(values = washed_palette) +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.border = element_rect(linewidth = 1, colour = "black"),
        axis.title = element_text(size = 32),
        axis.text = element_text(size = 12, colour = "black"))

############################################ EGG SAMPLES ####
#################################### 16S ####
#################### ALPHA DIV ####
eggs_microbiome_alpha <- alpha_div_16S_meta[which(alpha_div_16S_meta$matrix=="egg"),]
#### TEMPORAL ####
ggplot(eggs_microbiome_alpha, aes(x= Timepoint_Cat, y= Observed, fill = Timepoint_Cat, colour = Timepoint_Cat)) +
  theme_bw() +
  labs(y= "Observed ASVs") +
  geom_boxplot(alpha = 0.5, linewidth = 1) + geom_point(size = 4, shape = 18) +
  scale_x_discrete(limits = c("egg in nest","arrival at WPZ","1 month incubation")) +
  scale_y_continuous(expand = c(0.001,0,0.06,0)) +
  scale_fill_manual(values = egg_timepoint_palette) +
  scale_colour_manual(values = egg_timepoint_palette) +
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", linewidth = 1),
        panel.grid.major.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 36),
        axis.text.y = element_text(size = 16, colour = "black"),
        axis.ticks.y = element_line(colour = "black", size = 0.75))

pairwise.wilcox.test(eggs_microbiome_alpha$Observed, eggs_microbiome_alpha$Timepoint_Cat, p.adjust.method = "BH")

ggplot(eggs_microbiome_alpha, aes(x= Timepoint_Cat, y= Shannon, fill = Timepoint_Cat, colour = Timepoint_Cat)) +
  theme_bw() +
  labs(y= "Shannon Index") +
  geom_boxplot(alpha = 0.5, linewidth = 1) + geom_point(size = 4, shape = 18) +
  scale_x_discrete(limits = c("egg in nest","arrival at WPZ","1 month incubation")) +
  scale_y_continuous(expand = c(0.001,0,0.06,0)) +
  scale_fill_manual(values = egg_timepoint_palette) +
  scale_colour_manual(values = egg_timepoint_palette) +
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", linewidth = 1),
        panel.grid.major.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 36),
        axis.text.y = element_text(size = 16, colour = "black"),
        axis.ticks.y = element_line(colour = "black", size = 0.75))

pairwise.wilcox.test(eggs_microbiome_alpha$Shannon, eggs_microbiome_alpha$Timepoint_Cat, p.adjust.method = "BH")

ggplot(eggs_microbiome_alpha, aes(x= Timepoint_Cat, y= PD, fill = Timepoint_Cat, colour = Timepoint_Cat)) +
  theme_bw() +
  labs(y= "Faith's PD") +
  geom_boxplot(alpha = 0.5, linewidth = 1) + geom_point(size = 4, shape = 18) +
  scale_x_discrete(limits = c("egg in nest","arrival at WPZ","1 month incubation")) +
  scale_y_continuous(expand = c(0.001,0,0.06,0)) +
  scale_fill_manual(values = egg_timepoint_palette) +
  scale_colour_manual(values = egg_timepoint_palette) +
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", linewidth = 1),
        panel.grid.major.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 36),
        axis.text.y = element_text(size = 16, colour = "black"),
        axis.ticks.y = element_line(colour = "black", size = 0.75))

pairwise.wilcox.test(eggs_microbiome_alpha$PD, eggs_microbiome_alpha$Timepoint_Cat, p.adjust.method = "BH")

#### WASHED V UNWASHED ####
eggs_microbiome_alpha_onemonth <- eggs_microbiome_alpha[which(eggs_microbiome_alpha$Timepoint_Cat=="1 month incubation"),]

ggplot(eggs_microbiome_alpha_onemonth, aes(x= washed_unwashed, y= Observed, fill = washed_unwashed, colour = washed_unwashed)) +
  theme_bw() +
  labs(y= "Observed ASVs") +
  geom_boxplot(alpha = 0.5, linewidth = 1) + geom_point(size = 4, shape = 18) +
  scale_y_continuous(expand = c(0.040,0,0.06,0)) +
  scale_fill_manual(values = washed_palette) +
  scale_colour_manual(values = washed_palette) +
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", linewidth = 1),
        panel.grid.major.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 36),
        axis.text.y = element_text(size = 16, colour = "black"),
        axis.ticks.y = element_line(colour = "black", size = 0.75))

kruskal.test(eggs_microbiome_alpha_onemonth$Observed, eggs_microbiome_alpha_onemonth$washed_unwashed) #sig

ggplot(eggs_microbiome_alpha_onemonth, aes(x= washed_unwashed, y= Shannon, fill = washed_unwashed, colour = washed_unwashed)) +
  theme_bw() +
  labs(y= "Shannon Index") +
  geom_boxplot(alpha = 0.5, linewidth = 1) + geom_point(size = 4, shape = 18) +
  scale_y_continuous(expand = c(0.040,0,0.06,0)) +
  scale_fill_manual(values = washed_palette) +
  scale_colour_manual(values = washed_palette) +
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", linewidth = 1),
        panel.grid.major.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 36),
        axis.text.y = element_text(size = 16, colour = "black"),
        axis.ticks.y = element_line(colour = "black", size = 0.75))

kruskal.test(eggs_microbiome_alpha_onemonth$Shannon, eggs_microbiome_alpha_onemonth$washed_unwashed) #sig

ggplot(eggs_microbiome_alpha_onemonth, aes(x= washed_unwashed, y= PD, fill = washed_unwashed, colour = washed_unwashed)) +
  theme_bw() +
  labs(y= "Faith's PD") +
  geom_boxplot(alpha = 0.5, linewidth = 1) + geom_point(size = 4, shape = 18) +
  scale_y_continuous(expand = c(0.040,0,0.06,0)) +
  scale_fill_manual(values = washed_palette) +
  scale_colour_manual(values = washed_palette) +
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", linewidth = 1),
        panel.grid.major.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 36),
        axis.text.y = element_text(size = 16, colour = "black"),
        axis.ticks.y = element_line(colour = "black", size = 0.75))

kruskal.test(eggs_microbiome_alpha_onemonth$PD, eggs_microbiome_alpha_onemonth$washed_unwashed) # NS (p =)

#################### BETA DIV ####
eggs_microbiome.css <- subset_samples(microbiome.css, matrix=="egg")
eggs_microbiome.css <- prune_taxa(taxa_sums(eggs_microbiome.css) > 0, eggs_microbiome.css)
eggs_microbiome.gunifrac <- gunifrac(eggs_microbiome.css)
#### PCoA washed v unwashed ####
eggs_microbiome.ord <- ordinate(eggs_microbiome.css, method = "PCoA", distance = eggs_microbiome.gunifrac)
plot_ordination(eggs_microbiome.css, eggs_microbiome.ord, color = "washed_unwashed") +
  labs(x= "Axis 1 (37.8% variation explained)", y= "Axis 2 (11.5% variation explained)") +
  facet_wrap(~factor(matrix_timepoint_cat, levels = c("egg-egg in nest","egg-arrival at WPZ","egg-1 month incubation")), nrow = 1) + theme_bw() +
  stat_ellipse(geom = "polygon", alpha = 0.25, aes(fill= washed_unwashed), lty = 2, linewidth = 1) +
  geom_point(size = 6.5) +
  geom_text(aes(label = nest_cat), colour = "white", fontface = "bold") +
  scale_colour_manual(values = washed_palette) +
  scale_fill_manual(values = washed_palette) +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.border = element_rect(linewidth = 1, colour = "black"),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 12, colour = "black"))

#### DENDROGRAM TIMEPOINTS ####
egg_timepoint_palette <- c("#e14f00","#ff7800","#ffbb7f")
eggs_microbiome.hclust <- hclust(eggs_microbiome.gunifrac, method = "ward.D2")
eggs_microbiome.dendro <- as.dendrogram(eggs_microbiome.hclust)
eggs_microbiome.dendro.data <- dendro_data(eggs_microbiome.dendro, type = "rectangle")
eggs_microbiome.metadata_for_dendro <- as_tibble(eggs_microbiome.css@sam_data)
eggs_microbiome.dendro.data$labels <- eggs_microbiome.dendro.data$labels %>%
  left_join(eggs_microbiome.metadata_for_dendro, by = c("label" = "SampleID"))

ggplot(eggs_microbiome.dendro.data$segments) +
  theme_minimal() +
  labs(y= "Ward's Distance") +
  geom_segment(aes(x=x,y=y,xend=xend,yend=yend)) +
  geom_point(data = eggs_microbiome.dendro.data$labels, 
             aes(x=x,y=y, colour = Timepoint_Cat, fill= Timepoint_Cat, shape = washed_unwashed, size = washed_unwashed),
             stroke =1.5, position = position_nudge(y=-0.12)) +
  geom_text(data = eggs_microbiome.dendro.data$labels, aes(x=x,y=y, label = nest_cat), 
            position = position_nudge(y=-0.11), fontface = "bold", size = 5) +
  scale_shape_manual(values = c(18,20)) +
  scale_size_manual(values = c(10,11)) +
  scale_x_discrete(expand = c(0.01,0,0.01,0)) +
  scale_fill_manual(values = egg_timepoint_palette) +
  scale_colour_manual(values = egg_timepoint_palette) +
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

#### FAMILY RA FOR DENDRO ####
eggs_microbiome_ra <- transform_sample_counts(eggs_microbiome.css, function(x) {x/sum(x)} * 100)
eggs_microbiome_family_ra <- tax_glom(eggs_microbiome_ra, taxrank = "Family", NArm = F)
eggs_microbiome_family_ra_filt <- merge_low_abundance(eggs_microbiome_family_ra, threshold = 0.2) # 76 taxa
eggs_microbiome_family_palette <- distinctColorPalette(76)
write.csv(otu_table(eggs_microbiome_family_ra_filt),"eggs_family_otus.csv")
write.csv(tax_table(eggs_microbiome_family_ra_filt),"eggs_family_taxa.csv")


eggs_microbiome_family_ra_filt_melt <- psmelt(eggs_microbiome_family_ra_filt)

eggs_microbiome_sample_order <- eggs_microbiome.dendro.data$labels$label
egg_microbiome_family_palette <- c("#A53DF2",	"#988CE9",	"#ECA148",	"#EE97DD",	"#60ED5B",	"#F6C2EC",	"#DAA28F",	"#E531BE",	"seagreen3",	"#9DF4EE",	"#96AF59",	"#F45EE9",	"#E1BB96",	"#A0F2D4",	"#C5A0CF",	"#85CDBE",	"#59A1C6",	"#CAEC78",	"#A2D0D6",	"#528DEA",	"#82ABDE",	"#8EC6E9",	"#E695C7",	"#86E2A2",	"#DB74EF",	"#E8E631",	"#A532CE",	"#B9CE45",	"#E6CFF5",	"#B086A1",	"#B56063",	"#DBE376",	"#6470E5",	"#A45BEA",	"#EF82EA",	"#CF4BE1",	"#90BAAB",	"#45CDDD",	"#EAECB9",	"#BA37A5",	"#F34175",	"#E4A6C7",	"#AD9CE1",	"#66B636",	"#F03846",	"#4A45EC",	"#7FF434",	"#6EBB8B",	"#4EA7EF",	"#EE7736",	"#3BE67C",	"#F4F5E5",	"#6E20DC",	"#ED6376",	"#7FD2DB",	"#AFBACF",	"mediumorchid3",	"#D8F1B0",	"#7CF295",	"#CA84E4",	"#EA8DA4",	"#506268",	"#ED936E",	"#BC417F",	"#518097",	"#B94F35",	"#B8BAB2",	"#D4F1F2",	"#457747",	"#AFF33A",	"#48B5BA",	"#F3D9EC",	"#65F1BD",	"#83EF73",	"#EAD7C3",	"grey88")

ggplot(eggs_microbiome_family_ra_filt_melt, aes(x= SampleID, y= Abundance, fill = Family)) +
  theme_minimal() +
  labs(y= "Relative Abundance (%)") +
  geom_bar(stat = "identity", colour = "black") +
  scale_fill_manual(values = egg_microbiome_family_palette) +
  scale_x_discrete(limits = eggs_microbiome_sample_order, expand = c(0.009,0,0.015,0)) +
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

#### ABUNDANT FAMILIES AFTER TIME AT WPZ ####
## Streptomycetaceae ####
eggs_Streptomycetaceae_family <- eggs_microbiome_family_ra_filt_melt[which(eggs_microbiome_family_ra_filt_melt$Family=="Streptomycetaceae"),]

ggplot(eggs_Streptomycetaceae_family, aes(x= Timepoint_Cat, y= Abundance, fill = Timepoint_Cat, colour = Timepoint_Cat)) +
  theme_bw() +
  labs(y= "Streptomycetaceae RA (%)") +
  geom_boxplot(alpha = 0.5, linewidth = 1) + geom_point(size = 4, shape = 18) +
  scale_x_discrete(limits = c("egg in nest","arrival at WPZ","1 month incubation")) +
  scale_y_continuous(expand = c(0.001,0,0.06,0)) +
  scale_fill_manual(values = egg_timepoint_palette) +
  scale_colour_manual(values = egg_timepoint_palette) +
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", linewidth = 1),
        panel.grid.major.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 36),
        axis.text.y = element_text(size = 16, colour = "black"),
        axis.ticks.y = element_line(colour = "black", size = 0.75))

pairwise.wilcox.test(eggs_Streptomycetaceae_family$Abundance, eggs_Streptomycetaceae_family$Timepoint_Cat) # sig. more in 1 month (p < 0.001)

# washed v unwashed
eggs_Streptomycetaceae_family_oneMonth <- eggs_Streptomycetaceae_family[which(eggs_Streptomycetaceae_family$Timepoint_Cat=="1 month incubation"),]

ggplot(eggs_Streptomycetaceae_family_oneMonth, aes(x= washed_unwashed, y= Abundance, fill = washed_unwashed, colour = washed_unwashed)) +
  theme_bw() +
  labs(y= "Streptomycetaceae RA (%)") +
  geom_boxplot(alpha = 0.5, linewidth = 1) + geom_point(size = 4, shape = 18) +
  scale_y_continuous(expand = c(0.040,0,0.06,0)) +
  scale_fill_manual(values = washed_palette) +
  scale_colour_manual(values = washed_palette) +
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", linewidth = 1),
        panel.grid.major.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 36),
        axis.text.y = element_text(size = 16, colour = "black"),
        axis.ticks.y = element_line(colour = "black", size = 0.75))

kruskal.test(eggs_Streptomycetaceae_family_oneMonth$Abundance, eggs_Streptomycetaceae_family_oneMonth$washed_unwashed) # sig. p = 0.007

## Pseudonocardiaceae ####
eggs_Pseudonocardiaceae_family <- eggs_microbiome_family_ra_filt_melt[which(eggs_microbiome_family_ra_filt_melt$Family=="Pseudonocardiaceae"),]

ggplot(eggs_Pseudonocardiaceae_family, aes(x= Timepoint_Cat, y= Abundance, fill = Timepoint_Cat, colour = Timepoint_Cat)) +
  theme_bw() +
  labs(y= "Pseudonocardiaceae RA (%)") +
  geom_boxplot(alpha = 0.5, linewidth = 1) + geom_point(size = 4, shape = 18) +
  scale_x_discrete(limits = c("egg in nest","arrival at WPZ","1 month incubation")) +
  scale_y_continuous(expand = c(0.001,0,0.69,0)) +
  scale_fill_manual(values = egg_timepoint_palette) +
  scale_colour_manual(values = egg_timepoint_palette) +
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", linewidth = 1),
        panel.grid.major.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 36),
        axis.text.y = element_text(size = 16, colour = "black"),
        axis.ticks.y = element_line(colour = "black", size = 0.75))

pairwise.wilcox.test(eggs_Pseudonocardiaceae_family$Abundance, eggs_Pseudonocardiaceae_family$Timepoint_Cat) # sig. more in 1 month (p < 0.001)

# washed v unwashed
eggs_Pseudonocardiaceae_family_oneMonth <- eggs_Pseudonocardiaceae_family[which(eggs_Pseudonocardiaceae_family$Timepoint_Cat=="1 month incubation"),]

ggplot(eggs_Pseudonocardiaceae_family_oneMonth, aes(x= washed_unwashed, y= Abundance, fill = washed_unwashed, colour = washed_unwashed)) +
  theme_bw() +
  labs(y= "Pseudonocardiaceae RA (%)") +
  geom_boxplot(alpha = 0.5, linewidth = 1) + geom_point(size = 4, shape = 18) +
  scale_y_continuous(expand = c(0.001,0,0.69,0)) +
  scale_fill_manual(values = washed_palette) +
  scale_colour_manual(values = washed_palette) +
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", linewidth = 1),
        panel.grid.major.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 36),
        axis.text.y = element_text(size = 16, colour = "black"),
        axis.ticks.y = element_line(colour = "black", size = 0.75))

kruskal.test(eggs_Pseudonocardiaceae_family_oneMonth$Abundance, eggs_Pseudonocardiaceae_family_oneMonth$washed_unwashed) # sig. p < 0.001


#################################### ITS ####
eggs_ITS.css <- subset_samples(ITS.css, matrix=="egg")

#################### ALPHA DIV ####
eggs_ITS_alpha <- alpha_div_ITS_meta[which(alpha_div_ITS_meta$matrix=="egg"),]
#### TEMPORAL ####
ggplot(eggs_ITS_alpha, aes(x= Timepoint_Cat, y= Observed, fill = Timepoint_Cat, colour = Timepoint_Cat)) +
  theme_bw() +
  labs(y= "Observed ASVs") +
  geom_boxplot(alpha = 0.5, linewidth = 1) + geom_point(size = 4, shape = 18) +
  scale_x_discrete(limits = c("egg in nest","arrival at WPZ","1 month incubation")) +
  scale_y_continuous(expand = c(0.001,0,0.06,0)) +
  scale_fill_manual(values = egg_timepoint_palette) +
  scale_colour_manual(values = egg_timepoint_palette) +
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", linewidth = 1),
        panel.grid.major.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 36),
        axis.text.y = element_text(size = 16, colour = "black"),
        axis.ticks.y = element_line(colour = "black", size = 0.75))

pairwise.wilcox.test(eggs_ITS_alpha$Observed, eggs_ITS_alpha$Timepoint_Cat, p.adjust.method = "BH")

ggplot(eggs_ITS_alpha, aes(x= Timepoint_Cat, y= Shannon, fill = Timepoint_Cat, colour = Timepoint_Cat)) +
  theme_bw() +
  labs(y= "Shannon Index") +
  geom_boxplot(alpha = 0.5, linewidth = 1) + geom_point(size = 4, shape = 18) +
  scale_x_discrete(limits = c("egg in nest","arrival at WPZ","1 month incubation")) +
  scale_y_continuous(expand = c(0.001,0,0.06,0)) +
  scale_fill_manual(values = egg_timepoint_palette) +
  scale_colour_manual(values = egg_timepoint_palette) +
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", linewidth = 1),
        panel.grid.major.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 36),
        axis.text.y = element_text(size = 16, colour = "black"),
        axis.ticks.y = element_line(colour = "black", size = 0.75))

pairwise.wilcox.test(eggs_ITS_alpha$Shannon, eggs_ITS_alpha$Timepoint_Cat, p.adjust.method = "BH")

ggplot(eggs_ITS_alpha, aes(x= Timepoint_Cat, y= PD, fill = Timepoint_Cat, colour = Timepoint_Cat)) +
  theme_bw() +
  labs(y= "Faith's PD") +
  geom_boxplot(alpha = 0.5, linewidth = 1) + geom_point(size = 4, shape = 18) +
  scale_x_discrete(limits = c("egg in nest","arrival at WPZ","1 month incubation")) +
  scale_y_continuous(expand = c(0.001,0,0.06,0)) +
  scale_fill_manual(values = egg_timepoint_palette) +
  scale_colour_manual(values = egg_timepoint_palette) +
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", linewidth = 1),
        panel.grid.major.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 36),
        axis.text.y = element_text(size = 16, colour = "black"),
        axis.ticks.y = element_line(colour = "black", size = 0.75))

pairwise.wilcox.test(eggs_ITS_alpha$PD, eggs_ITS_alpha$Timepoint_Cat, p.adjust.method = "BH")

#### WASHED v UNWASHED ####
eggs_ITS_alpha_onemonth <- eggs_ITS_alpha[which(eggs_ITS_alpha$Timepoint_Cat=="1 month incubation"),]

ggplot(eggs_ITS_alpha_onemonth, aes(x= washed_unwashed, y= Observed, fill = washed_unwashed, colour = washed_unwashed)) +
  theme_bw() +
  labs(y= "Observed ASVs") +
  geom_boxplot(alpha = 0.5, linewidth = 1) + geom_point(size = 4, shape = 18) +
  scale_y_continuous(expand = c(0.040,0,0.06,0)) +
  scale_fill_manual(values = washed_palette) +
  scale_colour_manual(values = washed_palette) +
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", linewidth = 1),
        panel.grid.major.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 36),
        axis.text.y = element_text(size = 16, colour = "black"),
        axis.ticks.y = element_line(colour = "black", size = 0.75))

ggplot(eggs_ITS_alpha_onemonth, aes(x= washed_unwashed, y= Shannon, fill = washed_unwashed, colour = washed_unwashed)) +
  theme_bw() +
  labs(y= "Shannon Index") +
  geom_boxplot(alpha = 0.5, linewidth = 1) + geom_point(size = 4, shape = 18) +
  scale_y_continuous(expand = c(0.040,0,0.06,0)) +
  scale_fill_manual(values = washed_palette) +
  scale_colour_manual(values = washed_palette) +
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", linewidth = 1),
        panel.grid.major.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 36),
        axis.text.y = element_text(size = 16, colour = "black"),
        axis.ticks.y = element_line(colour = "black", size = 0.75))

ggplot(eggs_ITS_alpha_onemonth, aes(x= washed_unwashed, y= PD, fill = washed_unwashed, colour = washed_unwashed)) +
  theme_bw() +
  labs(y= "Observed ASVs") +
  geom_boxplot(alpha = 0.5, linewidth = 1) + geom_point(size = 4, shape = 18) +
  scale_y_continuous(expand = c(0.040,0,0.06,0)) +
  scale_fill_manual(values = washed_palette) +
  scale_colour_manual(values = washed_palette) +
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", linewidth = 1),
        panel.grid.major.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 36),
        axis.text.y = element_text(size = 16, colour = "black"),
        axis.ticks.y = element_line(colour = "black", size = 0.75))

#################### BETA DIV ####
eggs_ITS.css <- prune_taxa(taxa_sums(eggs_ITS.css) > 0, eggs_ITS.css)
eggs_ITS.gunifrac <- gunifrac(eggs_ITS.css)
eggs_ITS.ord <- ordinate(eggs_ITS.css, method = "PCoA", distance = eggs_ITS.gunifrac)

#### PCoA washed v unwashed ####
plot_ordination(eggs_ITS.css, eggs_ITS.ord, color = "washed_unwashed") +
  labs(x= "Axis 1 (21.1% variation explained)", y= "Axis 2 (10.7% variation explained)") +
  facet_wrap(~factor(matrix_timepoint_cat, levels = c("egg-egg in nest","egg-arrival at WPZ","egg-1 month incubation")), nrow = 1) + theme_bw() +
  stat_ellipse(geom = "polygon", alpha = 0.25, aes(fill= washed_unwashed), lty = 2, linewidth = 1) +
  geom_point(size = 6.5) +
  geom_text(aes(label = nest_cat), colour = "white", fontface = "bold") +
  scale_colour_manual(values = washed_palette) +
  scale_fill_manual(values = washed_palette) +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.border = element_rect(linewidth = 1, colour = "black"),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 12, colour = "black"))

#### DENDROGRAM TIMEPOINTS ####
eggs_ITS.hclust <- hclust(eggs_ITS.gunifrac, method = "ward.D2")
eggs_ITS.dendro <- as.dendrogram(eggs_ITS.hclust)
eggs_ITS.dendro.data <- dendro_data(eggs_ITS.dendro, type = "rectangle")
eggs_ITS.metadata_for_dendro <- as_tibble(eggs_ITS.css@sam_data)
eggs_ITS.dendro.data$labels <- eggs_ITS.dendro.data$labels %>%
  left_join(eggs_ITS.metadata_for_dendro, by = c("label" = "SampleID"))

ggplot(eggs_ITS.dendro.data$segments) +
  theme_minimal() +
  labs(y= "Ward's Distance") +
  geom_segment(aes(x=x,y=y,xend=xend,yend=yend)) +
  geom_point(data = eggs_ITS.dendro.data$labels, 
             aes(x=x,y=y, colour = Timepoint_Cat, fill= Timepoint_Cat, shape = washed_unwashed, size = washed_unwashed),
             stroke =1.5, position = position_nudge(y=-0.12)) +
  geom_text(data = eggs_ITS.dendro.data$labels, aes(x=x,y=y, label = nest_cat), 
            position = position_nudge(y=-0.11), fontface = "bold", size = 5) +
  scale_shape_manual(values = c(18,20)) +
  scale_size_manual(values = c(10,11)) +
  scale_x_discrete(expand = c(0.01,0,0.01,0)) +
  scale_fill_manual(values = egg_timepoint_palette) +
  scale_colour_manual(values = egg_timepoint_palette) +
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

#### FAMILY RA FOR DENDRO ####
eggs_ITS_ra <- transform_sample_counts(eggs_ITS.css, function(x) {x/sum(x)} * 100)
eggs_ITS_family_ra <- tax_glom(eggs_ITS_ra, taxrank = "Family", NArm = F) #323 families
eggs_ITS_family_ra_filt <- merge_low_abundance(eggs_ITS_family_ra, threshold = 0.2) # 35 taxa
write.csv(otu_table(eggs_ITS_family_ra_filt),"eggs_family_otus.csv")
write.csv(tax_table(eggs_ITS_family_ra_filt),"eggs_family_taxa.csv")


eggs_ITS_family_ra_filt_melt <- psmelt(eggs_ITS_family_ra_filt)

eggs_ITS_sample_order <- eggs_ITS.dendro.data$labels$label
egg_ITS_family_palette <- c("#EEF5E1",	"#85D8EC",	"#E86078",	"#B89CD2",	"#DBBBCA",	"#51ECD3",	"#7E8FDD",	"#E8AEBD",	"#E7E631",	"#6C54CE",	"#83EFF0",	"#A24FD6",	"#D8F0AE",	"#C0F19E",	"#E1A6C8",	"#DCBFF3",	"#D5478D",	"#518FA6",	"#A9819A",	"#8CEDD5",	"#C9E7C7",	"#98D342",	"#EDEA96",	"#5A7EEB",	"#45C5E4",	"#395964",	"#E381F2",	"#F1B8E9",	"#B7AECA",	"#EB906A",	"#9094BC",	"#E435ED",	"#52815F",	"#96B5C8",	"grey88")

ggplot(eggs_ITS_family_ra_filt_melt, aes(x= SampleID, y= Abundance, fill = Family)) +
  theme_minimal() +
  labs(y= "Relative Abundance (%)") +
  geom_bar(stat = "identity", colour = "black") +
  scale_fill_manual(values = egg_ITS_family_palette) +
  scale_x_discrete(limits = eggs_ITS_sample_order, expand = c(0.009,0,0.015,0)) +
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

#### ABUNDANT FAMILIES AFTER TIME AT WPZ ####
## Aspergillaceae ####
eggs_Aspergillaceae_family <- eggs_ITS_family_ra_filt_melt[which(eggs_ITS_family_ra_filt_melt$Family=="Aspergillaceae"),]

ggplot(eggs_Aspergillaceae_family, aes(x= Timepoint_Cat, y= Abundance, fill = Timepoint_Cat, colour = Timepoint_Cat)) +
  theme_bw() +
  labs(y= "Aspergillaceae RA (%)") +
  geom_boxplot(alpha = 0.5, linewidth = 1) + geom_point(size = 4, shape = 18) +
  scale_x_discrete(limits = c("egg in nest","arrival at WPZ","1 month incubation")) +
  scale_y_continuous(expand = c(0.028,0,0.014,0)) +
  scale_fill_manual(values = egg_timepoint_palette) +
  scale_colour_manual(values = egg_timepoint_palette) +
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", linewidth = 1),
        panel.grid.major.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 36),
        axis.text.y = element_text(size = 16, colour = "black"),
        axis.ticks.y = element_line(colour = "black", size = 0.75))

pairwise.wilcox.test(eggs_Aspergillaceae_family$Abundance, eggs_Aspergillaceae_family$Timepoint_Cat) # sig. more in 1 month (p < 0.001)

# washed v unwashed
eggs_Aspergillaceae_family_oneMonth <- eggs_Aspergillaceae_family[which(eggs_Aspergillaceae_family$Timepoint_Cat=="1 month incubation"),]

ggplot(eggs_Aspergillaceae_family_oneMonth, aes(x= washed_unwashed, y= Abundance, fill = washed_unwashed, colour = washed_unwashed)) +
  theme_bw() +
  labs(y= "Aspergillaceae RA (%)") +
  geom_boxplot(alpha = 0.5, linewidth = 1) + geom_point(size = 4, shape = 18) +
  scale_y_continuous(expand = c(0.324,0,0.018,0)) +
  scale_fill_manual(values = washed_palette) +
  scale_colour_manual(values = washed_palette) +
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", linewidth = 1),
        panel.grid.major.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 36),
        axis.text.y = element_text(size = 16, colour = "black"),
        axis.ticks.y = element_line(colour = "black", size = 0.75))

kruskal.test(eggs_Aspergillaceae_family_oneMonth$Abundance, eggs_Aspergillaceae_family_oneMonth$washed_unwashed) # NS

## Didymellaceae ####
eggs_Didymellaceae_family <- eggs_ITS_family_ra_filt_melt[which(eggs_ITS_family_ra_filt_melt$Family=="Didymellaceae"),]

ggplot(eggs_Didymellaceae_family, aes(x= Timepoint_Cat, y= Abundance, fill = Timepoint_Cat, colour = Timepoint_Cat)) +
  theme_bw() +
  labs(y= "Didymellaceae RA (%)") +
  geom_boxplot(alpha = 0.5, linewidth = 1) + geom_point(size = 4, shape = 18) +
  scale_x_discrete(limits = c("egg in nest","arrival at WPZ","1 month incubation")) +
  scale_y_continuous(expand = c(0.001,0,0.612,0)) +
  scale_fill_manual(values = egg_timepoint_palette) +
  scale_colour_manual(values = egg_timepoint_palette) +
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", linewidth = 1),
        panel.grid.major.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 36),
        axis.text.y = element_text(size = 16, colour = "black"),
        axis.ticks.y = element_line(colour = "black", size = 0.75))

pairwise.wilcox.test(eggs_Didymellaceae_family$Abundance, eggs_Didymellaceae_family$Timepoint_Cat) # sig. more in 1 month (p < 0.001)

# washed v unwashed
eggs_Didymellaceae_family_oneMonth <- eggs_Didymellaceae_family[which(eggs_Didymellaceae_family$Timepoint_Cat=="1 month incubation"),]

ggplot(eggs_Didymellaceae_family_oneMonth, aes(x= washed_unwashed, y= Abundance, fill = washed_unwashed, colour = washed_unwashed)) +
  theme_bw() +
  labs(y= "Didymellaceae RA (%)") +
  geom_boxplot(alpha = 0.5, linewidth = 1) + geom_point(size = 4, shape = 18) +
  scale_y_continuous(expand = c(0.001,0,2.35,0)) +
  scale_fill_manual(values = washed_palette) +
  scale_colour_manual(values = washed_palette) +
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", linewidth = 1),
        panel.grid.major.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 36),
        axis.text.y = element_text(size = 16, colour = "black"),
        axis.ticks.y = element_line(colour = "black", size = 0.75))

kruskal.test(eggs_Didymellaceae_family_oneMonth$Abundance, eggs_Didymellaceae_family_oneMonth$washed_unwashed) # sig. p = 0.007

## Didymosphaeriaceae ####
eggs_Didymosphaeriaceae_family <- eggs_ITS_family_ra_filt_melt[which(eggs_ITS_family_ra_filt_melt$Family=="Didymosphaeriaceae"),]

ggplot(eggs_Didymosphaeriaceae_family, aes(x= Timepoint_Cat, y= Abundance, fill = Timepoint_Cat, colour = Timepoint_Cat)) +
  theme_bw() +
  labs(y= "Didymosphaeriaceae RA (%)") +
  geom_boxplot(alpha = 0.5, linewidth = 1) + geom_point(size = 4, shape = 18) +
  scale_x_discrete(limits = c("egg in nest","arrival at WPZ","1 month incubation")) +
  scale_y_continuous(expand = c(0.001,0,1.18,0)) +
  scale_fill_manual(values = egg_timepoint_palette) +
  scale_colour_manual(values = egg_timepoint_palette) +
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", linewidth = 1),
        panel.grid.major.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 36),
        axis.text.y = element_text(size = 16, colour = "black"),
        axis.ticks.y = element_line(colour = "black", size = 0.75))

pairwise.wilcox.test(eggs_Didymosphaeriaceae_family$Abundance, eggs_Didymosphaeriaceae_family$Timepoint_Cat) # sig. more in 1 month (p < 0.001)

# washed v unwashed
eggs_Didymosphaeriaceae_family_oneMonth <- eggs_Didymosphaeriaceae_family[which(eggs_Didymosphaeriaceae_family$Timepoint_Cat=="1 month incubation"),]

ggplot(eggs_Didymosphaeriaceae_family_oneMonth, aes(x= washed_unwashed, y= Abundance, fill = washed_unwashed, colour = washed_unwashed)) +
  theme_bw() +
  labs(y= "Didymosphaeriaceae RA (%)") +
  geom_boxplot(alpha = 0.5, linewidth = 1) + geom_point(size = 4, shape = 18) +
  scale_y_continuous(expand = c(0.001,0,4.13,0)) +
  scale_fill_manual(values = washed_palette) +
  scale_colour_manual(values = washed_palette) +
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", linewidth = 1),
        panel.grid.major.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 36),
        axis.text.y = element_text(size = 16, colour = "black"),
        axis.ticks.y = element_line(colour = "black", size = 0.75))

kruskal.test(eggs_Didymosphaeriaceae_family_oneMonth$Abundance, eggs_Didymosphaeriaceae_family_oneMonth$washed_unwashed) # sig. p = 0.007



#################################### LOOK FOR EMYDOMYCES ####
ITS_ra_species <- tax_glom(ITS_ra, taxrank = "Species", NArm = F)
ITS_ra_genus <- tax_glom(ITS_ra, taxrank = "Genus", NArm = F)

emydomyces_genus <- subset_taxa(ITS_ra_genus, Genus=="Emydomyces")
emydomyces_genus_melt <- psmelt(emydomyces_genus)

write.csv(otu_table(emydomyces_genus),"emydomyces_genus_otus.csv")
write.csv(sample_data(emydomyces_genus),"emydomyces_genus_metadata.csv")

## EGGS ####
egg_emydomyces_genus_melt <- emydomyces_genus_melt[which(emydomyces_genus_melt$matrix=="egg"),]

ggplot(egg_emydomyces_genus_melt, aes(x= Timepoint_Cat, y= Abundance, fill = Timepoint_Cat, colour = Timepoint_Cat)) +
  theme_bw() +
  labs(y= expression(~italic(Emydomyces)~" RA (%)")) +
  geom_boxplot(alpha = 0.5, linewidth = 1) + geom_point(size = 4, shape = 18) +
  scale_x_discrete(limits = c("egg in nest","arrival at WPZ","1 month incubation")) +
  scale_y_continuous(expand = c(0.028,0,0.014,0)) +
  scale_fill_manual(values = egg_timepoint_palette) +
  scale_colour_manual(values = egg_timepoint_palette) +
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", linewidth = 1),
        panel.grid.major.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 36),
        axis.text.y = element_text(size = 16, colour = "black"),
        axis.ticks.y = element_line(colour = "black", size = 0.75))

kruskal.test(egg_emydomyces_genus_melt$Abundance, egg_emydomyces_genus_melt$Timepoint_Cat) #NS

egg_emydomyces_genus_melt_onemonth <- egg_emydomyces_genus_melt[which(egg_emydomyces_genus_melt$Timepoint_Cat=="1 month incubation"),]

ggplot(egg_emydomyces_genus_melt_onemonth, aes(x= washed_unwashed, y= Abundance, fill = washed_unwashed, colour = washed_unwashed)) +
  theme_bw() +
  labs(y= expression(~italic(Emydomyces)~" RA (%)")) +
  geom_boxplot(alpha = 0.5, linewidth = 1) + geom_point(size = 4, shape = 18) +
  scale_y_continuous(expand = c(0.001,0,0.13,0)) +
  scale_fill_manual(values = washed_palette) +
  scale_colour_manual(values = washed_palette) +
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", linewidth = 1),
        panel.grid.major.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 36),
        axis.text.y = element_text(size = 16, colour = "black"),
        axis.ticks.y = element_line(colour = "black", size = 0.75))

kruskal.test(egg_emydomyces_genus_melt_onemonth$Abundance,egg_emydomyces_genus_melt_onemonth$washed_unwashed) # NS, p = 0.9

#### TURTLES ####
## CARAPACE ####
carapace_emydomyces_genus_melt <- emydomyces_genus_melt[which(emydomyces_genus_melt$matrix_specific=="turtle_carapace"),]

ggplot(carapace_emydomyces_genus_melt, aes(x= washed_unwashed, y= Abundance, fill = washed_unwashed, colour = washed_unwashed)) +
  facet_wrap(~factor(Timepoint_Cat, levels = c("hatching of egg","3 months post hatch","6 months post hatch","1 month post move outside","release to wild")), nrow = 1) + theme_bw() +
  labs(y= expression(~italic(Emydomyces)~" RA (%)")) +
  geom_boxplot(alpha = 0.5, linewidth = 1) + geom_point(size = 4, shape = 18) +
  #scale_x_discrete(limits = c("egg in nest","arrival at WPZ","1 month incubation")) +
  scale_y_continuous(expand = c(0.028,0,0.014,0)) +
  scale_fill_manual(values = washed_palette) +
  scale_colour_manual(values = washed_palette) +
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", linewidth = 1),
        panel.grid.major.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 36),
        axis.text.y = element_text(size = 16, colour = "black"),
        axis.ticks.y = element_line(colour = "black", size = 0.75))

kruskal.test(egg_emydomyces_genus_melt$Abundance, egg_emydomyces_genus_melt$Timepoint_Cat) #NS

## PLASTRON ####
plastron_emydomyces_genus_melt <- emydomyces_genus_melt[which(emydomyces_genus_melt$matrix_specific=="turtle_plastron"),]

ggplot(plastron_emydomyces_genus_melt, aes(x= washed_unwashed, y= Abundance, fill = washed_unwashed, colour = washed_unwashed)) +
  facet_wrap(~factor(Timepoint_Cat, levels = c("hatching of egg","3 months post hatch","6 months post hatch","1 month post move outside","release to wild")), nrow = 1) + theme_bw() +
  labs(y= expression(~italic(Emydomyces)~" RA (%)")) +
  geom_boxplot(alpha = 0.5, linewidth = 1) + geom_point(size = 4, shape = 18) +
  #scale_x_discrete(limits = c("egg in nest","arrival at WPZ","1 month incubation")) +
  scale_y_continuous(expand = c(0.028,0,0.014,0)) +
  scale_fill_manual(values = washed_palette) +
  scale_colour_manual(values = washed_palette) +
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", linewidth = 1),
        panel.grid.major.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 36),
        axis.text.y = element_text(size = 16, colour = "black"),
        axis.ticks.y = element_line(colour = "black", size = 0.75))

kruskal.test(egg_emydomyces_genus_melt$Abundance, egg_emydomyces_genus_melt$Timepoint_Cat) #NS

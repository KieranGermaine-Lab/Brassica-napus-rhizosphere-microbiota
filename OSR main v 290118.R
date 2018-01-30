
# Ridhdhi Rathore
# Oilseed rape Rhizosphere/Bulk soil microbiome
# 29/01/2018
#
# useful information:
# http://joey711.github.io/phyloseq-demo/phyloseq-demo.html
# http://joey711.github.io/phyloseq/Example-Data.html


#These initial commands are required to clean-up the memory and start a new session
rm(list=ls())
dev.off()

#set working directory
setwd("C:/R script_Chap.2/")

#1st time Phyloseq installation
#source("http://bioconductor.org/biocLite.R")
#biocLite("phyloseq")
#biocLite("DESeq2")
#biocLite("VennDiagram")
#biocLite("PMCMR")
#biocLite("ComplexHeatmap")

#load the required packages (all included in Phyloseq, but is necessary to invoke them)
library("phyloseq")
library("DESeq2")
library("ggplot2")
library("vegan")
library ("ape")
library("VennDiagram")
library("PMCMR")
library("ComplexHeatmap")
sessionInfo()

#############################################################

#############################################################
#import the count matrix and the desing file

#OTU table this file has been generated using USEARCH v8 and QIIME 1.9.0. In the OTU ids, OTU abundance information has been removed
dat_info <- read.delim("Ridhdhi_OTU_table.txt", skip=1, sep = "\t", header=T, row.names= 1, blank.lines.skip = FALSE)

#inspect the file using the command dim
dim(dat_info)

#extract the total number of reads clustered at OTU 97% identiy (the number underneath the Bn identifier represents the total number of reads clustered for that sample) 
OTU_97_reads <- sort(colSums(dat_info[, 1:36]))
OTU_97_reads
#total reads
OTU_97_reads_sum <- sum(colSums(dat_info[, 1:36]))
OTU_97_reads_sum

#design file
design <- read.delim("Ridhdhi_Mapping_file_r180317.txt", sep = "\t", header=TRUE, row.names=1)
design

#remove chloroplast e mitochondria OTUs from the original dataset
Chloroplast <- dat_info[grepl("Chloroplast", dat_info$taxonomy), ]
dim(Chloroplast)
Chloroplast[1:10, ]
#write.table(Chloroplast, file="chloroplast.txt", sep="\t")

mitochondria <- dat_info[grepl("mitochondria", dat_info$taxonomy), ]
dim(mitochondria)
mitochondria[1:3, ]
#write.table(mitochondria, file="mitochondria.txt", sep="\t")

#Filter plant-derived OTUs from the OTU table
noPlants <- setdiff(rownames(dat_info), c(rownames(Chloroplast), rownames(mitochondria)))
#inspect the results
length(rownames(dat_info))
length(noPlants)

#save this piece of information for qiime purposes. We need to create a a new OTU table in QIIME to generate the taxa tables
#write(noPlants, "Ridhdhi_noPlant_OTUs_id.txt")

#Generate a new OTU table depleted with chloroplast and mitochondria OTUs
dat_info_noPlants <- dat_info[noPlants, ]
#write.table(dat_info_noPlants, file="Ridhdhi_dat_info_noPlants.txt", sep="\t")

#create a new count matrix without OTUs assigned to Choloplast and Mitochondria
dat_count <- dat_info[, rownames(design)]
dat_count_noplants <- dat_info[noPlants, rownames(design)]
dim(dat_count_noplants)
#write.table(dat_count_noplants, file="Ridhdhi_dat_count_noPlants.txt", sep="\t")

##########################################################################################################
#generate a phyloseq object
#########################################

#a) The OTU table counts
Ridhdhi_OTU <- otu_table(dat_count_noplants, taxa_are_rows=TRUE)

#b) The taxonomy information
#Note this is a new file generated in excell from the output of the command of lines 93-97
#it is a tab-delimited file with 8 columns, the header names are: OTU id; "Kingdom", "Phylum",  "Class",   "Order",   "Family",  "Genus",   "Species",
#and the empty cells are filled with the term 'Unassigned'
Ridhdhi_taxa_ordered <- read.delim ("Ridhdhi_dat_tax_noPlants_ordered.txt", sep = "\t", row.names=1, header=T, blank.lines.skip = FALSE)
Ridhdhi_taxa <- tax_table(as.matrix(Ridhdhi_taxa_ordered))
dim(Ridhdhi_taxa)

#c) The mapping file 
Ridhdhi_map <- sample_data(design)

#d) The phylogenetic tree 
Ridhdhi_tree <- read_tree("Ridhdhi_OTU_table_noPlants_pynast_tree.tre")
#check whether the tree is rooted
is.rooted(Ridhdhi_tree)

#Root the tree
#Identify unique taxa
unique_OTUs <- unique(Ridhdhi_taxa[,2])
unique_OTUs

#OK we could use a Planctomycetes as an outgroup 
tax <- as.data.frame(tax_table(Ridhdhi_taxa))
outgroup <- tax[grepl("p__Planctomycetes", tax$Phylum), ]
dim(outgroup)
outgroup 
#OK now we can identify the top abundant OTU of this group
sort(rowSums(dat_count_noplants[rownames(outgroup), ]))

#With 701 reads, OTUs281 is relatively abundant, let's pick this OTUs as an outgroup
newRoot = c("OTU281")
Ridhdhi_tree <- root(Ridhdhi_tree,newRoot,resolve.root = TRUE)

#let's check it now
is.rooted(Ridhdhi_tree)

#merge the files and create the phyloseq object
Ridhdhi_data_phyloseq <- merge_phyloseq(Ridhdhi_OTU, Ridhdhi_taxa, Ridhdhi_map,  Ridhdhi_tree)

#inspect the generated data
Ridhdhi_data_phyloseq
sum(colSums(otu_table(Ridhdhi_data_phyloseq)))
dim(dat_count_noplants)
sum(colSums(dat_count_noplants))

######################################################################################################
#Top 10 Phyla relative abundance (Figure 3)
################################################

#This table was generated in QIIME w/o plant-derived OTUs and biological replicates are already averaged according to the levels of the factor "Description" in the mapping file
dat_info_taxa_Phylum <- read.delim("Description_otu_table_L2.txt", skip=1, sep = "\t", header=T, row.names= 1, blank.lines.skip = FALSE)
colnames(dat_info_taxa_Phylum)

# transform the data in % for visualisation
dat_count_taxa_Phylum <- dat_info_taxa_Phylum
dat_norm_Phylum <- ((dat_count_taxa_Phylum)/colSums(dat_count_taxa_Phylum,na=T)) * 100 
dat_norm_Phylum[1:5, ]

# determine the average % of reads for each phylum
Phylum_mean_sorted <- dat_norm_Phylum[(order(-rowSums(dat_norm_Phylum))), ] 

#Calculate the contribution of the top 10 Phyla to the total dataset
Phylum_mean_topRank <- Phylum_mean_sorted[1:10, ]
Phylum_mean_topRank 
colSums(Phylum_mean_topRank)

#write.table(Phylum_mean_topRank, file="Phylum_mean_topRank.txt", sep="\t")

#overall
mean(colSums(Phylum_mean_topRank))

Phylum_mean_topRank_10 <- as.matrix(Phylum_mean_topRank[1:10, ]) 
#first we need to arrange the samples in a coherent way (i.e. according to the experiment)
colnames(Phylum_mean_topRank_10)

Phylum_mean_topRank_10_samples <- c( "VegBkConventional", "VegRzConventional", "VegBkConservation", "VegRzConservation", "FlBkConventional", "FlRzConventional", "FlBkConservation", "FlRzConservation", "HvBkConventional", "HvRzConventional", "HvBkConservation", "HvRzConservation")

#now we will use the order of samples we have generated above to sort the count matrix
Phylum_mean_topRank_10_ordered <- Phylum_mean_topRank_10[ ,Phylum_mean_topRank_10_samples] 
colnames(Phylum_mean_topRank_10_ordered)

#Inspect the generated files
Phylum_mean_topRank_10_ordered[1:5, ]

#Now we can plot them
# Stacked Bar Plot with Colors and Legend
barplot(Phylum_mean_topRank_10_ordered, main="Phylum Distribution",
        xlab="Samples", ylab = "Reads %", ylim = c(0,100), col=c("darkmagenta", "gold", "green3", "darkslateblue", "pink", "red", "lightblue1", "olivedrab", "deeppink", "gray"), beside=FALSE,   legend = rownames(Phylum_mean_topRank_10_ordered))

#Due to size limits the legend covers part of the graph, save the graph as .eps file 
#barplot_no legend: you can save this as image and then reconstruct the leged using the previous figure. Note: the order of samples can be inferred from the command on line 185

barplot(Phylum_mean_topRank_10_ordered, main="Phylum Distribution",
        xlab="Samples", ylab = "Reads %", ylim = c(0,100), col=c("darkmagenta", "gold", "green3", "darkslateblue", "pink", "red", "lightblue1", "olivedrab", "deeppink", "gray"), beside=FALSE)

##############################################################################################################################################################################################################################################
#Bacterial Family Distribution (Figure 5)
############################################

#This table was generated in QIIME w/o plant-derived OTUs and biological replicates are already averaged according to the levels of the factor "Description" in the mapping file
dat_info_taxa_Phylum <- read.delim("Family_discription.txt", skip=1, sep = "\t", header=T, row.names= 1, blank.lines.skip = FALSE)
colnames(dat_info_taxa_Phylum)

Family_toprank <- as.matrix(dat_info_taxa_Phylum[1: 18, ])
colnames(Family_toprank)

Family_toprank_samples <- c("BK_VG_CT", "BK_VG_ST", "BK_FL_CT", "BK_FL_ST", "BK_HV_CT", "BK_HV_ST", "RZ_VG_CT", "RZ_VG_ST", "RZ_FL_CT", "RZ_FL_ST", "RZ_HV_CT", "RZ_HV_ST")

Family_toprank_ordered <- Family_toprank[ ,Family_toprank_samples]
colnames(Family_toprank_ordered)

Family_toprank_ordered[1:5, ]

barplot(Family_toprank_ordered, main="Phylum Distribution",
        xlab="Samples", ylab = "Relative abundance %", ylim = c(0,50), col=c("darkblue", "yellow", "green3", "red1", "cyan", "blue", "greenyellow", "dimgray", "gray", "deeppink1", "forestgreen", "lightgoldenrod1", "darkred", "lightsalmon", "black", "mediumorchid1", "dodgerblue2", "gray98"), beside=FALSE, legend = rownames(Family_toprank_ordered))

barplot(Family_toprank_ordered, main="Phylum Distribution",
        xlab="Samples", ylab = "Relative abundance %", ylim = c(0,45), col=c("darkblue", "yellow", "green3", "red1", "cyan", "blue", "greenyellow", "dimgray", "gray", "deeppink1", "forestgreen", "lightgoldenrod1", "darkred", "lightsalmon", "black", "mediumorchid1", "dodgerblue2", "gray98"), beside=FALSE)

######################################################################################################
#Alpha diversity (Figure 1)
#########################################################

#import the rarefied OTU counts (note file name counts2.txt this file has been generated in excel and includes the #OTU ID as header of the first column)
dat_count_rare <- read.delim("Ridhdhi_data_phyloseq_rare_table2.txt", sep = "\t", header=T, row.names= 1, blank.lines.skip = FALSE)

#inspect the generated file
dim(dat_count_rare)

#generate a new phyloseq object which will contain only the rarefied counts and the design file (only these two pieces of information are required for alphadiversity calculation)
Ridhdhi_OTU_rare <- otu_table(dat_count_rare, taxa_are_rows=TRUE)
Ridhdhi_data_rare_phyloseq <- merge_phyloseq(Ridhdhi_OTU_rare , Ridhdhi_map)

#Inspect the generated file
Ridhdhi_data_rare_phyloseq

#number of reads per sample: original dataset different reads per sample
sample_sums(Ridhdhi_data_phyloseq)

#number of reads per sample: rarefied dataset, they are all the same
sample_sums(Ridhdhi_data_rare_phyloseq)

#Index calculation
Ridhdhi_alpha_rare <-  estimate_richness(Ridhdhi_data_rare_phyloseq, measures = c("Observed", "Chao1", "Shannon"))

#generate a new dataframe for data visualisation and statistical analysis
design2 <- design[rownames(Ridhdhi_alpha_rare), ]

Ridhdhi_alpha_rare_info <- cbind(design2, Ridhdhi_alpha_rare)

#check the new dataset: it contains both description of the samples and alpha diversity indices 
Ridhdhi_alpha_rare_info
#write.table(Ridhdhi_alpha_rare_info, file="Ridhdhi_alpha_rare_info.txt", sep="\t")

#generate a box plot for data visualisation
#re-order the factors
Ridhdhi_alpha_rare_info$ Description <- ordered(Ridhdhi_alpha_rare_info$Description, levels=c("VegBkConventional", "VegRzConventional", "VegBkConservation", "VegRzConservation", "FlBkConventional", "FlRzConventional", "FlBkConservation", "FlRzConservation", "HvBkConventional", "HvRzConventional", "HvBkConservation", "HvRzConservation"))

with(Ridhdhi_alpha_rare_info, boxplot(Observed  ~ Description, xlab = "Compartments", ylab = "Observed OTUs", main = "Observed OTUs",
                                      col=(c("olivedrab4", "olivedrab1", "olivedrab4", "olivedrab1", "orange", "yellow", "orange", "yellow", "tan4", "khaki1", "tan4", "khaki1"))))

with(Ridhdhi_alpha_rare_info, boxplot(Chao1  ~ Description, xlab = "Compartments", ylab = "Chao1", main = "Chao1",
                                      col=(c("olivedrab4", "olivedrab1", "olivedrab4", "olivedrab1", "orange", "yellow", "orange", "yellow", "tan4", "khaki1", "tan4", "khaki1"))))

with(Ridhdhi_alpha_rare_info, boxplot(Shannon  ~ Description, xlab = "Compartments", ylab = "Shannon", main = "Shannon",
                                      col=(c("olivedrab4", "olivedrab1", "olivedrab4", "olivedrab1", "orange", "yellow", "orange", "yellow", "tan4", "khaki1", "tan4", "khaki1"))))


#inspect the normal distribution of the datasets
shapiro.test(Ridhdhi_alpha_rare_info$Observed)
shapiro.test(Ridhdhi_alpha_rare_info$Chao1)
shapiro.test(Ridhdhi_alpha_rare_info$Shannon)

#Observed
Observed_stats <- aov(Ridhdhi_alpha_rare_info$Observed ~ TimePoint * Tillage * Compartments, data = Ridhdhi_alpha_rare_info)
summary(Observed_stats)
Tukey_HSD_Observed <- TukeyHSD(Observed_stats)
summary(Tukey_HSD_Observed)

#Chao1
Chao1_stats <- aov(Ridhdhi_alpha_rare_info$Chao1 ~ TimePoint * Tillage * Compartments, data = Ridhdhi_alpha_rare_info)
summary(Chao1_stats)
Tukey_HSD_Chao1 <- TukeyHSD(Chao1_stats)
summary(Tukey_HSD_Chao1)

#Shannon
Shannon_stats<- aov(Ridhdhi_alpha_rare_info$Shannon ~ TimePoint * Tillage * Compartments, data = Ridhdhi_alpha_rare_info)
summary(Shannon_stats)
Tukey_HSD_Shannon <- TukeyHSD(Shannon_stats)
summary(Tukey_HSD_Shannon)

###################################################################################################################################################################################################
# Betadiversity calculation (Figure 2, Supplementary figure S1, S2, and S3)
#############################################################################

#Transform the count in relative abundance
Ridhdhi_data_phyloseq_prop <- transform_sample_counts(Ridhdhi_data_phyloseq,  function(x) 1e+06 * x/sum(x))

##############################
#Rhizosphere
##############

Ridhdhi_data_phyloseq_prop_Rhizosphere <- subset_samples(Ridhdhi_data_phyloseq_prop, Compartments=="Rhizosphere")

# CAP TimePoint effect
# Bray distance
Ridhdhi_data_phyloseq_prop_Rhizosphere.CAP_bray1 <- ordinate(Ridhdhi_data_phyloseq_prop_Rhizosphere, "CAP", "bray", ~TimePoint)
plot_ordination(Ridhdhi_data_phyloseq_prop_Rhizosphere, Ridhdhi_data_phyloseq_prop_Rhizosphere.CAP_bray1, color = "TimePoint")
p2=plot_ordination(Ridhdhi_data_phyloseq_prop_Rhizosphere, Ridhdhi_data_phyloseq_prop_Rhizosphere.CAP_bray1 , type="samples", color = "TimePoint", shape = "Tillage")
p2 = p2 + geom_point(size = 5, alpha = 0.75)
p2 = p2 + scale_colour_manual(values = c("gold2", "tan4", "green4")) 
p2 + ggtitle("CAP~TimePoint, Rhizosphere, bray distance")

# Weighted UniFrac distance
Ridhdhi_data_phyloseq_prop_Rhizosphere.CAP_wunifrac1 <- ordinate(Ridhdhi_data_phyloseq_prop_Rhizosphere, "CAP", "unifrac", weighted =TRUE, ~TimePoint)
plot_ordination(Ridhdhi_data_phyloseq_prop_Rhizosphere, Ridhdhi_data_phyloseq_prop_Rhizosphere.CAP_wunifrac1, color = "TimePoint")
p2=plot_ordination(Ridhdhi_data_phyloseq_prop_Rhizosphere, Ridhdhi_data_phyloseq_prop_Rhizosphere.CAP_wunifrac1 , type="samples", color = "TimePoint", shape = "Tillage")
p2 = p2 + scale_colour_manual(values = c("gold2", "tan4", "green4")) 
p2 = p2 + geom_point(size = 5, alpha = 0.75)
p2 + ggtitle("CAP~TimePoint, Rhizosphere, wunifrac distance")

# CAP Tillage effect
# Bray distance
Ridhdhi_data_phyloseq_prop_Rhizosphere.CAP_bray2 <- ordinate(Ridhdhi_data_phyloseq_prop_Rhizosphere, "CAP", "bray", ~Tillage)
plot_ordination(Ridhdhi_data_phyloseq_prop_Rhizosphere, Ridhdhi_data_phyloseq_prop_Rhizosphere.CAP_bray2, color = "Tillage")
p2=plot_ordination(Ridhdhi_data_phyloseq_prop_Rhizosphere, Ridhdhi_data_phyloseq_prop_Rhizosphere.CAP_bray2 , type="samples", color = "Tillage", shape = "TimePoint")
p2 = p2 + geom_point(size = 5, alpha = 0.75)
p2 + ggtitle("CAP~Tillage, Rhizosphere, bray distance")

# Wighted UniFrac
# Weighted UniFrac distance
Ridhdhi_data_phyloseq_prop_Rhizosphere.CAP_wunifrac2 <- ordinate(Ridhdhi_data_phyloseq_prop_Rhizosphere, "CAP", "unifrac", weighted =TRUE, ~Tillage)
plot_ordination(Ridhdhi_data_phyloseq_prop_Rhizosphere, Ridhdhi_data_phyloseq_prop_Rhizosphere.CAP_wunifrac2, color = "Tillage")
p2=plot_ordination(Ridhdhi_data_phyloseq_prop_Rhizosphere, Ridhdhi_data_phyloseq_prop_Rhizosphere.CAP_wunifrac2 , type="samples", color = "Tillage", shape = "TimePoint")
p2 = p2 + geom_point(size = 5, alpha = 0.75)
p2 + ggtitle("CAP~Tillage, Rhizosphere, wunifrac distance")

anova(Ridhdhi_data_phyloseq_prop_Rhizosphere.CAP_bray1, permutations = 5000)
anova(Ridhdhi_data_phyloseq_prop_Rhizosphere.CAP_wunifrac1, permutations = 5000)
anova(Ridhdhi_data_phyloseq_prop_Rhizosphere.CAP_bray2, permutations = 5000)
anova(Ridhdhi_data_phyloseq_prop_Rhizosphere.CAP_wunifrac2, permutations = 5000)

###################################################################
#Bulk soil
#########################

Ridhdhi_data_phyloseq_prop_Bulk <- subset_samples(Ridhdhi_data_phyloseq_prop, Compartments=="Bulk soil")

# CAP TimePoint effect
# Bray distance
Ridhdhi_data_phyloseq_prop_Bulk.CAP_bray1 <- ordinate(Ridhdhi_data_phyloseq_prop_Bulk, "CAP", "bray", ~TimePoint)
plot_ordination(Ridhdhi_data_phyloseq_prop_Bulk, Ridhdhi_data_phyloseq_prop_Bulk.CAP_bray1, color = "TimePoint")
p2=plot_ordination(Ridhdhi_data_phyloseq_prop_Bulk, Ridhdhi_data_phyloseq_prop_Bulk.CAP_bray1 , type="samples", color = "TimePoint", shape = "Tillage")
p2 = p2 + geom_point(size = 5, alpha = 0.75)
p2 = p2 + scale_colour_manual(values = c("gold2", "tan3", "green4")) 
p2 + ggtitle("CAP~TimePoint, Bulk soil, bray distance")

# Wighted UniFrac
Ridhdhi_data_phyloseq_prop_Bulk.CAP_wunifrac1 <- ordinate(Ridhdhi_data_phyloseq_prop_Bulk, "CAP", "unifrac", weighted =TRUE, ~TimePoint)
plot_ordination(Ridhdhi_data_phyloseq_prop_Bulk, Ridhdhi_data_phyloseq_prop_Bulk.CAP_wunifrac1, color = "TimePoint")
p2=plot_ordination(Ridhdhi_data_phyloseq_prop_Bulk, Ridhdhi_data_phyloseq_prop_Bulk.CAP_wunifrac1 , type="samples", color = "TimePoint", shape = "Tillage")
p2 = p2 + geom_point(size = 5, alpha = 0.75)
p2 = p2 + scale_colour_manual(values = c("gold2", "tan4", "green4")) 
p2 + ggtitle("CAP~TimePoint, Bulk soil, wunifrac distance")

# CAP Tillage effect
# Bray distance
Ridhdhi_data_phyloseq_prop_Bulk.CAP_bray2 <- ordinate(Ridhdhi_data_phyloseq_prop_Bulk, "CAP", "bray", ~Tillage)
plot_ordination(Ridhdhi_data_phyloseq_prop_Bulk, Ridhdhi_data_phyloseq_prop_Bulk.CAP_bray2, color = "Tillage")
p2=plot_ordination(Ridhdhi_data_phyloseq_prop_Bulk, Ridhdhi_data_phyloseq_prop_Bulk.CAP_bray2 , type="samples", color = "Tillage", shape = "TimePoint")
p2 = p2 + geom_point(size = 5, alpha = 0.75)
p2 + ggtitle("CAP~Tillage, Bulk soil, bray distance")

# Wighted UniFrac
Ridhdhi_data_phyloseq_prop_Bulk.CAP_wunifrac2 <- ordinate(Ridhdhi_data_phyloseq_prop_Bulk, "CAP", "unifrac", weighted =TRUE, ~Tillage)
plot_ordination(Ridhdhi_data_phyloseq_prop_Bulk, Ridhdhi_data_phyloseq_prop_Bulk.CAP_wunifrac2, color = "Tillage")
p2=plot_ordination(Ridhdhi_data_phyloseq_prop_Bulk, Ridhdhi_data_phyloseq_prop_Bulk.CAP_wunifrac2 , type="samples", color = "Tillage", shape = "TimePoint")
p2 = p2 + geom_point(size = 5, alpha = 0.75)
p2 + ggtitle("CAP~Tillage, Bulk soil, wunifrac distance")

anova(Ridhdhi_data_phyloseq_prop_Bulk.CAP_bray1, permutations = 5000)
anova(Ridhdhi_data_phyloseq_prop_Bulk.CAP_wunifrac1, permutations = 5000)
anova(Ridhdhi_data_phyloseq_prop_Bulk.CAP_bray1, permutations = 5000)
anova(Ridhdhi_data_phyloseq_prop_Bulk.CAP_wunifrac1, permutations = 5000)

###################################################################
# Conservation tillage (ST)
#########################

Ridhdhi_data_phyloseq_prop_Conservation <- subset_samples(Ridhdhi_data_phyloseq_prop, Tillage=="Conservation")

# CAP TimePoint effect
# Bray distance
Ridhdhi_data_phyloseq_prop_Conservation.CAP_bray1 <- ordinate(Ridhdhi_data_phyloseq_prop_Conservation, "CAP", "bray", ~TimePoint)
plot_ordination(Ridhdhi_data_phyloseq_prop_Conservation, Ridhdhi_data_phyloseq_prop_Conservation.CAP_bray1, color = "TimePoint")
p2=plot_ordination(Ridhdhi_data_phyloseq_prop_Conservation, Ridhdhi_data_phyloseq_prop_Conservation.CAP_bray1, type="samples", color = "TimePoint", shape = "Compartments")
p2 = p2 + geom_point(size = 5, alpha = 0.75)
p2 = p2 + scale_colour_manual(values = c("gold2", "tan4", "green4")) 
p2 + ggtitle("CAP~TimePoint, Conservation tillage, bray distance")

# Wighted UniFrac
Ridhdhi_data_phyloseq_prop_Conservation.CAP_wunifrac1 <- ordinate(Ridhdhi_data_phyloseq_prop_Conservation, "CAP", "unifrac", weighted =TRUE, ~TimePoint)
plot_ordination(Ridhdhi_data_phyloseq_prop_Conservation, Ridhdhi_data_phyloseq_prop_Conservation.CAP_wunifrac1, color = "TimePoint")
p2=plot_ordination(Ridhdhi_data_phyloseq_prop_Conservation, Ridhdhi_data_phyloseq_prop_Conservation.CAP_wunifrac1, type="samples", color = "TimePoint", shape = "Compartments")
p2 = p2 + geom_point(size = 5, alpha = 0.75)
p2 = p2 + scale_colour_manual(values = c("gold2", "tan4", "green4")) 
p2 + ggtitle("CAP~TimePoint, Conservation tillage, wunifrac distance")

# CAP Compartment effect
# Bray distance
Ridhdhi_data_phyloseq_prop_Conservation.CAP_bray2 <- ordinate(Ridhdhi_data_phyloseq_prop_Conservation, "CAP", "bray", ~Compartments)
plot_ordination(Ridhdhi_data_phyloseq_prop_Conservation, Ridhdhi_data_phyloseq_prop_Conservation.CAP_bray2, color = "Compartments")
p2=plot_ordination(Ridhdhi_data_phyloseq_prop_Conservation, Ridhdhi_data_phyloseq_prop_Conservation.CAP_bray2 , type="samples", color = "Compartments", shape = "TimePoint")
p2 = p2 + geom_point(size = 5, alpha = 0.75)
p2 = p2 + scale_colour_manual(values = c("tan4", "tan1")) 
p2 + ggtitle("CAP~Compartments, Conservation tillage, bray distance")

# Wighted UniFrac
Ridhdhi_data_phyloseq_prop_Conservation.CAP_wunifrac2 <- ordinate(Ridhdhi_data_phyloseq_prop_Conservation, "CAP", "unifrac", weighted =TRUE, ~Compartments)
plot_ordination(Ridhdhi_data_phyloseq_prop_Conservation, Ridhdhi_data_phyloseq_prop_Conservation.CAP_wunifrac2, color = "Compartments")
p2=plot_ordination(Ridhdhi_data_phyloseq_prop_Conservation, Ridhdhi_data_phyloseq_prop_Conservation.CAP_wunifrac2, type="samples", color = "Compartments", shape = "TimePoint")
p2 = p2 + geom_point(size = 5, alpha = 0.75)
p2 = p2 + scale_colour_manual(values = c("tan4", "tan1")) 
p2 + ggtitle("CAP~Compartments, Conservation tillage, wunifrac distance")

anova(Ridhdhi_data_phyloseq_prop_Conservation.CAP_bray1, permutations = 5000)
anova(Ridhdhi_data_phyloseq_prop_Conservation.CAP_wunifrac1, permutations = 5000)
anova(Ridhdhi_data_phyloseq_prop_Conservation.CAP_bray2, permutations = 5000)
anova(Ridhdhi_data_phyloseq_prop_Conservation.CAP_wunifrac2, permutations = 5000)

###################################################################
# Conventional tillage (CT)
#########################
Ridhdhi_data_phyloseq_prop_Conventional <- subset_samples(Ridhdhi_data_phyloseq_prop, Tillage=="Conventional")

# CAP TimePoint effect
# Bray distance
Ridhdhi_data_phyloseq_prop_Conventional.CAP_bray1 <- ordinate(Ridhdhi_data_phyloseq_prop_Conventional, "CAP", "bray", ~TimePoint)
plot_ordination(Ridhdhi_data_phyloseq_prop_Conventional, Ridhdhi_data_phyloseq_prop_Conventional.CAP_bray1, color = "TimePoint")
p2=plot_ordination(Ridhdhi_data_phyloseq_prop_Conventional, Ridhdhi_data_phyloseq_prop_Conventional.CAP_bray1, type="samples", color = "TimePoint", shape = "Compartments")
p2 = p2 + geom_point(size = 5, alpha = 0.75)
p2 = p2 + scale_colour_manual(values = c("gold2", "tan4", "green4")) 
p2 + ggtitle("CAP~TimePoint, Conventional tillage, bray distance")

# Wighted UniFrac
Ridhdhi_data_phyloseq_prop_Conventional.CAP_wunifrac1 <- ordinate(Ridhdhi_data_phyloseq_prop_Conventional, "CAP", "wunifrac", weighted =TRUE, ~TimePoint)
plot_ordination(Ridhdhi_data_phyloseq_prop_Conventional, Ridhdhi_data_phyloseq_prop_Conventional.CAP_wunifrac1, color = "TimePoint")
p2=plot_ordination(Ridhdhi_data_phyloseq_prop_Conventional, Ridhdhi_data_phyloseq_prop_Conventional.CAP_wunifrac1 , type="samples", color = "TimePoint", shape = "Compartments")
p2 = p2 + geom_point(size = 5, alpha = 0.75)
p2 = p2 + scale_colour_manual(values = c("gold2", "tan4", "green4")) 
p2 + ggtitle("CAP~TimePoint, Conventional tillage, wunifrac distance")

# CAP Compartment effect
# Bray distance
Ridhdhi_data_phyloseq_prop_Conventional.CAP_bray2 <- ordinate(Ridhdhi_data_phyloseq_prop_Conventional, "CAP", "bray", ~Compartments)
plot_ordination(Ridhdhi_data_phyloseq_prop_Conventional, Ridhdhi_data_phyloseq_prop_Conventional.CAP_bray2, color = "Compartments")
p2=plot_ordination(Ridhdhi_data_phyloseq_prop_Conventional, Ridhdhi_data_phyloseq_prop_Conventional.CAP_bray2, type="samples", color = "Compartments", shape = "TimePoint")
p2 = p2 + geom_point(size = 5, alpha = 0.75)
p2 = p2 + scale_colour_manual(values = c("tan4", "tan1")) 
p2 + ggtitle("CAP~Compartments, Conventional tillage, bray distance")

# Wighted UniFrac
Ridhdhi_data_phyloseq_prop_Conventional.CAP_wunifrac2 <- ordinate(Ridhdhi_data_phyloseq_prop_Conventional, "CAP", "wunifrac", weighted =TRUE, ~TimePoint)
plot_ordination(Ridhdhi_data_phyloseq_prop_Conventional, Ridhdhi_data_phyloseq_prop_Conventional.CAP_wunifrac2, color = "Compartments")
p2=plot_ordination(Ridhdhi_data_phyloseq_prop_Conventional, Ridhdhi_data_phyloseq_prop_Conventional.CAP_wunifrac2, type="samples", color = "Compartments", shape = "TimePoint")
p2 = p2 + geom_point(size = 5, alpha = 0.75)
p2 = p2 + scale_colour_manual(values = c("tan4", "tan1")) 
p2 + ggtitle("CAP~Compartments, Conventional tillage, wunifrac distance")

anova(Ridhdhi_data_phyloseq_prop_Conventional.CAP_bray1, permutations = 5000)
anova(Ridhdhi_data_phyloseq_prop_Conventional.CAP_wunifrac1, permutations = 5000)
anova(Ridhdhi_data_phyloseq_prop_Conventional.CAP_bray2, permutations = 5000)
anova(Ridhdhi_data_phyloseq_prop_Conventional.CAP_wunifrac2, permutations = 5000)

#################################################################################################################################################################################
#permutational analysis of variance on dissimilarity matrices (Table - 2)
##########################################################################

#extract the dissimilarity matrices

BC1 <- phyloseq::distance(Ridhdhi_data_phyloseq_prop_Conservation, "bray")
WU1 <- phyloseq::distance(Ridhdhi_data_phyloseq_prop_Conservation, "unifrac", weighted= TRUE)

BC2 <- phyloseq::distance(Ridhdhi_data_phyloseq_prop_Conventional, "bray")
WU2 <- phyloseq::distance(Ridhdhi_data_phyloseq_prop_Conventional, "unifrac", weighted= TRUE)

BC3 <- phyloseq::distance(Ridhdhi_data_phyloseq_prop_Bulk, "bray")
WU3 <- phyloseq::distance(Ridhdhi_data_phyloseq_prop_Bulk, "unifrac", weighted= TRUE)

BC4 <- phyloseq::distance(Ridhdhi_data_phyloseq_prop_Rhizosphere, "bray")
WU4 <- phyloseq::distance(Ridhdhi_data_phyloseq_prop_Rhizosphere, "unifrac", weighted= TRUE)

Map_ST = read.table("Map_ST.txt", header = T, sep = "\t")
Map_CT = read.table("Map_CT.txt", header = T, sep = "\t")
Map_bulk = read.table("Map_bulk.txt", header = T, sep = "\t")
Map_rhizo = read.table("Map_rhizo.txt", header = T, sep = "\t")

#BC distance

adonis(BC1 ~ Map_ST$TimePoint*Map_ST$Compartments, data = Map_ST, permutations = 5000)
adonis(WU1 ~ Map_ST$TimePoint*Map_ST$Compartments, data = Map_ST, permutations = 5000)

adonis(BC2 ~ Map_CT$TimePoint*Map_CT$Compartments, data = Map_CT, permutations = 5000)
adonis(WU2 ~ Map_CT$TimePoint*Map_CT$Compartments, data = Map_CT, permutations = 5000)

adonis(BC3 ~ Map_bulk$TimePoint*Map_bulk$Tillage, data = Map_bulk, permutations = 5000)
adonis(WU3 ~ Map_bulk$TimePoint*Map_bulk$Tillage, data = Map_bulk, permutations = 5000)

adonis(BC2 ~ Map_rhizo$TimePoint*Map_rhizo$Tillage, data = Map_rhizo, permutations = 5000)
adonis(WU2 ~ Map_rhizo$TimePoint*Map_rhizo$Tillage, data = Map_rhizo, permutations = 5000)


#################################################################################
#create a deseq object
################################

#extract count data and 
Ridhdhi_OTU_counts_integer <- otu_table(Ridhdhi_data_phyloseq)
countData = as.data.frame(Ridhdhi_OTU_counts_integer)

#the design file containing sample information
colData = design

#construct a DESeq dataset combining count data and sample information
Ridhdhi_cds <- DESeqDataSetFromMatrix(countData =countData, colData=colData, design= ~ Description)

#execute the differential count analysis with the function DESeq 
Ridhdhi_cds_test <- DESeq(Ridhdhi_cds, fitType="local", betaPrior=FALSE)
plotDispEsts(Ridhdhi_cds_test)

##################################################################################
# Define the OTUs significantly enriched in TimePoints corrected for tillage 
###################################################################################

# TimePoint effect
# Rhizosphere 
# Conservation tillage ST

RZ_V_F_ST <- results(Ridhdhi_cds_test, contrast = c("Description", "VegRzConservation", "FlRzConservation"))
RZ_F_H_ST <- results(Ridhdhi_cds_test, contrast = c("Description", "FlRzConservation", "HvRzConservation"))
RZ_H_V_ST <- results(Ridhdhi_cds_test, contrast = c("Description", "HvRzConservation", "VegRzConservation"))

results(Ridhdhi_cds_test, contrast = c("Description", "VegRzConservation", "FlRzConservation"))

#1
##########
RZ_V_F_ST_FDR_005 <- RZ_V_F_ST[(rownames(RZ_V_F_ST)[which(RZ_V_F_ST$padj < 0.05)]), ]
RZ_V1_ST_enriched <- RZ_V_F_ST[(rownames(RZ_V_F_ST)[which(RZ_V_F_ST$log2FoldChange < 0)]), ]
RZ_V1_ST_enriched_FDR005 <- intersect(rownames(RZ_V_F_ST_FDR_005), rownames(RZ_V1_ST_enriched))
length(RZ_V1_ST_enriched_FDR005)
RZ_V1_ST_enriched_taxa <- Ridhdhi_taxa_ordered[RZ_V1_ST_enriched_FDR005, ]
RZ_V1_ST_enriched_FDR005_taxa <- cbind(as.data.frame(RZ_V_F_ST[RZ_V1_ST_enriched_FDR005, ]), RZ_V1_ST_enriched_taxa)
RZ_F1_ST_enriched <- RZ_V_F_ST[(rownames(RZ_V_F_ST)[which(RZ_V_F_ST$log2FoldChange > 0)]), ]
RZ_F1_ST_enriched_FDR005 <- intersect(rownames(RZ_V_F_ST_FDR_005), rownames(RZ_F1_ST_enriched))
length(RZ_F1_ST_enriched_FDR005)
RZ_F1_ST_enriched_taxa <- Ridhdhi_taxa_ordered[RZ_F1_ST_enriched_FDR005, ]
RZ_F1_ST_enriched_FDR005_taxa <- cbind(as.data.frame(RZ_V_F_ST[RZ_F1_ST_enriched_FDR005, ]), RZ_F1_ST_enriched_taxa)

#2
#########
RZ_F_H_ST_FDR_005 <- RZ_F_H_ST[(rownames(RZ_F_H_ST)[which(RZ_F_H_ST$padj < 0.05)]), ]
RZ_F2_ST_enriched <- RZ_F_H_ST[(rownames(RZ_F_H_ST)[which(RZ_F_H_ST$log2FoldChange < 0)]), ]
RZ_F2_ST_enriched_FDR005 <- intersect(rownames(RZ_F_H_ST_FDR_005), rownames(RZ_F2_ST_enriched))
length(RZ_F2_ST_enriched_FDR005)
RZ_F2_ST_enriched_taxa <- Ridhdhi_taxa_ordered[RZ_F2_ST_enriched_FDR005, ]
RZ_F2_ST_enriched_FDR005_taxa <- cbind(as.data.frame(RZ_F_H_ST[RZ_F2_ST_enriched_FDR005, ]), RZ_F2_ST_enriched_taxa)
#write.table(RZ_F2_ST_enriched_FDR005_taxa, file="RZ_F2_ST_enriched_FDR005_taxa.txt", sep="\t")
RZ_H1_ST_enriched <- RZ_F_H_ST[(rownames(RZ_F_H_ST)[which(RZ_F_H_ST$log2FoldChange > 0)]), ]
RZ_H1_ST_enriched_FDR005 <- intersect(rownames(RZ_F_H_ST_HDR_005), rownames(RZ_H1_ST_enriched))
length(RZ_H1_ST_enriched_FDR005)
RZ_H1_ST_enriched_taxa <- Ridhdhi_taxa_ordered[RZ_H1_ST_enriched_FDR005, ]
RZ_H1_ST_enriched_FDR005_taxa <- cbind(as.data.frame(RZ_F_H_ST[RZ_H1_ST_enriched_FDR005, ]), RZ_H1_ST_enriched_taxa)
#write.table(RZ_H1_ST_enriched_FDR005_taxa, file="RZ_H1_ST_enriched_FDR005_taxa.txt", sep="\t")

#3
#########
RZ_H_V_ST_FDR_005 <- RZ_H_V_ST[(rownames(RZ_H_V_ST)[which(RZ_H_V_ST$padj < 0.05)]), ]
RZ_H2_ST_enriched <- RZ_H_V_ST[(rownames(RZ_H_V_ST)[which(RZ_H_V_ST$log2FoldChange < 0)]), ]
RZ_H2_ST_enriched_FDR005 <- intersect(rownames(RZ_H_V_ST_FDR_005), rownames(RZ_H2_ST_enriched))
length(RZ_H2_ST_enriched_FDR005)
RZ_H2_ST_enriched_taxa <- Ridhdhi_taxa_ordered[RZ_H2_ST_enriched_FDR005, ]
RZ_H2_ST_enriched_FDR005_taxa <- cbind(as.data.frame(RZ_H_V_ST[RZ_H2_ST_enriched_FDR005, ]), RZ_H2_ST_enriched_taxa)
#write.table(RZ_H2_ST_enriched_FDR005_taxa, file="RZ_H2_ST_enriched_FDR005_taxa.txt", sep="\t")
RZ_V2_ST_enriched <- RZ_H_V_ST[(rownames(RZ_H_V_ST)[which(RZ_H_V_ST$log2FoldChange > 0)]), ]
RZ_V2_ST_enriched_FDR005 <- intersect(rownames(RZ_H_V_ST_FDR_005), rownames(RZ_V2_ST_enriched))
length(RZ_V2_ST_enriched_FDR005) 
RZ_V2_ST_enriched_taxa <- Ridhdhi_taxa_ordered[RZ_V2_ST_enriched_FDR005, ]
RZ_V2_ST_enriched_FDR005_taxa <- cbind(as.data.frame(RZ_H_V_ST[RZ_V2_ST_enriched_FDR005, ]), RZ_V2_ST_enriched_taxa)
#write.table(RZ_V2_ST_enriched_FDR005_taxa, file="RZ_V2_ST_enriched_FDR005_taxa.txt", sep="\t")

dev.off()
par(mfrow=c(3,1))
plotMA(RZ_V_F_ST, alpha = 0.05, main="Rhizosphere, Conservation tillage (ST), Rosette vs Flowering")
plotMA(RZ_F_H_ST, alpha = 0.05, main="Rhizosphere, Conservation tillage (ST), Flowering vs Harvesting")
plotMA(RZ_H_V_ST, alpha = 0.05, main="Rhizosphere, Conservation tillage (ST), Harvesting vs Rosette")


########################################################################################################################
# TimePoint effect 
# Rhizosphere 
# Conventional tillage CT

RZ_V_F_CT <- results(Ridhdhi_cds_test, contrast = c("Description", "VegRzConventional", "FlRzConventional"))
RZ_F_H_CT <- results(Ridhdhi_cds_test, contrast = c("Description", "FlRzConventional", "HvRzConventional"))
RZ_H_V_CT <- results(Ridhdhi_cds_test, contrast = c("Description", "HvRzConventional", "VegRzConventional"))

#1
##########
RZ_V_F_CT_FDR_005 <- RZ_V_F_CT[(rownames(RZ_V_F_CT)[which(RZ_V_F_CT$padj < 0.05)]), ]
RZ_V1_CT_enriched <- RZ_V_F_CT[(rownames(RZ_V_F_CT)[which(RZ_V_F_CT$log2FoldChange < 0)]), ]
RZ_V1_CT_enriched_FDR005 <- intersect(rownames(RZ_V_F_CT_FDR_005), rownames(RZ_V1_CT_enriched))
length(RZ_V1_CT_enriched_FDR005)
RZ_V1_CT_enriched_taxa <- Ridhdhi_taxa_ordered[RZ_V1_CT_enriched_FDR005, ]
RZ_V1_CT_enriched_FDR005_taxa <- cbind(as.data.frame(RZ_V_F_CT[RZ_V1_CT_enriched_FDR005, ]), RZ_V1_CT_enriched_taxa)
RZ_F1_CT_enriched <- RZ_V_F_CT[(rownames(RZ_V_F_CT)[which(RZ_V_F_CT$log2FoldChange > 0)]), ]
RZ_F1_CT_enriched_FDR005 <- intersect(rownames(RZ_V_F_CT_FDR_005), rownames(RZ_F1_CT_enriched))
length(RZ_F1_CT_enriched_FDR005)
RZ_F1_CT_enriched_taxa <- Ridhdhi_taxa_ordered[RZ_F1_CT_enriched_FDR005, ]
RZ_F1_CT_enriched_FDR005_taxa <- cbind(as.data.frame(RZ_V_F_CT[RZ_F1_CT_enriched_FDR005, ]), RZ_F1_CT_enriched_taxa)

#2
#########
RZ_F_H_CT_FDR_005 <- RZ_F_H_CT[(rownames(RZ_F_H_CT)[which(RZ_F_H_CT$padj < 0.05)]), ]
RZ_F2_CT_enriched <- RZ_F_H_CT[(rownames(RZ_F_H_CT)[which(RZ_F_H_CT$log2FoldChange < 0)]), ]
RZ_F2_CT_enriched_FDR005 <- intersect(rownames(RZ_F_H_CT_FDR_005), rownames(RZ_F2_CT_enriched))
length(RZ_F2_CT_enriched_FDR005)
RZ_F2_CT_enriched_taxa <- Ridhdhi_taxa_ordered[RZ_F2_CT_enriched_FDR005, ]
RZ_F2_CT_enriched_FDR005_taxa <- cbind(as.data.frame(RZ_F_H_CT[RZ_F2_CT_enriched_FDR005, ]), RZ_F2_CT_enriched_taxa)
RZ_H1_CT_enriched <- RZ_F_H_CT[(rownames(RZ_F_H_CT)[which(RZ_F_H_CT$log2FoldChange > 0)]), ]
RZ_H1_CT_enriched_FDR005 <- intersect(rownames(RZ_F_H_CT_FDR_005), rownames(RZ_H1_CT_enriched))
length(RZ_H1_CT_enriched_FDR005)
RZ_H1_CT_enriched_taxa <- Ridhdhi_taxa_ordered[RZ_H1_CT_enriched_FDR005, ]
RZ_H1_CT_enriched_FDR005_taxa <- cbind(as.data.frame(RZ_F_H_CT[RZ_H1_CT_enriched_FDR005, ]), RZ_H1_CT_enriched_taxa)

#3
#########
RZ_H_V_CT_FDR_005 <- RZ_H_V_CT[(rownames(RZ_H_V_CT)[which(RZ_H_V_CT$padj < 0.05)]), ]
RZ_H2_CT_enriched <- RZ_H_V_CT[(rownames(RZ_H_V_CT)[which(RZ_H_V_CT$log2FoldChange < 0)]), ]
RZ_H2_CT_enriched_FDR005 <- intersect(rownames(RZ_H_V_CT_FDR_005), rownames(RZ_H2_CT_enriched))
length(RZ_H2_CT_enriched_FDR005)
RZ_H2_CT_enriched_taxa <- Ridhdhi_taxa_ordered[RZ_H2_CT_enriched_FDR005, ]
RZ_H2_CT_enriched_FDR005_taxa <- cbind(as.data.frame(RZ_H_V_CT[RZ_H2_CT_enriched_FDR005, ]), RZ_H2_CT_enriched_taxa)
RZ_V2_CT_enriched <- RZ_H_V_CT[(rownames(RZ_H_V_CT)[which(RZ_H_V_CT$log2FoldChange > 0)]), ]
RZ_V2_CT_enriched_FDR005 <- intersect(rownames(RZ_H_V_CT_FDR_005), rownames(RZ_V2_CT_enriched))
length(RZ_V2_CT_enriched_FDR005)
RZ_V2_CT_enriched_taxa <- Ridhdhi_taxa_ordered[RZ_V2_CT_enriched_FDR005, ]
RZ_V2_CT_enriched_FDR005_taxa <- cbind(as.data.frame(RZ_H_V_CT[RZ_V2_CT_enriched_FDR005, ]), RZ_V2_CT_enriched_taxa)

dev.off()
par(mfrow=c(3,1))
plotMA(RZ_V_F_CT, alpha = 0.05, main="Rhizosphere, Conventional tillage (CT), Rosette vs Flowering")
plotMA(RZ_F_H_CT, alpha = 0.05, main="Rhizosphere, Conventional tillage (CT), Flowering vs Harvesting")
plotMA(RZ_H_V_CT, alpha = 0.05, main="Rhizosphere, Conventional tillage (CT), Harvesting vs Rosette")

######################################################################################################################################
# TimePoint effect 
# Bulk soil 
# Conservation tillage ST

BK_V_F_ST <- results(Ridhdhi_cds_test, contrast = c("Description", "VegBkConservation", "FlBkConservation"))
BK_F_H_ST <- results(Ridhdhi_cds_test, contrast = c("Description", "FlBkConservation", "HvBkConservation"))
BK_H_V_ST <- results(Ridhdhi_cds_test, contrast = c("Description", "HvBkConservation", "VegBkConservation"))

#1
##########
BK_V_F_ST_FDR_005 <- BK_V_F_ST[(rownames(BK_V_F_ST)[which(BK_V_F_ST$padj < 0.05)]), ]
BK_V1_ST_enriched <- BK_V_F_ST[(rownames(BK_V_F_ST)[which(BK_V_F_ST$log2FoldChange < 0)]), ]
BK_V1_ST_enriched_FDR005 <- intersect(rownames(BK_V_F_ST_FDR_005), rownames(BK_V1_ST_enriched))
length(BK_V1_ST_enriched_FDR005)
BK_V1_ST_enriched_taxa <- Ridhdhi_taxa_ordered[BK_V1_ST_enriched_FDR005, ]
BK_V1_ST_enriched_FDR005_taxa <- cbind(as.data.frame(BK_V_F_ST[BK_V1_ST_enriched_FDR005, ]), BK_V1_ST_enriched_taxa)
BK_F1_ST_enriched <- BK_V_F_ST[(rownames(BK_V_F_ST)[which(BK_V_F_ST$log2FoldChange > 0)]), ]
BK_F1_ST_enriched_FDR005 <- intersect(rownames(BK_V_F_ST_FDR_005), rownames(BK_F1_ST_enriched))
length(BK_F1_ST_enriched_FDR005)
BK_F1_ST_enriched_taxa <- Ridhdhi_taxa_ordered[BK_F1_ST_enriched_FDR005, ]
BK_F1_ST_enriched_FDR005_taxa <- cbind(as.data.frame(BK_V_F_ST[BK_F1_ST_enriched_FDR005, ]), BK_F1_ST_enriched_taxa)

#2
#########
BK_F_H_ST_FDR_005 <- BK_F_H_ST[(rownames(BK_F_H_ST)[which(BK_F_H_ST$padj < 0.05)]), ]
BK_F2_ST_enriched <- BK_F_H_ST[(rownames(BK_F_H_ST)[which(BK_F_H_ST$log2FoldChange < 0)]), ]
BK_F2_ST_enriched_FDR005 <- intersect(rownames(BK_F_H_ST_FDR_005), rownames(BK_F2_ST_enriched))
length(BK_F2_ST_enriched_FDR005)
BK_F2_ST_enriched_taxa <- Ridhdhi_taxa_ordered[BK_F2_ST_enriched_FDR005, ]
BK_F2_ST_enriched_FDR005_taxa <- cbind(as.data.frame(BK_F_H_ST[BK_F2_ST_enriched_FDR005, ]), BK_F2_ST_enriched_taxa)
BK_H1_ST_enriched <- BK_F_H_ST[(rownames(BK_F_H_ST)[which(BK_F_H_ST$log2FoldChange > 0)]), ]
BK_H1_ST_enriched_FDR005 <- intersect(rownames(BK_F_H_ST_FDR_005), rownames(BK_H1_ST_enriched))
length(BK_H1_ST_enriched_FDR005)
BK_H1_ST_enriched_taxa <- Ridhdhi_taxa_ordered[BK_H1_ST_enriched_FDR005, ]
BK_H1_ST_enriched_FDR005_taxa <- cbind(as.data.frame(BK_F_H_ST[BK_H1_ST_enriched_FDR005, ]), BK_H1_ST_enriched_taxa)

#3
#########
BK_H_V_ST_FDR_005 <- BK_H_V_ST[(rownames(BK_H_V_ST)[which(BK_H_V_ST$padj < 0.05)]), ]
BK_H2_ST_enriched <- BK_H_V_ST[(rownames(BK_H_V_ST)[which(BK_H_V_ST$log2FoldChange < 0)]), ]
BK_H2_ST_enriched_FDR005 <- intersect(rownames(BK_H_V_ST_FDR_005), rownames(BK_H2_ST_enriched))
length(BK_H2_ST_enriched_FDR005)
BK_H2_ST_enriched_taxa <- Ridhdhi_taxa_ordered[BK_H2_ST_enriched_FDR005, ]
BK_H2_ST_enriched_FDR005_taxa <- cbind(as.data.frame(BK_H_V_ST[BK_H2_ST_enriched_FDR005, ]), BK_H2_ST_enriched_taxa)
BK_V2_ST_enriched <- BK_H_V_ST[(rownames(BK_H_V_ST)[which(BK_H_V_ST$log2FoldChange > 0)]), ]
BK_V2_ST_enriched_FDR005 <- intersect(rownames(BK_H_V_ST_FDR_005), rownames(BK_V2_ST_enriched))
length(BK_V2_ST_enriched_FDR005) 
BK_V2_ST_enriched_taxa <- Ridhdhi_taxa_ordered[BK_V2_ST_enriched_FDR005, ]
BK_V2_ST_enriched_FDR005_taxa <- cbind(as.data.frame(BK_H_V_ST[BK_V2_ST_enriched_FDR005, ]), BK_V2_ST_enriched_taxa)

dev.off()
par(mfrow=c(3,1))
plotMA(Bk_V_F_ST, alpha = 0.05, main="Bulk soil, Conservation tillage (CT), Rosette vs Flowering")
plotMA(Bk_F_H_ST, alpha = 0.05, main="Bulk soil, Conservation tillage (CT), Flowering vs Harvesting")
plotMA(Bk_H_V_ST, alpha = 0.05, main="Bulk soil, Conservation tillage (CT), Harvesting vs Rosette")

################################################################################################################################################
# TimePoint effect
# Bulk soil 
# Conventional tillage CT

BK_V_F_CT <- results(Ridhdhi_cds_test, contrast = c("Description", "VegBkConventional", "FlBkConventional"))
BK_F_H_CT <- results(Ridhdhi_cds_test, contrast = c("Description", "FlBkConventional", "HvBkConventional"))
BK_H_V_CT <- results(Ridhdhi_cds_test, contrast = c("Description", "HvBkConventional", "VegBkConventional"))

#1
##########
BK_V_F_CT_FDR_005 <- BK_V_F_CT[(rownames(BK_V_F_CT)[which(BK_V_F_CT$padj < 0.05)]), ]
BK_V1_CT_enriched <- BK_V_F_CT[(rownames(BK_V_F_CT)[which(BK_V_F_CT$log2FoldChange < 0)]), ]
BK_V1_CT_enriched_FDR005 <- intersect(rownames(BK_V_F_CT_FDR_005), rownames(BK_V1_CT_enriched))
length(BK_V1_CT_enriched_FDR005)
BK_V1_CT_enriched_taxa <- Ridhdhi_taxa_ordered[BK_V1_CT_enriched_FDR005, ]
BK_V1_CT_enriched_FDR005_taxa <- cbind(as.data.frame(BK_V_F_CT[BK_V1_CT_enriched_FDR005, ]), BK_V1_CT_enriched_taxa)
BK_F1_CT_enriched <- BK_V_F_CT[(rownames(BK_V_F_CT)[which(BK_V_F_CT$log2FoldChange > 0)]), ]
BK_F1_CT_enriched_FDR005 <- intersect(rownames(BK_V_F_CT_FDR_005), rownames(BK_F1_CT_enriched))
length(BK_F1_CT_enriched_FDR005)
BK_F1_CT_enriched_taxa <- Ridhdhi_taxa_ordered[BK_F1_CT_enriched_FDR005, ]
BK_F1_CT_enriched_FDR005_taxa <- cbind(as.data.frame(BK_V_F_CT[BK_F1_CT_enriched_FDR005, ]), BK_F1_CT_enriched_taxa)

#2
#########
BK_F_H_CT_FDR_005 <- BK_F_H_CT[(rownames(BK_F_H_CT)[which(BK_F_H_CT$padj < 0.05)]), ]
BK_F2_CT_enriched <- BK_F_H_CT[(rownames(BK_F_H_CT)[which(BK_F_H_CT$log2FoldChange < 0)]), ]
BK_F2_CT_enriched_FDR005 <- intersect(rownames(BK_F_H_CT_FDR_005), rownames(BK_F2_CT_enriched))
length(BK_F2_CT_enriched_FDR005)
BK_F2_CT_enriched_taxa <- Ridhdhi_taxa_ordered[BK_F2_CT_enriched_FDR005, ]
BK_F2_CT_enriched_FDR005_taxa <- cbind(as.data.frame(BK_F_H_CT[BK_F2_CT_enriched_FDR005, ]), BK_F2_CT_enriched_taxa)
BK_H1_CT_enriched <- BK_F_H_CT[(rownames(BK_F_H_CT)[which(BK_F_H_CT$log2FoldChange > 0)]), ]
BK_H1_CT_enriched_FDR005 <- intersect(rownames(BK_F_H_CT_FDR_005), rownames(BK_H1_CT_enriched))
length(BK_H1_CT_enriched_FDR005)
BK_H1_CT_enriched_taxa <- Ridhdhi_taxa_ordered[BK_H1_CT_enriched_FDR005, ]
BK_H1_CT_enriched_FDR005_taxa <- cbind(as.data.frame(BK_F_H_CT[BK_H1_CT_enriched_FDR005, ]), BK_H1_CT_enriched_taxa)

#3
#########
BK_H_V_CT_FDR_005 <- BK_H_V_CT[(rownames(BK_H_V_CT)[which(BK_H_V_CT$padj < 0.05)]), ]
BK_H2_CT_enriched <- BK_H_V_CT[(rownames(BK_H_V_CT)[which(BK_H_V_CT$log2FoldChange < 0)]), ]
BK_H2_CT_enriched_FDR005 <- intersect(rownames(BK_H_V_CT_FDR_005), rownames(BK_H2_CT_enriched))
length(BK_H2_CT_enriched_FDR005)
BK_H2_CT_enriched_taxa <- Ridhdhi_taxa_ordered[BK_H2_CT_enriched_FDR005, ]
BK_H2_CT_enriched_FDR005_taxa <- cbind(as.data.frame(BK_H_V_CT[BK_H2_CT_enriched_FDR005, ]), BK_H2_CT_enriched_taxa)
BK_V2_CT_enriched <- BK_H_V_CT[(rownames(BK_H_V_CT)[which(BK_H_V_CT$log2FoldChange > 0)]), ]
BK_V2_CT_enriched_FDR005 <- intersect(rownames(BK_H_V_CT_FDR_005), rownames(BK_V2_CT_enriched))
length(BK_V2_CT_enriched_FDR005)
BK_V2_CT_enriched_taxa <- Ridhdhi_taxa_ordered[BK_V2_CT_enriched_FDR005, ]
BK_V2_CT_enriched_FDR005_taxa <- cbind(as.data.frame(BK_H_V_CT[BK_V2_CT_enriched_FDR005, ]), BK_V2_CT_enriched_taxa)

dev.off()
par(mfrow=c(3,1))
plotMA(Bk_V_F_CT, alpha = 0.05, main="Bulk soil, Conventional tillage (CT), Rosette vs Flowering")
plotMA(Bk_F_H_CT, alpha = 0.05, main="Bulk soil, Conventional tillage (CT), Flowering vs Harvesting")
plotMA(Bk_H_V_CT, alpha = 0.05, main="Bulk soil, Conventional tillage (CT), Harvesting vs Rosette")

#######################################################################################################################################
# Tillage effect (CT vs ST)
# Rhizosphere

RZ_V_CT_ST <- results(Ridhdhi_cds_test, contrast = c("Description", "VegRzConventional", "VegRzConservation"))
RZ_F_CT_ST <- results(Ridhdhi_cds_test, contrast = c("Description", "FlRzConventional", "FlRzConservation"))
RZ_H_CT_ST <- results(Ridhdhi_cds_test, contrast = c("Description", "HvRzConventional", "HvRzConservation"))

#1 Vegetative stage
##########
RZ_V_CT_ST_FDR_005 <- RZ_V_CT_ST[(rownames(RZ_V_CT_ST)[which(RZ_V_CT_ST$padj < 0.05)]), ]
RZ_V_CT_enriched <- RZ_V_CT_ST[(rownames(RZ_V_CT_ST)[which(RZ_V_CT_ST$log2FoldChange < 0)]), ]
RZ_V_CT_enriched_FDR005 <- intersect(rownames(RZ_V_CT_ST_FDR_005), rownames(RZ_V_CT_enriched))
length(RZ_V_CT_enriched_FDR005)
RZ_V_CT_enriched_taxa <- Ridhdhi_taxa_ordered[RZ_V_CT_enriched_FDR005, ]
RZ_V_CT_enriched_FDR005_taxa <- cbind(as.data.frame(RZ_V_CT_ST[RZ_V_CT_enriched_FDR005, ]), RZ_V_CT_enriched_taxa)
RZ_V_ST_enriched <- RZ_V_CT_ST[(rownames(RZ_V_CT_ST)[which(RZ_V_CT_ST$log2FoldChange > 0)]), ]
RZ_V_ST_enriched_FDR005 <- intersect(rownames(RZ_V_CT_ST_FDR_005), rownames(RZ_V_ST_enriched))
length(RZ_V_ST_enriched_FDR005)
RZ_V_ST_enriched_taxa <- Ridhdhi_taxa_ordered[RZ_V_ST_enriched_FDR005, ]
RZ_V_ST_enriched_FDR005_taxa <- cbind(as.data.frame(RZ_V_CT_ST[RZ_V_ST_enriched_FDR005, ]), RZ_V_ST_enriched_taxa)

#2 Flowering stage
#########
RZ_F_CT_ST_FDR_005 <- RZ_F_CT_ST[(rownames(RZ_F_CT_ST)[which(RZ_F_CT_ST$padj < 0.05)]), ]
RZ_F_CT_enriched <- RZ_F_CT_ST[(rownames(RZ_F_CT_ST)[which(RZ_F_CT_ST$log2FoldChange < 0)]), ]
RZ_F_CT_enriched_FDR005 <- intersect(rownames(RZ_F_CT_ST_FDR_005), rownames(RZ_F_CT_enriched))
length(RZ_F_CT_enriched_FDR005)
RZ_F_CT_enriched_taxa <- Ridhdhi_taxa_ordered[RZ_F_CT_enriched_FDR005, ]
RZ_F_CT_enriched_FDR005_taxa <- cbind(as.data.frame(RZ_F_CT_ST[RZ_F_CT_enriched_FDR005, ]), RZ_F_CT_enriched_taxa)
RZ_F_ST_enriched <- RZ_F_CT_ST[(rownames(RZ_F_CT_ST)[which(RZ_F_CT_ST$log2FoldChange > 0)]), ]
RZ_F_ST_enriched_FDR005 <- intersect(rownames(RZ_F_CT_ST_FDR_005), rownames(RZ_F_ST_enriched))
length(RZ_F_ST_enriched_FDR005)
RZ_F_ST_enriched_taxa <- Ridhdhi_taxa_ordered[RZ_F_ST_enriched_FDR005, ]
RZ_F_ST_enriched_FDR005_taxa <- cbind(as.data.frame(RZ_F_CT_ST[RZ_F_ST_enriched_FDR005, ]), RZ_F_ST_enriched_taxa)

#3 Harvesting stage
##########
RZ_H_CT_ST_FDR_005 <- RZ_H_CT_ST[(rownames(RZ_H_CT_ST)[which(RZ_H_CT_ST$padj < 0.05)]), ]
RZ_H_CT_enriched <- RZ_H_CT_ST[(rownames(RZ_H_CT_ST)[which(RZ_H_CT_ST$log2FoldChange < 0)]), ]
RZ_H_CT_enriched_FDR005 <- intersect(rownames(RZ_H_CT_ST_FDR_005), rownames(RZ_H_CT_enriched))
length(RZ_H_CT_enriched_FDR005)
RZ_H_CT_enriched_taxa <- Ridhdhi_taxa_ordered[RZ_H_CT_enriched_FDR005, ]
RZ_H_CT_enriched_FDR005_taxa <- cbind(as.data.frame(RZ_H_CT_ST[RZ_H_CT_enriched_FDR005, ]), RZ_H_CT_enriched_taxa)
RZ_H_ST_enriched <- RZ_H_CT_ST[(rownames(RZ_H_CT_ST)[which(RZ_H_CT_ST$log2FoldChange > 0)]), ]
RZ_H_ST_enriched_FDR005 <- intersect(rownames(RZ_H_CT_ST_FDR_005), rownames(RZ_H_ST_enriched))
length(RZ_H_ST_enriched_FDR005)
RZ_H_ST_enriched_taxa <- Ridhdhi_taxa_ordered[RZ_H_ST_enriched_FDR005, ]
RZ_H_ST_enriched_FDR005_taxa <- cbind(as.data.frame(RZ_H_CT_ST[RZ_H_ST_enriched_FDR005, ]), RZ_H_ST_enriched_taxa)

dev.off()
par(mfrow=c(3,1))
plotMA(RZ_V_CT_ST, alpha = 0.05, main="Rhizosphere, Rosette stage, CT vs ST")
plotMA(RZ_F_CT_ST, alpha = 0.05, main="Rhizosphere, Flowering stage, CT vs ST")
plotMA(RZ_H_CT_ST, alpha = 0.05, main="Rhizosphere, Harvesting stage, CT vs ST")

#######################################################################################################################################
# Tillage effect (CT vs ST)
# Bulk soil 

BK_V_CT_ST <- results(Ridhdhi_cds_test, contrast = c("Description", "VegBkConventional", "VegBkConservation"))
BK_F_CT_ST <- results(Ridhdhi_cds_test, contrast = c("Description", "FlBkConventional", "FlBkConservation"))
BK_H_CT_ST <- results(Ridhdhi_cds_test, contrast = c("Description", "HvBkConventional", "HvBkConservation"))

#1 Vegetative stage
##########
BK_V_CT_ST_FDR_005 <- BK_V_CT_ST[(rownames(BK_V_CT_ST)[which(BK_V_CT_ST$padj < 0.05)]), ]
BK_V_CT_enriched <- BK_V_CT_ST[(rownames(BK_V_CT_ST)[which(BK_V_CT_ST$log2FoldChange < 0)]), ]
BK_V_CT_enriched_FDR005 <- intersect(rownames(BK_V_CT_ST_FDR_005), rownames(BK_V_CT_enriched))
length(BK_V_CT_enriched_FDR005)
BK_V_CT_enriched_taxa <- Ridhdhi_taxa_ordered[BK_V_CT_enriched_FDR005, ]
BK_V_CT_enriched_FDR005_taxa <- cbind(as.data.frame(BK_V_CT_ST[BK_V_CT_enriched_FDR005, ]), BK_V_CT_enriched_taxa)
BK_V_ST_enriched <- BK_V_CT_ST[(rownames(BK_V_CT_ST)[which(BK_V_CT_ST$log2FoldChange > 0)]), ]
BK_V_ST_enriched_FDR005 <- intersect(rownames(BK_V_CT_ST_FDR_005), rownames(BK_V_ST_enriched))
length(BK_V_ST_enriched_FDR005)
BK_V_ST_enriched_taxa <- Ridhdhi_taxa_ordered[BK_V_ST_enriched_FDR005, ]
BK_V_ST_enriched_FDR005_taxa <- cbind(as.data.frame(BK_V_CT_ST[BK_V_ST_enriched_FDR005, ]), BK_V_ST_enriched_taxa)

#2 Flowering stage
#########
BK_F_CT_ST_FDR_005 <- BK_F_CT_ST[(rownames(BK_F_CT_ST)[which(BK_F_CT_ST$padj < 0.05)]), ]
BK_F_CT_enriched <- BK_F_CT_ST[(rownames(BK_F_CT_ST)[which(BK_F_CT_ST$log2FoldChange < 0)]), ]
BK_F_CT_enriched_FDR005 <- intersect(rownames(BK_F_CT_ST_FDR_005), rownames(BK_F_CT_enriched))
length(BK_F_CT_enriched_FDR005)
BK_F_CT_enriched_taxa <- Ridhdhi_taxa_ordered[BK_F_CT_enriched_FDR005, ]
BK_F_CT_enriched_FDR005_taxa <- cbind(as.data.frame(BK_F_CT_ST[BK_F_CT_enriched_FDR005, ]), BK_F_CT_enriched_taxa)
BK_F_ST_enriched <- BK_F_CT_ST[(rownames(BK_F_CT_ST)[which(BK_F_CT_ST$log2FoldChange > 0)]), ]
BK_F_ST_enriched_FDR005 <- intersect(rownames(BK_F_CT_ST_FDR_005), rownames(BK_F_ST_enriched))
length(BK_F_ST_enriched_FDR005)
BK_F_ST_enriched_taxa <- Ridhdhi_taxa_ordered[BK_F_ST_enriched_FDR005, ]
BK_F_ST_enriched_FDR005_taxa <- cbind(as.data.frame(BK_F_CT_ST[BK_F_ST_enriched_FDR005, ]), BK_F_ST_enriched_taxa)

#3 Harvesting stage
##########
BK_H_CT_ST_FDR_005 <- BK_H_CT_ST[(rownames(BK_H_CT_ST)[which(BK_H_CT_ST$padj < 0.05)]), ]
BK_H_CT_enriched <- BK_H_CT_ST[(rownames(BK_H_CT_ST)[which(BK_H_CT_ST$log2FoldChange < 0)]), ]
BK_H_CT_enriched_FDR005 <- intersect(rownames(BK_H_CT_ST_FDR_005), rownames(BK_H_CT_enriched))
length(BK_H_CT_enriched_FDR005)
BK_H_CT_enriched_taxa <- Ridhdhi_taxa_ordered[BK_H_CT_enriched_FDR005, ]
BK_H_CT_enriched_FDR005_taxa <- cbind(as.data.frame(BK_H_CT_ST[BK_H_CT_enriched_FDR005, ]), BK_H_CT_enriched_taxa)
BK_H_ST_enriched <- BK_H_CT_ST[(rownames(BK_H_CT_ST)[which(BK_H_CT_ST$log2FoldChange > 0)]), ]
BK_H_ST_enriched_FDR005 <- intersect(rownames(BK_H_CT_ST_FDR_005), rownames(BK_H_ST_enriched))
length(BK_H_ST_enriched_FDR005)
BK_H_ST_enriched_taxa <- Ridhdhi_taxa_ordered[BK_H_ST_enriched_FDR005, ]
BK_H_ST_enriched_FDR005_taxa <- cbind(as.data.frame(BK_H_CT_ST[BK_H_ST_enriched_FDR005, ]), BK_H_ST_enriched_taxa)

dev.off()
par(mfrow=c(3,1))
plotMA(Bk_V_CT_ST, alpha = 0.05, main="Bulk soil, Rosette stage, CT vs ST")
plotMA(Bk_F_CT_ST, alpha = 0.05, main="Bulk soil, Flowering stage, CT vs ST")
plotMA(Bk_H_CT_ST, alpha = 0.05, main="Bulk soil, Harvesting stage, CT vs ST")

##########################################################################################################################

# Compartment effect (Bulk vs Rhizosphere) 
# Conservation tillage ST

ST_V_BK_RZ <- results(Ridhdhi_cds_test, contrast = c("Description", "VegBkConservation", "VegRzConservation"))
ST_F_BK_RZ <- results(Ridhdhi_cds_test, contrast = c("Description", "FlBkConservation", "FlRzConservation"))
ST_H_BK_RZ <- results(Ridhdhi_cds_test, contrast = c("Description", "HvBkConservation", "HvRzConservation"))

#1 Vegetative stage
##########
ST_V_BK_RZ_FDR_005 <- ST_V_BK_RZ[(rownames(ST_V_BK_RZ)[which(ST_V_BK_RZ$padj < 0.05)]), ]
ST_V_BK_enriched <- ST_V_BK_RZ[(rownames(ST_V_BK_RZ)[which(ST_V_BK_RZ$log2FoldChange < 0)]), ]
ST_V_BK_enriched_FDR005 <- intersect(rownames(ST_V_BK_RZ_FDR_005), rownames(ST_V_BK_enriched))
length(ST_V_BK_enriched_FDR005)
ST_V_BK_enriched_taxa <- Ridhdhi_taxa_ordered[ST_V_BK_enriched_FDR005, ]
ST_V_BK_enriched_FDR005_taxa <- cbind(as.data.frame(ST_V_BK_RZ[ST_V_BK_enriched_FDR005, ]), ST_V_BK_enriched_taxa)
ST_V_RZ_enriched <- ST_V_BK_RZ[(rownames(ST_V_BK_RZ)[which(ST_V_BK_RZ$log2FoldChange > 0)]), ]
ST_V_RZ_enriched_FDR005 <- intersect(rownames(ST_V_BK_RZ_FDR_005), rownames(ST_V_RZ_enriched))
length(ST_V_RZ_enriched_FDR005)
ST_V_RZ_enriched_taxa <- Ridhdhi_taxa_ordered[ST_V_RZ_enriched_FDR005, ]
ST_V_RZ_enriched_FDR005_taxa <- cbind(as.data.frame(ST_V_BK_RZ[ST_V_RZ_enriched_FDR005, ]), ST_V_RZ_enriched_taxa)

#2 Flowering stage
##########
ST_F_BK_RZ_FDR_005 <- ST_F_BK_RZ[(rownames(ST_F_BK_RZ)[which(ST_F_BK_RZ$padj < 0.05)]), ]
ST_F_BK_enriched <- ST_F_BK_RZ[(rownames(ST_F_BK_RZ)[which(ST_F_BK_RZ$log2FoldChange < 0)]), ]
ST_F_BK_enriched_FDR005 <- intersect(rownames(ST_F_BK_RZ_FDR_005), rownames(ST_F_BK_enriched))
length(ST_F_BK_enriched_FDR005)
ST_F_BK_enriched_taxa <- Ridhdhi_taxa_ordered[ST_F_BK_enriched_FDR005, ]
ST_F_BK_enriched_FDR005_taxa <- cbind(as.data.frame(ST_F_BK_RZ[ST_F_BK_enriched_FDR005, ]), ST_F_BK_enriched_taxa)
ST_F_RZ_enriched <- ST_F_BK_RZ[(rownames(ST_F_BK_RZ)[which(ST_F_BK_RZ$log2FoldChange > 0)]), ]
ST_F_RZ_enriched_FDR005 <- intersect(rownames(ST_F_BK_RZ_FDR_005), rownames(ST_F_RZ_enriched))
length(ST_F_RZ_enriched_FDR005)
ST_F_RZ_enriched_taxa <- Ridhdhi_taxa_ordered[ST_F_RZ_enriched_FDR005, ]
ST_F_RZ_enriched_FDR005_taxa <- cbind(as.data.frame(ST_F_BK_RZ[ST_F_RZ_enriched_FDR005, ]), ST_F_RZ_enriched_taxa)

#3 Harvesting stage
##########
ST_H_BK_RZ_FDR_005 <- ST_H_BK_RZ[(rownames(ST_H_BK_RZ)[which(ST_H_BK_RZ$padj < 0.05)]), ]
ST_H_BK_enriched <- ST_H_BK_RZ[(rownames(ST_H_BK_RZ)[which(ST_H_BK_RZ$log2FoldChange < 0)]), ]
ST_H_BK_enriched_FDR005 <- intersect(rownames(ST_H_BK_RZ_FDR_005), rownames(ST_H_BK_enriched))
length(ST_H_BK_enriched_FDR005)
ST_H_BK_enriched_taxa <- Ridhdhi_taxa_ordered[ST_H_BK_enriched_FDR005, ]
ST_H_BK_enriched_FDR005_taxa <- cbind(as.data.frame(ST_H_BK_RZ[ST_H_BK_enriched_FDR005, ]), ST_H_BK_enriched_taxa)
#write.table(ST_H_BK_enriched_FDR005_taxa, file="ST_H_BK_enriched_FDR005_taxa.txt", sep="\t")
ST_H_RZ_enriched <- ST_H_BK_RZ[(rownames(ST_H_BK_RZ)[which(ST_H_BK_RZ$log2FoldChange > 0)]), ]
ST_H_RZ_enriched_FDR005 <- intersect(rownames(ST_H_BK_RZ_FDR_005), rownames(ST_H_RZ_enriched))
length(ST_H_RZ_enriched_FDR005)
ST_H_RZ_enriched_taxa <- Ridhdhi_taxa_ordered[ST_H_RZ_enriched_FDR005, ]
ST_H_RZ_enriched_FDR005_taxa <- cbind(as.data.frame(ST_H_BK_RZ[ST_H_RZ_enriched_FDR005, ]), ST_H_RZ_enriched_taxa)

dev.off()
par(mfrow=c(3,1))
plotMA(ST_V_BK_RZ, alpha = 0.05, main="Conservation tillage (ST), Rosette stage, Bulk soil vs Rhizosphere")
plotMA(ST_F_BK_RZ, alpha = 0.05, main="Conservation tillage (ST), Flowering stage, Bulk soil vs Rhizosphere")
plotMA(ST_H_BK_RZ, alpha = 0.05, main="Conservation tillage (ST), Harvesting stage, Bulk soil vs Rhizosphere")

####################################################################################################################
# Compartment effect (Bulk vs Rhizosphere)
# Conventional tillage CT

CT_V_BK_RZ <- results(Ridhdhi_cds_test, contrast = c("Description", "VegBkConventional", "VegRzConventional"))
CT_F_BK_RZ <- results(Ridhdhi_cds_test, contrast = c("Description", "FlBkConventional", "FlRzConventional"))
CT_H_BK_RZ <- results(Ridhdhi_cds_test, contrast = c("Description", "HvBkConventional", "HvRzConventional"))

#1 Vegetative stage
##########
CT_V_BK_RZ_FDR_005 <- CT_V_BK_RZ[(rownames(CT_V_BK_RZ)[which(CT_V_BK_RZ$padj < 0.05)]), ]
CT_V_BK_enriched <- CT_V_BK_RZ[(rownames(CT_V_BK_RZ)[which(CT_V_BK_RZ$log2FoldChange < 0)]), ]
CT_V_BK_enriched_FDR005 <- intersect(rownames(CT_V_BK_RZ_FDR_005), rownames(CT_V_BK_enriched))
length(CT_V_BK_enriched_FDR005)
CT_V_BK_enriched_taxa <- Ridhdhi_taxa_ordered[CT_V_BK_enriched_FDR005, ]
CT_V_BK_enriched_FDR005_taxa <- cbind(as.data.frame(CT_V_BK_RZ[CT_V_BK_enriched_FDR005, ]), CT_V_BK_enriched_taxa)
CT_V_RZ_enriched <- CT_V_BK_RZ[(rownames(CT_V_BK_RZ)[which(CT_V_BK_RZ$log2FoldChange > 0)]), ]
CT_V_RZ_enriched_FDR005 <- intersect(rownames(CT_V_BK_RZ_FDR_005), rownames(CT_V_RZ_enriched))
length(CT_V_RZ_enriched_FDR005)
CT_V_RZ_enriched_taxa <- Ridhdhi_taxa_ordered[CT_V_RZ_enriched_FDR005, ]
CT_V_RZ_enriched_FDR005_taxa <- cbind(as.data.frame(CT_V_BK_RZ[CT_V_RZ_enriched_FDR005, ]), CT_V_RZ_enriched_taxa)

#2 Flowering stage
##########
CT_F_BK_RZ_FDR_005 <- CT_F_BK_RZ[(rownames(CT_F_BK_RZ)[which(CT_F_BK_RZ$padj < 0.05)]), ]
CT_F_BK_enriched <- CT_F_BK_RZ[(rownames(CT_F_BK_RZ)[which(CT_F_BK_RZ$log2FoldChange < 0)]), ]
CT_F_BK_enriched_FDR005 <- intersect(rownames(CT_F_BK_RZ_FDR_005), rownames(CT_F_BK_enriched))
length(CT_F_BK_enriched_FDR005)
CT_F_BK_enriched_taxa <- Ridhdhi_taxa_ordered[CT_F_BK_enriched_FDR005, ]
CT_F_BK_enriched_FDR005_taxa <- cbind(as.data.frame(CT_F_BK_RZ[CT_F_BK_enriched_FDR005, ]), CT_F_BK_enriched_taxa)
CT_F_RZ_enriched <- CT_F_BK_RZ[(rownames(CT_F_BK_RZ)[which(CT_F_BK_RZ$log2FoldChange > 0)]), ]
CT_F_RZ_enriched_FDR005 <- intersect(rownames(CT_F_BK_RZ_FDR_005), rownames(CT_F_RZ_enriched))
length(CT_F_RZ_enriched_FDR005)
CT_F_RZ_enriched_taxa <- Ridhdhi_taxa_ordered[CT_F_RZ_enriched_FDR005, ]
CT_F_RZ_enriched_FDR005_taxa <- cbind(as.data.frame(CT_F_BK_RZ[CT_F_RZ_enriched_FDR005, ]), CT_F_RZ_enriched_taxa)

#3 Harvesting stage
##########
CT_H_BK_RZ_FDR_005 <- CT_H_BK_RZ[(rownames(CT_H_BK_RZ)[which(CT_H_BK_RZ$padj < 0.05)]), ]
CT_H_BK_enriched <- CT_H_BK_RZ[(rownames(CT_H_BK_RZ)[which(CT_H_BK_RZ$log2FoldChange < 0)]), ]
CT_H_BK_enriched_FDR005 <- intersect(rownames(CT_H_BK_RZ_FDR_005), rownames(CT_H_BK_enriched))
length(CT_H_BK_enriched_FDR005)
CT_H_BK_enriched_taxa <- Ridhdhi_taxa_ordered[CT_H_BK_enriched_FDR005, ]
CT_H_BK_enriched_FDR005_taxa <- cbind(as.data.frame(CT_H_BK_RZ[CT_H_BK_enriched_FDR005, ]), CT_H_BK_enriched_taxa)
CT_H_RZ_enriched <- CT_H_BK_RZ[(rownames(CT_H_BK_RZ)[which(CT_H_BK_RZ$log2FoldChange > 0)]), ]
CT_H_RZ_enriched_FDR005 <- intersect(rownames(CT_H_BK_RZ_FDR_005), rownames(CT_H_RZ_enriched))
length(CT_H_RZ_enriched_FDR005)
CT_H_RZ_enriched_taxa <- Ridhdhi_taxa_ordered[CT_H_RZ_enriched_FDR005, ]
CT_H_RZ_enriched_FDR005_taxa <- cbind(as.data.frame(CT_H_BK_RZ[CT_H_RZ_enriched_FDR005, ]), CT_H_RZ_enriched_taxa)

dev.off()
par(mfrow=c(3,1))
plotMA(CT_V_BK_RZ, alpha = 0.05, main="Conventional tillage (CT), Rosette stage, Bulk soil vs Rhizosphere")
plotMA(CT_F_BK_RZ, alpha = 0.05, main="Conventional tillage (CT), Flowering stage, Bulk soil vs Rhizosphere")
plotMA(CT_H_BK_RZ, alpha = 0.05, main="Conventional tillage (CT), Harvesting stage, Bulk soil vs Rhizosphere")

##################################################################################################################################################################
#Relative abundance of bacterial Phyla (Figure 4)
##############################################################################################################


Acidobacteria <- read.delim("Acidobacteria_R_B.txt", sep = "\t", header = TRUE)
attach(Acidobacteria)
f = ordered(Description, levels=c("B_V_CT", "R_V_CT", "B_V_ST", "R_V_ST", "B_F_CT", "R_F_CT", "B_F_ST", "R_F_ST", "B_H_CT", "R_H_CT", "B_H_ST", "R_H_ST"))
boxplot(Abundance ~ f, xlab = "TimePoint", ylab = "Abundance", main = "Acidobacteria", col=(c("olivedrab4", "olivedrab1", "olivedrab4", "olivedrab1", "orange", "yellow", "orange", "yellow", "tan4", "khaki1", "tan4", "khaki1")))
anova(Acidobacteria)

Actinobacteria <- read.delim("Actinobacteria_R_B.txt", sep = "\t", header = TRUE)
attach(Actinobacteria)
f1 = ordered(Description, levels=c("B_V_CT", "R_V_CT", "B_V_ST", "R_V_ST", "B_F_CT", "R_F_CT", "B_F_ST", "R_F_ST", "B_H_CT", "R_H_CT", "B_H_ST", "R_H_ST"))
boxplot(Abundance ~ f1, xlab = "TimePoint", ylab = "Abundance", main = "Actinobacteria", col=(c("olivedrab4", "olivedrab1", "olivedrab4", "olivedrab1", "orange", "yellow", "orange", "yellow", "tan4", "khaki1", "tan4", "khaki1")))
anova(Actinobacteria)

Bacteroidetes <- read.delim("Bacteroidetes_R_B.txt", sep = "\t", header = TRUE)
attach(Bacteroidetes)
f2 = ordered(Description, levels=c("B_V_CT", "R_V_CT", "B_V_ST", "R_V_ST", "B_F_CT", "R_F_CT", "B_F_ST", "R_F_ST", "B_H_CT", "R_H_CT", "B_H_ST", "R_H_ST"))
boxplot(Abundance ~ f2, xlab = "TimePoint", ylab = "Abundance", main = "Bacteroidetes", col=(c("olivedrab4", "olivedrab1", "olivedrab4", "olivedrab1", "orange", "yellow", "orange", "yellow", "tan4", "khaki1", "tan4", "khaki1")))
anova(Bacteroidetes)

Firmicutes <- read.delim("Firmicutes_R_B.txt", sep = "\t", header = TRUE)
attach(Firmicutes)
f3 = ordered(Description, levels=c("B_V_CT", "R_V_CT", "B_V_ST", "R_V_ST", "B_F_CT", "R_F_CT", "B_F_ST", "R_F_ST", "B_H_CT", "R_H_CT", "B_H_ST", "R_H_ST"))
boxplot(Abundance ~ f3, xlab = "TimePoint", ylab = "Abundance", main = "Firmicutes", col=(c("olivedrab4", "olivedrab1", "olivedrab4", "olivedrab1", "orange", "yellow", "orange", "yellow", "tan4", "khaki1", "tan4", "khaki1")))
anova(Firmicutes)

Proteobacteria <- read.delim("Proteobacteria_R_B.txt", sep = "\t", header = TRUE)
attach(Proteobacteria)
f4 = ordered(Description, levels=c("B_V_CT", "R_V_CT", "B_V_ST", "R_V_ST", "B_F_CT", "R_F_CT", "B_F_ST", "R_F_ST", "B_H_CT", "R_H_CT", "B_H_ST", "R_H_ST"))
boxplot(Abundance ~ f4, xlab = "TimePoint", ylab = "Abundance", main = "Proteobacteria", col=(c("olivedrab4", "olivedrab1", "olivedrab4", "olivedrab1", "orange", "yellow", "orange", "yellow", "tan4", "khaki1", "tan4", "khaki1")))
anova(Proteobacteria)

Chloroflexi <- read.delim("Chloroflexi_R_B.txt", sep = "\t", header = TRUE)
attach(Chloroflexi)
f5 = ordered(Description, levels=c("B_V_CT", "R_V_CT", "B_V_ST", "R_V_ST", "B_F_CT", "R_F_CT", "B_F_ST", "R_F_ST", "B_H_CT", "R_H_CT", "B_H_ST", "R_H_ST"))
boxplot(Abundance ~ f5, xlab = "TimePoint", ylab = "Abundance", main = "Chloroflexi", col=(c("olivedrab4", "olivedrab1", "olivedrab4", "olivedrab1", "orange", "yellow", "orange", "yellow", "tan4", "khaki1", "tan4", "khaki1")))
anova(Chloroflexi)

TM7 <- read.delim("TM7_R_B.txt", sep = "\t", header = TRUE)
attach(TM7)
f6 = ordered(Description, levels=c("B_V_CT", "R_V_CT", "B_V_ST", "R_V_ST", "B_F_CT", "R_F_CT", "B_F_ST", "R_F_ST", "B_H_CT", "R_H_CT", "B_H_ST", "R_H_ST"))
boxplot(Abundance ~ f6, xlab = "TimePoint", ylab = "Abundance", main = "TM7", col=(c("olivedrab4", "olivedrab1", "olivedrab4", "olivedrab1", "orange", "yellow", "orange", "yellow", "tan4", "khaki1", "tan4", "khaki1")))
anova(TM7)

Verrucomicrobia <- read.delim("Verrucomicrobia_R_B.txt", sep = "\t", header = TRUE)
attach(Verrucomicrobia)
f7 = ordered(Description, levels=c("B_V_CT", "R_V_CT", "B_V_ST", "R_V_ST", "B_F_CT", "R_F_CT", "B_F_ST", "R_F_ST", "B_H_CT", "R_H_CT", "B_H_ST", "R_H_ST"))
boxplot(Abundance ~ f7, xlab = "TimePoint", ylab = "Abundance", main = "Verrucomicrobia", col=(c("olivedrab4", "olivedrab1", "olivedrab4", "olivedrab1", "orange", "yellow", "orange", "yellow", "tan4", "khaki1", "tan4", "khaki1")))
anova(Verrucomicrobia)

Planctomycetes <- read.delim("Planctomycetes_R_B.txt", sep = "\t", header = TRUE)
attach(Planctomycetes)
f8 = ordered(Description, levels=c("B_V_CT", "R_V_CT", "B_V_ST", "R_V_ST", "B_F_CT", "R_F_CT", "B_F_ST", "R_F_ST", "B_H_CT", "R_H_CT", "B_H_ST", "R_H_ST"))
boxplot(Abundance ~ f8, xlab = "TimePoint", ylab = "Abundance", main = "Planctomycetes", col=(c("olivedrab4", "olivedrab1", "olivedrab4", "olivedrab1", "orange", "yellow", "orange", "yellow", "tan4", "khaki1", "tan4", "khaki1")))
anova(Planctomycetes)

Cynobacteria <- read.delim("Cynobacteria_R_B.txt", sep = "\t", header = TRUE)
attach(Cynobacteria)
f9 = ordered(Description, levels=c("B_V_CT", "R_V_CT", "B_V_ST", "R_V_ST", "B_F_CT", "R_F_CT", "B_F_ST", "R_F_ST", "B_H_CT", "R_H_CT", "B_H_ST", "R_H_ST"))
boxplot(Abundance ~ f9, xlab = "TimePoint", ylab = "Abundance", main = "Cynobacteria", col=(c("olivedrab4", "olivedrab1", "olivedrab4", "olivedrab1", "orange", "yellow", "orange", "yellow", "tan4", "khaki1", "tan4", "khaki1")))
anova(Cynobacteria)

Gemmatimonadetes <- read.delim("Gemmatimonadetes_R_B.txt", sep = "\t", header = TRUE)
attach(Gemmatimonadetes)
f10 = ordered(Description, levels=c("B_V_CT", "R_V_CT", "B_V_ST", "R_V_ST", "B_F_CT", "R_F_CT", "B_F_ST", "R_F_ST", "B_H_CT", "R_H_CT", "B_H_ST", "R_H_ST"))
boxplot(Abundance ~ f10, xlab = "TimePoint", ylab = "Abundance", main = "Gemmatimonadetes", col=(c("olivedrab4", "olivedrab1", "olivedrab4", "olivedrab1", "orange", "yellow", "orange", "yellow", "tan4", "khaki1", "tan4", "khaki1")))
anova(Gemmatimonadetes)

##############################################################################################################
# End
####################################################################################################################################################################

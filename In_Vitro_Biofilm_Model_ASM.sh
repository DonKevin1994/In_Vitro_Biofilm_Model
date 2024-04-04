#!/bin/bash

conda activate qiime2-2022.8
cd /mnt/c/Users/Kevin/QIIME_workflow/Peter_V/decontam/Peter_V_Sequencing_Analysis

#27/05/2021

#Split the ASM samples from the merged table 

qiime feature-table filter-samples \
--i-table Peter_V_merged_table.qza \
--m-metadata-file Peter_V_merged_map.txt \
--p-where "BarcodeSequence='TBD'" \
--o-filtered-table Peter_V_ASM_table.qza

qiime feature-table summarize \
--i-table Peter_V_ASM_table.qza \
--o-visualization Peter_V_ASM_table.qzv \
--m-sample-metadata-file Peter_V_merged_map.txt  

#Split the SHI samples from the merged table 

qiime feature-table filter-samples \
--i-table Peter_V_merged_table.qza \
--m-metadata-file Peter_V_merged_map.txt \
--p-where "BarcodeSequence='TBD'" \
--p-exclude-ids \
--o-filtered-table Peter_V_SHI_table.qza

qiime feature-table summarize \
--i-table Peter_V_SHI_table.qza \
--o-visualization Peter_V_SHI_table.qzv \
--m-sample-metadata-file Peter_V_merged_map.txt 

######### ASM SAMPLES ###########

#Run core-metrics 

# Run core metrics. Sampling depth is critical here.  Look at features in 
#Peter_V_ASM_table.qzv to see where to subsample to maximize sequences but limit number of samples excluded
# In this data, 900 sequences includes most samples 

  qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted_tree.qza \
  --i-table Peter_V_ASM_table.qza \
  --p-sampling-depth 770 \
  --m-metadata-file Peter_V_merged_map.txt \
  --output-dir Peter_V_ASM_CoreDiversity_Results

# Run stats to see if metadata explains any beta diversity characteristics
# This is for beta diversity 

  #qiime diversity beta-group-significance \
  #--i-distance-matrix Peter_V_ASM_CoreDiversity_Results/unweighted_unifrac_distance_matrix.qza \
  #--m-metadata-column PlaqueDayInProcess \
  #--m-metadata-file Peter_V_merged_map.txt \
  #--o-visualization Peter_V_ASM_CoreDiversity_Results/unweighted_unifrac_significance_PlaqueDayInProcess.qzv

  qiime diversity beta-group-significance \
  --i-distance-matrix Peter_V_ASM_CoreDiversity_Results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-column ControlvsBiological \
  --m-metadata-file Peter_V_merged_map.txt \
  --o-visualization Peter_V_ASM_CoreDiversity_Results/unweighted_unifrac_significance_ControlvsBiological.qzv

# ALPHA RAREFACTION
# Look at alpha rarefaction of your data

qiime diversity alpha-rarefaction \
  --i-table Peter_V_ASM_table.qza \
  --i-phylogeny rooted_tree.qza \
  --p-max-depth 770 \
  --m-metadata-file Peter_V_merged_map.txt \
  --o-visualization Peter_V_ASM_CoreDiversity_Results/alpha_rarefaction.qzv

# Run stats to see if metadata explains any alpha or beta diversity characteristics
# This is for alpha diversity

  qiime diversity alpha \
  --i-table Peter_V_ASM_table.qza \
  --p-metric observed_features \
  --o-alpha-diversity observed_features_vector.qza

  ^use this only if there is no observed_features_vector or observed_otus_vector in the CoreDiversity_Results folder 

  qiime diversity alpha-group-significance \
  --i-alpha-diversity Peter_V_ASM_CoreDiversity_Results/observed_features_vector.qza \
  --m-metadata-file Peter_V_merged_map.txt \
  --o-visualization Peter_V_ASM_CoreDiversity_Results/observed_features_significance.qzv

############################
# Taxonomic analysis of your data (first classify features and then visualise)

  qiime feature-classifier classify-sklearn \
  --i-classifier gg-13-8-99-515-806-nb-classifier.qza \
  --i-reads merged_repseqs.qza \
  --o-classification taxonomy.qza

qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy.qzv

# Building the plots to visualize

  qiime taxa barplot \
  --i-table Peter_V_ASM_table.qza \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file Peter_V_merged_map.txt \
  --o-visualization Peter_V_ASM_Taxonomy_BarPlots.qzv

########################

# Filtering metadata column - ControlvsBiological_ASM

qiime feature-table filter-samples \
--i-table Peter_V_merged_table.qza \
--m-metadata-file Peter_V_merged_map.txt \
--p-where "ControlvsBiological_ASM='na'" \
--p-exclude-ids \
--o-filtered-table Peter_V_ControlvsBiological_ASM_table.qza

qiime feature-table summarize \
--i-table Peter_V_ControlvsBiological_ASM_table.qza \
--o-visualization Peter_V_ControlvsBiological_ASM_table.qzv \
--m-sample-metadata-file Peter_V_merged_map.txt 

# Ran Core metrics

qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted_tree.qza \
  --i-table Peter_V_ControlvsBiological_ASM_table.qza \
  --p-sampling-depth 770 \
  --m-metadata-file Peter_V_merged_map.txt \
  --output-dir Peter_V_ControlvsBiological_ASM_CoreDiversity_Results

# Run stats to see if metadata explains any beta diversity characteristics
# This is for beta diversity 
# ControlvsBiological_ASM gives has different types of controls and biological samples (samples are not just labelled 
#Control or Biological)

  qiime diversity beta-group-significance \
  --i-distance-matrix Peter_V_ControlvsBiological_ASM_CoreDiversity_Results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-column ControlvsBiological_ASM \
  --m-metadata-file Peter_V_merged_map.txt \
  --o-visualization Peter_V_ControlvsBiological_ASM_CoreDiversity_Results/unweighted_unifrac_significance_ControlvsBiological_ASM.qzv

# ControlvsBiological_ASM_General metadata column has samples labelled just as Control or Biological

 qiime diversity beta-group-significance \
  --i-distance-matrix Peter_V_ControlvsBiological_ASM_CoreDiversity_Results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-column ControlvsBiological_ASM_General \
  --m-metadata-file Peter_V_merged_map.txt \
  --o-visualization Peter_V_ControlvsBiological_ASM_CoreDiversity_Results/unweighted_unifrac_significance_ControlvsBiological_ASM_General.qzv

# ALPHA RAREFACTION
# Look at alpha rarefaction of your data

qiime diversity alpha-rarefaction \
  --i-table Peter_V_ControlvsBiological_ASM_table.qza \
  --i-phylogeny rooted_tree.qza \
  --p-max-depth 770 \
  --m-metadata-file Peter_V_merged_map.txt \
  --o-visualization Peter_V_ControlvsBiological_ASM_CoreDiversity_Results/alpha_rarefaction.qzv

# Run stats to see if metadata explains any alpha or beta diversity characteristics
# This is for alpha diversity

qiime diversity alpha-group-significance \
  --i-alpha-diversity Peter_V_ControlvsBiological_ASM_CoreDiversity_Results/observed_features_vector.qza \
  --m-metadata-file Peter_V_merged_map.txt \
  --o-visualization Peter_V_ControlvsBiological_ASM_CoreDiversity_Results/observed_features_significance.qzv

############################
# Taxonomic analysis of your data (first classify features and then visualise)

  #qiime feature-classifier classify-sklearn \
  #--i-classifier gg-13-8-99-515-806-nb-classifier.qza \
  #--i-reads merged_repseqs.qza \
  #--o-classification taxonomy.qza

#qiime metadata tabulate \
  #--m-input-file taxonomy.qza \
  #--o-visualization taxonomy.qzv

# Building the plots to visualize

  qiime taxa barplot \
  --i-table Peter_V_ControlvsBiological_ASM_table.qza \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file Peter_V_merged_map.txt \
  --o-visualization Peter_V_ControlvsBiological_ASM_Taxonomy_BarPlots.qzv 

# Filtering table to remove all control samples 

  qiime feature-table filter-samples \
  --i-table Peter_V_ASM_table.qza \
  --m-metadata-file Peter_V_merged_map.txt \
  --p-where "ControlvsBiological='Control'" \
  --p-exclude-ids \
  --o-filtered-table Peter_V_ASM_table_Biologicalsamplesonly \

  qiime feature-table summarize \
  --i-table Peter_V_ASM_table_Biologicalsamplesonly.qza \
  --o-visualization Peter_V_ASM_table_Biologicalsamplesonly.qzv \
  --m-sample-metadata-file Peter_V_merged_map.txt

   qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted_tree.qza \
  --i-table Peter_V_ASM_table_Biologicalsamplesonly.qza \
  --p-sampling-depth 900 \
  --m-metadata-file Peter_V_merged_map.txt \
  --output-dir Peter_V_ASM_CoreDiversity_Results_Biologicalsamplesonly

  qiime diversity beta-group-significance \
  --i-distance-matrix Peter_V_ASM_CoreDiversity_Results_Biologicalsamplesonly/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-column PlankvsPlaq \
  --m-metadata-file Peter_V_merged_map.txt \
  --o-visualization Peter_V_ASM_CoreDiversity_Results_Biologicalsamplesonly/unweighted_unifrac_significance_PlankvsPlaq.qzv

  qiime diversity alpha-group-significance \
  --i-alpha-diversity Peter_V_ASM_CoreDiversity_Results_Biologicalsamplesonly/observed_features_vector.qza \
  --m-metadata-file Peter_V_merged_map.txt \
  --o-visualization Peter_V_ASM_CoreDiversity_Results_Biologicalsamplesonly/observed_features_significance.qzv

  ######## Decontam for ASM EBC1 ###########

  #Ran on Rstudio Server 2022.07.2 (R script was uploaded separately on to Github)

  # Installing packages 

library(tidyverse)
library(readr)
library(phyloseq)
library(ggplot2)
library(decontam)
library(scales)
library(qiime2R)

# Importing QIIME2 files for Decontam 

metadata<-read_tsv("Peter_V_merged_map.txt")

ASVs <- read_qza("Peter_V_ASM_table.qza")

greengenes_taxonomy <- read_qza("taxonomy.qza")

insertion_tree <- read_qza("rooted_tree.qza")

# Separating the headers for the greengenes taxonomy file 

taxtable<-greengenes_taxonomy$data %>% as_tibble() %>% separate(Taxon, sep=";", c("Kingdom","Phylum","Class","Order","Family","Genus","Species")) 


# Creating a phyloseq object  

physeq<-phyloseq(otu_table(ASVs$data, taxa_are_rows = T), phy_tree(insertion_tree$data), tax_table(as.data.frame(taxtable) %>% select(-Confidence) %>% column_to_rownames("Feature.ID") %>% as.matrix()), sample_data(metadata %>% as.data.frame() %>% column_to_rownames("sampleid")))


# Creating a dataframe of the phyloseq object 

Peter_V_EBC1_df <- as.data.frame(sample_data(physeq))


# Now to look at the library size 

Peter_V_EBC1_df$LibrarySize<-sample_sums(physeq)


# We want to visualise the library size to make sure most of the EBCs have a lower seq count (just do EBC1's first, do core metrics. Do FCControls only from the 
table you created after you did the EBC1's and run core metrics)

Peter_V_EBC1_df <- Peter_V_EBC1_df[order(Peter_V_EBC1_df$LibrarySize),]
Peter_V_EBC1_df$Index <- seq(nrow(Peter_V_EBC1_df))
ggplot(data=Peter_V_EBC1_df, aes(x=Index, y=LibrarySize, color=EBC1vsFCControlvsBiological)) +geom_point()


# Now to identify the prevalence of contaminants 
# After investigating the prevalence plot, we identified a threshold of 0.58 to be best for this dataset 

sample_data(physeq)$is.neg <- sample_data(physeq)$EBC1vsFCControlvsBiological == "EBC1"

Peter_V_EBC1_contamdf.prev <- isContaminant(physeq, method="prevalence", neg="is.neg", threshold=0.58)

ggplot(data = Peter_V_EBC1_contamdf.prev, aes(x=p)) +
  geom_histogram(binwidth = 0.01) +
  labs(x = 'decontam Score', y='Number of species')

#write frequency table to file 

write.table(Peter_V_EBC1_contamdf.prev, file = "/mnt/c/Users/Kevin/QIIME_workflow/Peter_V/decontam/Peter_V_Sequencing_Analysis/frequency_table_Peter_V_ASM_EBC1.csv",
            
            sep=",", quote = FALSE, col.names=TRUE, row.names=TRUE)

# Now to look at the frequencies 

#FALSE  TRUE 
 540    42


# We want to create a plot to look at the presence/absence of features in samples and controls (EBC1's)
## This includes making p-a objects of the controls (EBC1's) and samples and a new dataframe of prevalence for all samples 

physeq.pa <- transform_sample_counts(physeq, function(abund) 1*(abund>0))

physeq.pa.controls <- prune_samples(sample_data(physeq.pa)$EBC1vsFCControlvsBiological == "EBC1", physeq.pa)

physeq.pa.samples <- prune_samples(sample_data(physeq.pa)$EBC1vsFCControlvsBiological == "Biological", physeq.pa)

Peter_V_EBC1_df.pa <- data.frame(pa.pos=taxa_sums(physeq.pa.samples), pa.neg=taxa_sums(physeq.pa.controls), contaminant=Peter_V_EBC1_contamdf.prev$contaminant)

ggplot(data=Peter_V_EBC1_df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() + xlab("Prevalence (Negative Controls)") + ylab("Prevalence (Samples)") 

# Now we want to remove the contaminant ASVs listed below 

5f42f7e9e0d71fb8eb76daf377e8bdaa
d6da8ab6a0c7f969da0b1c93a085cda9
402bd5ac1951d19fbd3f312c7e5f1548
8dd93b6bc59f35b38d75d22f5cf66b23
bc4ae61d84daab04b0e888d22f7d3f69
d4d20c8a4b0b24df88825a93233ab689
a1999b4f10c91d6f8efa3c9d73f77fa4
6ebf4f21e7698424e1e3683df98e46fa
f2ad4d1da17d195fabd630ac869642cb
64487332c86a1167bd17614bf16f0b04
7a8d29c59b803baaed9cc1f04ce0dc33
51b77aac70949ef74778c977216dfa73
d56536a97c18adccf18eb08958d93345
d284b92d7518959fe5d52aa45b687394
f1da0f084f3e5b7957d11d01abe7013d
3f6c45001fd40e8b6f9d0558539fd820
8e63cbca108332aa618f7a75160529b3
a06f5fd3db9f925457a537135871f616
fad083926a3ee4a61302bf4336a15a53
e3a4a62adc4e83c274978806164b474a
a3fd29484f59d9bb1389beeef0eb80c7
f65114ea452e68521c4d87e46119613b
c75334a3dac048a23b4020109b8d0595
08a1598eb1ea4dca23db1135552cd5a7
0e832b78c3ed7217d516ae2f40426bef
54f25a2b09b2c1a1517c8b5472eca178
fefd20753a9f85375759b4750fc21d20
0d4b343802645f544e320bde8567ec75
ea28f843f0962949a7b82a70b75476fd
9f74473636e6874853d88dfcbcecff49
c5adecb440a6bac201ab573ccd2d6dda
1a0bea5ccf9d532d75136168ddb54654
47ad35356a9bfec68416d32e4f039021
16a283d6a5cd2c328067ce230b8f316b
1efc4603f7ff0b5563a445542d5ef209
e834390460ccf697f792be45c813a701
33b601a7a8295e7a6c8f413504129220
8ae84ad7043388535d2e24cb63c57526
59c4a678634e50b3fda861310820aadd
1334c4fecc102315458e181593b4a34d
34a957d43e36c8b937c6bfe54f010bb1
3df77157a242066dfb4eb9c19062436a

# As there are issues with biom format and the latest version of R, we will remove these manually 
## To do this we will create a .txt file with the column names "FeatureID" and "Frequency"
## Insert the contaminant ASVs under the "FeatureID" header
## Now we save that and move back into QIIME2 to remove these contaminants

# Assigning taxa 
# Import qiime taxonomy file (.qza)

library(qiime2R)
taxonomy <- read_qza("taxonomy.qza")
taxonomy<-parse_taxonomy(taxonomy$data)

# Convert row names to columns in dataframes, so they can be merged

library(data.table)
taxonomy <- tibble::rownames_to_column(taxonomy, "featureID")
Peter_V_EBC1_contamdf.prev <- tibble::rownames_to_column(Peter_V_EBC1_contamdf.prev, "featureID")

# Join tables to find taxonomy of contaminants
contaminants_taxonomy <- dplyr::semi_join(taxonomy, Peter_V_EBC1_contamdf.prev, by = "featureID")

# Now you can view the contaminant taxa and save this as a file to open in Excel if you want.

write.table(contaminants_taxonomy, file = "/mnt/c/Users/Kevin/QIIME_workflow/Peter_V/decontam/Peter_V_Sequencing_Analysis/frequency_table_Peter_V_ASM_EBC1_Taxa.csv",
            
            sep=",", quote = FALSE, col.names=TRUE, row.names=TRUE)

#Remove False feature IDs using excel
#Copy pasted the table containing all the assigned taxa (for both contaminants and non-contaminants)
#to my contaminants excel file, went to home > conditional formatting > highlight cell rules > 
#duplicates values > ok (which coloured all the contaminant feature ids). After that, I custom sorted
#based on the feature id column and the cell colour
#file saved as 'ContaminantsToRemove_ASM_EBC1.txt'

# We have to remove the contaminants from the .txt file created 

qiime feature-table filter-features \
--i-table Peter_V_ASM_table.qza \
--p-exclude-ids \
--m-metadata-file ContaminantsToRemove_ASM_EBC1.txt \
--o-filtered-table decontam_table_ASM_EBC1.qza

qiime feature-table summarize \
--i-table decontam_table_ASM_EBC1.qza \
--o-visualization decontam_table_ASM_EBC1.qzv \
--m-sample-metadata-file Peter_V_merged_map.txt 

############## REMOVING SINGELTONS AND SAMPLES WITH LOW READS ############ 

# Removing singletons that have been created from removing contaminants before samples with low reads 

qiime feature-table filter-features \
--i-table decontam_table_ASM_EBC1.qza \
--p-min-frequency 2 \
--o-filtered-table final_Peter_V_table_ASM_EBC1.qza

qiime feature-table summarize \
--i-table final_Peter_V_table_ASM_EBC1.qza \
--o-visualization final_Peter_V_table_ASM_EBC1.qzv \
--m-sample-metadata-file Peter_V_merged_map.txt  

# Now that we have our final table, we can perform diversity analyses, compositional analyses, and statistical analyses 

# Look at features in final_Peter_V_table_ASM_EBC1.qzv to see where to subsample to maximize sequences but limit number of samples excluded
# Look at the unweighted unifrac emperor plot and compare it with the same plot before decontam filteration (contaminant removal).
#This will show you the impact of decontam filtration
#Sampling depth 740 

qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted_tree.qza \
  --i-table final_Peter_V_table_ASM_EBC1.qza \
  --p-sampling-depth 740 \
  --m-metadata-file Peter_V_merged_map.txt \
  --output-dir Peter_V_ASM_EBC1_CoreDiversity_Results

# Run stats to see if metadata explains any beta diversity characteristics
# This is for beta diversity 

  qiime diversity beta-group-significance \
  --i-distance-matrix Peter_V_ASM_EBC1_CoreDiversity_Results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-column EBC1vsFCControlvsBiological \
  --m-metadata-file Peter_V_merged_map.txt \
  --o-visualization Peter_V_ASM_EBC1_CoreDiversity_Results/unweighted_unifrac_significance_EBC1vsFCControlvsBiological.qzv

# Only run beta diversity for ControlvsBiological column if you need to 

  qiime diversity beta-group-significance \
  --i-distance-matrix Peter_V_ASM_EBC1_CoreDiversity_Results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-column ControlvsBiological \
  --m-metadata-file Peter_V_merged_map.txt \
  --o-visualization Peter_V_ASM_EBC1_CoreDiversity_Results/unweighted_unifrac_significance_ControlvsBiological.qzv

# ALPHA RAREFACTION
# Look at alpha rarefaction of your data

qiime diversity alpha-rarefaction \
  --i-table decontam_table_ASM_EBC1.qza \
  --i-phylogeny rooted_tree.qza \
  --p-max-depth 740 \
  --m-metadata-file Peter_V_merged_map.txt \
  --o-visualization Peter_V_ASM_EBC1_CoreDiversity_Results/alpha_rarefaction.qzv

# Run stats to see if metadata explains any alpha or beta diversity characteristics
# This is for alpha diversity

   qiime diversity alpha-group-significance \
  --i-alpha-diversity Peter_V_ASM_EBC1_CoreDiversity_Results/observed_features_vector.qza \
  --m-metadata-file Peter_V_merged_map.txt \
  --o-visualization Peter_V_ASM_EBC1_CoreDiversity_Results/observed_features_significance.qzv

######## Decontam for ASM FCControl ###########

# Installing packages 

library(tidyverse)
library(readr)
library(phyloseq)
library(ggplot2)
library(decontam)
library(scales)
library(qiime2R)

# Importing QIIME2 files for Decontam 

metadata<-read_tsv("Peter_V_merged_map.txt")

ASVs <- read_qza("decontam_table_ASM_EBC1.qza")

greengenes_taxonomy <- read_qza("taxonomy.qza")

insertion_tree <- read_qza("rooted_tree.qza")

# Separating the headers for the greengenes taxonomy file 

taxtable<-greengenes_taxonomy$data %>% as_tibble() %>% separate(Taxon, sep=";", c("Kingdom","Phylum","Class","Order","Family","Genus","Species")) 


# Creating a phyloseq object  

physeq<-phyloseq(otu_table(ASVs$data, taxa_are_rows = T), phy_tree(insertion_tree$data), tax_table(as.data.frame(taxtable) %>% select(-Confidence) %>% column_to_rownames("Feature.ID") %>% as.matrix()), sample_data(metadata %>% as.data.frame() %>% column_to_rownames("sampleid")))


# Creating a dataframe of the phyloseq object 

Peter_V_FCControl_df <- as.data.frame(sample_data(physeq))


# Now to look at the library size 

Peter_V_FCControl_df$LibrarySize<-sample_sums(physeq)


# We want to visualise the library size to make sure most of the EBCs have a lower seq count

Peter_V_FCControl_df <- Peter_V_FCControl_df[order(Peter_V_FCControl_df$LibrarySize),]
Peter_V_FCControl_df$Index <- seq(nrow(Peter_V_FCControl_df))
ggplot(data=Peter_V_FCControl_df, aes(x=Index, y=LibrarySize, color=EBC1vsFCControlvsBiological)) +geom_point()


# Now to identify the prevalence of contaminants 
# After investigating the prevalence plot, we identified a threshold of 0.58 to be best for this dataset 

sample_data(physeq)$is.neg <- sample_data(physeq)$EBC1vsFCControlvsBiological == "FCControl"

Peter_V_FCControl_contamdf.prev <- isContaminant(physeq, method="prevalence", neg="is.neg", threshold=0.58)

ggplot(data = Peter_V_FCControl_contamdf.prev, aes(x=p)) +
  geom_histogram(binwidth = 0.01) +
  labs(x = 'decontam Score', y='Number of species')

#write frequency table to file 

write.table(Peter_V_FCControl_contamdf.prev, file = "/mnt/c/Users/Kevin/QIIME_workflow/Peter_V/decontam/Peter_V_Sequencing_Analysis/frequency_table_Peter_V_ASM_FCControl.csv",
            
            sep=",", quote = FALSE, col.names=TRUE, row.names=TRUE)

# Now to look at the frequencies 

#FALSE  TRUE 
505    35


# We want to create a plot to look at the presence/absence of features in samples and controls (FCControl)
## This includes making p-a objects of the controls (FCControl) and samples and a new dataframe of prevalence for all samples 

physeq.pa <- transform_sample_counts(physeq, function(abund) 1*(abund>0))

physeq.pa.controls <- prune_samples(sample_data(physeq.pa)$EBC1vsFCControlvsBiological == "FCControl", physeq.pa)

physeq.pa.samples <- prune_samples(sample_data(physeq.pa)$EBC1vsFCControlvsBiological == "Biological", physeq.pa)

Peter_V_FCControl_df.pa <- data.frame(pa.pos=taxa_sums(physeq.pa.samples), pa.neg=taxa_sums(physeq.pa.controls), contaminant=Peter_V_FCControl_contamdf.prev$contaminant)

ggplot(data=Peter_V_FCControl_df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() + xlab("Prevalence (Negative Controls)") + ylab("Prevalence (Samples)") 

# Now we want to remove the contaminant ASVs listed below 

a148bd855fdf35b23b20751aeb334ab8
93870814f567e00303722ba3fa36d6ea
1536bfe63a374c76fa0a4de8b927f947
053c4d5b8717b60387f777efae437b40
e66860dbaac1010e5b69f9780d84b3b2
68bb0065fa02ba7b52acf3f4d833895b
6e8986f8088b452f964e7968bb8bca87
caf54b4995a78b51266cb0e4aed4deae
2ea79eaaa0b23b2a7d2cf470cbaf72ea
a1dabd1d0db4a5aacc72a9e57e5c2f0d
088764d8589b20c4461df3ce11980d03
a930b17afabc7d14691776ff8b71efc2
fdcd6808ef8269653d25dce4a55a025d
505e1db623887f714a204ca9f6c7e33e
0374d69a06fe46ad16ded615ae85090d
8951cc9103140b71bd50c0932e096502
169f191613e3cc1f6a12b99c72948aca
aa5e660e6ac75fcc7aa5182aa2979650
73d8c2a8aca31125af37a01b15659445
614c1fdeb84fa86567a6cd31afb6ec1b
60e4d1ebe5228cf27ec0881e1d9f1933
59b93443b5560fef2e51483de0a22e4d
8505d535cd85cd48fbbe158406992f60
2d5bd471a9291e60c0e070a25bfa2074
38e3b7abf4e08699d7a1d1965f177086
3cd11adc295e84455609b45db54fed9a
ba96c9f5431c4bbaf2b9d348435c7f4e
8b1b161673921734fede6930a726fd71
3e9448d3db2cec00781e9b42ed85f8b1
0ceea5e52cad97a12c3a98d094754cd1
7db357a4adb452de8504d8ca1ddb3f99
38182ed483c50ea819ae212f3833d773
ffd7393b54add92dafea3368794ccb2c
4431e66562ad16df2201d5148364269f
465c45bc4b0ed05e986a5c3737cbc36c

# As there are issues with biom format and the latest version of R, we will remove these manually 
## To do this we will create a .txt file with the column names "FeatureID" and "Frequency"
## Insert the contaminant ASVs under the "FeatureID" header
## Now we save that and move back into QIIME2 to remove these contaminants

#Assigning taxa to FCControls (both true and false)

# Import qiime taxonomy file (.qza)
library(qiime2R)
taxonomy <- read_qza("taxonomy.qza")
taxonomy<-parse_taxonomy(taxonomy$data)

# Convert row names to columns in dataframes, so they can be merged
library(data.table)
taxonomy <- tibble::rownames_to_column(taxonomy, "featureID")
Peter_V_FCControl_contamdf.prev <- tibble::rownames_to_column(Peter_V_FCControl_contamdf.prev, "featureID")

# Join tables to find taxonomy of contaminants
contaminants_taxonomy <- dplyr::semi_join(taxonomy, Peter_V_FCControl_contamdf.prev, by = "featureID")

# Now you can view the contaminant taxa and save this as a file to open in Excel if you want.

write.table(contaminants_taxonomy, file = "/mnt/c/Users/Kevin/QIIME_workflow/Peter_V/decontam/Peter_V_Sequencing_Analysis/frequency_table_FCControl_AllTaxa.csv",
            
            sep=",", quote = FALSE, col.names=TRUE, row.names=TRUE)

#Remove False feature IDs using excel
#Copy pasted the table containing all the assigned taxa (for both contaminants and non-contaminants)
#to my contaminants excel file, went to home > conditional formatting > highlight cell rules > 
#duplicates values > ok (which coloured all the contaminant feature ids). After that, I custom sorted
#based on the feature id column and the cell colour
#file saved as 'ContaminantsToRemove_ASM_FCControl.txt'

# We have to remove the contaminants from the .txt file created 

qiime feature-table filter-features \
--i-table decontam_table_ASM_EBC1.qza \
--p-exclude-ids \
--m-metadata-file ContaminantsToRemove_ASM_FCControl.txt \
--o-filtered-table decontam_table_ASM_FCControl.qza

qiime feature-table summarize \
--i-table decontam_table_ASM_FCControl.qza \
--o-visualization decontam_table_ASM_FCControl.qzv \
--m-sample-metadata-file Peter_V_merged_map.txt 

############## REMOVING SINGELTONS AND SAMPLES WITH LOW READS ############ 

# Removing singletons that have been created from removing contaminants before samples with low reads 

qiime feature-table filter-features \
--i-table decontam_table_ASM_FCControl.qza \
--p-min-frequency 2 \
--o-filtered-table final_Peter_V_table_ASM_FCControl.qza

qiime feature-table summarize \
--i-table final_Peter_V_table_ASM_FCControl.qza \
--o-visualization final_Peter_V_table_ASM_FCControl.qzv \
--m-sample-metadata-file Peter_V_merged_map.txt  

# Now that we have our final table, we can perform diversity analyses, compositional analyses, and statistical analyses 

# Look at features in final_Peter_V_table_ASM_FCControl.qzv to see where to subsample to maximize sequences but limit number of samples excluded
# Look at the unweighted unifrac emperor plot and compare it with the same plot before decontam filteration (contaminant removal).
#This will show you the impact of decontam filtration
#Sampling depth 827

qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted_tree.qza \
  --i-table final_Peter_V_table_ASM_FCControl.qza \
  --p-sampling-depth 827 \
  --m-metadata-file Peter_V_merged_map.txt \
  --output-dir Peter_V_ASM_FCControl_CoreDiversity_Results

# Run stats to see if metadata explains any beta diversity characteristics
# This is for beta diversity 

  qiime diversity beta-group-significance \
  --i-distance-matrix Peter_V_ASM_FCControl_CoreDiversity_Results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-column EBC1vsFCControlvsBiological \
  --m-metadata-file Peter_V_merged_map.txt \
  --o-visualization Peter_V_ASM_FCControl_CoreDiversity_Results/unweighted_unifrac_significance_EBC1vsFCControlvsBiological.qzv

# Only run beta diversity for ControlvsBiological column if you need to 

  qiime diversity beta-group-significance \
  --i-distance-matrix Peter_V_ASM_FCControl_CoreDiversity_Results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-column ControlvsBiological \
  --m-metadata-file Peter_V_merged_map.txt \
  --o-visualization Peter_V_ASM_FCControl_CoreDiversity_Results/unweighted_unifrac_significance_ControlvsBiological.qzv

  #New metadata column for controls made 

   qiime diversity beta-group-significance \
  --i-distance-matrix Peter_V_ASM_FCControl_CoreDiversity_Results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-column New_ControlvsBiological \
  --m-metadata-file Peter_V_merged_map.txt \
  --o-visualization Peter_V_ASM_FCControl_CoreDiversity_Results/unweighted_unifrac_significance_New_ControlvsBiological.qzv


# ALPHA RAREFACTION
# Look at alpha rarefaction of your data

qiime diversity alpha-rarefaction \
  --i-table decontam_table_ASM_FCControl.qza \
  --i-phylogeny rooted_tree.qza \
  --p-max-depth 827 \
  --m-metadata-file Peter_V_merged_map.txt \
  --o-visualization Peter_V_ASM_FCControl_CoreDiversity_Results/alpha_rarefaction_FCControl.qzv

# Run stats to see if metadata explains any alpha or beta diversity characteristics
# This is for alpha diversity

   qiime diversity alpha-group-significance \
  --i-alpha-diversity Peter_V_ASM_FCControl_CoreDiversity_Results/observed_features_vector.qza \
  --m-metadata-file Peter_V_merged_map.txt \
  --o-visualization Peter_V_ASM_FCControl_CoreDiversity_Results/observed_features_significance.qzv

#Filtered samples on your metadata column which is "ControlvsBiological_ASM_General"

qiime feature-table filter-samples \
--i-table decontam_table_ASM_FCControl.qza \
--m-metadata-file Peter_V_merged_map.txt \
--p-where "ControlvsBiological_ASM_General='na'" \
--p-exclude-ids \
--o-filtered-table Peter_V__ControlvsBiological_ASM_General_table.qza

#Ran core metrics (you ran for the ASM_General column but there is not much difference between
#this column and ControlvsBiological_ASM)

qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted_tree.qza \
  --i-table Peter_V_ControlvsBiological_ASM_General_table.qza \
  --p-sampling-depth 827 \
  --m-metadata-file Peter_V_merged_map.txt \
  --output-dir Peter_V_ControlvsBiological_ASM_General_CoreDiversity_Results

#Ran beta diversity for columns "ControlvsBiological_ASM_General" 

qiime diversity beta-group-significance \
  --i-distance-matrix Peter_V_ControlvsBiological_ASM_General_CoreDiversity_Results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-column ControlvsBiological_ASM_General \
  --m-metadata-file Peter_V_merged_map.txt \
  --o-visualization Peter_V_ControlvsBiological_ASM_General_CoreDiversity_Results/unweighted_unifrac_significance_ControlvsBiological_ASM_General.qzv

# Filtering according to new metadata column "ASM_T0vsT14" to give ASM samples that are T0 (Collected and Innoculated) and T14
#This is just an overall column and DOES NOT go in to individual level (Individual level column is T0vsT14_ASM)
 
qiime feature-table filter-samples \
--i-table decontam_table_ASM_FCControl.qza \
--m-metadata-file Peter_V_merged_map.txt \
--p-where "ASM_T0vsT14='na'" \
--p-exclude-ids \
--o-filtered-table Peter_V_ASM_T0vsT14_table.qza

qiime feature-table summarize \
--i-table Peter_V_ASM_T0vsT14_table.qza \
--o-visualization Peter_V_ASM_T0vsT14_table.qzv \
--m-sample-metadata-file Peter_V_merged_map.txt  

# Running core-metrics. This folder contains only biological samples (No EBC1's and FCControls)
#Sampling depth 26608

qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted_tree.qza \
  --i-table Peter_V_ASM_T0vsT14_table.qza \
  --p-sampling-depth 26608 \
  --m-metadata-file Peter_V_merged_map.txt \
  --output-dir Peter_V_ASM_T0vsT14_CoreDiversity_Results

# Beta group significance (general level)

qiime diversity beta-group-significance \
  --i-distance-matrix Peter_V_ASM_T0vsT14_CoreDiversity_Results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-column ASM_T0vsT14 \
  --m-metadata-file Peter_V_merged_map.txt \
  --o-visualization Peter_V_ASM_T0vsT14_CoreDiversity_Results/unweighted_unifrac_significance_ASM_T0vsT14.qzv

#Individual level

qiime diversity beta-group-significance \
  --i-distance-matrix Peter_V_ASM_T0vsT14_CoreDiversity_Results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-column T0vsT14_ASM \
  --m-metadata-file Peter_V_merged_map.txt \
  --o-visualization Peter_V_ASM_T0vsT14_CoreDiversity_Results/unweighted_unifrac_significance_T0vsT14_ASM.qzv



# Alpha group significance 

qiime diversity alpha-group-significance \
  --i-alpha-diversity Peter_V_ASM_T0vsT14_CoreDiversity_Results/observed_features_vector.qza \
  --m-metadata-file Peter_V_merged_map.txt \
  --o-visualization Peter_V_ASM_T0vsT14_CoreDiversity_Results/observed_features_significance.qzv

#Alpha rarefaction

  qiime diversity alpha-rarefaction \
  --i-table Peter_V_ASM_T0vsT14_table.qza \
  --i-phylogeny rooted_tree.qza \
  --p-max-depth 26608 \
  --m-metadata-file Peter_V_merged_map.txt \
  --o-visualization Peter_V_ASM_T0vsT14_CoreDiversity_Results/alpha_rarefaction_ASM_T0vsT14.qzv

# Classify features and visualisation

qiime feature-classifier classify-sklearn \
  --i-classifier gg-13-8-99-515-806-nb-classifier.qza \
  --i-reads merged_repseqs.qza \
  --o-classification taxonomy.qza

# Create visualization

qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy.qzv

# Creating the plots to visualise

  qiime taxa barplot \
  --i-table Peter_V_ASM_T0vsT14_table.qza \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file Peter_V_merged_map.txt \
  --o-visualization Peter_V_ASM_T0vsT14_CoreDiversity_Results/ASM_T0vsT14_Taxonomy_BarPlots.qzv

#Preparing QIIME2 files for LEfSe

#ASM 
#Collapse table.qza to level 7

qiime taxa collapse \
  --i-table Peter_V_ASM_T0vsT14_table.qza \
  --i-taxonomy taxonomy.qza \
  --p-level 7 \
  --o-collapsed-table Peter-V-ASM-T0vsT14-table-L7.qza

#Calculate relative-frequency for the collapsed table (instead of counts you get relative abundance)

qiime feature-table relative-frequency \
--i-table Peter-V-ASM-T0vsT14-table-L7.qza \
--o-relative-frequency-table Peter-V-ASM-T0vsT14-frequency-table-L7.qza \
--output-dir Peter-V-ASM-T0vsT14-frequency/

#Export biom file

qiime tools export \
--input-path Peter-V-ASM-T0vsT14-frequency-table-L7.qza \
--output-path Peter-V-ASM-T0vsT14-frequency-L7/

#Convert biom to text file (for lefse comparison) 
#Had to make a copy of this file in Peter_v_Sequencing_Analysis folder from Peter-V-ASM-T0vsT14-frequency

biom convert \
-i Peter-V-ASM-T0vsT14-frequency-table-L7.biom \
-o Peter-V-ASM-T0vsT14-frequency-table-L7.txt \
--header-key “taxonomy” --to-tsv

#Running the above commands at for samples at an individual level 

#Collapse table.qza to level 7

qiime taxa collapse \
  --i-table Peter_V_T0vsT14_ASM_table.qza \
  --i-taxonomy taxonomy.qza \
  --p-level 7 \
  --o-collapsed-table Peter-V-T0vsT14-ASM-table-L7.qza

#Calculate relative-frequency for the collapsed table (instead of counts you get relative abundance)

qiime feature-table relative-frequency \
--i-table Peter-V-T0vsT14-ASM-table-L7.qza \
--o-relative-frequency-table Peter-V-T0vsT14-ASM-frequency-table-L7.qza \
--output-dir Peter-V-T0vsT14-ASM-frequency-L7/

#Export biom file

qiime tools export \
--input-path Peter-V-T0vsT14-ASM-frequency-table-L7.qza \
--output-path Peter-V-T0vsT14-ASM-frequency-L7/

#Convert biom to text file (for lefse comparison) 
#Had to make a copy of this file in Peter_v_Sequencing_Analysis folder from Peter-V-ASM-T0vsT14-frequency

biom convert \
-i Peter-V-T0vsT14-ASM-feature-table.biom \
-o Peter-V-T0vsT14-ASM-frequency-table-L7.txt \
--header-key “taxonomy” --to-tsv

#For level 2 
#Not as useful as hoped 

qiime taxa collapse \
  --i-table Peter_V_T0vsT14_ASM_table.qza \
  --i-taxonomy taxonomy.qza \
  --p-level 2 \
  --o-collapsed-table Peter-V-T0vsT14-ASM-table-L2.qza

#Calculate relative-frequency for the collapsed table (instead of counts you get relative abundance)

qiime feature-table relative-frequency \
--i-table Peter-V-T0vsT14-ASM-table-L2.qza \
--o-relative-frequency-table Peter-V-T0vsT14-ASM-frequency-table-L2.qza \
--output-dir Peter-V-T0vsT14-ASM-L2-frequency/

#Export biom file

qiime tools export \
--input-path Peter-V-T0vsT14-ASM-frequency-table-L2.qza \
--output-path Peter-V-T0vsT14-ASM-L2-frequency/

#Convert biom to text file (for lefse comparison) 
#Had to make a copy of this file in Peter_v_Sequencing_Analysis folder from Peter-V-ASM-T0vsT14-frequency

biom convert \
-i Peter-V-T0vsT14-ASM-L2-feature-table.biom \
-o Peter-V-T0vsT14-ASM-frequency-table-L2.txt \
--header-key “taxonomy” --to-tsv

#! /bin/bash

#Created a heatmap based based on individual patients over time 
#First done for ASM
#QIIMEView does not let you download the QZV file. So export the file to view it 
#Clustered on 'features'

qiime feature-table filter-samples \
--i-table Peter_V_ASM_T0vsT14_table.qza \
--m-metadata-file Peter_V_merged_map.txt \
--p-where "T0vsT14_ASM='na'" \
--p-exclude-ids \
--o-filtered-table Peter_V_T0vsT14_ASM_table.qza

qiime taxa collapse \
  --i-table Peter_V_T0vsT14_ASM_table.qza \
  --i-taxonomy taxonomy.qza \
  --p-level 7 \
  --o-collapsed-table Peter-V-T0vsT14-ASM-Individual-Level-table-L7.qza

qiime feature-table heatmap \
--i-table Peter-V-T0vsT14-ASM-Individual-Level-table-L7.qza \
--m-sample-metadata-file Peter_V_merged_map.txt \
--m-sample-metadata-column T0vsT14_ASM \
--p-cluster 'features' \
--o-visualization Peter-V-T0vsT14-ASM-heatmap-L7.qzv

qiime tools export \
  --input-path Peter-V-T0vsT14-ASM-heatmap-L7.qzv \
  --output-path exported-Peter-V-T0vsT14-ASM-heatmap-L7

#Clustered on 'samples'

qiime feature-table heatmap \
--i-table Peter-V-T0vsT14-ASM-Individual-Level-table-L7.qza \
--m-sample-metadata-file Peter_V_merged_map.txt \
--m-sample-metadata-column T0vsT14_ASM \
--p-cluster 'samples' \
--o-visualization Peter-V-T0vsT14-ASM-heatmap-samples-L7.qzv

qiime tools export \
  --input-path Peter-V-T0vsT14-ASM-heatmap-samples-L7.qzv \
  --output-path exported-Peter-V-T0vsT14-ASM-heatmap-samples-L7











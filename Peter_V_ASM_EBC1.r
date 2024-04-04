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
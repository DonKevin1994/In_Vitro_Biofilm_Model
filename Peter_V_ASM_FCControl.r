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

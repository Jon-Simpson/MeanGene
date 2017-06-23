# MeanGene
MeanGene is an unsupervised method for identifying groups of genes that are expressed differently between RNA-sequenced samples. The name MeanGene was chosen because K-Means is the engine behind the method. The idea was to find groups of meaningful genes (i.e. genes that all contribute to a certain function) that were expressed differently between groups of samples.

# Data
# 1) Samples
The initial data I explored was the gene counts from 9 RNA-sequenced dolphins from Charleston, SC.
Blood samples were collected from the dolphins in August 2013 in Charleston during the Dolphin Risk and Health Assessment.
The samples were taken to a lab where they were cultured. They were separated into groups and underwent in vitro exposure to varying levels of PFOS (perfluorooctane sulfonic acid). The purpose of the experiment was to study the immune response of the cells when exposed. PFOS is a highly stable compound found abundantly in the environment that is known to be immunotoxic. 9 samples were then sequenced for RNA. 

Each sequenced file, after genome alignment, consisisted of over 100,000 genes, however only around ~20,000 or so were significantly expressed.

# 2) Metadata
The metadata for the samples included sex, age (adult or juvenile), and a low or high label, depending on the experiental exposure to PFOS.



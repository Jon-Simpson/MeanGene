# MeanGene
MeanGene is a class written in Python with the purpose of implementing an unsupervised method for identifying patterns of gene expression given a DEseq output. DESeq is an R package for differential expression analysis of count data in high-throughput sequencing. More information on DESeq can be found [here](http://bioconductor.org/packages/release/bioc/html/DESeq.html). The name MeanGene was chosen because the K-Means algorithm is the engine behind the method. The aim was to find groups of associated genes (i.e. genes relating to cancer, etc.) that were expressed differently between groups of samples. 

# Manual
## Background
MeanGene was designed for clustering genetic expression data, and more specifically, count data that comes from a DESeq output. The reason for this is that the method was made for the use of researchers in the Bioinformatics Center for Genomic Medicine at the Medical University of South Carolina. However, only the initialization of meangene.py is specific to DESeq output, and with minor modifications, could be used for other inputs providing gene symbols and counts. Furthermore, The class could be used for other purposes other than genetic expression clustering. Any high dimensional data, where subsets of the dimensions are used as features for clustering is valid use of the tool. 

To further detail the aim of MeanGene, I'll explain the 3 types of data needed to optimize results:
## data
- **Count Data**: A dataframe of the count data for genes (represented as symbols), for as many samples as you may have
- **Metadata**: Such as phenotypic data. (for example; I knew the sex, age, and contamination level of a certain chemical for my samples). The more metadata the better, since the metadata allows for more opportunity to find correlations between patterns found and the phenotypes.
- **Subsets of genes**: A list of subsets of genes. This where MeanGene differentiates itself from other implementations of genetic pattern recognition. You supply MeanGene with many, many subsets of genes, and ideally the sets will be related to some biological function (i.e. breast cancer, lupus, influenza, etc). If a pattern is found among a set of genes, you have just discovered a function or disease that is differentially expressed between your samples. More on retrieving this information in a bit.

## Implementation
First initialize your MeanGene object with
```
meangene(DESeqDF, metadata)
```
Where DESeqDF is a dataframe of a DESeq output. 




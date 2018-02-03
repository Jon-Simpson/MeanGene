# MeanGene
MeanGene is a class written in Python that utilizes a combination of supervised and unsupervised learning methods for identifying patterns of gene expression given a DEseq output. DESeq is an R package for differential expression analysis of count data in high-throughput sequencing. More information on DESeq can be found [here](http://bioconductor.org/packages/release/bioc/html/DESeq.html). The name MeanGene is a play off the "K-Means" algorithm, as it is the engine behind the method. The aim of the class is to identify groups of genes that can explain the difference in expression among samples. 

# Manual
## Background
MeanGene was designed to cluster genetic expression count data coming from a DESeq output. Originally, the method was created for use by the researchers in the Bioinformatics Center for Genomic Medicine at the Medical University of South Carolina. However, only the initialization of meangene.py is specific to DESeq output, and with minor modifications, could be used for other inputs providing gene symbols and counts. Furthermore, The class could be used for other purposes other than genetic expression clustering. Any unsupervised problem with high dimensional data, where subsets of the dimensions are used as features for clustering, is valid use of the tool. 

To further detail the aim of MeanGene, I'll explain the 3 types of data needed to optimize results:
## data
- **Count Data**: A dataframe of the count data for each gene/sample. The dataframe should created from be a DESeq output, **trimmed of all columns excluding the count data**.
- **Metadata**: Such as phenotypic data. (for example; I knew the sex, age, and contamination level of a certain chemical for my samples). The more metadata the better, since the metadata allows for more opportunity to find correlations between patterns found and the phenotypes.
- **Subsets of genes**: A list of subsets of genes. This where MeanGene differentiates itself from other implementations of genetic pattern recognition. You supply MeanGene with many, many subsets of genes, and ideally the sets will be related to some biological function (i.e. breast cancer, lupus, influenza, etc). If a pattern is found among a set of genes, you have just discovered a function or disease that is differentially expressed between your samples. More on retrieving this information in a bit.

## Implementation
First initialize your MeanGene object with
```
meangene(DESeqDF, metadata)
```
Where DESeqDF is a dataframe of a DESeq output. 

## Methods
There are two methods ```runPCA(geneFunctions, components=2, subsets=None)``` and ```cluster(transformation=None)```.

```runPCA()``` performs principal component analysis on your features. Where 'geneFunctions' is an list of strings that describe each subset of genes given, 'components' is the number of principal components to actually use for the clustering algorithm (set to 2 by default), and 'subsets' is a 2-dimensional array of sets of genes.

```cluster()``` performs K-Means clustering on the subsets, using their principal components, if ```transformation=None```, then clustering is performed on the entire set of genes given from the DESeq (whole transcriptome).

## Analysis
After ```cluster()``` has been run, the meangene object will now have an attribute that highlights the analysis, which can be accessed using ```meangene.analysisDF```. The returned dataframe depicts the algorithm's results for each set of genes. There will be columns for the analysisID index number, the disease or function of the set of genes, the PCA object, an embedded dataframe containing grid coordinates for plotting, and a subsequent column for every metadata category that provides a number from 0 to 1 describing how well the data separated that category. A 1 indicates for a set of genes, that category was completely separated using those genes, meaning there is differential expression between the two groups. For example, if a 1 is given for a set of 45 genes relating to arthritis in the sex category, that means no cluster contained both males and females, suggesting that genes contributing the arthritis are likely expression differently between males and females. Whereas a 0 indicates perfect overlap between groups, where all types in a single category are represented equally in all clusters. 







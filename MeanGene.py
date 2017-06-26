import pandas as pd
import numpy as np
from sklearn import preprocessing as pp
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans


class meangene(object):
    def __init__(self, ExpressionFile, MetadataFile):

        # File1 is a DEseq expression (txt file) from OnRamp Bioinformatics pipeline
        # File2 is a metadata txt file
        # NOTE: Samples must be in the same order for both files & there can't be missing values in the metadata file

        self.data = ExpressionFile
        self.metadataDF = MetadataFile

        # Remove the genes with zero expressions across all samples
        zero_columns = self.data.mean() == 0
        self.data = self.data.drop(self.data.columns[zero_columns], axis=1)
        self.data.columns = map(str.upper, self.data.columns)

        # Generate sample IDs (Starting with 0 and going to n-1) and input as 1st column on metadat
        self.metadataDF.insert(0, 'SampleID', np.arange(len(self.metadataDF)))

        # Save data as attributes for accessibility
        self.data_head = self.data.head()

    def runPCA(self, geneFunctions, components=2, subsets=None):


        # PCA is robust for up to 10 Principal Components
        coordNames = ['PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10']
        sampleID = self.metadataDF.SampleID

        if subsets.empty != True:

            subsets = np.array(subsets)
            geneFunctions = np.array(geneFunctions)

            # Drop the sets of genes whose length is less than the number of 'components'
            lengths = np.array([len(i) for i in subsets])
            crop = np.where(lengths >= components)
            subsets = subsets[crop]
            geneFunctions = geneFunctions[crop]

            # Generate the number of PCA analyses to do (one for each set of genes)
            analysisIDs = list(range(len(subsets)))
            coords = []
            # Start building DataFrame
            self.analysisDF = pd.DataFrame({'AnalysisID': analysisIDs})
            pcaobject = []

            for subset in subsets:
                gene_indices = np.empty(0)
                for gene in subset:
                    inds = np.where(self.data.columns == gene)
                    inds = np.array(inds).reshape(-1, 1)
                    if len(inds) > 0:
                        means = []
                        for j in inds:
                            means.append(np.mean(DESeqDF.iloc[:, j].values))
                        gene_indices = np.append(gene_indices, inds[np.argmax(means)])
                gene_indices.flatten()
                subsetdata = self.data.iloc[:, gene_indices]

                if len(subsetdata.columns) >= components:
                    # Convert table into array format for normalization of data --> Converts all values to a 0-1 scale
                    data_arr = np.array(subsetdata)
                    norm_data = pp.normalize(data_arr)

                    # Reduce the dimensions down to n principal components
                    pca = PCA(n_components=components)
                    pca.fit(norm_data)
                    pcadata = pca.transform(norm_data)

                    pcaobject.append(pca)

                    # Assemble DataFrame for each subset of genes and append to coords
                    coordnames = coordNames[:components]
                    coordlists = []
                    for i in range(components):
                        coordlists.append(list(pcadata[:, i]))
                    zipped = list(zip(coordnames, coordlists))
                    df = pd.DataFrame(dict(zipped))
                    df.insert(0, "SampleID", sampleID)
                    coords.append(df)
                else:
                    pcaobject.append(None)
                    coords.append(None)

            # Add tp final analysisDF
            self.analysisDF.insert(1, 'Genes', geneFunctions)
            self.analysisDF.insert(2, 'PCAobject', pcaobject)
            self.analysisDF.insert(3, 'Coordinates', coords)

        else:
            # Convert table into array format for normalization of data --> Converts all values to a 0-1 scale
            data_arr = np.array(self.data)
            norm_data = pp.normalize(data_arr)

            # Reduce the dimensions down to n principal components
            pca = PCA(n_components=components)
            pca.fit(norm_data)
            pcadata = pca.transform(norm_data)

            coordnames = coordNames[:components]
            coordlists = []
            for i in range(components):
                coordlists.append(list(pcadata[:, i]))
            zipped = list(zip(coordnames, coordlists))
            df = pd.DataFrame(dict(zipped))
            df.insert(0, "SampleID", sampleID)
            self.analysisDF = pd.DataFrame({'AnalysisID': [0], 'Coordinates': [df]})
            self.analysisDF.insert(1, 'Genes', 'AllGenes')
            self.analysisDF.insert(2, 'PCAobject', pca)

    def cluster(self, transformation=None):

        if (transformation == "PCA"):

            if hasattr(self, 'analysisDF') != True:
                print("Perform PCA first!")

            else:
                metadatatype = self.metadataDF.columns[2:]
                metadataseps = []

                # Run K-Means on each subset PCA and add cluster-quality metric for each label type to DataFrame
                for i in range(len(self.analysisDF)):
                    seplevels = []
                    if self.analysisDF.iloc[i, 2] != None:
                        arr = np.array(self.analysisDF.iloc[i, 3].iloc[:, 1:]).astype(float)
                        n = len(self.analysisDF.iloc[:, 3][i].iloc[:])
                        labels = 0
                        # If samples are more than 15, find the optimal number of clusters (elbow method)
                        if n > 15:
                            inertias = []
                            pct_change = []
                            optimalK = 2
                            for i in range(1, n):
                                if i > 2:
                                    if pct_change[i - 3] < 0.2:
                                        optimalK = i - 2
                                        break
                                kmeans = KMeans(n_clusters=i)
                                kmeans.fit(arr)
                                inertias.append(kmeans.inertia_)
                                if i > 1:
                                    pct_change.append((inertias[i - 2] - inertias[i - 1]) / inertias[i - 2])

                            kmeans = KMeans(n_clusters=optimalK)
                            kmeans.fit(arr)
                            labels = list(kmeans.labels_)
                        # If samples are less than 10, use 2 clusters
                        else:
                            kmeans = KMeans(n_clusters=2)
                            kmeans.fit(arr)
                            labels = list(kmeans.labels_)

                        seplevels = []
                        for i in self.metadataDF.columns[2:]:
                            indices = self.metadataDF[i].dropna().index
                            labels2 = np.array(labels)[indices]
                            meta = self.metadataDF[i][indices]
                            df = pd.DataFrame({'Labels': labels2, 'Type': meta})
                            ct = pd.crosstab(df['Labels'], df['Type'])
                            clusterpurities = []
                            for j in range(len(ct.iloc[:, 0])):
                                clusterpurity = (max(ct.iloc[j]) / ct.iloc[j].sum())
                                clusterpurities.append(clusterpurity)
                            seplevels.append(sum(clusterpurities) / len(ct.iloc[:]))
                        metadataseps.append(seplevels)
                    else:
                        metadataseps.append([None] * len(self.metadataDF.columns[2:]))
                seps = np.array(metadataseps)
                for i in range(len(metadatatype)):
                    self.analysisDF[metadatatype[i]] = seps[:, i]

        else:
            # Runs clustering on using ALL gene expressions
            if hasattr(self, 'analysisDF'):
                self.analysisDF = self.analysisDF.drop(['Genes', 'PCAobject', 'Coordinates'], axis=1)
                self.analysisDF = self.analysisDF.drop(np.arange(1, len(self.analysisDF)), axis=0)
            else:
                self.analysisDF = pd.DataFrame({'AnalysisID': [0]})
            metadatatype = self.metadataDF.columns[2:]
            arr = np.array(self.data)
            inertias = []
            pct_change = []
            optimalK = 2
            n = len(self.metadataDF)
            labels = 0
            if n > 15:
                for i in range(1, n):
                    if i > 2:
                        if pct_change[i - 3] < 0.15:
                            optimalK = i - 2
                            break
                    kmeans = KMeans(n_clusters=i)
                    kmeans.fit(arr)
                    inertias.append(kmeans.inertia_)
                    if i > 1:
                        pct_change.append((inertias[i - 2] - inertias[i - 1]) / inertias[i - 2])

                kmeans = KMeans(n_clusters=optimalK)
                kmeans.fit(arr)
                labels = list(kmeans.labels_)
            else:
                kmeans = KMeans(n_clusters=2)
                kmeans.fit(arr)
                labels = list(kmeans.labels_)
            seplevels = []
            for i in range(len(self.metadataDF.columns[2:])):
                ls = list(self.metadataDF.iloc[:, i + 2])
                df = pd.DataFrame({'Labels': labels, 'Type': ls})
                ct = pd.crosstab(df['Labels'], df['Type'])
                clusterpurities = []
                for j in range(len(ct.iloc[:, 0])):
                    clusterpurity = (max(ct.iloc[j]) / ct.iloc[j].sum())
                    clusterpurities.append(clusterpurity)
                seplevels.append(sum(clusterpurities) / len(ct.iloc[:]))
            for i in range(len(metadatatype)):
                self.analysisDF[metadatatype[i]] = seplevels[i]

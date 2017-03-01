class meangene(object):
    def __init__(self, file1, file2):

        import pandas as pd
        import numpy as np

        # File1 is a DEseq expression txt file from OnRamp
        # File2 is a metadata txt file
        # NOTE: Samples must be in the same order for both files & there can't be missing values in the metadata file
        self.file1 = file1
        self.file2 = file2

        # Import the raw text file from the DEseq file into a DataFrame
        df = pd.read_table(file1, low_memory=False)

        # Import metadata file into a DataFrame
        metaDF = pd.read_table(file2)

        # Remove the annotation and metric columns & transpose the matrix
        df = df.drop(df.columns[[1, 2, 3, -1, -2, -3, -4, -5, -6]], axis=1).transpose()

        # Make the first row (gene symbols) the column names
        header = df.iloc[0]
        df = df.iloc[1:].rename(columns=header)

        # Remove the genes with zero expressions across all samples
        zero_columns = df.mean() == 0
        df = df.drop(df.columns[zero_columns], axis=1)

        # Generate sample IDs (Starting with 0 and going to n-1) and input as 1st column on metadata
        sampleID = []
        for i in range(len(metaDF.index)):
            sampleID.append(metaDF.index[i])
        metaDF.insert(0, 'SampleID', sampleID)

        # Save data as attributes for accessibility
        self.data = df
        self.data_head = df.head()
        self.metadata = metaDF

    def runPCA(self, components=2, subsets=None):

        import numpy as np
        import pandas as pd
        from sklearn import preprocessing as pp
        from sklearn.decomposition import PCA

        # PCA is robust for up to 10 Principal Components
        coordNames = ['PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10']
        sampleID = self.metadata.SampleID

        if subsets != None:

            # Generate the number of PCA analyses to do (one for each subset of genes)
            analysisIDs = list(range(len(subsets)))
            coords = []
            # Start building DataFrame
            self.analysisDF = pd.DataFrame({'AnalysisIDs': analysisIDs})

            for subset in subsets:
                gene_indices = np.empty(0)
                for gene in subset:
                    gene_indices = np.append(gene_indices, np.where(self.data.columns == gene))

                subsetdata = self.data.iloc[:, list(gene_indices)]

                # Convert table into array format for normalization of data --> Converts all values to a 0-1 scale
                data_arr = np.array(subsetdata)
                norm_data = pp.normalize(data_arr)

                # Reduce the dimensions down to n principal components
                pca = PCA(n_components=components)
                pca.fit(norm_data)
                pcadata = pca.transform(norm_data)

                # Assemble DataFrame for each subset of genes and append to coords
                coordnames = coordNames[:components]
                coordlists = []
                for i in range(components):
                    coordlists.append(list(pcadata[:, i]))
                zipped = list(zip(coordnames, coordlists))
                df = pd.DataFrame(dict(zipped))
                df.insert(0, "SampleID", sampleID)
                coords.append(df)

                # self.pc_variance = pca.explained_variance_ratio_
                # insert(1,"Coordinates", coords)

            # Construct final analysisDF
            self.analysisDF.insert(1, 'Coordinates', coords)

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
            self.analysisDF = pd.DataFrame({'AnalysisIDs': [0], 'Coordinates': [df]})

    def cluster(self, transformation=None, method='kmeans_classic', num_clusters=2):

        if method != 'kmeans_classic':
            print("Clustering method is not supported!")

        else:
            from sklearn.cluster import KMeans
            import numpy as np

            if (transformation == "PCA"):

                if hasattr(self, 'analysisDF') != True:
                    print("Perform PCA first!")

                else:
                    metadatatype = self.metadata.columns[2:]
                    metadataseps = []

                    # Run K-Means on each subset PCA and add cluster-quality metric for each label type to DataFrame
                    for i in range(len(self.analysisDF)):
                        seplevels = []
                        kmeans = KMeans(n_clusters=num_clusters)
                        arr = np.array(mg.analysisDF.iloc[:, 1][i].iloc[:, [1, 2]])
                        kmeans.fit(arr)
                        labels = list(kmeans.labels_)
                        seplevels = []
                        for i in range(len(self.metadata.columns[2:])):
                            ls = list(self.metadata.iloc[:, i + 2])
                            df = pd.DataFrame({'Labels': labels, 'Type': ls})
                            ct = pd.crosstab(df['Labels'], df['Type'])
                            clusterpurities = []
                            for j in range(len(ct.iloc[:, 0])):
                                clusterpurity = (max(ct.iloc[j]) / ct.iloc[j].sum())
                                clusterpurities.append(clusterpurity)
                            seplevels.append(sum(clusterpurities) / len(ct.iloc[:]))
                        metadataseps.append(seplevels)
                    seps = np.array(metadataseps)
                    for i in range(len(metadatatype)):
                        self.analysisDF[metadatatype[i]] = seps[:, i]



            else:
                # Runs clustering on using ALL gene expressions
                # Create the a kmeans object and fit it to the data
                metadatatype = self.metadata.columns[2:]
                kmeans = KMeans(n_clusters=num_clusters)
                kmeans.fit(self.data)
                # Get cluster centers and cluster labels
                labels = np.array(kmeans.labels_)
                seplevels = []
                for i in range(len(self.metadata.columns[2:])):
                    ls = list(self.metadata.iloc[:, i + 2])
                    df = pd.DataFrame({'Labels': labels, 'Type': ls})
                    ct = pd.crosstab(df['Labels'], df['Type'])
                    clusterpurities = []
                    for j in range(len(ct.iloc[:, 0])):
                        clusterpurity = (max(ct.iloc[j]) / ct.iloc[j].sum())
                        clusterpurities.append(clusterpurity)
                    seplevels.append(sum(clusterpurities) / len(ct.iloc[:]))
                for i in range(len(metadatatype)):
                    self.analysisDF[metadatatype[i]] = seplevels[i]

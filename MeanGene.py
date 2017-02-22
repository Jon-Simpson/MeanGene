class MeanGene(object):

    def __init__(self, file ):
        self.file = file

    def cluster(self):
        def import_file(file):

            import pandas as pd
            # Import the raw text file from OnRamp DEseq
            data = pd.read_table(file, low_memory=False)

            # Remove the annotation and metric columns & transpose the matrix
            data = data.drop(data.columns[[1, 2, 3, -1, -2, -3, -4, -5, -6]], axis=1).transpose()

            # Make the first row (gene symbols) the column names
            header = data.iloc[0]
            data = data.iloc[1:].rename(columns=header)

            # Remove the genes with zero expressions across all samples
            zero_columns = data.mean() == 0
            data = data.drop(data.columns[zero_columns], axis=1)
            return data

        def transform_data(data):

            import numpy as np

            from sklearn import preprocessing as pp
            data_arr = np.array(data)
            norm_data = pp.normalize(data_arr)

            from sklearn.decomposition import PCA

            pca = PCA(n_components=2)
            pca.fit(norm_data)
            finalData = pca.transform(norm_data)
            return finalData

        def kmeans(finalData):

            from sklearn.cluster import KMeans
            import numpy as np

            kmeans = KMeans(n_clusters=2)
            kmeans.fit(finalData)

            centroids = kmeans.cluster_centers_
            labels = np.array(kmeans.labels_)
            return centroids, labels

        def results(data, finalData, centroids, labels):

            group1_name = data.index[0].split(' ')[0]
            group2_name = data.index[-1].split(' ')[0]
            group1_indices = data.index.str.contains(group1_name)
            group2_indices = data.index.str.contains(group2_name)
            group1 = finalData[group1_indices]
            group2 = finalData[group2_indices]
            group1_labels = labels[group1_indices]
            group2_labels = labels[group2_indices]

            import matplotlib.pyplot as plt

            colors = ['b.', 'r.']

            print(group1_name, "Triangle")
            print(group2_name, "Square\n\nCluster 1: Blue\nCluster 2: Red\n")

            for i in range(len(group1)):
                plt.plot(group1[i, 0], group1[i, 1], colors[group1_labels[i]], marker="^", markersize=10)

            for i in range(len(group2)):
                plt.plot(group2[i, 0], group2[i, 1], colors[group2_labels[i]], marker="s", markersize=10)

            plt.title("Cluster Analysis using %s genes" % (len(data.iloc[0])))
            plt.xlabel('PC1')
            plt.ylabel('PC2')
            plt.show()

        data = import_file(self.file)
        finalData = transform_data(data)
        centroids, labels = kmeans(finalData)
        results(data, finalData, centroids, labels)

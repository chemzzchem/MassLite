import numpy as np
from scipy.cluster.hierarchy import linkage, dendrogram
import matplotlib.pyplot as plt

data = np.array([[1, 2, 3],
                 [4, 5, 6],
                 [7, 8, 9]])

data = np.array([[1, 2, 3],
                 [4, 5, 6],
                 [7, 8, 9],
                 [10, 11, 12]])

data = np.array([[1],
                 [4],
                 [8],
                 [10]])

data = np.array([[1],[2],[3]])
linkage_matrix = linkage(data, method='average')

plt.figure(figsize=(10, 6))
dendrogram(linkage_matrix)
plt.title('Hierarchical Clustering Dendrogram')
plt.xlabel('Sample Index')
plt.ylabel('Distance')
plt.show()

from scipy.cluster.hierarchy import fcluster

threshold_distance = 10  # Adjust this threshold as needed
clusters = fcluster(linkage_matrix, t=threshold_distance, criterion='distance')
print("Cluster assignments:", clusters)

print("done")

array([[   32.        ,   150.02671753,  8223.78125   ],
       [   43.        ,   150.02671753,  9096.1875    ],
       [   64.        ,   150.02671753, 10873.625     ],
       [   77.        ,   150.02671753, 11385.3125    ],
       [   72.        ,   150.02671753, 13448.875     ],
       [   51.        ,   150.02671753, 14756.        ],
       [   57.        ,   150.02671753, 20172.        ],
       [   59.        ,   150.02671753, 21444.375     ],
       [   37.        ,   150.02936893, 11403.3125    ],
       [   37.        ,   150.03452047, 12097.96875   ]])
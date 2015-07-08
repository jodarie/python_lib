import sys
import numpy as np
import matplotlib.pyplot as plt, matplotlib.gridspec as gridspec
from sklearn.cluster import DBSCAN
from sklearn.cluster import MeanShift, estimate_bandwidth
from sklearn import metrics


sys.path.append('/Users/Rowen/Documents/Codes/Python/')

# READ DATA

DIR = '/Users/Rowen/Documents/Projects/Euclid/challenge2/blind/svm_test/halo_test/'

train_fname = DIR + 'training_around'

train_data = np.genfromtxt(train_fname, dtype = 'S', unpack = True)

def extract_x(data):
    x0 = np.array(data[3, :], dtype = 'float')
    x1 = np.array(data[4, :], dtype = 'float')
    x2 = np.array(data[5, :], dtype = 'float')
    x3 = np.array(data[6, :], dtype = 'float')
    x4 = np.array(data[7, :], dtype = 'float')
    x5 = np.array(data[8, :], dtype = 'float')
    x6 = np.array(data[9, :], dtype = 'float')
    x7 = np.array(data[10, :], dtype = 'float')

    x45 = x4 - x5
    x67 = x6 - x7
    
    X = np.transpose(np.vstack((x45, x67)))
    return X

X = extract_x(train_data)

print 'Finished reading data...'

# Clustering

## DBSCAN

#eps = 0.2
#samples = 1

#db = DBSCAN(eps = eps, min_samples = samples).fit(X)

#core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
#core_samples_mask[db.core_sample_indices_] = True
#labels = db.labels_

# Number of clusters in labels, ignoring noise if present.
#n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)

#print('Estimated number of clusters: %d' % n_clusters_)

## MEAN SHIFT

bandwidth = estimate_bandwidth(X, quantile=0.4)

ms = MeanShift(bandwidth=bandwidth, bin_seeding=True)
ms.fit(X)
labels = ms.labels_
cluster_centers = ms.cluster_centers_

labels_unique = np.unique(labels)
n_clusters_ = len(labels_unique)

print("number of estimated clusters : %d" % n_clusters_)

#for i in range(len(labels)):
#    print train_data[12, i], labels[i]
    
print 'Finished clustering...'

# Plot result
import matplotlib.pyplot as plt
from itertools import cycle

plt.figure(1)
plt.clf()

colors = cycle('bgrcmykbgrcmykbgrcmykbgrcmyk')
for k, col in zip(range(n_clusters_), colors):
    my_members = labels == k
    cluster_center = cluster_centers[k]
    plt.plot(X[my_members, 0], X[my_members, 1], col + '.')
    plt.plot(cluster_center[0], cluster_center[1], 'o', markerfacecolor=col,
             markeredgecolor='k', markersize=14)
plt.title('Estimated number of clusters: %d' % n_clusters_)
plt.show()

# Print Labels

for i in range(len(labels)):
    print train_data[12, i], labels[i]

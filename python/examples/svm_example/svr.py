import sys
import numpy as np
import matplotlib.pyplot as plt, matplotlib.gridspec as gridspec
from sklearn import svm

sys.path.append('/Users/Rowen/Documents/Codes/Python/')

from functions import astro

# READ DATA

#DIR = '/Users/Rowen/Documents/Projects/Euclid/challenge2/blind/svm_test/sf_test/'
#DIR = '/Users/Rowen/Documents/Projects/Euclid/challenge2/blind/svm_test/fb_test/'
#DIR = '/Users/Rowen/Documents/Projects/Euclid/challenge2/blind/svm_test/halo_test/'
DIR = '/Users/Rowen/Documents/Projects/Euclid/challenge2/blind/pieces_7x7/sn_test/'

#train_fname = DIR + 'sf_train_data_mass'
#test_fname = DIR + 'sf_test_data_mass'

train_fname = DIR + 'temp1'

train_data = np.genfromtxt(train_fname, dtype = 'float', unpack = True)
#test_data = np.genfromtxt(test_fname, dtype = 'float', unpack = True)

def extract_xy(data):
    z = np.array(data[0, :], dtype = 'float')
    rich = np.log10(np.array(data[1, :], dtype = 'float'))
    sn = np.array(data[2, :], dtype = 'float')
    #sn2 = np.array(data[3, :], dtype = 'float')
    #size = np.array(data[6, :], dtype = 'float')
    #area = np.array(data[7, :], dtype = 'float')
    #mass = np.array(data[9, :], dtype = 'float')
    mass = np.array(data[3, :], dtype = 'float')
    X = np.transpose(np.vstack((z, rich, sn)))
    Y = mass
    return X, Y

X_train, Y_train = extract_xy(train_data)
#X_test, Y_test = extract_xy(test_data)

print 'Finished reading data...'

# SVM

kernel = 'rbf'
gamma = 1.0
C = 10.0

clf = svm.SVR(kernel = kernel, gamma = gamma, C = C)

# No Weights

clf.fit(X_train, Y_train)
print 'Score 1:', clf.score(X_train, Y_train)
#print 'Score 2:', clf.score(X_test, Y_test)

# With Weights

#clf.fit(X_train, Y_train, sample_weight = w_train)
#print 'Score 1:', clf.score(X_train, Y_train, sample_weight = w_train)
#print 'Score 2:', clf.score(X_test, Y_test, sample_weight = w_test)

print 'Finished fitting SVC...'

# PREDICT

#y_pred = clf.predict(X_test)
y_pred = clf.predict(X_train)

print 'Finished making preditions...'
    
# PLOT

def make_plot(x, y):
    fig = plt.figure(1)
    gs = gridspec.GridSpec(1, 1)
    ax = fig.add_subplot(gs[0])
    ax.scatter(x, y)
    ax.plot([13.0, 15.0], [13.0, 15.0], '--', c = 'r')
    ax.set_xlim(13.0, 15.0)
    ax.set_ylim(13.0, 15.0)
    ax.set_xlabel('True Mass')
    ax.set_ylabel('Estimated Mass')
      
    plt.show()
    
#make_plot(Y_test, y_pred)
make_plot(Y_train, y_pred)

print 'Finished producing plot...'

# OUTPUT

#out_file = test_fname + '_predictions'
out_file = train_fname + '_predictions'                                          
output = open(out_file,'w')

for i in range(len(y_pred)):
    #print>> output, X_test[i, 0], Y_test[i], y_pred[i]
    print>> output, X_train[i, 0], Y_train[i], y_pred[i]

print 'Output written to: ' + out_file

import sys
import numpy as np
import matplotlib.pyplot as plt, matplotlib.gridspec as gridspec
from sklearn import svm

sys.path.append('/Users/Rowen/Documents/Codes/Python/')

from functions import astro

# READ DATA

DIR = '/Users/Rowen/Documents/Projects/Euclid/challenge2/blind/svm_test/sf_test/'

train_fname = DIR + 'match_data_train'
test_fname = DIR + 'match_data_test'

train_data = np.genfromtxt(train_fname, dtype = 'S', unpack = True)
test_data = np.genfromtxt(test_fname, dtype = 'S', unpack = True)

def f_scale(value):
    return (value - min(value)) / (max(value) - min(value))

def extract_xy(data):
    z = f_scale(np.array(data[3, :], dtype = 'float'))
    rich = f_scale(np.array(data[4, :], dtype = 'float'))
    sn = f_scale(np.array(data[5, :], dtype = 'float'))
    #size = f_scale(np.array(data[6, :], dtype = 'float'))
    #area = f_scale(np.array(data[7, :], dtype = 'float'))
    match = np.array(data[6, :], dtype = 'int')
    #X = np.transpose(np.vstack((z, rich, sn, size, area)))
    X = np.transpose(np.vstack((z, rich, sn)))
    Y = match
    return X, Y

X_train, Y_train = extract_xy(train_data)
X_test, Y_test = extract_xy(test_data)

print 'Finished reading data...'

# SVM

kernel = 'rbf'
gamma = 10.0
C = 100

clf = svm.SVC(kernel = kernel, gamma = gamma, C = C)
clf.fit(X_train, Y_train)
print 'Train Score:', clf.score(X_train, Y_train)
print 'Test Score:', clf.score(X_test, Y_test)

#gamma_list = 10.0 ** np.arange(-5, 4)
#C_list = 10.0 ** np.arange(-2, 4)

#for gamma in gamma_list:
#   for c in C_list:
#      clf = svm.SVC(kernel = kernel, gamma = gamma, C = c)
#     clf.fit(X_train, Y_train)
#    score = clf.score(X_test, Y_test)
#   print 'Gamma: ', gamma, 'C: ', c, 'Score: ', score

print 'Finished fitting SVC...'

# PREDICT

y_pred = clf.predict(X_test)

print 'Finished making preditions...'
    
# PLOT

def make_plot(x, y):
    fig = plt.figure(1)
    gs = gridspec.GridSpec(1, 1)
    ax = fig.add_subplot(gs[0])
    
    ax.hist(Y_test, facecolor='green', alpha=0.5, label='Y')
    ax.hist(y_pred, facecolor='red', alpha=0.5, label='Y-Predicted')
    
    ax.legend()
      
    plt.show()
    
make_plot(Y_test, y_pred)

print 'Finished producing plot...'

# OUTPUT

out_file = test_fname + '_selection'                                          
output = open(out_file,'w')

for i in range(len(y_pred)):
    if y_pred[i] == 1:
        print>> output, test_data[0, i], test_data[1, i], test_data[2, i], test_data[3, i], test_data[4, i], test_data[5, i]

print 'Output written to ' + out_file + '...'

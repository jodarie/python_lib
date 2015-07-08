# Cluster Mass Neural Network Example

import numpy as np
import neurolab as nl

# Read Data

DIR = '/Users/Rowen/Documents/Projects/Euclid/challenge2/blind/svm_test/sf_test/'

train_fname = DIR + 'sf_train_data_mass'
test_fname = DIR + 'sf_test_data_mass'

train_data = np.genfromtxt(train_fname, dtype = 'S', unpack = True)
test_data = np.genfromtxt(test_fname, dtype = 'S', unpack = True)

# Random Sample

index_array = np.arange(train_data.shape[1])
np.random.shuffle(index_array)

train_data_random = train_data[:, index_array[:1000]]

def f_scale(value):
    return (value - min(value)) / (max(value) - min(value))

def extract_xy(data):
    z = np.array(data[3, :], dtype = 'float')
    rich = np.log(np.array(data[4, :], dtype = 'float'))
    sn = np.array(data[5, :], dtype = 'float')
    size = np.array(data[6, :], dtype = 'float')
    area = np.array(data[7, :], dtype = 'float')
    mass = np.array(data[9, :], dtype = 'float')
    X = np.transpose(np.vstack((f_scale(z), f_scale(rich), f_scale(sn))))
    Y = f_scale(mass).reshape(len(mass), 1)
    return X, Y

X_train, Y_train = extract_xy(train_data)
X_test, Y_test = extract_xy(test_data)

# Neural Net

arch = [10, 1]
epochs = 500
goal = 0.02
threshold = 100

def nn(X, Y, arch, epochs, goal):
    print 'Training Neural Network...'
    stuff = []
    for i in range(X.shape[1]):
        stuff.append([min(X[:, i]), max(X[:, i])])
    net = nl.net.newff(stuff, arch)
    error = net.train(X, Y, epochs = epochs, goal = goal, show = 1)
    return net, error

network, errors = nn(X_train, Y_train, arch, epochs, goal)

while errors[-1] > threshold:
    network, errors = nn(X_train, Y_train, arch, epochs, goal)

y_pred = network.sim(X_train)

# Plot results

import pylab as pl
pl.subplot(211)
pl.plot(errors)
pl.xlabel('Epoch number')
pl.ylabel('error (default SSE)')

pl.subplot(212)
pl.scatter(Y_train, y_pred)
pl.xlabel('Y-Test')
pl.ylabel('Y-NN')
pl.show()

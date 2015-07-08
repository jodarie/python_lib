"train a regression MLP"

import numpy as np
from math import sqrt
from pybrain.datasets.supervised import SupervisedDataSet as SDS
from pybrain.tools.shortcuts import buildNetwork
from pybrain.supervised.trainers import BackpropTrainer


DIR = '/Users/Rowen/Documents/Projects/Euclid/challenge2/calib/p_lambda/sf_test/'

train_fname = DIR + 'sf_train_data_mass'
test_fname = DIR + 'sf_test_data_mass'

hidden_size = 4
epochs = 600

train_data = np.genfromtxt(train_fname, dtype = 'S', unpack = True)
test_data = np.genfromtxt(test_fname, dtype = 'S', unpack = True)

def extract_xy(data):
    z = np.array(data[3, :], dtype = 'float')
    rich = np.log(np.array(data[4, :], dtype = 'float'))
    sn = np.array(data[5, :], dtype = 'float')
    size = np.array(data[6, :], dtype = 'float')
    area = np.array(data[7, :], dtype = 'float')
    mass = np.array(data[9, :], dtype = 'float')
    #X = np.transpose(np.vstack((z, rich, sn, size, area)))
    X = np.transpose(np.vstack((z, rich)))
    Y = mass
    return X, Y

X_train, Y_train = extract_xy(train_data)
X_test, Y_test = extract_xy(test_data)

Y_train = Y_train.reshape( -1, 1 )

input_size = X_train.shape[1]
target_size = Y_train.shape[1]

# prepare dataset

ds = SDS( input_size, target_size )
ds.setField( 'input', X_train )
ds.setField( 'target', Y_train )

# init and train

net = buildNetwork( input_size, hidden_size, target_size, bias = True )
trainer = BackpropTrainer( net,ds )

print "training for {} epochs...".format( epochs )

for i in range( epochs ):
	mse = trainer.train()
	rmse = sqrt( mse )
	print "training RMSE, epoch {}: {}".format( i + 1, rmse )
	

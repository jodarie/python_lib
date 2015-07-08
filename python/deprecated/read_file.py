import numpy as np

f2 = open('datas.dat', 'r')
lines = f2.readlines()
f2.close()

data = []

for line in lines:
    p = line.split()
    data.append(p)

print(data[0])

data2 = np.array(data)
data2 = data2.astype(float)

print(data2.shape)

x = data2[:,0]
y = data2[:,1]

print x

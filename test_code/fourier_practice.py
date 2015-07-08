import math
import extra_math as emath
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

L = 20
n = 128

x = np.array(np.linspace(-L/2, L/2, n + 1))[1:n + 1]

u = -np.cos(x)
ud = np.sin(x)

k = emath.k_val(n, L)

plt.plot(x, ud, label = 'True')
plt.plot(x, emath.fourier_derivative(u, k, 1), 'rx', label = 'Approximate')
plt.legend()
plt.title('sin(x)')
plt.show()


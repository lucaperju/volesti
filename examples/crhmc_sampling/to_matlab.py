from scipy.io import savemat
from scipy.io import loadmat
import numpy as np
import sys
import math

dim = int(sys.argv[1])

if int(sys.argv[2]) == 2:
    dim = math.ceil(math.sqrt(dim))
    dim += 1
    dim = dim * (dim - 2) + 1

matA = np.loadtxt('cpp_poly.txt', usecols = range(dim))

vecB = np.loadtxt('cpp_poly.txt', usecols = (dim), unpack = True)
data = {'Aineq': matA, 'bineq':vecB}

savemat('matlab_poly.mat', data)
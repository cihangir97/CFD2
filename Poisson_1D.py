import numpy as np
import math
from matplotlib import pyplot as plt
from scipy.sparse.linalg import spsolve

###
# Solves the 1d poisson equation
# with Dirichlet boundary condtion at x=0 of phi=0
# and Neumann boundary condtion at x=l
# where the first derivative wrt x equals 0 at this point
###

## Test function

def fx(x):
	y = np.sin(((3*np.pi*x)/2))
	return y

## Exact solution
def phix(x):
	y = (-4/(9*(np.pi)**2))*np.sin((3*np.pi*x)/2)
	return y

## Constants

l	= 1.0	# 1D domain length
N	= 5	# Grid points on [0,l]
h	= l/N	# Distance between points (constant)
x 	= np.linspace(0,l,N+1)

## Create linear system
# Tri diagonal matirx [-1 2 -1]
A =	(
		np.diag(np.full(N,-2)) + 
	 	np.diag(1*np.ones((N-1)),-1) + 
	 	np.diag(1*np.ones((N-1)),1)   
	)
# Implement neumann bc in matrix
A[N-1, N-2] = 2 # N-1 since python starts index 0 
# Divide everything by h squared
A = A / (h**2)
# right-hand side
fi = fx(x[1:])

## Solve the system
phii = np.zeros(N+1) # Add dirichelet bc.
phii[1:] = np.linalg.solve(A,fi)

## Compute errors

diff = phix(x) - phii 
eps = np.linalg.norm(diff) # L2 norm

plt.figure()
plt.title('Exact and approximate solution, N=%i'%N)
plt.grid()
plt.plot(x, phix(x), color='b', linewidth=1, label='Exact Solution')	# exact
plt.scatter(x, phii, color='r', label='2nd Order FD')			# approx 
plt.legend()

plt.figure()
plt.title('Error in the points, N=%i' %N)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.ylim(1.2*np.min(diff), 1.2*np.max(diff))
plt.scatter(x, diff, color='r')
plt.grid()
plt.show()

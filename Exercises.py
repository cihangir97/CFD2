from FD_functions import *

approx, exact, diff, eps, h  = Poisson_1D(10, False)

## Question 1, 2 and 3

approxN5, exactN5, diffN5, epsN5, h5  = Poisson_1D(5, False)
approxN10, exactN10, diffN10, epsN10, h10  = Poisson_1D(10, False)
approxN20, exactN20, diffN20, epsN20, h20  = Poisson_1D(20, False)
approxN40, exactN40, diffN40, epsN40, h40  = Poisson_1D(40, False)

factorN5toN10 = epsN5 / epsN10
factorN10toN20 = epsN10 / epsN20 
factorN20toN40 = epsN20 / epsN40 

averageFactor = (factorN5toN10 + factorN10toN20 + factorN20toN40)/(3)
print '############# Question 1 & 2 #############'
print 'Factor going from N=5 to N=10 is: %f' %factorN5toN10
print 'Factor going from N=10 to N=20 is: %f' %factorN10toN20
print 'Factor going from N=20 to N=40 is: %f' %factorN20toN40
print 'When the nodes are doubled, the L2 error reduces by %f' %averageFactor 

print '############# Question 3 #############'
print 'The Poisson equation contains a second order derivative.'
print 'The error is the largest where this 2nd order derivative is the largest.'
print 'In other words, the error is the largest where the slope of the function'
print 'changes the quickly.'
print 'So usualy at a maxima or minima, where the slope goes from positive to negative'
print '0r the other way around.'


## Question 4


Ns = 2**(np.arange(14)[4:])	# array with our Nvalues from 2**4 till specified

eps_hs	= [Poisson_1D(i,False)[3:] for i in Ns] # get all the L2 and h for different Ns
eps	= zip(*eps_hs)[0]			# all the L2 errors
hs	= zip(*eps_hs)[1]			# all the hs 

eps_log = np.log(eps)
hs_log = np.log(hs)

slope, intercept = np.polyfit(hs_log, eps_log, 1)

print 'The slope of this line equals %f' %slope
print 'Which means that this method has order %f convergence' %slope

plt.figure()
plt.title('log of the error agains the log of the mesh width, with slope = %f' %slope)
plt.xlabel('log(h)')
plt.ylabel('log(L2-error)')
plt.plot(hs_log, eps_log)
plt.gca().invert_xaxis()
plt.grid()

plt.show()

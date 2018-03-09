import numpy as np
import matplotlib.pyplot as plt
#%matplotlib inline 
%pylab
plt.rcParams['font.size']=20

# Generate synthetic data from a model.
# For simplicity, let us assume a LINEAR model y = m*x + b
# where we want to fit m and b
m_true = -0.9594
b_true = 4.294
N = 50

x = np.sort(10*np.random.rand(N))
yerr = 0.2 + 0.5*np.random.rand(N)
y = m_true*x + b_true
y += yerr * np.random.randn(N)

fig = plt.figure()
#fig.set_size_inches(12, 8)
plt.errorbar(x, y, yerr=yerr, fmt='.k')



# Now, let's setup some parameters that define the MCMC
ndim = 2
nwalkers = 500

# Initialize the chain
# Choice 1: chain uniformly distributed in the range of the parameters
pos_min = np.array([-5., 0.])
pos_max = np.array([5., 10.])
psize = pos_max - pos_min
pos = [pos_min + psize*np.random.rand(ndim) for i in range(nwalkers)]

# Visualize the initialization
import triangle
fig = triangle.corner(pos, labels=["$m$", "$b$"], extents=[[-5., 5.], [0., 10.]], 
        truths=[m_true, b_true])
fig.set_size_inches(10,10)
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.stats import norm 
import scipy.spatial as sp

class Soil(object):
    def __init__(self):
        self.c        = 1.
        self.phi      = 1.       # TAN(PHI)
        self.rho_cphi = -0.5        
  
class Prop(object):
    def __init__(self,m,s):
        self.mean    = m
        self.std     = s
        self.cov     = s/(m+1.e-12)
        self.mu_LN   = np.log(m+1.e-12/np.sqrt(1.+self.cov**2))
        self.sig_LN  = np.sqrt(np.log(1.+self.cov**2))

soil = Soil()
soil.c   = Prop(20,3)
soil.phi = Prop(0.57735026919,0.11547005383)

nx = 40
ny = 20
#x = np.linspace(-4, 4, 9)
x = np.linspace(-2, 2, nx)
y = np.linspace(-1, 1, ny)
tx = 20.
ty = 4.

N = nx*ny

#U = np.random.random((nx, ny)) 
Uc = np.random.normal(0,1,[N,1])
Up = np.random.normal(0,1,[N,1])

xx,yy = np.meshgrid(x,y)
#print ('xx=',xx)
#plt.plot(xx, yy, marker='.', color='k', linestyle='none')

X = np.array([xx.reshape(-1,)/tx,yy.reshape(-1,)/ty]).T

D = sp.distance_matrix(X,X)

#Dx = np.abs(xx-np.transpose(xx))  #delta x
#Dy = np.abs(yy-np.transpose(yy))  #delta y

#tau = np.sqrt((Dx/tx)**2 + (Dy/ty)**2) # ELLIPTIC EXPONENTIAL

tau = D

RR = np.exp(-2.*tau)
RRc= np.ones([2,2])
RRc[0,1] = soil.rho_cphi
RRc[1,0] = soil.rho_cphi

AA = np.linalg.cholesky(RR)
AAc  = np.linalg.cholesky(RRc)

Zc = AA.dot(Uc)
Zp   = AAc[1,0]*Zc + AAc[1,1]*AA.dot(Up)

c  = np.transpose(np.exp(soil.c.mu_LN + soil.c.sig_LN*Zc))


tphi = np.transpose(np.exp(soil.phi.mu_LN + soil.phi.sig_LN*Zp))
phi = np.rad2deg(tphi)

C = c.reshape(ny,nx)

plt.contourf(x, y, C)
#plt.imshow(c.reshape(ny,nx))
#plt.axis('off')
cbar = plt.colorbar(orientation='horizontal', fraction=.1)
#cbar.set_ticks(np.arange(0,50,5))
plt.show()
#plt.contourf(xx, yy, phi, cmap = 'Greys')
#plt.axis('off')
#plt.colorbar(orientation='horizontal', fraction=.1)
plt.imshow(C)


'''
fig, (ax1, ax2) = plt.subplots(2, 1)

im1 = ax1.contourf(xx, yy, c, cmap = 'Greys')
divider = make_axes_locatable(ax1)
cax1 = divider.append_axes("right", size="2%", pad=0.05)
fig.colorbar(im1, cax=cax1)
ax1.set_axis_off()

im2 = ax2.contourf(xx, yy, tphi, cmap = 'Greys')
divider = make_axes_locatable(ax2)
cax2 = divider.append_axes("right", size="2%", pad=0.05)
fig.colorbar(im2, cax=cax2)
ax2.set_axis_off()

fig.tight_layout()
plt.show()
'''



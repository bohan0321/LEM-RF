import matplotlib.pyplot as plt
plt.close('all')
f = plt.figure()

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.stats import norm 
import scipy.spatial as sp


class Res(object):
    def __init__(self):
        self.Bmcs


class Geometry(object):
    def __init__(self):
        self.H      = 10.       
        self.D      = 1.5         
        self.beta   = np.arctan(1.0)
        self.nslice = 5

class Soil(object):
    def __init__(self):
        self.c        = 1.
        self.phi      = 1.       # TAN(PHI)
        self.rho_cphi = -0.5
        self.gamma    = 20.
        self.t1       = 20.
        self.t2       = 4
        self.B        = 0.4
        self.fcor    = 'EExp'
        
class Prop(object):
    def __init__(self,m,s):
        self.mean    = m
        self.std     = s
        self.cov     = s/(m+1.e-12)
        self.mu_LN   = np.log(m+1.e-12/np.sqrt(1.+self.cov**2))
        self.sig_LN  = np.sqrt(np.log(1.+self.cov**2))        

class LEMMeth(object):
    def __init__(self):
        self.BS  = True 

class LEMType(object):
    def __init__(self):
        self.det  = True        ## Deterministic
        self.mcs1 = False       ## MonteCarlo
        self.mcs2 = False       ## MonteCarlo
    def count(self):
        self.n = 0
        if self.det == True:
            self.n+= 1
        if self.mcs1 == True:
            self.n+= 1
        if self.mcs2 == True:
            self.n+= 1     

class LEM(object):
    def __init__(self):
        self.Meth = LEMMeth()
        self.Type = LEMType()
        self.plot = False
        self.prnt = False
        self.nfld = 10
        self.nx   = 30
        self.ny   = 30
        self.nr   = 16
        self.rx = np.array([-10,10])
        self.ry = np.array([10.01,20.]) 

class RF(object):
    def __init__(self):
        self.nx   = 5
        self.ny   = 5
        self.domainx = np.array([-20,10])
        self.domainy = np.array([-5,10]) 
        
        
def radius(Xc, Yc, geom):
    xtop = geom.H/np.tan(geom.beta)
    Rmax = Yc + (geom.D - 1)*geom.H
    Rmin1 = np.sqrt((xtop+Xc)**2 + (Yc-geom.H)**2)
    Rmin2 = np.sqrt(Xc**2 + Yc**2)
    Rmin = max(Rmin1, Rmin2)
    return Rmin, Rmax   

'''
def dim(Xc, Yc, R, geom):
    h = Yc - geom.H
    if h == 0:
        a1 = 0
    else:
        a1 = np.arcsin(h/R)
    a2 = np.pi-np.arcsin(Yc/R)
    return h, a1, a2
'''

def CreateR(xrf,yrf,soil):
    
    xxrf,yyrf = np.meshgrid(xrf,yrf)
    #  
    tx = soil.t1
    ty = soil.t2
    #
    Xrf = np.array([xxrf.reshape(-1,)/tx,yyrf.reshape(-1,)/ty]).T
    D = sp.distance_matrix(Xrf,Xrf)
    tau = D
    #
    RR = np.exp(-2.*tau)
    
    return RR


def LEM_plot(Xc, Yc, R, geom):    

    h,a1,a2 = dim(Xc, Yc, R, geom)          # 
    a = np.linspace(a1,a2,geom.nslice+1)    # hoeken van alle slice intersections
    x = Xc+R*np.cos(np.pi+a)                # coordinaten [..]
    y = Yc+R*np.sin(np.pi+a)

    plt.plot([-geom.H/np.tan(geom.beta)-10,-geom.H/np.tan(geom.beta),0,10],[geom.H,geom.H,0,0],'b')
    if geom.nslice<25:
        plt.plot(x,y,'-|r')
        ht = geom.H-y                           # centre height of slice
        ht[x>geom.H/np.tan(-geom.beta)] = -y[x>geom.H/np.tan(-geom.beta)]+x[x>geom.H/np.tan(-geom.beta)]*np.tan(-geom.beta)
        ht[x>0.0] = -y[x>0.0]
        
        plt.plot([x,x],[y,y+ht],'k',linewidth=0.3)
                
    else:
        plt.plot(x,y,'-r')            
    plt.plot([x[0],Xc,x[-1]],[y[0],Yc,y[-1]],':k',linewidth=0.5)
    plt.plot([-geom.H/np.tan(geom.beta)-10,10],[-(geom.D-1.)*geom.H,-(geom.D-1.)*geom.H],'k')
    plt.plot(Xc,Yc,'+k')
    plt.axis('equal')
    #plt.xlim([-20,10])
    plt.ylim([-(geom.D-.9)*geom.H,(geom.D+2)*geom.H])


def LEM_S(Xc,Yc,R,geom,soil,LEM,U):
    #
    res = np.zeros(LEM.Type.n)
    ires = -1
    #
    h,a1,a2 = dim(Xc, Yc, R, geom)          # 
    a = np.linspace(a1,a2,geom.nslice+1)    # hoeken van alle slice intersections
    x = Xc+R*np.cos(np.pi+a)                # coordinaten [..]
    ##
    ## GENERAL 
    ## 
    ac = a[0:geom.nslice]+np.diff(a)/2      # hoeken van alle slice centres
    T = np.pi/2-ac                          # angles between slices and horizontal
    xc = Xc+R*np.cos(np.pi+ac)              # coordinates of slice centre of slide
    yc = Yc+R*np.sin(np.pi+ac)
    #    
    h = geom.H-yc                           # centre height of slice
    h[xc>geom.H/np.tan(-geom.beta)] = -yc[xc>geom.H/np.tan(-geom.beta)]+xc[xc>geom.H/np.tan(-geom.beta)]*np.tan(-geom.beta)
    h[xc>0.0] = -yc[xc>0.0]

    ## WEIGHTS AND DRIVING MOMENTS
    ##
    b = np.diff(x)          # width of slices
    w = h*b*soil.gamma      # weight of slices 
    W = w*(xc-Xc)           # moments of slices around centre of rotation [Xc,Yc]
    L = R*(a2-a1)           # length of sliding surface
    l = L/geom.nslice            # length of slice sliding surface

    if LEM.Meth.BS == True:
        
        if LEM.Type.det == True:
            ires+=1
            F = 1.0
            conv = False
            itr  = 0
            c    = soil.c.mean
            tphi = soil.phi.mean
            while conv == False: 
                itr+=1
                S = (c+(h*soil.gamma)*tphi)*l*R/((1.+np.tan(T)*tphi/F))
                FS = -np.sum(S)/np.sum(W)
                if np.abs(F-FS)<0.000001:
                    conv = True
                else:
                    F = FS
                if itr>50:
                    print('NO CONVERGENCE IN BISHOP ITERATION')
            res[ires] = FS
    return res            



soil = Soil()
LEM = LEM()
RF = RF()
geom = Geometry()
geom.beta = np.tan(np.deg2rad(30))





soil.c   = Prop(20,3)
soil.phi = Prop(0.57735026919,0.11547005383)

LEM.Type.count()

nx = LEM.nx
ny = LEM.ny
nr = LEM.nr


Xc   = np.linspace(LEM.rx[0], LEM.rx[1], LEM.nx)
Yc   = np.linspace(LEM.ry[0], LEM.ry[1], LEM.ny)


xrf  = np.linspace(RF.domainx[0], RF.domainx[1], RF.nx)
yrf  = np.linspace(RF.domainy[0], RF.domainy[1], RF.ny)


N = RF.nx*RF.ny

Uc = np.random.normal(0,1,[N,1])
Up = np.random.normal(0,1,[N,1])

    ##
    ## CALCULATE R
    ##
RR = CreateR(xrf,yrf,soil)
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

C = c.reshape(RF.nx,RF.ny)
tPHI = tphi.reshape(RF.nx,RF.ny)

xxrf,yyrf = np.meshgrid(xrf,yrf)



#plt.scatter(xxrf, yyrf, c=C, s=30)
#plt.colorbar(orientation='vertical', shrink=1, pad=0.01)

#plt.show()

Rmin = np.zeros([LEM.nx, LEM.ny])
Rmax = np.zeros([LEM.nx, LEM.ny])
#Xc   = np.linspace(LEM.rx[0], LEM.rx[1], LEM.nx)
#Yc   = np.linspace(LEM.ry[0], LEM.ry[1], LEM.ny)
RES  = np.zeros([LEM.nx, LEM.ny, LEM.nr, LEM.Type.n])
yy   = np.zeros([LEM.nx, LEM.ny, LEM.nr])
rr   = np.zeros([LEM.nx, LEM.ny, LEM.nr])
Xc = -1.7241379310344822
Yc = 19.65551724137931
R = 19.730991096024326

h = Yc - geom.H
a1 = np.arcsin(h/R)
a2 = np.pi-np.arcsin(Yc/R)

a = np.linspace(a1,a2,geom.nslice+1)    # hoeken van alle slice intersections
x = Xc+R*np.cos(np.pi+a)                # coordinaten [..]
y = Yc+R*np.sin(np.pi+a)


plt.plot(Xc,Yc,'+k')    #center of circle

plt.plot(x,y,'-|r')     #slice boundary red
ht = geom.H-y                           # centre height of slice
ht[x>geom.H/np.tan(-geom.beta)] = -y[x>geom.H/np.tan(-geom.beta)]+x[x>geom.H/np.tan(-geom.beta)]*np.tan(-geom.beta)
ht[x>0.0] = -y[x>0.0]
        
plt.plot([x,x],[y,y+ht],'k',linewidth=0.3)      #slice boundary vertical grey
#
plt.plot([x[0],Xc,x[-1]],[y[0],Yc,y[-1]],':k',linewidth=0.5)    #radius dash grey

#plt.axis('equal')
#plt.xlim([-20,10])
#plt.ylim([-(geom.D-.9)*geom.H,(geom.D+2)*geom.H])

#
plt.plot([-geom.H/np.tan(geom.beta)-10,-geom.H/np.tan(geom.beta),0,10],[geom.H,geom.H,0,0],'b') #upper boundary
plt.plot([-geom.H/np.tan(geom.beta)-10,10],[-(geom.D-1.)*geom.H,-(geom.D-1.)*geom.H],'k') #lower boundary
#
plt.scatter(xxrf, yyrf, c=C, s=30)

plt.show()

#f.savefig("foo.pdf", bbox_inches='tight')


'''
nt = 1

for t in range(nt):
    
    if LEM.Type.mcs2 == True:
        U = np.random.normal(0,1,[2*geom.nslice,LEM.nfld])
    else:
        U = []
    
    tt = time.time()
   
    for i in range(len(Xc)):
#        print(i,time.time() - tt)
        for j in range(len(Yc)):
            Rmin, Rmax = radius(Xc[i], Yc[j], geom)
            R = np.linspace(Rmin, 0.8*Rmin+0.2*Rmax, nr)
            for k in range(nr):
                yy[i,j,k] = Yc[j]
                rr[i,j,k] = R[k]
                RES[i,j,k,:] = LEM_S(Xc[i],Yc[j],R[k],geom,soil,LEM,U)

    elapsed = time.time() - tt
    print('elapsed time : ',elapsed)
    
#    F[1,2,3]=0
    
    x0 = LEM.rx[0] - (LEM.rx[1]-LEM.rx[0])/2/nx
    x1 = LEM.rx[1] + (LEM.rx[1]-LEM.rx[0])/2/nx
    y0 = LEM.ry[0] - (LEM.ry[1]-LEM.ry[0])/2/ny
    y1 = LEM.ry[1] + (LEM.ry[1]-LEM.ry[0])/2/ny
    ##
    ## PLOT MINIMUM OF FS_det
    ##
    iplt = 0
    if LEM.Type.det == True:
        iplt+=1
        FS = RES[:,:,:,iplt-1]
        if LEM.Type.n>1:
            plt.subplot(2,2,iplt)
        Fmin = np.amin(FS, 2)
        
        n = np.unravel_index(np.argmin(FS),(nx,ny,nr))
        Xcc = Xc[n[0]]
        Ycc = Yc[n[1]]
        Rmin,Rmax = radius(Xcc,Ycc,geom)
        R  = np.linspace(Rmin,Rmax,nr)
        R = R[n[2]]

        LEM_plot(Xcc, Ycc, R, geom,)    
        print('DTRMN: [Xc,Yc,Rc,FS] = %8.4f,%8.4f,%8.4f,%8.4f '%(Xcc,Ycc,R,np.min(FS)))
        
        plt.scatter(xxrf, yyrf, c=C, s=30)
        plt.imshow(np.transpose(Fmin), extent=[x0,x1,y1,y0])
#        plt.axis([-20,10,-5,10])
        cbar = plt.colorbar()
        cbar.set_label('Factor of Safety', rotation=90)
        plt.show()
'''      



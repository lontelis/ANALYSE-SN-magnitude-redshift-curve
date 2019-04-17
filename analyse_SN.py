import cosmology,fitting
from numpy import *
from pylab import *

# Read data
totoSN = loadtxt('data/supernovaUnion2011.txt',usecols=[1,2,3,4])
z_unsorted  = totoSN[:,0]
mu          = totoSN[:,1]
dmu         = totoSN[:,2]

permutter = sorted(range(len(z_unsorted)), key=lambda k: z_unsorted[k])
z  = z_unsorted[permutter]
mu0 =mu[permutter]
dmu0=dmu[permutter]

#Define model
def model2(x,pars): return log10( cosmology.get_dist(z,type='dl',params=[pars[0],pars[1],-1,0])/0.7)*5.0+25.

#Do fit
res = fitting.dothefit(z,mu0,dmu0,[.3,.7],functname=model2,method='minuit')


# plot
figure(1,figsize=(15,10)),clf()
errorbar(z,mu0,yerr=dmu0,fmt='.', label='SN Union2.1 Suzuki et Al. 2011') #,label='Riess et Al. 1998 High-Z Templates')
xscale('linear'),ylabel('$\mu(z)$',size=25),xlabel('$z$',size=25)
plot(z,model2(z,res[1]),'r-',label='$\mu(z)=5\log_{10} d_{L}(z)+25$, $(\Omega_{m},\Omega_{\Lambda})$=(%0.3f $\pm$ %0.3f,%0.3f $\pm$ %0.3f)'%(res[1][0],res[2][0],res[1][1],res[2][1]) )
legend(loc=4,frameon=False,numpoints=1,fontsize=25)
draw(),show()

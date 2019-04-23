import cosmology,fitting
from numpy import *
from pylab import *


totoSN = loadtxt('data/supernovaUnion2011.txt',usecols=[1,2,3,4])
z_unsorted  = totoSN[:,0]
mu          = totoSN[:,1]
dmu         = totoSN[:,2]
dmu_sys     = totoSN[:,3]

# sort the data according to the redshift, z.
permutter = sorted(range(len(z_unsorted)), key=lambda k: z_unsorted[k])
z  = z_unsorted[permutter] # redshift
mu0 =mu[permutter]      # magnitude
dmu0=dmu[permutter]     # statistical error on the magnitude, use
dmu0_sys=dmu_sys[permutter] # systematic error on the magnitude, no use

# Define your model
def model(x,pars): 
    # Om = pars[0] matter density ratio , free parameter
    # OL = pars[1] dark energy density ratio, free parameter
    # mu = 5*log10( dL(z,Om,OL) ) + 25
    # where dL is the luminosity distance in units of Mpc
    return log10( cosmology.get_dist(z,type='dl',params=[pars[0],pars[1],-1,0])/0.7)*5.0+25.

# Perform the fit, with the iminuit
res = fitting.dothefit(z,mu0,dmu0,[.3,.7],functname=model2,method='minuit')

# res[1] contains the best fitted parameters
# res[2] contains the best fitted standard deviation of the parameters, i.e. the error
# res[3] contains the best fitted covariance matrix of the parameters
# res[4] contains the best chi2 estimate of the fit
# res[5] constains the number of degrees of freedom (ndf) of the fit., ndf = Number of data points - number of fitted parameters

#%matplotlib inline
figure(1,figsize=(15,10)),clf()
errorbar(z,mu0,yerr=dmu0,fmt='.', label='SN Union2.1 Suzuki et Al. 2011') #,label='Riess et Al. 1998 High-Z Templates')
xscale('linear'),ylabel('$\mu(z)$',size=25),xlabel('$z$',size=25)
plot(z,model2(z,res[1]),'r-',label='$\mu(z)=5\log_{10} d_{L}(z)+25$ \n $(\Omega_{m},\Omega_{\Lambda})$=(%0.3f $\pm$ %0.3f,%0.3f $\pm$ %0.3f)'%(res[1][0],res[2][0],res[1][1],res[2][1]) )
legend(loc=4,frameon=False,numpoints=1,fontsize=25)

import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
import matplotlib as mpl
plt.style.use('classic')
import sys
import numpy as np
import os
from scipy.optimize import curve_fit
from matplotlib import cm
import matplotlib.colors as colors
import matplotlib.cm as cm
from functions import count_number_of_files
from functions import return_stuff 

params = {'legend.fontsize': 'medium',
'figure.figsize': (5, 5),
'axes.labelsize': 'x-large',
'axes.titlesize':'x-large',
'xtick.labelsize':'medium',
'ytick.labelsize':'medium'}
pylab.rcParams.update(params)
from matplotlib import rc
rc('font',**{'family':'Times New Roman'})

Mach_number=6
s0=0.65
var=2.*s0
stretch=2.
sigma=var**.5
s_k=1.5*s0

def k_formula(x):
    ll= 1/3*(1/2-1/np.pi*np.arctan((x-s_k)/2))
    return 1./ll
N_sample=200000
N_wait=1000
tau_Eddy=1
t_ff=.7*tau_Eddy
t_ff=0.7# in units of tau eddy
x=np.zeros(N_sample)
#fig,ax=plt.subplots(1,1,figsize=(5.,5.))
fig = plt.figure()
ax = fig.add_axes([0.15,0.15,0.7,.8])


x=np.random.normal(s0,sigma,N_sample)
dt0=1e-3#/tau_Eddy  # timestep of the simulation in T_eddy units
total_time=0.30# in tau_Eddy time unit 
N_steps=int(total_time/dt0)
print(N_steps)
cm_subsection =np.linspace(0, 1, N_steps) 
colors = [cm.tab20c_r(x) for x in cm_subsection]
colors = [cm.jet(x) for x in cm_subsection]
w=np.ones(len(x))
time=dt0 # to avoid division by zero, we set the starting time of the simulation to dt0
#while time<total_time-10*dt0:
def model(x,slope,offset):
    return slope*x+offset
pp=0
for step in range(0,N_steps):
    x_before=np.copy(x)
    n=np.random.normal(0,1,len(x))
    k=k_formula(x)
    D=2*k*var
    dt=dt0
    ind=(1>=(np.exp(-x)*(t_ff/time)**2)).nonzero()[0] #Which cells are affected by self-gravity
    if ind.size==0: #If there are no such cell, continue with a simple Langevin without gravity.
        x=x+n*(D*dt)**.5-(k*(x-s_k))*dt
    else:
        x=x+n*(D*dt)**.5-(k*(x-s_k))*dt #move all cells according to Langevin without gravity
        #Then add the self-gravity to the subset identified in ind=(1> xxx) line above.
        x[ind]+=dt/t_ff* (np.exp(x[ind]))**.5 * 3.*np.pi/2.
    #remove cells whose density is infinity from the list.
    including_ind=((np.logical_not(np.isnan(x))) & (x!=np.inf) & (x!=-np.inf)).nonzero()[0]
    x_before= x_before[including_ind]
    x= x[including_ind]
    values,bins=np.histogram(x,bins=50,normed=True,range=(-5,18))
    nonzeros=(values!=0).nonzero()[0]
    values=values[nonzeros]
    bins=bins[:-1][nonzeros]+np.diff(bins)[0]/2.
    values=values*np.exp(-bins)
    log10_bins=np.log10(np.exp(bins))
    ax.plot(log10_bins,np.log10(values/values.max()),c=colors[pp],lw=1,alpha=1)
    yy=np.log10(values)
    fit_range=((log10_bins>1)&(log10_bins<3)).nonzero()[0]
    xx=log10_bins[fit_range]
    yy=yy[fit_range]
    if yy.size<2:
        line=line=str(time)+str(' ')+str(0)+str("\n")
    else:
        popt, pcov = curve_fit(model, xx,yy)
        line=str(time)+str(' ')+str(popt[0])+str("\n")
    print(line)
    time+=dt
    pp+=1
	
xx=np.logspace(-1,10)
yy=2*xx**-1.695
ax.plot(np.log10(xx),np.log10(yy),c='k')
xx=np.logspace(-3,10)
yy=2e-5*xx**-1
#yy=2e-5*xx**(-var)
ax.plot(np.log10(xx),np.log10(yy),c='k',ls='--',lw=.001)
ax.text(.6,.9,r'$t_{ff,0}=0.7\,\tau_{\rm eddy}$',transform=ax.transAxes,fontsize=13,color='k')

f=np.loadtxt('K11_t_000.dat')
y_numbers=f[:,1]+.7
x_numbers=f[:,0]
ax.plot(x_numbers,y_numbers,ls='dashed',lw=2,label=r'$\rm K11, t=0$',c='k')#color=cm.jet(0.)

f=np.loadtxt('K11_t_018.dat')
y_numbers=f[:,1]+.7
x_numbers=f[:,0]
ax.plot(x_numbers,y_numbers,label=r'$\rm K11, t=0.18\tau_{Eddy}$',ls='dashed',lw=2,c='k')#color=cm.jet(0.18/.3),lw=2)
f=np.loadtxt('K11_t_030.dat')
y_numbers=(np.log10(f[:,2])+1.824)
ax.plot(np.log10(f[:,1]),y_numbers-y_numbers.max(),label=r'$\rm K11, t=0.3\tau_{Eddy}$',ls='--',lw=2,c='k')#,color=cm.jet(1.),lw=2)

ax.set_xlabel(r'$\rm log_{10}(\rho/\rho_{0})$',fontsize=20)
ax.set_ylabel(r'$\rm log_{10}<P_V>$',fontsize=20)
ax.set_xlim(-3,10)
ax.set_ylim(-16,0.1)


plt.legend(loc=3,frameon=False)
plt.savefig('Markov_vs_K11.pdf')

import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
import matplotlib
import sys
import numpy as np
import os
from functions import count_number_of_files
from functions import return_stuff 
from datetime  import datetime
print(datetime.now())
params = {'legend.fontsize': 10,
'figure.figsize': (10, 3),
'axes.labelsize': 'x-large',
'axes.titlesize':'x-large',
'xtick.labelsize':'large',
'ytick.labelsize':'large'}
pylab.rcParams.update(params)
matplotlib.rc('font',family='Times New Roman')
from matplotlib import rc
#rc('font',**{'family':'serif','serif':['Times']})
#rc('font',**{'family':'Times New Roman'})
rc('text', usetex=True)
N_sample=50000
N_wait=100

cc=1/2
def k_formula(tau_eddy,x,cc):
    ll= tau_eddy/3*(cc-1/np.pi*np.arctan((x-s_k)/2))
    return 1./ll
##############################################MAKING PDF PLOTS############################

#based on exact Mach numbers and average dt in the files, for Mach=2,3 and 6 in order
all_frac_lists=[
        [0.0959019253325,0.191434335434,0.382800291933,0.766521946909],
        [0.105386869842,0.210330581051,0.41973860318,0.8388270774],
        [0.098192528125,0.19664193314,0.393473397727,0.785058613636]
        ]
#for plots 
y_min_list=[-.5,-3,-3]
y_max_list=[.5,4,5]
x_min_list=[-2,-3,-4]
x_max_list=[2,3,4]
Mach_numbers=[2,3,6]
s0s=[.27,.52,.88] #s0 values from the MW PDFs of the simulations
stretch=2.
k_p=100. # never used
# loop over A_s : for A_s in [1.5,2] or for A_s in np.linspace(1,2,20) which is 20 numbers between 1 and 2
#for A_s in [3/2.,5./3.,7./4.,2.]:
#for k_n in [0.,0.5,1.]:
#for Mach_6_kp in [200,250,300]:
for A_s in [1.5]:
    for Mach_6_kn in [37.5]:
        kns=[Mach_6_kn*2.03/6.25,Mach_6_kn*3.17/6.25,Mach_6_kn]
        tau_eddys=[.5/2.03,.5/3.17,.5/6.25]
        fig,(ax1,ax2,ax3)=plt.subplots(1,3,sharex=True, sharey=True)
        pp=0
        for ax in [ax1,ax2,ax3]:
            Mach_number=Mach_numbers[pp]
            k_n=kns[pp]
            s0=s0s[pp]
            s_k=A_s*s0
            s_pdf=s_k
            sigma=(2*s0)**.5
            var=2*s0
            tau_eddy=tau_eddys[pp]
            x_initial=np.random.normal(s0,sigma,N_sample)
            values,bins=np.histogram(x_initial,bins=15)
            x_bins=bins[:-1]+np.diff(bins)[0]/2.
            dt=1e-3*tau_eddy
            Nsteps=int(1*tau_eddy/dt)
            for trials in range(0,1):
                print (trials)
                x=np.ones(N_sample)*s0
                x_initial=np.ones(N_sample)*s0
                values,bins=np.histogram(x,bins=201,range=(-5,5))
                x_bins=bins[:-1]+np.diff(bins)[0]
                save_pdf=[]
                ww=1
                for steps in range(0,int(Nsteps)+N_wait):
                    x_begin=np.copy(x)
                    n=np.random.normal(0,1,len(x))
                    #k=k_formula(k_n,k_p,x)
                    k=k_formula(tau_eddy,x,cc)
                    D=2*k*var
                    x=x+n*(D*dt)**.5-(k*(x-s_pdf))*dt
                    values,bins=np.histogram(x,bins=201,range=(-5,5),normed=True)
                    save_pdf.append(values/values.max())
                save_pdf=np.array(save_pdf)#*np.exp(x_bins)
                x_bins=bins[:-1]+np.diff(bins)[0]/2.
                mean_pdf=np.mean(save_pdf[N_wait:,:],0)
            #    ax.plot(x_bins,mean_pdf/mean_pdf.max(),c='r',lw=2)
            xx=np.linspace(-5,5,201)
            save_pdf=[]
            Nfiles=count_number_of_files(Mach_number,1)
            print('Nfiles',Nfiles)
            dir=return_stuff(Mach_number,mode='Fiducial')[0]
            for i in range(24,Nfiles):
                number = str(i).zfill(4)
                filename='s_s2_1_'+number+'.dat'
                f=np.loadtxt('../../Data/'+dir+'/'+filename)
                VW_pdf=f[7:7+201]
                MW_pdf=np.exp(xx)*VW_pdf
                save_pdf.append(MW_pdf/MW_pdf.max())
            save_pdf=np.array(save_pdf)
            mean_pdf=np.mean(save_pdf,0)
            ax.plot(x_bins,mean_pdf,c='k',lw=1)
            upper=np.percentile(save_pdf,84,axis=0)
            lower=np.percentile(save_pdf,16,axis=0)
            upper=np.percentile(save_pdf,97,axis=0)
            lower=np.percentile(save_pdf,2,axis=0)
            ax.fill_between(x_bins,lower,upper,facecolor='grey',interpolate=True,alpha=.7,edgecolor='none')

            x_initial=np.random.normal(s0,sigma,int(1e7))
            values,bins=np.histogram(x_initial,bins=201,range=(-5,5))
            x_bins=bins[:-1]+np.diff(bins)[0]/2.

            ax.axvline(s0,color='k',linestyle='--')

            ax.set_xlabel(r'$s$',fontsize=15)
            if ax==ax1:
                ax.set_ylabel(r'$\rm <P_M>$',fontsize=15)
            ax.text(.7,.8,r'$\mathcal{M}=%d$'%Mach_number,transform=ax.transAxes,fontsize=15,color='k')
            plt.tight_layout()
            ax.set_yscale('linear')
            ax.set_ylim(0,1.05)
            ax.set_xlim(-3,5)
            pp+=1
        fig.subplots_adjust(hspace=0.00,wspace=0.0)
        plt.setp([a.get_yticklabels() for a in fig.axes[1::]], visible=False)
        #plt.setp([a.get_yticklabels() for a in fig.axes], visible=True)
        ax.set_yticklabels(['',0.2,0.4,0.6,0.8,1])
        plt.savefig('row_plot_PDFs_only_sim.pdf')

##########################################PLOTTING ds-s #####################################################

        N_sample=200000
        fig,(ax1,ax2,ax3)=plt.subplots(1,3,sharex=True, sharey=True)
        pp=0
        for ax in [ax1,ax2,ax3]:
            Mach_number=Mach_numbers[pp]
            k_n=kns[pp]
            s0=s0s[pp]
            frac_list=all_frac_lists[pp]
            s_k=A_s*s0
            s_pdf=s_k
            sigma=(2*s0)**.5
            var=2*s0
            tau_eddy=tau_eddys[pp]
            dt=1e-2*tau_eddy
            x_initial=np.random.normal(s0,sigma,N_sample)
            values,bins=np.histogram(x_initial,bins=15)
            x_bins=bins[:-1]+np.diff(bins)[0]/2.
            x=x_initial
            #this portion is to run MCMC from initial Gaussian distribution for a short time####
            for steps in range(0,20):
                dt_large=(1e-4*tau_eddy)
                n=np.random.normal(0,1,len(x))
                #k=k_formula(k_n,k_p,x)
                k=k_formula(tau_eddy,x,cc)
                D=2*k*var
                x=x+n*(D*dt_large)**.5-(k*(x-s_k))*dt_large
            c_list=['b','r','g','k']#color list for each line corresponding to specific Mach number
            for color,frac in zip(c_list,frac_list):
                Nsteps=frac*tau_eddy/dt 
                print(frac,Nsteps)
                save_x=[]
                for steps in range(0,int(Nsteps)):
                    n=np.random.normal(0,1,len(x))
                    #k=k_formula(k_n,k_p,x)
                    k=k_formula(tau_eddy,x,cc)
                    D=(2*k*var)
                    x=x+n*(D*dt)**.5-(k*(x-s_pdf))*dt
                    save_x.append(x)#saves each step's pdf
                save_x=np.array(save_x)#turning the list into an array
                s=np.linspace(-5,5,201)
                s_old=save_x[0,:]#the first pdf we started with
                s_now=save_x[-1,:]#the last pdf we ended up with
                ds=np.zeros(201)
                dsbin=0.05
                for i,s_bin in enumerate(s):
                    index=((s_now>s_bin)&(s_now<s_bin+0.05)).nonzero()[0]
                    if len(index)<1:continue
                    mean_s_before=np.mean(s_old[index])
                    ds[i]=s_bin+dsbin/2.-mean_s_before

            ###
            #    ax.plot(s,ds,c=color)#,label=r't=$%.2f \tau_e$'%frac,c=color)
            ###

            #plotting the ds-s from the simulation data
            nstart_list=[25,12,6,3]
            N_files=count_number_of_files(Mach_number,[1,2,4,8])
            A_list=['1','2','4','8']
            c_list=['b','r','g','k']
            dir='../../Data/'+return_stuff(Mach_number,mode='Fiducial')[0]+'/'
            for Nfiles,nstart,A,color,frac in zip(N_files,nstart_list,A_list,c_list,frac_list):
                S1_sum=0
                DeltaS1_sum_eul=0
                DeltaS1_sum_lag=0
                x=np.linspace(-5,5,201)
                deltat=[]
                for i in range(nstart,int(Nfiles)):
                    number = str(i).zfill(4)
                    filename='s_s2_'+A+'_'+number+'.dat'
                    if os.path.exists(dir+filename)==False:
                        print('No files'),sys.exit()
                    f=np.loadtxt(dir+filename)
                    deltat.append(f[2])
                    S1=f[7:7+201]
                    S1_sum+=S1
                    x_S1max=x[(S1==S1.max()).nonzero()[0]]
                    DeltaS1_sum_lag+=f[7+201*8:7+201*9]
                t_frac=np.mean(np.array(deltat))/tau_eddy
                print(t_frac)
                ax.plot(x,(DeltaS1_sum_lag/S1_sum),c=color,\
                        label=r'$\Delta t=%.1f\tau_{\rm eddy}$'%frac,lw=1,ls='solid')

            #End of plotting the ds-s from the simulation data

            ax.text(.7,.1,r'$\mathcal{M}=%d$'%Mach_number,transform=ax.transAxes,fontsize=15,color='k')
            y_min=y_min_list[pp]
            y_max=y_max_list[pp]
            x_min=x_min_list[pp]
            x_max=x_max_list[pp]
            ax.set_xlim(x_min,x_max)
            ax.set_ylim(y_min,y_max)
            ax.set_xlim(-2.5,3.)
            ax.set_ylim(-3.,3.5)
            if ax==ax1:
                ax.legend(loc=2,frameon=False)
            if ax==ax1:
                ax.set_ylabel(r'$\Delta s$',fontsize=15)
            #if ax==ax3:
                #ax.tick_params(axis='y', which='both', labelleft='off', labelright='off')
                #ax.spines["left"].axis.axes.tick_params(direction="out")
                #ax.spines["right"].axis.axes.tick_params(direction="in")
                #ax.set_yticks([-2,-1,0,1,2])
            ax.set_xlabel(r'$s$',fontsize=15)
            pp+=1
        ax.set_yticks([-2,-1,0,1,2])
        ax.set_xticklabels(['',-2,-1,0,1,2,''])
        #ax.tick_params(labeltop=True, labelright=True)
        plt.tight_layout()
        fig.subplots_adjust(hspace=0,wspace=.0)
        #fig.savefig('April_19_row_final_ds_s_As_%.2f'%A_s+'_kn_%.2f'%k_n+'_kp_%.2f'%Mach_6_kn+'.pdf')
        fig.savefig('row_plot_ds_s_only_sim.pdf')


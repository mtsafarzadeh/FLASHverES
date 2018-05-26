import numpy as np
import os

n_sigma=2.
def return_stuff(Mach_number,mode='Fiducial'):
    if Mach_number==2 :
        dir='Mach2'
        Mid_point=106
        fit_low_bnd=Mid_point#-n_sigma*sigma_s/0.05+30 #0.05 is spacing in x
        fit_up_bnd=176#Mid_point+n_sigma*sigma_s/0.05+15
    if Mach_number==3 :
        if mode=='Fiducial':
            dir='Mach3'
            fit_low_bnd=120
        if mode=='short':
            dir='Mach3short'
            fit_low_bnd=120
        if mode=='veryshort':
            dir='Mach3VeryShort'
            fit_low_bnd=120
        if mode=='LowRes':
            dir='Mach3_256'
            fit_low_bnd=120
        if mode=='MoreModes':
            dir='Mach3_moremodes'
            fit_low_bnd=120
        Mid_point=120
        fit_up_bnd=200#Mid_point+n_sigma*sigma_s/0.05+15
    if Mach_number==5 :
        dir='Mach5'
        Mid_point=116
        fit_low_bnd=Mid_point#-n_sigma*sigma_s/0.05+80
        fit_up_bnd=200#Mid_point+n_sigma*sigma_s/0.05+60
    if Mach_number==4:
        dir='Mach4'
        Mid_point=113
        fit_low_bnd=Mid_point#-n_sigma*sigma_s/0.05+40
        fit_up_bnd=200#Mid_point+n_sigma*sigma_s/0.05+20
    if Mach_number==6:
        dir='Mach6'
        Mid_point=130
        fit_low_bnd=Mid_point#-n_sigma*sigma_s/0.05+40
        fit_up_bnd=200#Mid_point+n_sigma*sigma_s/0.05+20
    return dir,fit_low_bnd,fit_up_bnd

def count_number_of_files(Mach_number,A_list,mode='Fiducial'):
	dir=return_stuff(Mach_number,mode)[0]
	A_list=np.array(A_list)
	if A_list.size>1:
		all_counts=np.zeros(A_list.size)
		pq=0
		for A in A_list: 
			count=0
			for i in range(0,1000):
				number = str(i).zfill(4)
				filename='s_s2_'+str(A)+'_'+number+'.dat'
				if os.path.exists('../../Data/'+dir+'/'+filename)==True: count+=1
				if os.path.exists('../../Data/'+dir+'/'+filename)==False: break
			all_counts[pq]=count
			pq+=1
		return (all_counts)
	else:
		count=0
		for i in range(0,1000):
			number = str(i).zfill(4)
			filename='s_s2_'+str(A_list)+'_'+number+'.dat'
			if os.path.exists('../../Data/'+dir+'/'+filename)==True: count+=1
			if os.path.exists('../../Data/'+dir+'/'+filename)==False: break
		return (count)


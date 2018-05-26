rho = loaddata('noh_mhd_2d_hdf5_chk_0024', 'density',DOUBLE=1, XCOORDS=xc, YCOORDS=yc, TIME=tc)
bx = loaddata('noh_mhd_2d_hdf5_chk_0024', 'magx',DOUBLE=1, XCOORDS=xc, YCOORDS=yc, TIME=tc)
by = loaddata('noh_mhd_2d_hdf5_chk_0024', 'magy',DOUBLE=1, XCOORDS=xc, YCOORDS=yc, TIME=tc)
babs = sqrt(bx*bx+by*by)
bz = loaddata('noh_mhd_2d_hdf5_chk_0024', 'magz',DOUBLE=1, XCOORDS=xc, YCOORDS=yc, TIME=tc)
pr = loaddata('noh_mhd_2d_hdf5_chk_0024', 'pres',DOUBLE=1, XCOORDS=xc, YCOORDS=yc, TIME=tc)
vx = loaddata('noh_mhd_2d_hdf5_chk_0024', 'velx',DOUBLE=1, XCOORDS=xc, YCOORDS=yc, TIME=tc)
vy = loaddata('noh_mhd_2d_hdf5_chk_0024', 'vely',DOUBLE=1, XCOORDS=xc, YCOORDS=yc, TIME=tc)
vabs = sqrt(vx*vx+vy*vy)

window,0
plot,xc(0:*),rho(0:*,5)*1.0e-5,psym=1,symsize=0.5,chars=2,xrange=[0.0,0.6],XTITLE='Radius',YTITLE='Density'

window,1
plot,xc(0:*),vx(0:*,5)*1.0e7,psym=1,symsize=0.5,chars=2,xrange=[0.0,0.6],XTITLE='Radius',YTITLE='Velocity'

window,2
fact = sqrt(1.0e9*4.0*3.1415926)
plot,xc(0:*),bz(0:*,5)*fact,psym=1,symsize=0.5,chars=2,xrange=[0.0,0.6],XTITLE='Radius',YTITLE='Bphi'

window,3
fact = 1.0e9
plot,xc(0:*),pr(0:*,5)*fact,psym=1,symsize=0.5,chars=2,xrange=[0.0,0.6],XTITLE='Radius',YTITLE='Pressure'

window,4
fact = sqrt(1.0e9*4.0*3.1415926)
plot,xc(0:*),abs(by(0:*,5))*fact,psym=1,symsize=0.5,chars=2,xrange=[0.0,0.6],XTITLE='Radius',YTITLE='Bz'

























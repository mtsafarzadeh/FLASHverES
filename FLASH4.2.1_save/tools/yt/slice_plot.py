from yt.mods import *
pf = load("sedov_hdf5_chk_0003") # 3D Sedov simulation
ctr = [0.5,0.5,0.5] # is the actual domain center, can also use pf.domain_center
slc = SlicePlot(pf, 2, ["Density","Temperature"], center=ctr, origin='center-window')
slc.save("raw") # saves the raw plots without changing anything
slc.annotate_velocity() # Add velocity vectors
slc.save("vels")
slc.annotate_grids() # Add block boundaries
slc.save("blocks")
slc.set_width(0.5, "unitary") # sets the domain width to half the
slc.save("half")

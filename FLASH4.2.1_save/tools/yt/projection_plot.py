from yt.mods import *
pf = load("sloshing_low_res_hdf5_plt_cnt_0400") # sloshing simulation
ctr = [0.0,0.0,0.0] # is the actual domain center, can also use pf.domain_center
projx = ProjectionPlot(pf, 0, ["Temperature"], weight_field="Density", center=ctr,
                       width=500./pf["kpc"])
projx.set_log("Temperature", False) # Un-log the temperature field.
projx.save("proj_x") # saves the plot
projy = ProjectionPlot(pf, 1, ["Temperature"], weight_field="Density", center=ctr,
                       width=500./pf["kpc"])
projy.set_log("Temperature", False)
projy.save("proj_y")
projz = ProjectionPlot(pf, 2, ["Temperature"], weight_field="Density", center=ctr,
                       width=500./pf["kpc"])
projz.set_log("Temperature", False)
projz.save("proj_z")

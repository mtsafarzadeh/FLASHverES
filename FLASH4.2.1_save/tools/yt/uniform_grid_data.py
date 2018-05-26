from yt.mods import *
import h5py

# The DustCollapse problem was run with lrefine_max = 4.
pf = load("dust_collapse_hdf5_chk_0003")
max_level = pf.h.max_level # Determine the maximum refinement level from the hierarchy
le = [3.5e8,3.5e8,3.5e8] # Set the left edge of the region 
dims = [32, 32, 32] # Set the number of cells along each direction

# The right edge of the region is then re = le + (nx*dx,ny*dy,nz*dz), where the 
# (dx, dy, dz) are the cell widths at refinement level max_level. We don't need
# to compute this, however. 

cube = pf.h.covering_grid(max_level, le, dims, fields=["dens","temp","pres"])

f = h5py.File("%s_uniform_grid.h5" % pf, "w") # A writeable HDF5 file

# Write the datasets to the file
f.create_dataset("/dens", data=cube["dens"])
f.create_dataset("/temp", data=cube["temp"])
f.create_dataset("/pres", data=cube["pres"])

f.close() # close the file

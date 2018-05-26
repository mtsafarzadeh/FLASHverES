from yt.mods import *

pf = load("sloshing_low_res_hdf5_plt_cnt_0400")

# Define a sphere object with a center at the domain center and a radius of 100 kpc.
# Note that we can specify a tuple (radius, unit) instead of calculating the radius
# in cgs units
sp = pf.h.sphere(pf.domain_center, (100.0, "kpc"))

# Compute the mass-weighted temperature in the sphere
t_mw = sp.quantities["WeightedAverageQuantity"]("Temperature", "CellMass")

# Compute the total gas mass in the sphere
m_gas = sp.quantities["TotalQuantity"]("CellMassMsun")

# Compute the angular momentum vector of the sphere
L_vec = sp.quantities["AngularMomentumVector"]()

print t_mw, m_gas, L_vec

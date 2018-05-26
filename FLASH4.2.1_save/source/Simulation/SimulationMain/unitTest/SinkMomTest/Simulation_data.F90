

module Simulation_data
  
  implicit none
  
  integer,save  :: sim_globalMe
  
  real, save    :: sim_xMin, sim_xMax, sim_yMin, sim_yMax, sim_zMin, sim_zMax
  real, save    :: sim_xcenter, sim_ycenter, sim_zcenter, sim_radius, sim_dens, sim_cs
  real, save    :: sim_vx, sim_vy, sim_vz
  real, save    :: sim_sink_x, sim_sink_y, sim_sink_z, sim_sink_vx, sim_sink_vy, sim_sink_vz, sim_sink_mass
  
end module Simulation_data

subroutine sim_xToXi(x, t, n, xi)
  use Simulation_data, ONLY: sim_xi0, sim_alpha, sim_Q, sim_toffset, &
       sim_condTemperatureExponent
  implicit none
  real, intent(in) :: x, t, n
  real, intent(out) :: xi

!!$  n = sim_condTemperatureExponent

  xi = x / ( (sim_alpha * sim_Q**n * (t+sim_toffset))**(1.0/(n+2)) )
end subroutine sim_xToXi

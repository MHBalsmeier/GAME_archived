! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/MHBalsmeier/game

module radiation
  
  use iso_c_binding,
  use mo_gas_optics_rrtmgp, only: ty_gas_optics_rrtmgp
  use mo_rte_kind,          only: wp
  
  implicit none
  
  type(ty_gas_optics_rrtmgp) :: k_dist_sw, k_dist_lw

  character(len=3), dimension(8) :: active_gases = (/ &
      'H2O', 'CO2', 'O3 ', 'N2O', &
      'CO ', 'CH4', 'O2 ', 'N2 ' &
   /)
  
  character(len=128) :: rrtmgp_coefficients_file_sw, rrtmgp_coefficients_file_lw
  
  contains
  
  subroutine calc_radiative_flux_convergence(mass_densities, temperature_gas, radiation_tendency, &
  no_of_scalars, no_of_constituents) &
  bind(c, name = "calc_radiative_flux_convergence")
    
    integer, intent(in)              ::                    no_of_scalars
    integer, intent(in)              ::                    no_of_constituents
    real(8), intent(in)              :: mass_densities    &
    (no_of_constituents*no_of_scalars)
    real(8), intent(in)              :: temperature_gas   (no_of_scalars)
    real(8), intent(inout)           :: radiation_tendency(no_of_scalars)
    
    ! local variables
    integer                          :: ji
    
    do ji=1,no_of_scalars
      radiation_tendency(ji) = 0.000
    enddo
    
  end subroutine calc_radiative_flux_convergence

end module radiation

! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/MHBalsmeier/game

module radiation
  
  use iso_c_binding,
  use mo_rte_kind,           only: wp
  use mo_gas_optics_rrtmgp,  only: ty_gas_optics_rrtmgp
  use mo_rte_sw,             only: rte_sw
  use mo_rte_lw,             only: rte_lw
  
  implicit none
  
  type(ty_gas_optics_rrtmgp) :: k_dist_sw, k_dist_lw

  character(len=3), dimension(8) :: active_gases = (/ &
      'H2O', 'CO2', 'O3 ', 'N2O', &
      'CO ', 'CH4', 'O2 ', 'N2 ' &
   /)
  
  character(len=128) :: rrtmgp_coefficients_file_sw, rrtmgp_coefficients_file_lw
  
  contains
  
  subroutine calc_radiative_flux_convergence(latitude_scalar, longitude_scalar, &
  mass_densities, temperature_gas, radiation_tendency, &
  no_of_scalars, no_of_layers, no_of_constituents) &
  bind(c, name = "calc_radiative_flux_convergence")
    
    integer, intent(in)              ::                    no_of_scalars
    integer, intent(in)              ::                    no_of_layers
    integer, intent(in)              ::                    no_of_constituents
    real(8), intent(in)              :: latitude_scalar    (no_of_scalars/no_of_layers)
    real(8), intent(in)              :: longitude_scalar   (no_of_scalars/no_of_layers)
    real(8), intent(in)              :: mass_densities    &
    (no_of_constituents*no_of_scalars)
    real(8), intent(in)              :: temperature_gas   (no_of_scalars)
    real(8), intent(inout)           :: radiation_tendency(no_of_scalars)
    
    ! local variables
    ! solar zenith angle
    real(8)	                         :: mu_0(no_of_scalars/no_of_layers)
    integer                          :: ji, no_of_scalars_h
    real(8)                          :: time
    
    ! the number of scalars on every layer
    no_of_scalars_h = no_of_scalars/no_of_layers
    
    ! calculating the zenith angle
    do ji=1,no_of_scalars_h
      mu_0(ji) = zenith(latitude_scalar(ji), longitude_scalar(ji), time)
    enddo
    
    ! the result
    do ji=1,no_of_scalars
      radiation_tendency(ji) = 0.000
    enddo
    
  end subroutine calc_radiative_flux_convergence
  
  real(8) function zenith(lat, lon, t)
  
  ! calculates the cosine of the zenith angle at a given
  ! point and time
  
    real(8), intent(in)              :: lat
    real(8), intent(in)              :: lon
    real(8), intent(in)              :: t
  
    zenith = 0.
  
  end function zenith

end module radiation






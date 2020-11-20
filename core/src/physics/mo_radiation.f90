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
  no_of_scalars, no_of_layers, no_of_constituents, time_coord) &
  bind(c, name = "calc_radiative_flux_convergence")
    
    integer, intent(in)              ::                    no_of_scalars
    integer, intent(in)              ::                    no_of_layers
    integer, intent(in)              ::                    no_of_constituents
    real(8)                          :: time_coord
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
    
    ! the number of scalars on every layer
    no_of_scalars_h = no_of_scalars/no_of_layers
    
    ! calculating the zenith angle
    do ji=1,no_of_scalars_h
      mu_0(ji) = coszenith(latitude_scalar(ji), longitude_scalar(ji), time_coord)
    enddo
    
    ! the result
    do ji=1,no_of_scalars
      radiation_tendency(ji) = 0.000
    enddo
    
  end subroutine calc_radiative_flux_convergence
  
  real(8) function coszenith(lat, lon, t)
  
    ! calculates the cosine of the zenith angle at a given
    ! point and time
  
  	! the coordinates of the place we look at
    real(8), intent(in)              :: lat
    real(8), intent(in)              :: lon
    ! the unix time stamp of the time
    real(8), intent(in)              :: t
    
    ! local variables
    real(8)                          :: normal_vector_rel2_earth(3)
    real(8)                          :: normal_vector_rel2_sun  (3)
    real(8)                          :: sun_2_earth             (3)
    ! obliquity of the earth's axis
    real(8)                          :: obliquity
    ! rotation speed of the earth
    real(8)                          :: omega
    ! revolution speed of the earth around the sun
    real(8)                          :: omega_rev
    ! a reference time
    real(8)                          :: t_0
    ! a transformed time
    real(8)                          :: t_transformed
    ! the rotation angle of the earth
    real(8)                          :: rot_angle
    ! At the time t_0, the earth has been at an angle phi_0_earth_rotation
    ! around itself and at an angle phi_0_earth_around_sun around the sun.
    real(8)                          :: phi_0_earth_around_sun
    real(8)                          :: phi_0_earth_rotation
    real(8)                          :: trans_earth2sun         (3,3)
    
    omega                  = 7.292115d-5
    omega_rev              = 1.99099d-7
    obliquity              = 0.409092592843564
    
    ! refer to https://www.esrl.noaa.gov/gmd/grad/solcalc/azel.html
    ! Unix time coordinate of 2019-Dec-20, 12:00 UTC
    t_0                    = 1576843200.0
    ! this is a winter solstice
    phi_0_earth_around_sun = 0.
    phi_0_earth_rotation   = 0.
    
    ! transformation of the time coordinate
    t_transformed = t - t_0
    
    rot_angle = omega*t_transformed - phi_0_earth_rotation
    
    ! the normal vector of the place we look at in earth fixed coordinates
    normal_vector_rel2_earth(1) = cos(lat)*cos(lon)
    normal_vector_rel2_earth(2) = cos(lat)*sin(lon)
    normal_vector_rel2_earth(3) = sin(lat)
    
    ! the x vector of the earth fixed coordinate system in solar coordinates
    trans_earth2sun(1,1) = -cos(rot_angle)*cos(obliquity)
    trans_earth2sun(2,1) = -sin(rot_angle)
    trans_earth2sun(3,1) = cos(rot_angle)*sin(obliquity)
    ! the y vector of the earth fixed coordinate system in solar coordinates
    trans_earth2sun(1,2) = sin(rot_angle)*cos(obliquity)
    trans_earth2sun(2,2) = -cos(rot_angle)
    trans_earth2sun(3,2) = -sin(rot_angle)*sin(obliquity)
    ! the z vector of the earth fixed coordinate system in solar coordinates
    trans_earth2sun(1,3) = sin(obliquity)
    trans_earth2sun(2,3) = 0
    trans_earth2sun(3,3) = cos(obliquity)
    
    ! transforming the normal vector of the place to solar coordinates
    normal_vector_rel2_sun      = matmul(trans_earth2sun, normal_vector_rel2_earth)
    
    sun_2_earth             (1) = cos(omega_rev*t_transformed + phi_0_earth_around_sun)
    sun_2_earth             (2) = sin(omega_rev*t_transformed + phi_0_earth_around_sun)
    sun_2_earth             (3) = 0
    
    ! the result
    coszenith = DOT_PRODUCT(normal_vector_rel2_earth, -sun_2_earth)
    
    ! the night case
    if (coszenith < 0) then
      coszenith = 0
    endif
  
  end function coszenith

end module radiation
















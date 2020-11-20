! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/MHBalsmeier/game

module radiation
  
  use iso_c_binding,
  use mo_rte_kind,           only: wp
  use mo_gas_optics_rrtmgp,  only: ty_gas_optics_rrtmgp
  use mo_gas_concentrations, only: ty_gas_concs
  use mo_fluxes_byband,      only: ty_fluxes_byband
  use mo_rte_sw,             only: rte_sw
  use mo_rte_lw,             only: rte_lw
  
  implicit none
  
  ! the number of bands in the short wave region
  integer, parameter                 :: no_of_sw_bands = 14
  ! the number of bands in the long wave region
  integer, parameter                 :: no_of_lw_bands = 16
  ! the number of g points is the short wave region
  integer                            :: no_of_sw_g_points
  ! the number of g points is the long wave region
  integer                            :: no_of_lw_g_points
  
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
    ! the gas concentrations
    type(ty_gas_concs)               :: gas_concentrations
    ! number of scalars per layer (number of columns)
    integer                          :: no_of_scalars_h
    ! number of points where it is day
    integer                          :: no_day_points
    ! number of points where it is night
    integer                          :: no_night_points
    ! loop indices
    integer                          :: ji, j_day, j_night
    ! the indices of columns where it is day
    integer                          :: day_indices(no_of_scalars/no_of_layers)
    ! the indices of columns where it is night
    integer                          :: night_indices(no_of_scalars/no_of_layers)
    ! the resulting clear sky fluxes
    type(ty_fluxes_byband)           :: fluxes_clearsky_day
    real(8)                          :: surface_emissivity(no_of_lw_bands, no_of_scalars/no_of_layers)
    real(8)                          :: albedo_dir        (no_of_lw_bands, no_of_scalars/no_of_layers)
    real(8)                          :: albedo_dif        (no_of_lw_bands, no_of_scalars/no_of_layers)
    
    ! calculation of the number of columns
    no_of_scalars_h = no_of_scalars/no_of_layers
    
    ! set the surface emissivity to one
    surface_emissivity(:,:) = 1.
    
    ! set the surface albedos to 0.5
    albedo_dir        (:,:) = 0.5
    albedo_dif        (:,:) = 0.5
    
    ! calculating the zenith angle, and counting day and night points
    j_day = 0
    j_night = 0
    do ji=1,no_of_scalars_h
      mu_0(ji) = coszenith(latitude_scalar(ji), longitude_scalar(ji), time_coord)
      if (mu_0(ji) > 0) then
        j_day   = j_day + 1
        day_indices(j_day)     = ji
      else
        j_night = j_night + 1
        night_indices(j_night) = ji
      endif
    enddo
    
    no_day_points   = j_day
    no_night_points = j_night
    
    no_of_sw_g_points = k_dist_sw%get_ngpt()
    no_of_lw_g_points = k_dist_lw%get_ngpt()
    
    ! calculate shorwave radiative fluxes
    ! rte_sw(k_dist_sw, gas_concentrations, albedo_dir, albedo_dif, mu_0)
    
    ! calculate longwave radiative fluxes
    ! rte_lw(k_dist_lw, gas_concentrations, surface_emissivity(:,:))
    
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
















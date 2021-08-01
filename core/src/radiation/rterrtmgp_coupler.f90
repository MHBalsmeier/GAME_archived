! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/AUN4GFD/game

! This file calculates radiative fluxes.

module radiation
  
  use iso_c_binding,
  use mo_rte_kind,           only: wp
  use mo_rrtmgp_util_string, only: lower_case
  use mo_gas_optics_rrtmgp,  only: ty_gas_optics_rrtmgp
  use mo_load_coefficients,  only: load_and_init
  use mo_gas_concentrations, only: ty_gas_concs
  use mo_fluxes_byband,      only: ty_fluxes_broadband
  use mo_source_functions,   only: ty_source_func_lw
  use mo_rte_sw,             only: rte_sw
  use mo_rte_lw,             only: rte_lw
  use mo_optical_props,      only: ty_optical_props_1scl, &
                                   ty_optical_props_2str
  
  implicit none
  
  ! the number of bands in the short wave region
  integer, parameter                 :: no_of_sw_bands = 14
  ! the number of bands in the long wave region
  integer, parameter                 :: no_of_lw_bands = 16
  ! the number of g points is the short wave region
  integer                            :: no_of_sw_g_points
  ! the number of g points is the long wave region
  integer                            :: no_of_lw_g_points
  ! used for C interoperability
  integer                            :: zero = 0
  integer                            :: one = 1
  ! the gas concentrations (object holding all information on the composition
  ! of the gas phase)
  type(ty_gas_concs)                 :: gas_concentrations_sw
  type(ty_gas_concs)                 :: gas_concentrations_lw
  
  type(ty_gas_optics_rrtmgp)         :: k_dist_sw, k_dist_lw

  character(len = 3), dimension(8) :: active_gases =  (/ &
   "N2 ", "O2 ", "CH4", "O3 ", "CO2", "H2O", "N2O", "CO " &
   /)
  
  character(len = *), parameter      :: rrtmgp_coefficients_file_sw =  &
  ! insert the name of the short wave data file here
  "/home/max/code/rte-rrtmgp/rrtmgp/data/rrtmgp-data-sw-g224-2018-12-04.nc"
  character(len = *), parameter      :: rrtmgp_coefficients_file_lw =  &
  ! insert the name of the long wave data file here
  "/home/max/code/rte-rrtmgp/rrtmgp/data/rrtmgp-data-lw-g256-2018-12-04.nc"
  ! the gases in lowercase
  character(len = 32), dimension(size(active_gases)) :: gases_lowercase
  
  ! interface to C functions
  interface
    real(C_DOUBLE) function specific_gas_constants(gas_number) bind(c, name = "specific_gas_constants")
      use, intrinsic::iso_c_binding
      implicit none
      integer(C_INT), value :: gas_number
    end function specific_gas_constants
    real(C_DOUBLE) function molar_fraction_in_dry_air(gas_number) bind(c, name = "molar_fraction_in_dry_air")
      use, intrinsic::iso_c_binding
      implicit none
      integer(C_INT), value :: gas_number
    end function molar_fraction_in_dry_air
    real(C_DOUBLE) function calc_o3_vmr(z_height) bind(c, name = "calc_o3_vmr")
      use, intrinsic::iso_c_binding
      implicit none
      real(C_DOUBLE), value :: z_height
    end function calc_o3_vmr
  end interface
  
  contains
  
  subroutine radiation_init() &
  bind(c, name =  "radiation_init")
    
    ! This is called only once, in the beginning.
    
    ! local variables
    ! loop index
    integer                          :: ji
    
    ! getting the number of g points in both spectral regions
    no_of_sw_g_points =  k_dist_sw%get_ngpt()
    no_of_lw_g_points =  k_dist_lw%get_ngpt()
    
    ! formatting the gas names
    do ji =  1,size(active_gases)
      gases_lowercase(ji) =  trim(lower_case(active_gases(ji)))
    end do
    ! here, the names of the gases are written to the gas_concentrations object
    call handle_error(gas_concentrations_sw%init(gases_lowercase))
    call handle_error(gas_concentrations_lw%init(gases_lowercase))
    
    ! loading the short wave radiation properties
    call load_and_init(k_dist_sw, rrtmgp_coefficients_file_sw, gas_concentrations_sw)
    ! loading the long wave radiation properties
    call load_and_init(k_dist_lw, rrtmgp_coefficients_file_lw, gas_concentrations_lw)
    
  end subroutine radiation_init
  
  subroutine calc_radiative_flux_convergence(latitude_scalar, longitude_scalar, &
  z_scalar, z_vector, &
  mass_densities, temperature_gas, radiation_tendency, &
  temp_sfc, sfc_sw_in, sfc_lw_out, &
  no_of_scalars, no_of_layers, no_of_constituents, no_of_condensed_constituents, &
  time_coord) &
  
  ! This is the function that is called by the dynamical core. The dycore hands over
  ! the thermodynamic state as well as meta data (time stamp, coordinates) and gets
  ! back radiative fux convergences in W/m^3.
  
  bind(c, name =  "calc_radiative_flux_convergence")
    
    ! the number of scalar points of the model grid
    integer, intent(in)              ::                    no_of_scalars
    ! the number of layers of the model grid
    integer, intent(in)              ::                    no_of_layers
    ! the number of constituents of the model atmosphere
    integer, intent(in)              ::                    no_of_constituents
    ! the numer of condensed constituents of the model atmosphere
    integer, intent(in)              ::                    no_of_condensed_constituents
    ! the time coordinate (UTC time stamp)
    real(8)                          :: time_coord
    ! the latitude coordinates of the scalar data points
    real(8), intent(in)              :: latitude_scalar    (no_of_scalars/no_of_layers)
    ! the longitude coordinates of the scalar data points
    real(8), intent(in)              :: longitude_scalar   (no_of_scalars/no_of_layers)
    ! the vertical positions of the scalar data points
    real(8), intent(in)              :: z_scalar           (no_of_scalars)
    ! the vertical positions of the vector data points
    real(8), intent(in)              :: z_vector           (no_of_scalars)
    ! the mass densities of the model atmosphere
    real(8), intent(in)              :: mass_densities    &
    (no_of_constituents*no_of_scalars)
    ! the temperature of the model atmosphere
    real(8), intent(in)              :: temperature_gas   (no_of_scalars)
    ! the result (in W/m^3)
    real(8), intent(inout)           :: radiation_tendency(no_of_scalars)
    ! surface temperature
    real(8), intent(in)              :: temp_sfc          (no_of_scalars/no_of_layers)
    ! surface shortwave in
    real(8), intent(inout)           :: sfc_sw_in         (no_of_scalars/no_of_layers)
    ! surface longwave out
    real(8), intent(inout)           :: sfc_lw_out        (no_of_scalars/no_of_layers)
    
    ! local variables
    ! solar zenith angle
    real(8)                          :: mu_0(no_of_scalars/no_of_layers)
    ! number of points where it is day
    integer                          :: no_of_day_points
    ! loop indices
    integer                          :: ji, j_day, jk
    ! the indices of columns where it is day
    integer                          :: day_indices(no_of_scalars/no_of_layers)
    ! number of scalars per layer (number of columns)
    integer                          :: no_of_scalars_h
    ! the resulting fluxes
    type(ty_fluxes_broadband)        :: fluxes, fluxes_day
    ! short wave optical properties
    type(ty_optical_props_2str)      :: optical_props_sw
    ! long wave optical properties
    type(ty_optical_props_1scl)      :: optical_props_lw
    ! top of atmosphere short wave flux
    real(wp), dimension(:,:), allocatable          :: toa_flux ! no_of_day_points, no_of_sw_g_points
    ! long wave source function
    type(ty_source_func_lw)          :: sources_lw
    ! the surface emissivity
    real(8)                          :: surface_emissivity(no_of_lw_bands, no_of_scalars/no_of_layers)
    ! surface albedo for direct radiation
    real(8)                          :: albedo_dir        (no_of_sw_bands, no_of_scalars/no_of_layers)
    ! surface albedo for diffusive radiation
    real(8)                          :: albedo_dif        (no_of_sw_bands, no_of_scalars/no_of_layers)
    ! surface albedo for direct radiation (day points only)
    real(8)                          :: albedo_dir_day    (no_of_sw_bands, no_of_scalars/no_of_layers)
    ! surface albedo for diffusive radiation (day points only)
    real(8)                          :: albedo_dif_day    (no_of_sw_bands, no_of_scalars/no_of_layers)
    ! solar zenith angle (day points only)
    real(8)                          :: mu_0_day(no_of_scalars/no_of_layers)
    ! temperature at the surface (day points only)
    real(8)                          :: temp_sfc_day(no_of_scalars/no_of_layers)
    ! reformatted temperature field
    real(8)                          :: temperature_rad           (no_of_scalars/no_of_layers, no_of_layers)
    ! reformatted pressure field
    real(8)                          :: pressure_rad              (no_of_scalars/no_of_layers, no_of_layers)
    ! pressure at cell interfaces
    real(8)                          :: pressure_interface_rad    (no_of_scalars/no_of_layers, no_of_layers+1)
    ! temperature at cell interfaces
    real(8)                          :: temperature_interface_rad (no_of_scalars/no_of_layers, no_of_layers+1)
    ! temperature at cells restricted to day points
    real(8)                          :: temperature_rad_day       (no_of_scalars/no_of_layers, no_of_layers)
    ! pressure at cells restricted to day points
    real(8)                          :: pressure_rad_day          (no_of_scalars/no_of_layers, no_of_layers)
    ! pressure at cell interfaces restricted to day points
    real(8)                          :: pressure_interface_rad_day(no_of_scalars/no_of_layers, no_of_layers+1)
    ! scale height of the atmosphere
    real(8), parameter               :: scale_height = 8.e3_wp
    
    ! calculation of the number of columns
    no_of_scalars_h =  no_of_scalars/no_of_layers
    
    ! set the surface emissivity (a longwave property) to a standard value
    surface_emissivity(:,:) =  0.98_wp
    
    ! setting the surface albedos to 0.12 (compare Zdunkowski, Trautmann & Bott:
    ! Radiation in the Atmosphere, 2007, p. 444)
    albedo_dir        (:,:) =  0.12_wp
    albedo_dif        (:,:) =  0.12_wp
    
    ! reformatting the thermodynamical state for RTE+RRTMGP
    do ji = 1,no_of_scalars_h
      do jk = 1,no_of_layers
        temperature_rad(ji,jk) = temperature_gas((jk-1)*no_of_scalars_h+ji)
        ! the pressure is diagnozed here, using the equation of state for ideal gases
        pressure_rad(ji,jk) = specific_gas_constants(0) &
        *mass_densities(no_of_condensed_constituents*no_of_scalars &
        + (jk-1)*no_of_scalars_h+ji)*temperature_rad(ji,jk)
      enddo
    enddo
    
    ! moving the temperature into the allowed area
    do ji = 1,no_of_scalars_h
      do jk = 1,no_of_layers
        if (temperature_rad(ji,jk) > k_dist_lw%get_temp_max()) then
          temperature_rad(ji,jk) = k_dist_lw%get_temp_max()
        endif
        if (temperature_rad(ji,jk) < k_dist_lw%get_temp_min()) then
          temperature_rad(ji,jk) = k_dist_lw%get_temp_min()
        endif
        if (temperature_rad(ji,jk) > k_dist_sw%get_temp_max()) then
          temperature_rad(ji,jk) = k_dist_sw%get_temp_max()
        endif
        if (temperature_rad(ji,jk) < k_dist_sw%get_temp_min()) then
          temperature_rad(ji,jk) = k_dist_sw%get_temp_min()
        endif
      enddo
    enddo
    
    ! the properties at cell interfaces
    do ji = 1,no_of_scalars_h
      do jk = 1,no_of_layers+1
        ! values at TOA
        if (jk == 1) then
          ! temperature at TOA (linear extrapolation)
          ! the value in the highest layer
          temperature_interface_rad(ji,jk) = temperature_rad(ji,jk) &
          ! the gradient
          ! delta T
          + (temperature_rad(ji,jk) - temperature_rad(ji,jk+1))/     &
          ! delta z
          (z_scalar(ji+(jk-1)*no_of_scalars_h)-z_scalar(ji+jk*no_of_scalars_h)) &
          ! times delta_z
          *(z_vector(ji)-z_scalar(ji+(jk-1)*no_of_scalars_h))
          ! pressure at TOA
          ! here, the barometric height formula is used
          pressure_interface_rad   (ji,jk) =  pressure_rad   (ji,jk) &
          *EXP(-(z_vector(ji)-z_scalar(ji+(jk-1)*no_of_scalars_h))/scale_height)
        ! values at the surface
        elseif (jk == no_of_layers+1) then
          ! temperature at the surfac
          ! the value in the lowest layer
          temperature_interface_rad(ji,jk) = temp_sfc(ji)
          ! surface pressure
          pressure_interface_rad   (ji,jk) =  pressure_rad   (ji,jk-1) &
          *EXP(-(z_vector(no_of_layers*no_of_scalars_h+ji) &
          -z_scalar(ji+(jk-2)*no_of_scalars_h))/scale_height)
        else
          ! just the arithmetic mean
          temperature_interface_rad(ji,jk) =  0.5_wp*(temperature_rad(ji,jk-1)+temperature_rad(ji,jk))
          pressure_interface_rad   (ji,jk) =  0.5_wp*(pressure_rad   (ji,jk-1)+pressure_rad   (ji,jk))
        endif
      enddo
    enddo
    
    ! calculating the zenith angle, and counting day and night points
    j_day =  0
    do ji = 1,no_of_scalars_h
      mu_0(ji) =  coszenith(latitude_scalar(ji), longitude_scalar(ji), time_coord)
      ! it should be > 0, but this would lead to problems with the slicing procedure
      if (mu_0(ji) >= 0) then
        j_day  = j_day + 1
        day_indices(j_day)    = ji
      endif
    enddo
    
    no_of_day_points = j_day
    
    ! now we start the actual radiation calculation
    ! clearing the radiation tendency (it still contains the results of the previous call
    ! from the dycore)
    do ji = 1,no_of_scalars
      radiation_tendency(ji) = 0._wp
    enddo
    
    ! short wave first
    ! filling up the arrays restricted to day points
    do j_day = 1,no_of_day_points
      temperature_rad_day(j_day,:)       = temperature_rad(day_indices(j_day),:)
      pressure_rad_day(j_day,:)          = pressure_rad(day_indices(j_day),:)
      pressure_interface_rad_day(j_day,:)= pressure_interface_rad(day_indices(j_day),:)
      mu_0_day(j_day)                    = mu_0(day_indices(j_day))
      temp_sfc_day(j_day)                = temp_sfc(day_indices(j_day))
      albedo_dir_day(:,j_day)            = albedo_dir(:,day_indices(j_day))  
      albedo_dif_day(:,j_day)            = albedo_dif(:,day_indices(j_day))   
    end do
    
    ! setting the volume mixing ratios of the gases for the short wave calculation
    call set_vol_mix_ratios(mass_densities, .true., no_of_day_points, no_of_scalars_h, &
    no_of_layers, no_of_scalars, no_of_condensed_constituents, day_indices, z_scalar)
    
    ! initializing the short wave fluxes
    call init_fluxes(fluxes_day, no_of_day_points, no_of_layers+1, no_of_sw_bands)
    
    ! allocating the short wave optical properties
    call handle_error(optical_props_sw%alloc_2str(no_of_day_points, no_of_layers, k_dist_sw))
    
    ! allocating the TOA flux
    allocate(toa_flux(no_of_day_points, k_dist_sw%get_ngpt()))
    
    ! setting the short wave optical properties
    call handle_error(k_dist_sw%gas_optics(pressure_rad_day(1:no_of_day_points,:),           &
                                           pressure_interface_rad_day(1:no_of_day_points,:), &
                                           temperature_rad_day(1:no_of_day_points,:),        &
                                           gas_concentrations_sw,                            &
                                           optical_props_sw,                                 &
                                           toa_flux))
    
    ! calculate shortwave radiative fluxes (only the day points are handed over
    ! for efficiency)
    call handle_error(rte_sw(optical_props_sw,                     &
                             .true.,                               &
                             mu_0_day(1:no_of_day_points),         &
                             toa_flux,                             &
                             albedo_dir_day(:,1:no_of_day_points), &
                             albedo_dif_day(:,1:no_of_day_points), &
                             fluxes_day))
    
    ! short wave result (in Wm^-3)
    call calc_power_density(.true., no_of_scalars, &
    no_of_layers, no_of_scalars_h, no_of_day_points, day_indices, &
    fluxes_day, z_vector, radiation_tendency)
    
    ! saving the surface shortwave inward radiative flux density
    do ji = 1,no_of_day_points
      sfc_sw_in(day_indices(ji)) = fluxes_day%flux_dn(ji,no_of_layers+1) &
                                 - fluxes_day%flux_up(ji,no_of_layers+1)
    enddo
    
    ! freeing the short wave fluxes
    call free_fluxes(fluxes_day)
    
    ! now long wave
    ! setting the volume mixing ratios of the gases for the long wave calculation
    call set_vol_mix_ratios(mass_densities, .false., no_of_day_points, no_of_scalars_h, &
    no_of_layers, no_of_scalars, no_of_condensed_constituents, day_indices, z_scalar)
    
    ! initializing the long wave fluxes
    call init_fluxes(fluxes, no_of_scalars_h, no_of_layers+1, no_of_lw_bands)
    
    ! allocating the long wave optical properties
    call handle_error(optical_props_lw%alloc_1scl(no_of_scalars_h, no_of_layers, k_dist_lw))
    
    ! allocating the long wave source function
    call handle_error(sources_lw%alloc(no_of_scalars_h, no_of_layers, k_dist_lw))
    
    ! setting the long wave optical properties
    call handle_error(k_dist_lw%gas_optics(pressure_rad,                      &
                                           pressure_interface_rad,            &
                                           temperature_rad,                   &
                                           temp_sfc,                          &
                                           gas_concentrations_lw,             &
                                           optical_props_lw,                  &
                                           sources_lw,                        &
                                           tlev = temperature_interface_rad))
    
    ! calculate longwave radiative fluxes
    call handle_error(rte_lw(optical_props_lw,   &
                             .true.,             &
                             sources_lw,         &
                             surface_emissivity, &
                             fluxes))
   
    ! add long wave result (in Wm^-3)
    call calc_power_density(.false., no_of_scalars,               &
    no_of_layers, no_of_scalars_h, no_of_day_points, day_indices, &
    fluxes, z_vector, radiation_tendency)
    
    ! saving the surface longwave outward radiative flux density
    do ji = 1,no_of_scalars_h
      sfc_lw_out(ji) = fluxes%flux_up(ji,no_of_layers+1) &
                     - fluxes%flux_dn(ji,no_of_layers+1)
    enddo
    
    ! freeing the long wave fluxes
    call free_fluxes(fluxes)
    
  end subroutine calc_radiative_flux_convergence
    
  subroutine calc_power_density(day_only, no_of_scalars,        &
  no_of_layers, no_of_scalars_h, no_of_day_points, day_indices, &
  fluxes, z_vector, radiation_tendency)
  
    ! this is essentially the negative vertical divergence operator
    
    ! true for short wave calculations (for efficiency)
    logical, intent(in)              :: day_only
    ! as usual
    integer, intent(in)              :: no_of_scalars
    ! as usual
    integer, intent(in)              :: no_of_layers
    ! as usual
    integer, intent(in)              :: no_of_scalars_h
    ! as usual
    integer, intent(in)              :: no_of_day_points
    ! the indices of the columns where it is day
    integer, intent(in)              :: day_indices(no_of_scalars_h)
    type(ty_fluxes_broadband), intent(in):: fluxes
    ! as usual
    real(8), intent(in)              :: z_vector(no_of_scalars + no_of_scalars_h)
    ! the result (in W/m^3)
    real(8), intent(inout)           :: radiation_tendency(no_of_scalars)
  
    ! local variables
    ! the layer index
    integer                          :: ji
    ! the index of the relevant column
    integer                          :: j_column
    ! the horizontal index
    integer                          :: jk
    ! the number of columns taken into account
    integer                          :: no_of_relevant_columns
    
    if (day_only) then
      no_of_relevant_columns = no_of_day_points
    else
      no_of_relevant_columns = no_of_scalars_h
    endif
  
    ! loop over all columns
    do j_column = 1,no_of_relevant_columns
      ! loop over all layers
      do ji = 1,no_of_layers
        ! finding the relevant horizontal index
        if (day_only) then
          jk = day_indices(j_column)
        else
          jk = j_column
        endif
        radiation_tendency((ji-1)*no_of_scalars_h+jk) = &
        ! this function is called four times, therefore we need to
        ! add up the tendencies
        radiation_tendency((ji-1)*no_of_scalars_h+jk) + &
        ! this is a sum of four fluxes
        ( &
        ! upward flux (going in)
        fluxes%flux_up  (j_column,ji+1)                 &
        ! upward flux (going out)
        - fluxes%flux_up(j_column,ji  )                 &
        ! downward flux (going in)
        + fluxes%flux_dn(j_column,ji)                   &
        ! downward flux (going out)
        - fluxes%flux_dn(j_column,ji+1))                &
        ! dividing by the column thickness (the shallow atmosphere
        ! approximation is made at this point)
        /(z_vector((ji-1)*no_of_scalars_h+jk) - z_vector(ji*no_of_scalars_h+jk))
      enddo
    enddo
  
  end subroutine calc_power_density
  
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
    
    omega                 = 7.292115e-5_wp
    omega_rev             = 1.99099e-7_wp
    obliquity             = 0.409092592843564_wp
    
    ! refer to https://www.esrl.noaa.gov/gmd/grad/solcalc/azel.html
    ! Unix time coordinate of 2019-Dec-20, 12:00 UTC
    t_0                   = 1576843200.0_wp
    ! this is a winter solstice
    phi_0_earth_around_sun = 0._wp
    phi_0_earth_rotation  = 0._wp
    
    ! transformation of the time coordinate
    t_transformed =  t - t_0
    
    rot_angle =  omega*t_transformed - phi_0_earth_rotation
    
    ! the normal vector of the place we look at in earth fixed coordinates
    normal_vector_rel2_earth(1) =  cos(lat)*cos(lon)
    normal_vector_rel2_earth(2) =  cos(lat)*sin(lon)
    normal_vector_rel2_earth(3) =  sin(lat)
    
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
    trans_earth2sun(2,3) = 0._wp
    trans_earth2sun(3,3) = cos(obliquity)
    
    ! transforming the normal vector of the place to solar coordinates
    normal_vector_rel2_sun = matmul(trans_earth2sun, normal_vector_rel2_earth)
    
    sun_2_earth             (1) = cos(omega_rev*t_transformed + phi_0_earth_around_sun)
    sun_2_earth             (2) = sin(omega_rev*t_transformed + phi_0_earth_around_sun)
    sun_2_earth             (3) = 0._wp
    
    ! the result
    coszenith = DOT_PRODUCT(normal_vector_rel2_sun, -sun_2_earth)
    
    ! the night case
    if (coszenith < 0._wp) then
      coszenith = 0._wp
    endif
  
  end function coszenith
  
  subroutine set_vol_mix_ratios(mass_densities, sw_bool, no_of_day_points, no_of_scalars_h, &
  no_of_layers, no_of_scalars, no_of_condensed_constituents, day_indices, z_scalar)
    
    ! computes volume mixing ratios out of the model variables
    
    ! mass densities of the constituents
    real(8), intent(in)              :: mass_densities(:)
    ! short wave switch
    logical, intent(in)              :: sw_bool
    ! as usual
    integer, intent(in)              :: no_of_day_points
    ! as usual
    integer, intent(in)              :: no_of_scalars_h
    ! as usual
    integer, intent(in)              :: no_of_layers
    ! as usual
    integer, intent(in)              :: no_of_scalars
    ! as usual
    integer, intent(in)              :: no_of_condensed_constituents
    ! the indices of the points where it is day
    integer, intent(in)              :: day_indices(no_of_scalars/no_of_layers)
    ! z coordinates of scalar data points
    real(8), intent(in)              :: z_scalar(no_of_scalars)
    
    ! the volume mixing ratio of a gas
    real(8)                          :: vol_mix_ratio(no_of_scalars_h, no_of_layers)
    ! loop indices
    integer                          :: ji,jk,jl
    
    ! setting the volume mixing ratios of the gases
    do ji = 1,size(active_gases)
      ! the default
      vol_mix_ratio(:,:) = 0.0_wp
      select case (gases_lowercase(ji))
        ! reading the VMRs from atmostracers library
        case("n2")
          vol_mix_ratio(:,:) = molar_fraction_in_dry_air(2)
        case("o2")
          vol_mix_ratio(:,:) = molar_fraction_in_dry_air(3)
        case("ch4")
          vol_mix_ratio(:,:) = molar_fraction_in_dry_air(8)
        case("o3")
          if (sw_bool) then
            do jk=1,no_of_day_points
              do jl=1,no_of_layers
                vol_mix_ratio(jk,jl) = calc_o3_vmr(z_scalar(day_indices(jk)+(jl-1)*no_of_scalars_h))
              enddo
            enddo
          else
            do jk=1,no_of_scalars_h
              do jl=1,no_of_layers
                vol_mix_ratio(jk,jl) = calc_o3_vmr(z_scalar(jk+(jl-1)*no_of_scalars_h))
              enddo
            enddo
          endif
        case("co2")
          vol_mix_ratio(:,:) = molar_fraction_in_dry_air(5)
        case("co")
          vol_mix_ratio(:,:) = molar_fraction_in_dry_air(9)
        case("n2o")
          vol_mix_ratio(:,:) = molar_fraction_in_dry_air(11)
        case("h2o")
          ! no_of_condensed_constituents > 0 is equivalent to the presence of water in the model atmosphere
          ! in the short wave case, only the day points matter
          if (sw_bool .and. no_of_condensed_constituents > 0) then
            do jk=1,no_of_day_points
              do jl=1,no_of_layers
                vol_mix_ratio(jk,jl) = & 
                mass_densities((no_of_condensed_constituents+1)*no_of_scalars+day_indices(jk)+(jl-1)*no_of_scalars_h) &
                *specific_gas_constants(1)/ &
                (mass_densities(no_of_condensed_constituents*no_of_scalars+day_indices(jk)+(jl-1)*no_of_scalars_h) &
                *specific_gas_constants(0))
              enddo
            enddo
          ! in the long wave case, all points matter
          elseif (no_of_condensed_constituents > 0) then
            do jk=1,no_of_scalars_h
              do jl=1,no_of_layers
                vol_mix_ratio(jk,jl) = & 
                mass_densities((no_of_condensed_constituents+1)*no_of_scalars+jk+(jl-1)*no_of_scalars_h) &
                *specific_gas_constants(1)/ &
                (mass_densities(no_of_condensed_constituents*no_of_scalars+jk+(jl-1)*no_of_scalars_h) &
                *specific_gas_constants(0))
              enddo
            enddo
          endif
        end select
      ! finally setting the VMRs to the gas_concentrations objects
      if (sw_bool) then
        call handle_error(gas_concentrations_sw%set_vmr(gases_lowercase(ji), vol_mix_ratio(1:no_of_day_points,:)))
      else
        call handle_error(gas_concentrations_lw%set_vmr(gases_lowercase(ji), vol_mix_ratio(:,:)))
      endif
    enddo ! ji
  
  end subroutine set_vol_mix_ratios
  
  subroutine init_fluxes(fluxes, n_hor, n_vert, n_bands)
  
    ! initializing a flux object
    ! the fluxes to initialize
    type(ty_fluxes_broadband), intent(inout) :: fluxes
    ! the number of columns
    integer, intent(in)                   :: n_hor
    ! the number of levels
    integer, intent(in)                   :: n_vert
    ! the number of bads
    integer, intent(in)                   :: n_bands
 	
 	! broad band fluxes
    allocate(fluxes%flux_up (n_hor, n_vert))
    allocate(fluxes%flux_dn (n_hor, n_vert))
    allocate(fluxes%flux_net(n_hor, n_vert))
    
    call reset_fluxes(fluxes)
    
  end subroutine init_fluxes
  
  subroutine reset_fluxes(fluxes)

    ! resets all fluxes to zero

    type(ty_fluxes_broadband), intent(inout) :: fluxes

    ! reset broadband fluxes
    fluxes%flux_up(:,:) =  0._wp
    fluxes%flux_dn(:,:) =  0._wp
    fluxes%flux_net(:,:) =  0._wp
    if (associated(fluxes%flux_dn_dir)) fluxes%flux_dn_dir(:,:) =  0._wp

  end subroutine reset_fluxes
  
  subroutine free_fluxes(fluxes)
  
    ! freeing a flux object
    ! the fluxes to free
    type(ty_fluxes_broadband), intent(inout) :: fluxes
    
    if (associated(fluxes%flux_up)) deallocate(fluxes%flux_up)
    if (associated(fluxes%flux_dn)) deallocate(fluxes%flux_dn)
    if (associated(fluxes%flux_net)) deallocate(fluxes%flux_net)
    if (associated(fluxes%flux_dn_dir)) deallocate(fluxes%flux_dn_dir)
  
  end subroutine free_fluxes
  
  subroutine handle_error(error_message)
  
    character(len = *), intent(in) :: error_message
    
    ! write the error message if its real length is larger than zero
    if (len(trim(error_message)) > 0) then
      write(*,*) error_message
    endif
  
  end subroutine handle_error
  
end module radiation
















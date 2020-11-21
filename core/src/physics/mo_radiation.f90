! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/MHBalsmeier/game

module radiation
  
  use iso_c_binding,
  use mo_rte_kind,           only: wp
  use mo_rrtmgp_util_string, only: lower_case
  use mo_gas_optics_rrtmgp,  only: ty_gas_optics_rrtmgp
  use mo_load_coefficients,  only: load_and_init
  use mo_gas_concentrations, only: ty_gas_concs
  use mo_fluxes_byband,      only: ty_fluxes_byband
  use mo_rrtmgp_clr_all_sky, only: rte_sw
  use mo_rrtmgp_clr_all_sky, only: rte_lw
  use mo_optical_props,      only: ty_optical_props_2str, &
                                   ty_optical_props_1scl
  
  implicit none
  
  ! the number of bands in the short wave region
  integer, parameter                 :: no_of_sw_bands = 14
  ! the number of bands in the long wave region
  integer, parameter                 :: no_of_lw_bands = 16
  ! the number of g points is the short wave region
  integer                            :: no_of_sw_g_points
  ! the number of g points is the long wave region
  integer                            :: no_of_lw_g_points
  ! the gas concentrations (object holding all information on the composition
  ! of the gas phase)
  type(ty_gas_concs)                 :: gas_concentrations_sw
  type(ty_gas_concs)                 :: gas_concentrations
  
  type(ty_gas_optics_rrtmgp)         :: k_dist_sw, k_dist_lw

  character(len=3), dimension(7) :: active_gases = (/ &
   'N2 ', 'O2 ', 'CH4', 'O3 ', 'CO2', 'H2O', 'N2O' &
   /)
  
  character(len=*), parameter      :: rrtmgp_coefficients_file_sw = &
  ! insert the name of the short wave data file here
  '/home/max/code/rte-rrtmgp/rrtmgp/data/rrtmgp-data-sw-g224-2018-12-04.nc'
  character(len=*), parameter      :: rrtmgp_coefficients_file_lw = &
  ! insert the name of the long wave data file here
  '/home/max/code/rte-rrtmgp/rrtmgp/data/rrtmgp-data-lw-g256-2018-12-04.nc'
  ! the gases in lowercase
  character(len=32), dimension(size(active_gases)) :: gases_lowercase
  
  ! interface to C function
  interface
    real(8) function specific_gas_constants(gas_number) bind(c, name = "specific_gas_constants")
      integer :: gas_number
    end function specific_gas_constants
  end interface
    
  contains
  
  subroutine radiation_init() &
  bind(c, name = "radiation_init")
    
    ! This is called only once, in the beginning.
    
    ! local variables
    ! loop index
    integer                          :: ji
    
    no_of_sw_g_points = k_dist_sw%get_ngpt()
    no_of_lw_g_points = k_dist_lw%get_ngpt()
    
    ! formatting the gas names
    do ji = 1,size(active_gases)
      gases_lowercase(ji) = trim(lower_case(active_gases(ji)))
    end do
    ! here, the names of the gases are written to the gas_concentrations object
    call handle_error(gas_concentrations_sw%init(gases_lowercase))
    call handle_error(gas_concentrations%init(gases_lowercase))
    
    ! loading the short wave radiation properties
    call load_and_init(k_dist_sw, rrtmgp_coefficients_file_sw, gas_concentrations_sw)
    ! loading the long wave radiation properties
    call load_and_init(k_dist_lw, rrtmgp_coefficients_file_lw, gas_concentrations)
    
  end subroutine radiation_init
  
  subroutine calc_radiative_flux_convergence(latitude_scalar, longitude_scalar, &
  volume, area, &
  mass_densities, temperature_gas, radiation_tendency, &
  no_of_scalars, no_of_vectors, no_of_vectors_per_layer, &
  no_of_layers, no_of_constituents, no_of_condensed_constituents, &
  time_coord) &
  bind(c, name = "calc_radiative_flux_convergence")
    
    integer, intent(in)              ::                    no_of_scalars
    integer, intent(in)              ::                    no_of_vectors
    integer, intent(in)              ::                    no_of_vectors_per_layer
    integer, intent(in)              ::                    no_of_layers
    integer, intent(in)              ::                    no_of_constituents
    integer, intent(in)              ::                    no_of_condensed_constituents
    real(8)                          :: time_coord
    real(8), intent(in)              :: latitude_scalar    (no_of_scalars/no_of_layers)
    real(8), intent(in)              :: longitude_scalar   (no_of_scalars/no_of_layers)
    real(8), intent(in)              :: volume             (no_of_scalars)
    real(8), intent(in)              :: area               (no_of_scalars)
    real(8), intent(in)              :: mass_densities    &
    (no_of_constituents*no_of_scalars)
    real(8), intent(in)              :: temperature_gas   (no_of_scalars)
    real(8), intent(inout)           :: radiation_tendency(no_of_scalars)
    
    ! local variables
    ! solar zenith angle
    real(8)	                         :: mu_0(no_of_scalars/no_of_layers)
    ! number of points where it is day
    integer                          :: no_of_day_points
    ! number of points where it is night
    integer                          :: no_of_night_points
    ! loop indices
    integer                          :: ji, j_day, j_night, jk
    ! the indices of columns where it is day
    integer                          :: day_indices(no_of_scalars/no_of_layers)
    ! the indices of columns where it is night
    integer                          :: night_indices(no_of_scalars/no_of_layers)
    ! number of scalars per layer (number of columns)
    integer                          :: no_of_scalars_h
    ! the resulting clear sky fluxes
    type(ty_fluxes_byband)           :: fluxes_clearsky, fluxes_clearsky_day
    ! the resulting all sky fluxes
    type(ty_fluxes_byband)           :: fluxes_allsky, fluxes_allsky_day
    real(8)                          :: surface_emissivity(no_of_lw_bands, no_of_scalars/no_of_layers)
    ! surface albedo for direct radiation
    real(8)                          :: albedo_dir        (no_of_lw_bands, no_of_scalars/no_of_layers)
    ! surface albedo for diffusive radiation
    real(8)                          :: albedo_dif        (no_of_lw_bands, no_of_scalars/no_of_layers)
    ! temperature at cells
    real(8)                          :: temperature_rad           (no_of_scalars/no_of_layers, no_of_layers)
    ! pressure at cells
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
    ! the volume mixing ratio of a gas
    real(8)                          :: vol_mix_ratio             (no_of_scalars/no_of_layers, no_of_layers+1)
    ! cloud optical properties
    type(ty_optical_props_2str)      :: cloud_optics_sw
    type(ty_optical_props_1scl)      :: cloud_optics_lw
    
    ! calculation of the number of columns
    no_of_scalars_h = no_of_scalars/no_of_layers
    
    ! set the surface emissivity (a longwave property) to one
    surface_emissivity(:,:) = 1._wp
    
    ! set the surface albedos to 0.5
    albedo_dir        (:,:) = 0.5_wp
    albedo_dif        (:,:) = 0.5_wp
    
    ! reformatting the thermodynamical state for RTE+RRTMGP
    do ji=1,no_of_scalars_h
      do jk=1,no_of_layers
        temperature_rad(ji,jk) = temperature_gas((jk-1)*no_of_scalars_h+ji)
        ! the pressure is diagnozed here, using the equation of state for ideal gases
        pressure_rad(ji,jk)    = specific_gas_constants(0) &
        *mass_densities(no_of_condensed_constituents*no_of_scalars &
        + (jk-1)*no_of_scalars_h+ji)*temperature_rad(ji,jk)
      enddo
    enddo
    ! the properties at cell interfaces
    do ji=1,no_of_scalars_h
      do jk=1,no_of_layers+1
        if (jk==1) then
          temperature_interface_rad(ji,jk) = temperature_rad(ji,jk)
          pressure_interface_rad   (ji,jk) = pressure_rad   (ji,jk)
        elseif (jk==no_of_layers+1) then
          temperature_interface_rad(ji,jk) = temperature_rad(ji,jk-1)
          pressure_interface_rad   (ji,jk) = pressure_rad   (ji,jk-1)
        else
          temperature_interface_rad(ji,jk) = 0.5*(temperature_rad(ji,jk-1)+temperature_rad(ji,jk))
          pressure_interface_rad   (ji,jk) = 0.5*(pressure_rad   (ji,jk-1)+pressure_rad   (ji,jk))
        endif
      enddo
    enddo
    
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
    
    no_of_day_points   = j_day
    no_of_night_points = j_night
    
    ! filling up the arrays restricted to day points
    do j_day = 1,no_of_day_points
      temperature_rad_day(j_day,:)        = temperature_rad(day_indices(j_day),:)
      pressure_rad_day(j_day,:)           = pressure_rad(day_indices(j_day),:)
      pressure_interface_rad_day(j_day,:) = pressure_interface_rad(day_indices(j_day),:)
    end do
    
    ! setting the volume mixing ratios of the gases for the short wave calculation
    do ji=1,size(active_gases)
      ! the default
      vol_mix_ratio(:,:)=0.0_wp
      select case (gases_lowercase(ji))
        case('n2')
          vol_mix_ratio(:,:)=0.8_wp
        case('o2')
          vol_mix_ratio(:,:)=0.2_wp
      end select
      call handle_error(gas_concentrations_sw%set_vmr(gases_lowercase(ji), vol_mix_ratio &
      (1:no_of_day_points,1:no_of_layers+1)))
    enddo
   
    ! initializing the short wave fluxes
    call init_fluxes(fluxes_clearsky_day, no_of_day_points, no_of_layers+1, no_of_sw_bands)
    call init_fluxes(fluxes_allsky_day, no_of_day_points, no_of_layers+1, no_of_sw_bands)
    
    ! calculate shortwave radiative fluxes (only the day points are handed over
    ! for efficiency)
    write(*,*) "a"
    call handle_error(rte_sw(k_dist_sw, gas_concentrations_sw, pressure_rad_day, &
    temperature_rad_day, pressure_interface_rad_day, &
    mu_0, albedo_dir, albedo_dif, cloud_optics_sw, fluxes_allsky_day, fluxes_clearsky_day))
    write(*,*) "a"
    
    ! clearing the radiation tendency
    do ji=1,no_of_scalars
      radiation_tendency(ji)=0._wp
    enddo
    ! short wave result (in Wm^-3)
    ! clear sky
    call calc_power_density(.TRUE., no_of_scalars, no_of_vectors, &
    no_of_layers, no_of_scalars_h, no_of_vectors_per_layer, no_of_day_points, day_indices, &
    fluxes_clearsky_day, volume, area, radiation_tendency)
    ! all sky
    call calc_power_density(.TRUE., no_of_scalars, no_of_vectors, &
    no_of_layers, no_of_scalars_h, no_of_vectors_per_layer, no_of_day_points, day_indices, &
    fluxes_allsky_day, volume, area, radiation_tendency)
    
    ! freeing the short wave fluxes
    call free_fluxes(fluxes_clearsky_day)
    call free_fluxes(fluxes_allsky_day)
    
    ! setting the volume mixing ratios of the gases for the long wave calculation
    do ji=1,size(active_gases)
      ! the default
      vol_mix_ratio(:,:)=0.0_wp
      select case (gases_lowercase(ji))
        case('n2')
          vol_mix_ratio(:,:)=0.8_wp
        case('o2')
          vol_mix_ratio(:,:)=0.2_wp
      end select
      call handle_error(gas_concentrations%set_vmr(gases_lowercase(ji), vol_mix_ratio &
      (1:no_of_scalars_h,1:no_of_layers+1)))
    enddo
    
    ! initializing the long wave fluxes
    call init_fluxes(fluxes_clearsky, no_of_scalars_h, no_of_layers+1, no_of_lw_bands)
    call init_fluxes(fluxes_allsky, no_of_scalars_h, no_of_layers+1, no_of_lw_bands)
    ! calculate longwave radiative fluxes
    write(*,*) "a"
    call handle_error(rte_lw(k_dist_lw, gas_concentrations, pressure_rad, &
    temperature_rad, pressure_interface_rad, temperature_interface_rad(:,no_of_layers+1), &
    surface_emissivity(:,:), cloud_optics_lw, fluxes_allsky, fluxes_clearsky, &
    t_lev=temperature_interface_rad(:,:)))
    write(*,*) "a"
   
    ! add long wave result (in Wm^-3)
    ! clear sky
    call calc_power_density(.FALSE., no_of_scalars, no_of_vectors, &
    no_of_layers, no_of_scalars_h, no_of_vectors_per_layer, no_of_day_points, day_indices, &
    fluxes_clearsky, volume, area, radiation_tendency)
    ! all sky
    call calc_power_density(.FALSE., no_of_scalars, no_of_vectors, &
    no_of_layers, no_of_scalars_h, no_of_vectors_per_layer, no_of_day_points, day_indices, &
    fluxes_allsky, volume, area, radiation_tendency)
    
    ! freeing the long wave fluxes
    call free_fluxes(fluxes_clearsky)
    call free_fluxes(fluxes_allsky)
    
  end subroutine calc_radiative_flux_convergence
    
  subroutine calc_power_density(day_only, no_of_scalars, no_of_vectors, &
  no_of_layers, no_of_scalars_h, no_of_vectors_per_layer, no_of_day_points, day_indices, &
  fluxes, volume, area, radiation_tendency)
  
    ! this is essentially the negative vertical divergence operator
    
    ! true for short wave calculations (for efficiency)
    logical, intent(in)              :: day_only
    integer, intent(in)              :: no_of_scalars
    integer, intent(in)              :: no_of_vectors
    integer, intent(in)              :: no_of_layers
    integer, intent(in)              :: no_of_scalars_h
    integer, intent(in)              :: no_of_vectors_per_layer
    integer, intent(in)              :: no_of_day_points
    integer, intent(in)              :: day_indices(no_of_day_points)
    type(ty_fluxes_byband), intent(in):: fluxes
    real(8), intent(in)              :: volume(no_of_scalars)
    real(8), intent(in)              :: area(no_of_vectors)
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
      no_of_relevant_columns=size(day_indices)
    else
      no_of_relevant_columns=no_of_scalars_h
    endif
  
  
    do ji=1,no_of_layers
      do j_column=1,no_of_relevant_columns
        if (day_only) then
          jk=day_indices(j_column)
        else
          jk=j_column
        endif
        radiation_tendency((ji-1)*no_of_scalars_h+jk) = &
        radiation_tendency((ji-1)*no_of_scalars_h+jk) +&
        ! this is a sum of four fluxes
        ( &
        ! upward flux (going in)
        fluxes%flux_up  (j_column,ji+1)*area((ji    *no_of_vectors_per_layer)+jk) &
        ! upward flux (going out)
        - fluxes%flux_up(j_column,ji  )*area(((ji-1)*no_of_vectors_per_layer)+jk) &
        ! downward flux (going in)
        + fluxes%flux_dn(j_column,ji)  *area(((ji-1)*no_of_vectors_per_layer)+jk) &
        ! downward flux (going out)
        - fluxes%flux_dn(j_column,ji+1)*area((ji    *no_of_vectors_per_layer)+jk) &
        )/volume((ji-1)*no_of_scalars_h+jk)
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
  
  subroutine init_fluxes(fluxes, n_hor, n_vert, n_bands)
  
    ! initializing a flux object
    ! the fluxes to initialize
    type(ty_fluxes_byband), intent(inout) :: fluxes
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
  
 	! band-by-band fluxes
    allocate(fluxes%bnd_flux_up (n_hor, n_vert, n_bands))
    allocate(fluxes%bnd_flux_dn (n_hor, n_vert, n_bands))
    allocate(fluxes%bnd_flux_net(n_hor, n_vert, n_bands))
    
    call reset_fluxes(fluxes)
    
  end subroutine init_fluxes
  
  subroutine reset_fluxes(fluxes)

    type(ty_fluxes_byband), intent(inout) :: fluxes

    ! reset broadband fluxes
    fluxes%flux_up(:,:) = 0._wp
    fluxes%flux_dn(:,:) = 0._wp
    fluxes%flux_net(:,:) = 0._wp
    if (associated(fluxes%flux_dn_dir)) fluxes%flux_dn_dir(:,:) = 0._wp

    ! reset band-by-band fluxes
    fluxes%bnd_flux_up(:,:,:) = 0._wp
    fluxes%bnd_flux_dn(:,:,:) = 0._wp
    fluxes%bnd_flux_net(:,:,:) = 0._wp
    if (associated(fluxes%bnd_flux_dn_dir)) fluxes%bnd_flux_dn_dir(:,:,:) = 0._wp

  end subroutine reset_fluxes
  
  subroutine free_fluxes(fluxes)
  
    ! freeing a flux object
    ! the fluxes to free
    type(ty_fluxes_byband), intent(inout) :: fluxes
    
    if (associated(fluxes%flux_up)) deallocate(fluxes%flux_up)
    if (associated(fluxes%flux_dn)) deallocate(fluxes%flux_dn)
    if (associated(fluxes%flux_net)) deallocate(fluxes%flux_net)
    if (associated(fluxes%flux_dn_dir)) deallocate(fluxes%flux_dn_dir)
    if (associated(fluxes%bnd_flux_up)) deallocate(fluxes%bnd_flux_up)
    if (associated(fluxes%bnd_flux_dn)) deallocate(fluxes%bnd_flux_dn)
    if (associated(fluxes%bnd_flux_net)) deallocate(fluxes%bnd_flux_net)
    if (associated(fluxes%bnd_flux_dn_dir)) deallocate(fluxes%bnd_flux_dn_dir)
  
  end subroutine free_fluxes
  
  subroutine handle_error(error_message)
  
    character(len=*), intent(in) :: error_message
    
    if (len(trim(error_message)) > 0) then
      write(*,*) error_message
    endif
  
  end subroutine handle_error

end module radiation
















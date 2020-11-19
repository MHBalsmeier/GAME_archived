! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/MHBalsmeier/game

module radiation
  
  use mo_gas_optics_rrtmgp, only: ty_gas_optics_rrtmgp
  use mo_rte_kind, only: wp
  
  implicit none
  
  contains
  
  subroutine calc_radiative_flux_convergence() &
  bind(c, name = "calc_radiative_flux_convergence")
  
    
    
  end subroutine calc_radiative_flux_convergence

end module radiation

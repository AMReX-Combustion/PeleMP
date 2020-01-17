subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(c)

  use probdata_module
  use amrex_error_module
  implicit none

  integer :: init, namlen
  integer :: name(namlen)
  double precision :: problo(3), probhi(3)

  integer untin,i

  namelist /fortin/ u_ref, p_ref, T_ref, Y_O2, Y_N2
    
  ! Build "probin" filename -- the name of file containing fortin namelist.
  integer, parameter :: maxlen = 256
  character probin*(maxlen)

  if (namlen .gt. maxlen) then
     call amrex_error('probin file name too long')
  end if

  do i = 1, namlen
     probin(i:i) = char(name(i))
  end do

  ! Read namelists
  untin = 9
  open(untin,file=probin(1:namlen),form='formatted',status='old')
  read(untin,fortin)
  close(unit=untin)

end subroutine amrex_probinit


! ::: -----------------------------------------------------------
! ::: This routine is called at problem setup time and is used
! ::: to initialize data on each grid.
! :::
! ::: NOTE:  all arrays have one cell of ghost zones surrounding
! :::        the grid interior.  Values in these cells need not
! :::        be set here.
! :::
! ::: INPUTS/OUTPUTS:
! :::
! ::: level     => amr level of grid
! ::: time      => time at which to init data
! ::: lo,hi     => index limits of grid interior (cell centered)
! ::: nstate    => number of state components.  You should know
! :::		   this already!
! ::: state     <=  Scalar array
! ::: delta     => cell size
! ::: xlo,xhi   => physical locations of lower left and upper
! :::              right hand corner of grid.  (does not include
! :::		   ghost region).
! ::: -----------------------------------------------------------
  subroutine pc_initdata(level,time,lo,hi,nvar, &
       state,state_lo,state_hi, &
       delta,xlo,xhi) bind(C, name="pc_initdata")

  use probdata_module
  use network, only: nspecies, naux
  use chemistry_module, only : get_species_index
  use eos_type_module
  use meth_params_module, only : URHO, UMX, UMY, UMZ, &
       UEDEN, UEINT, UFS, UTEMP, small_temp
  use amrex_constants_module, only: ZERO, HALF, M_PI
  !use extern_probin_module, only: const_viscosity, const_bulk_viscosity, const_conductivity, const_diffusivity
  use eos_module
  
  implicit none

  integer :: level, nvar
  integer :: lo(3), hi(3)
  integer :: state_lo(3),state_hi(3)
  double precision :: xlo(3), xhi(3), time, delta(3)
  double precision :: state(state_lo(1):state_hi(1), &
                            state_lo(2):state_hi(2), &
                            state_lo(3):state_hi(3),nvar)

  integer :: i,j,k
  double precision :: x,y,z,rho,u,v,w,p,eint
  double precision :: p0, T0
  integer, parameter :: out_unit=20
  
  type(eos_t) :: eos_state
  
  integer :: iN2, iO2

  iN2 = get_species_index("N2")
  iO2 = get_species_index("O2")
  call build(eos_state)

  ! Initial pressure and temperature
  p0 = p_ref ! [erg cm^-3]
  T0 = T_ref

  ! Set the equation of state variables
  eos_state % p = p0
  eos_state % T = T0
  eos_state % massfrac    = 0.d0
  eos_state % massfrac(iO2) = Y_O2
  eos_state % massfrac(iN2) = Y_N2


  call eos_tp(eos_state)
  write(*,*) 'Initial T', eos_state % T
  write(*,*) 'Initial p ', eos_state % p
  write(*,*) 'Initial rho ', eos_state % rho
  write(*,*) 'Initial e ',  eos_state % e
  state(:,:,:,UTEMP) = T0

  do k = lo(3), hi(3)
     z = (k+HALF)*delta(3)

     do j = lo(2), hi(2)
        y = xlo(2) + delta(2)*(dble(j-lo(2)) + HALF)

        do i = lo(1), hi(1)
           x = xlo(1) + delta(1)*(dble(i-lo(1)) + HALF)

           u = u_ref
           v = 0.d0
           w = 0.d0

           eos_state % p   =  p0

           ! Call EOS by specifying the temperature and pressure
           call eos_tp(eos_state)
           rho  = eos_state % rho
           eint = eos_state % e

           ! Fill the states
           state(i,j,k,URHO)            = rho
           state(i,j,k,UFS:UFS+nspecies-1) = rho * eos_state % massfrac(1:nspecies)
           state(i,j,k,UMX)             = rho * u
           state(i,j,k,UMY)             = rho * v
           state(i,j,k,UMZ)             = rho * w
           state(i,j,k,UEINT)           = rho * eint
           state(i,j,k,UEDEN)           = rho * (eint + HALF * (u**2 + v**2 + w**2))

        enddo
     enddo
  enddo

end subroutine pc_initdata


subroutine pc_prob_close() &
     bind(C, name="pc_prob_close")

end subroutine pc_prob_close

module pc_prob_module

  implicit none

  private

  public :: amrex_probinit, pc_initdata, pc_prob_close

contains
  
  subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(C, name = "amrex_probinit")
    
    use probdata_module
    !use bl_error_module
    use turbinflow_module
    implicit none

    integer :: init, namlen
    integer :: name(namlen)
    double precision :: problo(3), probhi(3)
    character flct_file*(256)

    integer untin,i

    namelist /fortin/ p_ref, u_ref, T_ref, vel_jet, vel_coflow, jet_width, decay_factor, temp_jet, temp_coflow
    namelist /flctin/ flct_file

    ! Build "probin" filename -- the name of file containing fortin namelist.
    integer, parameter :: maxlen = 256
    character probin*(maxlen)

    if (namlen .gt. maxlen) then
   !   call bl_error('probin file name too long')
    end if

    do i = 1, namlen
       probin(i:i) = char(name(i))
    end do

    ! set namelist defaults here
    p_ref = 1013000.0d0
    u_ref = 10000.0d0
    T_ref = 300.0d0

    vel_jet = 2500.0d0 !m/s
    vel_coflow = -1000.0d0 !m/s
    jet_width = 0.02d0 !m
    decay_factor = 0.1d0 !non-dimensional

    temp_jet = 470d0 !K
    temp_coflow = 700d0 !K

    ! Read namelists
    untin = 9
    open(untin,file=probin(1:namlen),form='formatted',status='old')
    read(untin,fortin)
    read(untin,flctin)
    close(unit=untin)

    !Initialize the turbinflow module
    if(.not.turbinflow_initialized) then
      call init_turbinflow(flct_file, .true.)
    endif

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
  ! ::: nvar      => number of state components.
  ! ::: state     <= scalar array
  ! ::: dx        => cell size
  ! ::: xlo, xhi  => physical locations of lower left and upper
  ! :::              right hand corner of grid.  (does not include
  ! :::		   ghost region).
  ! ::: -----------------------------------------------------------

  subroutine pc_initdata(level,time,lo,hi,nvar, &
       state,state_lo,state_hi, &
       dx,xlo,xhi) bind(C, name="pc_initdata")
    use eos_type_module
    use probdata_module, only : p_ref, T_ref, vel_jet, vel_coflow, jet_width, decay_factor, temp_jet, temp_coflow
    use network, only: nspec, naux !, molec_wt
    use meth_params_module, only : URHO, UMX, UMY, UMZ, UEINT, UEDEN, UTEMP, UFS
    use prob_params_module, only : problo, probhi
    use eos_module
    use turbinflow_module
    use chemistry_module, only : get_species_index

    implicit none

    integer :: level, nvar
    integer :: lo(3), hi(3)
    integer :: state_lo(3), state_hi(3)
    double precision :: xlo(3), xhi(3), time, dx(3)
    double precision :: state(state_lo(1):state_hi(1), &
         state_lo(2):state_hi(2), &
         state_lo(3):state_hi(3),nvar)

    ! local variables
    integer :: i, j, k, n
    double precision :: x(lo(1):hi(1)), y(lo(2):hi(2)), vfluc(lo(1):hi(1),lo(2):hi(2),3)
    double precision :: z, jet_prof, jet_loc_up, jet_loc_down, delta, yy, wmean, vel_tot(3)
    integer :: iN2, iO2
    double precision vflucsx,vflucsy,vflucsz

    type (eos_t) :: eos_state


    iN2 = get_species_index("N2")
    iO2 = get_species_index("O2")
    call build(eos_state)

    !Assuming jet is always centered around x = 0.0
    jet_loc_up = 0.5*jet_width
    jet_loc_down = -0.5*jet_width
    delta = decay_factor*0.5*jet_width

    eos_state % massfrac = 0.d0
    eos_state % massfrac(iO2) = 0.15d0
    eos_state % massfrac(iN2) = 0.85d0
    eos_state % p = p_ref
    eos_state % T = temp_jet

    do i = lo(1), hi(1)
       x(i) = problo(1) + dx(1)*(dble(i) + HALF)
    enddo

    do j = lo(2), hi(2)
       y(j) = problo(2) + dx(2)*(dble(j) + HALF)
    enddo

    !Loop over z planes, and initialize velocity to that read from turb file
    vflucsx = 0d0
    vflucsy = 0d0
    vflucsz = 0d0
    do k = lo(3), hi(3)
       z = problo(3) + dx(3)*(dble(k) + HALF)
       call get_turbvelocity(lo(1),lo(2),hi(1),hi(2),x,y,z,vfluc)

       do j = lo(2), hi(2)
          
          do i = lo(1), hi(1)

             !Compute jet profile to impose on (mean+fluctuation)
             jet_prof = tanh((x(i)-jet_loc_down)/delta) - tanh((x(i)-jet_loc_up)/delta)

             wmean = vel_jet*0.5d0*jet_prof + vel_coflow*(1.0d0 - 0.5d0*jet_prof)

             vflucsx = vflucsx+jet_prof*vfluc(i,j,1)
             vflucsy = vflucsy+jet_prof*vfluc(i,j,2)
             vflucsz = vflucsz+jet_prof*vfluc(i,j,3)
             vel_tot(1) = jet_prof*vfluc(i,j,1)
             vel_tot(2) = jet_prof*vfluc(i,j,2)
             vel_tot(3) = wmean + jet_prof*vfluc(i,j,3)
             
             eos_state % p = p_ref
             eos_state % T = temp_jet*0.5d0*jet_prof + temp_coflow*(1.0d0 - 0.5d0*jet_prof)           
             call eos_tp(eos_state)
             state(i,j,k,URHO )  = eos_state % rho

             state(i,j,k,UMX  )  = vel_tot(1) * eos_state % rho
             state(i,j,k,UMY  ) = vel_tot(2) * eos_state % rho
             state(i,j,k,UMZ  ) = vel_tot(3) * eos_state % rho 
             state(i,j,k,UEINT) = eos_state % rho  *  eos_state % e
             state(i,j,k,UEDEN) = eos_state % rho * (eos_state % e+0.5d0 *(vel_tot(1)**2+vel_tot(2)**2+vel_tot(3)**2))
             state(i,j,k,UTEMP) = eos_state % T
             do n=1, nspec
                state(i,j,k,UFS+n-1) = eos_state % rho  *  eos_state % massfrac(n)
             end do

          end do
       end do
    end do
  
  call destroy(eos_state)

  end subroutine pc_initdata

  subroutine pc_prob_close() &
       bind(C, name="pc_prob_close")
  end subroutine pc_prob_close

end module pc_prob_module

! This code is part of RRTM for GCM Applications - Parallel (RRTMGP)
!
! Contacts: Robert Pincus and Eli Mlawer
! email:  rrtmgp@aer.com
!
! Copyright 2015-2018,  Atmospheric and Environmental Research and
! Regents of the University of Colorado.  All right reserved.
!
! Use and duplication is permitted under the terms of the
!    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
! -------------------------------------------------------------------------------------------------
!
! Example program to demonstrate the calculation of shortwave radiative fluxes in clear, aerosol-free skies.
!   The example files come from the Radiative Forcing MIP (https://www.earthsystemcog.org/projects/rfmip/)
!   The large problem (1800 profiles) is divided into blocks
!
! Program is invoked as rrtmgp_rfmip_sw [block_size input_file  coefficient_file upflux_file downflux_file]
!   All arguments are optional but need to be specified in order.
!
! -------------------------------------------------------------------------------------------------
!
! Error checking: Procedures in rte+rrtmgp return strings which are empty if no errors occured
!   Check the incoming string, print it out and stop execution if non-empty
!
subroutine stop_on_err(error_msg)
  use iso_fortran_env, only : error_unit
  character(len=*), intent(in) :: error_msg

  if(error_msg /= "") then
    write (error_unit,*) trim(error_msg)
    write (error_unit,*) "rrtmg_rfmip_lw stopping"
    stop
  end if
end subroutine stop_on_err
! -------------------------------------------------------------------------------------------------
!
! Main program
!
! -------------------------------------------------------------------------------------------------
program rrtmg_rfmip_lw
  ! --------------------------------------------------
  !
  ! Modules for working with rte and rrtmgp
  !
  ! Working precision for real variables
  !
  use mo_rte_kind,           only: wp
  use mo_rfmip_io,           only: read_size, read_and_block_pt, &
    read_and_block_gases_ty, unblock_and_write, &
    read_and_block_lw_bc, read_and_block_sw_bc, read_kdist_gas_names
  use mo_gas_concentrations, only: ty_gas_concs, get_subset_range
  use rrtmg_lw_rad,          only: rrtmg_lw
  use rrtmg_lw_init,         only: rrtmg_lw_ini

  implicit none
  ! --------------------------------------------------
  !
  ! Local variables
  !
  character(len=132) :: rfmip_file, &
    kdist_file = 'coefficients_lw.nc', &
    flxdn_file = 'rld_template.nc', flxup_file = 'rlu_template.nc'
  integer :: nargs, ncol, nlay, nexp, nblocks, block_size, &
    dumInt=0, ngpt=16
  logical :: top_at_1
  integer :: b, icol, igpt, nband=1
  character(len=6) :: block_size_char

  character(len=32 ), dimension(:), allocatable :: gases_to_use

  ! block_size, nlay, nblocks
  real(wp), dimension(:,:,:), allocatable :: &
    p_lay, p_lev, t_lay, t_lev, dum3D

  ! block_size, nlay (nlev for HR)
  real(wp), dimension(:,:), allocatable :: &
    swuflx , swdflx, swhr, swuflxc, swdflxc, swhrc, dum2D

  ! gas concentrations (which differ by experiment, and some are 
  ! global averages)
  real(wp), dimension(:,:), allocatable :: &
    h2o, co2, o3, n2o, co, ch4, o2, n2, &
    cfc11, cfc12, cfc22, ccl4, &
    reicmcl, relqmcl, hr, uflxc, dflxc, hrc, &
    duflx_dt, duflxc_dt

  real(wp), dimension(:,:,:), allocatable :: flux_up, flux_dn, &
    lay_src, lev_src_dec, lev_src_inc, &
    cldfmcl, taucmcl, ciwpmcl, clwpmcl, tauaer
  real(wp), dimension(:,:), allocatable :: & ! block_size, nblocks
    emiss_sfc, t_sfc, conc

  type(ty_gas_concs), dimension(:), allocatable  :: gas_conc_array
  character(len=128) :: error_msg

  real(wp), dimension(:,:  ), allocatable :: sfc_emis ! block_size, nblocks

  rfmip_file = &
    'multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-0-4_none.nc'

  nargs = command_argument_count()
  if(nargs >= 2) call get_command_argument(2, rfmip_file)
  if(nargs >= 3) call get_command_argument(3, flxup_file)
  if(nargs >= 4) call get_command_argument(4, flxdn_file)

  ! How big is the problem? 
  ! Does it fit into blocks of the size we've specified?
  call read_size(rfmip_file, ncol, nlay, nexp)
  if(nargs >= 1) then
    call get_command_argument(1, block_size_char)
    read(block_size_char, '(i6)') block_size
  else
    block_size = ncol
  end if

  if(mod(ncol*nexp, block_size) /= 0 ) &
    call stop_on_err("# of columns doesn't fit evenly into blocks.")
  nblocks = (ncol*nexp)/block_size
  print *, "Doing ",  nblocks, "blocks of size ", block_size

  ! read in RFMIP specs
  call read_and_block_pt(rfmip_file, block_size, p_lay, p_lev, &
    t_lay, t_lev)

  ! RRTMG needs surface T
  call read_and_block_lw_bc(rfmip_file, block_size, emiss_sfc, t_sfc)

  ! are we going surface to TOA or TOA to surface?
  top_at_1 = p_lay(1, 1, 1) < p_lay(1, nlay, 1)

  ! RRTMGP won't run with pressure less than its minimum. The top 
  ! level in the RFMIP file is set to 10^-3 Pa. Here we pretend the 
  ! layer is just a bit less deep. This introduces an error but shows 
  ! input sanitizing. Pernak: just going with 10**-3
  if(top_at_1) then
    p_lev(:,1,:) = 1e-3
  else
    p_lev(:,nlay+1,:) = 1e-3
  end if

  ! Names of gases known to the k-distribution.
  ! Which gases will be included in the calculation?
  ! By default we'll use all the gases the k-distribution can handle, 
  ! but we could provide variants i.e. using equivalent concentrations 
  ! per RFMIP
  call read_kdist_gas_names(kdist_file, gases_to_use)
  print *, "Radiation calculation uses gases "
  print *, "  ", &
    [(trim(gases_to_use(b)) // " ", b = 1, size(gases_to_use))]

  ! grab concentrations
  call read_and_block_gases_ty(rfmip_file, block_size, &
    gases_to_use, gas_conc_array)

  allocate(h2o(block_size, nlay), co2(block_size, nlay), &
    o3(block_size, nlay), n2o(block_size, nlay), &
    co(block_size, nlay), ch4(block_size, nlay), &
    o2(block_size, nlay), n2(block_size, nlay), &
    cfc11(block_size, nlay), cfc12(block_size, nlay), &
    cfc22(block_size, nlay), ccl4(block_size, nlay))

  allocate(dum2D(block_size, nlay+1), dum3D(block_size, nlay, nband))
  allocate(flux_up(block_size, nlay+1, nblocks), &
           flux_dn(block_size, nlay+1, nblocks))

  allocate(cldfmcl(ngpt, block_size, nlay), &
    taucmcl(ngpt, block_size, nlay), ciwpmcl(ngpt, block_size, nlay), &
    clwpmcl(ngpt, block_size, nlay), tauaer(block_size, nlay, ngpt))

  allocate(emiss_sfc(block_size, nblocks), &
    t_sfc(block_size, nblocks), &
    reicmcl(block_size, nlay), relqmcl(block_size, nlay), &
    hr(block_size, nlay), uflxc(block_size, nlay+1), &
    dflxc(block_size, nlay+1), hrc(block_size, nlay), &
    duflx_dt(block_size, nlay+1), duflxc_dt(block_size, nlay+1))

  dum2D(:,:) = 0._wp
  dum3D(:,:,:) = 0._wp
  cfc11(:,:) = 0._wp 
  cfc12(:,:) = 0._wp
  cfc22(:,:) = 0._wp
  ccl4(:,:) = 0._wp
  reicmcl(:,:) = 0._wp
  relqmcl(:,:) = 0._wp
  hr(:,:) = 0._wp
  uflxc(:,:) = 0._wp
  dflxc(:,:) = 0._wp
  hrc(:,:) = 0._wp
  duflx_dt(:,:) = 0._wp
  duflxc_dt(:,:) = 0._wp
  cldfmcl(:,:,:) = 0._wp
  taucmcl(:,:,:) = 0._wp
  ciwpmcl(:,:,:) = 0._wp
  clwpmcl(:,:,:) = 0._wp
  tauaer(:,:,:) = 0._wp

  ! Loop over blocks
  call rrtmg_sw_ini(1004.64_wp)
  do b = 1, nblocks
    if (b .ne. 1) then
      cycle
    endif

    error_msg = gas_conc_array(b)%get_vmr('h2o', h2o(:,:))
    error_msg = gas_conc_array(b)%get_vmr('co2', co2(:,:))
    error_msg = gas_conc_array(b)%get_vmr('o3', o3(:,:))
    error_msg = gas_conc_array(b)%get_vmr('n2o', n2o(:,:))
    error_msg = gas_conc_array(b)%get_vmr('co', co(:,:))
    error_msg = gas_conc_array(b)%get_vmr('ch4', ch4(:,:))
    error_msg = gas_conc_array(b)%get_vmr('o2', o2(:,:))
    error_msg = gas_conc_array(b)%get_vmr('n2', n2(:,:))

    call rrtmg_lw_ini(1004.64_wp)
    call rrtmg_lw(block_size, nlay, 0, 0, &
      p_lay(:,:,b)/100.0, p_lev(:,:,b)/100.0, &
      t_lay(:,:,b), t_lev(:,:,b), t_sfc(:,b), & 
      h2o, o3, co2, ch4, n2o, o2, &
      cfc11, cfc12, cfc22, ccl4, emiss_sfc, &
      0, 0, 0, cldfmcl, &
      taucmcl ,ciwpmcl ,clwpmcl ,reicmcl ,relqmcl, tauaer, &
      flux_up(:,:,b), flux_dn(:,:,b), hr, uflxc, dflxc, hrc, &
      duflx_dt, duflxc_dt)

!  do icol = 1, block_size
!    if(.not. usecol(icol)) then
!      flux_up(icol,:,b)  = 0._wp
!      flux_dn(icol,:,b)  = 0._wp
!    end if
!  end do

  end do

  !call unblock_and_write(trim(flxup_file), 'rsu', flux_up)
  !call unblock_and_write(trim(flxdn_file), 'rsd', flux_dn)
end program rrtmg_rfmip_lw

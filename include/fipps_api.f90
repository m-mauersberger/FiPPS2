!
!  This program developed in FORTRAN is a Finite Element solver for linear-
!  static analyses as well as linearized stability analyses. It is inherently
!  coupled to the open-source panel method APAME for providing fluid-structure-
!  interaction capabilites.
!    
!  Copyright (C) 2024 TUD Dresden University of Technology
! 
!  This file is part of FiPPS².
!
!  FiPPS² is free software: you can redistribute and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!
!  FiPPS² is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
! =================================================================================================
!
!> @brief
!> Programmierschnittstelle (API) zur Modelldefinition ohne *.fipps-Dateien
!>
!> @details
!> Das Modul fipps_api ersetzt die Datei-basierte Modelldefinition (control.fipps
!> und die zugehoerigen *.fipps-Dateien, die von input_tf eingelesen werden)
!> durch Subroutinen mit Input-Argumenten.
!>
!> Verwendung (nach call init_default(fesim)):
!>
!>   call init_default(fesim)
!>   call fipps_api_control (fesim, sol=1, outputVTK=.true., err)
!>   call fipps_api_nodes   (fesim, cids, coords, err)
!>   call fipps_api_mat1s   (fesim, mids, yms, sms, nus, rhos, aths, trefs, ges, fids, err)
!>   call fipps_api_pbeams  (fesim, pids, mids, aas, i11s, i22s, i12s, its, t1s, t2s, angles, 'rad', nsms, err)
!>   call fipps_api_beam2s  (fesim, pids, nids, xis, n0s, err)
!>   call fipps_api_spc1s   (fesim, sids, dofs, n1s, thrus, nns, err)
!>   call fipps_api_forces  (fesim, lids, nids, facs, nis, err)
!>   call fipps_api_loads   (fesim, lcids, sfacis, lides, err)
!>   call fipps_api_subcases(fesim, scids, spcaddids, loadids, mpcaddids, &
!>                           & skipBuckling, upgeom, upstress, outputs, readApameInputs, err)
!>   call fipps_api_check   (fesim, err)
!>   call init_values(fesim)
!>
!> Konventionen:
!> - Jede Section muss (wie die zugehoerige *.fipps-Datei) genau einmal
!>   bereitgestellt werden (err = 3, falls schon vorhanden).
!> - Jeder Aufruf setzt das zugehoerige is_<section>-Flag von control.fipps.
!> - Die Reihenfolge der Elementsections muss beam2 -> quad8 -> lsolid20
!>   sein (wie in input_tf), damit die globalen Element-IDs (eid) identisch
!>   mit der Datei-Variante sind.
!> - fipps_api_subcases muss vor fipps_api_lam20s aufgerufen werden.
!>
!> Fehlercodes (err):
!>   0 - Erfolg
!>   1 - ungueltige Argumente (z.B. Anzahl <= 0)
!>   2 - ungueltiger Winkeltyp ('deg'/'rad' erwartet)
!>   3 - Section wurde bereits bereitgestellt
!>   4 - Pruefbedingung nicht erfuellt (z.B. lam20 ohne vorherige subcases)
!
module fipps_api

use fesimulation_typen

implicit none

private

public :: fipps_api_control
public :: fipps_api_nodes
public :: fipps_api_loads
public :: fipps_api_temperatures
public :: fipps_api_beam2temps
public :: fipps_api_quad8temps
public :: fipps_api_lsolid20temps
public :: fipps_api_p2loads
public :: fipps_api_spcadds
public :: fipps_api_mpcadds
public :: fipps_api_beam2s
public :: fipps_api_quad8s
public :: fipps_api_lsolid20s
public :: fipps_api_forces
public :: fipps_api_moments
public :: fipps_api_p8loads
public :: fipps_api_p20loads
public :: fipps_api_aeroload2ds
public :: fipps_api_aeroload3ds
public :: fipps_api_mat1s
public :: fipps_api_mat8s
public :: fipps_api_mat20s
public :: fipps_api_failure_tresca
public :: fipps_api_failure_mises
public :: fipps_api_failure_maxpstress
public :: fipps_api_failure_puck
public :: fipps_api_failure_hill
public :: fipps_api_failure_norris
public :: fipps_api_failure_fibre
public :: fipps_api_failure_maxstrain
public :: fipps_api_failure_cuntze
public :: fipps_api_failure_maxstrain3d
public :: fipps_api_failure_tsaiwu3d
public :: fipps_api_pshells
public :: fipps_api_pcomps
public :: fipps_api_pbeams
public :: fipps_api_plsolids
public :: fipps_api_spc1s
public :: fipps_api_mpcs
public :: fipps_api_coords
public :: fipps_api_subcases
public :: fipps_api_lam8s
public :: fipps_api_lam20s
public :: fipps_api_couplings
public :: fipps_api_contact_node_beam2
public :: fipps_api_contact_node_quad8
public :: fipps_api_contact_node_lsolid20
public :: fipps_api_aero_coupling_2d
public :: fipps_api_aero_coupling_3d
public :: fipps_api_check

contains

  !
  ! =================================================================================================
  !
  !> Ersetzt control.fipps: Rechen- und Ausgabesteuerung (var=value-Eintraege
  !> und Ausgabeflagge-Woerter). Alle Argumente sind optional, nicht
  !> uebergebene Werte bleiben unveraendert.
  !
  subroutine fipps_api_control(fesim, sol, numEigVal, blocksize, shellResPos, &
                                & xmin, xmax, ymin, ymax, zmin, zmax, aeroDisplEps, &
                                & outputVTK, outputUser, outputShort, outputKoopt, &
                                & outputAdviLa, outputOptitube, outputGlawi, &
                                & outputHyMoWi, outputElemCoord, outputBoundCond, &
                                & outputApamePressures, skipFailed, calculateTSE, &
                                & calculateElementalTSE, calculateReactForce, &
                                & globalReactForce, err)

    type(fe_simulation), intent(inout)          :: fesim
    integer,            intent(in), optional    :: sol, numEigVal, blocksize, shellResPos
    double precision,   intent(in), optional    :: xmin, xmax, ymin, ymax, zmin, zmax, aeroDisplEps
    logical,            intent(in), optional    :: outputVTK, outputUser, outputShort, outputKoopt
    logical,            intent(in), optional    :: outputAdviLa, outputOptitube, outputGlawi
    logical,            intent(in), optional    :: outputHyMoWi, outputElemCoord, outputBoundCond
    logical,            intent(in), optional    :: outputApamePressures
    logical,            intent(in), optional    :: skipFailed, calculateTSE
    logical,            intent(in), optional    :: calculateElementalTSE, calculateReactForce
    logical,            intent(in), optional    :: globalReactForce
    integer,            intent(out)             :: err

    err = 0

    if (present(sol))               fesim%sol                 = sol
    if (present(numEigVal))         fesim%numEigVal           = numEigVal
    if (present(blocksize))         fesim%blocksize           = blocksize
    if (present(shellResPos))       fesim%shellResPos         = shellResPos
    if (present(xmin))              fesim%ausgabe%xmin        = xmin
    if (present(xmax))              fesim%ausgabe%xmax        = xmax
    if (present(ymin))              fesim%ausgabe%ymin        = ymin
    if (present(ymax))              fesim%ausgabe%ymax        = ymax
    if (present(zmin))              fesim%ausgabe%zmin        = zmin
    if (present(zmax))              fesim%ausgabe%zmax        = zmax
    if (present(aeroDisplEps))      fesim%aeroDisplEps        = aeroDisplEps

    if (present(outputVTK))         fesim%ausgabe%outputVTK        = outputVTK
    if (present(outputUser))        fesim%ausgabe%outputUser       = outputUser
    if (present(outputShort))       fesim%ausgabe%outputShort      = outputShort
    if (present(outputKoopt))       fesim%ausgabe%outputKoopt      = outputKoopt
    if (present(outputAdviLa))      fesim%ausgabe%outputAdviLa     = outputAdviLa
    if (present(outputOptitube))    fesim%ausgabe%outputOptitube   = outputOptitube
    if (present(outputGlawi))       fesim%ausgabe%outputGlawi      = outputGlawi
    if (present(outputHyMoWi))      fesim%ausgabe%outputHyMoWi     = outputHyMoWi
    if (present(outputElemCoord))   fesim%ausgabe%outputElemCoord  = outputElemCoord
    if (present(outputBoundCond))   fesim%ausgabe%outputBoundCond  = outputBoundCond
    if (present(outputApamePressures)) fesim%ausgabe%outputApamePressures = outputApamePressures

    if (present(skipFailed))        fesim%skipFailed           = skipFailed
    if (present(calculateTSE))      fesim%calculateTSE         = calculateTSE
    if (present(calculateElementalTSE)) fesim%calculateElementalTSE = calculateElementalTSE
    if (present(calculateReactForce))   fesim%calculateReactForce   = calculateReactForce
    if (present(globalReactForce))      fesim%globalReactForce      = globalReactForce

  end subroutine fipps_api_control

  !
  ! =================================================================================================
  !
  !> Ersetzt nodes.fipps: Zeilenformat "cid x y z" (nid wird automatisch
  !> 1..n vergeben, wie im Datei-Einlesen).
  !
  subroutine fipps_api_nodes(fesim, cids, coords, err)

    type(fe_simulation), intent(inout)                     :: fesim
    integer,            intent(in), dimension(:)           :: cids
    double precision,   intent(in), dimension(:,:)         :: coords ! (3, n)
    integer,            intent(out)                        :: err

    integer :: n, ii

    err = 0
    n = size(cids,1)

    if (n .le. 0) then
      err = 1
      return
    end if
    if (allocated(fesim%knoten%nodes)) then
      err = 3
      return
    end if

    allocate(fesim%knoten%nodes(1:n))
    do ii = 1, n
      fesim%knoten%nodes(ii)%nid    = ii
      fesim%knoten%nodes(ii)%cid    = cids(ii)
      fesim%knoten%nodes(ii)%coords = coords(1:3, ii)
    end do
    fesim%num_nodes = n
    fesim%is_node   = .true.

  end subroutine fipps_api_nodes

  !
  ! =================================================================================================
  !
  !> Ersetzt loads.fipps: Zeilenformat "lcid sfaci lidi" (sfac=1 wie im
  !> Datei-Einlesen).
  !
  subroutine fipps_api_loads(fesim, lcids, sfacis, lides, err)

    type(fe_simulation), intent(inout)                     :: fesim
    integer,            intent(in), dimension(:)           :: lcids
    double precision,   intent(in), dimension(:)           :: sfacis
    integer,            intent(in), dimension(:)           :: lides
    integer,            intent(out)                        :: err

    integer :: n, ii

    err = 0
    n = size(lcids,1)

    if (n .le. 0) then
      err = 1
      return
    end if
    if (allocated(fesim%lasten%loads)) then
      err = 3
      return
    end if

    allocate(fesim%lasten%loads(1:n))
    do ii = 1, n
      fesim%lasten%loads(ii)%sfac  = 1.d0
      fesim%lasten%loads(ii)%lcid  = lcids(ii)
      fesim%lasten%loads(ii)%sfaci = sfacis(ii)
      fesim%lasten%loads(ii)%lidi  = lides(ii)
    end do
    fesim%is_load = .true.

  end subroutine fipps_api_loads

  !
  ! =================================================================================================
  !
  !> Ersetzt temperatures.fipps: Zeilenformat "lid nid temp"
  !
  subroutine fipps_api_temperatures(fesim, lids, nids, temps, err)

    type(fe_simulation), intent(inout)             :: fesim
    integer,            intent(in), dimension(:)   :: lids
    integer,            intent(in), dimension(:)   :: nids
    double precision,   intent(in), dimension(:)   :: temps
    integer,            intent(out)                :: err

    integer :: n, ii

    err = 0
    n = size(lids,1)

    if (n .le. 0) then
      err = 1
      return
    end if
    if (allocated(fesim%lasten%nodetemps)) then
      err = 3
      return
    end if

    allocate(fesim%lasten%nodetemps(1:n))
    do ii = 1, n
      fesim%lasten%nodetemps(ii)%lid  = lids(ii)
      fesim%lasten%nodetemps(ii)%nid  = nids(ii)
      fesim%lasten%nodetemps(ii)%temp = temps(ii)
    end do
    fesim%is_temperature = .true.

  end subroutine fipps_api_temperatures

  !
  ! =================================================================================================
  !
  !> Ersetzt beam2temps.fipps: Zeilenformat "lid eid temp"
  !
  subroutine fipps_api_beam2temps(fesim, lids, eids, temps, err)

    type(fe_simulation), intent(inout)             :: fesim
    integer,            intent(in), dimension(:)   :: lids
    integer,            intent(in), dimension(:)   :: eids
    double precision,   intent(in), dimension(:)   :: temps
    integer,            intent(out)                :: err

    integer :: n, ii

    err = 0
    n = size(lids,1)

    if (n .le. 0) then
      err = 1
      return
    end if
    if (allocated(fesim%lasten%beam2temps)) then
      err = 3
      return
    end if

    allocate(fesim%lasten%beam2temps(1:n))
    do ii = 1, n
      fesim%lasten%beam2temps(ii)%lid  = lids(ii)
      fesim%lasten%beam2temps(ii)%eid  = eids(ii)
      fesim%lasten%beam2temps(ii)%temp = temps(ii)
    end do
    fesim%is_beam2temp = .true.

  end subroutine fipps_api_beam2temps

  !
  ! =================================================================================================
  !
  !> Ersetzt quad8temps.fipps: Zeilenformat "lid eid temp"
  !
  subroutine fipps_api_quad8temps(fesim, lids, eids, temps, err)

    type(fe_simulation), intent(inout)             :: fesim
    integer,            intent(in), dimension(:)   :: lids
    integer,            intent(in), dimension(:)   :: eids
    double precision,   intent(in), dimension(:)   :: temps
    integer,            intent(out)                :: err

    integer :: n, ii

    err = 0
    n = size(lids,1)

    if (n .le. 0) then
      err = 1
      return
    end if
    if (allocated(fesim%lasten%quad8temps)) then
      err = 3
      return
    end if

    allocate(fesim%lasten%quad8temps(1:n))
    do ii = 1, n
      fesim%lasten%quad8temps(ii)%lid  = lids(ii)
      fesim%lasten%quad8temps(ii)%eid  = eids(ii)
      fesim%lasten%quad8temps(ii)%temp = temps(ii)
    end do
    fesim%is_quad8temp = .true.

  end subroutine fipps_api_quad8temps

  !
  ! =================================================================================================
  !
  !> Ersetzt lsolid20temps.fipps: Zeilenformat "lid eid temp"
  !
  subroutine fipps_api_lsolid20temps(fesim, lids, eids, temps, err)

    type(fe_simulation), intent(inout)             :: fesim
    integer,            intent(in), dimension(:)   :: lids
    integer,            intent(in), dimension(:)   :: eids
    double precision,   intent(in), dimension(:)   :: temps
    integer,            intent(out)                :: err

    integer :: n, ii

    err = 0
    n = size(lids,1)

    if (n .le. 0) then
      err = 1
      return
    end if
    if (allocated(fesim%lasten%lsolid20temps)) then
      err = 3
      return
    end if

    allocate(fesim%lasten%lsolid20temps(1:n))
    do ii = 1, n
      fesim%lasten%lsolid20temps(ii)%lid  = lids(ii)
      fesim%lasten%lsolid20temps(ii)%eid  = eids(ii)
      fesim%lasten%lsolid20temps(ii)%temp = temps(ii)
    end do
    fesim%is_lsolid20temp = .true.

  end subroutine fipps_api_lsolid20temps

  !
  ! =================================================================================================
  !
  !> Ersetzt p2loads.fipps: Zeilenformat
  !> "lid eid1 dir pi1 pi2 thru eid2"
  !
  subroutine fipps_api_p2loads(fesim, lids, eid1s, dirs, pis, thrus, eid2s, err)

    type(fe_simulation), intent(inout)             :: fesim
    integer,            intent(in), dimension(:)   :: lids
    integer,            intent(in), dimension(:)   :: eid1s
    integer,            intent(in), dimension(:)   :: dirs
    double precision,   intent(in), dimension(:,:) :: pis    ! (2, n)
    logical,            intent(in), dimension(:)   :: thrus
    integer,            intent(in), dimension(:)   :: eid2s
    integer,            intent(out)                :: err

    integer :: n, ii

    err = 0
    n = size(lids,1)

    if (n .le. 0) then
      err = 1
      return
    end if
    if (allocated(fesim%lasten%p2loads)) then
      err = 3
      return
    end if

    allocate(fesim%lasten%p2loads(1:n))
    do ii = 1, n
      fesim%lasten%p2loads(ii)%lid  = lids(ii)
      fesim%lasten%p2loads(ii)%eid1 = eid1s(ii)
      fesim%lasten%p2loads(ii)%dir  = dirs(ii)
      fesim%lasten%p2loads(ii)%pi   = pis(1:2, ii)
      fesim%lasten%p2loads(ii)%thru = thrus(ii)
      fesim%lasten%p2loads(ii)%eid2 = eid2s(ii)
    end do
    fesim%is_p2load = .true.

  end subroutine fipps_api_p2loads

  !
  ! =================================================================================================
  !
  !> Ersetzt spcadd.fipps: Zeilenformat "scid sid"
  !
  subroutine fipps_api_spcadds(fesim, scids, sids, err)

    type(fe_simulation), intent(inout)             :: fesim
    integer,            intent(in), dimension(:)   :: scids
    integer,            intent(in), dimension(:)   :: sids
    integer,            intent(out)                :: err

    integer :: n, ii

    err = 0
    n = size(scids,1)

    if (n .le. 0) then
      err = 1
      return
    end if
    if (allocated(fesim%randbedingungen%spcadds)) then
      err = 3
      return
    end if

    allocate(fesim%randbedingungen%spcadds(1:n))
    do ii = 1, n
      fesim%randbedingungen%spcadds(ii)%scid = scids(ii)
      fesim%randbedingungen%spcadds(ii)%sid  = sids(ii)
    end do
    fesim%is_spcadd = .true.

  end subroutine fipps_api_spcadds

  !
  ! =================================================================================================
  !
  !> Ersetzt mpcadd.fipps: Zeilenformat "scid sid"
  !
  subroutine fipps_api_mpcadds(fesim, scids, sids, err)

    type(fe_simulation), intent(inout)             :: fesim
    integer,            intent(in), dimension(:)   :: scids
    integer,            intent(in), dimension(:)   :: sids
    integer,            intent(out)                :: err

    integer :: n, ii

    err = 0
    n = size(scids,1)

    if (n .le. 0) then
      err = 1
      return
    end if
    if (allocated(fesim%randbedingungen%mpcadds)) then
      err = 3
      return
    end if

    allocate(fesim%randbedingungen%mpcadds(1:n))
    do ii = 1, n
      fesim%randbedingungen%mpcadds(ii)%scid = scids(ii)
      fesim%randbedingungen%mpcadds(ii)%sid  = sids(ii)
    end do
    fesim%is_mpcadd = .true.

  end subroutine fipps_api_mpcadds

  !
  ! =================================================================================================
  !
  ! Elementsektionen
  !
  ! Die globale Element-ID (eid) wird wie in input_tf fortlauf ueber alle
  ! Elementtypen vergeben. Dafuer muessen die Sektionen in der Reihenfolge
  ! beam2 -> quad8 -> lsolid20 aufgerufen werden.
  !
  subroutine get_element_base(fesim, base)

    type(fe_simulation), intent(in) :: fesim
    integer,            intent(out) :: base

    base = 0
    if (allocated(fesim%elemente%beam2s))    base = base + size(fesim%elemente%beam2s,1)
    if (allocated(fesim%elemente%quad8s))    base = base + size(fesim%elemente%quad8s,1)
    if (allocated(fesim%elemente%lsolid20s)) base = base + size(fesim%elemente%lsolid20s,1)

  end subroutine get_element_base

  !
  ! =================================================================================================
  !
  !> Ersetzt beam2.fipps: Zeilenformat
  !> "pid nid1 nid2 xi1 xi2 xi3 n0"
  !
  subroutine fipps_api_beam2s(fesim, pids, nids, xis, n0s, err)

    type(fe_simulation), intent(inout)             :: fesim
    integer,            intent(in), dimension(:)   :: pids
    integer,            intent(in), dimension(:,:) :: nids   ! (2, n)
    double precision,   intent(in), dimension(:,:) :: xis    ! (3, n)
    integer,            intent(in), dimension(:)   :: n0s
    integer,            intent(out)                :: err

    integer :: n, ii, base

    err = 0
    n = size(pids,1)

    if (n .le. 0) then
      err = 1
      return
    end if
    if (allocated(fesim%elemente%beam2s)) then
      err = 3
      return
    end if

    call get_element_base(fesim, base)
    allocate(fesim%elemente%beam2s(1:n))
    do ii = 1, n
      fesim%elemente%beam2s(ii)%eid  = base + ii
      fesim%elemente%beam2s(ii)%pid  = pids(ii)
      fesim%elemente%beam2s(ii)%nids = nids(1:2, ii)
      fesim%elemente%beam2s(ii)%xi   = xis(1:3, ii)
      fesim%elemente%beam2s(ii)%n0   = n0s(ii)
    end do
    fesim%is_beam2 = .true.

  end subroutine fipps_api_beam2s

  !
  ! =================================================================================================
  !
  !> Ersetzt quad8.fipps: Zeilenformat
  !> "pid nid1..nid8 theta" (offset=0 wie im Datei-Einlesen)
  !
  subroutine fipps_api_quad8s(fesim, pids, nids, thetas, err)

    type(fe_simulation), intent(inout)             :: fesim
    integer,            intent(in), dimension(:)   :: pids
    integer,            intent(in), dimension(:,:) :: nids   ! (8, n)
    double precision,   intent(in), dimension(:)   :: thetas
    integer,            intent(out)                :: err

    integer :: n, ii, base

    err = 0
    n = size(pids,1)

    if (n .le. 0) then
      err = 1
      return
    end if
    if (allocated(fesim%elemente%quad8s)) then
      err = 3
      return
    end if

    call get_element_base(fesim, base)
    allocate(fesim%elemente%quad8s(1:n))
    do ii = 1, n
      fesim%elemente%quad8s(ii)%eid    = base + ii
      fesim%elemente%quad8s(ii)%pid    = pids(ii)
      fesim%elemente%quad8s(ii)%nids   = nids(1:8, ii)
      fesim%elemente%quad8s(ii)%theta  = thetas(ii)
      fesim%elemente%quad8s(ii)%offset = 0.d0
    end do
    fesim%is_quad8 = .true.

  end subroutine fipps_api_quad8s

  !
  ! =================================================================================================
  !
  !> Ersetzt lsolid20.fipps: Zeilenformat "pid nid1..nid20"
  !
  subroutine fipps_api_lsolid20s(fesim, pids, nids, err)

    type(fe_simulation), intent(inout)             :: fesim
    integer,            intent(in), dimension(:)   :: pids
    integer,            intent(in), dimension(:,:) :: nids   ! (20, n)
    integer,            intent(out)                :: err

    integer :: n, ii, base

    err = 0
    n = size(pids,1)

    if (n .le. 0) then
      err = 1
      return
    end if
    if (allocated(fesim%elemente%lsolid20s)) then
      err = 3
      return
    end if

    call get_element_base(fesim, base)
    allocate(fesim%elemente%lsolid20s(1:n))
    do ii = 1, n
      fesim%elemente%lsolid20s(ii)%eid  = base + ii
      fesim%elemente%lsolid20s(ii)%pid  = pids(ii)
      fesim%elemente%lsolid20s(ii)%nids = nids(1:20, ii)
    end do
    fesim%is_lsolid20 = .true.

  end subroutine fipps_api_lsolid20s

  !
  ! =================================================================================================
  !
  !> Ersetzt forces.fipps: Zeilenformat "lid nid fac ni1 ni2 ni3"
  !
  subroutine fipps_api_forces(fesim, lids, nids, facs, nis, err)

    type(fe_simulation), intent(inout)             :: fesim
    integer,            intent(in), dimension(:)   :: lids
    integer,            intent(in), dimension(:)   :: nids
    double precision,   intent(in), dimension(:)   :: facs
    double precision,   intent(in), dimension(:,:) :: nis    ! (3, n)
    integer,            intent(out)                :: err

    integer :: n, ii

    err = 0
    n = size(lids,1)

    if (n .le. 0) then
      err = 1
      return
    end if
    if (allocated(fesim%lasten%forces)) then
      err = 3
      return
    end if

    allocate(fesim%lasten%forces(1:n))
    do ii = 1, n
      fesim%lasten%forces(ii)%lid  = lids(ii)
      fesim%lasten%forces(ii)%nid  = nids(ii)
      fesim%lasten%forces(ii)%fac  = facs(ii)
      fesim%lasten%forces(ii)%ni   = nis(1:3, ii)
    end do
    fesim%is_force = .true.

  end subroutine fipps_api_forces

  !
  ! =================================================================================================
  !
  !> Ersetzt moments.fipps: Zeilenformat "lid nid fac ni1 ni2 ni3"
  !
  subroutine fipps_api_moments(fesim, lids, nids, facs, nis, err)

    type(fe_simulation), intent(inout)             :: fesim
    integer,            intent(in), dimension(:)   :: lids
    integer,            intent(in), dimension(:)   :: nids
    double precision,   intent(in), dimension(:)   :: facs
    double precision,   intent(in), dimension(:,:) :: nis    ! (3, n)
    integer,            intent(out)                :: err

    integer :: n, ii

    err = 0
    n = size(lids,1)

    if (n .le. 0) then
      err = 1
      return
    end if
    if (allocated(fesim%lasten%moments)) then
      err = 3
      return
    end if

    allocate(fesim%lasten%moments(1:n))
    do ii = 1, n
      fesim%lasten%moments(ii)%lid  = lids(ii)
      fesim%lasten%moments(ii)%nid  = nids(ii)
      fesim%lasten%moments(ii)%fac  = facs(ii)
      fesim%lasten%moments(ii)%ni   = nis(1:3, ii)
    end do
    fesim%is_moment = .true.

  end subroutine fipps_api_moments

  !
  ! =================================================================================================
  !
  !> Ersetzt p8loads.fipps: Zeilenformat
  !> "lid eid1 cid pi1..pi4 thru eid2"
  !
  subroutine fipps_api_p8loads(fesim, lids, eid1s, cids, pis, thrus, eid2s, err)

    type(fe_simulation), intent(inout)             :: fesim
    integer,            intent(in), dimension(:)   :: lids
    integer,            intent(in), dimension(:)   :: eid1s
    integer,            intent(in), dimension(:)   :: cids
    double precision,   intent(in), dimension(:,:) :: pis    ! (4, n)
    logical,            intent(in), dimension(:)   :: thrus
    integer,            intent(in), dimension(:)   :: eid2s
    integer,            intent(out)                :: err

    integer :: n, ii

    err = 0
    n = size(lids,1)

    if (n .le. 0) then
      err = 1
      return
    end if
    if (allocated(fesim%lasten%p8loads)) then
      err = 3
      return
    end if

    allocate(fesim%lasten%p8loads(1:n))
    do ii = 1, n
      fesim%lasten%p8loads(ii)%lid  = lids(ii)
      fesim%lasten%p8loads(ii)%eid1 = eid1s(ii)
      fesim%lasten%p8loads(ii)%cid  = cids(ii)
      fesim%lasten%p8loads(ii)%pi   = pis(1:4, ii)
      fesim%lasten%p8loads(ii)%thru = thrus(ii)
      fesim%lasten%p8loads(ii)%eid2 = eid2s(ii)
    end do
    fesim%is_p8load = .true.

  end subroutine fipps_api_p8loads

  !
  ! =================================================================================================
  !
  !> Ersetzt p20loads.fipps: Zeilenformat "lid eid1 p surf thru eid2"
  !
  subroutine fipps_api_p20loads(fesim, lids, eid1s, ps, surfs, thrus, eid2s, err)

    type(fe_simulation), intent(inout)             :: fesim
    integer,            intent(in), dimension(:)   :: lids
    integer,            intent(in), dimension(:)   :: eid1s
    double precision,   intent(in), dimension(:)   :: ps
    integer,            intent(in), dimension(:)   :: surfs
    logical,            intent(in), dimension(:)   :: thrus
    integer,            intent(in), dimension(:)   :: eid2s
    integer,            intent(out)                :: err

    integer :: n, ii

    err = 0
    n = size(lids,1)

    if (n .le. 0) then
      err = 1
      return
    end if
    if (allocated(fesim%lasten%p20loads)) then
      err = 3
      return
    end if

    allocate(fesim%lasten%p20loads(1:n))
    do ii = 1, n
      fesim%lasten%p20loads(ii)%lid  = lids(ii)
      fesim%lasten%p20loads(ii)%eid1 = eid1s(ii)
      fesim%lasten%p20loads(ii)%p    = ps(ii)
      fesim%lasten%p20loads(ii)%surf = surfs(ii)
      fesim%lasten%p20loads(ii)%thru = thrus(ii)
      fesim%lasten%p20loads(ii)%eid2 = eid2s(ii)
    end do
    fesim%is_p20load = .true.

  end subroutine fipps_api_p20loads

  !
  ! =================================================================================================
  !
  !> Ersetzt aeroload2ds.fipps: Zeilenformat "lid mthd dfac"
  !> (dfac wird wie im Datei-Einlesen auf (0,1] zurueckgesetzt)
  !
  subroutine fipps_api_aeroload2ds(fesim, lids, mthds, dfacs, err)

    type(fe_simulation), intent(inout)             :: fesim
    integer,            intent(in), dimension(:)   :: lids
    integer,            intent(in), dimension(:)   :: mthds
    double precision,   intent(in), dimension(:)   :: dfacs
    integer,            intent(out)                :: err

    integer :: n, ii

    err = 0
    n = size(lids,1)

    if (n .le. 0) then
      err = 1
      return
    end if
    if (allocated(fesim%lasten%aeroload2ds)) then
      err = 3
      return
    end if

    allocate(fesim%lasten%aeroload2ds(1:n))
    do ii = 1, n
      fesim%lasten%aeroload2ds(ii)%lid  = lids(ii)
      fesim%lasten%aeroload2ds(ii)%mthd = mthds(ii)
      fesim%lasten%aeroload2ds(ii)%dfac = dfacs(ii)
      if (fesim%lasten%aeroload2ds(ii)%dfac .le. 0.d0 .or. &
          & fesim%lasten%aeroload2ds(ii)%dfac .gt. 1.d0) then
        fesim%lasten%aeroload2ds(ii)%dfac = 1.d0
      end if
    end do
    fesim%is_aeroload2d = .true.

  end subroutine fipps_api_aeroload2ds

  !
  ! =================================================================================================
  !
  !> Ersetzt aeroload3ds.fipps: Zeilenformat "lid mthd dfac"
  !> (dfac wird wie im Datei-Einlesen auf (0,1] zurueckgesetzt)
  !
  subroutine fipps_api_aeroload3ds(fesim, lids, mthds, dfacs, err)

    type(fe_simulation), intent(inout)             :: fesim
    integer,            intent(in), dimension(:)   :: lids
    integer,            intent(in), dimension(:)   :: mthds
    double precision,   intent(in), dimension(:)   :: dfacs
    integer,            intent(out)                :: err

    integer :: n, ii

    err = 0
    n = size(lids,1)

    if (n .le. 0) then
      err = 1
      return
    end if
    if (allocated(fesim%lasten%aeroload3ds)) then
      err = 3
      return
    end if

    allocate(fesim%lasten%aeroload3ds(1:n))
    do ii = 1, n
      fesim%lasten%aeroload3ds(ii)%lid  = lids(ii)
      fesim%lasten%aeroload3ds(ii)%mthd = mthds(ii)
      fesim%lasten%aeroload3ds(ii)%dfac = dfacs(ii)
      if (fesim%lasten%aeroload3ds(ii)%dfac .le. 0.d0 .or. &
          & fesim%lasten%aeroload3ds(ii)%dfac .gt. 1.d0) then
        fesim%lasten%aeroload3ds(ii)%dfac = 1.d0
      end if
    end do
    fesim%is_aeroload3d = .true.

  end subroutine fipps_api_aeroload3ds

  !
  ! =================================================================================================
  !
  ! Materialsektionen (fid-Werte < 1 werden wie im Datei-Einlesen auf 0 gesetzt)
  !
  subroutine fipps_api_mat1s(fesim, mids, yms, sms, nus, rhos, aths, trefs, ges, fids, err)

    type(fe_simulation), intent(inout)             :: fesim
    integer,            intent(in), dimension(:)   :: mids
    double precision,   intent(in), dimension(:)   :: yms, sms, nus, rhos, aths, trefs, ges
    integer,            intent(in), dimension(:,:) :: fids   ! (4, n)
    integer,            intent(out)                :: err

    integer :: n, ii, jj

    err = 0
    n = size(mids,1)

    if (n .le. 0) then
      err = 1
      return
    end if
    if (allocated(fesim%materialien%mat1s)) then
      err = 3
      return
    end if

    allocate(fesim%materialien%mat1s(1:n))
    do ii = 1, n
      fesim%materialien%mat1s(ii)%mid  = mids(ii)
      fesim%materialien%mat1s(ii)%ym   = yms(ii)
      fesim%materialien%mat1s(ii)%sm   = sms(ii)
      fesim%materialien%mat1s(ii)%nu   = nus(ii)
      fesim%materialien%mat1s(ii)%rho  = rhos(ii)
      fesim%materialien%mat1s(ii)%ath  = aths(ii)
      fesim%materialien%mat1s(ii)%tref = trefs(ii)
      fesim%materialien%mat1s(ii)%ge   = ges(ii)
      do jj = 1, 4
        fesim%materialien%mat1s(ii)%fid(jj) = fids(jj, ii)
        if (fesim%materialien%mat1s(ii)%fid(jj) .lt. 1) &
          & fesim%materialien%mat1s(ii)%fid(jj) = 0
      end do
    end do
    fesim%is_mat1 = .true.

  end subroutine fipps_api_mat1s

  !
  ! =================================================================================================
  !
  subroutine fipps_api_mat8s(fesim, mids, ym11s, ym22s, nu12s, sm12s, sm13s, sm23s, &
                             & rhos, ath11s, ath22s, trefs, ges, fids, err)

    type(fe_simulation), intent(inout)             :: fesim
    integer,            intent(in), dimension(:)   :: mids
    double precision,   intent(in), dimension(:)   :: ym11s, ym22s, nu12s, sm12s, sm13s, sm23s
    double precision,   intent(in), dimension(:)   :: rhos, ath11s, ath22s, trefs, ges
    integer,            intent(in), dimension(:,:) :: fids   ! (4, n)
    integer,            intent(out)                :: err

    integer :: n, ii, jj

    err = 0
    n = size(mids,1)

    if (n .le. 0) then
      err = 1
      return
    end if
    if (allocated(fesim%materialien%mat8s)) then
      err = 3
      return
    end if

    allocate(fesim%materialien%mat8s(1:n))
    do ii = 1, n
      fesim%materialien%mat8s(ii)%mid   = mids(ii)
      fesim%materialien%mat8s(ii)%ym11  = ym11s(ii)
      fesim%materialien%mat8s(ii)%ym22  = ym22s(ii)
      fesim%materialien%mat8s(ii)%nu12  = nu12s(ii)
      fesim%materialien%mat8s(ii)%sm12  = sm12s(ii)
      fesim%materialien%mat8s(ii)%sm13  = sm13s(ii)
      fesim%materialien%mat8s(ii)%sm23  = sm23s(ii)
      fesim%materialien%mat8s(ii)%rho   = rhos(ii)
      fesim%materialien%mat8s(ii)%ath11 = ath11s(ii)
      fesim%materialien%mat8s(ii)%ath22 = ath22s(ii)
      fesim%materialien%mat8s(ii)%tref  = trefs(ii)
      fesim%materialien%mat8s(ii)%ge    = ges(ii)
      do jj = 1, 4
        fesim%materialien%mat8s(ii)%fid(jj) = fids(jj, ii)
        if (fesim%materialien%mat8s(ii)%fid(jj) .lt. 1) &
          & fesim%materialien%mat8s(ii)%fid(jj) = 0
      end do
    end do
    fesim%is_mat8 = .true.

  end subroutine fipps_api_mat8s

  !
  ! =================================================================================================
  !
  subroutine fipps_api_mat20s(fesim, mids, ym11s, ym22s, ym33s, nu12s, nu13s, nu23s, &
                              & sm12s, sm13s, sm23s, ath11s, ath22s, ath33s, rhos, fids, err)

    type(fe_simulation), intent(inout)             :: fesim
    integer,            intent(in), dimension(:)   :: mids
    double precision,   intent(in), dimension(:)   :: ym11s, ym22s, ym33s, nu12s, nu13s, nu23s
    double precision,   intent(in), dimension(:)   :: sm12s, sm13s, sm23s, ath11s, ath22s, ath33s, rhos
    integer,            intent(in), dimension(:,:) :: fids   ! (4, n)
    integer,            intent(out)                :: err

    integer :: n, ii, jj

    err = 0
    n = size(mids,1)

    if (n .le. 0) then
      err = 1
      return
    end if
    if (allocated(fesim%materialien%mat20s)) then
      err = 3
      return
    end if

    allocate(fesim%materialien%mat20s(1:n))
    do ii = 1, n
      fesim%materialien%mat20s(ii)%mid   = mids(ii)
      fesim%materialien%mat20s(ii)%ym11  = ym11s(ii)
      fesim%materialien%mat20s(ii)%ym22  = ym22s(ii)
      fesim%materialien%mat20s(ii)%ym33  = ym33s(ii)
      fesim%materialien%mat20s(ii)%nu12  = nu12s(ii)
      fesim%materialien%mat20s(ii)%nu13  = nu13s(ii)
      fesim%materialien%mat20s(ii)%nu23  = nu23s(ii)
      fesim%materialien%mat20s(ii)%sm12  = sm12s(ii)
      fesim%materialien%mat20s(ii)%sm13  = sm13s(ii)
      fesim%materialien%mat20s(ii)%sm23  = sm23s(ii)
      fesim%materialien%mat20s(ii)%ath11 = ath11s(ii)
      fesim%materialien%mat20s(ii)%ath22 = ath22s(ii)
      fesim%materialien%mat20s(ii)%ath33 = ath33s(ii)
      fesim%materialien%mat20s(ii)%rho   = rhos(ii)
      do jj = 1, 4
        fesim%materialien%mat20s(ii)%fid(jj) = fids(jj, ii)
        if (fesim%materialien%mat20s(ii)%fid(jj) .lt. 1) &
          & fesim%materialien%mat20s(ii)%fid(jj) = 0
      end do
    end do
    fesim%is_mat20 = .true.

  end subroutine fipps_api_mat20s

  !
  ! =================================================================================================
  !
  !> Ersetzt failure.fipps: je Versagenskriterium eine eigene Sektion
  !> (die Datei verwendet Zeilen mit Typname, z.B. "  tresca  fid  ys").
  !> Jeder Aufruf setzt is_failure.
  !
  subroutine fipps_api_failure_tresca(fesim, fids, ysv, err)

    type(fe_simulation), intent(inout)             :: fesim
    integer,            intent(in), dimension(:)   :: fids
    double precision,   intent(in), dimension(:)   :: ysv
    integer,            intent(out)                :: err

    integer :: n, ii

    err = 0
    n = size(fids,1)

    if (n .le. 0) then
      err = 1
      return
    end if

    allocate(fesim%versagenskriterien%failTrescas(1:n))
    do ii = 1, n
      fesim%versagenskriterien%failTrescas(ii)%fid = fids(ii)
      fesim%versagenskriterien%failTrescas(ii)%ys  = ysv(ii)
    end do
    fesim%is_failure = .true.

  end subroutine fipps_api_failure_tresca

  !
  ! =================================================================================================
  !
  subroutine fipps_api_failure_mises(fesim, fids, ysv, err)

    type(fe_simulation), intent(inout)             :: fesim
    integer,            intent(in), dimension(:)   :: fids
    double precision,   intent(in), dimension(:)   :: ysv
    integer,            intent(out)                :: err

    integer :: n, ii

    err = 0
    n = size(fids,1)

    if (n .le. 0) then
      err = 1
      return
    end if

    allocate(fesim%versagenskriterien%failMises(1:n))
    do ii = 1, n
      fesim%versagenskriterien%failMises(ii)%fid = fids(ii)
      fesim%versagenskriterien%failMises(ii)%ys  = ysv(ii)
    end do
    fesim%is_failure = .true.

  end subroutine fipps_api_failure_mises

  !
  ! =================================================================================================
  !
  subroutine fipps_api_failure_maxpstress(fesim, fids, ysv, yscv, err)

    type(fe_simulation), intent(inout)             :: fesim
    integer,            intent(in), dimension(:)   :: fids
    double precision,   intent(in), dimension(:)   :: ysv
    double precision,   intent(in), dimension(:)   :: yscv
    integer,            intent(out)                :: err

    integer :: n, ii

    err = 0
    n = size(fids,1)

    if (n .le. 0) then
      err = 1
      return
    end if

    allocate(fesim%versagenskriterien%failMaxprincstresses(1:n))
    do ii = 1, n
      fesim%versagenskriterien%failMaxprincstresses(ii)%fid  = fids(ii)
      fesim%versagenskriterien%failMaxprincstresses(ii)%ys   = ysv(ii)
      fesim%versagenskriterien%failMaxprincstresses(ii)%ysC  = yscv(ii)
    end do
    fesim%is_failure = .true.

  end subroutine fipps_api_failure_maxpstress

  !
  ! =================================================================================================
  !
  subroutine fipps_api_failure_puck(fesim, fids, rparTens, rparComs, rnorTens, rnorComs, &
                                    & rshears, pspds, pspzs, a0s, lambdamins, err)

    type(fe_simulation), intent(inout)             :: fesim
    integer,            intent(in), dimension(:)   :: fids
    double precision,   intent(in), dimension(:)   :: rparTens, rparComs, rnorTens, rnorComs
    double precision,   intent(in), dimension(:)   :: rshears, pspds, pspzs, a0s, lambdamins
    integer,            intent(out)                :: err

    integer :: n, ii

    err = 0
    n = size(fids,1)

    if (n .le. 0) then
      err = 1
      return
    end if

    allocate(fesim%versagenskriterien%failPucks(1:n))
    do ii = 1, n
      fesim%versagenskriterien%failPucks(ii)%fid       = fids(ii)
      fesim%versagenskriterien%failPucks(ii)%RParTen   = rparTens(ii)
      fesim%versagenskriterien%failPucks(ii)%RParCom   = rparComs(ii)
      fesim%versagenskriterien%failPucks(ii)%RNorTen   = rnorTens(ii)
      fesim%versagenskriterien%failPucks(ii)%RNorCom   = rnorComs(ii)
      fesim%versagenskriterien%failPucks(ii)%RShear    = rshears(ii)
      fesim%versagenskriterien%failPucks(ii)%Pspd      = pspds(ii)
      fesim%versagenskriterien%failPucks(ii)%Pspz      = pspzs(ii)
      fesim%versagenskriterien%failPucks(ii)%a0        = a0s(ii)
      fesim%versagenskriterien%failPucks(ii)%lambdamin = lambdamins(ii)
    end do
    fesim%is_failure = .true.

  end subroutine fipps_api_failure_puck

  !
  ! =================================================================================================
  !
  subroutine fipps_api_failure_hill(fesim, fids, rparTens, rparComs, rnorTens, rnorComs, &
                                    & rshears, f12stars, err)

    type(fe_simulation), intent(inout)             :: fesim
    integer,            intent(in), dimension(:)   :: fids
    double precision,   intent(in), dimension(:)   :: rparTens, rparComs, rnorTens, rnorComs
    double precision,   intent(in), dimension(:)   :: rshears, f12stars
    integer,            intent(out)                :: err

    integer :: n, ii

    err = 0
    n = size(fids,1)

    if (n .le. 0) then
      err = 1
      return
    end if

    allocate(fesim%versagenskriterien%failHills(1:n))
    do ii = 1, n
      fesim%versagenskriterien%failHills(ii)%fid      = fids(ii)
      fesim%versagenskriterien%failHills(ii)%RParTen  = rparTens(ii)
      fesim%versagenskriterien%failHills(ii)%RParCom  = rparComs(ii)
      fesim%versagenskriterien%failHills(ii)%RNorTen  = rnorTens(ii)
      fesim%versagenskriterien%failHills(ii)%RNorCom  = rnorComs(ii)
      fesim%versagenskriterien%failHills(ii)%RShear   = rshears(ii)
      fesim%versagenskriterien%failHills(ii)%F12star  = f12stars(ii)
    end do
    fesim%is_failure = .true.

  end subroutine fipps_api_failure_hill

  !
  ! =================================================================================================
  !
  subroutine fipps_api_failure_norris(fesim, fids, rpars, rnorms, rshears, err)

    type(fe_simulation), intent(inout)             :: fesim
    integer,            intent(in), dimension(:)   :: fids
    double precision,   intent(in), dimension(:)   :: rpars, rnorms, rshears
    integer,            intent(out)                :: err

    integer :: n, ii

    err = 0
    n = size(fids,1)

    if (n .le. 0) then
      err = 1
      return
    end if

    allocate(fesim%versagenskriterien%failNorris(1:n))
    do ii = 1, n
      fesim%versagenskriterien%failNorris(ii)%fid    = fids(ii)
      fesim%versagenskriterien%failNorris(ii)%RPar   = rpars(ii)
      fesim%versagenskriterien%failNorris(ii)%RNor   = rnorms(ii)
      fesim%versagenskriterien%failNorris(ii)%RShear = rshears(ii)
    end do
    fesim%is_failure = .true.

  end subroutine fipps_api_failure_norris

  !
  ! =================================================================================================
  !
  subroutine fipps_api_failure_fibre(fesim, fids, rparTens, rparComs, err)

    type(fe_simulation), intent(inout)             :: fesim
    integer,            intent(in), dimension(:)   :: fids
    double precision,   intent(in), dimension(:)   :: rparTens, rparComs
    integer,            intent(out)                :: err

    integer :: n, ii

    err = 0
    n = size(fids,1)

    if (n .le. 0) then
      err = 1
      return
    end if

    allocate(fesim%versagenskriterien%failFibres(1:n))
    do ii = 1, n
      fesim%versagenskriterien%failFibres(ii)%fid     = fids(ii)
      fesim%versagenskriterien%failFibres(ii)%RParTen = rparTens(ii)
      fesim%versagenskriterien%failFibres(ii)%RParCom = rparComs(ii)
    end do
    fesim%is_failure = .true.

  end subroutine fipps_api_failure_fibre

  !
  ! =================================================================================================
  !
  subroutine fipps_api_failure_maxstrain(fesim, fids, epsParTens, epsParComs, epsNorTens, &
                                         & epsNorComs, epsShears, useGlobals, err)

    type(fe_simulation), intent(inout)             :: fesim
    integer,            intent(in), dimension(:)   :: fids
    double precision,   intent(in), dimension(:)   :: epsParTens, epsParComs, epsNorTens
    double precision,   intent(in), dimension(:)   :: epsNorComs, epsShears
    logical,            intent(in), dimension(:)   :: useGlobals
    integer,            intent(out)                :: err

    integer :: n, ii

    err = 0
    n = size(fids,1)

    if (n .le. 0) then
      err = 1
      return
    end if

    allocate(fesim%versagenskriterien%failMaxStrains(1:n))
    do ii = 1, n
      fesim%versagenskriterien%failMaxStrains(ii)%fid       = fids(ii)
      fesim%versagenskriterien%failMaxStrains(ii)%epsParTen = epsParTens(ii)
      fesim%versagenskriterien%failMaxStrains(ii)%epsParCom = epsParComs(ii)
      fesim%versagenskriterien%failMaxStrains(ii)%epsNorTen = epsNorTens(ii)
      fesim%versagenskriterien%failMaxStrains(ii)%epsNorCom = epsNorComs(ii)
      fesim%versagenskriterien%failMaxStrains(ii)%epsShear  = epsShears(ii)
      fesim%versagenskriterien%failMaxStrains(ii)%useGlobal = useGlobals(ii)
    end do
    fesim%is_failure = .true.

  end subroutine fipps_api_failure_maxstrain

  !
  ! =================================================================================================
  !
  subroutine fipps_api_failure_cuntze(fesim, fids, rparTens, rparComs, rnorTens, rnorComs, &
                                      & rshears, muNorPars, ms, err)

    type(fe_simulation), intent(inout)             :: fesim
    integer,            intent(in), dimension(:)   :: fids
    double precision,   intent(in), dimension(:)   :: rparTens, rparComs, rnorTens, rnorComs
    double precision,   intent(in), dimension(:)   :: rshears, muNorPars, ms
    integer,            intent(out)                :: err

    integer :: n, ii

    err = 0
    n = size(fids,1)

    if (n .le. 0) then
      err = 1
      return
    end if

    allocate(fesim%versagenskriterien%failCuntzes(1:n))
    do ii = 1, n
      fesim%versagenskriterien%failCuntzes(ii)%fid      = fids(ii)
      fesim%versagenskriterien%failCuntzes(ii)%RParTen  = rparTens(ii)
      fesim%versagenskriterien%failCuntzes(ii)%RParCom  = rparComs(ii)
      fesim%versagenskriterien%failCuntzes(ii)%RNorTen  = rnorTens(ii)
      fesim%versagenskriterien%failCuntzes(ii)%RNorCom  = rnorComs(ii)
      fesim%versagenskriterien%failCuntzes(ii)%RShear   = rshears(ii)
      fesim%versagenskriterien%failCuntzes(ii)%muNorPar = muNorPars(ii)
      fesim%versagenskriterien%failCuntzes(ii)%m        = ms(ii)
    end do
    fesim%is_failure = .true.

  end subroutine fipps_api_failure_cuntze

  !
  ! =================================================================================================
  !
  subroutine fipps_api_failure_maxstrain3d(fesim, fids, eps11Tens, eps11Coms, eps22Tens, &
                                           & eps22Coms, eps33Tens, eps33Coms, eps12Shears, &
                                           & eps13Shears, eps23Shears, useGlobals, err)

    type(fe_simulation), intent(inout)             :: fesim
    integer,            intent(in), dimension(:)   :: fids
    double precision,   intent(in), dimension(:)   :: eps11Tens, eps11Coms, eps22Tens, eps22Coms
    double precision,   intent(in), dimension(:)   :: eps33Tens, eps33Coms, eps12Shears
    double precision,   intent(in), dimension(:)   :: eps13Shears, eps23Shears
    logical,            intent(in), dimension(:)   :: useGlobals
    integer,            intent(out)                :: err

    integer :: n, ii

    err = 0
    n = size(fids,1)

    if (n .le. 0) then
      err = 1
      return
    end if

    allocate(fesim%versagenskriterien%failMaxStrain3Ds(1:n))
    do ii = 1, n
      fesim%versagenskriterien%failMaxStrain3Ds(ii)%fid        = fids(ii)
      fesim%versagenskriterien%failMaxStrain3Ds(ii)%eps11Ten   = eps11Tens(ii)
      fesim%versagenskriterien%failMaxStrain3Ds(ii)%eps11Com   = eps11Coms(ii)
      fesim%versagenskriterien%failMaxStrain3Ds(ii)%eps22Ten   = eps22Tens(ii)
      fesim%versagenskriterien%failMaxStrain3Ds(ii)%eps22Com   = eps22Coms(ii)
      fesim%versagenskriterien%failMaxStrain3Ds(ii)%eps33Ten   = eps33Tens(ii)
      fesim%versagenskriterien%failMaxStrain3Ds(ii)%eps33Com   = eps33Coms(ii)
      fesim%versagenskriterien%failMaxStrain3Ds(ii)%eps12Shear = eps12Shears(ii)
      fesim%versagenskriterien%failMaxStrain3Ds(ii)%eps13Shear = eps13Shears(ii)
      fesim%versagenskriterien%failMaxStrain3Ds(ii)%eps23Shear = eps23Shears(ii)
      fesim%versagenskriterien%failMaxStrain3Ds(ii)%useGlobal  = useGlobals(ii)
    end do
    fesim%is_failure = .true.

  end subroutine fipps_api_failure_maxstrain3d

  !
  ! =================================================================================================
  !
  subroutine fipps_api_failure_tsaiwu3d(fesim, fids, r11Tens, r11Coms, r22Tens, r22Coms, &
                                        & r33Tens, r33Coms, r12Shears, r13Shears, r23Shears, &
                                        & coupl12s, coupl13s, coupl23s, err)

    type(fe_simulation), intent(inout)             :: fesim
    integer,            intent(in), dimension(:)   :: fids
    double precision,   intent(in), dimension(:)   :: r11Tens, r11Coms, r22Tens, r22Coms
    double precision,   intent(in), dimension(:)   :: r33Tens, r33Coms, r12Shears, r13Shears
    double precision,   intent(in), dimension(:)   :: r23Shears, coupl12s, coupl13s, coupl23s
    integer,            intent(out)                :: err

    integer :: n, ii

    err = 0
    n = size(fids,1)

    if (n .le. 0) then
      err = 1
      return
    end if

    allocate(fesim%versagenskriterien%failTsaiWu3Ds(1:n))
    do ii = 1, n
      fesim%versagenskriterien%failTsaiWu3Ds(ii)%fid       = fids(ii)
      fesim%versagenskriterien%failTsaiWu3Ds(ii)%R11Ten    = r11Tens(ii)
      fesim%versagenskriterien%failTsaiWu3Ds(ii)%R11Com    = r11Coms(ii)
      fesim%versagenskriterien%failTsaiWu3Ds(ii)%R22Ten    = r22Tens(ii)
      fesim%versagenskriterien%failTsaiWu3Ds(ii)%R22Com    = r22Coms(ii)
      fesim%versagenskriterien%failTsaiWu3Ds(ii)%R33Ten    = r33Tens(ii)
      fesim%versagenskriterien%failTsaiWu3Ds(ii)%R33Com    = r33Coms(ii)
      fesim%versagenskriterien%failTsaiWu3Ds(ii)%R12Shear  = r12Shears(ii)
      fesim%versagenskriterien%failTsaiWu3Ds(ii)%R13Shear  = r13Shears(ii)
      fesim%versagenskriterien%failTsaiWu3Ds(ii)%R23Shear  = r23Shears(ii)
      fesim%versagenskriterien%failTsaiWu3Ds(ii)%coupl12   = coupl12s(ii)
      fesim%versagenskriterien%failTsaiWu3Ds(ii)%coupl13   = coupl13s(ii)
      fesim%versagenskriterien%failTsaiWu3Ds(ii)%coupl23   = coupl23s(ii)
    end do
    fesim%is_failure = .true.

  end subroutine fipps_api_failure_tsaiwu3d

  !
  ! =================================================================================================
  !
  !> Ersetzt pshell.fipps: Zeilenformat
  !> "pid mid1 mt mid2 bmr mid3 tst nsm z1 z2 mid4"
  !
  subroutine fipps_api_pshells(fesim, pids, mid1s, mts, mid2s, bmrs, mid3s, tstes, &
                               & nsms, z1s, z2s, mid4s, err)

    type(fe_simulation), intent(inout)             :: fesim
    integer,            intent(in), dimension(:)   :: pids
    integer,            intent(in), dimension(:)   :: mid1s, mid2s, mid3s, mid4s
    double precision,   intent(in), dimension(:)   :: mts, bmrs, tstes, nsms, z1s, z2s
    integer,            intent(out)                :: err

    integer :: n, ii

    err = 0
    n = size(pids,1)

    if (n .le. 0) then
      err = 1
      return
    end if
    if (allocated(fesim%eigenschaften%pshells)) then
      err = 3
      return
    end if

    allocate(fesim%eigenschaften%pshells(1:n))
    do ii = 1, n
      fesim%eigenschaften%pshells(ii)%pid  = pids(ii)
      fesim%eigenschaften%pshells(ii)%mid1 = mid1s(ii)
      fesim%eigenschaften%pshells(ii)%mt   = mts(ii)
      fesim%eigenschaften%pshells(ii)%mid2 = mid2s(ii)
      fesim%eigenschaften%pshells(ii)%bmr  = bmrs(ii)
      fesim%eigenschaften%pshells(ii)%mid3 = mid3s(ii)
      fesim%eigenschaften%pshells(ii)%tst  = tstes(ii)
      fesim%eigenschaften%pshells(ii)%nsm  = nsms(ii)
      fesim%eigenschaften%pshells(ii)%z1   = z1s(ii)
      fesim%eigenschaften%pshells(ii)%z2   = z2s(ii)
      fesim%eigenschaften%pshells(ii)%mid4 = mid4s(ii)
    end do
    fesim%is_pshell = .true.

  end subroutine fipps_api_pshells

  !
  ! =================================================================================================
  !
  !> Ersetzt pcomp.fipps: Zeilenformat "pid lamid offset nsm sb ft"
  !> (ft: 2 Zeichen, lay wird wie im Datei-Einlesen auf 0 gesetzt und
  !> spaeter aus den lam8-Karten in init_values bestimmt)
  !
  subroutine fipps_api_pcomps(fesim, pids, lamids, offsets, nsms, sbs, fts, err)

    type(fe_simulation), intent(inout)                   :: fesim
    integer,            intent(in), dimension(:)         :: pids
    integer,            intent(in), dimension(:)         :: lamids
    double precision,   intent(in), dimension(:)         :: offsets, nsms, sbs
    character(2),       intent(in), dimension(:)         :: fts
    integer,            intent(out)                      :: err

    integer :: n, ii

    err = 0
    n = size(pids,1)

    if (n .le. 0) then
      err = 1
      return
    end if
    if (allocated(fesim%eigenschaften%pcomps)) then
      err = 3
      return
    end if

    allocate(fesim%eigenschaften%pcomps(1:n))
    do ii = 1, n
      fesim%eigenschaften%pcomps(ii)%pid    = pids(ii)
      fesim%eigenschaften%pcomps(ii)%lamid  = lamids(ii)
      fesim%eigenschaften%pcomps(ii)%offset = offsets(ii)
      fesim%eigenschaften%pcomps(ii)%nsm    = nsms(ii)
      fesim%eigenschaften%pcomps(ii)%sb     = sbs(ii)
      fesim%eigenschaften%pcomps(ii)%ft     = fts(ii)
      fesim%eigenschaften%pcomps(ii)%lay    = 0
    end do
    fesim%is_pcomp = .true.

  end subroutine fipps_api_pcomps

  !
  ! =================================================================================================
  !
  !> Ersetzt pbeam.fipps: Zeilenformat
  !> "pid mid AA I11 I22 I12 It t1 t2 angle [deg|rad] nsm"
  !
  subroutine fipps_api_pbeams(fesim, pids, mids, aas, i11s, i22s, i12s, its, t1s, t2s, &
                              & angles, angle_unit, nsms, err)

    type(fe_simulation), intent(inout)          :: fesim
    integer,            intent(in), dimension(:) :: pids
    integer,            intent(in), dimension(:) :: mids
    double precision,   intent(in), dimension(:) :: aas, i11s, i22s, i12s, its
    double precision,   intent(in), dimension(:) :: t1s, t2s
    double precision,   intent(in), dimension(:) :: angles
    character(3),       intent(in)              :: angle_unit ! 'deg' oder 'rad'
    double precision,   intent(in), dimension(:) :: nsms
    integer,            intent(out)             :: err

    integer :: n, ii

    err = 0
    n = size(pids,1)

    if (n .le. 0) then
      err = 1
      return
    end if
    if (allocated(fesim%eigenschaften%pbeams)) then
      err = 3
      return
    end if
    if ((angle_unit .ne. 'deg') .and. (angle_unit .ne. 'rad')) then
      err = 2
      return
    end if

    allocate(fesim%eigenschaften%pbeams(1:n))
    do ii = 1, n
      fesim%eigenschaften%pbeams(ii)%pid   = pids(ii)
      fesim%eigenschaften%pbeams(ii)%mid   = mids(ii)
      fesim%eigenschaften%pbeams(ii)%AA    = aas(ii)
      fesim%eigenschaften%pbeams(ii)%I11   = i11s(ii)
      fesim%eigenschaften%pbeams(ii)%I22   = i22s(ii)
      fesim%eigenschaften%pbeams(ii)%I12   = i12s(ii)
      fesim%eigenschaften%pbeams(ii)%It    = its(ii)
      fesim%eigenschaften%pbeams(ii)%t1    = t1s(ii)
      fesim%eigenschaften%pbeams(ii)%t2    = t2s(ii)
      fesim%eigenschaften%pbeams(ii)%angle = angles(ii)
      fesim%eigenschaften%pbeams(ii)%nsm   = nsms(ii)
      if (angle_unit == 'deg') then
        fesim%eigenschaften%pbeams(ii)%angle = &
          & fesim%eigenschaften%pbeams(ii)%angle/180.d0*acos(-1.d0)
      end if
    end do
    fesim%is_pbeam = .true.

  end subroutine fipps_api_pbeams

  !
  ! =================================================================================================
  !
  !> Ersetzt plsolid.fipps: Zeilenformat "pid lamid cid resLay globOut"
  !
  subroutine fipps_api_plsolids(fesim, pids, lamids, cids, resLays, globOuts, err)

    type(fe_simulation), intent(inout)             :: fesim
    integer,            intent(in), dimension(:)   :: pids
    integer,            intent(in), dimension(:)   :: lamids
    integer,            intent(in), dimension(:)   :: cids
    integer,            intent(in), dimension(:)   :: resLays
    logical,            intent(in), dimension(:)   :: globOuts
    integer,            intent(out)                :: err

    integer :: n, ii

    err = 0
    n = size(pids,1)

    if (n .le. 0) then
      err = 1
      return
    end if
    if (allocated(fesim%eigenschaften%plsolids)) then
      err = 3
      return
    end if

    allocate(fesim%eigenschaften%plsolids(1:n))
    do ii = 1, n
      fesim%eigenschaften%plsolids(ii)%pid     = pids(ii)
      fesim%eigenschaften%plsolids(ii)%lamid   = lamids(ii)
      fesim%eigenschaften%plsolids(ii)%cid     = cids(ii)
      fesim%eigenschaften%plsolids(ii)%resLay  = resLays(ii)
      fesim%eigenschaften%plsolids(ii)%globOut = globOuts(ii)
      fesim%eigenschaften%plsolids(ii)%lay     = 0
    end do
    fesim%is_plsolid = .true.

  end subroutine fipps_api_plsolids

  !
  ! =================================================================================================
  !
  !> Ersetzt spc1.fipps: Zeilenformat "sid dof n1 thru nn"
  !
  subroutine fipps_api_spc1s(fesim, sids, dofs, n1s, thrus, nns, err)

    type(fe_simulation), intent(inout)             :: fesim
    integer,            intent(in), dimension(:)   :: sids
    integer,            intent(in), dimension(:)   :: dofs
    integer,            intent(in), dimension(:)   :: n1s
    logical,            intent(in), dimension(:)   :: thrus
    integer,            intent(in), dimension(:)   :: nns
    integer,            intent(out)                :: err

    integer :: n, ii

    err = 0
    n = size(sids,1)

    if (n .le. 0) then
      err = 1
      return
    end if
    if (allocated(fesim%randbedingungen%spc1s)) then
      err = 3
      return
    end if

    allocate(fesim%randbedingungen%spc1s(1:n))
    do ii = 1, n
      fesim%randbedingungen%spc1s(ii)%sid  = sids(ii)
      fesim%randbedingungen%spc1s(ii)%dof  = dofs(ii)
      fesim%randbedingungen%spc1s(ii)%n1   = n1s(ii)
      fesim%randbedingungen%spc1s(ii)%thru = thrus(ii)
      fesim%randbedingungen%spc1s(ii)%nn   = nns(ii)
    end do
    fesim%is_spc1 = .true.

  end subroutine fipps_api_spc1s

  !
  ! =================================================================================================
  !
  !> Ersetzt mpc.fipps. Die Datei verwendet pro MPC zwei Blöcke
  !> (1 Zeile abhaengiger DOF + n_master Zeilen unabhaengiger DOFs).
  !> Im API werden die unabhaengigen DOFs als Matrix (nmax, n) uebergeben;
  !> fuer den MPC ii sind die Spalten 1:num_masters(ii) gueltig.
  !> sid wird wie im Datei-Einlesen 1..n vergeben, mpc_type = 0.
  !
  subroutine fipps_api_mpcs(fesim, num_masters, dep_nids, dep_dofs, dep_facs, &
                            & master_nids, master_dofs, master_facs, err)

    type(fe_simulation), intent(inout)                     :: fesim
    integer,            intent(in), dimension(:)           :: num_masters
    integer,            intent(in), dimension(:)           :: dep_nids
    integer,            intent(in), dimension(:)           :: dep_dofs
    double precision,   intent(in), dimension(:)           :: dep_facs
    integer,            intent(in), dimension(:,:)         :: master_nids   ! (nmax, n)
    integer,            intent(in), dimension(:,:)         :: master_dofs   ! (nmax, n)
    double precision,   intent(in), dimension(:,:)         :: master_facs   ! (nmax, n)
    integer,            intent(out)                        :: err

    integer :: n, nmax, ii, jj, nm

    err = 0
    n = size(num_masters,1)

    if (n .le. 0) then
      err = 1
      return
    end if
    if (allocated(fesim%randbedingungen%mpcs)) then
      err = 3
      return
    end if
    nmax = size(master_nids,1)
    do ii = 1, n
      if (num_masters(ii) .lt. 0 .or. num_masters(ii) .gt. nmax) then
        err = 1
        return
      end if
    end do

    allocate(fesim%randbedingungen%mpcs(1:n))
    do ii = 1, n
      fesim%randbedingungen%mpcs(ii)%sid      = ii
      fesim%randbedingungen%mpcs(ii)%mpc_type = 0
      fesim%randbedingungen%mpcs(ii)%dependend%nid  = dep_nids(ii)
      fesim%randbedingungen%mpcs(ii)%dependend%dof  = dep_dofs(ii)
      fesim%randbedingungen%mpcs(ii)%dependend%fac  = dep_facs(ii)
      nm = num_masters(ii)
      allocate(fesim%randbedingungen%mpcs(ii)%independend(1:nm))
      do jj = 1, nm
        fesim%randbedingungen%mpcs(ii)%independend(jj)%nid = master_nids(jj, ii)
        fesim%randbedingungen%mpcs(ii)%independend(jj)%dof = master_dofs(jj, ii)
        fesim%randbedingungen%mpcs(ii)%independend(jj)%fac = master_facs(jj, ii)
      end do
    end do
    fesim%is_mpc = .true.

  end subroutine fipps_api_mpcs

  !
  ! =================================================================================================
  !
  !> Ersetzt coord.fipps: Zeilenformat "cid" gefolgt von den 3 Spalten
  !> (jeweils 3 Komponenten) der Transformationssmatrix. Die Spalten werden
  !> wie im Datei-Einlesen normiert.
  !
  subroutine fipps_api_coords(fesim, cids, transMats, err)

    type(fe_simulation), intent(inout)             :: fesim
    integer,            intent(in), dimension(:)   :: cids
    double precision,   intent(in), dimension(:,:,:) :: transMats ! (3,3,n): Spalten = X-, Y-, Z-Vektor
    integer,            intent(out)                :: err

    integer :: n, ii
    double precision :: norm1, norm2, norm3

    err = 0
    n = size(cids,1)

    if (n .le. 0) then
      err = 1
      return
    end if
    if (allocated(fesim%koordinatensysteme%coords)) then
      err = 3
      return
    end if

    allocate(fesim%koordinatensysteme%coords(1:n))
    do ii = 1, n
      fesim%koordinatensysteme%coords(ii)%cid      = cids(ii)
      fesim%koordinatensysteme%coords(ii)%transMat = transMats(1:3, 1:3, ii)
      norm1 = sqrt(dot_product(transMats(1:3, 1, ii), transMats(1:3, 1, ii)))
      norm2 = sqrt(dot_product(transMats(1:3, 2, ii), transMats(1:3, 2, ii)))
      norm3 = sqrt(dot_product(transMats(1:3, 3, ii), transMats(1:3, 3, ii)))
      fesim%koordinatensysteme%coords(ii)%transMat(1:3,1) = transMats(1:3, 1, ii)/norm1
      fesim%koordinatensysteme%coords(ii)%transMat(1:3,2) = transMats(1:3, 2, ii)/norm2
      fesim%koordinatensysteme%coords(ii)%transMat(1:3,3) = transMats(1:3, 3, ii)/norm3
    end do
    fesim%is_coord = .true.

  end subroutine fipps_api_coords

  !
  ! =================================================================================================
  !
  !> Ersetzt subcase.fipps: Zeilenformat
  !> "scid spcaddid loadid mpcaddid skipBuckling upgeom upstress output readApameInput"
  !> (upmats wird wie im Datei-Einlesen von init_values bestimmt)
  !
  subroutine fipps_api_subcases(fesim, scids, spcaddids, loadids, mpcaddids, &
                                & skipBucklings, upgeoms, upstresses, outputs, &
                                & readApameInputs, err)

    type(fe_simulation), intent(inout)             :: fesim
    integer,            intent(in), dimension(:)   :: scids
    integer,            intent(in), dimension(:)   :: spcaddids
    integer,            intent(in), dimension(:)   :: loadids
    integer,            intent(in), dimension(:)   :: mpcaddids
    logical,            intent(in), dimension(:)   :: skipBucklings
    logical,            intent(in), dimension(:)   :: upgeoms
    logical,            intent(in), dimension(:)   :: upstresses
    logical,            intent(in), dimension(:)   :: outputs
    logical,            intent(in), dimension(:)   :: readApameInputs
    integer,            intent(out)                :: err

    integer :: n, ii

    err = 0
    n = size(scids,1)

    if (n .le. 0) then
      err = 1
      return
    end if
    if (allocated(fesim%lasten%subcases)) then
      err = 3
      return
    end if

    allocate(fesim%lasten%subcases(1:n))
    do ii = 1, n
      fesim%lasten%subcases(ii)%scid           = scids(ii)
      fesim%lasten%subcases(ii)%spcaddid       = spcaddids(ii)
      fesim%lasten%subcases(ii)%loadid         = loadids(ii)
      fesim%lasten%subcases(ii)%mpcaddid       = mpcaddids(ii)
      fesim%lasten%subcases(ii)%skipBuckling   = skipBucklings(ii)
      fesim%lasten%subcases(ii)%upgeom         = upgeoms(ii)
      fesim%lasten%subcases(ii)%upstress       = upstresses(ii)
      fesim%lasten%subcases(ii)%output         = outputs(ii)
      fesim%lasten%subcases(ii)%readApameInput = readApameInputs(ii)
      if (fesim%lasten%subcases(ii)%upgeom .eq. .true. .or. &
          & fesim%lasten%subcases(ii)%upstress .eq. .true.) then
        fesim%is_multistep = .true.
      end if
    end do
    fesim%is_subcase = .true.

  end subroutine fipps_api_subcases

  !
  ! =================================================================================================
  !
  !> Ersetzt lam8.fipps: Zeilenformat
  !> "lamid plyid mat8id th angle [deg|rad]"
  !
  subroutine fipps_api_lam8s(fesim, lamids, plyids, mat8ids, ths, angles, angle_unit, err)

    type(fe_simulation), intent(inout)          :: fesim
    integer,            intent(in), dimension(:) :: lamids
    integer,            intent(in), dimension(:) :: plyids
    integer,            intent(in), dimension(:) :: mat8ids
    double precision,   intent(in), dimension(:) :: ths
    double precision,   intent(in), dimension(:) :: angles
    character(3),       intent(in)              :: angle_unit ! 'deg' oder 'rad'
    integer,            intent(out)             :: err

    integer :: n, ii

    err = 0
    n = size(lamids,1)

    if (n .le. 0) then
      err = 1
      return
    end if
    if (allocated(fesim%laminate%lam8s)) then
      err = 3
      return
    end if
    if ((angle_unit .ne. 'deg') .and. (angle_unit .ne. 'rad')) then
      err = 2
      return
    end if

    allocate(fesim%laminate%lam8s(1:n))
    do ii = 1, n
      fesim%laminate%lam8s(ii)%lamid  = lamids(ii)
      fesim%laminate%lam8s(ii)%plyid  = plyids(ii)
      fesim%laminate%lam8s(ii)%mat8id = mat8ids(ii)
      fesim%laminate%lam8s(ii)%th     = ths(ii)
      fesim%laminate%lam8s(ii)%angle  = angles(ii)
      if (angle_unit == 'deg') then
        fesim%laminate%lam8s(ii)%angle = &
          & fesim%laminate%lam8s(ii)%angle/180.d0*acos(-1.d0)
      end if
    end do
    fesim%is_lam8 = .true.

  end subroutine fipps_api_lam8s

  !
  ! =================================================================================================
  !
  !> Ersetzt lam20.fipps. Die Datei erlaubt (im Multistep) pro Lage einen
  !> Material-20-Bezug pro Subcase-Zeile; das wird hier durch mat20ids(nsub, n)
  !> abgebildet: im Nicht-Multistep-Fall nsub=1, im Multistep-Fall nsub =
  !> Anzahl der bereitgestellten Subcase-Zeilen.
  !> fipps_api_subcases muss vorher aufgerufen worden sein (wie in der Datei).
  !
  subroutine fipps_api_lam20s(fesim, lamids, plyids, mat20ids, ths, angles, angle_unit, &
                              & nops, err)

    type(fe_simulation), intent(inout)          :: fesim
    integer,            intent(in), dimension(:) :: lamids
    integer,            intent(in), dimension(:) :: plyids
    integer,            intent(in), dimension(:,:) :: mat20ids ! (nsub, n)
    double precision,   intent(in), dimension(:) :: ths
    double precision,   intent(in), dimension(:) :: angles
    character(3),       intent(in)              :: angle_unit ! 'deg' oder 'rad'
    integer,            intent(in), dimension(:) :: nops
    integer,            intent(out)             :: err

    integer :: n, nsub, ii, jj

    err = 0
    n = size(lamids,1)

    if (n .le. 0) then
      err = 1
      return
    end if
    if (allocated(fesim%laminate%lam20s)) then
      err = 3
      return
    end if
    if ((angle_unit .ne. 'deg') .and. (angle_unit .ne. 'rad')) then
      err = 2
      return
    end if
    if (fesim%is_multistep .eqv. .true.) then
      if (.not. allocated(fesim%lasten%subcases)) then
        err = 4
        return
      end if
      nsub = size(fesim%lasten%subcases,1)
    else
      nsub = 1
    end if
    if (size(mat20ids,1) .lt. nsub) then
      err = 1
      return
    end if

    allocate(fesim%laminate%lam20s(1:n))
    do ii = 1, n
      fesim%laminate%lam20s(ii)%lamid = lamids(ii)
      fesim%laminate%lam20s(ii)%plyid = plyids(ii)
      allocate(fesim%laminate%lam20s(ii)%mat20id(1:nsub))
      do jj = 1, nsub
        fesim%laminate%lam20s(ii)%mat20id(jj) = mat20ids(jj, ii)
      end do
      fesim%laminate%lam20s(ii)%th    = ths(ii)
      fesim%laminate%lam20s(ii)%angle = angles(ii)
      fesim%laminate%lam20s(ii)%nop   = nops(ii)
      if (angle_unit == 'deg') then
        fesim%laminate%lam20s(ii)%angle = &
          & fesim%laminate%lam20s(ii)%angle/180.d0*acos(-1.d0)
      end if
    end do
    fesim%is_lam20 = .true.

  end subroutine fipps_api_lam20s

  !
  ! =================================================================================================
  !
  !> Ersetzt couplings.fipps: Zeilenformat "cpsid dof nid"
  !
  subroutine fipps_api_couplings(fesim, cpsids, dofs, nids, err)

    type(fe_simulation), intent(inout)             :: fesim
    integer,            intent(in), dimension(:)   :: cpsids
    integer,            intent(in), dimension(:)   :: dofs
    integer,            intent(in), dimension(:)   :: nids
    integer,            intent(out)                :: err

    integer :: n, ii

    err = 0
    n = size(cpsids,1)

    if (n .le. 0) then
      err = 1
      return
    end if
    if (allocated(fesim%randbedingungen%couplings)) then
      err = 3
      return
    end if

    allocate(fesim%randbedingungen%couplings(1:n))
    do ii = 1, n
      fesim%randbedingungen%couplings(ii)%cpsid = cpsids(ii)
      fesim%randbedingungen%couplings(ii)%dof   = dofs(ii)
      fesim%randbedingungen%couplings(ii)%nid   = nids(ii)
    end do
    fesim%is_coupling = .true.

  end subroutine fipps_api_couplings

  !
  ! =================================================================================================
  !
  !> Ersetzt contact_node_beam2.fipps: Zeilenformat "beam2ID nodeID xi".
  !> Die Zeilen werden wie im Datei-Einlesen nach Elementen zusammengefasst
  !> (Elementreihenfolge = Reihenfolge des ersten Auftretens).
  !
  subroutine fipps_api_contact_node_beam2(fesim, eids, nids, xis, err)

    type(fe_simulation), intent(inout)             :: fesim
    integer,            intent(in), dimension(:)   :: eids
    integer,            intent(in), dimension(:)   :: nids
    double precision,   intent(in), dimension(:)   :: xis
    integer,            intent(out)                :: err

    logical,    allocatable, dimension(:)          :: seen
    integer     :: n, ii, jj, jj2, conNum, anz, nelem

    err = 0
    n = size(eids,1)

    if (n .le. 0) then
      err = 1
      return
    end if
    if (allocated(fesim%randbedingungen%contact_node_beam2)) then
      err = 3
      return
    end if
    if (.not. allocated(fesim%elemente%beam2s)) then
      err = 4
      return
    end if
    nelem = size(fesim%elemente%beam2s,1)

    allocate(seen(1:nelem))
    seen = .false.

    do ii = 1, n
      if (eids(ii) .lt. 1 .or. eids(ii) .gt. nelem) then
        deallocate(seen)
        err = 1
        return
      end if
      if (.not. seen(eids(ii))) then
        seen(eids(ii)) = .true.
      end if
    end do

    conNum = 0
    do ii = 1, nelem
      if (seen(ii)) conNum = conNum + 1
    end do
    if (conNum .eq. 0) then
      deallocate(seen)
      err = 1
      return
    end if
    allocate(fesim%randbedingungen%contact_node_beam2(1:conNum))
    conNum = 0
    do ii = 1, nelem
      if (.not. seen(ii)) cycle
      anz = 0
      do jj = 1, n
        if (eids(jj) .eq. ii) anz = anz + 1
      end do
      conNum = conNum + 1
      fesim%randbedingungen%contact_node_beam2(conNum)%beam2ID = ii
      allocate(fesim%randbedingungen%contact_node_beam2(conNum)%nodeIDs(anz))
      allocate(fesim%randbedingungen%contact_node_beam2(conNum)%xi(anz))
      jj = 0
      do jj2 = 1, n
        if (eids(jj2) .eq. ii) then
          jj = jj + 1
          fesim%randbedingungen%contact_node_beam2(conNum)%nodeIDs(jj) = nids(jj2)
          fesim%randbedingungen%contact_node_beam2(conNum)%xi(jj)      = xis(jj2)
        end if
      end do
    end do

    deallocate(seen)
    fesim%is_contact_node_beam2 = .true.

  end subroutine fipps_api_contact_node_beam2

  !
  ! =================================================================================================
  !
  !> Ersetzt contact_node_quad8.fipps: Zeilenformat "quad8ID nodeID".
  !> Die Zeilen werden wie im Datei-Einlesen nach Elementen zusammengefasst.
  !
  subroutine fipps_api_contact_node_quad8(fesim, eids, nids, err)

    type(fe_simulation), intent(inout)             :: fesim
    integer,            intent(in), dimension(:)   :: eids
    integer,            intent(in), dimension(:)   :: nids
    integer,            intent(out)                :: err

    logical,    allocatable, dimension(:)          :: seen
    integer     :: n, ii, jj, jj2, conNum, anz, nelem

    err = 0
    n = size(eids,1)

    if (n .le. 0) then
      err = 1
      return
    end if
    if (allocated(fesim%randbedingungen%contact_node_quad8)) then
      err = 3
      return
    end if
    if (.not. allocated(fesim%elemente%quad8s)) then
      err = 4
      return
    end if
    nelem = size(fesim%elemente%quad8s,1)

    allocate(seen(1:nelem))
    seen = .false.

    do ii = 1, n
      if (eids(ii) .lt. 1 .or. eids(ii) .gt. nelem) then
        deallocate(seen)
        err = 1
        return
      end if
      seen(eids(ii)) = .true.
    end do

    conNum = 0
    do ii = 1, nelem
      if (seen(ii)) conNum = conNum + 1
    end do
    if (conNum .eq. 0) then
      deallocate(seen)
      err = 1
      return
    end if
    allocate(fesim%randbedingungen%contact_node_quad8(1:conNum))
    conNum = 0
    do ii = 1, nelem
      if (.not. seen(ii)) cycle
      anz = 0
      do jj = 1, n
        if (eids(jj) .eq. ii) anz = anz + 1
      end do
      conNum = conNum + 1
      fesim%randbedingungen%contact_node_quad8(conNum)%quad8ID = ii
      allocate(fesim%randbedingungen%contact_node_quad8(conNum)%nodeIDs(anz))
      jj = 0
      do jj2 = 1, n
        if (eids(jj2) .eq. ii) then
          jj = jj + 1
          fesim%randbedingungen%contact_node_quad8(conNum)%nodeIDs(jj) = nids(jj2)
        end if
      end do
    end do

    deallocate(seen)
    fesim%is_contact_node_quad8 = .true.

  end subroutine fipps_api_contact_node_quad8

  !
  ! =================================================================================================
  !
  !> Ersetzt contact_node_lsolid20.fipps: Zeilenformat "lsolid20ID nodeID".
  !> Die Zeilen werden wie im Datei-Einlesen nach Elementen zusammengefasst.
  !
  subroutine fipps_api_contact_node_lsolid20(fesim, eids, nids, err)

    type(fe_simulation), intent(inout)             :: fesim
    integer,            intent(in), dimension(:)   :: eids
    integer,            intent(in), dimension(:)   :: nids
    integer,            intent(out)                :: err

    logical,    allocatable, dimension(:)          :: seen
    integer     :: n, ii, jj, jj2, conNum, anz, nelem

    err = 0
    n = size(eids,1)

    if (n .le. 0) then
      err = 1
      return
    end if
    if (allocated(fesim%randbedingungen%contact_node_lsolid20)) then
      err = 3
      return
    end if
    if (.not. allocated(fesim%elemente%lsolid20s)) then
      err = 4
      return
    end if
    nelem = size(fesim%elemente%lsolid20s,1)

    allocate(seen(1:nelem))
    seen = .false.

    do ii = 1, n
      if (eids(ii) .lt. 1 .or. eids(ii) .gt. nelem) then
        deallocate(seen)
        err = 1
        return
      end if
      seen(eids(ii)) = .true.
    end do

    conNum = 0
    do ii = 1, nelem
      if (seen(ii)) conNum = conNum + 1
    end do
    if (conNum .eq. 0) then
      deallocate(seen)
      err = 1
      return
    end if
    allocate(fesim%randbedingungen%contact_node_lsolid20(1:conNum))
    conNum = 0
    do ii = 1, nelem
      if (.not. seen(ii)) cycle
      anz = 0
      do jj = 1, n
        if (eids(jj) .eq. ii) anz = anz + 1
      end do
      conNum = conNum + 1
      fesim%randbedingungen%contact_node_lsolid20(conNum)%lsolid20ID = ii
      allocate(fesim%randbedingungen%contact_node_lsolid20(conNum)%nodeIDs(anz))
      jj = 0
      do jj2 = 1, n
        if (eids(jj2) .eq. ii) then
          jj = jj + 1
          fesim%randbedingungen%contact_node_lsolid20(conNum)%nodeIDs(jj) = nids(jj2)
        end if
      end do
    end do

    deallocate(seen)
    fesim%is_contact_node_lsolid20 = .true.

  end subroutine fipps_api_contact_node_lsolid20

  !
  ! =================================================================================================
  !
  !> Ersetzt aeroelem2structnode2d.fipps und structelem2aeronode2d.fipps.
  !> 1. Block (Datei aeroelem2structnode2d.fipps): nodeID, elemID, xi
  !> (Knoten am Aero-Element), danach die betroffenen Struktur-Beam2-Elemente.
  !> 2. Block (Datei structelem2aeronode2d.fipps): nodeID, elemID, xi
  !> (Knoten am Struktur-Element).
  !
  subroutine fipps_api_aero_coupling_2d(fesim, &
                                        & aeroNodeIDs, aeroElemIDs, aeroXis, &
                                        & structBeam2IDs, &
                                        & structNodeIDs, structElemIDs, structXis, err)

    type(fe_simulation), intent(inout)             :: fesim
    integer,            intent(in), dimension(:)   :: aeroNodeIDs
    integer,            intent(in), dimension(:)   :: aeroElemIDs
    double precision,   intent(in), dimension(:)   :: aeroXis
    integer,            intent(in), dimension(:)   :: structBeam2IDs
    integer,            intent(in), dimension(:)   :: structNodeIDs
    integer,            intent(in), dimension(:)   :: structElemIDs
    double precision,   intent(in), dimension(:)   :: structXis
    integer,            intent(out)                :: err

    integer :: n1, m, n3, ii

    err = 0
    n1 = size(aeroNodeIDs,1)
    m  = size(structBeam2IDs,1)
    n3 = size(structNodeIDs,1)

    if (n1 .le. 0 .or. m .le. 0 .or. n3 .le. 0) then
      err = 1
      return
    end if
    if (allocated(fesim%internals%aeroElem2structNode)) then
      err = 3
      return
    end if

    allocate(fesim%internals%aeroElem2structNode(1:n1))
    do ii = 1, n1
      fesim%internals%aeroElem2structNode(ii)%nodeID = aeroNodeIDs(ii)
      fesim%internals%aeroElem2structNode(ii)%elemID = aeroElemIDs(ii)
      fesim%internals%aeroElem2structNode(ii)%xi     = aeroXis(ii)
    end do

    allocate(fesim%internals%structBeam2IDs(1:m))
    fesim%internals%structBeam2IDs(:) = structBeam2IDs(:)

    allocate(fesim%internals%structElem2aeroNode(1:n3))
    do ii = 1, n3
      fesim%internals%structElem2aeroNode(ii)%nodeID = structNodeIDs(ii)
      fesim%internals%structElem2aeroNode(ii)%elemID = structElemIDs(ii)
      fesim%internals%structElem2aeroNode(ii)%xi     = structXis(ii)
    end do

  end subroutine fipps_api_aero_coupling_2d

  !
  ! =================================================================================================
  !
  !> Ersetzt aeroelem2structnode3d.fipps und structelem2aeronode3d.fipps.
  !> Wie fipps_api_aero_coupling_2d, zusaezlich eta-Koordinaten; betroffene
  !> Struktur-Quad8-Elemente.
  !
  subroutine fipps_api_aero_coupling_3d(fesim, &
                                        & aeroNodeIDs, aeroElemIDs, aeroXis, aeroEtas, &
                                        & structQuad8IDs, &
                                        & structNodeIDs, structElemIDs, structXis, structEtas, err)

    type(fe_simulation), intent(inout)             :: fesim
    integer,            intent(in), dimension(:)   :: aeroNodeIDs
    integer,            intent(in), dimension(:)   :: aeroElemIDs
    double precision,   intent(in), dimension(:)   :: aeroXis
    double precision,   intent(in), dimension(:)   :: aeroEtas
    integer,            intent(in), dimension(:)   :: structQuad8IDs
    integer,            intent(in), dimension(:)   :: structNodeIDs
    integer,            intent(in), dimension(:)   :: structElemIDs
    double precision,   intent(in), dimension(:)   :: structXis
    double precision,   intent(in), dimension(:)   :: structEtas
    integer,            intent(out)                :: err

    integer :: n1, m, n3, ii

    err = 0
    n1 = size(aeroNodeIDs,1)
    m  = size(structQuad8IDs,1)
    n3 = size(structNodeIDs,1)

    if (n1 .le. 0 .or. m .le. 0 .or. n3 .le. 0) then
      err = 1
      return
    end if
    if (allocated(fesim%internals%aeroElem2structNode)) then
      err = 3
      return
    end if

    allocate(fesim%internals%aeroElem2structNode(1:n1))
    do ii = 1, n1
      fesim%internals%aeroElem2structNode(ii)%nodeID = aeroNodeIDs(ii)
      fesim%internals%aeroElem2structNode(ii)%elemID = aeroElemIDs(ii)
      fesim%internals%aeroElem2structNode(ii)%xi     = aeroXis(ii)
      fesim%internals%aeroElem2structNode(ii)%eta    = aeroEtas(ii)
    end do

    allocate(fesim%internals%structQuad8IDs(1:m))
    fesim%internals%structQuad8IDs(:) = structQuad8IDs(:)

    allocate(fesim%internals%structElem2aeroNode(1:n3))
    do ii = 1, n3
      fesim%internals%structElem2aeroNode(ii)%nodeID = structNodeIDs(ii)
      fesim%internals%structElem2aeroNode(ii)%elemID = structElemIDs(ii)
      fesim%internals%structElem2aeroNode(ii)%xi     = structXis(ii)
      fesim%internals%structElem2aeroNode(ii)%eta    = structEtas(ii)
    end do

  end subroutine fipps_api_aero_coupling_3d

  !
  ! =================================================================================================
  !
  !> Konsistenzpruefung, die input_tf nach dem Einlesen der control.fipps-
  !> Woerter durchfuehrt (Ausschlussbedingungen und Pflichtsections).
  !> err = 0, wenn alles konsistent ist.
  !
  subroutine fipps_api_check(fesim, err)

    type(fe_simulation), intent(in) :: fesim
    integer,            intent(out) :: err

    err = 0

    if ((fesim%is_aeroload2d .eq. .true.) .and. (fesim%is_aeroload3d .eq. .true.)) then
      write(*,*) "Keywords 'aeroload2d' and 'aeroload3d' may not"
      write(*,*) "be used in combination."
      err = 1
      return
    end if

    if ((fesim%is_temperature .eq. .true.) .and. (fesim%is_beam2temp .eq. .true.)) then
      write(*,*) "Keywords 'temperature' and 'beam2temp' may not"
      write(*,*) "be used in combination."
      err = 1
      return
    end if

    if ((fesim%is_temperature .eq. .true.) .and. (fesim%is_quad8temp .eq. .true.)) then
      write(*,*) "Keywords 'temperature' and 'quad8temp' may not"
      write(*,*) "be used in combination."
      err = 1
      return
    end if

    if ((fesim%is_temperature .eq. .true.) .and. (fesim%is_lsolid20temp .eq. .true.)) then
      write(*,*) "Keywords 'temperature' and 'lsolid20temp' may not"
      write(*,*) "be used in combination."
      err = 1
      return
    end if

    if (fesim%is_node .eqv. .false.) then
      write(*,*) 'No nodes defined (fipps_api_nodes not called).'
      err = 1
      return
    end if

    if (fesim%is_load .eqv. .false.) then
      write(*,*) 'No loads defined (fipps_api_loads not called).'
      err = 1
      return
    end if

  end subroutine fipps_api_check

end module fipps_api

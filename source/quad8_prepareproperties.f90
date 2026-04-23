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
!  FiPPS² is free software: you can redistribute it and/or modify
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
!> Preperation of the 8-node shell element properties
!
!> @details
!> Preperation of the 8-node shell element properties to avoid computational effort during the
!> calculation of needed elemental data (striffness matrix a.s.o.)
!
!> @author Andreas Hauffe, TU Dresden, wiss. Mitarbeiter, 29.06.2010
!
!> $Id: quad8_prepareproperties.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
subroutine quad8_prepareProperties(fesim)
! =================================================================================================
!
! Use
!
  use konstanten
  use fesimulation_typen
!
! =================================================================================================
!
  implicit none
!
! =================================================================================================
!
! Include
!

!
! =================================================================================================
!
! Data types
!
! Input
!
! Output
!
!
! Input + Output
!
  type(fe_simulation) :: fesim
!
! Internal
!
  integer :: elem
  integer :: ii, kk, ll
  integer :: propid, propm1
  integer :: int_pid, propType
  integer :: err_code = 0
  double precision, dimension(8,3) :: node_coords
!
! =================================================================================================
!
! Initialisation
!
  propm1   = 0
  propType = 0
  int_pid  = 0
!
! =================================================================================================
!
! Calculation
!

 do ii = 1,size(fesim%eigenschaften%pcomps,1)
   ! allocate arrays dependent on number of layers
   allocate(fesim%eigenschaften%pcomps(ii)%C(6,6,fesim%eigenschaften%pcomps(ii)%lay))
   allocate(fesim%eigenschaften%pcomps(ii)%T(6,6,fesim%eigenschaften%pcomps(ii)%lay))
   allocate(fesim%eigenschaften%pcomps(ii)%TInv(6,6,fesim%eigenschaften%pcomps(ii)%lay))
   allocate(fesim%eigenschaften%pcomps(ii)%lth(fesim%eigenschaften%pcomps(ii)%lay))
   allocate(fesim%eigenschaften%pcomps(ii)%ath(6,fesim%eigenschaften%pcomps(ii)%lay))
   allocate(fesim%eigenschaften%pcomps(ii)%nop(fesim%eigenschaften%pcomps(ii)%lay))
   allocate(fesim%eigenschaften%pcomps(ii)%angle(fesim%eigenschaften%pcomps(ii)%lay))
   allocate(fesim%eigenschaften%pcomps(ii)%intMatID(fesim%eigenschaften%pcomps(ii)%lay))
 end do


  ! Loop over all flat (reduced) quad8 elements
  do elem=1,size(fesim%elemente%quad8s,1)

    propid = fesim%elemente%quad8s(elem)%pid

    if (propid /= propm1) then

      propType = 0
      int_pid  = 0

      ! isotropic material

      if (fesim%is_pshell == .true.) then
        do kk=1,size(fesim%eigenschaften%pshells,1)

          if (propid == fesim%eigenschaften%pshells(kk)%pid) then

            propType = 1
            int_pid  = kk

            do ll=1,size(fesim%materialien%mat1s,1)
              if (fesim%eigenschaften%pshells(kk)%mid1 == fesim%materialien%mat1s(ll)%mid) then
                fesim%eigenschaften%pshells(kk)%intMat1ID = ll
                fesim%eigenschaften%pshells(kk)%areaWeight = fesim%eigenschaften%pshells(kk)%mt*fesim%materialien%mat1s(ll)%rho
                call mat1_calc_missing_value(fesim,ll,fesim%materialien%mat1s(ll)%ym,fesim%materialien%mat1s(ll)%sm,fesim%materialien%mat1s(ll)%nu)
                exit
              end if
            end do

            fesim%eigenschaften%pshells(kk)%weight = 0.d0

          end if
        end do
      end if

      ! layerwise othotropic material

      if (fesim%is_pcomp == .true.) then
        do kk=1,size(fesim%eigenschaften%pcomps,1)

          if (propid == fesim%eigenschaften%pcomps(kk)%pid) then

            ! Wenn zwei Properties gefunden wurden
            if (propType .NE. 0) then
              write(*,*) 'More than one property found for Quad8 element ', fesim%elemente%quad8s(elem)%eid
              err_code = 1
              goto 9999
            end if

            propType = 2
            int_pid  = kk

            call mat8(fesim,kk,fesim%eigenschaften%pcomps(kk)%lay,fesim%eigenschaften%pcomps(kk)%C,  &
                fesim%eigenschaften%pcomps(kk)%T,fesim%eigenschaften%pcomps(kk)%TInv,fesim%eigenschaften%pcomps(kk)%lth,  &
                fesim%eigenschaften%pcomps(kk)%ath,fesim%eigenschaften%pcomps(kk)%nop,fesim%eigenschaften%pcomps(kk)%angle,  &
                fesim%eigenschaften%pcomps(kk)%tth,fesim%eigenschaften%pcomps(kk)%penaltyStiffness,  &
                fesim%eigenschaften%pcomps(kk)%intMatID,fesim%eigenschaften%pcomps(kk)%areaWeight)

            fesim%eigenschaften%pcomps(kk)%weight = 0.d0

          end if
        end do
      end if

      ! Wenn keine Property gefunden wurde
      if (propType .EQ. 0) then
        write(*,*) 'No Property found for Quad8 element ', fesim%elemente%quad8s(elem)%eid
        err_code = 2
        goto 9999
      end if

      propm1 = propid

    end if

    fesim%elemente%quad8s(elem)%propType = propType
    fesim%elemente%quad8s(elem)%int_pid  = int_pid


    do kk = 1,8
      node_coords(kk,1:3) = fesim%knoten%nodes(fesim%elemente%quad8s(elem)%nids(kk))%coords(1:3)
    end do

    call quad8_area(node_coords,fesim%elemente%quad8s(elem)%area)

  end do

!
! =================================================================================================
!
! Error handling
!
9999 continue

  if (err_code /= 0) then

    write(*,*)                      'An error occured in subroutine'
    write(*,*)                      'quad8_prepareProperties'
    write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
    write(*,*)                      'exit program '
    stop

  end if

end subroutine quad8_prepareProperties

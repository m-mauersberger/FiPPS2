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
!> @file fipps_main.f90
!> Hauptprogramm (Vorschlag: schlankes Programm mit Bibliothek-tauglichen
!> Subroutinen fipps_initialize, fipps_subcase und fipps_finalize)
!
!> @author  Martin Rädel, TU Dresden, Diplomarbeit
!> @date  30.06.2010
!
! =================================================================================================
PROGRAM FiPPS
!
  use globale_variablen, ONLY : shortoutFile
  use fesimulation_typen
  use pre_assemble_types
  use konstanten
!
!
#include "petsc/finclude/petscksp.h"
#include "slepc/finclude/slepceps.h"

  use petscksp
  use slepceps

IMPLICIT NONE

#if defined (PETSC_HAVE_MPIUNI)
    MPIUNI_FInt MPI_INTEGER,MPI_LOGICAL
#endif

! Interfaces
interface
    subroutine aero_coupling(Ausgabe, fesim, scloop, cdi)
        use fesimulation_typen
        implicit none
        logical, intent(in)                                      :: Ausgabe
        type(fe_simulation), intent(in)                          :: fesim
        integer, intent(in)                                      :: scloop
        double precision, allocatable, dimension(:), intent(out) :: cdi
    end subroutine aero_coupling

    subroutine readpressure(fesim, pAeroElem, F_total_aero)
        use fesimulation_typen
        implicit none
        type(fe_simulation)                                 :: fesim
        double precision, dimension(:,:), intent(in)        :: pAeroElem
        double precision, dimension(:), intent(in)          :: F_total_aero
    end subroutine readpressure
    
    subroutine sol1_output_glawi2(fesim,Utot,scloop,cdi,aeroConverged)
        use fesimulation_typen
        implicit none
        type(fe_simulation)                                 :: fesim
        double precision, dimension(fesim%num_dof)          :: Utot
        integer,intent(in)                                  :: scloop
        double precision, dimension(:)                      :: cdi
        logical                                             :: aeroConverged
    end subroutine sol1_output_glawi2
    
    subroutine sol1_output_hymowi2(fesim,Utot,scloop,aeroConverged)
        use fesimulation_typen
        implicit none  
        type(fe_simulation), intent(in)                         :: fesim
        double precision, dimension(fesim%num_dof), intent(in)  :: Utot
        integer,intent(in)                                      :: scloop
        logical, intent(in)                                     :: aeroConverged
    end subroutine sol1_output_hymowi2
    
    subroutine sol1_output_control(fesim,Utot,scloop,cdi,aeroConverged)
        use fesimulation_typen
        implicit none
        type(fe_simulation), intent(inout)                      :: fesim
        double precision, dimension(fesim%num_dof), intent(in)  :: Utot
        integer, intent(in)                                     :: scloop
        double precision, dimension(:), allocatable, intent(in) :: cdi
        logical, intent(in)                                     :: aeroConverged
    end subroutine sol1_output_control
    
    subroutine sol1_reaction_forces(fesim, Utot, Ftot, rank)
        use fesimulation_typen
        implicit none
        type(fe_simulation)                                     :: fesim
        double precision, dimension(:), allocatable, intent(in) :: Utot
        double precision, dimension(:), allocatable, intent(in) :: Ftot
        PetscMPIInt, intent(in)                                 :: rank
    end subroutine sol1_reaction_forces
end interface

#include "../include/version.include"

  double precision, dimension(:), allocatable           :: Utot         ! Lösungsvektor mit allen Freiheiten
  integer                                               :: start, ende, ctr, ctm ! Variablen zur Zeitmessung
  integer                                               :: startC, endeC, ctrC, ctmC
  integer                                               :: startS, endeS, ctrS, ctmS
  integer                                               :: startI1, endeI1, ctrI1, ctmI1
  integer                                               :: startI2, endeI2, ctrI2, ctmI2
  integer                                               :: startL, endeL, ctrL, ctmL
  double precision                                      :: compTime
  double precision                                      :: solvTime
  double precision                                      :: solvEigTime
  double precision                                      :: stiffTime
  double precision                                      :: initTime
  double precision                                      :: lvecTime
  double precision, dimension(:), allocatable           :: eigvals
  double precision, dimension(:,:), allocatable         :: eigenvectors
  character(18)                                         :: fname
! !
! !
  Mat                                                   :: KaaS         ! Steifigkeitsmatrix ohne gesperrte Freiheiten als Sparse-Matrix
  Mat                                                   :: KgaaS        ! geometrische Steifigkeitsmatrix ohne gesperrte Freiheiten als Sparse-Matrix
  Vec                                                   :: UaS          ! Verschiebungsvektor als MPI-Sparse-Vektor
  Vec                                                   :: UaSseq       ! Verschiebungsvektor als Seq-Sparse-Vektor zur Ausgabe und Berechnung der Loesung
  Vec                                                   :: FaS          ! Kraftvektor ohne gesperrte Freiheiten als Sparse-Vektor
! !      
  EPS                                                   :: eps
  KSP                                                   :: ksp
  PetscErrorCode                                        :: ierr
  PetscMPIInt                                           :: rank, mpisize
  PetscInt                                              :: nconv
  VecScatter                                            :: scatter
  PetscInt                                              :: major
  PetscInt                                              :: minor
  PetscInt                                              :: subminor
  
  type(fe_simulation)                                   :: fesim

  integer                                               :: scloop

  double precision, allocatable, dimension(:)           :: cdi
  double precision                                      :: displ, maxDispl, oldMaxDispl
  integer                                               :: ii
  logical                                               :: aeroConverged = .false.
  logical                                               :: moved = .false.
  character(3)                                          :: majorChar
  character(3)                                          :: minorChar
  character(3)                                          :: subminorChar
  character(198)                                        :: buf

  double precision, dimension(:), allocatable           :: Fout         ! Kraftvektor mit allen Freiheiten, wobei alle Kraftewerte im globalen System sind
!
! =================================================================================================
!
! Driver: initialisierung, Multistepschleife, abschliessende Rueckraeumung
!
! Die enthaltenen Subroutinen fipps_initialize, fipps_subcase und
! fipps_finalize sind Kandidaten, in einem nachfolgenden Schritt in eigene
! Quelldateien (und damit in die Bibliothek libfipps.a) ausgegliedert zu werden.
!
  call fipps_initialize()

  do scloop = 1,fesim%num_subcases
      call fipps_subcase(scloop)
  end do

  call fipps_finalize()

contains

  !
  ! =================================================================================================
  !
  ! Initialisation (SLEPc/MPI, Banner, Modell lesen, Initialisierung, Aufteilung)
  !
  subroutine fipps_initialize()
!
    call slepcInitialize(PETSC_NULL_CHARACTER,ierr); CHKERRA(ierr)
    call MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr); CHKERRA(ierr)
    call MPI_Comm_size(PETSC_COMM_WORLD,mpisize,ierr)

    if (rank == 0) then
      write(*,*) ''
      write(*,*) ''
      write(*,'(a)') '************************************************************************************************************************'
      write(*,'(a)') '*'
      ! FiPPS2
      write(*,'(a)') '*  FiPPS2 ' // version
      ! PETSc
      call PetscGetVersionNumber(major,minor,subminor,PETSC_NULL_INTEGER,ierr); CHKERRA(ierr)
      write(majorChar,'(I3)') major
      write(minorChar,'(I3)') minor
      write(subminorChar,'(I3)') subminor
      write(*,'(a)') '*  PETSc  ' // trim(adjustl(majorChar)) // '.' // trim(adjustl(minorChar)) // '.' // trim(adjustl(subminorChar))
      ! SLEPc
      call SlepcGetVersionNumber(major,minor,subminor,PETSC_NULL_INTEGER,ierr); CHKERRA(ierr)
      write(majorChar,'(I3)') major
      write(minorChar,'(I3)') minor
      write(subminorChar,'(I3)') subminor
      write(*,'(a)') '*  SLEPc  ' // trim(adjustl(majorChar)) // '.' // trim(adjustl(minorChar)) // '.' // trim(adjustl(subminorChar))
      ! MKL
      call mkl_get_version_string(buf)
      write(*,'(a)') '*  ' // TRIM(buf)
      write(*,'(a)') '*'
      write(*,'(a)') '************************************************************************************************************************'
      write(*,*) ''
      write(*,*) ''
    end if  
  
    if (rank == 0) then
        compTime = 0.d0
        solvTime = 0.d0
        solvEigTime = 0.d0
        stiffTime = 0.d0
        initTime = 0.d0
        lvecTime = 0.d0
        maxDispl = 0.d0
    end if
    
    ! Initialisieren der Defaultwerte
    call init_default(fesim)

    !
    ! =================================================================================================
    !
    ! Einlesen der Dateien des Problems
    !
    if (rank == 0) then
    
      if (textoutput .eq. .true.) write(*,*) 'START - Lese Modell'
      call input_tf (fesim)                                                       ! Einlesen der Inputdateien
      if (textoutput .eq. .true.) write(*,*) 'ENDE  - Lese Modell'

      call SYSTEM_CLOCK(startI1,ctrI1,ctmI1)
          
      if (textoutput .eq. .true.) write(*,*) 'START - Initialisiere Variablen'
      call init_values (fesim)
      if (textoutput .eq. .true.) write(*,*) 'ENDE  - Initialisiere Variablen'

      ! Der Aufruf der internen MPCs sollte vor dem Herausschreiben des Netzes
      ! erfolgen, damit eventuelle Anpassungen an z.B. Kontenkoordinaten entahlten sind
      if (textoutput .eq. .true.) write(*,*) 'START - Erzeuge interne MPCs (z.B. Kontakte)'
      call impc_control(fesim)                                                                              
      if (textoutput .eq. .true.) write(*,*) 'ENDE  - Erzeuge interne MPCs (z.B. Kontakte)'

      call SYSTEM_CLOCK(endeI1,ctrI1,ctmI1)
      initTime  = initTime + (endeI1-startI1)/DBLE(ctrI1)
      
    end if

    if (mpisize .GT. 1) then
        call bcast_fesim(fesim)
    end if
    
    if (textoutput .eq. .true. .and. rank .eq. 0) write(*,*) 'START - Bestimme MPI-Interna'
    call init_elem_proc_dist(fesim, mpisize)
    if (textoutput .eq. .true. .and. rank .eq. 0) write(*,*) 'ENDE  - Bestimme MPI-Interna'

  end subroutine fipps_initialize

  !
  ! =================================================================================================
  !
  ! Ein Multistep-Durchlauf: Matrizen, Aero-Kopplung, Last, statische Loesung,
  ! Ausgabe, Biege-Eigenwertproblem
  !
  subroutine fipps_subcase(scloop)
!
    integer, intent(in)                                       :: scloop
!
    if (rank == 0 .and. textoutput .eq. .true.) then
        write(*,*)
        write(*,*) 'SUBCASE ', scloop, ' von ', fesim%num_subcases
    end if
    
    if (aeroConverged .eq. .true.) then
      if (fesim%lasten%subcases(scloop)%readApameInput .ne. .true.) then
        return   ! entspricht 'cycle' im Umlaufprogramm main.f90
      else
        aeroConverged = .false.
        moved = .false.
        if (allocated(fesim%internals%aeroDispl) .eq. .true.) deallocate(fesim%internals%aeroDispl)
        oldMaxDispl = 0.d0
        maxDispl = 0.d0
      end if
    end if
    
    if (fesim%lasten%subcases(scloop)%upmats .eq. .true.) then
        if (rank == 0) then
            allocate ( fesim%internals%num_dof_vec(fesim%num_dof) )
            
            !
            ! =================================================================================================
            !
            ! Integrate homogeneous boundary conditions
            !
            ! only for the case of constant boundary conditions and varying loadvectors for subcases
             
            call SYSTEM_CLOCK(startC,ctrC,ctmC)
            call SYSTEM_CLOCK(startI2,ctrI2,ctmI2)

            if (textoutput .eq. .true.) write(*,*) 'START - Erstelle Indexvektor'
            call boundary_varload_control (fesim,scloop)
            if (textoutput .eq. .true.) write(*,*) 'ENDE  - Erstelle Indexvektor'

        end if
         
        call bcast_internals(fesim%internals,fesim%num_dof)
         
        if (scloop .GT. 1) then
            if (fesim%is_lsolid20 == .true.) then
                if (rank .eq. 0) then
                    if (textoutput .eq. .true.) write(*,*) 'START - Aktualisiere LSolid20 Eigenschaften'
                    call lsolid20_prepareProperties(fesim,scloop)
                    if (textoutput .eq. .true.) write(*,*) 'ENDE  - Aktualisiere LSolid20 Eigenschaften'
                end if
                 
                WRITE(*,*) 'Verteile Eigenschaften'
                 
                call bcast_eigenschaften(fesim%eigenschaften,.false.,.false.,.false.,.true.)
                 
            end if
        end if

        !________________________________________________________________________________________
        !
        ! Init Matrices
        !
        !----------------------------------------------------------------------------------------
        if (textoutput .eq. .true. .and. rank .eq. 0) write(*,*) 'START - Initialisiere Matrizen'
        call init_matrices(fesim, KaaS, KgaaS, rank)
        if (textoutput .eq. .true. .and. rank .eq. 0)  write(*,*) 'ENDE  - Initialisiere Matrizen'
        !
        !
        ! =================================================================================================
        !
        ! Compute loadvector(s) and stiffness matrices
        !
        if (rank == 0) then
            call SYSTEM_CLOCK(endeI2,ctrI2,ctmI2)
            initTime  = initTime + (endeI2-startI2)/DBLE(ctrI2)
            call SYSTEM_CLOCK(startS,ctrS,ctmS)
        end if

        if (textoutput .eq. .true. .and. rank .eq. 0) write(*,*) 'START - Erstelle Steifigkeitsmatrix'
        call sol1_get_stiffness_matrices(fesim, KaaS)
        if (textoutput .eq. .true. .and. rank .eq. 0) write(*,*) 'ENDE  - Erstelle Steifigkeitsmatrix'
         
        if (rank == 0) then
            call SYSTEM_CLOCK(endeS,ctrS,ctmS)
            stiffTime = stiffTime + (endeS-startS)/DBLE(ctrS)
        end if

        if (textoutput .eq. .true. .and. rank .eq. 0) write(*,*) 'Kaas assemble'
        call MatAssemblyBegin(KaaS,MAT_FINAL_ASSEMBLY,ierr); CHKERRA(ierr)
        call MatAssemblyEnd  (KaaS,MAT_FINAL_ASSEMBLY,ierr); CHKERRA(ierr)
        if (textoutput .eq. .true. .and. rank .eq. 0) write(*,*) 'Kaas assemble fertig'
         
    end if
    
    !
    ! =================================================================================================
    !
    ! Herausschreiben des Netzes als VTK
    !
    if (rank == 0 .and. fesim%ausgabe%outputVTK .eq. .true.) then
        if (textoutput .eq. .true.) write(*,*) 'START - Schreibe Modell als VTK'
        write(fname, '(A10,I4.4,A4)') 'output_sc_', scloop, '.vtk'

        call vtkoutprep(fname)                                                                                    ! VTK-Datei öffnen und vorbereiten
        call vtkoutmesh(fesim) 
        call vtkoutelementprop(fesim)                                                                             ! Netz herausschreiben
        if (fesim%ausgabe%outputElemCoord .eq. .true.) call vtkoutelemcoord(fesim)                                ! Herausschreiben der Elementkoordinatensysteme
        if (fesim%ausgabe%outputBoundCond .eq. .true.) call vtkoutfixeddispl(fesim)                               ! Herausschreiben der Randbedingungen (Kräfte, Momente, Nullverschiebungen, Nullrotationen)
         
        if (textoutput .eq. .true.) write(*,*) 'ENDE  - Schreibe Modell als VTK'
    end if

    if (rank == 0) call SYSTEM_CLOCK(startL,ctrL,ctmL)
    
    if (textoutput .eq. .true. .and. rank .eq. 0) write(*,*) 'START - Initialisiere Vektoren'
    call init_vectors (fesim%internals%dim_dof, FaS, UaS)
    if (textoutput .eq. .true. .and. rank .eq. 0) write(*,*) 'ENDE  - Initialisiere Vektoren'
         
    if (rank == 0) then

        call aero_coupling(textoutput, fesim, scloop, cdi)
    
        if (textoutput .eq. .true.) write(*,*) 'START - Erstelle Lastvektor'
        allocate(Fout(fesim%num_dof))
        call load_get_loadvector(fesim, FaS, Fout, scloop)
        if (fesim%ausgabe%outputUser .ne. .true. .and. fesim%calculateReactForce .ne. .true.) then
          deallocate(Fout)
        end if
        if (textoutput .eq. .true.) write(*,*) 'ENDE  - Erstelle Lastvektor'

        if (textoutput .eq. .true.) then
            write(*,'(A41,I10)')   'Elementanzahl                           :', fesim%num_elements
            write(*,'(A41,I10)')   'Knotenanzahl                            :', size(fesim%knoten%nodes,1)
            write(*,'(A41,I10)')   'maximale Anzahl an Freiheiten           :', fesim%num_dof
            write(*,'(A41,I10)')   'Anzahl an Freiheiten                    :', fesim%internals%dim_dof
        endif
    
    end if

    if (rank == 0 .and. ((fesim%is_aeroload2d .eq. .true.) .or. (fesim%is_aeroload3d .eq. .true.))) then
        if (textoutput .eq. .true.) write(*,*) 'START - Entferne Aero-Lasten'
        call load_remove_aeropressure(fesim)
        if (textoutput .eq. .true.) write(*,*) 'ENDE  - Entferne Aero-Lasten'
    end if
    
    if (rank == 0) then
        call SYSTEM_CLOCK(endeL,ctrL,ctmL) 
        lvecTime  = lvecTime + (endeL-startL)/DBLE(ctrL)
    end if

    if (textoutput .eq. .true. .and. rank .eq. 0) write(*,*) 'FaS assemble'
    call VecAssemblyBegin(FaS,ierr); CHKERRA(ierr);
    call VecAssemblyEnd(FaS,ierr); CHKERRA(ierr);
    if (textoutput .eq. .true. .and. rank .eq. 0) write(*,*) 'FaS assemble fertig'

    ! =================================================================================================
    !
    ! Static solution
    !
    !  call write_matrix(KaaS,'Kaas_mpc')
    if (rank == 0) call SYSTEM_CLOCK(start,ctr,ctm)
    call sol1_solve(fesim,ksp,KaaS,UaS,FaS,rank,scloop)
    if (rank == 0) then
        call SYSTEM_CLOCK(ende,ctr,ctm)
        solvTime  = solvTime + (ende-start)/DBLE(ctr)
    end if
    !
    ! Calculation of result values requiring stiffness matrix
    if (fesim%calculateTSE .EQV. .TRUE.) then
        call sol1_calculate_tse(fesim, KaaS, UaS, rank)
    end if

    if (fesim%sol == 1) then
        if (scloop .LT. fesim%num_subcases) then
            if (fesim%lasten%subcases(scloop+1)%upmats .eq. .true.) then
                call MatDestroy(KaaS,ierr); CHKERRA(ierr)
                call KSPDestroy(ksp,ierr); CHKERRA(ierr)
            end if
        else
            call MatDestroy(KaaS,ierr); CHKERRA(ierr)
            call KSPDestroy(ksp,ierr); CHKERRA(ierr)
        end if
#ifndef mod_slepc
    else
        call KSPDestroy(ksp,ierr); CHKERRA(ierr)
#endif
    end if
    
    if (textoutput .eq. .true. .and. rank .eq. 0) write(*,*) 'statische Loesung fertig'

    if (rank == 0) then
        call SYSTEM_CLOCK(endeC,ctrC,ctmC)
        compTime  = compTime + (endeC-startC)/DBLE(ctrC)
        startC = endeC
    end if

    if (textoutput .eq. .true. .and. rank .eq. 0) write(*,*) '  Loesungsvektor sammeln'

    call VecCreate(PETSC_COMM_WORLD,UaSseq,ierr); CHKERRA(ierr)
    call VecSetSizes(UaSseq,PETSC_DECIDE,fesim%internals%dim_dof,ierr); CHKERRA(ierr)
    call VecSetFromOptions(UaSseq,ierr); CHKERRA(ierr)
    
    call VecScatterCreateToZero(UaS, scatter, UaSseq, ierr); CHKERRA(ierr)
    
    call VecScatterBegin(scatter, UaS, UaSseq, INSERT_VALUES, SCATTER_FORWARD, ierr); CHKERRA(ierr)
    call VecScatterEnd  (scatter, UaS, UaSseq, INSERT_VALUES, SCATTER_FORWARD, ierr); CHKERRA(ierr)
    
    call VecScatterDestroy(scatter, ierr); CHKERRA(ierr)
    
    if (textoutput .eq. .true. .and. rank .eq. 0) write(*,*) '  Loesungsvektor sammeln fertig'

    if (rank == 0) then
    
        ! Neusortieren der Loesung
        allocate(Utot(fesim%num_dof))
    
        call sol1_get_Utot(fesim, Utot, UaSseq)
         
        if (textoutput .eq. .true.) write(*,*) 'Verschiebungsvektor fertig'         
 
        if ((fesim%is_aeroload2d .eq. .true.) .or. (fesim%is_aeroload3d .eq. .true.)) then
            call printdisplacements(fesim,scloop,Utot)
            oldMaxDispl = maxDispl
            maxDispl = 0.d0
            do ii =  1, fesim%num_nodes
                displ = SQRT(Utot((ii-1)*ndof+1)**2 + Utot((ii-1)*ndof+2)**2 + Utot((ii-1)*ndof+3)**2)
                if (displ .gt. maxDispl) then
                    maxDispl = displ
                end if
            end do
            write(*,'(A34,F8.3,A2,F8.3)') '  maximale Verschiebungsänderung: ', ABS((maxDispl-oldMaxDispl)/maxDispl)*100.d0, ' %', maxDispl-oldMaxDispl
            if (ABS((maxDispl-oldMaxDispl)/maxDispl) .lt. fesim%aeroDisplEps) then
              aeroConverged = .true.
            end if
        end if
        !
        ! =================================================================================================
        !
        ! Output static solution
         
        call sol1_output_control(fesim,Utot,scloop,cdi,aeroConverged)
        !
        if (fesim%ausgabe%outputUser .eq. .true.) then
            call sol1_output_user(fesim,Utot,Fout,scloop,cdi,aeroConverged)
        end if
    end if
    
    if (fesim%calculateReactForce .EQ. .TRUE.) then
        call sol1_reaction_forces(fesim, Utot, Fout, rank)
    end if
         
    if (rank == 0) then
        if (fesim%ausgabe%outputUser .eq. .true. .or. fesim%calculateReactForce .eq. .true.) then
            deallocate(Fout)
        end if
         
        if (scloop .LT. fesim%num_subcases) then
            if (fesim%lasten%subcases(scloop+1)%upgeom .eq. .true.) then
                call update_geometry(fesim,Utot)
            end if
        end if
         
        if (fesim%sol == 1) deallocate(Utot)
    end if
! 
    call VecDestroy(UaSseq,ierr); CHKERRA(ierr)
    call VecDestroy(FaS,ierr); CHKERRA(ierr)
    call VecDestroy(UaS,ierr); CHKERRA(ierr)
    
    if (scloop .LT. fesim%num_subcases) then
        if (fesim%lasten%subcases(scloop+1)%upgeom .eq. .true.) then
            call bcast_knoten(fesim%knoten,.true.)
        end if
    end if

    !
    ! =================================================================================================
    !
    ! Buckling Calculations
    !

    if (fesim%sol == 2) then
    
      if (fesim%is_multistep .eq. .true.) then
          write(*,*) 'Im Falle einer Multisteprechnung sind Stabilitätsbetrachtungen nicht möglich'
          stop
      end if
      
      call MPI_Bcast (fesim%internals%failed, 1, MPI_LOGICAL, 0, PETSC_COMM_WORLD, ierr); CHKERRA(ierr)

      if ((fesim%skipFailed .EQ. .TRUE. .AND. fesim%internals%failed .EQ. .TRUE.) .OR. fesim%lasten%subcases(scloop)%skipBuckling .EQ. .TRUE.) then
        if (fesim%ausgabe%outputShort .eq. .true.) then
            if (rank == 0) then
                OPEN(shortoutFile, file = 'output_FEM.txt', STATUS='OLD', POSITION='APPEND')
                write(shortoutFile,'(A8,I3,A27,E25.18)') 'Subcase ', scloop, ', Eigenwert             :  ', 1.d300
                CLOSE(shortoutFile)
            end if
        end if
      else
        !
        ! =================================================================================================
        !
        ! Compute element stiffness matrices and assemble it to total stiffness matrix
        !
        !call MatZeroEntries(KgaaS,ierr); CHKERRA(ierr)
        if (fesim%lasten%subcases(scloop)%upmats .ne. .true.) then
            if (textoutput .eq. .true. .and. rank .eq. 0) WRITE(*,*) 'Reset KGeo'
            if (scloop .gt. 1) call sol2_reset_kgeo(fesim%internals%dim_dof, KgaaS)
            if (textoutput .eq. .true. .and. rank .eq. 0) WRITE(*,*) 'Reset KGeo fertig'
        end if

        ! Create KgaaS
        if (textoutput .eq. .true. .and. rank .eq. 0) WRITE(*,*) 'START - Erstelle geometrische Steifigkeitsmatrix'
        if (rank == 0) call sol2_get_geometric_stiffness_matrix(fesim, scloop, KgaaS, Utot(:))
        if (textoutput .eq. .true. .and. rank .eq. 0) WRITE(*,*) 'ENDE  - Erstelle geometrische Steifigkeitsmatrix'
         
        if (textoutput .eq. .true. .and. rank .eq. 0) WRITE(*,*) 'START - MPI Assemblieren der geometrischen Steifigkeitsmatrix'
        call MatAssemblyBegin(KgaaS,MAT_FINAL_ASSEMBLY,ierr); CHKERRA(ierr)
        call MatAssemblyEnd  (KgaaS,MAT_FINAL_ASSEMBLY,ierr); CHKERRA(ierr)
        if (textoutput .eq. .true. .and. rank .eq. 0) WRITE(*,*) 'ENDE  - MPI Assemblieren der geometrischen Steifigkeitsmatrix'
         
        !call sol2_solve_condensed_problem (rank,fesim,ksp,KaaS,KgaaS,scloop)
        !STOP
         
        if (rank == 0) call SYSTEM_CLOCK(start,ctr,ctm)
        call sol2_solve(eps, ksp, KaaS, KgaaS, fesim%numEigVal)
        if (rank == 0) then
            call SYSTEM_CLOCK(ende,ctr,ctm)
            solvEigTime  = solvEigTime + (ende-start)/DBLE(ctr)
        end if
        if (textoutput .eq. .true. .and. rank .eq. 0) write(*,*) 'Eigenwertproblem geloest'
        call EPSGetConverged(eps,nconv,ierr); CHKERRA(ierr)

        if (rank == 0) then  
            allocate(eigvals(nconv))
            allocate(eigenvectors(fesim%internals%dim_dof,nconv))
        else
            allocate(eigvals(1))
            allocate(eigenvectors(1,1))
        end if
         
        if (textoutput .eq. .true. .and. rank .eq. 0) WRITE(*,*) 'START - Extrahieren der Eigenwerte'
        call sol2_extract_eigenvalues(eps, eigvals, eigenvectors, nconv, fesim%internals%dim_dof, rank)
        if (textoutput .eq. .true. .and. rank .eq. 0) WRITE(*,*) 'ENDE  - Extrahieren der Eigenwerte'

        if (rank == 0) call sol2_output_control(fesim, scloop, eigvals, eigenvectors, nconv, aeroConverged)

        deallocate(eigvals)
        deallocate(eigenvectors)
        call EPSDestroy(eps,ierr); CHKERRA(ierr)
        call MatDestroy(KgaaS,ierr); CHKERRA(ierr)
      
      end if
         
      if (rank == 0) deallocate(Utot)
  
      if (scloop .LT. fesim%num_subcases) then
        if (fesim%lasten%subcases(scloop+1)%upmats .eq. .true.) then
          call MatDestroy(KaaS,ierr); CHKERRA(ierr)
       
          if (allocated(d_nzz_kgeo) .eq. .true.) deallocate(d_nzz_kgeo)
          if (allocated(o_nzz_kgeo) .eq. .true.) deallocate(o_nzz_kgeo)

#ifdef mod_slepc                
          call KSPDestroy(ksp,ierr); CHKERRA(ierr)
#endif
        end if
      else
        call MatDestroy(KaaS,ierr); CHKERRA(ierr)
         
        if (allocated(d_nzz_kgeo) .eq. .true.) deallocate(d_nzz_kgeo)
        if (allocated(o_nzz_kgeo) .eq. .true.) deallocate(o_nzz_kgeo)

#ifdef mod_slepc            
        call KSPDestroy(ksp,ierr); CHKERRA(ierr)
#endif
      end if
!         
    end if
#ifdef video      
    if (rank == 0) call vtkoutend()
    if (aeroConverged .eq. .true. .and. moved .eq. .false.) then
!       moved = .true.
      call system('mv ' // TRIM(fname) // ' final/.')
    end if
#endif

    if (scloop .LT. fesim%num_subcases) then
        if (fesim%lasten%subcases(scloop+1)%upmats .eq. .true.) then
            deallocate(fesim%internals%num_dof_vec)
        end if
    else
        deallocate(fesim%internals%num_dof_vec)
    end if

  end subroutine fipps_subcase

  !
  ! =================================================================================================
  !
  ! Abschluss: Rechenzeitausgabe, Rueckraeumung, Finalisierung von SLEPc
  !
  subroutine fipps_finalize()
!
    if (rank == 0) then
      call writeCompTime(compTime,solvTime,solvEigTime,initTime,stiffTime,lvecTime)
    end if
       
    
    if (allocated(d_nzz_kgeo) .eq. .true.) deallocate(d_nzz_kgeo)
    if (allocated(o_nzz_kgeo) .eq. .true.) deallocate(o_nzz_kgeo)

    call free_globals()
    call free_mem_fesim(fesim)
  
    call SlepcFinalize(ierr); CHKERRA(ierr)

  end subroutine fipps_finalize

! call SYSTEM_CLOCK(ende,ctr,ctm)
! 
! WRITE(*,*) 'Verwendete Zeit (s): ', (ende-start)/dble(ctr)
END PROGRAM FiPPS

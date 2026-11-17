! *****************************COPYRIGHT*******************************
!
! Copyright (c) 2008, Regents of the University of California
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are
! met:
!
!     * Redistributions of source code must retain the above copyright
!       notice, this list of conditions and the following disclaimer.
!     * Redistributions in binary form must reproduce the above
!       copyright notice, this list of conditions and the following
!       disclaimer in the documentation and/or other materials provided
!       with the distribution.
!     * Neither the name of the University of California, Irvine nor the
!       names of its contributors may be used to endorse or promote
!       products derived from this software without specific prior
!       written permission.
!
!       THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
!       IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!       TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
!       PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
!       OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
!       EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
!       PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
!       PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
!       LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
!       NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
!       SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
! *****************************COPYRIGHT*******************************
!
!  Description:
!    Driver for fully implicit ODE integrator.
!    Part of the ASAD chemical solver.
!
!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
!   Called from ASAD_CDRIVE
!
!
!     MJPDRIV  - Driver for MJP fully implicit ODE integrator.
!
!     Michael Prather            Earth System Science
!     Oliver Wild                University of California, Irvine
!
!     ASAD: mjpdriv              Version: mjpdriv.f 1.0 04/17/97
!
!     Purpose.
!     --------
!     To organise the integration of the chemical rate equations using
!     the MJP implicit integrator.
!
!     Interface
!     ---------
!     Called from chemistry driver routine *cdrive*.
!
!     This routine assumes that all the bi-,tri-, phot- and het-
!     reaction rates have been computed prior to this routine.
!     It is also assumed that the species array, y, has been set
!     from the array, f, passed to ASAD, and that constant species
!     have been set. This can be done by calling the routine fyinit.
!
!     Method.
!     -------
!     This routine calls the MJP integrator once for each gridpoint
!     of the one-dimensional arrays passed to it.
!     If convergence isn't achieved, the time step is halved, and the
!     integrator called again - this is continued until either
!     convergence is achieved or the minimum time step length is
!     encountered (currently 1.E-05 seconds).
!
!     Local variables
!     ---------------
!     ncsteps_initial -  Stores number of basic chemical steps 'ncsteps'
!     cdt_initial     -  Stores basic chemical time step length 'ctd'
!     f_initial       -  Stores family concentrations at beginning of call
!     exit_code       -  Error code from the integrator:
!                        0 = successful return
!                        1 = negatives encountered
!                        2 = convergence failure after 'nrsteps' iterations
!                        3 = convergence failure due to divergence - 'NaN's
!                        4 = convergence failure (as '2') but set debugging
!
!  Code Description:
!    Language:  FORTRAN 90
!
! ######################################################################
!
MODULE asad_spmjpdriv_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'ASAD_SPMJPDRIV_MOD'

CONTAINS

SUBROUTINE asad_spmjpdriv(ix,jy,nlev,n_points)

USE asad_mod, ONLY: cdt, f, jpcspf, jpspec, ltrig,                             &
                    ncsteps, ncsteps_factor, nitfg, speci, specf, y
USE ukca_config_specification_mod, ONLY: ukca_config
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
USE ereport_mod, ONLY: ereport
USE umPrintMgr, ONLY: umMessage, umPrint, PrintStatus, PrStatus_Oper

USE ukca_um_legacy_mod,  ONLY: mype

USE errormessagelength_mod, ONLY: errormessagelength

USE asad_diffun_mod, ONLY: asad_diffun
USE asad_spimpmjp_mod, ONLY: asad_spimpmjp
USE asad_ftoy_mod, ONLY: asad_ftoy
USE asad_range_mod, ONLY: asad_range_print, asad_range_reset,                  &
                          asad_range_track, asad_range_print_tracked
USE asad_range_mod, ONLY: asad_range_print_conv, asad_range_accumulate,        &
                          asad_range_accumulate_rates

IMPLICIT NONE

! Subroutine interface
INTEGER, INTENT(IN) :: n_points     ! Number of points in chunk
INTEGER, INTENT(IN) :: ix           ! i counter
INTEGER, INTENT(IN) :: jy           ! j counter
INTEGER, INTENT(IN) :: nlev         ! Number of model levels

! Local variables
INTEGER, PARAMETER :: max_redo=128  ! Max times for halving TS, was 16
INTEGER :: exit_code                ! Convergence exit code
INTEGER :: ncsteps_initial          ! Initial number of chemistry steps
INTEGER :: iredo                    ! Number of iterations for convergence
INTEGER :: js                       ! Species counter
INTEGER :: iter                     ! Timestep iteration for current window
INTEGER :: location                 ! Array index for logging
INTEGER :: solver_iter              ! Nonlinear solver iteration count
INTEGER :: num_iter                 ! Total nonlinear solver iteration count
INTEGER :: jit                      ! Iteration for asad_ftoy
INTEGER :: n_calls                  ! Number of calls of asad_spimpmjp
INTEGER :: n_fail                   ! Number of failed calls
INTEGER :: n_fail2                  ! Failed calls with exit code 2
INTEGER :: n_fail3                  ! Failed calls with exit code 3
LOGICAL :: l_noconv                 ! .TRUE. if solver gave up
LOGICAL :: l_event                  ! .TRUE. if this call gets the per-call lines

INTEGER :: errcode                  ! Variable passed to ereport
LOGICAL :: not_first_call = .FALSE.

CHARACTER(LEN=errormessagelength) :: cmessage

REAL :: cdt_initial                 ! Initial chemistry timestep
REAL :: f_initial(n_points,jpcspf)  ! Saved f array from previous solver call

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ASAD_SPMJPDRIV'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Stash initial values
ncsteps_initial = ncsteps
cdt_initial = cdt
iredo = ncsteps_factor
f_initial(1:n_points,:)=f(1:n_points,:)

ltrig=.FALSE.

IF (ukca_config%l_ukca_asad_full) THEN
  location = 1
ELSE IF (ukca_config%l_ukca_asad_columns) THEN
  ! mapping to theta_field
  location = ix + ((jy - 1)*ukca_config%row_length)
ELSE
  location = nlev
END IF

! ASADTRACK: start of tracking for this solver call (y at entry)
CALL asad_range_reset()
CALL asad_range_track(y(1:n_points,:))

! Start iterations here.
num_iter = 0
n_calls = 0
n_fail = 0
n_fail2 = 0
n_fail3 = 0
l_noconv = .FALSE.
iter = 1
DO WHILE (iter <= iredo)
  CALL asad_spimpmjp(exit_code, ix, jy, nlev, n_points, location, solver_iter)
  n_calls = n_calls + 1
  num_iter = num_iter + solver_iter

  IF (exit_code == 0) THEN
    ! Solver convergence
    iter = iter + 1
  ELSE
    n_fail = n_fail + 1
    IF (exit_code == 2) n_fail2 = n_fail2 + 1
    IF (exit_code == 3) n_fail3 = n_fail3 + 1
    f(1:n_points,:)=f_initial(1:n_points,:)

    IF (exit_code == 4) THEN
      ! Debug slow convergence systems - switch this on in 'spimpmjp'
      IF (ltrig) THEN
        errcode=1
        cmessage='Slow-converging system, Set printstatus for Jacobian debug'
        DO js=1,jpspec
          WRITE(umMessage,'(a4,i6,a12,2e14.5,i12)') 'y: ',js,speci(js),        &
              MAXVAL(y(:,js)), MINVAL(y(:,js)),SIZE(y(:,js))
          CALL umPrint(umMessage,src='asad_spmjpdriv')
        END DO
        CALL ereport('ASAD_SPMJPDRIV',errcode,cmessage)
      END IF

      ltrig=.TRUE.
    ELSE

      ! Reset for failed convergence
      ncsteps = ncsteps*2
      cdt = cdt/2.0
      iredo = iredo*2

      IF (ukca_config%l_ukca_debug_asad) THEN
        ! Added extra print statements here for verbosity
        WRITE(umMessage,"('ASAD: failed to converge at location  = ',I0)")     &
              location
        CALL umPrint(umMessage,src='asad_spmjpdriv')
        WRITE(umMessage,'(A,I0,A,I0,A,E18.8)')                                 &
              'ASAD: halving timestep: ncsteps = ', ncsteps,                   &
              ' iredo = ', iredo, ' cdt = ', cdt
        CALL umPrint(umMessage,src='asad_spmjpdriv')
      END IF

      IF (cdt < 1.0e-05) THEN
        errcode=2
        cmessage=' Time step now too short'
        CALL ereport('ASAD_SPMJPDRIV',errcode,cmessage)
      END IF

      ! Drop out if too many successive halvings fail
      IF (iredo >= max_redo) THEN
        IF (printstatus >= prstatus_oper) THEN
          WRITE(umMessage,"(' Resetting array after',i4,' iterations')") iredo
          CALL umPrint(umMessage,src='asad_spmjpdriv')
          WRITE(umMessage,"('NO CONVERGENCE location: ',i4,' pe: ',i4)")       &
              location,mype
          CALL umPrint(umMessage,src='asad_spmjpdriv')
        END IF
        l_noconv = .TRUE.
        EXIT
      END IF

    END IF

    ! Call asad_ftoy with jit = 0 to reinitialise y array
    jit = 0
    CALL asad_ftoy( not_first_call, nitfg, jit, n_points, ix, jy, nlev )
    CALL asad_diffun( n_points )
    CALL asad_range_track(y(1:n_points,:))
    iter = 1

  END IF
END DO

IF (iredo > 2) THEN
  WRITE(umMessage,"('   No. iterations =',i2)") iredo
  CALL umPrint(umMessage,src='asad_spmjpdriv')
END IF

IF (ukca_config%l_ukca_debug_asad) THEN
  ! Select which print statement to write depending on how asad solver is called
  IF (ukca_config%l_ukca_asad_full) THEN
    WRITE(umMessage, "('Iterations in spmjpdriv = ',I0)") num_iter
  ELSE IF (ukca_config%l_ukca_asad_columns) THEN
    WRITE(umMessage,                                                           &
    "('Iterations in spmjpdriv = ',I0,' ix = ',I0,' jy = ',I0)")               &
    num_iter, ix, jy
  ELSE
    WRITE(umMessage,"('Iterations in spmjpdriv = ',I0,' k = ',I0)")            &
    num_iter, nlev
  END IF
  CALL umPrint(umMessage,src='asad_spmjpdriv')
END IF

! ASADRANGE: min and max of y (state at end of solver call), one line per
! species. Draft for the GC_Faster range study.
CALL asad_range_accumulate(y(1:n_points,:), speci, nlev, num_iter, n_fail,     &
                           n_fail2, n_fail3, l_noconv, l_event)
CALL asad_range_accumulate_rates(speci, specf, nlev, n_points)
IF (l_event) THEN
  CALL asad_range_print('y', y(1:n_points,:), speci, ix, jy, nlev)
  CALL asad_range_print_tracked('y', speci, ix, jy, nlev)
  CALL asad_range_print_conv(ix, jy, nlev, n_points, ncsteps_factor, iredo,    &
                             n_calls, n_fail, num_iter, exit_code, l_noconv,   &
                             n_fail2, n_fail3)
END IF
! end of added lines for GC_Faster project

ncsteps = ncsteps_initial
cdt = cdt_initial

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE asad_spmjpdriv

END MODULE asad_spmjpdriv_mod

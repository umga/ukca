! Draft for the GC_Faster project's range study.
!
! This code prints the smallest and largest value of an array.
! One line for each column (a species or a reaction).
! Each line has: the word ASADRANGE (to grep), the time step number
! and the date.
!
MODULE asad_range_mod

USE umPrintMgr,         ONLY: umMessage, umPrint
USE ukca_um_legacy_mod, ONLY: mype
USE, INTRINSIC :: iso_fortran_env, ONLY: real32

IMPLICIT NONE
PRIVATE

PUBLIC :: asad_range_set_time, asad_range_print
PUBLIC :: asad_range_reset, asad_range_track, asad_range_print_tracked
PUBLIC :: asad_range_print_conv

INTERFACE asad_range_print
  MODULE PROCEDURE asad_range_print_2d, asad_range_print_1d
END INTERFACE

! Top-level time information, set once per UKCA chemistry call (by the
! master thread, before the OpenMP region that calls the solver)
INTEGER, SAVE :: step_no = -1
INTEGER, SAVE :: ymdh(4) = -1

CHARACTER(LEN=*), PARAMETER :: fmt_line =                                      &
  '(A,5(1X,I0),1X,A,1X,I0,1X,A,2(1X,ES15.6E3),4(1X,I0),2(1X,ES15.6E3))'

! Running min and max of y during one solver call (one copy per thread)
INTEGER, PARAMETER :: maxsp = 400
REAL, PARAMETER    :: tiny32 = TINY(1.0_real32)
REAL, PARAMETER    :: huge32 = HUGE(1.0_real32)
REAL, PARAMETER    :: big    = HUGE(1.0)

REAL, SAVE    :: trk_min(maxsp)
REAL, SAVE    :: trk_max(maxsp)
INTEGER, SAVE :: trk_nzero(maxsp)
INTEGER, SAVE :: trk_nsmall(maxsp)
INTEGER, SAVE :: trk_nbig(maxsp)
INTEGER, SAVE :: trk_nneg(maxsp)
INTEGER, SAVE :: trk_nsnap
!$OMP THREADPRIVATE(trk_min, trk_max, trk_nzero, trk_nsmall, trk_nbig,         &
!$OMP               trk_nneg, trk_nsnap)

CHARACTER(LEN=*), PARAMETER :: fmt_trk =                                       &
  '(A,5(1X,I0),1X,A,1X,I0,1X,A,2(1X,ES15.6E3),4(1X,I0),2(1X,ES15.6E3),5(1X,I0))'

CONTAINS

!------------------------------------------------------------------------------
SUBROUTINE asad_range_set_time(timestep_number, year, month, day, hour)

IMPLICIT NONE

INTEGER, INTENT(IN) :: timestep_number
INTEGER, INTENT(IN) :: year
INTEGER, INTENT(IN) :: month
INTEGER, INTENT(IN) :: day
INTEGER, INTENT(IN) :: hour

step_no = timestep_number
ymdh(:) = [year, month, day, hour]

END SUBROUTINE asad_range_set_time

!------------------------------------------------------------------------------
! arr(n_points, ncol): one line per column. labels(ncol) gives a name per
! column (species name); if absent the label is '-' and the column number
! is the index (e.g. reaction number).
SUBROUTINE asad_range_print_2d(array_name, arr, labels, ix, jy, nlev)

IMPLICIT NONE

CHARACTER(LEN=*), INTENT(IN)           :: array_name
REAL,             INTENT(IN)           :: arr(:,:)
CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: labels(:)
INTEGER,          INTENT(IN)           :: ix
INTEGER,          INTENT(IN)           :: jy
INTEGER,          INTENT(IN)           :: nlev

INTEGER                :: j
CHARACTER(LEN=20)      :: lab

DO j = 1, SIZE(arr, 2)
  lab = '-'
  IF (PRESENT(labels)) lab = labels(j)
  WRITE(umMessage, fmt_line) 'ASADRANGE', step_no, ymdh,                       &
        TRIM(array_name), j, TRIM(lab),                                        &
        MINVAL(arr(:,j)), MAXVAL(arr(:,j)), mype, ix, jy, nlev,                &
        REAL(MINVAL(arr(:,j)), KIND=real32),                                   &
        REAL(MAXVAL(arr(:,j)), KIND=real32)
  CALL umPrint(umMessage, src='asad_range_mod')
END DO

END SUBROUTINE asad_range_print_2d

!------------------------------------------------------------------------------
! arr(n_points): one line (column number 1, label '-').
SUBROUTINE asad_range_print_1d(array_name, arr, ix, jy, nlev)

IMPLICIT NONE

CHARACTER(LEN=*), INTENT(IN) :: array_name
REAL,             INTENT(IN) :: arr(:)
INTEGER,          INTENT(IN) :: ix
INTEGER,          INTENT(IN) :: jy
INTEGER,          INTENT(IN) :: nlev

WRITE(umMessage, fmt_line) 'ASADRANGE', step_no, ymdh,                         &
      TRIM(array_name), 1, '-',                                                &
      MINVAL(arr), MAXVAL(arr), mype, ix, jy, nlev,                            &
      REAL(MINVAL(arr), KIND=real32), REAL(MAXVAL(arr), KIND=real32)
CALL umPrint(umMessage, src='asad_range_mod')

END SUBROUTINE asad_range_print_1d

!------------------------------------------------------------------------------
! Call at the start of every solver call: sets min, max and counters to start.
SUBROUTINE asad_range_reset()

IMPLICIT NONE

trk_min(:)    = big
trk_max(:)    = -big
trk_nzero(:)  = 0
trk_nsmall(:) = 0
trk_nbig(:)   = 0
trk_nneg(:)   = 0
trk_nsnap     = 0

END SUBROUTINE asad_range_reset

!------------------------------------------------------------------------------
! Call every time y has been recalculated: updates min, max and counters.
! The counters count values (points) over all snapshots of one solver call.
SUBROUTINE asad_range_track(arr)

IMPLICIT NONE

REAL, INTENT(IN) :: arr(:,:)

INTEGER :: j

DO j = 1, SIZE(arr, 2)
  trk_min(j) = MIN(trk_min(j), MINVAL(arr(:,j)))
  trk_max(j) = MAX(trk_max(j), MAXVAL(arr(:,j)))
  trk_nzero(j)  = trk_nzero(j)  + COUNT(arr(:,j) == 0.0)
  trk_nsmall(j) = trk_nsmall(j) +                                              &
                  COUNT(arr(:,j) > 0.0 .AND. arr(:,j) < tiny32)
  trk_nbig(j)   = trk_nbig(j)   + COUNT(arr(:,j) > huge32)
  trk_nneg(j)   = trk_nneg(j)   + COUNT(arr(:,j) < 0.0)
END DO
trk_nsnap = trk_nsnap + 1

END SUBROUTINE asad_range_track

!------------------------------------------------------------------------------
! Call at the end of the solver call: one ASADTRACK line per species.
SUBROUTINE asad_range_print_tracked(array_name, labels, ix, jy, nlev)

IMPLICIT NONE

CHARACTER(LEN=*), INTENT(IN) :: array_name
CHARACTER(LEN=*), INTENT(IN) :: labels(:)
INTEGER,          INTENT(IN) :: ix
INTEGER,          INTENT(IN) :: jy
INTEGER,          INTENT(IN) :: nlev

INTEGER :: j

DO j = 1, SIZE(labels)
  WRITE(umMessage, fmt_trk) 'ASADTRACK', step_no, ymdh,                        &
        TRIM(array_name), j, TRIM(labels(j)),                                  &
        trk_min(j), trk_max(j), mype, ix, jy, nlev,                            &
        REAL(trk_min(j), KIND=real32), REAL(trk_max(j), KIND=real32),          &
        trk_nsnap, trk_nzero(j), trk_nsmall(j), trk_nbig(j), trk_nneg(j)
  CALL umPrint(umMessage, src='asad_range_mod')
END DO

END SUBROUTINE asad_range_print_tracked

!------------------------------------------------------------------------------
! Call at the end of the solver call: one ASADCONV line per call.
! iredo0 = starting number of sub-steps, iredo = final number (it doubles at
! every halving of the time step), n_calls = calls of asad_spimpmjp,
! n_fail = failed calls, num_iter = iterations of all calls added up,
! exit_code = code of the last call (0 = converged), the last number is 1 if
! the solver gave up (no convergence after the maximum number of halvings).
SUBROUTINE asad_range_print_conv(ix, jy, nlev, n_points, iredo0, iredo,       &
                                 n_calls, n_fail, num_iter, exit_code,         &
                                 l_noconv, n_fail2, n_fail3)

IMPLICIT NONE

INTEGER, INTENT(IN) :: ix
INTEGER, INTENT(IN) :: jy
INTEGER, INTENT(IN) :: nlev
INTEGER, INTENT(IN) :: n_points
INTEGER, INTENT(IN) :: iredo0
INTEGER, INTENT(IN) :: iredo
INTEGER, INTENT(IN) :: n_calls
INTEGER, INTENT(IN) :: n_fail
INTEGER, INTENT(IN) :: num_iter
INTEGER, INTENT(IN) :: exit_code
LOGICAL, INTENT(IN) :: l_noconv
INTEGER, INTENT(IN) :: n_fail2
INTEGER, INTENT(IN) :: n_fail3

CHARACTER(LEN=*), PARAMETER :: fmt_conv = '(A,19(1X,I0))'

WRITE(umMessage, fmt_conv) 'ASADCONV', step_no, ymdh,                          &
      mype, ix, jy, nlev, n_points, iredo0, iredo, n_calls, n_fail,            &
      num_iter, exit_code, MERGE(1, 0, l_noconv), n_fail2, n_fail3
CALL umPrint(umMessage, src='asad_range_mod')

END SUBROUTINE asad_range_print_conv

END MODULE asad_range_mod

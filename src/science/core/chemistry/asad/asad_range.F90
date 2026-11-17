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
PUBLIC :: asad_range_accumulate, asad_range_flush

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
! Tables by level and species, kept over a period (shared by the threads:
! each model level is handled by one thread at a time)
INTEGER, PARAMETER :: maxlev = 150
INTEGER, PARAMETER :: period_steps = 72   ! 10 days = 720 at 72 steps per day; 72 for a one-day run
REAL, PARAMETER    :: event_big = 1.0e29
REAL, PARAMETER    :: tiny_r = TINY(1.0)

REAL, SAVE    :: acc_emin(maxlev,maxsp)
REAL, SAVE    :: acc_emax(maxlev,maxsp)
REAL, SAVE    :: acc_tmin(maxlev,maxsp)
REAL, SAVE    :: acc_tmax(maxlev,maxsp)
REAL, SAVE    :: acc_sumlog(maxlev,maxsp)
REAL, SAVE    :: acc_sumlog2(maxlev,maxsp)
INTEGER, SAVE :: acc_nzero(maxlev,maxsp)
INTEGER, SAVE :: acc_nsmall(maxlev,maxsp)
INTEGER, SAVE :: acc_nbig(maxlev,maxsp)
INTEGER, SAVE :: acc_nneg(maxlev,maxsp)
INTEGER, SAVE :: acc_npos(maxlev,maxsp)
INTEGER, SAVE :: lev_ncalls(maxlev)
INTEGER, SAVE :: lev_npts(maxlev)
INTEGER, SAVE :: lev_ntrack(maxlev)
INTEGER, SAVE :: lev_iter(maxlev)
INTEGER, SAVE :: lev_nfailcalls(maxlev)
INTEGER, SAVE :: lev_nf2(maxlev)
INTEGER, SAVE :: lev_nf3(maxlev)
INTEGER, SAVE :: lev_nnoconv(maxlev)
CHARACTER(LEN=20), SAVE :: acc_label(maxsp) = ' '
INTEGER, SAVE :: acc_nsp = 0
LOGICAL, SAVE :: acc_init = .FALSE.

CHARACTER(LEN=*), PARAMETER :: fmt_lev =                                       &
  '(A,5(1X,I0),1X,A,2(1X,I0),4(1X,ES15.6E3),5(1X,I0),2(1X,ES24.14E3))'
CHARACTER(LEN=*), PARAMETER :: fmt_lvl = '(A,15(1X,I0))'

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
IF (.NOT. acc_init) THEN
  CALL asad_range_acc_reset()
  acc_init = .TRUE.
END IF

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
!------------------------------------------------------------------------------
! Set the tables to start values (called once, and after each print).
SUBROUTINE asad_range_acc_reset()

IMPLICIT NONE

acc_emin(:,:)    = big
acc_emax(:,:)    = -big
acc_tmin(:,:)    = big
acc_tmax(:,:)    = -big
acc_sumlog(:,:)  = 0.0
acc_sumlog2(:,:) = 0.0
acc_nzero(:,:)   = 0
acc_nsmall(:,:)  = 0
acc_nbig(:,:)    = 0
acc_nneg(:,:)    = 0
acc_npos(:,:)    = 0
lev_ncalls(:)    = 0
lev_npts(:)      = 0
lev_ntrack(:)    = 0
lev_iter(:)      = 0
lev_nfailcalls(:) = 0
lev_nf2(:)       = 0
lev_nf3(:)       = 0
lev_nnoconv(:)   = 0

END SUBROUTINE asad_range_acc_reset

!------------------------------------------------------------------------------
! Call at the end of every solver call (after the last asad_range_track):
! adds this call to the tables of its model level. l_event is .TRUE. for a
! call that had a failed attempt, gave up, or reached event_big: only those
! calls get the old per-call lines.
SUBROUTINE asad_range_accumulate(arr, labels, nlev, num_iter, n_fail,          &
                                 n_fail2, n_fail3, l_noconv, l_event)

IMPLICIT NONE

REAL,             INTENT(IN)  :: arr(:,:)
CHARACTER(LEN=*), INTENT(IN)  :: labels(:)
INTEGER,          INTENT(IN)  :: nlev
INTEGER,          INTENT(IN)  :: num_iter
INTEGER,          INTENT(IN)  :: n_fail
INTEGER,          INTENT(IN)  :: n_fail2
INTEGER,          INTENT(IN)  :: n_fail3
LOGICAL,          INTENT(IN)  :: l_noconv
LOGICAL,          INTENT(OUT) :: l_event

INTEGER :: j
INTEGER :: nsp
INTEGER :: npts

nsp  = MIN(SIZE(arr, 2), maxsp)
npts = SIZE(arr, 1)
l_event = (n_fail > 0) .OR. l_noconv .OR. ANY(trk_max(1:nsp) > event_big)

IF (nlev < 1 .OR. nlev > maxlev) RETURN

IF (acc_nsp == 0) THEN
!$OMP CRITICAL (asad_range_labels)
  IF (acc_nsp == 0) THEN
    DO j = 1, nsp
      acc_label(j) = labels(j)
    END DO
    acc_nsp = nsp
  END IF
!$OMP END CRITICAL (asad_range_labels)
END IF

DO j = 1, nsp
  acc_emin(nlev,j) = MIN(acc_emin(nlev,j), MINVAL(arr(:,j)))
  acc_emax(nlev,j) = MAX(acc_emax(nlev,j), MAXVAL(arr(:,j)))
  acc_tmin(nlev,j) = MIN(acc_tmin(nlev,j), trk_min(j))
  acc_tmax(nlev,j) = MAX(acc_tmax(nlev,j), trk_max(j))
  acc_nzero(nlev,j)  = acc_nzero(nlev,j)  + trk_nzero(j)
  acc_nsmall(nlev,j) = acc_nsmall(nlev,j) + trk_nsmall(j)
  acc_nbig(nlev,j)   = acc_nbig(nlev,j)   + trk_nbig(j)
  acc_nneg(nlev,j)   = acc_nneg(nlev,j)   + trk_nneg(j)
  acc_npos(nlev,j)   = acc_npos(nlev,j)   + COUNT(arr(:,j) > 0.0)
  acc_sumlog(nlev,j) = acc_sumlog(nlev,j) +                                    &
        SUM(LOG10(MAX(arr(:,j), tiny_r)), MASK=(arr(:,j) > 0.0))
  acc_sumlog2(nlev,j) = acc_sumlog2(nlev,j) +                                  &
        SUM(LOG10(MAX(arr(:,j), tiny_r))**2, MASK=(arr(:,j) > 0.0))
END DO

lev_ncalls(nlev) = lev_ncalls(nlev) + 1
lev_npts(nlev)   = lev_npts(nlev) + npts
lev_ntrack(nlev) = lev_ntrack(nlev) + trk_nsnap * npts
lev_iter(nlev)   = lev_iter(nlev) + num_iter
IF (n_fail > 0) lev_nfailcalls(nlev) = lev_nfailcalls(nlev) + 1
lev_nf2(nlev)    = lev_nf2(nlev) + n_fail2
lev_nf3(nlev)    = lev_nf3(nlev) + n_fail3
IF (l_noconv) lev_nnoconv(nlev) = lev_nnoconv(nlev) + 1

END SUBROUTINE asad_range_accumulate

!------------------------------------------------------------------------------
! Call once per chemistry step after the solver (master thread, outside the
! OpenMP region): prints the tables when the period is over and resets them.
SUBROUTINE asad_range_flush()

IMPLICIT NONE

INTEGER :: lev
INTEGER :: j

IF (MOD(step_no, period_steps) /= 0) RETURN

DO lev = 1, maxlev
  IF (lev_ncalls(lev) == 0) CYCLE
  WRITE(umMessage, fmt_lvl) 'ASADLVL', step_no, ymdh, lev, mype,               &
        lev_ncalls(lev), lev_npts(lev), lev_ntrack(lev), lev_iter(lev),        &
        lev_nfailcalls(lev), lev_nf2(lev), lev_nf3(lev), lev_nnoconv(lev)
  CALL umPrint(umMessage, src='asad_range_mod')
  DO j = 1, acc_nsp
    WRITE(umMessage, fmt_lev) 'ASADLEV', step_no, ymdh, TRIM(acc_label(j)),    &
          lev, mype, acc_emin(lev,j), acc_emax(lev,j),                         &
          acc_tmin(lev,j), acc_tmax(lev,j),                                    &
          acc_nzero(lev,j), acc_nsmall(lev,j), acc_nbig(lev,j),                &
          acc_nneg(lev,j), acc_npos(lev,j),                                    &
          acc_sumlog(lev,j), acc_sumlog2(lev,j)
    CALL umPrint(umMessage, src='asad_range_mod')
  END DO
END DO

CALL asad_range_acc_reset()

END SUBROUTINE asad_range_flush

END MODULE asad_range_mod

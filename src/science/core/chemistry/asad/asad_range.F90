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
USE asad_mod, ONLY: jpbk, jptk, jppj, jphk, nbrkx, ntrkx, nprkx, nhrkx,        &
                    spb, spt, spj, sph

IMPLICIT NONE
PRIVATE

PUBLIC :: asad_range_set_time, asad_range_print
PUBLIC :: asad_range_reset, asad_range_track, asad_range_print_tracked
PUBLIC :: asad_range_print_conv
PUBLIC :: asad_range_accumulate, asad_range_flush
PUBLIC :: asad_range_track_rates, asad_range_accumulate_rates
PUBLIC :: asad_range_accumulate_rk

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

! prod (1), slos (2) and fdot (3) tracked through the solver call, one copy
! per OpenMP thread (like trk_*)
INTEGER, PARAMETER :: nrt = 3
REAL, SAVE    :: rt_min(maxsp,nrt)
REAL, SAVE    :: rt_max(maxsp,nrt)
INTEGER, SAVE :: rt_nzero(maxsp,nrt)
INTEGER, SAVE :: rt_nsmall(maxsp,nrt)
INTEGER, SAVE :: rt_nbig(maxsp,nrt)
INTEGER, SAVE :: rt_nneg(maxsp,nrt)
INTEGER, SAVE :: rt_n(nrt)
INTEGER, SAVE :: rt_nsnap
!$OMP THREADPRIVATE(rt_min, rt_max, rt_nzero, rt_nsmall, rt_nbig, rt_nneg,     &
!$OMP               rt_n, rt_nsnap)

! Tables by level for the three arrays, kept over a period (like acc_*)
REAL, SAVE    :: acr_tmin(maxlev,maxsp,nrt)
REAL, SAVE    :: acr_tmax(maxlev,maxsp,nrt)
INTEGER, SAVE :: acr_nzero(maxlev,maxsp,nrt)
INTEGER, SAVE :: acr_nsmall(maxlev,maxsp,nrt)
INTEGER, SAVE :: acr_nbig(maxlev,maxsp,nrt)
INTEGER, SAVE :: acr_nneg(maxlev,maxsp,nrt)
INTEGER, SAVE :: lev_nrate(maxlev)
CHARACTER(LEN=20), SAVE :: acr_label(maxsp,nrt) = ' '
INTEGER, SAVE :: acr_nsp(nrt) = 0
CHARACTER(LEN=7), PARAMETER :: rt_tag(nrt) =                                   &
  [CHARACTER(LEN=7) :: 'ASADPRD', 'ASADSLS', 'ASADFDT']
CHARACTER(LEN=*), PARAMETER :: fmt_rt =                                        &
  '(A,5(1X,I0),1X,A,2(1X,I0),2(1X,ES15.6E3),5(1X,I0))'

! rk (rate coefficients) by level and reaction, kept over a period
INTEGER, PARAMETER :: maxrx = 1500
REAL, SAVE    :: rk_tmin(maxlev,maxrx)
REAL, SAVE    :: rk_tmax(maxlev,maxrx)
INTEGER, SAVE :: rk_nzero(maxlev,maxrx)
INTEGER, SAVE :: rk_nsmall(maxlev,maxrx)
INTEGER, SAVE :: rk_nbig(maxlev,maxrx)
INTEGER, SAVE :: rk_nneg(maxlev,maxrx)
INTEGER, SAVE :: rk_nval(maxlev)
CHARACTER(LEN=64), SAVE :: rk_label(maxrx) = ' '
INTEGER, SAVE :: rk_nrx = 0

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
rt_min(:,:)    = big
rt_max(:,:)    = -big
rt_nzero(:,:)  = 0
rt_nsmall(:,:) = 0
rt_nbig(:,:)   = 0
rt_nneg(:,:)   = 0
rt_n(:)        = 0
rt_nsnap       = 0

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
acr_tmin(:,:,:)   = big
acr_tmax(:,:,:)   = -big
acr_nzero(:,:,:)  = 0
acr_nsmall(:,:,:) = 0
acr_nbig(:,:,:)   = 0
acr_nneg(:,:,:)   = 0
lev_nrate(:)      = 0
rk_tmin(:,:)      = big
rk_tmax(:,:)      = -big
rk_nzero(:,:)     = 0
rk_nsmall(:,:)    = 0
rk_nbig(:,:)      = 0
rk_nneg(:,:)      = 0
rk_nval(:)        = 0

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
INTEGER :: id

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
  DO id = 1, nrt
    DO j = 1, acr_nsp(id)
      WRITE(umMessage, fmt_rt) rt_tag(id), step_no, ymdh,                      &
            TRIM(acr_label(j,id)), lev, mype,                                  &
            acr_tmin(lev,j,id), acr_tmax(lev,j,id),                            &
            acr_nzero(lev,j,id), acr_nsmall(lev,j,id), acr_nbig(lev,j,id),     &
            acr_nneg(lev,j,id), lev_nrate(lev)
      CALL umPrint(umMessage, src='asad_range_mod')
    END DO
  END DO
  DO j = 1, rk_nrx
    IF (rk_label(j) == ' ') CYCLE
    WRITE(umMessage, fmt_rt) 'ASADRK ', step_no, ymdh,                         &
          TRIM(rk_label(j)), lev, mype,                                        &
          rk_tmin(lev,j), rk_tmax(lev,j),                                      &
          rk_nzero(lev,j), rk_nsmall(lev,j), rk_nbig(lev,j),                   &
          rk_nneg(lev,j), rk_nval(lev)
    CALL umPrint(umMessage, src='asad_range_mod')
  END DO
END DO

CALL asad_range_acc_reset()

END SUBROUTINE asad_range_flush

!------------------------------------------------------------------------------
! Track one rate array (id = 1 prod, 2 slos, 3 fdot) for one iteration.
! Counts use the absolute value; min and max keep the sign.
SUBROUTINE asad_range_track_one(id, arr)

IMPLICIT NONE

INTEGER, INTENT(IN) :: id
REAL,    INTENT(IN) :: arr(:,:)

INTEGER :: j
INTEGER :: nsp

nsp = MIN(SIZE(arr, 2), maxsp)
rt_n(id) = nsp
DO j = 1, nsp
  rt_min(j,id) = MIN(rt_min(j,id), MINVAL(arr(:,j)))
  rt_max(j,id) = MAX(rt_max(j,id), MAXVAL(arr(:,j)))
  rt_nzero(j,id)  = rt_nzero(j,id)  + COUNT(arr(:,j) == 0.0)
  rt_nsmall(j,id) = rt_nsmall(j,id) +                                          &
        COUNT(ABS(arr(:,j)) > 0.0 .AND. ABS(arr(:,j)) < tiny32)
  rt_nbig(j,id)   = rt_nbig(j,id)   + COUNT(ABS(arr(:,j)) > huge32)
  rt_nneg(j,id)   = rt_nneg(j,id)   + COUNT(arr(:,j) < 0.0)
END DO

END SUBROUTINE asad_range_track_one

!------------------------------------------------------------------------------
! Called at the end of asad_diffun: prod, slos and fdot of this iteration.
SUBROUTINE asad_range_track_rates(p_arr, l_arr, d_arr)

IMPLICIT NONE

REAL, INTENT(IN) :: p_arr(:,:)
REAL, INTENT(IN) :: l_arr(:,:)
REAL, INTENT(IN) :: d_arr(:,:)

CALL asad_range_track_one(1, p_arr)
CALL asad_range_track_one(2, l_arr)
CALL asad_range_track_one(3, d_arr)
rt_nsnap = rt_nsnap + 1

END SUBROUTINE asad_range_track_rates

!------------------------------------------------------------------------------
! Called at the end of every solver call (after asad_range_accumulate): adds
! the tracked prod, slos and fdot of this call to the tables of its level.
SUBROUTINE asad_range_accumulate_rates(labels_s, labels_f, nlev, npts)

IMPLICIT NONE

CHARACTER(LEN=*), INTENT(IN) :: labels_s(:)
CHARACTER(LEN=*), INTENT(IN) :: labels_f(:)
INTEGER,          INTENT(IN) :: nlev
INTEGER,          INTENT(IN) :: npts

INTEGER :: id
INTEGER :: j
INTEGER :: nn

IF (nlev < 1 .OR. nlev > maxlev) RETURN
IF (rt_nsnap == 0) RETURN

DO id = 1, nrt
  IF (id < 3) THEN
    nn = MIN(rt_n(id), SIZE(labels_s))
  ELSE
    nn = MIN(rt_n(id), SIZE(labels_f))
  END IF
  IF (acr_nsp(id) == 0) THEN
!$OMP CRITICAL (asad_range_labels_rt)
    IF (acr_nsp(id) == 0) THEN
      DO j = 1, nn
        IF (id < 3) THEN
          acr_label(j,id) = labels_s(j)
        ELSE
          acr_label(j,id) = labels_f(j)
        END IF
      END DO
      acr_nsp(id) = nn
    END IF
!$OMP END CRITICAL (asad_range_labels_rt)
  END IF
  DO j = 1, nn
    acr_tmin(nlev,j,id) = MIN(acr_tmin(nlev,j,id), rt_min(j,id))
    acr_tmax(nlev,j,id) = MAX(acr_tmax(nlev,j,id), rt_max(j,id))
    acr_nzero(nlev,j,id)  = acr_nzero(nlev,j,id)  + rt_nzero(j,id)
    acr_nsmall(nlev,j,id) = acr_nsmall(nlev,j,id) + rt_nsmall(j,id)
    acr_nbig(nlev,j,id)   = acr_nbig(nlev,j,id)   + rt_nbig(j,id)
    acr_nneg(nlev,j,id)   = acr_nneg(nlev,j,id)   + rt_nneg(j,id)
  END DO
END DO

lev_nrate(nlev) = lev_nrate(nlev) + rt_nsnap * npts

END SUBROUTINE asad_range_accumulate_rates
!------------------------------------------------------------------------------
! Label of one reaction: <jr>:<type>:<r1>+<r2>><p1>+<p2>
SUBROUTINE asad_range_set_rk_label(jr, ktype, s1, s2, s3, s4)

IMPLICIT NONE

INTEGER,          INTENT(IN) :: jr
CHARACTER(LEN=1), INTENT(IN) :: ktype
CHARACTER(LEN=*), INTENT(IN) :: s1
CHARACTER(LEN=*), INTENT(IN) :: s2
CHARACTER(LEN=*), INTENT(IN) :: s3
CHARACTER(LEN=*), INTENT(IN) :: s4

CHARACTER(LEN=8) :: idx

IF (jr < 1 .OR. jr > maxrx) RETURN
WRITE(idx, '(I0)') jr
rk_label(jr) = TRIM(idx)//':'//ktype//':'//TRIM(s1)//'+'//TRIM(s2)//           &
               '>'//TRIM(s3)//'+'//TRIM(s4)

END SUBROUTINE asad_range_set_rk_label

!------------------------------------------------------------------------------
! Called once per solver call at the entry of asad_spmjpdriv: min, max and
! counts of the rate coefficients rk(:,reaction) of this call, added to the
! tables of its level. rk is computed before the integrator (asad_cdrive)
! and does not change during the iterations.
SUBROUTINE asad_range_accumulate_rk(rk_arr, nlev, npts)

IMPLICIT NONE

REAL,    INTENT(IN) :: rk_arr(:,:)
INTEGER, INTENT(IN) :: nlev
INTEGER, INTENT(IN) :: npts

INTEGER :: j
INTEGER :: nn

IF (nlev < 1 .OR. nlev > maxlev) RETURN
nn = MIN(SIZE(rk_arr, 2), maxrx)

IF (rk_nrx == 0) THEN
!$OMP CRITICAL (asad_range_labels_rk)
  IF (rk_nrx == 0) THEN
    rk_label(:) = ' '
    DO j = 1, jpbk
      CALL asad_range_set_rk_label(nbrkx(j), 'B', spb(j,1), spb(j,2),          &
                                   spb(j,3), spb(j,4))
    END DO
    DO j = 1, jptk
      CALL asad_range_set_rk_label(ntrkx(j), 'T', spt(j,1), spt(j,2),          &
                                   spt(j,3), spt(j,4))
    END DO
    DO j = 1, jppj
      CALL asad_range_set_rk_label(nprkx(j), 'J', spj(j,1), spj(j,2),          &
                                   spj(j,3), spj(j,4))
    END DO
    DO j = 1, jphk
      CALL asad_range_set_rk_label(nhrkx(j), 'H', sph(j,1), sph(j,2),          &
                                   sph(j,3), sph(j,4))
    END DO
    rk_nrx = nn
  END IF
!$OMP END CRITICAL (asad_range_labels_rk)
END IF

DO j = 1, nn
  rk_tmin(nlev,j) = MIN(rk_tmin(nlev,j), MINVAL(rk_arr(1:npts,j)))
  rk_tmax(nlev,j) = MAX(rk_tmax(nlev,j), MAXVAL(rk_arr(1:npts,j)))
  rk_nzero(nlev,j)  = rk_nzero(nlev,j)  + COUNT(rk_arr(1:npts,j) == 0.0)
  rk_nsmall(nlev,j) = rk_nsmall(nlev,j) +                                      &
        COUNT(ABS(rk_arr(1:npts,j)) > 0.0 .AND. ABS(rk_arr(1:npts,j)) < tiny32)
  rk_nbig(nlev,j)   = rk_nbig(nlev,j)   + COUNT(ABS(rk_arr(1:npts,j)) > huge32)
  rk_nneg(nlev,j)   = rk_nneg(nlev,j)   + COUNT(rk_arr(1:npts,j) < 0.0)
END DO
rk_nval(nlev) = rk_nval(nlev) + npts

END SUBROUTINE asad_range_accumulate_rk

END MODULE asad_range_mod

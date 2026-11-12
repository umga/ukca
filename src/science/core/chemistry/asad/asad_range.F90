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

IMPLICIT NONE
PRIVATE

PUBLIC :: asad_range_set_time, asad_range_print

INTERFACE asad_range_print
  MODULE PROCEDURE asad_range_print_2d, asad_range_print_1d
END INTERFACE

! Top-level time information, set once per UKCA chemistry call (by the
! master thread, before the OpenMP region that calls the solver)
INTEGER, SAVE :: step_no = -1
INTEGER, SAVE :: ymdh(4) = -1

CHARACTER(LEN=*), PARAMETER :: fmt_line =                                      &
  '(A,5(1X,I0),1X,A,1X,I0,1X,A,2(1X,ES15.6E3),4(1X,I0))'

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
        MINVAL(arr(:,j)), MAXVAL(arr(:,j)), mype, ix, jy, nlev
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
      MINVAL(arr), MAXVAL(arr), mype, ix, jy, nlev
CALL umPrint(umMessage, src='asad_range_mod')

END SUBROUTINE asad_range_print_1d

END MODULE asad_range_mod

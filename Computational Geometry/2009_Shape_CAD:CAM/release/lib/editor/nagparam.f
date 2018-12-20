C************************************************************************
C *									*
C			Copyright (C) 1992 by
C	Massachusetts Institute of Technology, Cambridge, MA
C			 All rights reserved
C *									*
C ************************************************************************


      SUBROUTINE NAGPARAM(I)
      CHARACTER*40 STRING
      INTEGER I

      IF (I.EQ.0) THEN
         STRING = 'Nolist'
      ELSE IF (I.EQ.1) THEN
         STRING = 'Minor Print Level = 0'
C         STRING = 'Minor Print Level = 10'
      ELSE IF (I.EQ.2) THEN
         STRING = 'Major Print Level = 0'
C         STRING = 'Major Print Level = 20'
      ELSE IF (I.EQ.3) THEN
         STRING = 'Optimality tolerance = 1D-5'
      ELSE IF (I.EQ.4) THEN
         STRING = 'Linesearch tolerance = 1D-5'
      ELSE IF (I.EQ.5) THEN
         STRING = 'Function precision = 1D-6'
      ELSE IF (I.EQ.6) THEN
         STRING = 'Non Feas Tol = 1D-6'
      ELSE IF (I.EQ.7) THEN
         STRING = 'Verify Level = 0'
      ELSE IF (I.EQ.8) THEN
         STRING = 'Difference Interval = 1D-6'
      ELSE IF (I.EQ.9) THEN
         STRING = 'Central difference Interval = 1D-6'
      ELSE IF (I.EQ.10) THEN
         STRING = 'Derivative Level = 3'
      END IF

      CALL E04UEF(STRING)

      RETURN
      END

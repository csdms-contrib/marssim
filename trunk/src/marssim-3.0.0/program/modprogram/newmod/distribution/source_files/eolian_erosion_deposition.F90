    !     ####################################################################################################
    !     this is one of several source files for the marssim landform evolution model
    !     copyright (c) 2009 alan d. howard
    !     developer can be contacted by ah6p`virginia.edu and department of environmental sciences, p.o. box 400123,
    !                university of virginia, charlottesville, va 22904-4123
    !     this program is free software; you can redistribute it and/or modify it under the terms of the gnu general public license
    !       as published by the free software foundation; either version 2 of the
    !       license, or (at your option) any later version.
    !     this program is distributed in the hope that it will be useful, but without any warranty;
    !        without even the implied warranty of merchantability or fitness for a particular purpose. see the gnu
    !        general public license for more details.
    !      you should have received a copy of the gnu general public license along with this program; if not, write to
    !        the free software foundation, inc., 51 franklin street, fifth floor, boston, ma 02110-1301 usa.
    !        a web link:   http://www.gnu.org/licenses/gpl-2.0.txt
    !     ####################################################################################################
    !       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      *********************************************************************
    !       Do eolian deposion and erosion based upon an "exposure" index which is
    !       positive if the location is on a relative high, zero on a level plain,
    !       and negative in depressions or valleys
    !       For examples of use and model presentation see:
    !          forsberg-taylor, n. k., howard, a. d., and craddock, r. a., 2004,
    !          j. geophys. res., 109, e05002, doi:10.1029/2004je002242
    !       and
    !          luo, w., and howard, a. d., 2008, j. geophysical research, 113,
    !          e05002,doi:10.1029/2007je002981
    !       Two variant subroutines to define relative exposure are exposure and totalexposure.
    !         In the exposure subroutine only cells visible to the local cell are used in computing
    !           the weighted exposure
    !         In the totalexposure subroutine all cells within the calculation window are used to
    !           compute exposure
    !       this routine is not mass-conservative.  It assumes that deposition is primarily by
    !       suspension from imported sediment and eroded sediment is exported as suspended
    !       material.  Therefore it is not suitable for modeling eolian deposition and
    !       erosion by saltation.
    !       If default_eolian_process is true, deposition occurs normal to the surface otherwise
    !         it is vertically-directed (without the inverse cosine correction)
    !       If eolian erosion occurs and fluvial and slope modeling is also used, then
    !        the rate of erosion is dependent upon whether bedrock or regolith occurs at the surface
    !   MODIFIES: ERODE_SLOPE, ELEVATION, CUMULATIVE_EOLIAN_CHANGE, IS_ROCK_SURFACE
    !             SEDIMENT_BASE
    !   CALLS: TOTAL_EXPOSURE, EXPOSURE
    !      *********************************************************************
    SUBROUTINE DO_EOLIAN_CHANGE()
        USE ERODE_GLOBALS
        USE EOLIAN_GLOBALS
        IMPLICIT NONE
        INTEGER :: I,J
        REAL (8) :: EXPOSE,ABSCHANGE,RATEMULT
        REAL (8) :: TEMP1,EOLRATE,MAXVAL,MINVAL
        REAL (8) :: MINEOLIAN,MAXEOLIAN,ETEMP,MINEXPOSERATE,MAXEXPOSERATE
        REAL (8) :: MINEXPOSE,MAXEXPOSE,AVGEXPOSE,EXPOSESSQ
        MAXVAL=-1.0E+25
        MINVAL=-MAXVAL
        MINEXPOSE=1.0E+25
        MAXEXPOSE=-MINEXPOSE
        AVGEXPOSE=0.0
        EXPOSESSQ=0.0
        MINEXPOSERATE=1.0E+25
        MAXEXPOSERATE=-MINEXPOSERATE
        MINEOLIAN=MINVAL
        MAXEOLIAN=MAXVAL
        L310: DO  J=1,MY
            M310: DO  I=1,MX
                IF (USE_TOTAL_EXPOSURE) THEN
                    CALL TOTAL_EXPOSURE(I,J,EXPOSE)
                ELSE
                    CALL EXPOSURE(I,J,EXPOSE)
                ENDIF
                AVGEXPOSE=AVGEXPOSE+EXPOSE
                EXPOSESSQ=EXPOSESSQ+EXPOSE**2
                TEMP1=EOLIAN_CONSTANT_3*(EXPOSE-EXPOSURE_50_PERCENT)
                EOLRATE=EOLIAN_CONSTANT_2+EOLIAN_CONSTANT_1*TEMP1/DSQRT(1.0+TEMP1**2)
                IF (FLUVIAL_AND_SLOPE_MODELING) THEN
                    IF (IS_ROCK_SURFACE(I,J)) THEN
                        IF (EOLRATE < 0.0) THEN
                            ERODE_SLOPE(I,J)=EOLRATE/REGOLITH_CRITICAL_SHEAR_FACTOR
                        ELSE
                            ERODE_SLOPE(I,J)=EOLRATE
                        ENDIF
                    ELSE
                        IF (DEFAULT_EOLIAN_PROCESS) THEN
                            ERODE_SLOPE(I,J)=DSQRT(D8_GRADIENT(I,J)**2+1.0)*EOLRATE
                        ELSE
                            ERODE_SLOPE(I,J)=EOLRATE
                        ENDIF
                    ENDIF
                ELSE
                    IF (DEFAULT_EOLIAN_PROCESS) THEN
                        ERODE_SLOPE(I,J)=DSQRT(D8_GRADIENT(I,J)**2+1.0)*EOLRATE
                    ELSE
                        ERODE_SLOPE(I,J)=EOLRATE
                    ENDIF
                ENDIF
                IF (.NOT.FLUVIAL_AND_SLOPE_MODELING) THEN
                    IF ((ERODE_SLOPE(I,J) < 0.0).AND.(ELEVATION(I,J) < INITIAL_ELEVATION(I,J)))&
                    THEN
                        ERODE_SLOPE(I,J)=ERODE_SLOPE(I,J)/REGOLITH_CRITICAL_SHEAR_FACTOR
                    ENDIF
                ENDIF
                IF (EXPOSE > MAXEXPOSE) THEN
                    MAXEXPOSE=EXPOSE
                    MAXEXPOSERATE=ERODE_SLOPE(I,J)
                ENDIF
                IF (EXPOSE < MINEXPOSE) THEN
                    MINEXPOSE=EXPOSE
                    MINEXPOSERATE=ERODE_SLOPE(I,J)
                ENDIF
                IF (ERODE_SLOPE(I,J) < MINEOLIAN) THEN
                    MINEOLIAN=ERODE_SLOPE(I,J)
                ENDIF
                IF (ERODE_SLOPE(I,J) > MAXEOLIAN) THEN
                    MAXEOLIAN=ERODE_SLOPE(I,J)
                ENDIF
            ENDDO M310
        ENDDO L310
        MAXIMUM_ELEVATION_CHANGE=0.0
        DO  J=1,MY
            DO  I=1,MX
                ABSCHANGE=DABS(ERODE_SLOPE(I,J))
                IF (ABSCHANGE > MAXIMUM_ELEVATION_CHANGE) MAXIMUM_ELEVATION_CHANGE=ABSCHANGE
            ENDDO
        ENDDO
        RATEMULT=EOLIAN_TIME_INCREMENT
        EXPOSESSQ=DSQRT((EXPOSESSQ-AVGEXPOSE**2/(MX*MY))/&
        (MX*MY-1))
        AVGEXPOSE=AVGEXPOSE/(MX*MY)
        WRITE(OUTHIST,11310) MINEOLIAN,MAXEOLIAN,EOLIAN_TIME_INCREMENT
        WRITE(*,11310) MINEOLIAN,MAXEOLIAN,EOLIAN_TIME_INCREMENT
        11310         FORMAT(' MINEOLIAN=',G12.5,' MAXEOLIAN=',G12.5, &
        ' EOLIAN_TIME_INCREMENT=',G12.5)
        WRITE(OUTHIST,21310) MINEXPOSE,MINEXPOSERATE, &
        MAXEXPOSE,MAXEXPOSERATE,  &
        AVGEXPOSE,EXPOSESSQ
        WRITE(*,21310) MINEXPOSE,MINEXPOSERATE, &
        MAXEXPOSE,MAXEXPOSERATE, &
        AVGEXPOSE,EXPOSESSQ
        21310         FORMAT('MINEXPOSE=',G12.5,' MINEXPOSERATE=',G12.5,/,&
        ' MAXEXPOSE=',G12.5,' MAXEXPOSERATE=',G12.5,        &
        ' AVG=',G12.5,' SD=',G12.5)
        DO  J=1,MY
            DO  I=1,MX
                ETEMP=ERODE_SLOPE(I,J)*RATEMULT
                ELEVATION(I,J)=ELEVATION(I,J)+ETEMP
                CUMULATIVE_EOLIAN_CHANGE(I,J)=CUMULATIVE_EOLIAN_CHANGE(I,J)+ETEMP
                ERODE_SLOPE(I,J)=0.0
            ENDDO
        ENDDO
        IF (FLUVIAL_AND_SLOPE_MODELING) THEN
            DO  J=1,MY
                DO  I=1,MX
                    IF (ELEVATION(I,J) > SEDIMENT_BASE(I,J)) THEN
                        IS_ROCK_SURFACE(I,J)=.FALSE.
                    ELSE
                        IS_ROCK_SURFACE(I,J)=.TRUE.
                        SEDIMENT_BASE(I,J)=ELEVATION(I,J)
                    ENDIF
                ENDDO
            ENDDO
        ENDIF
        RETURN
    END
    !       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE EXPOSURE(I,J,EXPOSE)
        USE ERODE_GLOBALS
        USE EOLIAN_GLOBALS
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: I,J
        INTEGER :: K,ILOC,JLOC
        REAL (8) :: ECOMP,ELOC,WLOC,LSLOPE,WEIGHTXUM,WEIGHTDUM
        REAL (8) :: SLOPE1SUM,SLOPE2SUM,SLOPE3SUM,SLOPE4SUM
        REAL (8) :: SLOPE5SUM,SLOPE6SUM,SLOPE7SUM,SLOPE8SUM, MAXSLOPE
        REAL (8), INTENT(OUT) :: EXPOSE
        LOGICAL DOCALC
        DOCALC=.TRUE.
        ECOMP=ELEVATION(I,J)
        ILOC=I+1
        IF (ILOC > MX) THEN
            IF(IS_X_PERIODIC) THEN
                ILOC=1
            ELSE
                SLOPE1SUM=0.0
                DOCALC=.FALSE.
            ENDIF
        ENDIF
        IF (DOCALC) THEN
            ELOC=ELEVATION(ILOC,J)
            WLOC=WEIGHTX(1)
            LSLOPE=ELOC-ECOMP
            MAXSLOPE=LSLOPE
            SLOPE1SUM=LSLOPE*WLOC
            WEIGHTXUM=WLOC
            L100: DO K=2,MAXIMUM_WEIGHT_DISTANCE
                ILOC=I+K
                IF (ILOC > MX) THEN
                    IF (IS_X_PERIODIC) THEN
                        ILOC=ILOC-MX
                    ELSE
                        EXIT L100
                    ENDIF
                ENDIF
                ELOC=ELEVATION(ILOC,J)
                LSLOPE=(ELOC-ECOMP)/K
                IF (LSLOPE > MAXSLOPE) THEN
                    MAXSLOPE=LSLOPE
                    WLOC=WEIGHTX(K)
                    SLOPE1SUM=SLOPE1SUM+LSLOPE*WLOC
                    WEIGHTXUM=WEIGHTXUM+WLOC
                ENDIF
            ENDDO L100
            SLOPE1SUM=SLOPE1SUM/WEIGHTXUM
        ELSE
            DOCALC=.TRUE.
        ENDIF
        ILOC=I-1
        IF (ILOC < 1) THEN
            IF(IS_X_PERIODIC) THEN
                ILOC=MX
            ELSE
                SLOPE2SUM=0.0
                DOCALC=.FALSE.
            ENDIF
        ENDIF
        IF (DOCALC) THEN
            ELOC=ELEVATION(ILOC,J)
            WLOC=WEIGHTX(1)
            LSLOPE=ELOC-ECOMP
            MAXSLOPE=LSLOPE
            SLOPE2SUM=LSLOPE*WLOC
            WEIGHTXUM=WLOC
            L110: DO K=2,MAXIMUM_WEIGHT_DISTANCE
                ILOC=I-K
                IF (ILOC < 1) THEN
                    IF (IS_X_PERIODIC) THEN
                        ILOC=ILOC+MX
                    ELSE
                        EXIT L110
                    ENDIF
                ENDIF
                ELOC=ELEVATION(ILOC,J)
                LSLOPE=(ELOC-ECOMP)/K
                IF (LSLOPE > MAXSLOPE) THEN
                    MAXSLOPE=LSLOPE
                    WLOC=WEIGHTX(K)
                    SLOPE2SUM=SLOPE2SUM+LSLOPE*WLOC
                    WEIGHTXUM=WEIGHTXUM+WLOC
                ENDIF
            ENDDO L110
            SLOPE2SUM=SLOPE2SUM/WEIGHTXUM
        ELSE
            DOCALC=.TRUE.
        ENDIF
        JLOC=J+1
        IF (JLOC > MY) THEN
            IF(IS_Y_PERIODIC) THEN
                JLOC=1
            ELSE
                SLOPE3SUM=0.0
                DOCALC=.FALSE.
            ENDIF
        ENDIF
        IF (DOCALC) THEN
            ELOC=ELEVATION(I,JLOC)
            WLOC=WEIGHTX(1)
            LSLOPE=ELOC-ECOMP
            MAXSLOPE=LSLOPE
            SLOPE3SUM=LSLOPE*WLOC
            WEIGHTXUM=WLOC
            L120: DO K=2,MAXIMUM_WEIGHT_DISTANCE
                JLOC=J+K
                IF (JLOC > MY) THEN
                    IF (IS_Y_PERIODIC) THEN
                        JLOC=JLOC-MY
                    ELSE
                        EXIT L120
                    ENDIF
                ENDIF
                ELOC=ELEVATION(I,JLOC)
                LSLOPE=(ELOC-ECOMP)/K
                IF (LSLOPE > MAXSLOPE) THEN
                    MAXSLOPE=LSLOPE
                    WLOC=WEIGHTX(K)
                    SLOPE3SUM=SLOPE3SUM+LSLOPE*WLOC
                    WEIGHTXUM=WEIGHTXUM+WLOC
                ENDIF
            ENDDO L120
            SLOPE3SUM=SLOPE3SUM/WEIGHTXUM
        ELSE
            DOCALC=.TRUE.
        ENDIF
        JLOC=J-1
        IF (JLOC < 1) THEN
            IF(IS_Y_PERIODIC) THEN
                JLOC=MY
            ELSE
                SLOPE4SUM=0.0
                DOCALC=.FALSE.
            ENDIF
        ENDIF
        IF (DOCALC) THEN
            ELOC=ELEVATION(I,JLOC)
            WLOC=WEIGHTX(1)
            LSLOPE=ELOC-ECOMP
            MAXSLOPE=LSLOPE
            SLOPE4SUM=LSLOPE*WLOC
            WEIGHTXUM=WLOC
            L130: DO K=2,MAXIMUM_WEIGHT_DISTANCE
                JLOC=J-K
                IF (JLOC < 1) THEN
                    IF (IS_Y_PERIODIC) THEN
                        JLOC=JLOC+MY
                    ELSE
                        EXIT L130
                    ENDIF
                ENDIF
                ELOC=ELEVATION(I,JLOC)
                LSLOPE=(ELOC-ECOMP)/K
                IF (LSLOPE > MAXSLOPE) THEN
                    MAXSLOPE=LSLOPE
                    WLOC=WEIGHTX(K)
                    SLOPE4SUM=SLOPE4SUM+LSLOPE*WLOC
                    WEIGHTXUM=WEIGHTXUM+WLOC
                ENDIF
            ENDDO L130
            SLOPE4SUM=SLOPE4SUM/WEIGHTXUM
        ELSE
            DOCALC=.TRUE.
        ENDIF
        ILOC=I+1
        JLOC=J+1
        IF (ILOC > MX) THEN
            IF(IS_X_PERIODIC) THEN
                ILOC=1
            ELSE
                SLOPE5SUM=0.0
                DOCALC=.FALSE.
            ENDIF
        ENDIF
        IF (JLOC > MY) THEN
            IF(IS_Y_PERIODIC) THEN
                JLOC=1
            ELSE
                SLOPE5SUM=0.0
                DOCALC=.FALSE.
            ENDIF
        ENDIF
        IF (DOCALC) THEN
            ELOC=ELEVATION(ILOC,JLOC)
            WLOC=WEIGHTD(1)
            LSLOPE=(ELOC-ECOMP)
            MAXSLOPE=LSLOPE
            SLOPE5SUM=LSLOPE*WLOC
            WEIGHTDUM=WLOC
            L140: DO K=2,MAXIMUM_WEIGHT_DISTANCE
                ILOC=I+K
                JLOC=J+K
                IF (ILOC > MX) THEN
                    IF (IS_X_PERIODIC) THEN
                        ILOC=ILOC-MX
                    ELSE
                        EXIT L140
                    ENDIF
                ENDIF
                IF (JLOC > MY) THEN
                    IF (IS_Y_PERIODIC) THEN
                        JLOC=JLOC-MY
                    ELSE
                        EXIT L140
                    ENDIF
                ENDIF
                ELOC=ELEVATION(ILOC,JLOC)
                LSLOPE=(ELOC-ECOMP)/K
                IF (LSLOPE > MAXSLOPE) THEN
                    MAXSLOPE=LSLOPE
                    WLOC=WEIGHTD(K)
                    SLOPE5SUM=SLOPE5SUM+LSLOPE*WLOC
                    WEIGHTDUM=WEIGHTDUM+WLOC
                ENDIF
            ENDDO L140
            SLOPE5SUM=SLOPE5SUM/WEIGHTDUM
        ELSE
            DOCALC=.TRUE.
        ENDIF
        JLOC=J-1
        ILOC=I-1
        IF (ILOC < 1) THEN
            IF(IS_X_PERIODIC) THEN
                ILOC=MX
            ELSE
                SLOPE6SUM=0.0
                DOCALC=.FALSE.
            ENDIF
        ENDIF
        IF (JLOC < 1) THEN
            IF(IS_Y_PERIODIC) THEN
                JLOC=MY
            ELSE
                SLOPE6SUM=0.0
                DOCALC=.FALSE.
            ENDIF
        ENDIF
        IF (DOCALC) THEN
            ELOC=ELEVATION(ILOC,JLOC)
            WLOC=WEIGHTD(1)
            LSLOPE=ELOC-ECOMP
            MAXSLOPE=LSLOPE
            SLOPE6SUM=LSLOPE*WLOC
            WEIGHTDUM=WLOC
            L150: DO K=2,MAXIMUM_WEIGHT_DISTANCE
                ILOC=I-K
                JLOC=J-K
                IF (ILOC < 1) THEN
                    IF (IS_X_PERIODIC) THEN
                        ILOC=ILOC+MX
                    ELSE
                        EXIT L150
                    ENDIF
                ENDIF
                IF (JLOC < 1) THEN
                    IF (IS_Y_PERIODIC) THEN
                        JLOC=JLOC+MY
                    ELSE
                        EXIT L150
                    ENDIF
                ENDIF
                ELOC=ELEVATION(ILOC,JLOC)
                LSLOPE=(ELOC-ECOMP)/K
                IF (LSLOPE > MAXSLOPE) THEN
                    MAXSLOPE=LSLOPE
                    WLOC=WEIGHTD(K)
                    SLOPE6SUM=SLOPE6SUM+LSLOPE*WLOC
                    WEIGHTDUM=WEIGHTDUM+WLOC
                ENDIF
            ENDDO L150
            SLOPE6SUM=SLOPE6SUM/WEIGHTDUM
        ELSE
            DOCALC=.TRUE.
        ENDIF
        JLOC=J+1
        ILOC=I-1
        IF (JLOC > MY) THEN
            IF(IS_Y_PERIODIC) THEN
                JLOC=1
            ELSE
                SLOPE7SUM=0.0
                DOCALC=.FALSE.
            ENDIF
        ENDIF
        IF (ILOC < 1) THEN
            IF(IS_X_PERIODIC) THEN
                ILOC=MX
            ELSE
                SLOPE7SUM=0.0
                DOCALC=.FALSE.
            ENDIF
        ENDIF
        IF (DOCALC) THEN
            ELOC=ELEVATION(ILOC,JLOC)
            WLOC=WEIGHTD(1)
            LSLOPE=ELOC-ECOMP
            MAXSLOPE=LSLOPE
            SLOPE7SUM=LSLOPE*WLOC
            WEIGHTDUM=WLOC
            L160: DO  K=2,MAXIMUM_WEIGHT_DISTANCE
                JLOC=J+K
                ILOC=I-K
                IF (JLOC > MY) THEN
                    IF (IS_Y_PERIODIC) THEN
                        JLOC=JLOC-MY
                    ELSE
                        EXIT L160
                    ENDIF
                ENDIF
                IF (ILOC < 1) THEN
                    IF (IS_X_PERIODIC) THEN
                        ILOC=ILOC+MX
                    ELSE
                        EXIT L160
                    ENDIF
                ENDIF
                ELOC=ELEVATION(ILOC,JLOC)
                LSLOPE=(ELOC-ECOMP)/K
                IF (LSLOPE > MAXSLOPE) THEN
                    MAXSLOPE=LSLOPE
                    WLOC=WEIGHTD(K)
                    SLOPE7SUM=SLOPE7SUM+LSLOPE*WLOC
                    WEIGHTDUM=WEIGHTDUM+WLOC
                ENDIF
            ENDDO L160
            SLOPE7SUM=SLOPE7SUM/WEIGHTDUM
        ELSE
            DOCALC=.TRUE.
        ENDIF
        JLOC=J-1
        ILOC=I+1
        IF (JLOC < 1) THEN
            IF(IS_Y_PERIODIC) THEN
                JLOC=MY
            ELSE
                SLOPE8SUM=0.0
                DOCALC=.FALSE.
            ENDIF
        ENDIF
        IF (ILOC > MX) THEN
            IF(IS_X_PERIODIC) THEN
                ILOC=1
            ELSE
                SLOPE8SUM=0.0
                DOCALC=.FALSE.
            ENDIF
        ENDIF
        IF (DOCALC) THEN
            ELOC=ELEVATION(ILOC,JLOC)
            WLOC=WEIGHTD(1)
            LSLOPE=ELOC-ECOMP
            MAXSLOPE=LSLOPE
            SLOPE8SUM=LSLOPE*WLOC
            WEIGHTDUM=WLOC
            L170: DO K=2,MAXIMUM_WEIGHT_DISTANCE
                JLOC=J-K
                ILOC=I+K
                IF (JLOC < 1) THEN
                    IF (IS_Y_PERIODIC) THEN
                        JLOC=JLOC+MY
                    ELSE
                        EXIT L170
                    ENDIF
                ENDIF
                IF (ILOC > MX) THEN
                    IF (IS_X_PERIODIC) THEN
                        ILOC=ILOC-MX
                    ELSE
                        EXIT L170
                    ENDIF
                ENDIF
                ELOC=ELEVATION(ILOC,JLOC)
                LSLOPE=(ELOC-ECOMP)/K
                IF (LSLOPE > MAXSLOPE) THEN
                    MAXSLOPE=LSLOPE
                    WLOC=WEIGHTD(K)
                    SLOPE8SUM=SLOPE8SUM+LSLOPE*WLOC
                    WEIGHTDUM=WEIGHTDUM+WLOC
                ENDIF
            ENDDO L170
            SLOPE8SUM=SLOPE8SUM/WEIGHTDUM
        ELSE
            DOCALC=.TRUE.
        ENDIF
        EXPOSE=(SLOPE1SUM+SLOPE2SUM+SLOPE3SUM+SLOPE4SUM+SLOPE5SUM+ &
        SLOPE6SUM+SLOPE7SUM+SLOPE8SUM)/(8.0*CELL_SIZE)
        RETURN
    END

    !       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE TOTAL_EXPOSURE(I,J,EXPOSE)
        USE ERODE_GLOBALS
        USE EOLIAN_GLOBALS
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: I,J
        INTEGER :: K,ILOC,JLOC
        REAL (8) :: ECOMP,ELOC,WLOC,LSLOPE,WEIGHTXUM,WEIGHTDUM
        REAL (8) :: SLOPE1SUM,SLOPE2SUM,SLOPE3SUM,SLOPE4SUM
        REAL (8) :: SLOPE5SUM,SLOPE6SUM,SLOPE7SUM,SLOPE8SUM
        REAL (8), INTENT(OUT) :: EXPOSE
        LOGICAL DOCALC
        DOCALC=.TRUE.
        ECOMP=ELEVATION(I,J)
        ILOC=I+1
        IF (ILOC > MX) THEN
            IF(IS_X_PERIODIC) THEN
                ILOC=1
            ELSE
                SLOPE1SUM=0.0
                DOCALC=.FALSE.
            ENDIF
        ENDIF
        IF (DOCALC) THEN
            ELOC=ELEVATION(ILOC,J)
            WLOC=WEIGHTX(1)
            LSLOPE=ELOC-ECOMP
            SLOPE1SUM=LSLOPE*WLOC
            WEIGHTXUM=WLOC
            L100: DO K=2,MAXIMUM_WEIGHT_DISTANCE
                ILOC=I+K
                IF (ILOC > MX) THEN
                    IF (IS_X_PERIODIC) THEN
                        ILOC=ILOC-MX
                    ELSE
                        EXIT L100
                    ENDIF
                ENDIF
                ELOC=ELEVATION(ILOC,J)
                LSLOPE=(ELOC-ECOMP)/K
                WLOC=WEIGHTX(K)
                SLOPE1SUM=SLOPE1SUM+LSLOPE*WLOC
                WEIGHTXUM=WEIGHTXUM+WLOC
            ENDDO L100
            SLOPE1SUM=SLOPE1SUM/WEIGHTXUM
        ELSE
            DOCALC=.TRUE.
        ENDIF
        ILOC=I-1
        IF (ILOC < 1) THEN
            IF(IS_X_PERIODIC) THEN
                ILOC=MX
            ELSE
                SLOPE2SUM=0.0
                DOCALC=.FALSE.
            ENDIF
        ENDIF
        IF (DOCALC) THEN
            ELOC=ELEVATION(ILOC,J)
            WLOC=WEIGHTX(1)
            LSLOPE=ELOC-ECOMP
            SLOPE2SUM=LSLOPE*WLOC
            WEIGHTXUM=WLOC
            L110: DO K=2,MAXIMUM_WEIGHT_DISTANCE
                ILOC=I-K
                IF (ILOC < 1) THEN
                    IF (IS_X_PERIODIC) THEN
                        ILOC=ILOC+MX
                    ELSE
                        EXIT L110
                    ENDIF
                ENDIF
                ELOC=ELEVATION(ILOC,J)
                LSLOPE=(ELOC-ECOMP)/K
                WLOC=WEIGHTX(K)
                SLOPE2SUM=SLOPE2SUM+LSLOPE*WLOC
                WEIGHTXUM=WEIGHTXUM+WLOC
            ENDDO L110
            SLOPE2SUM=SLOPE2SUM/WEIGHTXUM
        ELSE
            DOCALC=.TRUE.
        ENDIF
        JLOC=J+1
        IF (JLOC > MY) THEN
            IF(IS_Y_PERIODIC) THEN
                JLOC=1
            ELSE
                SLOPE3SUM=0.0
                DOCALC=.FALSE.
            ENDIF
        ENDIF
        IF (DOCALC) THEN
            ELOC=ELEVATION(I,JLOC)
            WLOC=WEIGHTX(1)
            LSLOPE=ELOC-ECOMP
            SLOPE3SUM=LSLOPE*WLOC
            WEIGHTXUM=WLOC
            L120: DO K=2,MAXIMUM_WEIGHT_DISTANCE
                JLOC=J+K
                IF (JLOC > MY) THEN
                    IF (IS_Y_PERIODIC) THEN
                        JLOC=JLOC-MY
                    ELSE
                        EXIT L120
                    ENDIF
                ENDIF
                ELOC=ELEVATION(I,JLOC)
                LSLOPE=(ELOC-ECOMP)/K
                WLOC=WEIGHTX(K)
                SLOPE3SUM=SLOPE3SUM+LSLOPE*WLOC
                WEIGHTXUM=WEIGHTXUM+WLOC
            ENDDO L120
            SLOPE3SUM=SLOPE3SUM/WEIGHTXUM
        ELSE
            DOCALC=.TRUE.
        ENDIF
        JLOC=J-1
        IF (JLOC < 1) THEN
            IF(IS_Y_PERIODIC) THEN
                JLOC=MY
            ELSE
                SLOPE4SUM=0.0
                DOCALC=.FALSE.
            ENDIF
        ENDIF
        IF (DOCALC) THEN
            ELOC=ELEVATION(I,JLOC)
            WLOC=WEIGHTX(1)
            LSLOPE=ELOC-ECOMP
            SLOPE4SUM=LSLOPE*WLOC
            WEIGHTXUM=WLOC
            L130: DO K=2,MAXIMUM_WEIGHT_DISTANCE
                JLOC=J-K
                IF (JLOC < 1) THEN
                    IF (IS_Y_PERIODIC) THEN
                        JLOC=JLOC+MY
                    ELSE
                        EXIT L130
                    ENDIF
                ENDIF
                ELOC=ELEVATION(I,JLOC)
                LSLOPE=(ELOC-ECOMP)/K
                WLOC=WEIGHTX(K)
                SLOPE4SUM=SLOPE4SUM+LSLOPE*WLOC
                WEIGHTXUM=WEIGHTXUM+WLOC
            ENDDO L130
            SLOPE4SUM=SLOPE4SUM/WEIGHTXUM
        ELSE
            DOCALC=.TRUE.
        ENDIF
        ILOC=I+1
        JLOC=J+1
        IF (ILOC > MX) THEN
            IF(IS_X_PERIODIC) THEN
                ILOC=1
            ELSE
                SLOPE5SUM=0.0
                DOCALC=.FALSE.
            ENDIF
        ENDIF
        IF (JLOC > MY) THEN
            IF(IS_Y_PERIODIC) THEN
                JLOC=1
            ELSE
                SLOPE5SUM=0.0
                DOCALC=.FALSE.
            ENDIF
        ENDIF
        IF (DOCALC) THEN
            ELOC=ELEVATION(ILOC,JLOC)
            WLOC=WEIGHTD(1)
            LSLOPE=(ELOC-ECOMP)
            SLOPE5SUM=LSLOPE*WLOC
            WEIGHTDUM=WLOC
            L140: DO K=2,MAXIMUM_WEIGHT_DISTANCE
                ILOC=I+K
                JLOC=J+K
                IF (ILOC > MX) THEN
                    IF (IS_X_PERIODIC) THEN
                        ILOC=ILOC-MX
                    ELSE
                        EXIT L140
                    ENDIF
                ENDIF
                IF (JLOC > MY) THEN
                    IF (IS_Y_PERIODIC) THEN
                        JLOC=JLOC-MY
                    ELSE
                        EXIT L140
                    ENDIF
                ENDIF
                ELOC=ELEVATION(ILOC,JLOC)
                LSLOPE=(ELOC-ECOMP)/K
                WLOC=WEIGHTD(K)
                SLOPE5SUM=SLOPE5SUM+LSLOPE*WLOC
                WEIGHTDUM=WEIGHTDUM+WLOC
            ENDDO L140
            SLOPE5SUM=SLOPE5SUM/WEIGHTDUM
        ELSE
            DOCALC=.TRUE.
        ENDIF
        JLOC=J-1
        ILOC=I-1
        IF (ILOC < 1) THEN
            IF(IS_X_PERIODIC) THEN
                ILOC=MX
            ELSE
                SLOPE6SUM=0.0
                DOCALC=.FALSE.
            ENDIF
        ENDIF
        IF (JLOC < 1) THEN
            IF(IS_Y_PERIODIC) THEN
                JLOC=MY
            ELSE
                SLOPE6SUM=0.0
                DOCALC=.FALSE.
            ENDIF
        ENDIF
        IF (DOCALC) THEN
            ELOC=ELEVATION(ILOC,JLOC)
            WLOC=WEIGHTD(1)
            LSLOPE=ELOC-ECOMP
            SLOPE6SUM=LSLOPE*WLOC
            WEIGHTDUM=WLOC
            L150: DO K=2,MAXIMUM_WEIGHT_DISTANCE
                ILOC=I-K
                JLOC=J-K
                IF (ILOC < 1) THEN
                    IF (IS_X_PERIODIC) THEN
                        ILOC=ILOC+MX
                    ELSE
                        EXIT L150
                    ENDIF
                ENDIF
                IF (JLOC < 1) THEN
                    IF (IS_Y_PERIODIC) THEN
                        JLOC=JLOC+MY
                    ELSE
                        EXIT L150
                    ENDIF
                ENDIF
                ELOC=ELEVATION(ILOC,JLOC)
                LSLOPE=(ELOC-ECOMP)/K
                WLOC=WEIGHTD(K)
                SLOPE6SUM=SLOPE6SUM+LSLOPE*WLOC
                WEIGHTDUM=WEIGHTDUM+WLOC
            ENDDO L150
            SLOPE6SUM=SLOPE6SUM/WEIGHTDUM
        ELSE
            DOCALC=.TRUE.
        ENDIF
        JLOC=J+1
        ILOC=I-1
        IF (JLOC > MY) THEN
            IF(IS_Y_PERIODIC) THEN
                JLOC=1
            ELSE
                SLOPE7SUM=0.0
                DOCALC=.FALSE.
            ENDIF
        ENDIF
        IF (ILOC < 1) THEN
            IF(IS_X_PERIODIC) THEN
                ILOC=MX
            ELSE
                SLOPE7SUM=0.0
                DOCALC=.FALSE.
            ENDIF
        ENDIF
        IF (DOCALC) THEN
            ELOC=ELEVATION(ILOC,JLOC)
            WLOC=WEIGHTD(1)
            LSLOPE=ELOC-ECOMP
            SLOPE7SUM=LSLOPE*WLOC
            WEIGHTDUM=WLOC
            L160: DO  K=2,MAXIMUM_WEIGHT_DISTANCE
                JLOC=J+K
                ILOC=I-K
                IF (JLOC > MY) THEN
                    IF (IS_Y_PERIODIC) THEN
                        JLOC=JLOC-MY
                    ELSE
                        EXIT L160
                    ENDIF
                ENDIF
                IF (ILOC < 1) THEN
                    IF (IS_X_PERIODIC) THEN
                        ILOC=ILOC+MX
                    ELSE
                        EXIT L160
                    ENDIF
                ENDIF
                ELOC=ELEVATION(ILOC,JLOC)
                LSLOPE=(ELOC-ECOMP)/K
                WLOC=WEIGHTD(K)
                SLOPE7SUM=SLOPE7SUM+LSLOPE*WLOC
                WEIGHTDUM=WEIGHTDUM+WLOC
            ENDDO L160
            SLOPE7SUM=SLOPE7SUM/WEIGHTDUM
        ELSE
            DOCALC=.TRUE.
        ENDIF
        JLOC=J-1
        ILOC=I+1
        IF (JLOC < 1) THEN
            IF(IS_Y_PERIODIC) THEN
                JLOC=MY
            ELSE
                SLOPE8SUM=0.0
                DOCALC=.FALSE.
            ENDIF
        ENDIF
        IF (ILOC > MX) THEN
            IF(IS_X_PERIODIC) THEN
                ILOC=1
            ELSE
                SLOPE8SUM=0.0
                DOCALC=.FALSE.
            ENDIF
        ENDIF
        IF (DOCALC) THEN
            ELOC=ELEVATION(ILOC,JLOC)
            WLOC=WEIGHTD(1)
            LSLOPE=ELOC-ECOMP
            SLOPE8SUM=LSLOPE*WLOC
            WEIGHTDUM=WLOC
            L170: DO K=2,MAXIMUM_WEIGHT_DISTANCE
                JLOC=J-K
                ILOC=I+K
                IF (JLOC < 1) THEN
                    IF (IS_Y_PERIODIC) THEN
                        JLOC=JLOC+MY
                    ELSE
                    ENDIF
                ENDIF
                IF (ILOC > MX) THEN
                    IF (IS_X_PERIODIC) THEN
                        ILOC=ILOC-MX
                    ELSE
                        EXIT L170
                    ENDIF
                ENDIF
                ELOC=ELEVATION(ILOC,JLOC)
                LSLOPE=(ELOC-ECOMP)/K
                WLOC=WEIGHTD(K)
                SLOPE8SUM=SLOPE8SUM+LSLOPE*WLOC
                WEIGHTDUM=WEIGHTDUM+WLOC
            ENDDO L170
            SLOPE8SUM=SLOPE8SUM/WEIGHTDUM
        ELSE
            DOCALC=.TRUE.
        ENDIF
        EXPOSE=(SLOPE1SUM+SLOPE2SUM+SLOPE3SUM+SLOPE4SUM+SLOPE5SUM+ &
        SLOPE6SUM+SLOPE7SUM+SLOPE8SUM)/(8.0*CELL_SIZE)
        RETURN
    END

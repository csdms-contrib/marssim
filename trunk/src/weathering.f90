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
    !       ********************************************************************
    !        calculate topgraphic divergence
    !    ACCESSES: MX,MY,IS_X_PERIODIC,IS_Y_PERIODIC,ELEVATION
    !    MODIFIES: DIVERGENCE
    !       ********************************************************************
    SUBROUTINE CALCULATE_DIVERGENCE()
        USE ERODE_GLOBALS
        IMPLICIT NONE
        INTEGER :: I,J,IE,IW,JN,JS
        DO  J=1,MY
            DO  I=1,MX
                IW=I-1
                IE=I+1
                IF (IS_X_PERIODIC) THEN
                    IF (IW < 1) IW=MX
                    IF (IE > MX) IE=1
                ELSE
                    IF (IW < 1) IW=2
                    IF (IE > MX) IE=MX-1
                ENDIF
                JN=J-1
                JS=J+1
                IF (IS_Y_PERIODIC) THEN
                    IF (JN < 1) JN=MY
                    IF (JS > MY) JS=1
                ELSE
                    IF (JN < 1) JN=2
                    IF (JS > MY) JS=MY-1
                ENDIF
                !      ********************************************************************
                !      Divergence based on 9 point formula with reflection at top
                !         and bottom boundaries and periodic left&right boundaries
                !      ********************************************************************
                DIVERGENCE(I,J)=(ELEVATION(IW,JN)+ELEVATION(IE,JN)+ELEVATION(IW,JS)+ELEVATION(IE,JS)+ &
                4.0*(ELEVATION(IE,J)+ELEVATION(IW,J)+ELEVATION(I,JN)+ELEVATION(I,JS))-20.0*ELEVATION(I,J))/ &
                (6.0*D2X)
            ENDDO
        ENDDO
        RETURN
    END
    !       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE DO_WEATHERING()
        USE ERODE_GLOBALS
        USE ACCRETION_GLOBALS
        IMPLICIT NONE
        !      ********************************************************************
        !       This subroutine determines the weathering rate on bare bedrock
        !         and weathering rate and regolith thickness on regolith covered
        !         slopes
        !       For regolith-mantled slopes the matrix regolith(i,j) is the regolith thickness
        !       For bedrock slopes regolith(i,j) is less than zero and is the rate of weathering of rock into regolith
        !   ACCESSES: MX, MY, MYY, SUBMERGED, ROCK_WEATHERING_RATE, WEATHERING_TERM_2
        !             SEEPAGE_WEATHERING_SCALING,GROUNDWATER_FLUX,SEEPAGE_WEATHERING_EXPONENT
        !             CRITICAL_BEDROCK_GRADIENT,D8_GRADIENT, WEATHER_MULT,TWO_TERM_WEATHERING
        !             RESISTANT_SURFACE_LAYER,RAD_ERODE,RAD_CONST,RELATIVE_RESISTANCE 
        !             TIME_INCREMENT,IS_DEPRESSION,WEATHER_DECAY_RATE,SEEPAGE_WEATHERING
        !             D8_GRADIENT,USE_SOLAR_EROSION,DIVERGENCE,
        !   MODIFIES: REGOLITH,IS_ROCK_SURFACE 
        !   CALLS: ROCK_MASS_WASTING, FIND_DEPRESSION
        !      ********************************************************************
        REAL (8) :: DIVERGE,MAXRATE,WTERM,FAILTERM,ROCK_MASS_WASTING,MAXCOSFACT,RAD_ERODE
        REAL (8) :: COSFACTOR,AVGWTERM,MAXWTERM,XAVGSEEP,MAXSEEP,XSEEP,XTOTAL
        REAL (8) :: FRACTBEDROCK,MINREGOLITH,MAXREGOLITH,AVGREGOLITH
        EXTERNAL ROCK_MASS_WASTING, RAD_ERODE
        LOGICAL :: IS_DEPRESSION,WRITEDETAIL,SKIPBEDROCK
        INTEGER :: I,J
        WRITEDETAIL=.FALSE.
        FRACTBEDROCK=0.0
        MINREGOLITH=1.0E+25
        MAXREGOLITH=-1.0E+25
        AVGREGOLITH=0.0
        XTOTAL=0.0
        MAXSEEP=-1.0E+25
        MAXWTERM=-1.0E+25
        XAVGSEEP=0.0
        MAXSEEP=0.0
        AVGWTERM=0.0
        MAXWTERM=0.0
        MAXCOSFACT=10.0
        TOTALRADERODE=0.0
        L110: DO  J=1,MYY
            M110: DO  I=1,MX
                SKIPBEDROCK=.FALSE.
                COSFACTOR=MIN(MAXCOSFACT,DSQRT(1.0+D8_GRADIENT(I,J)**2))
                IF (SUBMERGED(I,J)) CYCLE M110
                IF (IS_ROCK_SURFACE(I,J)) THEN
                    IF (TWO_TERM_WEATHERING) THEN
                        MAXRATE = COSFACTOR*(ROCK_WEATHERING_RATE - WEATHERING_TERM_2)
                    ELSE
                        MAXRATE = COSFACTOR*ROCK_WEATHERING_RATE
                    ENDIF
                    IF (SEEPAGE_WEATHERING) THEN
                        IF (SEEPAGE_AVERAGING) THEN
                            XSEEP=SEEPAGE_WEATHERING_SCALING*FILTERED_GROUNDWATER_FLUX(I,J)**SEEPAGE_WEATHERING_EXPONENT
                        ELSE
                            XSEEP=SEEPAGE_WEATHERING_SCALING*GROUNDWATER_FLUX(I,J)**SEEPAGE_WEATHERING_EXPONENT
                        ENDIF
                        IF (XSEEP > 0.0) THEN
                            IS_ROCK_SURFACE(I,J)=.FALSE.
                            REGOLITH(I,J)=XSEEP* TIME_INCREMENT/RELATIVE_RESISTANCE(I,J)
                            SKIPBEDROCK=.TRUE.
                        ENDIF
                        XAVGSEEP=XAVGSEEP+XSEEP
                        XTOTAL=XTOTAL+1.0
                        IF (XSEEP > MAXSEEP) MAXSEEP=XSEEP
                        MAXRATE=MAXRATE+XSEEP
                    ENDIF
                    !     *********************************************************************
                    !      Following is for bare bedrock slopes
                    !      The first part calculates a divergence-dependent weathering term
                    !      (if weather_divergence is greater than zero
                    !     *********************************************************************
                    IF (.NOT.SKIPBEDROCK) THEN
                        IF (WEATHER_DIVERGENCE > 0.0) THEN
                            DIVERGE=DEXP(-WEATHER_DIVERGENCE*DIVERGENCE(I,J))
                        ELSE
                        !     *********************************************************************
                        !      Just set divergence to unity if there is no divergence dependency
                        !     *********************************************************************
                            DIVERGE=1.0
                        ENDIF
                        IF (CRITICAL_BEDROCK_GRADIENT > 0.0) THEN
                            FAILTERM=WEATHER_MULT*ROCK_MASS_WASTING(D8_GRADIENT(I,J))
                        ELSE
                            FAILTERM=0.0
                        ENDIF
                         !     *********************************************************************
                         !      Weathering rate is equal to a flat surface value possibly enhanced
                         !         by slope gradient (if critical_bedrock_gradient>0) and by surface convexity
                         !         (if weather_divergence>0)
                         !     *********************************************************************
                        IF (USE_SOLAR_EROSION) THEN
                            REGOLITH(I,J)=-RAD_CONST*(RAD_ERODE(I,J)+FAILTERM)
                            TOTALRADERODE=TOTALRADERODE+RAD_CONST*RAD_ERODE(I,J)
                        ELSE
                            REGOLITH(I,J)=-MAXRATE*(DIVERGE+FAILTERM)     &
                            /RELATIVE_RESISTANCE(I,J)
                            IF (RESISTANT_SURFACE_LAYER) REGOLITH(I,J)=REGOLITH(I,J)/RELATIVE_RESISTANCE(I,J)
                        ENDIF
                    ENDIF
                        
                ELSE
                    !     *********************************************************************
                    !      For regolith-mantled slopes calculate increase in regolith thickness
                    !         Rate law assumes negative exponential dependency on regolith
                    !         thickness with rock_weathering_rate being the maximum rate at zero
                    !         thickness.  Other functional dependencies could be used, including
                    !         a maximum rate at a finite thickness.
                    !     *********************************************************************
                    WTERM=COSFACTOR*ROCK_WEATHERING_RATE
                    IF (REGOLITH(I,J) <= 0.0) THEN
                        REGOLITH(I,J)=0.0
                        IF ((D8_GRADIENT(I,J) >= 0.0).AND.((ELEVATION(I,J)-SEDIMENT_BASE(I,J)) <= 0.0 )) THEN
                            IS_ROCK_SURFACE(I,J)=.TRUE.
                            IS_SEDIMENT_COVERED(I,J)=.FALSE.
                        ELSE
                            IS_ROCK_SURFACE(I,J)=.FALSE.
                        ENDIF
                    ELSE
                        WTERM = WTERM*DEXP(-WEATHER_DECAY_RATE*REGOLITH(I,J))
                    ENDIF
                    IF (SEEPAGE_WEATHERING) THEN
                        IF (SEEPAGE_AVERAGING) THEN
                            XSEEP=SEEPAGE_WEATHERING_SCALING*FILTERED_GROUNDWATER_FLUX(I,J)**SEEPAGE_WEATHERING_EXPONENT
                        ELSE
                            XSEEP=SEEPAGE_WEATHERING_SCALING*GROUNDWATER_FLUX(I,J)**SEEPAGE_WEATHERING_EXPONENT
                        ENDIF
                        ! EXPERIMENTAL
                        XSEEP=XSEEP*DEXP(-WEATHER_DECAY_RATE*REGOLITH(I,J))
                        ! EXPERIMENTAL
                        IF (XSEEP < 0.0) XSEEP=0.0
                        XAVGSEEP=XAVGSEEP+XSEEP
                        AVGWTERM=AVGWTERM+WTERM
                        XTOTAL=XTOTAL+1.0
                        IF (WTERM > MAXWTERM) MAXWTERM=WTERM
                        IF (XSEEP > MAXSEEP) MAXSEEP=XSEEP
                        WTERM=WTERM+XSEEP
                    ENDIF
                    IF (TWO_TERM_WEATHERING) THEN
                        WTERM = WTERM - WEATHERING_TERM_2*DEXP(-WEATHERING_DECAY_2*REGOLITH(I,J))
                    ENDIF
                    IF (RESISTANT_SURFACE_LAYER) WTERM=WTERM/RELATIVE_RESISTANCE(I,J)
                    REGOLITH(I,J)=REGOLITH(I,J) +  &
                    TIME_INCREMENT * WTERM/RELATIVE_RESISTANCE(I,J)
                ENDIF
                IF (IS_ROCK_SURFACE(I,J)) FRACTBEDROCK=FRACTBEDROCK+1.0
                IF (REGOLITH(I,J) > MAXREGOLITH) MAXREGOLITH=REGOLITH(I,J)
                IF (REGOLITH(I,J) < MINREGOLITH) MINREGOLITH=REGOLITH(I,J)
                AVGREGOLITH=AVGREGOLITH+REGOLITH(I,J)
            ENDDO M110
        ENDDO L110
        IF (USE_SOLAR_EROSION) THEN
            DO J=1,MY
                DO I=1,MX
                    CALL FIND_DEPRESSION(I,J,IS_DEPRESSION)
                    IF (IS_DEPRESSION) THEN
                        IS_ROCK_SURFACE(I,J)=.FALSE.
                    ENDIF
                ENDDO
            ENDDO
        ENDIF
        FRACTBEDROCK=FRACTBEDROCK/(MX*MY)
        AVGREGOLITH=AVGREGOLITH/(MX*MY)
        IF (WRITEDETAIL) WRITE(*,782) FRACTBEDROCK,MINREGOLITH,AVGREGOLITH,MAXREGOLITH,COSFACTOR
        782   FORMAT(' FRACTROCK=',G12.5,' MINREG=',G12.5,' AVGREG=',G12.5,/, &
        ' MAXREG=',G12.5,' COSFACTOR=',G12.5)
        IF (XTOTAL > 0.0) THEN
            XAVGSEEP=XAVGSEEP/XTOTAL
            AVGWTERM=AVGWTERM/XTOTAL
        ENDIF
        IF (WRITEDETAIL) THEN
            IF (SEEPAGE_WEATHERING) THEN
                WRITE(*,298) AVGWTERM,MAXWTERM,XAVGSEEP,MAXSEEP
                298 FORMAT ('AVG & MAX WTERM=',2G12.5,'AVG & MAX XSEEP=',2G12.5)
            ENDIF
            IF (USE_SOLAR_EROSION) THEN
                WRITE(*,844) TOTALRADERODE
                844 FORMAT(' TOTALRADERODE=',G12.5)
            ENDIF
        ENDIF
        RETURN
    END
    !       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !       ********************************************************************
    !        Find depressions
    !    ACCESSES:  IS_X_PERIODIC, IS_Y_PERIODIC, ELEVATION, MX, MY
    !       ********************************************************************
    SUBROUTINE FIND_DEPRESSION(I,J,IS_DEPRESSION)
        USE ERODE_GLOBALS
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: I,J
        INTEGER :: IL,IH,JL,JH,II,JJ,IC,JC
        REAL (8) :: ECOMP
        LOGICAL, INTENT(OUT) :: IS_DEPRESSION
        ECOMP=ELEVATION(I,J)
        IL=I-1
        IH=I+1
        JL=J-1
        JH=J+1
        IS_DEPRESSION=.TRUE.
        DO JJ=JL,JH
                JC=JJ
                IF (IS_Y_PERIODIC) THEN
                    IF (JC < 1) JC=MY
                    IF (JC > MY) JC=1
                ELSE
                    IF (JC < 1) CYCLE
                    IF (JC > MY) CYCLE
                ENDIF
            DO II=IL,IH
                IC=II
                IF (IS_X_PERIODIC) THEN
                    IF (IC < 1) IC=MX
                    IF (IC > MX) IC=1
                ELSE
                    IF (IC < 1) CYCLE
                    IF (IC > MY) CYCLE
                ENDIF
                IF ((IC == I).AND.(JC == J)) CYCLE
                IF (ELEVATION(IC,JC) <= ECOMP) IS_DEPRESSION=.FALSE.
            ENDDO
        ENDDO
        RETURN
    END

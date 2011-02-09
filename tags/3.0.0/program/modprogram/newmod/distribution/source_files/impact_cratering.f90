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
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE DO_IMPACT_CRATERING()
        USE ERODE_GLOBALS
        USE CRATER_GLOBALS
        IMPLICIT NONE
        !     **********************************************************************
        !       This routine creates an impact crater in a random location with a
        !         size determined by a production function and geometry governed by
        !         a number of input parameters.  the model is documeted in:
        !          forsberg-taylor, n. k., howard, a. d., and craddock, r. a., 2004,
        !            j. geophys. res., 109, e05002, doi:10.1029/2004je002242
        !         and
        !          howard, a. d, 2007, geomorphology, 91, 332-363
        !  CALLS: GET_CRATER_SIZE, FIND_IMPACT_SITE, FIND_REFERENCE_ELEVATION, CREATE_CRATER
        !     **********************************************************************
        COUNT2=0
        COUNT10=0
        !     **********************************************************************
        !
        !     ******** Generate a random crater size from the assumed population
        !              distribution
        !
        !     **********************************************************************
        CALL GET_CRATER_SIZE()
        !     **********************************************************************
        !
        !     ******** Where in the total playing field did it hit?
        !
        !     **********************************************************************
        CALL FIND_IMPACT_SITE()
        !     **********************************************************************
        !
        !     ******** Determine the average elevation of the impactor footprint
        !
        !     **********************************************************************
        CALL FIND_REFERENCE_ELEVATION()
        !     **********************************************************************
        !
        !     ******** Do the excavation and deposition
        !
        !     **********************************************************************
        CALL CREATE_CRATER()
        RETURN
    END
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    SUBROUTINE GET_CRATER_SIZE()
        USE ERODE_GLOBALS
        USE CRATER_GLOBALS
        IMPLICIT NONE
        REAL (4) :: TEMP
        REAL (8) :: TEMP1,DIAMRATIO
        REAL (8) :: RRAND
        EXTERNAL RRAND
        !     **********************************************************************
        !
        !     ******* Determine a random impactor diameter from a negative
        !             power function cumulative diameter size.  Also determine
        !             parameters governing excavation and deposition shape
        !             functions and the maximum range of the impactor footprint
        !  MODIFIES: DIAMETER, RADIUS, CRATER_DEPTH, RIM_HEIGHT,INTERIOR_SHAPE_EXPONENT
        !            EXTERIOR_SHAPE_EXPONENT, DIAMRATIO, MAXRANGE
        !
        !     **********************************************************************
        L110: DO
            TEMP=RRAND()
            IF (TEMP == 0.0) CYCLE L110
            DIAMETER=(CCON/(TEMP-BCON))**ALINV
            TEMP=RRAND()
            TEMP1=EXP(FREQUENCY_CUTOFF_SCALING*(SMALLEST_POSSIBLE_CRATER-DIAMETER))
            IF (TEMP < TEMP1) CYCLE L110
            EXIT L110
        ENDDO L110
        WRITE(*,7000) DIAMETER
        WRITE(OUTHIST,7000) DIAMETER
        7000  FORMAT(' CRATERING, DIAM=',G12.5)
        RADIUS=DIAMETER/2.0
        DIAMRATIO=DIAMETER/TRANSITION_DIAMETER
        IF (DIAMRATIO >= 1.0) THEN
            CRATER_DEPTH=LARGE_CRATER_DEPTH_SCALE*DIAMETER**LARGE_CRATER_DEPTH_EXPONENT
            RIM_HEIGHT=LARGE_CRATER_RIM_SCALE*DIAMETER**LARGE_CRATER_RIM_EXPONENT
            INTERIOR_SHAPE_EXPONENT=LARGE_CRATER_SHAPE_SCALE*DIAMETER**LARGE_CRATER_SHAPE_EXPONENT
        ELSE
            CRATER_DEPTH=SMALL_CRATER_DEPTH_SCALE*DIAMETER**SMALL_CRATER_DEPTH_EXPONENT
            RIM_HEIGHT=SMALL_CRATER_RIM_SCALE*DIAMETER**SMALL_CRATER_RIM_EXPONENT
            INTERIOR_SHAPE_EXPONENT=SMALL_CRATER_SHAPE_SCALE*DIAMETER**SMALL_CRATER_SHAPE_EXPONENT
        ENDIF
        EXTERIOR_SHAPE_EXPONENT=2.0-(RIM_HEIGHT/((RIM_HEIGHT-CRATER_DEPTH)/2.0+CRATER_DEPTH/(INTERIOR_SHAPE_EXPONENT+2)))
        MAXRANGE=DMIN1(MX*CELL_SIZE,RADIUS*(1.0/0.01)**(1.0/(EXTERIOR_SHAPE_EXPONENT-1.0)))
        !WRITE (*,100) MAXRANGE, CRATER_DEPTH, RIM_HEIGHT,INTERIOR_SHAPE_EXPONENT, EXTERIOR_SHAPE_EXPONENT
        100   FORMAT(' MAXRANGE=',G12.5,' CRATER_DEPTH=',G12.5,' RIM_HEIGHT=',G12.5, &
        ' INTERIOR_SHAPE_EXPONENT=',G12.5,' EXTERIOR_SHAPE_EXPONENT=',G12.5)
        RETURN
    END
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    SUBROUTINE FIND_MODIFICATION_RANGE()
        USE ERODE_GLOBALS
        USE CRATER_GLOBALS
    !  MODIFIES: XMIN, XMAX, YMIN, YMAX
        IMPLICIT NONE
        IF (IS_X_PERIODIC) THEN
            XMIN=0.0
            XMAX=MX*CELL_SIZE
        ELSE
            XMIN=-MAXRANGE
            XMAX=MX*CELL_SIZE+MAXRANGE
        ENDIF
        IF (IS_Y_PERIODIC) THEN
            YMIN=0.0
            YMAX=MY*CELL_SIZE
        ELSE
            YMIN=-MAXRANGE
            YMAX=MY*CELL_SIZE+MAXRANGE
        ENDIF
        RETURN
    END
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    SUBROUTINE FIND_IMPACT_SITE()
        USE ERODE_GLOBALS
        USE CRATER_GLOBALS
        IMPLICIT NONE
        REAL (4) :: XTEMP,YTEMP
        INTEGER :: ICNTR,JCNTR
        REAL (8) :: RRAND
        EXTERNAL RRAND
        !     **********************************************************************
        !
        !     ****** Find a random location within the virtual target domain defined
        !            by xmin,xmax,ymin,ymax and determine the matrix indices of the
        !            footprint
        !  MODIFIES: XCENTER, YCENTER, ICNTR, JCNTR, IMIN, IMAX, JMIN, JMAX
        !            IMINR, IMAXR, JMINR, JMAXR
        !     **********************************************************************
        XTEMP=RRAND()
        YTEMP=RRAND()
        XCENTER=XMIN+(XMAX-XMIN)*XTEMP
        YCENTER=YMIN+(YMAX-YMIN)*YTEMP
        ICNTR=XCENTER/CELL_SIZE+1
        JCNTR=YCENTER/CELL_SIZE+1
        IMIN=MAX0(1,INT(1+(XCENTER-MAXRANGE)/CELL_SIZE))
        IMAX=MIN0(MX,INT(1+(XCENTER+MAXRANGE)/CELL_SIZE))
        JMIN=MAX0(1,INT(1+(YCENTER-MAXRANGE)/CELL_SIZE))
        JMAX=MIN0(MY,INT(1+(YCENTER+MAXRANGE)/CELL_SIZE))
        !     **********************************************************************
        !
        !     ***** If periodic left and right boundaries are used, iminr and imaxr
        !           are the larger domain that allow boundary overlap
        !
        !     **********************************************************************
        IF (IS_X_PERIODIC) THEN
            IF (DO_EJECTA_WRAPAROUND) THEN
                IMINR=MAX0(-MX+1,INT(1+(XCENTER-MAXRANGE)/CELL_SIZE))
                IMAXR=MIN0(2*MX,INT(1+(XCENTER+MAXRANGE)/CELL_SIZE))
            ELSE
                IMINR=MAX0(ICNTR-MX/2+1,INT(1+(XCENTER-MAXRANGE)/CELL_SIZE))
                IMAXR=MIN0(ICNTR+MX/2,INT(1+(XCENTER+MAXRANGE)/CELL_SIZE))
            ENDIF
        ELSE
            IMINR=IMIN
            IMAXR=IMAX
        ENDIF
        IF (IS_Y_PERIODIC) THEN
            IF (DO_EJECTA_WRAPAROUND) THEN
                JMINR=MAX0(-MY+1,INT(1+(YCENTER-MAXRANGE)/CELL_SIZE))
                JMAXR=MIN0(2*MY,INT(1+(YCENTER+MAXRANGE)/CELL_SIZE))
            ELSE
                JMINR=MAX0(JCNTR-MY/2+1,INT(1+(YCENTER-MAXRANGE)/CELL_SIZE))
                JMAXR=MIN0(JCNTR+MY/2,INT(1+(YCENTER+MAXRANGE)/CELL_SIZE))
            ENDIF
        ELSE
            JMINR=JMIN
            JMAXR=JMAX
        ENDIF
        IF     &
        ((IMINR == -128).OR.(IMAXR == 256).OR.(JMINR == -128).OR. &
        (JMAXR == 256)) THEN
        ENDIF
        RETURN
    END
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    SUBROUTINE CREATE_CRATER()
        USE ERODE_GLOBALS
        USE CRATER_GLOBALS
        IMPLICIT NONE
        INTEGER :: I,J,ITEMP,JTEMP
        REAL (8) :: XCOMP,YCOMP,DCOMP,TEMP1,TEMP2,DEPOSITTHICK
        REAL (8) :: EUSE,LOGNORMAL_RANDOM_DEVIATE,TOTALVOL,BIASADD,TOTPOINTS
        REAL (8) :: LOCALGRAV,LOCALWORK
        REAL (8) :: NETCHANGE,NETPOINTS,DELELEV
        LOGICAL (1) :: NOTDONE(IMMX,JMMX),NEWSEDIMENT_BASE(IMMX,JMMX)
        !     **********************************************************************
        !
        !      ******  This subroutine does the actual work of excavating and
        !              depositing.  Small craters off of the target area are
        !              bypassed
        !  MODIFIES: ERODE_SLOPE, CFNE, REGOLITH, IS_ROCK_SURFACE, IS_SEDIMENT_COVERED
        !            ELEVATION, CUMULATIVE_CRATERING_CHANGE, SEDIMENT_BASE
        !            CRATERWORK, LOCALWORK, CRATERGRAV, LOCALGRAV
        !  CALLS: LOGNORMAL_RANDOM_DEVIATE
        !
        !     **********************************************************************
        IF ((IMAX >= IMIN).AND.(JMAX >= JMIN)) THEN
            TOTALVOL=0.0
            TOTPOINTS=0.0
            NETCHANGE=0.0
            NETPOINTS=0.0
            LOCALGRAV=0.0
            LOCALWORK=0.0
            IF (FLUVIAL_AND_SLOPE_MODELING) THEN
                WRITE(OUTCRATER,120) DIAMETER &
                ,XCENTER,YCENTER,ITERATION,PRESENT_TIME
                120   FORMAT(3(G12.5,' '),I5,' ',G12.5)
            ELSE
                WRITE(OUTCRATER,1120) DIAMETER &
                ,XCENTER,YCENTER,ITERATION
                1120   FORMAT(3(G12.5,' '),I5)
            ENDIF
            !     **********************************************************************
            !
            !      ****** If we're here we model the impact.  Update counters.
            !
            !     **********************************************************************
            ITOTHITS=ITOTHITS+1
            IF (DIAMETER > 2.0) COUNT2=COUNT2+1
            IF (DIAMETER > 10.0) COUNT10=COUNT10+1
            NOISECALL=.TRUE.
            L2100: DO  JTEMP=JMINR,JMAXR
                M2100: DO ITEMP=IMINR,IMAXR
                    I=ITEMP
                    J=JTEMP
                    !     **********************************************************************
                    !
                    !     ******   If the lateral boundaries are periodic, then allow the
                    !              boundary overlap
                    !
                    !     **********************************************************************
                    IF (IS_X_PERIODIC) THEN
                        IF (ITEMP < 1) I=ITEMP+MX
                        IF (ITEMP > MX) I=ITEMP-MX
                    ENDIF
                    IF (IS_Y_PERIODIC) THEN
                        IF (JTEMP < 1) J=JTEMP+MY
                        IF (JTEMP > MY) J=JTEMP-MY
                    ENDIF
                    ERODE_SLOPE(I,J)=0.0
                    NOTDONE(I,J)=.TRUE.
                    NEWSEDIMENT_BASE(I,J)=.FALSE.
                ENDDO M2100
            ENDDO L2100
            !     **********************************************************************
            !
            !     ****** Cycle across the footprint area on the target matrix
            !
            !     **********************************************************************
            L100: DO JTEMP=JMINR,JMAXR
                M100: DO ITEMP=IMINR,IMAXR
                    I=ITEMP
                    J=JTEMP
                    !     **********************************************************************
                    !
                    !     ******   If the lateral boundaries are periodic, then allow the
                    !              boundary overlap
                    !
                    !     **********************************************************************
                    IF (IS_X_PERIODIC) THEN
                        IF (ITEMP < 1) I=ITEMP+MX
                        IF (ITEMP > MX) I=ITEMP-MX
                    ENDIF
                    IF (IS_Y_PERIODIC) THEN
                        IF (JTEMP < 1) J=JTEMP+MY
                        IF (JTEMP > MY) J=JTEMP-MY
                    ENDIF
                    XCOMP=(ITEMP-1)*CELL_SIZE
                    YCOMP=(JTEMP-1)*CELL_SIZE
                    DCOMP=DSQRT((XCOMP-XCENTER)**2+(YCOMP-YCENTER)**2)
                    !     **********************************************************************
                    !
                    !     *****   dcomp is the radial  distance of the present location on
                    !             the footprint from the crater center.  We model the interior
                    !             and exterior of the crater separately
                    !
                    !     **********************************************************************
                    IF (DCOMP <= RADIUS) THEN
                        !     **********************************************************************
                        !
                        !     *****   Model the crater interior
                        !
                        !     **********************************************************************
                        !                     depositthick=(rim_height-crater_depth)+crater_depth*(dcomp/radius)**interior_shape_exponent
                        DEPOSITTHICK=CRATER_DEPTH*(DCOMP/RADIUS)**INTERIOR_SHAPE_EXPONENT
                        !     **********************************************************************
                        !
                        !     *****   Random variability of crater elevation is presently allowed
                        !             only in areas of net deposition.  Perhaps it also ought to
                        !             be allowed on the lower crater rim
                        !
                        !     **********************************************************************
                        IF (RANDOM_EJECTA_THICKNESS) THEN
                            DEPOSITTHICK=DEPOSITTHICK &
                            *LOGNORMAL_RANDOM_DEVIATE(EJECTA_THICKNESS_VARIABILITY)
                        ENDIF
                        DEPOSITTHICK=DEPOSITTHICK+(RIM_HEIGHT-CRATER_DEPTH)
                        !     **********************************************************************
                        !
                        !     *****   The amount of inheritance of the original topography increases
                        !             from zero at the crater center (where the average elevation of
                        !             the footprint is used) to inheritance_parameter at the rim
                        !
                        !     **********************************************************************
                        EUSE=INHERITANCE_PARAMETER*(DCOMP/RADIUS)**2
                        ERODE_SLOPE(I,J)=ERODE_SLOPE(I,J)+  &
                        (1.0-EUSE)*(ELEVBASE-ELEVATION(I,J))+ &
                        DEPOSITTHICK
                        IF (FLUVIAL_AND_SLOPE_MODELING) THEN
                            IS_SEDIMENT_COVERED(I,J)=.FALSE.
                            NEWSEDIMENT_BASE(I,J)=.TRUE.
                            CFNE(I,J)=0.0
                            IF (IS_REGOLITH_CRATER )  THEN
                                IS_ROCK_SURFACE(I,J)=.FALSE.
                                REGOLITH(I,J)=INITIAL_REGOLITH_THICKNESS
                            ELSE
                                IS_ROCK_SURFACE(I,J)=.TRUE.
                                REGOLITH(I,J)=-ROCK_WEATHERING_RATE
                            ENDIF
                        ENDIF
                    ELSE
                        !     **********************************************************************
                        !
                        !     ******  Model the crater exterior.  The fractional inheritance of the
                        !             original topography goes from inheritance_parameter at the rim to 100%
                        !             at the far margins of deposition
                        !
                        !     **********************************************************************
                        TEMP1=(RADIUS/DCOMP)**EXTERIOR_SHAPE_EXPONENT
                        TEMP2=DMIN1((1.0-INHERITANCE_PARAMETER),TEMP1)
                        DEPOSITTHICK=RIM_HEIGHT*TEMP1
                        !     **********************************************************************
                        !
                        !     ******  Determine the random variability of depositional amount.
                        !
                        !     **********************************************************************
                        IF (RANDOM_EJECTA_THICKNESS) THEN
                            DEPOSITTHICK=DEPOSITTHICK*LOGNORMAL_RANDOM_DEVIATE(EJECTA_THICKNESS_VARIABILITY)
                        ENDIF
                        ERODE_SLOPE(I,J)=ERODE_SLOPE(I,J)+(ELEVBASE-ELEVATION(I,J))*TEMP2  &
                        +DEPOSITTHICK
                        IF (FLUVIAL_AND_SLOPE_MODELING) THEN
                            IF (ABS(ERODE_SLOPE(I,J)) > 0.001) THEN
                                IS_SEDIMENT_COVERED(I,J)=.FALSE.
                                NEWSEDIMENT_BASE(I,J)=.TRUE.
                                IF (IS_REGOLITH_CRATER )  THEN
                                    IS_ROCK_SURFACE(I,J)=.FALSE.
                                    REGOLITH(I,J)=INITIAL_REGOLITH_THICKNESS
                                ELSE
                                    IS_ROCK_SURFACE(I,J)=.TRUE.
                                    REGOLITH(I,J)=-ROCK_WEATHERING_RATE
                                ENDIF
                                CFNE(I,J)=0.0
                            ENDIF
                        ENDIF
                    ENDIF
                ENDDO M100
            ENDDO L100
            L3100: DO JTEMP=JMINR,JMAXR
                M3100: DO ITEMP=IMINR,IMAXR
                    I=ITEMP
                    J=JTEMP
                    !     **********************************************************************
                    !
                    !     ******   If the lateral boundaries are periodic, then allow the
                    !              boundary overlap
                    !
                    !     **********************************************************************
                    IF (IS_X_PERIODIC) THEN
                        IF (ITEMP < 1) I=ITEMP+MX
                        IF (ITEMP > MX) I=ITEMP-MX
                    ENDIF
                    IF (IS_Y_PERIODIC) THEN
                        IF (JTEMP < 1) J=JTEMP+MY
                        IF (JTEMP > MY) J=JTEMP-MY
                    ENDIF
                    IF (NOTDONE(I,J)) THEN
                        TOTALVOL=TOTALVOL+ERODE_SLOPE(I,J)
                        TOTPOINTS=TOTPOINTS+1.0
                        NOTDONE(I,J)=.FALSE.
                    ENDIF
                ENDDO M3100
            ENDDO L3100
            BIASADD=-TOTALVOL/TOTPOINTS
            !            write(*,3110) totalvol,totpoints,biasadd
            3110      FORMAT(' OV=',G12.5,' T=',G12.5,' BA=',G12.5)
            L4100: DO  JTEMP=JMINR,JMAXR
                M4100: DO ITEMP=IMINR,IMAXR
                    I=ITEMP
                    J=JTEMP
                    !     **********************************************************************
                    !
                    !     ******   If the lateral boundaries are periodic, then allow the
                    !              boundary overlap
                    !
                    !     **********************************************************************
                    IF (IS_X_PERIODIC) THEN
                        IF (ITEMP < 1) I=ITEMP+MX
                        IF (ITEMP > MX) I=ITEMP-MX
                    ENDIF
                    IF (IS_Y_PERIODIC) THEN
                        IF (JTEMP < 1) J=JTEMP+MY
                        IF (JTEMP > MY) J=JTEMP-MY
                    ENDIF
                    NOTDONE(I,J)=.TRUE.
                ENDDO M4100
            ENDDO L4100
            L1100: DO  JTEMP=JMINR,JMAXR
                M1100: DO  ITEMP=IMINR,IMAXR
                    I=ITEMP
                    J=JTEMP
                    !     **********************************************************************
                    !
                    !     ******   If the lateral boundaries are periodic, then allow the
                    !              boundary overlap
                    !
                    !     **********************************************************************
                    IF (IS_X_PERIODIC) THEN
                        IF (ITEMP < 1) I=ITEMP+MX
                        IF (ITEMP > MX) I=ITEMP-MX
                    ENDIF
                    IF (IS_Y_PERIODIC) THEN
                        IF (JTEMP < 1) J=JTEMP+MY
                        IF (JTEMP > MY) J=JTEMP-MY
                    ENDIF
                    XCOMP=(ITEMP-1)*CELL_SIZE
                    YCOMP=(JTEMP-1)*CELL_SIZE
                    DCOMP=DSQRT((XCOMP-XCENTER)**2+(YCOMP-YCENTER)**2)
                    IF (NOTDONE(I,J)) THEN
                        DELELEV=+ERODE_SLOPE(I,J)+BIASADD
                        ELEVATION(I,J)=ELEVATION(I,J)+DELELEV
                        CRATERWORK=CRATERWORK+DCOMP*DELELEV
                        LOCALWORK=LOCALWORK+DCOMP*DELELEV
                        CRATERGRAV=CRATERGRAV-0.5*ABS(DELELEV)
                        LOCALGRAV=LOCALGRAV-0.5*ABS(DELELEV)
                        CUMULATIVE_CRATERING_CHANGE(I,J)=CUMULATIVE_CRATERING_CHANGE(I,J) &
                        + DELELEV
                        IF (FLUVIAL_AND_SLOPE_MODELING) THEN
                            IF (NEWSEDIMENT_BASE(I,J)) SEDIMENT_BASE(I,J)=ELEVATION(I,J)
                        ENDIF
                        NETCHANGE=NETCHANGE+ERODE_SLOPE(I,J)+BIASADD
                        NETPOINTS=NETPOINTS+1.0
                        NOTDONE(I,J)=.FALSE.
                    ENDIF
                ENDDO M1100
            ENDDO L1100
            NETCHANGE=NETCHANGE/NETPOINTS
        ELSE
            !WRITE(*,130)
            130   FORMAT(' MISS')
        ENDIF
        !WRITE(*,7552) LOCALWORK,LOCALGRAV
        7552  FORMAT(' LOCALWORK=',G12.5,' LOCALGRAV=',G12.5)
        RETURN
    END
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE FIND_REFERENCE_ELEVATION()
        USE ERODE_GLOBALS
        USE CRATER_GLOBALS
        IMPLICIT NONE
        REAL (8) :: TEMP,WEIGHT,XCOMP,YCOMP,DCOMP
        INTEGER :: I,J,ITEMP,JTEMP
        !     **********************************************************************
        !
        !     ***** Finds the average elevation of the footprint. The only trick
        !           is that the locations lying within the new crater interior are
        !           given greatest weighting and that points outside the crater rim
        !           are given a weighting in determining the average elevation in
        !           proportion to the amount of deposition that will occur
        !  MODIFIES: ELEVBASE
        !
        !     **********************************************************************
        ELEVBASE=0.0
        WEIGHT=0.0
        IF ((IMAX >= IMIN).AND.(JMAX >= JMIN)) THEN
            L100: DO  JTEMP=JMINR,JMAXR
                M100: DO ITEMP=IMINR,IMAXR
                    I=ITEMP
                    J=JTEMP
                    !     **********************************************************************
                    !
                    !     ******   If the lateral boundaries are periodic, then allow the
                    !              boundary overlap
                    !
                    !     **********************************************************************
                    IF (IS_X_PERIODIC) THEN
                        IF (ITEMP < 1) I=ITEMP+MX
                        IF (ITEMP > MX) I=ITEMP-MX
                    ENDIF
                    IF (IS_Y_PERIODIC) THEN
                        IF (JTEMP < 1) J=JTEMP+MY
                        IF (JTEMP > MY) J=JTEMP-MY
                    ENDIF
                    XCOMP=(I-1)*CELL_SIZE
                    YCOMP=(J-1)*CELL_SIZE
                    DCOMP=DSQRT((XCOMP-XCENTER)**2+(YCOMP-YCENTER)**2)
                    IF (DCOMP <= RADIUS) THEN
                        ELEVBASE=ELEVBASE+ELEVATION(I,J)
                        WEIGHT=WEIGHT+1.0
                    ELSE
                        TEMP=(RADIUS/DCOMP)**EXTERIOR_SHAPE_EXPONENT
                        ELEVBASE=ELEVBASE+ELEVATION(I,J)*TEMP
                        WEIGHT=WEIGHT+TEMP
                    ENDIF
                ENDDO M100
            ENDDO L100
            ELEVBASE=ELEVBASE/WEIGHT
        ENDIF
        RETURN
    END

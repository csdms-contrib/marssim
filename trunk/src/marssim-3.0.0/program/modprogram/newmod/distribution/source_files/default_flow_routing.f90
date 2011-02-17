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
    REAL (8) FUNCTION DISCHARGE_FROM_CELL(I,J)
        USE ERODE_GLOBALS
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: I,J
        REAL (8) :: TEMP1,TEMP2,TEMP3,TEMP4,RECHARGERATIO
        RECHARGERATIO=GROUNDWATER_RECHARGE_RATE*1.0
        !     **********************************************************************
        !      Determines the local specific runoff yield.  Can depend upon the local slope
        !         divergence and/or whether there is a bistable erosion rate
        !     **********************************************************************
        IF (DIVERGENCE_DEPENDENT_RUNOFF) THEN
            TEMP1 = RUNOFF_CONSTANT_3*(DIVERGENCE(I,J)-MEAN_CONVERGENCE_RUNOFF)
            TEMP1 = RUNOFF_CONSTANT_2+RUNOFF_CONSTANT_1*TEMP1/DSQRT(1.0+TEMP1**2)
        ELSE
            TEMP1 =1.0
        ENDIF
        TEMP1=DISCHARGE_SCALE_FACTOR*TEMP1
        IF (ACCELERATED_EROSION(I,J)) THEN
            DISCHARGE_FROM_CELL=BISTABLE_RUNOFF_FACTOR*TEMP1
        ELSE
            DISCHARGE_FROM_CELL=TEMP1
        ENDIF
        100   FORMAT(' DISCHARGE_FROM_CELL=',G12.5)
        !     **********************************************************************
        !     *****   New groundwater code - add in groundwater flow component
        !             If groundwater_flow_fraction=1 all surface flow from groundwater plus
        !             recharge (rainfall) on saturated ground (fixed_head true)
        !     **********************************************************************
        IF (MODEL_GROUNDWATER) THEN
            TEMP2=DMAX1(0.0D0,GROUNDWATER_FLUX(I,J))
            IF (FIXED_HEAD(I,J)) THEN
                TEMP2=TEMP2+RECHARGERATIO
            ELSE
                TEMP2=0.0
            ENDIF
            TEMP3=DISCHARGE_FROM_CELL*(1.0-GROUNDWATER_FLOW_FRACTION)
            TEMP4=TEMP2*GROUNDWATER_FLOW_FRACTION
            DISCHARGE_FROM_CELL=TEMP3+TEMP4
            AVGRUNOFF=AVGRUNOFF+TEMP3
            AVGSEEP=AVGSEEP+TEMP4
            IF (TEMP2 > 0.0) THEN
                AVGSEEPFRACT=AVGSEEPFRACT+TEMP4/DISCHARGE_FROM_CELL
                NFRACTAVG=NFRACTAVG+1.0
            ENDIF
            NSEEPAVG=NSEEPAVG+1.0
        ENDIF
        DISCHARGE_FROM_CELL=DISCHARGE_FROM_CELL*D2X
        RETURN
    END
    !       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE IS_IT_SUBMERGED()
        !      ************************************************************
        !        Determines if location i,j is submerged (submerged(i,j)=.true.)
        !           or not (submerged(i,j)=.false.)
        !   MODIFIES: SUBMERGED
        !      ***********************************************************
        USE ERODE_GLOBALS
        IMPLICIT NONE
        INTEGER :: I,J,NNN
        SUBMERGED=.FALSE.
        DO  J=1,MYY
            DO  I=1,MX
                IF (USE_AN_OCEAN) THEN
                    IF (ELEVATION(I,J) <= OCEAN_ELEVATION) SUBMERGED(I,J)=.TRUE.
                ENDIF
                IF(.NOT.COMPLETE_RUNOFF) CYCLE
                NNN=BASIN_NUMBER(I,J)
                IF (MODEL_LAKE_EVAPORATION) THEN
                    IF ((COMPLETE_RUNOFF).AND.(ENCLOSED(NNN)).AND.  &
                    (ELEVATION(I,J) < LAKE_SURFACE_ELEVATION(NNN))) THEN
                        SUBMERGED(I,J)=.TRUE.
                    ENDIF
                ELSE
                    IF ((COMPLETE_RUNOFF).AND.(ENCLOSED(NNN)).AND.  &
                    (ELEVATION(I,J) < LAKE_OUTLET_ELEVATION(NNN))) THEN
                        SUBMERGED(I,J)=.TRUE.
                    ENDIF
                ENDIF
                IF (DO_FLOW_BOUNDARIES) THEN
                ELSE
                    IF (COMPLETE_RUNOFF.AND.(.NOT.IS_Y_PERIODIC)) THEN
                        IF (ELEVATION(I,J) <= ELEVATION(I,MY)) SUBMERGED(I,J)=.TRUE.
                    ENDIF
                ENDIF
            ENDDO
        ENDDO
        RETURN
    END
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE DRAINAGE_BASIN_AREA_FLOW()
        USE ERODE_GLOBALS, CHANNEL_EROSION_RATE => CFNE
        USE AREA_GLOBALS
        USE LAKE_GLOBALS
        IMPLICIT NONE
        !      *******************************************************************
        !       Calculates drainage areas for each location via flow routing,
        !       including overflow from enclosed depressions
        !
        !       Also determines total sediment yield from channel erosion
        !
        !       Local vectors are defined in area_globals and lake_globals:
        !           (i=integer, l=logical, else real)
        !
        !           tobeadded(k)   - the cumulative drainage area of  enclosed drainage
        !                            basins upstream
        !       l   needtodo(k)    - working variable which indicates which enclosed
        !                            basins need to be considered for grouping with
        !                            neighboring basins
        !       l   outer(k) - true if the exit point from an enclosed basin is
        !                      lower than any exterior contiguous point outside the basin
        !           local_basin_discharge(k) - the total discharge from the enclosed basin
        !           qtoadd(k) - the cumulative contributed discharge from enclosed
        !                       drainage basins upstream
        !   MODIFIES: I_OUTFLOW, J_OUTFLOW, ENCLOSED, SEDIMENT_YIELD, FLOW_DIRECTION
        !             DISCHARGE, DRAINAGE_AREA, BASIN_NUMBER,LAKE_AREA,IS_EXIT_BASIN
        !             ROUTED_DRAINAGE_AREA, DOWNSTREAM_BASIN
        !             LAKE_OUTLET_ELEVATION, LAKE_SURFACE_ELEVATION,LAKE_AREA
        !             ERODETOADD,ERGRAVTOADD 
        !   CALLS:  DRAWLINE, DISCHARGE_FROM_CELL, BASIN_REPORT 
        !
        !      *******************************************************************
        INTEGER :: I,J,LOCAL_DIRECTION,I_OLD,J_OLD,I_NEW,J_NEW,K,L,KK,KKFROM
        INTEGER :: IFROM,JFROM,NUMBBASIN,NNN,KKK,III,JJJ
        INTEGER :: IOUTLET,JOUTLET,NCONTIG
        INTEGER :: II,JJ,ILL,JLL,IUU,JUU,ILOW,IHIGH,JLOW,JHIGH,KAROUND
        INTEGER :: KKKK,KEND,KKCOUNT,KL
        INTEGER :: KKSQUARE
        INTEGER :: NBOTTOMBASIN,NBOTTOMADD,CURRENT_BASIN 
        REAL (8) :: ADDAREA
        REAL (8) :: AAMAX,AAMIN,AAAVG,AASD
        LOGICAL :: CONTIGUOUS,ISFOUND,TOSEDCOVER
        LOGICAL :: TOLAKE,NEWBASIN,REDO,XCONTINUE,YCONTINUE
        REAL (8) :: XCOMP,SEDADD,RIVERADD
        REAL (8) :: WSUM,WNUM,DISCHARGE_FROM_CELL,QADD,ADDQ
        REAL (8) :: TBEDCONT,TSEDCONT,ELEVLOWEST
        REAL (8) :: THEAVG,THEMAX,THEMIN
        REAL (8) :: GLOBAL_MAXIMUM_DISCHARGE
        AAMAX=-1.0E+25
        ELEVLOWEST=1.0E+25
        AAMIN=-AAMAX
        AAAVG=0.0
        AASD=0.0
        NBOTTOMBASIN=0
        NBOTTOMADD=0
        GLOBAL_MAXIMUM_DISCHARGE=-1.0E+25
        AVGRUNOFF=0.0
        AVGSEEP=0.0
        AVGSEEPFRACT=0.0
        NSEEPAVG=0.0
        NFRACTAVG=0.0
        !      *******************************************************************
        !       no_lakes_exist is index of presence of any enclosed basins
        !      *******************************************************************
        NO_LAKES_EXIST=.TRUE.
        ABORT_SIMULATION=.FALSE.
        !      *******************************************************************
        !       set initial values
        !            current_basin is the current basin number
        !      *******************************************************************
        CURRENT_BASIN=0
        TOBEADDED=0
        DOWNSTREAM_BASIN=0
        I_OUTFLOW=0
        J_OUTFLOW=0
        BASIN_DRAINAGE_AREA=0
        LOCAL_BASIN_DISCHARGE=0.0
        QTOADD=0.0
        ENCLOSED =.FALSE.
        IS_EXIT_BASIN=.FALSE.
        TBEDCONT=0.0
        TSEDCONT=0.0
        !      *********************************************************************
        !        Initialize matrices
        !         "sediment_yield" will hold the total sediment supplied by bed erosion of bedrock
        !         channels upstream from the given location, except if a depression
        !         intervenes in the flow path, which zeros the sediment yield.  Sediment_yield is
        !         defifined for both alluvial and bedrock areas
        !         "sediment_yield" is used for downstream sediment routing modeling of
        !         erosion/deposition
        !      *********************************************************************
        BASIN_NUMBER=0
        SEDIMENT_YIELD=0.0
        DISCHARGE=0.0
        DRAINAGE_AREA=0.0
        IF (IS_Y_PERIODIC) THEN
            DO  J=1,MY
                DO  I=1,MX
                    IF (ELEVATION(I,J) < ELEVLOWEST) THEN
                        ELEVLOWEST=ELEVATION(I,J)
                        ILOWEST=I
                        JLOWEST=J
                    ENDIF
                ENDDO
            ENDDO
        ELSE
            ILOWEST=0
            JLOWEST=0
        ENDIF
        IF (IS_Y_PERIODIC) THEN
        ELSE
            DO  I=1,MX
                FLOW_DIRECTION(I,MY)=1
            ENDDO
        ENDIF
        !      *******************************************************************
        !       Main loop for determination of drainage areas
        !      *******************************************************************
        L110: DO  J=1,MYY
            1105  FORMAT('J=',I5)
            M110: DO  I=1,MX
                1106  FORMAT('I=',I5)
                !      *******************************************************************
                !       tolake is index of whether basin is enclosed
                !       newbasin is index of whether a new basin has been detected
                !         (either a new exit or a new enclosed depression)
                !       riveradd is contribution (if any) from rivers entering from outside
                !         the matrix
                !      *******************************************************************
                TOLAKE=.FALSE.
                NEWBASIN=.FALSE.
                RIVERADD=0.0
                I_OLD=I
                J_OLD=J
                !      *******************************************************************
                !       Increment drainage area for locally contributed drainage and for
                !         incoming streams
                !      *******************************************************************
                IF (HAVE_INFLUENT_RIVERS.AND.IS_INFLUENT_RIVER_LOCATION(I,J)) THEN
                    DO KL=1,NUMBER_OF_INFLUENT_RIVERS
                        IF((I == I_INFLUENT_RIVER(KL)).AND.(J == J_INFLUENT_RIVER(KL)))THEN
                            RIVERADD=INFLUENT_RIVER_DISCHARGE(KL)
                        ENDIF
                    ENDDO
                ENDIF
                DRAINAGE_AREA(I,J)=DRAINAGE_AREA(I,J)+D2X+RIVERADD
                QADD=DISCHARGE_FROM_CELL(I,J)+RIVERADD
                DISCHARGE(I,J)=DISCHARGE(I,J)+QADD
                L120: DO
                    !      *******************************************************************
                    !       Determine local drainage direction and check for bad value
                    !       Implicit loop for traversing downstream
                    !      *******************************************************************
                    LOCAL_DIRECTION = FLOW_DIRECTION(I_OLD,J_OLD)
                    IF (DO_FLOW_BOUNDARIES) THEN
                        IF (LOCAL_DIRECTION == 1) THEN
                            IF (BASIN_NUMBER(I_OLD,J_OLD) == 0) THEN
                                NEWBASIN=.TRUE.
                                TOLAKE=.FALSE.
                                NBOTTOMBASIN=NBOTTOMBASIN+1
                            ENDIF
                            EXIT L120
                        ENDIF
                    ENDIF
                    !      *******************************************************************
                    !       Check for enclosed basin [point has negative flow_direction(i,j)]
                    !      *******************************************************************
                    IF (LOCAL_DIRECTION < -1) THEN
                        IF ((I_OLD == ILOWEST).AND.(J_OLD == JLOWEST)) THEN
                        ELSE
                            TOLAKE=.TRUE.
                            NO_LAKES_EXIST=.FALSE.
                        ENDIF
                        !      *******************************************************************
                        !       If it is a depression then check if this is a new basin
                        !      *******************************************************************
                        IF (BASIN_NUMBER(I_OLD,J_OLD) == 0) THEN
                            NEWBASIN=.TRUE.
                        ENDIF
                        !      *******************************************************************
                        !       Branch if this is an enclosed basin
                        !      *******************************************************************
                        EXIT L120
                    ENDIF
                    !      *******************************************************************
                    !       If its not a depression then determine the next downstream location
                    !       Also change reference location to the next downstream point and
                    !          increment its drainage area and sediment yield
                    !       We increment in all areas except after depressions for sediment_yield but only in
                    !          bedrock channels for bederode
                    !      *******************************************************************
                    IF (LOCAL_DIRECTION > 0) THEN
                        J_NEW=J_OLD+DOWNSTREAM(LOCAL_DIRECTION,2)
                        I_NEW=I_OLD+DOWNSTREAM(LOCAL_DIRECTION,1)
                        IF (IS_X_PERIODIC) THEN
                            IF (I_NEW < 1) I_NEW = MX
                            IF (I_NEW > MX) I_NEW=1
                        ENDIF
                        IF (IS_Y_PERIODIC) THEN
                            IF (J_NEW < 1) J_NEW=MY
                            IF (J_NEW > MY) J_NEW=1
                        ENDIF
                        I_OLD=I_NEW
                        J_OLD=J_NEW
                        DRAINAGE_AREA(I_OLD,J_OLD)=DRAINAGE_AREA(I_OLD,J_OLD)+D2X+RIVERADD

                        DISCHARGE(I_OLD,J_OLD)=DISCHARGE(I_OLD,J_OLD)+QADD
                    ENDIF
                    !      *******************************************************************
                    !       Check for exit point if non-periodic lower boundary (j_old=my)
                    !            and also check for a new
                    !           basin, branching if at an exit
                    !      *******************************************************************
                    IF (DO_FLOW_BOUNDARIES) THEN
                        IF (FLOW_DIRECTION(I_OLD,J_OLD) == 1) THEN
                            IF (BASIN_NUMBER(I_OLD,J_OLD) == 0) THEN
                                NEWBASIN=.TRUE.
                                TOLAKE=.FALSE.
                                NBOTTOMBASIN=NBOTTOMBASIN+1
                                IS_EXIT_BASIN(CURRENT_BASIN+1)=.TRUE.
                            ENDIF
                            EXIT L120
                        ENDIF
                    ELSE
                        IF (IS_Y_PERIODIC) THEN
                            IF ((I_OLD == ILOWEST).AND.(J_OLD == JLOWEST)) THEN
                                IF (BASIN_NUMBER(I_OLD,J_OLD) == 0) NEWBASIN=.TRUE.
                                EXIT L120
                            ENDIF
                        ELSE
                            IF (J_OLD == MY) THEN
                                IF (BASIN_NUMBER(I_OLD,J_OLD) == 0) THEN
                                    NEWBASIN=.TRUE.
                                    TOLAKE=.FALSE.
                                    NBOTTOMBASIN=NBOTTOMBASIN+1
                                    IS_EXIT_BASIN(CURRENT_BASIN+1)=.TRUE.
                                ENDIF
                                EXIT L120
                            ENDIF
                        ENDIF
                    ENDIF
                    !      *******************************************************************
                    !       If not aborted and not an exit point, continue downstream,
                    !       incrementing the drainage area and sediment yield
                    !      *******************************************************************
                    IF (ABORT_SIMULATION) EXIT L120
                ENDDO L120
                !      *******************************************************************
                !       If we are here then we have reached the end of a basin - either
                !       an outlet or depression
                !       First check for aborted run and exit
                !      *******************************************************************
                IF (ABORT_SIMULATION) THEN
                    IABORTMIN=I_OLD-3
                    IABORTMAX=I_OLD+3
                    JABORTMIN=J_OLD-3
                    JABORTMAX=J_OLD+3
                    WRITE(OUTHIST,814) I_OLD,J_OLD,BASIN_NUMBER(I_OLD,J_OLD),DRAINAGE_AREA(I_OLD,J_OLD)
                    814         FORMAT(' *** INFINITE CYCLE AT ',2I6,' ***'  &
                    ,/,' AT BASIN NUMBER: ',I5,' AR=',G12.5)
                    RETURN
                ENDIF
                !      *******************************************************************
                !       If this is a new basin then increment the basin counter and assign
                !       the current_basin to the basin endpoint
                !      *******************************************************************
                IF (NEWBASIN) THEN
                    CURRENT_BASIN=CURRENT_BASIN+1
                    BASIN_NUMBER(I_OLD,J_OLD)=CURRENT_BASIN
                    IF (TOLAKE) THEN
                        ENCLOSED(CURRENT_BASIN)=.TRUE.
                        !      *******************************************************************
                        !      Look for surrounding points with same elevation - if points have
                        !          the same elevation then we make sure they are included in the same
                        !          basin number
                        !      Initially look only at surrounding 8 points but widen the search
                        !         if a contiguous point at the same elevation is found
                        !      *****************************************************************
                        KAROUND=1
                        IF ((FLOW_DIRECTION(I_OLD,J_OLD) < 1).OR.(D8_GRADIENT(I_OLD,J_OLD) <= 0.0)) THEN
                            L667: DO
                                ISFOUND=.FALSE.
                                IF (IS_X_PERIODIC) THEN
                                    ILOW=I_OLD-KAROUND
                                    IHIGH=I_OLD+KAROUND
                                    IF (ILOW < 1) ILOW=MX+ILOW
                                    IF (IHIGH > MX) IHIGH=IHIGH-MX
                                ELSE
                                    ILOW=MAX(1,I_OLD-KAROUND)
                                    IHIGH=MIN(MX,I_OLD+KAROUND)
                                ENDIF
                                IF (IS_Y_PERIODIC) THEN
                                    JLOW=J_OLD-KAROUND
                                    JHIGH=J_OLD+KAROUND
                                    IF (JLOW < 1) JLOW=MY+JLOW
                                    IF (JHIGH > MY) JHIGH=JHIGH-MY
                                ELSE
                                    JLOW=MAX(1,J_OLD-KAROUND)
                                    JHIGH=MIN(MYY,J_OLD+KAROUND)
                                ENDIF
                                II=ILOW
                                L665: DO
                                    JJ=JLOW
                                    L669: DO
                                        YCONTINUE=.TRUE.
                                        IF (((II == ILOW).OR.(II == IHIGH)).AND.  &
                                        ((JJ == JLOW).OR. &
                                        (JJ == JHIGH))) THEN
                                            !      *******************************************************************
                                            !       Ignore the present point and any exit points (jj=my)
                                            !      *******************************************************************
                                            IF ((II == I_OLD).AND.(JJ == J_OLD)) &
                                            YCONTINUE=.FALSE.
                                            IF (IS_Y_PERIODIC) THEN
                                                IF ((II == ILOWEST).AND.(JJ == JLOWEST)) &
                                                YCONTINUE=.FALSE. !    goto 666
                                            ELSE
                                                IF (JJ == MY) YCONTINUE=.FALSE. !    goto 666
                                            ENDIF
                                            IF (YCONTINUE) THEN
                                                IF(ELEVATION(II,JJ) == ELEVATION(I_OLD,J_OLD)) THEN
                                                    !      *******************************************************************
                                                    !       If we are looking at more than 8 points then things are a little
                                                    !       more subtle -- We need to make sure that the point at an equal
                                                    !       elevation is actually contiguous with an already identified basin
                                                    !       point - So we look around the identified point for one already
                                                    !       found
                                                    !      *******************************************************************
                                                    IF(KAROUND > 1) THEN
                                                        KKSQUARE=(KAROUND-1)**2
                                                        IF (IS_X_PERIODIC) THEN
                                                            ILL=II-1
                                                            IF (ILL < 1) ILL=ILL+MX
                                                            IUU=II+1
                                                            IF (IUU > MX) IUU=IUU-MX
                                                        ELSE
                                                            ILL=MAX(II-1,1)
                                                            IUU=MIN(II+1,MX)
                                                        ENDIF
                                                        IF (IS_Y_PERIODIC) THEN
                                                            JLL=JJ-1
                                                            IF (JLL < 1) JLL=JLL+MY
                                                            JUU=JJ+1
                                                            IF (JUU > MY) JUU=JUU-MY
                                                        ELSE
                                                            JLL=MAX(JJ-1,1)
                                                            JUU=MIN(JJ+1,MY-1)
                                                        ENDIF
                                                        CONTIGUOUS=.FALSE.
                                                        III=ILL
                                                        L664: DO
                                                            JJJ=JLL
                                                            L670: DO
                                                                XCONTINUE=.TRUE.
                                                                IF ((III == II).AND.(JJJ == JJ)) &
                                                                XCONTINUE=.FALSE.
                                                                IF ((((III-II)**2) <= KKSQUARE)  &
                                                                .AND.(((JJJ-JJ)**2) <= KKSQUARE)) &
                                                                XCONTINUE=.FALSE.
                                                                IF (IS_Y_PERIODIC) THEN
                                                                    IF ((II == ILOWEST).AND. &
                                                                    (JJ == JLOWEST)) &
                                                                    XCONTINUE=.FALSE.
                                                                ELSE
                                                                    IF (JJ == MY) &
                                                                    XCONTINUE=.FALSE.
                                                                ENDIF
                                                                IF (XCONTINUE) THEN
                                                                    IF (ELEVATION(III,JJJ) == ELEVATION(II,JJ)) THEN
                                                                        IF (  (((III-II)**2) <= 1).AND.(((JJJ-JJ)**2) <= 1)) THEN
                                                                            CONTIGUOUS=.TRUE.
                                                                        ELSE
                                                                            EXIT L667
                                                                        ENDIF
                                                                    ENDIF
                                                                ENDIF
                                                                IF (JJJ /= JUU) THEN
                                                                    JJJ=JJJ+1
                                                                    IF (IS_Y_PERIODIC) THEN
                                                                        IF (JJJ > MY) JJJ=1
                                                                    ENDIF
                                                                    CYCLE L670
                                                                ENDIF
                                                                EXIT L670
                                                            ENDDO L670
                                                            IF (III == IUU) EXIT L664
                                                            III=III+1
                                                            IF (IS_X_PERIODIC) THEN
                                                                IF (III > MX) III=1
                                                            ENDIF
                                                        ENDDO L664
                                                    ELSE
                                                        CONTIGUOUS=.TRUE.
                                                    ENDIF
                                                    !      *******************************************************************
                                                    !       If the equal-elevation point is contiguous it is given the same
                                                    !       basin number
                                                    !      *******************************************************************
                                                    IF (CONTIGUOUS) THEN
                                                        BASIN_NUMBER(II,JJ)=CURRENT_BASIN
                                                        NCONTIG=NCONTIG+1
                                                        IF ((.NOT.IS_Y_PERIODIC).AND. &
                                                        (JJ == MY-1)) THEN
                                                        ELSE
                                                            ISFOUND=.TRUE.
                                                        ENDIF
                                                    ENDIF
                                                ENDIF
                                            ENDIF
                                        ENDIF
                                        IF (JJ /= JHIGH) THEN
                                            JJ=JJ+1
                                            IF (IS_Y_PERIODIC) THEN
                                                IF (JJ > MY) JJ=1
                                            ENDIF
                                            CYCLE L669
                                        ENDIF
                                        EXIT L669
                                    ENDDO L669
                                    IF (II == IHIGH) EXIT L665
                                    II=II+1
                                    IF (IS_X_PERIODIC) THEN
                                        IF (II > MX) II=1
                                    ENDIF
                                ENDDO L665
                                KAROUND=KAROUND+1
                                IF (KAROUND > NMAX) THEN
                                    RETURN
                                ENDIF
                                !      *******************************************************************
                                !        We widen the search if a contiguous point at the same elevation
                                !        has been found to include points possibly contiguous to the new
                                !        point but not the original point
                                !      *******************************************************************
                                IF (.NOT.ISFOUND) EXIT L667
                            ENDDO L667
                        ENDIF
                        !      *******************************************************************
                        !        End of code relating to grouping surrounding points at the same
                        !            elevation
                        !      *******************************************************************
                    ENDIF
                ENDIF
                !      *******************************************************************
                !       We need to retrace our steps from the initial location, setting the
                !       basin number of all points along the path to that of the exit or sink point
                !      *******************************************************************
                NUMBBASIN=BASIN_NUMBER(I_OLD,J_OLD)
                I_OLD=I
                J_OLD=J
                L140: DO
                    !      *******************************************************************
                    !       End if we reach an already-identified basin location
                    !      *******************************************************************
                    IF (BASIN_NUMBER(I_OLD,J_OLD) /= 0) THEN
                        EXIT
                    ELSE
                        BASIN_NUMBER(I_OLD,J_OLD)=NUMBBASIN
                        LOCAL_DIRECTION = FLOW_DIRECTION(I_OLD,J_OLD)
                        J_NEW=J_OLD+DOWNSTREAM(LOCAL_DIRECTION,2)
                        I_NEW=I_OLD+DOWNSTREAM(LOCAL_DIRECTION,1)
                        IF (DO_FLOW_BOUNDARIES) THEN
                            IF (FLOW_DIRECTION(I_NEW,J_NEW) == 1) THEN
                                EXIT
                            ENDIF
                        ELSE
                            IF (IS_X_PERIODIC) THEN
                                IF (I_NEW < 1) I_NEW = MX
                                IF (I_NEW > MX) I_NEW=1
                            ENDIF
                            IF (IS_Y_PERIODIC) THEN
                                IF (J_NEW < 1) J_NEW=MY
                                IF (J_NEW > MY) J_NEW=1
                            ELSE
                                IF (J_NEW == MY) EXIT
                            ENDIF
                        ENDIF
                        I_OLD=I_NEW
                        J_OLD=J_NEW
                    ENDIF
                ENDDO L140
            ENDDO M110
        ENDDO L110
        !      *******************************************************************
        !       The initial drainage areas and basin numbers have now been determined
        !       Check to make sure that we have not exceeded the dimensioned
        !       number of basins
        !       Also we have determined the sediment yield from channel erosion
        !      *******************************************************************

        IF (CURRENT_BASIN > LMMX) THEN
            WRITE(OUTHIST,933) CURRENT_BASIN
            933     FORMAT(' TOO MANY BASINS - REDIMENSION.  CURRENT_BASIN=',I6)
            RETURN
        ENDIF
        !      *******************************************************************
        !       Look at boundary points and set values for those that have not
        !       received drainage from upstream to unity area and assign a basin number
        !      *******************************************************************
        IF (DO_FLOW_BOUNDARIES) THEN
            DO J=1,MY
                DO I=1,MX
                    IF (FLOW_DIRECTION(I,J) == 1) THEN
                        IF (DISCHARGE(I,J) == 0.0) THEN
                            DISCHARGE(I,J)=DISCHARGE_FROM_CELL(I,J)
                        ENDIF
                        IF (BASIN_NUMBER(I,J) == 0) THEN
                            CURRENT_BASIN=CURRENT_BASIN+1
                            BASIN_NUMBER(I,J)=CURRENT_BASIN
                            NBOTTOMADD=NBOTTOMADD+1
                            IS_EXIT_BASIN(CURRENT_BASIN)=.TRUE.
                        ENDIF
                    ENDIF
                ENDDO
            ENDDO
        ELSE
            IF (IS_Y_PERIODIC) THEN
                IF (DRAINAGE_AREA(ILOWEST,JLOWEST) == 0) THEN
                    DRAINAGE_AREA(ILOWEST,JLOWEST)=D2X
                    DISCHARGE(ILOWEST,JLOWEST)=DISCHARGE_FROM_CELL(ILOWEST,JLOWEST)
                ENDIF
                IF (BASIN_NUMBER(ILOWEST,JLOWEST) == 0) THEN
                    CURRENT_BASIN=CURRENT_BASIN+1
                    BASIN_NUMBER(ILOWEST,JLOWEST)=CURRENT_BASIN
                ENDIF
            ELSE
                DO  I=1,MX
                    IF (DRAINAGE_AREA(I,MY) == 0) THEN
                        DRAINAGE_AREA(I,MY)=D2X
                        DISCHARGE(I,MY)=DISCHARGE_FROM_CELL(I,MY)
                    ENDIF
                    IF (BASIN_NUMBER(I,MY) == 0) THEN
                        CURRENT_BASIN=CURRENT_BASIN+1
                        BASIN_NUMBER(I,MY)=CURRENT_BASIN
                        NBOTTOMADD=NBOTTOMADD+1
                        IS_EXIT_BASIN(CURRENT_BASIN)=.TRUE.
                    ENDIF
                    IF (.NOT.IS_EXIT_BASIN(BASIN_NUMBER(I,MY))) THEN
                        WRITE(*,8227) 
                        8227 FORMAT(' PROBLEM IN LOWER BOUNDARY BASIN IDENTIFICATION')
                        IS_EXIT_BASIN(BASIN_NUMBER(I,MY))=.TRUE.
                    ENDIF
                ENDDO
            ENDIF
        ENDIF
        TOTAL_BASINS=CURRENT_BASIN
        DO  I=2,MX-1
            DO  J=2,MYY

                IF (BASIN_NUMBER(I,J) == 0) THEN
                    DO  JJ=J-1,J+1
                        DO  II=I-1,I+1
                            WRITE (OUTHIST,4955) II,JJ,FLOW_DIRECTION(II,JJ),BASIN_NUMBER(II,JJ)
                            4955 FORMAT('NOBASIN, II=',I5,' JJ=',I5,' ID=',I4,' N=',I5)
                        ENDDO
                    ENDDO
                    WRITE (OUTHIST,4956)
                    4956  FORMAT(' ')
                ENDIF
            ENDDO
        ENDDO
        IF (.NOT.IS_Y_PERIODIC) THEN
            DO  K=1,TOTAL_BASINS
                KKK=0
                DO  J=1,MY
                    DO I=1,MX
                        IF (BASIN_NUMBER(I,J) == K) KKK=KKK+1
                    ENDDO
                ENDDO
                IF (KKK == 0) ENCLOSED(K)=.FALSE.
                IF ((.NOT.ENCLOSED(K)).AND.(.NOT.IS_EXIT_BASIN(K))) THEN
                    9221   FORMAT('BAD NOENCLOSE=',I7)
                ENDIF
            ENDDO
        ENDIF
        ROUTED_DISCHARGE=DISCHARGE
        !      *******************************************************************
        !       If there are enclosed basins (*and if runoff overflow is permitted*)
        !       then we identify those enclosed basins which need to be examined
        !       for grouping and flow routing to the outlets - skip all this
        !       otherwise
        !      *******************************************************************
        IF ((COMPLETE_RUNOFF).AND.(.NOT.NO_LAKES_EXIST)) THEN
            DO  K=1,TOTAL_BASINS
                IF (ENCLOSED(K)) THEN
                    NEEDTODO(K)=.TRUE.
                ELSE
                    NEEDTODO(K)=.FALSE.
                ENDIF
            ENDDO
            !      *******************************************************************
            !       We need to keep regrouping basins and identifying lake outlets and
            !       cumulative drainage areas in the case where there is enclosed basins
            !       and lake overflow - redo is an index of whether another level of
            !       basin consolidation is needed
            !      *******************************************************************
            L250: DO
                REDO=.FALSE.
                !      *******************************************************************
                !       If we still need to work in a basin (an outlet location has still
                !       not been determined) then we reset the outlet elevation and
                !       lake area
                !      *******************************************************************
                L150: DO  K=1,TOTAL_BASINS
                    IF (NEEDTODO(K)) THEN
                        LAKE_OUTLET_ELEVATION(K)=1.0E+25
                        !      *******************************************************************
                        !         Determine the present lake drainage area
                        !      *******************************************************************
                        BASIN_DRAINAGE_AREA(K)=0
                        LOCAL_BASIN_DISCHARGE(K)=0.0
                        L160: DO  J=1,MYY
                            M160: DO  I=1,MX
                                IF (BASIN_NUMBER(I,J) == K) THEN
                                    IF (FLOW_DIRECTION(I,J) < 0) THEN
                                        BASIN_DRAINAGE_AREA(K)=   &
                                        BASIN_DRAINAGE_AREA(K)+DRAINAGE_AREA(I,J)
                                        LOCAL_BASIN_DISCHARGE(K)=LOCAL_BASIN_DISCHARGE(K)+DISCHARGE(I,J)
                                    ENDIF
                                    !      *******************************************************************
                                    !       We need to look at all points on the periphery of the basin to
                                    !       determine if it is an outlet
                                    !      *******************************************************************
                                    L170: DO  L=2,9
                                        J_NEW=J+DOWNSTREAM(L,2)
                                        IF (IS_Y_PERIODIC) THEN
                                            IF (J_NEW < 1) J_NEW = MY
                                            IF (J_NEW > MY) J_NEW=1
                                        ELSE
                                            IF (J_NEW < 1) CYCLE L170
                                        ENDIF
                                        I_NEW=I+DOWNSTREAM(L,1)
                                        IF (IS_X_PERIODIC) THEN
                                            IF (I_NEW < 1) I_NEW = MX
                                            IF (I_NEW > MX) I_NEW=1
                                        ELSE
                                            IF ((I_NEW < 1).OR.(I_NEW > MX)) CYCLE L170
                                        ENDIF
                                        !      *******************************************************************
                                        !       Do the following if the present point is at a peripheral location
                                        !       Check to see if the periphery point or a contiguous point outside
                                        !         the basin is lower than any identified lake outlet
                                        !      *******************************************************************
                                        IF (BASIN_NUMBER(I_NEW,J_NEW) /= K) THEN
                                            XCOMP=MAX(ELEVATION(I_NEW,J_NEW),ELEVATION(I,J))
                                            !      *******************************************************************
                                            !       If this is a potential lake outlet then there are two cases:
                                            !         1.  The basin peripery is lower than any contiguous point outside
                                            !          the basin.  If this is the case (outer is true) the outlet is
                                            !             identified with the outside point and the outlet elevation is
                                            !             equal to that of the outside point.
                                            !         2.  The basin periphery is higher than the lowest contiguous point
                                            !          outside.  The outlet is the periphery point and the lake outlet
                                            !          elevation is its elevation.  In this case we need to record where
                                            !             the flow is going to and its direction (ifrom,jfrom,kfrom)
                                            !      *******************************************************************
                                            IF (XCOMP < LAKE_OUTLET_ELEVATION(K)) THEN
                                                IF (ELEVATION(I,J) <= ELEVATION(I_NEW,J_NEW)) THEN
                                                    OUTER(K)=.TRUE.
                                                    IOUTLET=I_NEW
                                                    JOUTLET=J_NEW
                                                    LAKE_OUTLET_ELEVATION(K)=ELEVATION(I_NEW,J_NEW)
                                                ELSE
                                                    OUTER(K)=.FALSE.
                                                    IOUTLET=I
                                                    JOUTLET=J
                                                    LAKE_OUTLET_ELEVATION(K)=ELEVATION(I,J)
                                                    IFROM=I_NEW
                                                    JFROM=J_NEW
                                                    KKFROM=L
                                                ENDIF
                                                !      *******************************************************************
                                                !       downstream_basin is the basin into which the present basin drains
                                                !      *******************************************************************
                                                DOWNSTREAM_BASIN(K)=BASIN_NUMBER(I_NEW,J_NEW)
                                            ENDIF
                                        ENDIF
                                    ENDDO L170
                                ENDIF
                            ENDDO M160
                        ENDDO L160
                        !      *******************************************************************
                        !       We have now identified the lowest outlet for the present basin
                        !       so we set its outlet location.  If the basin outlet is on the
                        !       periphery (not outer) then we need to redefine the drainage direction
                        !       and gradient of the peripery point for flow to the outlet rather
                        !       than into the present basin
                        !      *******************************************************************
                        I_OUTFLOW(K)=IOUTLET
                        J_OUTFLOW(K)=JOUTLET
                        IF (.NOT.OUTER(K)) THEN
                            FLOW_DIRECTION(IOUTLET,JOUTLET)=KKFROM
                            D8_GRADIENT(IOUTLET,JOUTLET)=(ELEVATION(IOUTLET,JOUTLET) &
                            -ELEVATION(IFROM,JFROM))/CELL_SIZE
                        ENDIF
                    ENDIF
                    !      *******************************************************************
                    !       We have now done this basin - we don't need to worry about it for
                    !       the time being.
                    !      *******************************************************************
                    NEEDTODO(K)=.FALSE.
                ENDDO L150
                !      *******************************************************************
                !       We have now completed one stage of basin agglomeration.  We need
                !       to check for the case where lake basins drain into each other
                !       if this is true then redo is set true and needtodo(k) is set true
                !       redo is false if each chain of enclosed basin outlets following downstream
                !       eventually leads to an external outlet
                !      *******************************************************************
                K=0
                KKK=0
                !      *******************************************************************
                !         Cycle through all basins
                !      *******************************************************************
                L260: DO
                    K=K+1
                    KKCOUNT=0
                    !      *******************************************************************
                    !         Exit from loop if we have looked at all basins
                    !      *******************************************************************
                    IF (K >= TOTAL_BASINS) EXIT L260
                    !      *******************************************************************
                    !         If we find an enclosed basin we need to check that tracing successive
                    !         outlets downstream eventually leads to an external outlet, otherwise
                    !         we need an additional stage of basin agglomeration
                    !      *******************************************************************
                    IF (ENCLOSED(K)) THEN
                        KK=K
                        IF ((K > LMMX).OR.(KK > LMMX).OR.(KKK > LMMX) &
                        .OR.(K == 0).OR.(KK == 0)) THEN
                            WRITE (*,4623) K,KK,KKK,DOWNSTREAM_BASIN(KK),  &
                            LMMX
                            4623  FORMAT(' ABORTING, K=',I7,' KK=',I7,' KKK=',I7, &
                            ' DOWNSTREAM_BASIN(KK)=',I7,' LMMX=',I7)
                            CALL BASIN_REPORT()
                            RETURN
                        ENDIF
                        L623: DO
                            KKK=DOWNSTREAM_BASIN(KK)
                            IF ((KKK > LMMX).OR.(KKK == 0)) THEN
                                WRITE (*,5623) KK,KKK
                                5623  FORMAT(' ABORTING, KK=',I7,' KKK=',I7)
                                CALL BASIN_REPORT()
                                RETURN
                            ENDIF

                            IF (DOWNSTREAM_BASIN(KKK) == K) EXIT L623
                            IF (.NOT.ENCLOSED(KKK)) CYCLE L260
                            IF (DO_FLOW_BOUNDARIES) THEN
                                IF ((I_OUTFLOW(KKK) == 1).OR.(I_OUTFLOW(KKK) == MX).OR. &
                                (J_OUTFLOW(KKK) == 1).OR.(J_OUTFLOW(KKK) == MY)) THEN
                                    CYCLE L260
                                ENDIF
                            ELSE
                                IF (IS_Y_PERIODIC) THEN
                                    IF ((I_OUTFLOW(KKK) == ILOWEST).AND.(J_OUTFLOW(KKK) == JLOWEST)) &
                                    CYCLE L260
                                ELSE
                                    IF (J_OUTFLOW(KKK) == MY) CYCLE L260 !    goto 260
                                ENDIF
                            ENDIF
                            KK=KKK
                            KKCOUNT=KKCOUNT+1
                            IF (KKCOUNT > 500) CYCLE L260
                        ENDDO L623
                        !      *******************************************************************
                        !       This is the case of mutually-draining basins
                        !      *******************************************************************
                        KK=K
                        KEND=KKK
                        REDO=.TRUE.
                        NEEDTODO(K)=.TRUE.
                        KKK=DOWNSTREAM_BASIN(KK)
                        !      *******************************************************************
                        !       If there are two mutually-draining basins, then they are arbitrarily
                        !       agglomerated by reassigning one basins points (kkk) to the other (k)
                        !      *******************************************************************
                        DO  J=1,MYY
                            DO  I=1,MX
                                IF (BASIN_NUMBER(I,J) == KKK) BASIN_NUMBER(I,J)=K
                            ENDDO
                        ENDDO
                        !      *******************************************************************
                        !       If other basins drain to the reassigned drainage basin, then their
                        !         exit basin have to be reassined to the agglomerated basin
                        !      *******************************************************************
                        L627: DO  KKKK=1,TOTAL_BASINS
                            IF (ENCLOSED(KKKK)) THEN
                                IF (KKKK == K) CYCLE L627
                                IF (KKKK == KKK) CYCLE L627
                                IF (DOWNSTREAM_BASIN(KKKK) == KKK) DOWNSTREAM_BASIN(KKKK)=K
                            ENDIF
                        ENDDO L627
                        !      *******************************************************************
                        !       We need to deactivate the various values of the deactivated
                        !          (agglomerated) basin (kkk)
                        !      *******************************************************************
                        ENCLOSED(KKK)=.FALSE.
                        BASIN_DRAINAGE_AREA(KKK)=0.0
                        LOCAL_BASIN_DISCHARGE(KKK)=0.0
                        I_OUTFLOW(KKK)=-1
                        J_OUTFLOW(KKK)=-1
                        KK=DOWNSTREAM_BASIN(KKK)
                        TOBEADDED(KKK)=0
                        DOWNSTREAM_BASIN(KKK)=0
                        !      *******************************************************************
                        !       Look at additional basins if we have not examined all of them
                        !      *******************************************************************
                        IF (KKK == KEND) CYCLE L260
                        KKK=KK
                    ENDIF
                    !      *******************************************************************
                    !       Go on to the next basin
                    !      *******************************************************************
                ENDDO L260
                !      *******************************************************************
                !       if we have not completely agglomerated drainage basins so that they
                !       flow to the exterior, then cycle back to do it again
                !      *******************************************************************
                IF (.NOT.REDO) THEN
                    EXIT L250
                ENDIF
            ENDDO L250
            !      *******************************************************************
            !       We now need to route flow downstream through all lakes, such that
            !       the drainage area downstream from the outlet includes the lake area
            !       and all drainage basins and upstream lakes draining into the given
            !       lake
            !      *******************************************************************
            L750: DO K=1,TOTAL_BASINS
                KK=K
                NNN=0
                !      *******************************************************************
                !         Check for enclosed drainage basins
                !      *******************************************************************
                L760: DO
                    IF (ENCLOSED(KK)) THEN
                        KKKK=1
                        IF (ENCLOSED(DOWNSTREAM_BASIN(KK))) THEN
                            KKK=1
                        ELSE
                            KKK=0
                        ENDIF
                    ELSE
                        KKKK=0
                    ENDIF
                    IF (ENCLOSED(KK)) THEN
                        !      *******************************************************************
                        !         If its outlet basin is also enclosed, then we need to continue
                        !         downstream.  tobeadded(k) is the cumulative added drainage from
                        !         enclosed basins upstream from the given basin
                        !      *******************************************************************
                        IF (ENCLOSED(DOWNSTREAM_BASIN(KK))) THEN
                            TOBEADDED(DOWNSTREAM_BASIN(KK))=  &
                            TOBEADDED(DOWNSTREAM_BASIN(KK))+BASIN_DRAINAGE_AREA(K)
                            QTOADD(DOWNSTREAM_BASIN(KK))=QTOADD(DOWNSTREAM_BASIN(KK)) &
                            +LOCAL_BASIN_DISCHARGE(K)
                            KK=DOWNSTREAM_BASIN(KK)
                            3444 FORMAT(' KK=',I7,' TO',I7,' ENC=',I5,' TO ENC=',I5,' N=',I7)
                            NNN=NNN+1
                            !      *******************************************************************
                            !       Check for infinite loops between basins (basins draining into one
                            !       another or through a loop of connecting basins).  Abort if this
                            !       occurs
                            !      *******************************************************************
                            IF (NNN > TOTAL_BASINS) THEN
                                WRITE(*,799) ITERATION,KK
                                WRITE(OUTHIST,799) ITERATION,KK
                                799                 FORMAT(' **** BASIN SCRAMBLE AT ITERATION',I5, &
                                ' AT BASIN',I5)
                                CALL BASIN_REPORT()
                                RETURN
                            ENDIF
                            CYCLE L760
                        ENDIF
                        GOTO 8443
                    ELSE
                        GOTO 8443
                    ENDIF
                ENDDO L760
                8443 CONTINUE
                IF (ENCLOSED(K)) THEN
                    IF(ENCLOSED(DOWNSTREAM_BASIN(K))) THEN
                        IF (ENCLOSED(DOWNSTREAM_BASIN(KK))) THEN
                            KKK=1
                        ELSE
                            KKK=0
                        ENDIF
                    ENDIF
                ENDIF
            ENDDO L750
            !      *******************************************************************
            !       Now that we have calculated contributions from upstream enclosed basins
            !         we need to calculate the total lake drainage area of each basin
            !      *******************************************************************
            DO  K=1,TOTAL_BASINS
                IF (ENCLOSED(K)) THEN
                    BASIN_DRAINAGE_AREA(K)=BASIN_DRAINAGE_AREA(K)+TOBEADDED(K)
                    LOCAL_BASIN_DISCHARGE(K)=LOCAL_BASIN_DISCHARGE(K)+QTOADD(K)
                ELSE
                    IF(.NOT.IS_Y_PERIODIC) THEN
                        LAKE_OUTLET_ELEVATION(K)=ELEVATION(1,MY)
                    ENDIF
                ENDIF
            ENDDO
            !      *******************************************************************
            !       We now go downstream from each enclosed basin outlet and increment the
            !       drainage area by the amount of drainage flowing through the outlet
            !      *******************************************************************
            L180: DO K=1,TOTAL_BASINS
                IF (ENCLOSED(K)) THEN
                    I_OLD=I_OUTFLOW(K)
                    J_OLD=J_OUTFLOW(K)
                    IF ((I_OLD == I_OUTFLOW(DOWNSTREAM_BASIN(K))).AND. &
                    (J_OLD == J_OUTFLOW(DOWNSTREAM_BASIN(K)))) &
                    CYCLE L180 !    go to 180
                    IF ((I_OLD < 1).OR.(I_OLD > MX).OR. &
                    (J_OLD < 1).OR.(J_OLD > MY)) THEN
                        WRITE(OUTHIST,377) K,I_OLD,J_OLD
                        377 FORMAT(' ***BASIN INTEGRITY PROBLEM AT BASIN' &
                        ,I5,', I_OLD & J_OLD=', &
                        2I5)
                    ENDIF
                    ADDAREA=BASIN_DRAINAGE_AREA(K)
                    ADDQ=LOCAL_BASIN_DISCHARGE(K)
                    IF(OUTER(K)) THEN
                        DRAINAGE_AREA(I_OLD,J_OLD)=DRAINAGE_AREA(I_OLD,J_OLD)+ADDAREA
                        DISCHARGE(I_OLD,J_OLD)=DISCHARGE(I_OLD,J_OLD)+ADDQ
                        ROUTED_DISCHARGE(I_OLD,J_OLD)=DISCHARGE(I_OLD,J_OLD)
                    ELSE
                        DRAINAGE_AREA(I_OLD,J_OLD)=ADDAREA
                        DISCHARGE(I_OLD,J_OLD)=ADDQ
                        ROUTED_DISCHARGE(I_OLD,J_OLD)=DISCHARGE(I_OLD,J_OLD)
                    ENDIF
                    IF (DO_FLOW_BOUNDARIES) THEN
                        IF (FLOW_DIRECTION(I_OLD,J_OLD) == 1) THEN
                            CYCLE L180
                        ENDIF
                    ELSE
                        IF (IS_Y_PERIODIC) THEN
                            IF ((I_OLD == ILOWEST).AND. &
                            (J_OLD == JLOWEST)) CYCLE L180
                        ELSE
                            IF (J_OLD == MY) CYCLE L180
                        ENDIF
                    ENDIF
                    LOCAL_DIRECTION = FLOW_DIRECTION(I_OLD,J_OLD)
                    IF (LOCAL_DIRECTION < 0) CYCLE L180
                    L190: DO
                        J_NEW=J_OLD+DOWNSTREAM(LOCAL_DIRECTION,2)
                        I_NEW=I_OLD+DOWNSTREAM(LOCAL_DIRECTION,1)
                        IF (IS_X_PERIODIC) THEN
                            IF (I_NEW < 1) I_NEW = MX
                            IF (I_NEW > MX) I_NEW=1
                        ENDIF
                        IF (IS_Y_PERIODIC) THEN
                            IF (J_NEW < 1) J_NEW=MY
                            IF (J_NEW > MY) J_NEW=1
                        ENDIF
                        I_OLD=I_NEW
                        J_OLD=J_NEW
                        IF (BASIN_NUMBER(I_OLD,J_OLD) /= DOWNSTREAM_BASIN(K)) THEN
                            WRITE(OUTHIST,823) K,DOWNSTREAM_BASIN(K),&
                            BASIN_NUMBER(I_OLD,J_OLD),I_OLD,J_OLD &
                            ,DRAINAGE_AREA(I_OLD,J_OLD),ADDAREA
                            823               FORMAT(' PROBLEM IN ADDAREA FROM BASIN',I5, &
                            ' TO BASIN' &
                            ,I5,/,' NOW IN BASIN',I5,' AT I=',I5,' J=',I5, &
                            ' AREA=',G12.5,' ADDAREA=',I5)
                            EXIT L190 !    goto 200
                            CALL BASIN_REPORT()
                            ABORT_SIMULATION=.TRUE.
                            RETURN
                        ENDIF
                        IF ((I_OLD == I_OUTFLOW(DOWNSTREAM_BASIN(K))).AND.(J_OLD == J_OUTFLOW(DOWNSTREAM_BASIN(K)))) &
                        EXIT L190
                        DRAINAGE_AREA(I_OLD,J_OLD)=DRAINAGE_AREA(I_OLD,J_OLD)+ADDAREA
                        DISCHARGE(I_OLD,J_OLD)=DISCHARGE(I_OLD,J_OLD)+ADDQ
                        ROUTED_DISCHARGE(I_OLD,J_OLD)=DISCHARGE(I_OLD,J_OLD)
                        LOCAL_DIRECTION = FLOW_DIRECTION(I_OLD,J_OLD)
                        IF (LOCAL_DIRECTION < 0) EXIT L190
                        IF (DO_FLOW_BOUNDARIES) THEN
                            IF (LOCAL_DIRECTION == 1) THEN
                                EXIT L190
                            ENDIF
                        ELSE
                            IF (IS_Y_PERIODIC) THEN
                                IF ((I_OLD == ILOWEST).AND. &
                                (J_OLD == JLOWEST)) EXIT L190
                            ELSE
                                IF (J_OLD == MY) EXIT L190
                            ENDIF
                        ENDIF
                    ENDDO L190
                ENDIF
            ENDDO L180
            !      *******************************************************************
            !       Finished with basin and drainage area determination
            !      *******************************************************************
        ENDIF
        ADDAREA=0
        ADDQ=0.0
        IF (.NOT.IS_Y_PERIODIC) THEN
            DO  I=1,MX
                ADDQ=ADDQ+DISCHARGE(I,MY)
                ADDAREA=ADDAREA+DRAINAGE_AREA(I,MY)
            ENDDO
            IF (ADDAREA < ((MX-1)*(MY-1)-50)) THEN
                WRITE (OUTHIST,889) ITERATION,ADDAREA
                889     FORMAT(' ITER=',I5,' TOTAL AREA = ',G12.5)
            ENDIF
        ENDIF
        DO J=1,MY
            DO I=1,MX
                IF (FLOW_DIRECTION(I,J) < 0) THEN
                    K=BASIN_NUMBER(I,J)
                    IF (ENCLOSED(K)) THEN
                        IF ((I_OUTFLOW(K) == I).AND.(J_OUTFLOW(K) == J)) THEN
                            KKK=0
                            DO JJ=1,MY
                                DO II=1,MX
                                    IF (BASIN_NUMBER(II,JJ) == K) KKK=KKK+1
                                ENDDO
                            ENDDO
                            WRITE(*,2890) I,J,K,KKK
                            2890  FORMAT('DEPRESSION IS OUTFLOW, I=',I5,' J=',I5, &
                            ' BN=',I6,' SIZE=',I6)
                        ENDIF
                    ELSE
                        !WRITE(*,2889) I,J,K
                        2889  FORMAT('NON-ENCLOSED DEPRESSION AT I=',I5,' J=',I5,' BN=',I6)
                    ENDIF
                ENDIF
            ENDDO
        ENDDO
        WSUM=0.0
        WNUM=0.0
        !      ********************************************************************
        !        If discharge is not directly proportional to drainage area, then
        !        correct for non-linearity
        !      ********************************************************************
        DO J=1,MY
            DO I=1,MX
                IF (DISCHARGE(I,J) > GLOBAL_MAXIMUM_DISCHARGE) GLOBAL_MAXIMUM_DISCHARGE=DISCHARGE(I,J)
            ENDDO
        ENDDO
        GLOBAL_MAXIMUM_DISCHARGE=GLOBAL_MAXIMUM_DISCHARGE/(D2X*DISCHARGE_SCALE_FACTOR)
        8549  FORMAT(' MAX. YIELD=',G12.5)
        IF (RESCALE_DISCHARGES)  THEN
            DISCHARGE=DISCHARGE**DISCHARGE_EXPONENT
        ENDIF
        DRAINAGE_AREA=DRAINAGE_AREA+AREAFACTOR
        IF (MODEL_GROUNDWATER) THEN
            IF (NSEEPAVG > 0.0) THEN
                AVGRUNOFF=AVGRUNOFF*D2X/NSEEPAVG
                AVGSEEP=AVGSEEP*D2X/NSEEPAVG
                IF (NFRACTAVG > 0.0) THEN
                    AVGSEEPFRACT=AVGSEEPFRACT/NFRACTAVG
                ENDIF
                WRITE(OUTHIST,2881) AVGRUNOFF,AVGSEEP,AVGSEEPFRACT,NSEEPAVG,NFRACTAVG
                2881  FORMAT(' AVG. RUNOFF=',G12.5,' AVG. SEEP FLOW=',G12.5,/, &
                ' AVG. FRACTION OF FLOW FROM SEEPAGE=',G12.5,' N=',G12.5,' N1=',G12.5)
            ENDIF
        ENDIF
        !     **********************************************************************
        !      Determine channel widths
        !     **********************************************************************
        CHANNEL_WIDTH=CHANNEL_WIDTH_CONSTANT*DISCHARGE**CHANNEL_WIDTH_EXPONENT
        !     **********************************************************************
        !      Determine sediment transport through each cell based upon elevation
        !         changes locally and upstream during last iteration.  Sediment is not
        !         routed through depressions
        !     **********************************************************************
        ERODETOADD=0.0
        ERGRAVTOADD=0.0
        !      *********************************************************************
        !       sedadd is the sediment load added to the stream flow by local bedrock
        !         erosion
        !       tsedcont sums up the sediment flux contributions to
        !         local sediment yield - used for debugging
        !         CHANNEL_EROSION_RATE holds the local bedrock channel erosion from the last iteration
        !      *********************************************************************
        !      *********************************************************************
        !        sediment_yield is accumulated downstream - it is assumed that all sediment in
        !         bedrock channels moves instantaneously downstream - i.e., moves much more
        !         rapidly than the timescale for vertical elevation changes
        !      *********************************************************************
        L710: DO  J=1,MYY
            M710: DO  I=1,MX
                I_OLD=I
                J_OLD=J
                IF (.NOT.IS_SEDIMENT_COVERED(I,J)) THEN
                    IF (CHANNEL_EROSION_RATE(I,J) > 0.0) THEN
                        SEDADD=CHANNEL_EROSION_RATE(I,J)*D2X
                    ELSE
                        SEDADD=0.0
                    ENDIF
                    !      *********************************************************************
                    !       If there is an incoming river from out of the matrix,
                    !         add in its sediment supply
                    !      *********************************************************************
                    IF (HAVE_INFLUENT_RIVERS.AND.IS_INFLUENT_RIVER_LOCATION(I,J)) THEN
                        DO KL=1,NUMBER_OF_INFLUENT_RIVERS
                            IF((I == I_INFLUENT_RIVER(KL)).AND.(J == J_INFLUENT_RIVER(KL)))THEN
                                SEDADD=SEDADD+INFLUENT_RIVER_SEDIMENT_LOAD(KL)
                            ENDIF
                        ENDDO
                    ENDIF
                    TSEDCONT=TSEDCONT+SEDADD
                    TOSEDCOVER=.FALSE.
                ELSE
                    CYCLE M710
                ENDIF
                SEDIMENT_YIELD(I,J)=SEDIMENT_YIELD(I,J)+SEDADD
                L720: DO
                    LOCAL_DIRECTION = FLOW_DIRECTION(I_OLD,J_OLD)
                    !     **********************************************************************
                    !      Branch if a drainage exit or depression
                    !     **********************************************************************
                    IF (LOCAL_DIRECTION < 2) THEN
                        !      *********************************************************************
                        !       Increments total work and transport variables
                        !      *********************************************************************
                        IF (DO_FLUVIAL_DETACHMENT) THEN
                            ERODETOADD=ERODETOADD+SEDADD*SQRT(REAL((I_OLD-I)**2+ &
                            (J_OLD-J)**2))
                            ERGRAVTOADD=ERGRAVTOADD+SEDADD*(ELEVATION(I,J)-ELEVATION(I_OLD,J_OLD))
                        ENDIF
                        EXIT L720
                    ENDIF
                    J_NEW=J_OLD+DOWNSTREAM(LOCAL_DIRECTION,2)
                    I_NEW=I_OLD+DOWNSTREAM(LOCAL_DIRECTION,1)
                    IF (IS_X_PERIODIC) THEN
                        IF (I_NEW < 1) I_NEW=MX
                        IF (I_NEW > MX) I_NEW=1
                    ENDIF
                    IF (IS_Y_PERIODIC) THEN
                        IF (J_NEW > MY) J_NEW=1
                        IF (J_NEW < 1) J_NEW=MY
                    ENDIF
                    IF (IS_SEDIMENT_COVERED(I_NEW,J_NEW)) THEN
                        !      *********************************************************************
                        !       Increments total work and transport variables
                        !      *********************************************************************
                        IF (DO_FLUVIAL_DETACHMENT.AND.(.NOT.TOSEDCOVER)) THEN
                            ERODETOADD=ERODETOADD+SEDADD*SQRT(REAL((I_OLD-I)**2+ &
                            (J_OLD-J)**2))
                            ERGRAVTOADD=ERGRAVTOADD+SEDADD*(ELEVATION(I,J)-ELEVATION(I_OLD,J_OLD))
                        ENDIF
                        TOSEDCOVER=.TRUE.
                    ELSE
                        IF (TOSEDCOVER) EXIT L720 !    goto 710
                    ENDIF
                    I_OLD=I_NEW
                    J_OLD=J_NEW
                    !     **********************************************************************
                    !      Otherwise add in contribution from upstream
                    !     **********************************************************************
                    SEDIMENT_YIELD(I_OLD,J_OLD)=SEDIMENT_YIELD(I_OLD,J_OLD) + SEDADD
                    IF (DO_FLOW_BOUNDARIES) THEN
                        IF (FLOW_DIRECTION(I_OLD,J_OLD) == 1) THEN
                            EXIT L720
                        ENDIF
                    ELSE
                        IF (IS_Y_PERIODIC) THEN
                            IF ((I_OLD == ILOWEST).AND.(J_OLD == JLOWEST)) EXIT L720
                        ELSE
                            IF (J_OLD == MY) EXIT L720
                        ENDIF
                    ENDIF
                ENDDO L720
            ENDDO M710
        ENDDO L710
        IF (COMPLETE_RUNOFF) LAKE_SURFACE_ELEVATION=LAKE_OUTLET_ELEVATION
        CALL IS_IT_SUBMERGED()
        LAKE_AREA=0.0
        DO K=1,TOTAL_BASINS
            IF (ENCLOSED(K)) THEN
                DO J=1,MY
                    DO I=1,MX
                        IF (BASIN_NUMBER(I,J) == K) THEN
                            IF (ELEVATION(I,J) <= LAKE_SURFACE_ELEVATION(K)) THEN
                                LAKE_AREA(K)=LAKE_AREA(K)+D2X
                            ENDIF
                        ENDIF
                    ENDDO
                ENDDO
            ENDIF
        ENDDO

        !      *********************************************************************
        !       Print out some info for debugging purposes
        !      *********************************************************************
        TSEDCONT=TSEDCONT*TIME_INCREMENT
        IF (DO_SEDIMENT_TRANSPORT) THEN
            IF(MOD((ITERATION+1),OUTPUT_PRINT_INTERVAL) == 0) THEN
                WRITE(OUTHIST,764) TSEDCONT
                764           FORMAT(      &
                ' BED EROSION CONTRIBUTION TO SEDIMENT YIELD=',G12.5)
                WRITE(OUTHIST,436)
                436         FORMAT(/,' SEDIMENT YIELD')
                CALL SUMMARIZE_MATRIX_DATA(SEDIMENT_YIELD,THEAVG,THEMAX,THEMIN)
            ENDIF
        ENDIF
        IF (COMPLETE_RUNOFF) THEN
            L8554: DO KK=1,TOTAL_BASINS
                IF (ENCLOSED(KK)) THEN
                    K=DOWNSTREAM_BASIN(KK)
                    KKK=K
                    IF (ENCLOSED(K)) THEN
                        II=I_OUTFLOW(K)
                        JJ=J_OUTFLOW(K)
                        I_OLD=II
                        J_OLD=JJ
                        LOCAL_DIRECTION=FLOW_DIRECTION(I_OLD,J_OLD)
                        IF (SUBMERGED(I_OLD,J_OLD).OR.(LOCAL_DIRECTION < 2)) CYCLE L8554
                        L8556:    DO
                            J_NEW=J_OLD+DOWNSTREAM(LOCAL_DIRECTION,2)
                            I_NEW=I_OLD+DOWNSTREAM(LOCAL_DIRECTION,1)
                            IF (IS_X_PERIODIC) THEN
                                IF (I_NEW < 1) I_NEW=MX
                                IF (I_NEW > MX) I_NEW=1
                            ENDIF
                            IF (IS_Y_PERIODIC) THEN
                                IF (J_NEW > MY) J_NEW=1
                                IF (J_NEW < 1) J_NEW=MY
                            ENDIF
                            LOCAL_DIRECTION=FLOW_DIRECTION(I_NEW,J_NEW)
                            IF (SUBMERGED(I_NEW,J_NEW).OR.(LOCAL_DIRECTION < 2) &
                            .OR.(BASIN_NUMBER(I_NEW,J_NEW) /= BASIN_NUMBER(I_OLD,J_OLD) )) &
                            THEN
                                K=BASIN_NUMBER(I_OLD,J_OLD)
                                IF (ENCLOSED(K)) THEN
                                    IF  ((II >= 1).AND.(II <= MX) &
                                    .AND.(I_OUTFLOW(DOWNSTREAM_BASIN(K)) >= 1) &
                                    .AND.(J_OUTFLOW(K) >= 1).AND.(I_OUTFLOW(K) <= MX).AND. &
                                    (J_OUTFLOW(K) <= MY).AND.(JJ >= 1).AND.(JJ <= MY)) THEN
                                        CALL DRAWLINE(II,JJ,I_OUTFLOW(K), &
                                        J_OUTFLOW(K),ROUTED_DISCHARGE(I_NEW,J_NEW))
                                        !                write(outhist,9054) k,downstream_basin(k),i_old,j_old,i_outflow(downstream_basin(k))&
                                        !                ,j_outflow(downstream_basin(k)),routed_discharge(i_outflow(k),j_outflow(k))
                                        9054 FORMAT('DL, K=',I5,' TO=',I5,' I_OLD=',I5,' J_OLD=',I5,' ITO=',I5, &
                                        ' JTO=',I5,' QE=',G12.5)
                                    ENDIF
                                ENDIF
                                EXIT L8556
                            ENDIF
                            I_OLD=I_NEW
                            J_OLD=J_NEW
                            IF (DO_FLOW_BOUNDARIES) THEN
                                IF (FLOW_DIRECTION(I_OLD,J_OLD) == 1) THEN
                                    EXIT L8556
                                ENDIF
                            ELSE
                                IF (IS_Y_PERIODIC) THEN
                                    IF ((I_OLD == ILOWEST).AND.(J_OLD == JLOWEST)) EXIT L8556
                                    IF (J_OLD == MY) EXIT L8556
                                ENDIF
                            ENDIF
                        ENDDO L8556
                    ENDIF
                ENDIF
            ENDDO L8554
        ENDIF
        999 RETURN
    END
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE BASIN_REPORT()
        USE ERODE_GLOBALS
        IMPLICIT NONE
        !      *******************************************************************
        !       Debugging report in case of improper agglomeration of drainage
        !       basins
        !      *******************************************************************
        INTEGER :: KKK,I,J,II,JJ
        REAL (8) :: TAREA,MAXAREA
        L934: DO KKK=1,TOTAL_BASINS
            TAREA=0
            MAXAREA=0
            DO  J=1,MY
                DO  I=1,MX
                    IF (BASIN_NUMBER(I,J) == KKK) THEN
                        TAREA=TAREA+D2X
                        IF (DRAINAGE_AREA(I,J) > MAXAREA) MAXAREA=DRAINAGE_AREA(I,J)
                    ENDIF
                ENDDO
            ENDDO
            IF(ENCLOSED(KKK)) THEN
                WRITE(OUTHIST,935)KKK,DOWNSTREAM_BASIN(KKK),I_OUTFLOW(KKK),J_OUTFLOW(KKK), &
                LAKE_OUTLET_ELEVATION(KKK),TAREA,MAXAREA,BASIN_DRAINAGE_AREA(KKK)
                935         FORMAT(I6,' ENCL., TO:',I5,' AT ',2I6,' ELEVATION=',E15.7,&
                ' AR=',G12.5, &
                ' MA=',G12.5,' LA=',G12.5)
                IF (DOWNSTREAM_BASIN(KKK) == 0) THEN
                    DO  II=1,MX
                        DO  JJ=1,MY
                            IF (BASIN_NUMBER(II,JJ) == KKK) THEN
                                WRITE(OUTHIST,4936) II,JJ
                                4936  FORMAT ('NULL AT I=',I5,' J=',I5)
                            ENDIF
                        ENDDO
                    ENDDO
                ENDIF
            ELSE
                IF(I_OUTFLOW(KKK) < 0) THEN
                    WRITE (OUTHIST,936) KKK
                    936             FORMAT(I6,' DECOMMISSIONED')
                ELSE
                    WRITE (OUTHIST,937) KKK,TAREA,MAXAREA
                    937             FORMAT(I6,' OUTLET, AR=',G12.5,' MA=',G12.5)
                ENDIF
            ENDIF
        ENDDO L934
        WRITE (OUTHIST, 942)
        942 FORMAT(' BASIN NUMBERS',/)
        DO  J=1,MY
            WRITE (OUTHIST,944)  J,(BASIN_NUMBER(I,J),I=1,MX)
            944     FORMAT(I5,10I5,12(/,'     ',10I5))
        ENDDO
        RETURN
    END
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


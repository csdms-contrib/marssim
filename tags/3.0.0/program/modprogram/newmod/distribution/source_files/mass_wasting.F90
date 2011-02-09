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
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE DO_MASS_WASTING()
        USE ERODE_GLOBALS
        IMPLICIT NONE
        !      ********************************************************************
        !        This subroutine models slope mass wasting and returns the net elevation
        !         change in ERODE_SLOPE(i,j)
        !   MODIFIES:  CFW, CFN, CFNE, CFNW, ERODE_SLOPE, IS_ROCK_SURFACE, NOMINAL_ERODE_SLOPE
        !              SLOPETOADD, SLGRAVTOADD
        !      ********************************************************************
        REAL (8) :: RAPID_CREEP
        EXTERNAL RAPID_CREEP
        INTEGER :: I,J,I_NEW,J_NEW,LOCAL_DIRECTION,I_OLD,J_OLD,IE,IW,IPATH,JS,JN
        REAL (8) :: GRADIENT,SUMSLOPE
        LOGICAL (1) :: REGROUTED(IMMX,JMMX)
        REAL (8) :: REGEMIN,REGEMAX,REGEAVG,REGENUM
        REAL (8) :: ROCEMIN,ROCEMAX,ROCEAVG,ROCENUM
        REAL (8) :: DMULT
        LOGICAL :: IREG, WRITEDETAIL
        WRITEDETAIL=.FALSE.
        ROCEMIN=1.0E+25
        REGEMIN=ROCEMIN
        ROCEMAX=-ROCEMIN
        REGEMAX=-ROCEMIN
        REGEAVG=0.0
        REGENUM=0.0
        ROCEAVG=0.0
        ROCENUM=0.0
        ROUTE_REGOLITH_OVER_ROCK=.TRUE.
        SLOPETOADD=0.0
        SLGRAVTOADD=0.0
        !      ********************************************************************
        !       This is the main loop determining rates of mass wasting
        !      ********************************************************************
        REGROUTED=.FALSE.
            DO  J = 2, MY
                DO  I = 2, MX-1
                    !      ********************************************************************
                    !       first determine mass wasting rate to the west, cfw(i,j)
                    !       initially use the absolute gradient
                    !      ********************************************************************
                    GRADIENT=ABS(ELEVATION(I,J)-ELEVATION(I-1,J))/CELL_SIZE
                    !      ********************************************************************
                    !       the mass wasting rate dependency on gradient is a function of the
                    !         input parameter powertouse
                    !      ********************************************************************
                    CFW(I,J)=SLOPE_DIFFUSIVITY*GRADIENT
                    !      ********************************************************************
                    !       if there is provision for a critical slope gradient then adjust
                    !         the mass wasting rate in accordance with critical_gradient_term, gradient,
                    !         slope_failure_diffusivity, and slope_gradient_exponent.  
                    !       Also make sure that the rate enhancement
                    !         does not exceed maximum_diffusivity_increase - this is done in subroutine failsafe
                    !      ********************************************************************
                    IF (USE_CRITICAL_SLOPE_CRADIENT) THEN
                        IF (USE_ROERING_MASS_WASTING) THEN
                            CFW(I,J)=CFW(I,J)*RAPID_CREEP(GRADIENT)
                        ELSE
                            CFW(I,J)=CFW(I,J)+RAPID_CREEP(GRADIENT)
                        ENDIF
                    ENDIF
                    !      ********************************************************************
                    !       Change the sign of cfw if the gradient is negative
                    !      ********************************************************************
                    IF (ELEVATION(I,J) < ELEVATION(I-1,J)) CFW(I,J)=-CFW(I,J)
                    !      ********************************************************************
                    !       Do the same for the other directions, including cfn, and if nine-
                    !          point divergence is used, cfne and cfnw -- also do it slightly
                    !         differently for matrix edges - taking into account that we have periodic
                    !         boundaries on the two lateral edges and a no-mass-flow boundary to the top
                    !      ********************************************************************
                    GRADIENT=ABS(ELEVATION(I,J)-ELEVATION(I+1,J-1))/CELL_SIZE
                    CFNE(I,J)=SLOPE_DIFFUSIVITY*GRADIENT
                    IF (USE_CRITICAL_SLOPE_CRADIENT) THEN
                        GRADIENT=GRADIENT*0.7071
                        IF (USE_ROERING_MASS_WASTING) THEN
                            CFNE(I,J)=CFNE(I,J)*RAPID_CREEP(GRADIENT)
                        ELSE
                            CFNE(I,J)=CFNE(I,J)+RAPID_CREEP(GRADIENT)
                        ENDIF
                    ENDIF
                    IF (ELEVATION(I,J) < ELEVATION(I+1,J-1)) CFNE(I,J)=-CFNE(I,J)
                    GRADIENT=ABS(ELEVATION(I,J)-ELEVATION(I,J-1))/CELL_SIZE
                    CFN(I,J)=SLOPE_DIFFUSIVITY*GRADIENT
                    IF (USE_CRITICAL_SLOPE_CRADIENT) THEN
                        IF (USE_ROERING_MASS_WASTING) THEN
                            CFN(I,J)=CFN(I,J)*RAPID_CREEP(GRADIENT)
                        ELSE
                            CFN(I,J)=CFN(I,J)+RAPID_CREEP(GRADIENT)
                        ENDIF
                    ENDIF
                    IF (ELEVATION(I,J) < ELEVATION(I,J-1)) CFN(I,J)=-CFN(I,J)
                    GRADIENT=ABS(ELEVATION(I,J)-ELEVATION(I-1,J-1))/CELL_SIZE
                    CFNW(I,J)=SLOPE_DIFFUSIVITY*GRADIENT
                    IF (USE_CRITICAL_SLOPE_CRADIENT) THEN
                        GRADIENT=GRADIENT*0.7071
                        IF (USE_ROERING_MASS_WASTING) THEN
                            CFNW(I,J)=CFNW(I,J)*RAPID_CREEP(GRADIENT)
                        ELSE
                            CFNW(I,J)=CFNW(I,J)+RAPID_CREEP(GRADIENT)
                        ENDIF
                    ENDIF
                    IF (ELEVATION(I,J) < ELEVATION(I-1,J-1)) CFNW(I,J)=-CFNW(I,J)
                ENDDO
            ENDDO
            DO  J = 2, MY
                I=1
                IF (IS_X_PERIODIC) THEN
                    IW=MX
                ELSE
                    IW=2
                ENDIF
                GRADIENT=ABS(ELEVATION(I,J)-ELEVATION(IW,J))/CELL_SIZE
                CFW(I,J)=SLOPE_DIFFUSIVITY*GRADIENT
                IF (USE_CRITICAL_SLOPE_CRADIENT) THEN
                    IF (USE_ROERING_MASS_WASTING) THEN
                        CFW(I,J)=CFW(I,J)*RAPID_CREEP(GRADIENT)
                    ELSE
                        CFW(I,J)=CFW(I,J)+RAPID_CREEP(GRADIENT)
                    ENDIF
                ENDIF
                IF (ELEVATION(I,J) < ELEVATION(IW,J)) CFW(I,J)=-CFW(I,J)
                GRADIENT=ABS(ELEVATION(I,J)-ELEVATION(I+1,J-1))/CELL_SIZE
                CFNE(I,J)=SLOPE_DIFFUSIVITY*GRADIENT
                IF (USE_CRITICAL_SLOPE_CRADIENT) THEN
                    GRADIENT=GRADIENT*0.7071
                    IF (USE_ROERING_MASS_WASTING) THEN
                        CFNE(I,J)=CFNE(I,J)*RAPID_CREEP(GRADIENT)
                    ELSE
                        CFNE(I,J)=CFNE(I,J)+RAPID_CREEP(GRADIENT)
                    ENDIF
                ENDIF
                IF (ELEVATION(I,J) < ELEVATION(I+1,J-1)) CFNE(I,J)=-CFNE(I,J)
                GRADIENT=ABS(ELEVATION(I,J)-ELEVATION(I,J-1))/CELL_SIZE
                CFN(I,J)=SLOPE_DIFFUSIVITY*GRADIENT
                IF (USE_CRITICAL_SLOPE_CRADIENT) THEN
                    IF (USE_ROERING_MASS_WASTING) THEN
                        CFN(I,J)=CFN(I,J)*RAPID_CREEP(GRADIENT)
                    ELSE
                        CFN(I,J)=CFN(I,J)+RAPID_CREEP(GRADIENT)
                    ENDIF
                ENDIF
                IF (ELEVATION(I,J) < ELEVATION(I,J-1)) CFN(I,J)=-CFN(I,J)
                GRADIENT=ABS(ELEVATION(I,J)-ELEVATION(IW,J-1))/CELL_SIZE
                CFNW(I,J)=SLOPE_DIFFUSIVITY*GRADIENT
                IF (USE_CRITICAL_SLOPE_CRADIENT) THEN
                    GRADIENT=GRADIENT*0.7071
                    IF (USE_ROERING_MASS_WASTING) THEN
                        CFNW(I,J)=CFNW(I,J)*RAPID_CREEP(GRADIENT)
                    ELSE
                        CFNW(I,J)=CFNW(I,J)+RAPID_CREEP(GRADIENT)
                    ENDIF
                ENDIF
                IF (ELEVATION(I,J) < ELEVATION(IW,J-1)) CFNW(I,J)=-CFNW(I,J)
                I=MX
                IF (IS_X_PERIODIC) THEN
                    IE=1
                ELSE
                    IE=MX-1
                ENDIF
                GRADIENT=ABS(ELEVATION(I,J)-ELEVATION(I-1,J))/CELL_SIZE
                CFW(I,J)=SLOPE_DIFFUSIVITY*GRADIENT
                IF (USE_CRITICAL_SLOPE_CRADIENT) THEN
                    IF (USE_ROERING_MASS_WASTING) THEN
                        CFW(I,J)=CFW(I,J)*RAPID_CREEP(GRADIENT)
                    ELSE
                        CFW(I,J)=CFW(I,J)+RAPID_CREEP(GRADIENT)
                    ENDIF
                ENDIF
                IF (ELEVATION(I,J) < ELEVATION(I-1,J)) CFW(I,J)=-CFW(I,J)
                GRADIENT=ABS(ELEVATION(I,J)-ELEVATION(IE,J-1))/CELL_SIZE
                CFNE(I,J)=SLOPE_DIFFUSIVITY*GRADIENT
                IF (USE_CRITICAL_SLOPE_CRADIENT) THEN
                    GRADIENT=GRADIENT*0.7071
                    IF (USE_ROERING_MASS_WASTING) THEN
                        CFNE(I,J)=CFNE(I,J)*RAPID_CREEP(GRADIENT)
                    ELSE
                        CFNE(I,J)=CFNE(I,J)+RAPID_CREEP(GRADIENT)
                    ENDIF
                ENDIF
                IF (ELEVATION(I,J) < ELEVATION(IE,J-1)) CFNE(I,J)=-CFNE(I,J)
                GRADIENT=ABS(ELEVATION(I,J)-ELEVATION(I,J-1))/CELL_SIZE
                CFN(I,J)=SLOPE_DIFFUSIVITY*GRADIENT
                IF (USE_CRITICAL_SLOPE_CRADIENT) THEN
                    IF (USE_ROERING_MASS_WASTING) THEN
                        CFN(I,J)=CFN(I,J)*RAPID_CREEP(GRADIENT)
                    ELSE
                        CFN(I,J)=CFN(I,J)+RAPID_CREEP(GRADIENT)
                    ENDIF
                ENDIF
                IF (ELEVATION(I,J) < ELEVATION(I,J-1)) CFN(I,J)=-CFN(I,J)
                GRADIENT=ABS(ELEVATION(I,J)-ELEVATION(I-1,J-1))/CELL_SIZE
                CFNW(I,J)=SLOPE_DIFFUSIVITY*GRADIENT
                IF (USE_CRITICAL_SLOPE_CRADIENT) THEN
                    GRADIENT=GRADIENT*0.7071
                    IF (USE_ROERING_MASS_WASTING) THEN
                        CFNW(I,J)=CFNW(I,J)*RAPID_CREEP(GRADIENT)
                    ELSE
                        CFNW(I,J)=CFNW(I,J)+RAPID_CREEP(GRADIENT)
                    ENDIF
                ENDIF
                IF (ELEVATION(I,J) < ELEVATION(I-1,J-1)) CFNW(I,J)=-CFNW(I,J)
            ENDDO
            J=1
            IF (IS_Y_PERIODIC) THEN
                JN=MY
            ELSE
                JN=2
            ENDIF
            DO  I = 2, MX-1
                GRADIENT=ABS(ELEVATION(I,J)-ELEVATION(I-1,J))/CELL_SIZE
                CFW(I,J)=SLOPE_DIFFUSIVITY*GRADIENT
                IF (USE_CRITICAL_SLOPE_CRADIENT) THEN
                    IF (USE_ROERING_MASS_WASTING) THEN
                        CFW(I,J)=CFW(I,J)*RAPID_CREEP(GRADIENT)
                    ELSE
                        CFW(I,J)=CFW(I,J)+RAPID_CREEP(GRADIENT)
                    ENDIF
                ENDIF
                IF (ELEVATION(I,J) < ELEVATION(I-1,J)) CFW(I,J)=-CFW(I,J)
                GRADIENT=ABS(ELEVATION(I,J)-ELEVATION(I+1,JN))/CELL_SIZE
                CFNE(I,J)=SLOPE_DIFFUSIVITY*GRADIENT
                IF (USE_CRITICAL_SLOPE_CRADIENT) THEN
                    GRADIENT=GRADIENT*0.7071
                    IF (USE_ROERING_MASS_WASTING) THEN
                        CFNE(I,J)=CFNE(I,J)*RAPID_CREEP(GRADIENT)
                    ELSE
                        CFNE(I,J)=CFNE(I,J)+RAPID_CREEP(GRADIENT)
                    ENDIF
                ENDIF
                IF (ELEVATION(I,J) < ELEVATION(I+1,JN)) CFNE(I,J)=-CFNE(I,J)
                GRADIENT=ABS(ELEVATION(I,J)-ELEVATION(I,JN))/CELL_SIZE
                CFN(I,J)=SLOPE_DIFFUSIVITY*GRADIENT
                IF (USE_CRITICAL_SLOPE_CRADIENT) THEN
                    IF (USE_ROERING_MASS_WASTING) THEN
                        CFN(I,J)=CFN(I,J)*RAPID_CREEP(GRADIENT)
                    ELSE
                        CFN(I,J)=CFN(I,J)+RAPID_CREEP(GRADIENT)
                    ENDIF
                ENDIF
                IF (ELEVATION(I,J) < ELEVATION(I,JN)) CFN(I,J)=-CFN(I,J)
                GRADIENT=ABS(ELEVATION(I,J)-ELEVATION(I-1,JN))/CELL_SIZE
                CFNW(I,J)=SLOPE_DIFFUSIVITY*GRADIENT
                IF (USE_CRITICAL_SLOPE_CRADIENT) THEN
                    GRADIENT=GRADIENT*0.7071
                    IF (USE_ROERING_MASS_WASTING) THEN
                        CFNW(I,J)=CFNW(I,J)*RAPID_CREEP(GRADIENT)
                    ELSE
                        CFNW(I,J)=CFNW(I,J)+RAPID_CREEP(GRADIENT)
                    ENDIF
                ENDIF
                IF (ELEVATION(I,J) < ELEVATION(I-1,JN)) CFNW(I,J)=-CFNW(I,J)
            ENDDO
            I=1
            IF (IS_X_PERIODIC) THEN
                IW=MX
            ELSE
                IW=2
            ENDIF
            GRADIENT=ABS(ELEVATION(I,J)-ELEVATION(IW,J))/CELL_SIZE
            CFW(I,J)=SLOPE_DIFFUSIVITY*GRADIENT
            IF (USE_CRITICAL_SLOPE_CRADIENT) THEN
                IF (USE_ROERING_MASS_WASTING) THEN
                    CFW(I,J)=CFW(I,J)*RAPID_CREEP(GRADIENT)
                ELSE
                    CFW(I,J)=CFW(I,J)+RAPID_CREEP(GRADIENT)
                ENDIF
            ENDIF
            IF (ELEVATION(I,J) < ELEVATION(IW,J)) CFW(I,J)=-CFW(I,J)
            GRADIENT=ABS(ELEVATION(I,J)-ELEVATION(I+1,JN))/CELL_SIZE
            CFNE(I,J)=SLOPE_DIFFUSIVITY*GRADIENT
            IF (USE_CRITICAL_SLOPE_CRADIENT) THEN
                GRADIENT=GRADIENT*0.7071
                IF (USE_ROERING_MASS_WASTING) THEN
                    CFNE(I,J)=CFNE(I,J)*RAPID_CREEP(GRADIENT)
                ELSE
                    CFNE(I,J)=CFNE(I,J)+RAPID_CREEP(GRADIENT)
                ENDIF
            ENDIF
            IF (ELEVATION(I,J) < ELEVATION(I+1,JN)) CFNE(I,J)=-CFNE(I,J)
            GRADIENT=ABS(ELEVATION(I,J)-ELEVATION(I,JN))/CELL_SIZE
            CFN(I,J)=SLOPE_DIFFUSIVITY*GRADIENT
            IF (USE_CRITICAL_SLOPE_CRADIENT) THEN
                IF (USE_ROERING_MASS_WASTING) THEN
                    CFN(I,J)=CFN(I,J)*RAPID_CREEP(GRADIENT)
                ELSE
                    CFN(I,J)=CFN(I,J)+RAPID_CREEP(GRADIENT)
                ENDIF
            ENDIF
            IF (ELEVATION(I,J) < ELEVATION(I,JN)) CFN(I,J)=-CFN(I,J)
            GRADIENT=ABS(ELEVATION(I,J)-ELEVATION(IW,JN))/CELL_SIZE
            CFNW(I,J)=SLOPE_DIFFUSIVITY*GRADIENT
            IF (USE_CRITICAL_SLOPE_CRADIENT) THEN
                GRADIENT=GRADIENT*0.7071
                IF (USE_ROERING_MASS_WASTING) THEN
                    CFNW(I,J)=CFNW(I,J)*RAPID_CREEP(GRADIENT)
                ELSE
                    CFNW(I,J)=CFNW(I,J)+RAPID_CREEP(GRADIENT)
                ENDIF
            ENDIF
            IF (ELEVATION(I,J) < ELEVATION(IW,JN)) CFNW(I,J)=-CFNW(I,J)
            I=MX
            IF (IS_X_PERIODIC) THEN
                IE=1
            ELSE
                IE=MX-1
            ENDIF
            GRADIENT=ABS(ELEVATION(I,J)-ELEVATION(I-1,J))/CELL_SIZE
            CFW(I,J)=SLOPE_DIFFUSIVITY*GRADIENT
            IF (USE_CRITICAL_SLOPE_CRADIENT) THEN
                IF (USE_ROERING_MASS_WASTING) THEN
                    CFW(I,J)=CFW(I,J)*RAPID_CREEP(GRADIENT)
                ELSE
                    CFW(I,J)=CFW(I,J)+RAPID_CREEP(GRADIENT)
                ENDIF
            ENDIF
            IF (ELEVATION(I,J) < ELEVATION(I-1,J)) CFW(I,J)=-CFW(I,J)
            GRADIENT=ABS(ELEVATION(I,J)-ELEVATION(IE,JN))/CELL_SIZE
            CFNE(I,J)=SLOPE_DIFFUSIVITY*GRADIENT
            IF (USE_CRITICAL_SLOPE_CRADIENT) THEN
                GRADIENT=GRADIENT*0.7071
                IF (USE_ROERING_MASS_WASTING) THEN
                    CFNE(I,J)=CFNE(I,J)*RAPID_CREEP(GRADIENT)
                ELSE
                    CFNE(I,J)=CFNE(I,J)+RAPID_CREEP(GRADIENT)
                ENDIF
            ENDIF
            IF (ELEVATION(I,J) < ELEVATION(IE,JN)) CFNE(I,J)=-CFNE(I,J)
            GRADIENT=ABS(ELEVATION(I,J)-ELEVATION(I,JN))/CELL_SIZE
            CFN(I,J)=SLOPE_DIFFUSIVITY*GRADIENT
            IF (USE_CRITICAL_SLOPE_CRADIENT) THEN
                IF (USE_ROERING_MASS_WASTING) THEN
                    CFN(I,J)=CFN(I,J)*RAPID_CREEP(GRADIENT)
                ELSE
                    CFN(I,J)=CFN(I,J)+RAPID_CREEP(GRADIENT)
                ENDIF
            ENDIF
            IF (ELEVATION(I,J) < ELEVATION(I,JN)) CFN(I,J)=-CFN(I,J)
            GRADIENT=ABS(ELEVATION(I,J)-ELEVATION(I-1,JN))/CELL_SIZE
            CFNW(I,J)=SLOPE_DIFFUSIVITY*GRADIENT
            IF (USE_CRITICAL_SLOPE_CRADIENT) THEN
                GRADIENT=GRADIENT*0.7071
                IF (USE_ROERING_MASS_WASTING) THEN
                    CFNW(I,J)=CFNW(I,J)*RAPID_CREEP(GRADIENT)
                ELSE
                    CFNW(I,J)=CFNW(I,J)+RAPID_CREEP(GRADIENT)
                ENDIF
            ENDIF
            IF (ELEVATION(I,J) < ELEVATION(I-1,JN)) CFNW(I,J)=-CFNW(I,J)
            IF (ROUTE_REGOLITH_OVER_ROCK)  THEN 
                !      ********************************************************************
                !       This loop
                !       routes mass wasting material downslope on steep bedrock slopes
                !      ********************************************************************
                L2086: DO  J=MYY,1,-1
                    M2086: DO  I=1,MX
                        !      ********************************************************************
                        !       We're only looking at steep bedrock slopes
                        !      ********************************************************************
                        IF (.NOT.IS_ROCK_SURFACE(I,J)) CYCLE M2086 !    goto 2086
                        !      ********************************************************************
                        !       If the direction of mass transport is towards the point i,j then
                        !       add in the contributions from the neighboring points, but not if the
                        !       surrounding point is bedrock
                        !      ********************************************************************
                        REGROUTED(I,J)=.TRUE.
                        IPATH=0
                        IW=I-1
                        IE=I+1
                        IF (IS_X_PERIODIC) THEN
                            IF (IW < 1) IW=MX
                            IF (IE > MX) IE=1
                        ENDIF
                        JN=J-1
                        JS=J+1
                        IF (IS_Y_PERIODIC) THEN
                            IF (JN < 1) JN=MY
                            IF (JS > MY) JS=1
                        ENDIF
                        IF ((J > 1).OR.IS_Y_PERIODIC) THEN
                            IF ((CFN(I,J) < 0.0).AND.(.NOT.IS_ROCK_SURFACE(I,JN))) &
                            ERODE_SLOPE(I,J)=ERODE_SLOPE(I,J)-CFN(I,J)
                        ENDIF
                        IF ((I > 1).OR.IS_X_PERIODIC) THEN
                            IF ((CFW(I,J) < 0.0).AND.(.NOT.IS_ROCK_SURFACE(IW,J))) &
                            ERODE_SLOPE(I,J)=ERODE_SLOPE(I,J)-CFW(I,J)
                        ENDIF
                        IF ((IS_Y_PERIODIC).OR.(J < MYY)) THEN
                            IF ((CFN(I,JS) > 0.0).AND.(.NOT.IS_ROCK_SURFACE(I,JS))) &
                            ERODE_SLOPE(I,J)=ERODE_SLOPE(I,J)+CFN(I,JS)
                        ENDIF
                        IF ((I < MX).OR.IS_X_PERIODIC) THEN
                            IF ((CFW(IE,J) > 0.0).AND.(.NOT.IS_ROCK_SURFACE(IE,J)))  &
                            ERODE_SLOPE(I,J)=ERODE_SLOPE(I,J)+CFW(IE,J)
                        ENDIF
                        ERODE_SLOPE(I,J)=ERODE_SLOPE(I,J)*CROSS_WEIGHTING
                        IF ((J > 1).OR.IS_Y_PERIODIC) THEN
                            IF ((I < MX).OR.IS_X_PERIODIC) THEN
                                IF ((CFNE(I,J) < 0.0).AND.(.NOT.IS_ROCK_SURFACE(IE,JN))) &
                                ERODE_SLOPE(I,J)=ERODE_SLOPE(I,J)-CFNE(I,J)*DIAGONAL_WEIGHTING
                            ENDIF
                            IF ((I > 1).OR.IS_X_PERIODIC) THEN
                                IF ((CFNW(I,J) < 0.0).AND.(.NOT.IS_ROCK_SURFACE(IW,JN))) &
                                ERODE_SLOPE(I,J)=ERODE_SLOPE(I,J)-CFNW(I,J)*DIAGONAL_WEIGHTING
                            ENDIF
                        ENDIF
                        IF ((IS_Y_PERIODIC).OR.(J < MYY)) THEN
                            IF ((I > 1).OR.IS_X_PERIODIC) THEN
                                IF ((CFNE(IW,JS) > 0.0).AND.(.NOT.IS_ROCK_SURFACE(IW,JS))) &
                                ERODE_SLOPE(I,J)=ERODE_SLOPE(I,J)+CFNE(IW,JS)*DIAGONAL_WEIGHTING
                            ENDIF
                            IF ((I < MX).OR.IS_X_PERIODIC) THEN
                                IF ((CFNW(IE,JS) > 0.0).AND.(.NOT.IS_ROCK_SURFACE(IE,JS))) &
                                ERODE_SLOPE(I,J)=ERODE_SLOPE(I,J)+CFNW(IE,JS)*DIAGONAL_WEIGHTING
                            ENDIF
                        ENDIF
                        ERODE_SLOPE(I,J)=ERODE_SLOPE(I,J)-REGOLITH(I,J)
                        !      ********************************************************************
                        !       We now need to route all mass wasting contributions from surrounding
                        !        points downslope to the nearest non-bedrock location.  We assume that
                        !        debris transport on bedrock slopes is "instantaneous"
                        !      ********************************************************************
                        I_OLD=I
                        J_OLD=J
                        L2087: DO
                            LOCAL_DIRECTION = FLOW_DIRECTION(I_OLD,J_OLD)
                            !      ********************************************************************
                            !       If we are in a depression we need to deposit, so branch
                            !      ********************************************************************
                            IF (LOCAL_DIRECTION < 2) THEN
                                IREG=.FALSE.
                                EXIT L2087
                            ENDIF
                            J_NEW=J_OLD+DOWNSTREAM(LOCAL_DIRECTION,2)
                            I_NEW=I_OLD+DOWNSTREAM(LOCAL_DIRECTION,1)
                            IF (IS_X_PERIODIC) THEN
                                IF (I_NEW < 1) I_NEW=MX
                                IF (I_NEW > MX) I_NEW=1
                            ENDIF
                            IPATH=IPATH+1
                            IF (IPATH > 300) THEN
                                ABORT_SIMULATION=.TRUE.
                            ENDIF
                            !      ********************************************************************
                            !       Forget it if we reach a non-eroding or border location
                            !      ********************************************************************
                            IF (IS_Y_PERIODIC) THEN
                                IF (J_NEW < 1) J_NEW=MY
                                IF (J_NEW > MY) J_NEW=1
                            ELSE
                                IF (J_NEW == MY) THEN
                                    ERODE_SLOPE(I,J)=REGOLITH(I,J)
                                    CYCLE M2086
                                ENDIF
                            ENDIF
                            !      ********************************************************************
                            !       Branch if we have reached a regolith-mantled slope
                            !      ********************************************************************
                            IREG=.FALSE.
                            IF ((.NOT.IS_ROCK_SURFACE(I_NEW,J_NEW))) THEN
                                IREG=.TRUE.
                                EXIT L2087
                            ENDIF
                            !      ********************************************************************
                            !       Otherwise continue marching downslope
                            !      ********************************************************************
                            I_OLD=I_NEW
                            J_OLD=J_NEW
                            IF (.NOT.ABORT_SIMULATION) CYCLE L2087
                            IF (ABORT_SIMULATION) THEN
                                IABORTMIN=I_OLD-3
                                IABORTMAX=I_OLD+3
                                JABORTMIN=J_OLD-3
                                JABORTMAX=J_OLD+3
                                WRITE(*,814) I_OLD,J_OLD,I,J
                                WRITE(OUTHIST,814) I_OLD,J_OLD,I,J
                                814         FORMAT(' *** INFINITE CYCLE AT ',2I6,' ***' &
                                ,/,' IN REGOLITH ROUTING STARTING FROM: ',2I6)
                                IF ((.NOT.IS_Y_PERIODIC).AND.(J_NEW == MY)) THEN
                                    CYCLE M2086
                                ELSE
                                    EXIT L2086
                                ENDIF
                            ENDIF
                        ENDDO L2087
                        IF (.NOT.IREG) THEN
                            I_NEW=I_OLD
                            J_NEW=J_OLD
                            IS_ROCK_SURFACE(I_NEW,J_NEW)=.FALSE.
                        ENDIF
                        !      ********************************************************************
                        !       We're at a regolith mantled slope, so we need to deposit our sediment,
                        !       which includes material mass wasted from sourrounding points and
                        !       the local contribution of weathered material
                        !      ********************************************************************
                        ERODE_SLOPE(I_NEW,J_NEW)=ERODE_SLOPE(I_NEW,J_NEW)+ERODE_SLOPE(I,J)
                        SLOPETOADD=SLOPETOADD+ERODE_SLOPE(I,J)*DABS(DFLOAT((I-I_NEW)**2) &
                        + DFLOAT((J-J_NEW)**2))
                        SLGRAVTOADD=SLGRAVTOADD+ERODE_SLOPE(I,J)*(ELEVATION(I,J)-ELEVATION(I_NEW,J_NEW))
                        !      ********************************************************************
                        !       Make sure that our local rate of erosion is zero for now
                        !      ********************************************************************
                        ERODE_SLOPE(I,J)=REGOLITH(I,J)
                    ENDDO M2086
                ENDDO L2086
            ENDIF
            !      ********************************************************************
            !       This is the main loop for regolith mantled slopes
            !      ********************************************************************
            L90: DO  J=MYY,1,-1
                M90: DO  I=1,MX
                    !      ********************************************************************
                    !       Include contributions from surrounding points if they are regolith
                    !         mantled and losses to surrounding points if they are either
                    !         regolith or bedrock
                    !      *******************************************************************
                    IE=I+1
                    IW=I-1
                    IF (IS_X_PERIODIC) THEN
                        IF (IE > MX) IE=1
                        IF (IW < 1) IW=MX
                    ENDIF
                    JN=J-1
                    JS=J+1
                    IF (IS_Y_PERIODIC) THEN
                        IF (JN < 1) JN=MY
                        IF (JS > MY) JS=1
                    ENDIF
                    !      ********************************************************************
                    !       Skip this if it is non-eroding border point or if it is bedrock
                    !      ********************************************************************
                    IF (.NOT.REGROUTED(I,J)) THEN
                        IF ((IE <= MX).OR.IS_X_PERIODIC) THEN
                            IF (CFW(IE,J) <= 0.0) THEN
                                ERODE_SLOPE(I,J)=ERODE_SLOPE(I,J)+CFW(IE,J)*CROSS_WEIGHTING
                            ELSE
                                IF (.NOT.REGROUTED(IE,J)) THEN
                                    ERODE_SLOPE(I,J)=ERODE_SLOPE(I,J)+CFW(IE,J)*CROSS_WEIGHTING
                                ENDIF
                            ENDIF
                        ENDIF
                        IF (CFN(I,JS) <= 0.0) THEN
                            ERODE_SLOPE(I,J)=ERODE_SLOPE(I,J)+CFN(I,JS)*CROSS_WEIGHTING
                        ELSE
                            IF (.NOT.REGROUTED(I,JS))  &
                            ERODE_SLOPE(I,J)=ERODE_SLOPE(I,J)+CFN(I,JS)*CROSS_WEIGHTING
                        ENDIF
                        IF (CFW(I,J) >= 0.0) THEN
                            ERODE_SLOPE(I,J)=ERODE_SLOPE(I,J)-CFW(I,J)*CROSS_WEIGHTING
                        ELSE
                            IF ((I > 1).OR.IS_X_PERIODIC) THEN
                                IF (.NOT.REGROUTED(IW,J)) THEN
                                    ERODE_SLOPE(I,J)=ERODE_SLOPE(I,J)-CFW(I,J)*CROSS_WEIGHTING
                                ENDIF
                            ENDIF
                        ENDIF
                        IF (CFN(I,J) >= 0.0) THEN
                            ERODE_SLOPE(I,J)=ERODE_SLOPE(I,J)-CFN(I,J)*CROSS_WEIGHTING
                        ELSE
                            IF ((J > 1).OR.(IS_Y_PERIODIC)) THEN
                                IF (.NOT.REGROUTED(I,JN)) THEN
                                    ERODE_SLOPE(I,J)=ERODE_SLOPE(I,J)-CFN(I,J)*CROSS_WEIGHTING
                                ENDIF
                            ENDIF
                        ENDIF
                        IF (CFNE(I,J) > 0.0) THEN
                            ERODE_SLOPE(I,J)=ERODE_SLOPE(I,J)-CFNE(I,J)*DIAGONAL_WEIGHTING
                        ELSE
                            IF ((J > 1).OR.(IS_Y_PERIODIC)) THEN
                                IF ((I < MX).OR.IS_X_PERIODIC) THEN
                                    IF (.NOT.REGROUTED(IE,JN)) THEN
                                        ERODE_SLOPE(I,J)=ERODE_SLOPE(I,J)-CFNE(I,J)*DIAGONAL_WEIGHTING
                                    ENDIF
                                ENDIF
                            ENDIF
                        ENDIF
                        IF (CFNW(I,J) > 0.0) THEN
                            ERODE_SLOPE(I,J)=ERODE_SLOPE(I,J)-CFNW(I,J)*DIAGONAL_WEIGHTING
                        ELSE
                            IF ((J > 1).OR.(IS_Y_PERIODIC)) THEN
                                IF ((I > 1).OR.IS_X_PERIODIC) THEN
                                    IF (.NOT.REGROUTED(IW,JN)) THEN
                                        ERODE_SLOPE(I,J)=ERODE_SLOPE(I,J)-CFNW(I,J)*DIAGONAL_WEIGHTING
                                    ENDIF
                                ENDIF
                            ENDIF
                        ENDIF
                        IF ((I > 1).OR.IS_X_PERIODIC) THEN
                            IF (CFNE(IW,JS) <= 0.0) THEN
                                ERODE_SLOPE(I,J)=ERODE_SLOPE(I,J)+CFNE(IW,JS)*DIAGONAL_WEIGHTING
                            ELSE
                                IF (.NOT.REGROUTED(IW,JS)) THEN
                                    ERODE_SLOPE(I,J)=ERODE_SLOPE(I,J)+CFNE(IW,JS)*DIAGONAL_WEIGHTING
                                ENDIF
                            ENDIF
                        ENDIF
                        IF ((I < MX).OR.IS_X_PERIODIC) THEN
                            IF (CFNW(IE,JS) <= 0.0) THEN
                                ERODE_SLOPE(I,J)=ERODE_SLOPE(I,J)+CFNW(IE,JS)*DIAGONAL_WEIGHTING
                            ELSE
                                IF (.NOT.REGROUTED(IE,JS)) THEN
                                    ERODE_SLOPE(I,J)=ERODE_SLOPE(I,J)+CFNW(IE,JS)*DIAGONAL_WEIGHTING
                                ENDIF
                            ENDIF
                        ENDIF
                        DMULT=CELL_SIZE
                        IF (USE_ROERING_MASS_WASTING) THEN
                            SLOPETOADD=SLOPETOADD+  &
                            (SLOPE_DIFFUSIVITY*D8_GRADIENT(I,J)*RAPID_CREEP(D8_GRADIENT(I,J)))*DMULT
                            SLGRAVTOADD=SLGRAVTOADD+  &
                            (SLOPE_DIFFUSIVITY*D8_GRADIENT(I,J)*RAPID_CREEP(D8_GRADIENT(I,J)))*DMULT*D8_GRADIENT(I,J)
                        ELSE
                            SLOPETOADD=SLOPETOADD+  &
                            (SLOPE_DIFFUSIVITY*D8_GRADIENT(I,J)+RAPID_CREEP(D8_GRADIENT(I,J)))*DMULT
                            SLGRAVTOADD=SLGRAVTOADD+  &
                            (SLOPE_DIFFUSIVITY*D8_GRADIENT(I,J)+RAPID_CREEP(D8_GRADIENT(I,J)))*DMULT*D8_GRADIENT(I,J)
                        ENDIF
                    ENDIF
                    NOMINAL_ERODE_SLOPE(I,J)=(CFN(I,JS)-CFN(I,J)-CFW(I,J))*CROSS_WEIGHTING  &
                    +(-CFNE(I,J)-CFNW(I,J))*DIAGONAL_WEIGHTING
                    IF ((I < MX).OR.IS_X_PERIODIC) THEN
                        NOMINAL_ERODE_SLOPE(I,J)=NOMINAL_ERODE_SLOPE(I,J)+CFW(IE,J)*CROSS_WEIGHTING &
                        + CFNW(IE,JS)*DIAGONAL_WEIGHTING
                    ENDIF
                    IF ((I > 1).OR.IS_X_PERIODIC) THEN
                        NOMINAL_ERODE_SLOPE(I,J)=NOMINAL_ERODE_SLOPE(I,J)+CFNE(IW,JS)*DIAGONAL_WEIGHTING
                    ENDIF
                ENDDO M90
            ENDDO L90
        DO  J=1,MY
            DO  I=1,MX
                SUMSLOPE=SUMSLOPE+ERODE_SLOPE(I,J)
            ENDDO
        ENDDO
        IF (WRITEDETAIL) WRITE(OUTHIST,4234) SUMSLOPE
        4234  FORMAT ('NET MASS WASTING ELEV CHANGE=',G12.5)
        DO J=1,MY
            DO I=1,MX
                IF (IS_ROCK_SURFACE(I,J))THEN
                    IF (ERODE_SLOPE(I,J) > ROCEMAX) ROCEMAX=ERODE_SLOPE(I,J)
                    IF (ERODE_SLOPE(I,J) < ROCEMIN) ROCEMIN=ERODE_SLOPE(I,J)
                    ROCEAVG=ROCEAVG+ERODE_SLOPE(I,J)
                    ROCENUM=ROCENUM+1.0
                ELSE
                    IF (ERODE_SLOPE(I,J) > REGEMAX) REGEMAX=ERODE_SLOPE(I,J)
                    IF (ERODE_SLOPE(I,J) < REGEMIN) REGEMIN=ERODE_SLOPE(I,J)
                    REGEAVG=REGEAVG+ERODE_SLOPE(I,J)
                    REGENUM=REGENUM+1.0
                ENDIF
            ENDDO
        ENDDO
        IF (ROCENUM > 0.0) ROCEAVG=ROCEAVG/ROCENUM
        IF (REGENUM > 0.0) REGEAVG=REGEAVG/REGENUM
        IF (WRITEDETAIL) WRITE(OUTHIST,14234) ROCEMIN,ROCEAVG,ROCEMAX,ROCENUM,REGEMIN, &
        REGEAVG,REGEMAX,REGENUM
        14234 FORMAT(' ROCK SLOPE EROSION: MIN=',G12.5,' AVG=',G12.5,' MAX=', &
        G12.5,' N=',G12.5,/, &
        ' REGOLITH SLOPE EROSION: MIN=',G12.5,' AVG=',G12.5,' MAX=', &
        G12.5,' N=',G12.5)
        RETURN
    END
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    REAL (8) FUNCTION RAPID_CREEP(GRADIENT)
        USE ERODE_GLOBALS
        !      ********************************************************************
        !        This is the look-up function for slope mass-wasting rate if there is a
        !         critical gradient for slope failure
        !      ********************************************************************
        IMPLICIT NONE
        REAL (8), INTENT(IN)  :: GRADIENT
        INTEGER ICAT
        IF (GRADIENT <= 0.0) THEN
            RAPID_CREEP=0.0
            RETURN
        ENDIF
        IF (GRADIENT >= GRADMAX) THEN
            RAPID_CREEP=FAILMAX
        ELSE
            ICAT=INT((100.0*GRADIENT)/GRADMAX)+1
            RAPID_CREEP= &
            RATECAT(ICAT)+(GRADIENT-GRADCAT(ICAT))*DELRATE(ICAT)
        ENDIF
        RETURN
    END
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    REAL (8) FUNCTION ROCK_MASS_WASTING(GRADIENT)
        USE ERODE_GLOBALS
        !      ********************************************************************
        !        This is the look-up function for bedrock mass-wasting rate if there is a
        !         critical gradient for slope failure
        !      ********************************************************************
        IMPLICIT NONE
        REAL (8), INTENT(IN) :: GRADIENT
        INTEGER ICAT
        IF (GRADIENT <= 0.0) THEN
           ROCK_MASS_WASTING=0.0
           RETURN
        ENDIF
        IF (GRADIENT >= RGRADMAX) THEN
            ROCK_MASS_WASTING=FAILMAX
        ELSE
            ICAT=INT((100.0*GRADIENT)/RGRADMAX)+1
            ROCK_MASS_WASTING= &
            RRATECAT(ICAT)+(GRADIENT-RGRADCAT(ICAT))*RDELRATE(ICAT)
        ENDIF
        RETURN
    END
    

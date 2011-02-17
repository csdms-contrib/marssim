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
    SUBROUTINE INITIALIZE_VARIABLES()
        USE ERODE_GLOBALS
        IMPLICIT NONE
        INTEGER :: I,J,K,IFLAG
        REAL (8) :: XXXX,XRATE,XINC,XGRAD,RAPID_CREEP,ROCK_MASS_WASTING,SUM
        EXTERNAL RAPID_CREEP
        EXTERNAL ROCK_MASS_WASTING
        !      *******************************************************************
        !       Set up the matrix im which determines relative x,y locations of
        !       surrounding points.  The matrix id indicates local drainage
        !       direction according to the following scheme, where 1 is present
        !       location and entries are cells towards which flow is directed
        !                         relative i
        !                         -1   0   1
        !                          ---------
        !                     -1 | 7   4   8  |
        !        relative j    0 | 3   1   5  | values of the matrix flow_direction(i,j)
        !                      1 | 6   2   9  |
        !
        !      flow_direction(i,j) is given a negative sign if the d8_gradient(i,j)<=0.0 (a depression)
        !
        !      flow_direction(i,j)=1 specifies an exit point
        !
        !       The matrix downstream(k,l) takes the value of flow_direction(i,j) as the k argument
        !       and gives relative i value of the downstream location for l=1 and relative j for l=2
        !  MODIFIES: DOWNSTREAM, WEIGHTS, MP, SQRTOFTWO, FORCE_SEDIMENT_CONSERVATION, PREVIOUS_TIME_INCREMENT
        !            CUMULEXCESS, DEPOSITWORK, ERODEWORK, ERODEGRAV, DEPOSITGRAV, CRATERGRAV
        !            CRATERWORK, LAVAWORK, LAVADEPTH, EOLIANDEPTH, ELAPSEDTIME, SLOPEWORK, SLOPEGRAV
        !            LASTSEDBIAS, SEDBIAS, SEDIMENT_DEPOSITED, CUMULATIVE_ELEVATION_CHANGE,
        !            CUMULATIVE_EOLIAN_CHANGE, CUMULATIVE_LAVA_CHANGE, CUMULATIVE_CRATERING_CHANGE
        !            FLEN, TIME_INCREMENT, ITERATION, JABORTMIN, IABORTMIN, IABORTMAX, JABORTMAX
        !            CROSS_WEIGHTING, DIAGONAL_WEIGHTING, ONE_SIXTH, MAXCRIT, NCRITS, D2X
        !            PREVIOUS_DISCHARGE_COEFFICIENT, GRADMAX,FAILMAX,GRADCAT,RATECAT,DELRATE
        !            RGRADCAT, RRATECAT, RDELRATE, RFAILMAX, STICKYFACTOR, STICKYPROB, STICKYCAT, STICKYDEL
        !            PRESENT_TIME, INUMLEVELS, OCEANNEXTTIME, PRESENT_TIME
        ! CALLS: FINDOCEAN_ELEVATION

        !      *******************************************************************
        DOWNSTREAM(1,1) =  0
        DOWNSTREAM(1,2) =  0
        DOWNSTREAM(2,1) =  0
        DOWNSTREAM(2,2) =  1
        DOWNSTREAM(3,1) = -1
        DOWNSTREAM(3,2) =  0
        DOWNSTREAM(4,1) =  0
        DOWNSTREAM(4,2) = -1
        DOWNSTREAM(5,1) =  1
        DOWNSTREAM(5,2) =  0
        DOWNSTREAM(6,1) = -1
        DOWNSTREAM(6,2) =  1
        DOWNSTREAM(7,1) = -1
        DOWNSTREAM(7,2) = -1
        DOWNSTREAM(8,1) =  1
        DOWNSTREAM(8,2) = -1
        DOWNSTREAM(9,1) =  1
        DOWNSTREAM(9,2) =  1
        MP = MX * MY
        SQRTOFTWO = SQRT(2.0)
        !       weights for printing shaded-relief images
        WEIGHTS(1,1,1)=4.0/4.0
        WEIGHTS(1,1,2)=0.0/4.0
        WEIGHTS(1,1,3)=0.0/4.0
        WEIGHTS(1,1,4)=0.0/4.0
        WEIGHTS(2,1,1)=2.0/4.0
        WEIGHTS(2,1,2)=2.0/4.0
        WEIGHTS(2,1,3)=0.0/4.0
        WEIGHTS(2,1,4)=0.0/4.0
        WEIGHTS(1,2,1)=2.0/4.0
        WEIGHTS(1,2,2)=0.0/4.0
        WEIGHTS(1,2,3)=2.0/4.0
        WEIGHTS(1,2,4)=0.0/4.0
        WEIGHTS(2,2,1)=1.0/4.0
        WEIGHTS(2,2,2)=1.0/4.0
        WEIGHTS(2,2,3)=1.0/4.0
        WEIGHTS(2,2,4)=1.0/4.0
        L1100: DO I=1, 2
            L1101: DO J= 1, 2
                SUM = 0.0
                L1110: DO K= 1,4
                    SUM = SUM+WEIGHTS(I,J,K)
                ENDDO L1110
                IF (SUM /= 1.0) THEN
                    WRITE(*,1120) I,J
                    1120              FORMAT(' BAD WEIGHTS AT I=',I5,' J=',I5)
                    STOP
                ENDIF
            ENDDO L1101
        ENDDO L1100
        FORCE_SEDIMENT_CONSERVATION=.FALSE.
        IF ((.NOT.IS_X_PERIODIC).OR.(.NOT.IS_Y_PERIODIC)) FORCE_SEDIMENT_CONSERVATION=.FALSE.
        IF (.NOT.DO_FLUVIAL_DETACHMENT) FORCE_SEDIMENT_CONSERVATION=.FALSE.
        IF (.NOT.DO_SEDIMENT_TRANSPORT) FORCE_SEDIMENT_CONSERVATION=.FALSE.
        IF (DO_ROCK_DEFORMATION) FORCE_SEDIMENT_CONSERVATION=.FALSE.
        IF (MODEL_EOLIAN_CHANGES) FORCE_SEDIMENT_CONSERVATION=.FALSE.
        IF (MODEL_LAVA_FLOWS) FORCE_SEDIMENT_CONSERVATION=.FALSE.
        IF (HAVE_INFLUENT_RIVERS) FORCE_SEDIMENT_CONSERVATION=.FALSE.
        IF (.NOT.DO_ALLUVIAL_SMOOTHING) FORCE_SEDIMENT_CONSERVATION=.FALSE.
        RSUFFIX='.RAW'
        RPREFIX='BSHADE'
        GRADCUT=0.0001
        XDIRECT='RESULTS'
        X1PREFIX='UNKN'
        X2PREFIX='DQS_'
        X3PREFIX='BDQS'
        IJDEBUG=0
        IF ( NEW_SIMULATION ) THEN
            IFILE1=48
            IFILE2=48
            IFILE3=48
            XFILE1=48
            XFILE2=48
            XFILE3=48
            RNUM1=CHAR(IFILE1)
            RNUM2=CHAR(IFILE2)
            RNUM3=CHAR(IFILE3)
            PREVIOUS_TIME_INCREMENT=0.0
            CUMULEXCESS=0.0
            DEPOSITWORK=0.0
            ERODEWORK=0.0
            ERODEGRAV=0.0
            DEPOSITGRAV=0.0
            CRATERGRAV=0.0
            CRATERWORK=0.0
            LAVAWORK=0.0
            LAVADEPTH=0.0
            EOLIANDEPTH=0.0
            ELAPSEDTIME=0.0
            SLOPEWORK=0.0
            SLOPEGRAV=0.0
            LASTSEDBIAS=0.0
            SEDBIAS=1.0
            IF (FLUVIAL_AND_SLOPE_MODELING) THEN
                SEDIMENT_DEPOSITED=.FALSE.
            ENDIF
            L355: DO  J=1,MYY
                L356: DO  I=1,MX
                    CUMULATIVE_ELEVATION_CHANGE(I,J)=0.0
                    CUMULATIVE_EOLIAN_CHANGE(I,J)=0.0
                    CUMULATIVE_LAVA_CHANGE(I,J)=0.0
                    CUMULATIVE_CRATERING_CHANGE(I,J)=0.0
                ENDDO L356
            ENDDO L355
        ENDIF
        !      *******************************************************************
        !       flen(k)is the distance to the surrounding matrix point for values
        !       of the argument k, which is the appropriate value of flow_direction(i,j)
        !      *******************************************************************
        L100: DO I=1,5
            FLEN(I)=1.0
        ENDDO L100
        L110: DO I=6,9
            FLEN(I)=1.0/SQRTOFTWO
        ENDDO L110
        NMAX = MAX(MX,MY)
        LMAX = NMAX*NMAX
        TIME_INCREMENT = 1.0E-5
        IF (.NOT.FLUVIAL_AND_SLOPE_MODELING)TIME_INCREMENT=1.0
        IF (NEW_SIMULATION) THEN
            ITERATION = 0
        ENDIF
        !      *******************************************************************
        !       Sets range of i,j for debugging printout
        !      *******************************************************************
        JABORTMIN=NMAX+1
        IABORTMIN=NMAX+1
        IABORTMAX=0
        JABORTMAX=0
        !      *******************************************************************
        !       cross_weighting, diagonal_weighting, and one_sixth are used in divergence
        !       calculation
        !      *******************************************************************
        CROSS_WEIGHTING = 2.0/(3.0*CELL_SIZE)
        DIAGONAL_WEIGHTING = 1.0/(6.0*CELL_SIZE)
        ONE_SIXTH = 1.0/6.0
        !      *******************************************************************
        !       Maxcrit is used in writnew.f90 to determine the channel network
        !      *******************************************************************
        MAXCRIT = 5.0
        !      *******************************************************************
        !       Ncrits is also used in writnew.f90
        !      *******************************************************************
        NCRITS = 50
        D2X=CELL_SIZE*CELL_SIZE
        PREVIOUS_DISCHARGE_COEFFICIENT=DISCHARGE_COEFFICIENT
        !      ********************************************************************
        !       this sets up a table-based look-up function (rapid_creep) for use when there is
        !         a critical value for slope failure - this avoids having to take powers
        !         at each iteration and each cell direction
        !      ********************************************************************
        IF (USE_CRITICAL_SLOPE_CRADIENT) THEN
            GRADMAX=CRITICAL_SLOPE_GRADIENT*(1.0-MAXIMUM_DIFFUSIVITY_INCREASE)**(1.0/SLOPE_GRADIENT_EXPONENT)
            FAILMAX=SLOPE_FAILURE_DIFFUSIVITY*(1.0/MAXIMUM_DIFFUSIVITY_INCREASE-1.0)
            WRITE(OUTHIST,135) GRADMAX,FAILMAX
            135 FORMAT(' FAILURE GRADIENT=',G12.5,' MAXIMUM RATE=',G12.5)
            L130: DO I=1,101
                GRADCAT(I)=(I-1)*GRADMAX/100.0
                IF (USE_ROERING_MASS_WASTING) THEN
                    RATECAT(I)=SLOPE_FAILURE_DIFFUSIVITY*(1.0/(1.0-CRITICAL_GRADIENT_TERM &
                    *GRADCAT(I)**SLOPE_GRADIENT_EXPONENT))
                ELSE
                    RATECAT(I)=SLOPE_FAILURE_DIFFUSIVITY*(1.0/(1.0-CRITICAL_GRADIENT_TERM &
                    *GRADCAT(I)**SLOPE_GRADIENT_EXPONENT)-1.0)
                ENDIF
                WRITE(OUTHIST,150)GRADCAT(I),RATECAT(I)
            ENDDO L130
            L120: DO  I=1,100
                DELRATE(I)=(RATECAT(I+1)-RATECAT(I))/(GRADCAT(I+1)-GRADCAT(I))
            ENDDO L120
            DELRATE(101)=0.0
            WRITE(OUTHIST,160)
            160 FORMAT(' TABLE OF GRADIENTS AND RATES')
            XXXX=1.5*GRADMAX
            XINC=XXXX/75.0
            XGRAD=0.0
            L140: DO I=1,75
                XRATE=RAPID_CREEP(XGRAD)
                XGRAD=XGRAD+XINC
                WRITE(OUTHIST,150) XGRAD,XRATE
                150 FORMAT(2G12.5)
            ENDDO L140
        ENDIF
        !      ********************************************************************
        !       This sets up a table-based look-up function (rapid_mass_wasting) for use when there is
        !         a critical value for bedrock slope failure - this avoids having to take powers
        !         at each iteration and each cell direction
        !      ********************************************************************
        IF (CRITICAL_BEDROCK_GRADIENT > 0.0) THEN
            XXXX=1.0/CRITICAL_BEDROCK_GRADIENT**SLOPE_GRADIENT_EXPONENT
            RGRADMAX=CRITICAL_BEDROCK_GRADIENT*(1.0-MAXIMUM_DIFFUSIVITY_INCREASE)**(1.0/SLOPE_GRADIENT_EXPONENT)
            FAILMAX=SLOPE_FAILURE_DIFFUSIVITY*(1.0/MAXIMUM_DIFFUSIVITY_INCREASE-1.0)
            WRITE(OUTHIST,148) RGRADMAX,FAILMAX
            148 FORMAT(' ROCK FAILURE GRADIENT=',G12.5,' MAXIMUM RATE=',G12.5)
            L180: DO I=1,101
                RGRADCAT(I)=(I-1)*RGRADMAX/100.0
                IF (USE_ROERING_MASS_WASTING)THEN
                    RRATECAT(I)=SLOPE_FAILURE_DIFFUSIVITY*(1.0/(1.0-XXXX &
                    *RGRADCAT(I)**SLOPE_GRADIENT_EXPONENT))
                ELSE
                    RRATECAT(I)=SLOPE_FAILURE_DIFFUSIVITY*(1.0/(1.0-XXXX &
                    *RGRADCAT(I)**SLOPE_GRADIENT_EXPONENT)-1.0)
                ENDIF
            ENDDO L180
            L190: DO I=1,100
                RDELRATE(I)=(RRATECAT(I+1)-RRATECAT(I))/ &
                (RGRADCAT(I+1)-RGRADCAT(I))
            ENDDO L190
            RDELRATE(101)=0.0
        ENDIF
        !     **********************************************************************
        !      Determines the lookup table for probability of drainage direction
        !         change on alluvial surfaces if sticky_sediment_routing is true.  Probabilities
        !         depend upon the parameter sticky_routing_critical_value
        !     **********************************************************************
        IF  (STICKY_SEDIMENT_ROUTING) THEN
            STICKYFACTOR = - DLOG(0.5D0)/(0.5*(STICKY_ROUTING_CRITICAL_VALUE-1.0))
            L201: DO I=1,11
                STICKYPROB(I)=1.0-EXP(-0.5*(I-1))
                STICKYCAT(I)=1.0+(I-1)/STICKYFACTOR
            ENDDO L201
            L202: DO I=1,10
                STICKYDEL(I)=(STICKYPROB(I+1)-STICKYPROB(I))/ &
                (STICKYCAT(I+1) -STICKYCAT(I))
            ENDDO L202
        ENDIF
        !     **********************************************************************
        !      If a variable ocean base level is used, this reads in the levels and times for recalculation
        !     **********************************************************************
        IF (VARIABLE_OCEAN_ELEVATION) THEN
            OPEN(67,FILE='OCEANLEVELS.DAT')
            I=1
            L6430: DO
                READ(67,*,IOSTAT=IFLAG) OCEAN_RECALCULATION_TIMES(I),OCEAN_LEVELS(I)
                IF (IFLAG < 0) EXIT L6430
                I=I+1
            ENDDO L6430
            !             goto 6430
            INUMLEVELS=I
            OCEANNEXTTIME=0.0
            WRITE(*,7440) INUMLEVELS
            7440     FORMAT(' NO OF OCEAN_RECALCULATION_TIMES=',I6)
            CLOSE(67)
            CALL FINDOCEAN_ELEVATION(PRESENT_TIME)
        ENDIF
        RETURN
    END

    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    REAL (8) FUNCTION NORMAL_RANDOM_DEVIATE()
        IMPLICIT NONE
        !     ********************************************************************
        !      Generates a normally disributed random deviate
        !     ********************************************************************
        REAL (8) :: V1, V2, SS
        REAL (4) :: TEMP1,TEMP2
        REAL (8) :: TEMP3, RRAND
        EXTERNAL RRAND
        L10: DO
            TEMP1=RRAND()
            TEMP2=RRAND()
            V1 = 2.0*TEMP1-1.0
            V2 = 2.0*TEMP2-1.0
            SS = V1*V1 + V2*V2
            IF (SS < 1.0) EXIT L10
        ENDDO L10
        TEMP3=-2.0*DLOG(SS)/SS
        IF (TEMP3 > 0.0) THEN
            SS = DSQRT(TEMP3)
        ELSE
            WRITE(*,150)
            150      FORMAT(' NEGATIVE ARGUMENT TO RANDOM DEVIATE')
            STOP
        ENDIF
        NORMAL_RANDOM_DEVIATE = V1*SS
        RETURN
    END
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    REAL (8) FUNCTION LOGNORMAL_RANDOM_DEVIATE1(NORMMEAN,NORMSD)
        IMPLICIT NONE
        !     ********************************************************************
        !       Generates a lognormal_random_deviately distributed random deviate
        !     ********************************************************************
        REAL (8), INTENT(IN) :: NORMSD,NORMMEAN
        REAL (8) ::  Y,NORMAL_RANDOM_DEVIATE
        !         write(*,200) normmean,normsd
        200   FORMAT(' NORMMEAN=',G12.5,' NORMSD=',G12.5)
        Y = NORMAL_RANDOM_DEVIATE()
        LOGNORMAL_RANDOM_DEVIATE1 = DEXP(NORMMEAN+NORMSD*Y)
        RETURN
    END
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    REAL (8) FUNCTION LOGNORMAL_RANDOM_DEVIATE(STDEV)
        IMPLICIT NONE
        !     ********************************************************************
        !       Generates a lognormally distributed random deviate
        !     ********************************************************************
        REAL (8), INTENT(IN) :: STDEV
        REAL (8) ::  Y,NORMAL_RANDOM_DEVIATE
        Y = NORMAL_RANDOM_DEVIATE()
        LOGNORMAL_RANDOM_DEVIATE = DEXP((1.0-0.5*STDEV)*Y*STDEV)
        RETURN
    END
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    REAL (8) FUNCTION EXPONENTIAL_DISTRIBUTION()
        IMPLICIT NONE
        !     ********************************************************************
        !       Generates a exponentially distributed random deviate
        !     ********************************************************************c
        REAL (8) :: STDEV
        REAL (4) :: Y 
        REAL (8) :: RRAND
        EXTERNAL RRAND
        !          write(*,300)
        300   FORMAT(' EXPONENTIAL_DISTRIBUTION')
        L100: DO
            Y=RRAND()
            IF (Y /= 0.0) EXIT L100
        ENDDO L100
        STDEV=Y
        EXPONENTIAL_DISTRIBUTION=-DLOG(STDEV)
        RETURN
    END
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    REAL(8) FUNCTION RRAND()
    REAL :: RANDNUMBER
        CALL RANDOM_NUMBER(RANDNUMBER)
        RRAND=RANDNUMBER
    RETURN
    END FUNCTION RRAND
    


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
    SUBROUTINE BOUNDARY_CONDITIONS()
        USE ERODE_GLOBALS
        IMPLICIT NONE
        INTEGER :: I,J
        REAL (8) :: UPLIFT
        !     **********************************************************************
        !       This subroutine applies boundary conditions at the end of the
        !         iteration.  The default is just a lowering of the
        !         lower matrix boundary.  But it could also include variable surface
        !         and subsurface uplift or depression through direct manipulation of
        !         elevation(i,j) and of the direct access resistance file or more indirectly
        !         through manipulation of the erosion_depth_index(i,j) file through, as shown below, a
        !         new matrix deformation(i,j)
        !   MODIFIES:  ELEVATION, SEDIMENT_BASE, DEFORMATION, NEWBASE
        !     **********************************************************************
        IF (.NOT.IS_Y_PERIODIC) THEN
            IF (MODEL_OCEAN_LEVEL) THEN
                DO I=1,MX
                    IF (ELEVATION(I,MY) > OCEAN_ELEVATION) THEN
                        ELEVATION(I,MY)=OCEAN_ELEVATION
                        IF (ELEVATION(I,MY) < SEDIMENT_BASE(I,MY)) SEDIMENT_BASE(I,MY)=ELEVATION(I,MY)
                        IS_ROCK_SURFACE(I,MY)=.FALSE.
                    ENDIF
                ENDDO
            ELSE
                IF (HORIZONTAL_LOWER_BOUNDARY) THEN
                    NEWBASE =TIME_INCREMENT * BOUNDARY_LOWERING_RATE
                    DO  I = 1 ,MX
                        ELEVATION(I,MY) = ELEVATION(I,MY) - NEWBASE
                    ENDDO
                ELSE
                    DO  I = 1, MX
                        IF (ERODING_LOWER_BOUNDARY(I)) THEN
                            ELEVATION(I,MY) = ELEVATION(I,MY) +TIME_INCREMENT * CFNW(I,MY-1)
                        ENDIF
                    ENDDO
                ENDIF
            ENDIF
        ENDIF
        !     **********************************************************************
        !      Here is an example of a spatially-varying uplift rate
        !          affecting both the surface elevation and the
        !         rock structural elevation
        !         to implement this I have:
        !           declared deformation as a matrix in erode.ins
        !           declared deformscale as a global variable
        !              ncluded it in the time-dependent
        !              list of simulation parameters included in the matrix simulation_parameters
        !              and read in from file inrates.dat
        !           Write out deformation as an output matrix and read it in again if
        !              the program is restarted (much as, say, sedbase is read in)
        !           Modified the subroutine finderodibility as shown below so as to
        !             have the entry into the 3-d resistance matrix affected by the
        !             cumulative deformation
        !     **********************************************************************
        IF (DO_ROCK_DEFORMATION) THEN
            L140:   DO  J=1,MY
                IF (J <= 40) THEN
                    UPLIFT=DEFORMSCALE *TIME_INCREMENT
                ELSE
                    !   Example of linearly varying deformation
                    IF (J < 50) THEN
                        UPLIFT=(5.0-J/10.0)*DEFORMSCALE *TIME_INCREMENT
                    ELSE
                        UPLIFT=0.0
                    ENDIF
                ENDIF
                L150:   DO  I=1,MX
                    !    Example of sinusoidal deformation
                    !                  uplift=dsin(aiparam*i+biparam)*dsin(ajparam*j+bjparam)
                    !         +          *time_increment * deformscale
                    143       FORMAT(2I5,3G12.5)
                    IF (USE_3D_SLOPE_RESISTANCE) THEN
                        DEFORMATION(I,J)=DEFORMATION(I,J)+UPLIFT
                    ENDIF
                    ELEVATION(I,J)=ELEVATION(I,J)+UPLIFT
                    IF (DO_SEDIMENT_TRANSPORT) SEDIMENT_BASE(I,J)=SEDIMENT_BASE(I,J)+UPLIFT
                ENDDO L150
            ENDDO L140
        ENDIF
        RETURN
    END
    !       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE DETERMINE_ERODIBILITY()
        USE ERODE_GLOBALS
        IMPLICIT NONE
        REAL (8) :: READ_ERODIBILITY
        !     ***********************************************************************
        !      Determines the erodibilty as a function of time and space.  The
        !      weatherability and the bedrock fluvial erodibilty are assumed to be
        !      proportional.  Creep rate and regolith fluvial erodibility are
        !      assumed not to vary with erodibility.
        !  MODIFIES: EROSION_DEPTH_INDEX, RELATIVE_RESISTANCE
        !  CALLS: READ_ERODIBILITY
        !     ***********************************************************************
        INTEGER :: I,J,K
        L100: DO  J=1,MY
            L101: DO  I=1,MX
                IF (IS_ROCK_SURFACE(I,J)) THEN
                    IF (DO_ROCK_DEFORMATION) THEN
                        K=((EREFERENCE-ELEVATION(I,J))+DEFORMATION(I,J))* &
                        VERTICAL_RESISTANCE_SCALING+1
                    ELSE
                        K=(EREFERENCE-ELEVATION(I,J))*VERTICAL_RESISTANCE_SCALING+1
                    ENDIF
                ELSE
                    IF (DO_ROCK_DEFORMATION) THEN
                        K=((EREFERENCE-ELEVATION(I,J))+DEFORMATION(I,J)&
                        - REGOLITH(I,J))*VERTICAL_RESISTANCE_SCALING+1
                    ELSE
                        K=(EREFERENCE-ELEVATION(I,J)-REGOLITH(I,J))*VERTICAL_RESISTANCE_SCALING+1
                    ENDIF
                ENDIF
                IF (K > MZ) K=MZ
                IF (K < 1) K=1
                IF (K /= EROSION_DEPTH_INDEX(I,J)) THEN
                    RELATIVE_RESISTANCE(I,J)=READ_ERODIBILITY(I,J,K)
                    EROSION_DEPTH_INDEX(I,J)=K
                ENDIF
            ENDDO L101
        ENDDO L100
        RETURN
    END
    !       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    REAL (8) FUNCTION READ_ERODIBILITY(I,J,K)
        !     ***********************************************************************
        !      Reads a value from the erodibility "cube" input file of random deviates
        !      and scales it to a lognormal distribution with mean unity and standard
        !      deviation diffusivity_variability
        !     ***********************************************************************
        USE ERODE_GLOBALS
        IMPLICIT NONE
        INTEGER I,J,K,IREC
        REAL (4) :: RRRR
        IREC=J+MY*(I-1)+MX*MY*(K-1)
        READ(INRESIST,REC=IREC) RRRR
        IF (SCALE_3D_ROCK_ERODIBILITY) THEN
            READ_ERODIBILITY=RANDMULT*DEXP(RRRR*SIGMANORM-SIGMASQ)
        ELSE
            READ_ERODIBILITY=RRRR
        ENDIF
        RETURN
    END
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE DETERMINE_EROSION_RATE()
        USE ERODE_GLOBALS
        IMPLICIT NONE
        !     **********************************************************************
        !      Eetermines the time-varying simulation parameters during each iteration
        !         first checks to see if we have exceeded the maximum time for the
        !         current erosion rate category.  if so, then set the new simulation
        !         parameters from the matrices erosion_rate_values and simulation_parameters
        !       %%%%%% needs to be updated with newer simumaltion parameters that
        !              might be temporally varying %%%%%%%%%
        !  MODIFIES: PARAMETER_CHANGE_INDEX
        !  CALLS: PRINT_RATE_STATISTICS
        !     **********************************************************************
        IF (PRESENT_TIME > TIMES_FOR_PARAMETER_CHANGES(PARAMETER_CHANGE_INDEX)) THEN
            IF (PARAMETER_CHANGE_INDEX > NUMBER_OF_PARAMETER_CHANGES) THEN
                 PARAMETER_CHANGE_INDEX=NUMBER_OF_PARAMETER_CHANGES
            ELSE
                WRITE(OUTHIST,100) PARAMETER_CHANGE_INDEX
                100  FORMAT(' NOW USING EROSION PARAMETERS FOR PARAMETER_CHANGE_INDEX=',I5)
                CALL PRINT_RATE_STATISTICS()
                BOUNDARY_LOWERING_RATE=EROSION_RATE_VALUES(PARAMETER_CHANGE_INDEX)
                SLOPE_DIFFUSIVITY=SIMULATION_PARAMETERS(PARAMETER_CHANGE_INDEX,1)
                BEDROCK_ERODIBILITY=SIMULATION_PARAMETERS(PARAMETER_CHANGE_INDEX,2)
                DISCHARGE_COEFFICIENT=SIMULATION_PARAMETERS(PARAMETER_CHANGE_INDEX,3)
                DETACHMENT_CRITICAL_SHEAR=SIMULATION_PARAMETERS(PARAMETER_CHANGE_INDEX,4)
                TRANSPORTFACTOR=SIMULATION_PARAMETERS(PARAMETER_CHANGE_INDEX,5)
                TRANSPORT_CRITICAL_DIM_SHEAR=SIMULATION_PARAMETERS(PARAMETER_CHANGE_INDEX,6)
                BISTABLE_CRITICAL_SHEAR=SIMULATION_PARAMETERS(PARAMETER_CHANGE_INDEX,7)
                LOW_EROSION_THRESHOLD=SIMULATION_PARAMETERS(PARAMETER_CHANGE_INDEX,8)
                HIGH_EROSION_THRESHOLD=SIMULATION_PARAMETERS(PARAMETER_CHANGE_INDEX,9)
                BISTABLE_RUNOFF_FACTOR=SIMULATION_PARAMETERS(PARAMETER_CHANGE_INDEX,10)
                DEFORMSCALE=SIMULATION_PARAMETERS(PARAMETER_CHANGE_INDEX,11)
                BISTABLE_BEDROCK_ERODIBILITY=SIMULATION_PARAMETERS(PARAMETER_CHANGE_INDEX,12)
                WRITE(OUTHIST,200) BOUNDARY_LOWERING_RATE,SLOPE_DIFFUSIVITY,BEDROCK_ERODIBILITY,DISCHARGE_COEFFICIENT &
                ,DETACHMENT_CRITICAL_SHEAR,TRANSPORTFACTOR,TRANSPORT_CRITICAL_DIM_SHEAR,BISTABLE_CRITICAL_SHEAR &
                ,LOW_EROSION_THRESHOLD,HIGH_EROSION_THRESHOLD,BISTABLE_RUNOFF_FACTOR,DEFORMSCALE,BISTABLE_BEDROCK_ERODIBILITY
                200 FORMAT(' BOUNDARY LOWERING RATE=',G12.5,/, &
                ' SLOPE DIFFUSIVITY=',G12.5,/, &
                ' BEDROCK ERODIBILITY=',G12.5,/, &                
                ' DISCHARGE COEFFICIENT=',G12.5,/, &                
                ' DETACHMENT CRITICAL SHEAR=',G12.5,/, &                
                ' TRANSPORT FACTOR=',G12.5,/, &                
                ' TRANSPORT CRITICAL DIMENSIONLESS SHEAR=',G12.5,/, &                
                ' BISTABLE CRITICAL SHEAR=',G12.5,/, &                
                ' LOW EROSION THRESHOLD=',G12.5,/, &                
                ' HIGH EROSION THRESHOLD=',G12.5,/, &                
                ' BISTABLE RUNOFF FACTOR=',G12.5,/, &                
                ' DEFORMSCALE=',G12.5,/, &                
                ' BISTABLE BEDROCK ERODIBILITY=',G12.5)
                DISCHARGE_COEFFICIENT=WATER_DENSITY*(MANNING/CHANNEL_WIDTH_CONSTANT)**0.6*GRAVITY * &
               (9.8/GRAVITY)**0.3
                PREVIOUS_DISCHARGE_COEFFICIENT=DISCHARGE_COEFFICIENT
                PREVIOUS_CRITICAL_SHEAR=DETACHMENT_CRITICAL_SHEAR
            ENDIF
            PARAMETER_CHANGE_INDEX=PARAMETER_CHANGE_INDEX+1
        ENDIF
        RETURN
    END

    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE MAKE_EVENT(EVENT_NUMBER)
        USE ERODE_GLOBALS
        IMPLICIT NONE
        REAL (8) :: ZPLANE,X1,X2,Y1,Y2
        REAL (8) :: AA,BB,DD,SS,SLINE,INTLINE,ZCOMP,RRAND
        INTEGER :: I,J,JCOMP
        EXTERNAL RRAND
        INTEGER, INTENT(IN) :: EVENT_NUMBER
        !     **********************************************************************
        !      Does user-defined modifications of the simulation at specified times
        !        during the simulations.  those times are read in from the 'events.prm'
        !        file.  this routine is called only once for each event
        !        the times can be accessed in the matrix event_times(event_number)
        !  MODIFIES: Whatever appropriate simulation variables
        !     **********************************************************************
        WRITE(*,100) EVENT_NUMBER
        WRITE(OUTHIST,100) EVENT_NUMBER
        100   FORMAT(' EVENT NUMBER ',I5,' HAS OCCURRED')
        !     An example: changing to universal runoff
    !    COMPLETE_RUNOFF=.TRUE.
    !    DO_ALLUVIAL_REROUTING=.FALSE.
    !    RETURN
    !END
    !     **********************************************************************
    !      some other possible event assigments:
    !           change erosion threshold:
    !               detachment_critical_shear=200.0
    !           add a resistant layer:
    !               resistant_surface_layer=.true.
    !               surface_layer_thickness=10.0
    !               surface_layer_resistance=10.0
    !               do j=1,my
    !               do i=1,mx
    !                   previous_elevation(i,j)=elevation(i,j)-surface_layer_thickness
    !               enddo
    !               enddo
    !           change the discharge to a lower, but more steady value
    !               discharge_constant=1.0e-05
    !               use_random_discharge=.false.
    !               flow_fraction=1.0
    !           increase downstream discharges relative to upstream values
    !               discharge_exponent=1.0
    !           reduce creep and rate of weathering
    !               slope_diffusivity=0.002
    !               rock_weathering_rate=0.0001
    !     **********************************************************************
    !        the following commented-out section has been used to model wave-cut erosion in a simulation
    !           of erosion of a coastal-plain landscape
    !     **********************************************************************

               ss=0.2
               if (event_number == 1) then
                  zplane=23.0
                  x1=cell_size*(mx/10)
                  x2=cell_size*mx
                  y1=cell_size*my
                  y2=cell_size*(my/20)
                else
                  zplane=8.0
                  x1=cell_size*(mx/2)
                  x2=cell_size*mx
                  y1=cell_size*my
                  y2=cell_size*(my/2)
                endif
                aa=sqrt((ss**2*(y1-y2)**2)/((y1-y2)**2+(x1-x2)**2))
                bb=sqrt(ss**2-aa**2)
                dd=-aa*x1-bb*y1-zplane
                sline=(y1-y2)/(x1-x2)
                intline=y1-sline*x1
                !write (*,200) aa,bb,dd,sline,intline
                200 format('aa=',g12.5,' bb=',g12.5,' dd=',g12.5,' sline=',g12.5, &
                ' inter=',g12.5)
                do j=1,my
                   do i=1,mx
                      jcomp=(sline*i*cell_size+intline)/cell_size
                      if (jcomp <= j) then
                         if (elevation(i,j) > zplane) then
                            elevation(i,j)=zplane+(rrand()-0.5)*0.1
                            is_rock_surface(i,j)=.false.
                            regolith(i,j)=initial_regolith_thickness
                            if (elevation(i,j) <= sediment_base(i,j)) then
                               sediment_base(i,j)=elevation(i,j)
                               is_sediment_covered(i,j)=.false.
                            endif
                         endif
                      else
                         zcomp=-aa*i*cell_size-bb*j*cell_size-dd
                         if (elevation(i,j) > zcomp) then
                            elevation(i,j)=zcomp+(rrand()-0.5)*0.1
                            is_rock_surface(i,j)=.false.
                            regolith(i,j)=initial_regolith_thickness
                            if (elevation(i,j) <= sediment_base(i,j)) then
                               sediment_base(i,j)=elevation(i,j)
                               is_sediment_covered(i,j)=.false.
                            endif
                         endif
                      endif
                   enddo
                enddo
               return
               end
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     SUBROUTINE FINDOCEAN_ELEVATION(TDUMMY)
        USE ERODE_GLOBALS
        IMPLICIT NONE
        REAL (8), INTENT(INOUT) :: TDUMMY
        REAL (8) :: EDUMMY,LASTOCEANLEVEL
        INTEGER :: I,J
        LOGICAL :: HASFOUND
    !  MODIFIES: TDUMMY, OCEAN_ELEVATION, PRESENT_TIME,MAXIMUMELEVATION, OCEANLSOPE, OCEANINT
    !            OCEANNEXTTIME
        IF (.NOT.VARIABLE_OCEAN_ELEVATION) RETURN
        LASTOCEANLEVEL=OCEAN_ELEVATION
        EDUMMY=-1.0E+25
        DO J=1,MY
            DO I=1,MX
                IF (ELEVATION(I,J) > EDUMMY) EDUMMY=ELEVATION(I,J)
            ENDDO
        ENDDO
        MAXIMUMELEVATION=EDUMMY
        L200: DO
            IF (TDUMMY >= OCEANNEXTTIME) THEN
                HASFOUND=.FALSE.
                DO I=1,INUMLEVELS-1
                    IF (TDUMMY <= OCEAN_RECALCULATION_TIMES(I+1)) THEN
                        HASFOUND=.TRUE.
                        EXIT
                    ENDIF
                ENDDO
                IF (HASFOUND) THEN
                    OCEANNEXTTIME=OCEAN_RECALCULATION_TIMES(I+1)
                    OCEANTSLOPE=(OCEAN_LEVELS(I+1)-OCEAN_LEVELS(I))/  &
                    (OCEAN_RECALCULATION_TIMES(I+1)-OCEAN_RECALCULATION_TIMES(I))
                    OCEANTINT=OCEAN_LEVELS(I)-OCEANTSLOPE*OCEAN_RECALCULATION_TIMES(I)
                ELSE
                    RETURN
                ENDIF

            ENDIF
            OCEAN_ELEVATION=OCEANTSLOPE*PRESENT_TIME+OCEANTINT
            IF (OCEAN_ELEVATION >= MAXIMUMELEVATION) THEN
                PRESENT_TIME=PRESENT_TIME+100.0
                TDUMMY=TDUMMY+100.0
                IF (PRESENT_TIME < MAXIMUM_SIMULATION_TIME) CYCLE L200
            ENDIF
            EXIT L200
        ENDDO L200
        RETURN
    END
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    LOGICAL (4) FUNCTION CHANGE_FLOW_DIRECTION(GRADIENT_RATIO)
        USE ERODE_GLOBALS
        IMPLICIT NONE
        REAL (8), INTENT(IN) :: GRADIENT_RATIO
        REAL (8) :: THEPROB, RRAND
        EXTERNAL RRAND
        INTEGER ICAT
        !     **********************************************************************
        !      This subroutine determines if the direction of flow on an alluvial
        !         surface remains the same or changes direction to a path with steeper
        !         gradient.  The decision is random with probability being an increasing
        !         function of the ratio of the new to the existing gradient and
        !         of the magnitude of the parameter sticky_routing_critical_value.  Uses a lookup
        !         function
        !     **********************************************************************
        IF (GRADIENT_RATIO <= 1.0) THEN
            CHANGE_FLOW_DIRECTION=.FALSE.
        ELSEIF (GRADIENT_RATIO > STICKYCAT(11)) THEN
            CHANGE_FLOW_DIRECTION=.TRUE.
        ELSE
            ICAT=INT((STICKYFACTOR*(GRADIENT_RATIO-1.0))+1.0)
            THEPROB=STICKYPROB(ICAT)+(GRADIENT_RATIO-STICKYCAT(ICAT))* &
            STICKYDEL(ICAT)
            IF (RRAND() > THEPROB) THEN
                CHANGE_FLOW_DIRECTION=.FALSE.
            ELSE
                CHANGE_FLOW_DIRECTION=.TRUE.
            ENDIF
        ENDIF
        RETURN
    END
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    



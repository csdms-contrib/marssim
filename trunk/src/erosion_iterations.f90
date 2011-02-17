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
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE SETUP_FLUVIAL_SLOPE_EROSION
        USE ERODE_GLOBALS
        IMPLICIT NONE
        REAL (8) :: LOGSD,NORMVAR,NORMMEAN,NORMSD,RRRR,LOGNORMAL_RANDOM_DEVIATE1,NORMAL_RANDOM_DEVIATE
        INTEGER :: INSTOP
        !      ********************************************************************
        !   MODIFIES: CFW, ABORT_SIMULATION,TEMPQCONSTANT,LASTQCONSTANT, DISCHARGE_SCALE_FACTOR
        !             EVAPORATION_SCALE, EVAPORATION_RATE
        !             
        !   CALLS: INITIALIZE_VARIABLES, CALCULATE_DIVERGENCE, EXPONENTIAL_HYDR_COND_GRNDWTR
        !          FIND_GROUNDWATER_FLUX, CONSTANT_HYDR_COND_GRNDWTR,DRAINAGE_BASIN_LAKE_FLOW
        !          DRAINAGE_BASIN_AREA_FLOW, DETERMINE_EROSION_RATE                 
        !
        !      ********************************************************************
        !      Set initial parameter values
        !      ********************************************************************
       
        CALL INITIALIZE_VARIABLES()
        !      ********************************************************************
        !      Abort_simulation is true if run is to be aborted
        !      ********************************************************************
        ABORT_SIMULATION = .FALSE.
        INSTOP=31
        !     **********************************************************************
        !      Determine topographic divergence
        !     **********************************************************************
        IF (FLUVIAL_AND_SLOPE_MODELING) CALL CALCULATE_DIVERGENCE()
        !      ********************************************************************
        !      Determine flow directions and gradients, and reset sed transport
        !         matrix
        !      ********************************************************************
        IF (FLUVIAL_AND_SLOPE_MODELING) THEN
            CFW=0.0
        ENDIF
        !      ********************************************************************
        !      Calculate groundwater flow
        !      ********************************************************************
        
        IF (MODEL_GROUNDWATER) THEN
            IF (EXPONENTIAL_PERMEABILITY_DECAY) THEN
                CALL EXPONENTIAL_HYDR_COND_GRNDWTR()
            ELSE
                CALL CONSTANT_HYDR_COND_GRNDWTR()
            ENDIF
            CALL FIND_GROUNDWATER_FLUX()
        ENDIF
        !      ********************************************************************
        !      Determine discharge scaling
        !      ********************************************************************
        
        IF (USE_RANDOM_DISCHARGE) THEN
            LOGSD=DISCHARGE_CONSTANT*DISCHARGE_COEFF_VARIATION*OCORRECT
            NORMVAR=DLOG((LOGSD**2+DISCHARGE_CONSTANT**2)/DISCHARGE_CONSTANT**2)
            NORMMEAN=DLOG(DISCHARGE_CONSTANT**2/DSQRT(DISCHARGE_CONSTANT**2+LOGSD**2))
            NORMSD=DSQRT(NORMVAR)
            WRITE(*,8211) NORMMEAN,NORMSD
            8211   FORMAT(' NORMMEAN=',G12.5,' NORMSD=',G12.5)
            RRRR=LOGNORMAL_RANDOM_DEVIATE1(NORMMEAN,NORMSD)
            TEMPQCONSTANT=LASTQCONSTANT*(1.0-OMEGA_WEIGHT)  &
            +OMEGA_WEIGHT*RRRR
            LASTQCONSTANT=TEMPQCONSTANT
        ELSE
            TEMPQCONSTANT=DISCHARGE_CONSTANT
            LASTQCONSTANT=DISCHARGE_CONSTANT
        ENDIF
        !      ********************************************************************
        !      Determine discharges and lake levels
        !      ********************************************************************
        DISCHARGE_SCALE_FACTOR=(TEMPQCONSTANT)**(1.0/DISCHARGE_EXPONENT)
        STEADY_DISCHARGE_FACTOR=DISCHARGE_CONSTANT/TEMPQCONSTANT
        CALL GRADIENT_AND_FLOW_DIRECTION()        
        IF (FLUVIAL_AND_SLOPE_MODELING) THEN
            IF (DO_SEDIMENT_TRANSPORT.OR.DO_FLUVIAL_DETACHMENT) THEN
                IF (MODEL_LAKE_EVAPORATION) THEN
                    EVAPORATION_SCALE=(EVAPORATION_MEAN+EVAPORATION_STANDARD_DEVIATION*NORMAL_RANDOM_DEVIATE())
                    IF (EVAPORATION_SCALE < 0.0) EVAPORATION_SCALE=0.0
                    WRITE(*,8298) EVAPORATION_SCALE
                    WRITE(OUTHIST,8298) EVAPORATION_SCALE
                    8298 FORMAT('EVAPORATION_RATE=',G12.5)
                    EVAPORATION_RATE=DISCHARGE_SCALE_FACTOR*EVAPORATION_SCALE
                    CALL DRAINAGE_BASIN_LAKE_FLOW()
                ELSE
                    CALL DRAINAGE_BASIN_AREA_FLOW()
                ENDIF
            ENDIF
        ENDIF
        !      ********************************************************************
        !      Stop simulation if something has gone wrong in subroutine drbasins
        !      ********************************************************************
        IF (ABORT_SIMULATION) THEN
            WRITE(OUTHIST,1777)
            WRITE(*,1777)
            1777    FORMAT(' ABORTING DUE TO DRBASINS ERROR')
            STOP
        ENDIF
        !     **********************************************************************
        !      Find present value of erosion rate and simulation parameters if they
        !         are time-varying
        !     **********************************************************************
        IF (VARIABLE_EROSION_RATE) CALL DETERMINE_EROSION_RATE()
        RETURN
    END
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE DO_FLUVIAL_AND_SLOPE()
        USE ERODE_GLOBALS
        IMPLICIT NONE
        INTEGER I,J,ILU,IUM
        REAL (8) :: LOGSD,NORMVAR,NORMMEAN,NORMSD,RRRR,LOGNORMAL_RANDOM_DEVIATE1,NORMAL_RANDOM_DEVIATE
        REAL (8) :: THEAVG,THEMIN,THEMAX
        INTEGER INSTOP,IDOSTOP
        INTEGER :: Reason

        !      ********************************************************************
        !      This routine is called each iteration to do fluvial and slope erosion
        !  MODIFIES:  ABORT_SIMULATION,TEMPQCONSTANT,LASTQCONSTANT,DISCHARGE_SCALE_FACTOR
        !             EVAPORATION_RATE, EVAPORATION_SCALE,, PREVIOUS_TIME_INCREMENT,CHANGECOUNT
        !              
        !  CALLS:  DETERMINE_ERODIBILITY, DETERMINE_EROSION_RATE, EXPONENTIAL_HYDR_COND_GRNDWTR
        !          CONSTANT_HYDR_COND_GRNDWTR, DRAINAGE_BASIN_LAKE_FLOW
        !          DRAINAGE_BASIN_AREA_FLOW, CALCULATE_DIVERGENCE, IS_IT_SUBMERGED
        !          PELAGIC_DEPOSIT, DO_THE_EROSION, PRINT_MORPHOMETRY, SUMMARIZE_CHANNELS
        !          CALCULATE_TOPO_CONVERGENCE, WRITE_FIRST_DATA_MATRIX, WRITE_ALLUVIAL_LOCATIONS
        !          WRITE_BEDROCK_LOCATIONS, WRITE_LAKE_INFO, WRITE_EROSION_DEPTH_INDEX
        !          WRITE_DEFORMATION, WRITE_SEDIMENT_BASE, WRITE_REGOLITH_THICKNESS
        !          WRITE_ROCK_RESISTANCE, PRINT_SIMULATION_INFORMATION, WRITE_ACCELERATED_EROSION_STATE
        !          WRITE_SHADED_RELIEF_IMAGE, WRITE_IMAGE, WRITE_GROUNDWATER_ELEVATION,
        !          WRITE_GROUNDWATER_FLOW, WRITE_SECOND_DATA_MATRIX

        !      ********************************************************************
        !      Abort_simulation is true if run is to be aborted
        !      ********************************************************************
        ABORT_SIMULATION = .FALSE.
        INSTOP=31
        !     **********************************************************************
        !      Find matrix of erodibilities if spatio-temporal variability of bedrock
        !             resistance is used
        !     **********************************************************************
        IF (USE_3D_SLOPE_RESISTANCE) THEN
            CALL DETERMINE_ERODIBILITY()
        ENDIF
        !     **********************************************************************
        !      Determine amount of bedrock weathering for both regolith-mantled and
        !             bedrock slopes
        !     **********************************************************************
        ICLAST=ICENT
        JCLAST=JCENT
        !     **********************************************************************
        !      Find present value of erosion rate and simulation parameters if they
        !         are time-varying
        !     **********************************************************************
        IF (VARIABLE_EROSION_RATE) CALL DETERMINE_EROSION_RATE()
        !      ********************************************************************
        !      Calculate groundwater flow
        !      ********************************************************************
        IF (MODEL_GROUNDWATER) THEN
            IF (MOD(ITERATION,SEEPAGE_ITERATION_INTERVAL) == 0) THEN
                IF (EXPONENTIAL_PERMEABILITY_DECAY) THEN
                    CALL EXPONENTIAL_HYDR_COND_GRNDWTR()
                ELSE
                    CALL CONSTANT_HYDR_COND_GRNDWTR()
                ENDIF
            ENDIF
        ENDIF
        !      ********************************************************************
        !      Determine discharge scaling
        !      ********************************************************************
        IF (MOD(ITERATION,RECALCULATE_GRADIENT_INTERVAL) == 0) THEN
            IF (USE_RANDOM_DISCHARGE) THEN
                LOGSD=DISCHARGE_CONSTANT*DISCHARGE_COEFF_VARIATION*OCORRECT
                NORMVAR=DLOG((LOGSD**2+DISCHARGE_CONSTANT**2)/DISCHARGE_CONSTANT**2)
                NORMMEAN=DLOG(DISCHARGE_CONSTANT**2/DSQRT(DISCHARGE_CONSTANT**2+LOGSD**2))
                NORMSD=DSQRT(NORMVAR)
                RRRR=LOGNORMAL_RANDOM_DEVIATE1(NORMMEAN,NORMSD)
                TEMPQCONSTANT=LASTQCONSTANT*(1.0-OMEGA_WEIGHT)  &
                +OMEGA_WEIGHT*RRRR
                LASTQCONSTANT=TEMPQCONSTANT
            ELSE
                TEMPQCONSTANT=DISCHARGE_CONSTANT
                LASTQCONSTANT=DISCHARGE_CONSTANT
            ENDIF
            DISCHARGE_SCALE_FACTOR=(TEMPQCONSTANT)**(1.0/DISCHARGE_EXPONENT)
            STEADY_DISCHARGE_FACTOR=DISCHARGE_CONSTANT/TEMPQCONSTANT
            4337      FORMAT('QS=',G12.5)
        !      ********************************************************************
        !      Determine discharges and lake levels
        !      ********************************************************************
            IF (FLUVIAL_AND_SLOPE_MODELING) THEN
                IF (DO_SEDIMENT_TRANSPORT.OR.DO_FLUVIAL_DETACHMENT) THEN
                    IF (MODEL_LAKE_EVAPORATION) THEN
                        IF (MOD(ITERATION,NCALCEVAP) == 0) THEN
                            IF(MOD(ITERATION,NCHANGEEVAP) == 0) THEN
                                EVAPORATION_SCALE=(EVAPORATION_MEAN+EVAPORATION_STANDARD_DEVIATION*NORMAL_RANDOM_DEVIATE())
                                IF (EVAPORATION_SCALE < 0.0) EVAPORATION_SCALE=0.0
                                WRITE(*,8298) EVAPORATION_SCALE
                                WRITE(OUTHIST,8298) EVAPORATION_SCALE
                                8298 FORMAT('EVAPORATION_RATE=',G12.5)
                            ENDIF
                            EVAPORATION_RATE=DISCHARGE_SCALE_FACTOR*EVAPORATION_SCALE
                            CALL DRAINAGE_BASIN_LAKE_FLOW()
                        ENDIF
                    ELSE
                        CALL DRAINAGE_BASIN_AREA_FLOW()
                    ENDIF
                ENDIF
            ELSE
                IF (DO_SEDIMENT_TRANSPORT.OR.DO_FLUVIAL_DETACHMENT) THEN
                    IF (MODEL_LAKE_EVAPORATION) THEN
                        IF(MOD(ITERATION,NCALCEVAP) == 0) THEN
                            IF(MOD(ITERATION,NCHANGEEVAP) == 0) THEN
                                EVAPORATION_SCALE=(EVAPORATION_MEAN+EVAPORATION_STANDARD_DEVIATION*NORMAL_RANDOM_DEVIATE())
                                IF (EVAPORATION_SCALE < 0.0) EVAPORATION_SCALE=0.0
                                WRITE(*,8298) EVAPORATION_SCALE
                                WRITE(OUTHIST,8298) EVAPORATION_SCALE
                            ENDIF
                            EVAPORATION_RATE=DISCHARGE_SCALE_FACTOR*EVAPORATION_SCALE
                            CALL DRAINAGE_BASIN_LAKE_FLOW()
                        ENDIF
                    ELSE
                        CALL DRAINAGE_BASIN_AREA_FLOW()
                    ENDIF
                ENDIF
            ENDIF
        ENDIF
        IF (MOD(ITERATION,DIVERGEINTERVAL) == 0) THEN
            IF (FLUVIAL_AND_SLOPE_MODELING) CALL CALCULATE_DIVERGENCE()
        ENDIF
        !      ********************************************************************
        !      Do the actual erosion
        !      ********************************************************************
        IF (MODEL_PELAGIC_DEPOSITION) PELAGIC_SEDIMENT_VOLUME=0.0
        CALL DO_THE_EROSION()
        IF (MODEL_PELAGIC_DEPOSITION) THEN
            CALL IS_IT_SUBMERGED()
            CALL PELAGIC_DEPOSIT()
        ENDIF
        PREVIOUS_TIME_INCREMENT=TIME_INCREMENT
        IF (.NOT.ABORT_SIMULATION) THEN
            !      ********************************************************************
            !      If it is time to calculate basin statistics do so
            !      ********************************************************************
            IF (MOD(ITERATION,OUTPUT_PRINT_INTERVAL) == 0) THEN
                IF (MOD(ITERATION,ELEVATION_PRINT_INTERVAL) /= 0) THEN
                    IF (DO_MORPHOMETRY) THEN
                        CALL PRINT_MORPHOMETRY()
                        CALL CALCULATE_TOPO_DIVERGENCE()
                        WRITE_OUT_CHANNELS=.FALSE.
                        CALL SUMMARIZE_CHANNELS()
                    ENDIF
                ENDIF
            ENDIF
            !      ********************************************************************
            !      Initialize simulation summary variables
            !      ********************************************************************
            IF (ITERATION < 2) THEN
                NUMBIDCHANGE=0.0
                SUMGRADCHANGE=0.0
                ABSGRADCHANGE=0.0
                CHANGECOUNT=0.0
            ENDIF
            !      ********************************************************************
            !      If it is time to summarize basin changes then do so
            !      ********************************************************************
            IF(MOD(ITERATION,WRITE_CHANGE_INTERVAL) == 0) THEN
                IF (CHANGECOUNT > 0.0) THEN
                    SUMGRADCHANGE=SUMGRADCHANGE/  &
                    (CHANGECOUNT*(MX-1)*(MY-1))
                    ABSGRADCHANGE=ABSGRADCHANGE/  &
                    (CHANGECOUNT*(MX-1)*(MY-1))
                    NUMBIDCHANGE=NUMBIDCHANGE/(CHANGECOUNT*(MX-1)*(MY-1))
                    WRITE(OUTRECORD,580) PRESENT_TIME,SUMGRADCHANGE, &
                    ABSGRADCHANGE, &
                    NUMBIDCHANGE
                    580             FORMAT(4(' ',E15.7))
                ENDIF
                NUMBIDCHANGE=0.0
                SUMGRADCHANGE=0.0
                ABSGRADCHANGE=0.0
                CHANGECOUNT=0.0
            ENDIF
            !      ********************************************************************
            !      Changecount is number of iterations over which changes are
            !        summarized
            !      ********************************************************************
            CHANGECOUNT=CHANGECOUNT+1.0
            !      ********************************************************************
            !      Check to see if operator intervention is requesting a stop
            !      ********************************************************************
            IF(MOD(ITERATION,10) == 0) THEN
                OPEN(INSTOP,FILE=trim(INPUT_DIRECTORY_PATH) // 'ERODE.STOP')
                READ(INSTOP,*,IOSTAT=Reason) IDOSTOP
                !Job running on the code were having issues on reading the ERODE.STOP
                !file for each iteration.IOSTAT catches the error and continues the execution
                CLOSE(INSTOP)
                IF (Reason == 0)  THEN
	                IF (IDOSTOP > 0) THEN
	                	PRINT *, "IN IF"
	                    IF(MOD(ITERATION,WRITE_CHANGE_INTERVAL) /= 0) THEN
	                        IF (CHANGECOUNT > 0.0) THEN
	                            SUMGRADCHANGE=SUMGRADCHANGE/  &
	                            (CHANGECOUNT*(MX-1)*(MY-1))
	                            ABSGRADCHANGE=ABSGRADCHANGE/  &
	                            (CHANGECOUNT*(MX-1)*(MY-1))
	                            NUMBIDCHANGE=NUMBIDCHANGE/(CHANGECOUNT* &
	                            (MX-1)*(MY-1))
	                            WRITE(OUTRECORD,580) PRESENT_TIME,SUMGRADCHANGE, &
	                            ABSGRADCHANGE,  &
	                            NUMBIDCHANGE
	                        ENDIF
	                    ENDIF
	                    IF(MOD(ITERATION,WRITE_CHANGE_INTERVAL) /= 0) THEN
	                        CALL WRITE_REPORT()
	                    ENDIF
	                    IF (MOD(ITERATION,OUTPUT_PRINT_INTERVAL) /= 0) THEN
	                        IF (MOD(ITERATION,ELEVATION_PRINT_INTERVAL) /= 0) THEN
	                            IF (FLUVIAL_AND_SLOPE_MODELING) CALL CHANNEL_PROPERTIES()
	                            CALL WRITE_FIRST_DATA_MATRIX()
	                            CALL WRITE_ALLUVIAL_LOCATIONS(OUTALLUV)
	                            CALL WRITE_BEDROCK_LOCATIONS(OUTROCK)
	                            CALL WRITE_LAKE_INFO(OUTLAKE)
	                            CALL WRITE_EROSION_DEPTH_INDEX(OUTEROSION_DEPTH_INDEX)
	                            CALL WRITE_DEFORMATION(OUTDEFORM)
	                            CALL WRITE_SEDIMENT_BASE()
	                            CALL WRITE_REGOLITH_THICKNESS()
	                            CALL WRITE_ROCK_RESISTANCE()
	                            CALL PRINT_SIMULATION_INFORMATION()
	                            CALL WRITE_ACCELERATED_EROSION_STATE(OUTHIGH)
	                            IF (MODEL_GROUNDWATER) THEN
	                                CALL WRITE_GROUNDWATER_ELEVATION(IOTEMP1)
	                                CALL WRITE_GROUNDWATER_FLOW(IOTEMP1)
	                            ENDIF
	                            IF (DO_MORPHOMETRY) THEN
	                                CALL PRINT_MORPHOMETRY()
	                                CALL CALCULATE_TOPO_DIVERGENCE()
	                                WRITE_OUT_CHANNELS=.TRUE.
	                                CALL SUMMARIZE_CHANNELS()
	                                CALL WRITE_SECOND_DATA_MATRIX()
	                            ENDIF
	                            IF (MOD(ITERATION,IMAGE_OUTPUT_INTERVAL) /= 0) THEN
	                                CALL WRITE_IMAGE()
	                                CALL WRITE_SHADED_RELIEF_IMAGE()
	                            ENDIF
	                        ENDIF
	                    ENDIF
	                    CLOSE(OUTHIST)
	                    CLOSE(OUTCHAN)
	                    CLOSE(OUTRECORD)
	                    CLOSE(OUTSUMMARY)
	                    CLOSE(OUTSOURCE)
	                    CLOSE(OUTREPORT)
	                    CLOSE(IOTEMP1)
	                    CLOSE(IOTEMP2)
	                    CLOSE(OUTCRATER)
	                    STOP
	                ENDIF
	        	ENDIF
            ENDIF
            !      ********************************************************************
            !      Recycle for another iteration if everything is ok and we are
            !           not done
            !      ********************************************************************
            IF ((MAXGRADIENT < 10000.0*CELL_SIZE) &
            .AND.(.NOT.ABORT_SIMULATION)) RETURN
        ENDIF
        IF (ABORT_SIMULATION) THEN
            !      ********************************************************************
            !      Something is wrong - print out summary information and abort
            !      ********************************************************************
            WRITE(OUTHIST,501)
            501     FORMAT('*** ABORTED DUE TO SINKHOLES ***')
            IF (FLUVIAL_AND_SLOPE_MODELING) CALL CHANNEL_PROPERTIES()
            CALL WRITE_FIRST_DATA_MATRIX()
            CALL WRITE_ALLUVIAL_LOCATIONS(OUTALLUV)
            CALL WRITE_BEDROCK_LOCATIONS(OUTROCK)
            CALL WRITE_LAKE_INFO(OUTLAKE)
            CALL WRITE_EROSION_DEPTH_INDEX(OUTEROSION_DEPTH_INDEX)
            CALL WRITE_DEFORMATION(OUTDEFORM)
            CALL WRITE_SEDIMENT_BASE()
            CALL WRITE_REGOLITH_THICKNESS()
            CALL WRITE_ROCK_RESISTANCE()
            CALL PRINT_SIMULATION_INFORMATION()
            CALL WRITE_ACCELERATED_EROSION_STATE(OUTHIGH)
            CALL WRITE_SHADED_RELIEF_IMAGE()
            CALL WRITE_IMAGE()
            IF (MODEL_GROUNDWATER) THEN
                CALL WRITE_GROUNDWATER_ELEVATION(IOTEMP1)
                CALL WRITE_GROUNDWATER_FLOW(IOTEMP1)
            ENDIF
            IF (DO_MORPHOMETRY) THEN
                CALL PRINT_MORPHOMETRY()
                WRITE_OUT_CHANNELS=.TRUE.
                CALL CALCULATE_TOPO_DIVERGENCE()
                CALL SUMMARIZE_CHANNELS()
                CALL WRITE_SECOND_DATA_MATRIX()
            ENDIF
            ILU=-1
            IUM=MX+1
            IABORTMIN = IABORTMIN-1
            IF (IABORTMIN < 1) THEN
                IUM=MX-(1-IABORTMIN)
                IABORTMIN=1
            ENDIF
            IABORTMAX = IABORTMAX+1
            IF (IABORTMAX  >  MX) THEN
                ILU=(IABORTMAX-MX)
                IABORTMAX=MX
            ENDIF
            JABORTMIN = JABORTMIN-1
            IF (JABORTMIN < 1) JABORTMIN=1
            JABORTMAX = JABORTMAX+1
            IF (JABORTMAX  >  MY) JABORTMAX =MY
            WRITE(OUTHIST,702)
            702     FORMAT('*****  AREAS ****')
            !             write(outhist,703)
            703     FORMAT('I=     ')
            WRITE(OUTHIST,704)(I,I=IUM,MX), &
            (I,I=IABORTMIN,IABORTMAX),(I,I=1,ILU)
            704     FORMAT(' I=    ',20I5,//)
            DO  J= JABORTMIN, JABORTMAX
                WRITE(OUTHIST,705)J,(DRAINAGE_AREA(I,J),I=IUM,MX),  &
                (DRAINAGE_AREA(I,J),I=IABORTMIN,IABORTMAX),(DRAINAGE_AREA(I,J),I=1,ILU)
                705         FORMAT(20I5)
            ENDDO
            WRITE(OUTHIST,706)
            706     FORMAT(//)
            WRITE(OUTHIST,712)
            712     FORMAT('*****  BASINS ****')
            !             write(outhist,713)
            713     FORMAT('I=     ')
            WRITE(OUTHIST,714)(I,I=IUM,MX),&
            (I,I=IABORTMIN,IABORTMAX),(I,I=1,ILU)
            714     FORMAT(' I=    ',20I5,//)
            DO  J= JABORTMIN, JABORTMAX
                WRITE(OUTHIST,715)J,(BASIN_NUMBER(I,J),I=IUM,MX), &
                (BASIN_NUMBER(I,J),I=IABORTMIN,IABORTMAX), &
                (BASIN_NUMBER(I,J),I=1,ILU)
                715         FORMAT(20I5)
            ENDDO
            WRITE(OUTHIST,716)
            716     FORMAT(//)
            WRITE(OUTHIST,502)
            502     FORMAT('*****  DIRECTIONS ****')
            503     FORMAT('I=     ')
            WRITE(OUTHIST,504)(I,I=IUM,MX), &
            (I,I=IABORTMIN,IABORTMAX),(I,I=1,ILU)
            504     FORMAT(' I=    ',20I5,//)
            DO  J= JABORTMIN, JABORTMAX
                WRITE(OUTHIST,505)J,(FLOW_DIRECTION(I,J),I=IUM,MX), &
                (FLOW_DIRECTION(I,J),I=IABORTMIN,IABORTMAX),(FLOW_DIRECTION(I,J),I=1,ILU)
                505         FORMAT(20I5)
            ENDDO
            WRITE(OUTHIST,506)
            506     FORMAT(//)
            WRITE(OUTHIST,507)
            507     FORMAT ('*****  ELEVATIONS ****')
            WRITE(OUTHIST,508)(I,I=IUM,MX), &
            (I,I=IABORTMIN,IABORTMAX),(I,I=1,ILU)
            508     FORMAT(' I= ',10I8,//)
            DO  J= JABORTMIN, JABORTMAX
                WRITE(OUTHIST,509)J,(ELEVATION(I,J),I=IUM,MX), &
                (ELEVATION(I,J),I=IABORTMIN,IABORTMAX),(ELEVATION(I,J),I=1,ILU)
                509         FORMAT(I4,6(G12.4),/,3(G12.4))
            ENDDO
            WRITE(OUTHIST,506)
            WRITE(OUTHIST,599)
            599     FORMAT ('*****  GRADIENTS ****')
            WRITE(OUTHIST,503)
            WRITE(OUTHIST,508)(I,I=IUM,MX),    &
            (I,I=IABORTMIN,IABORTMAX),(I,I=1,ILU)
            DO  J= JABORTMIN, JABORTMAX
                WRITE(OUTHIST,509)J,(D8_GRADIENT(I,J),I=IUM,MX), &
                (D8_GRADIENT(I,J),I=IABORTMIN,IABORTMAX),(D8_GRADIENT(I,J),I=1,ILU)
            ENDDO
            WRITE(OUTHIST,506)
            WRITE(OUTHIST,587)
            587     FORMAT ('*****  RESISTANCE ****')
            WRITE(OUTHIST,588)(I,I=IUM,MX), &
            (I,I=IABORTMIN,IABORTMAX),(I,I=1,ILU)
            588     FORMAT(' I= ',10I8,//)
            DO  J= JABORTMIN, JABORTMAX
                WRITE(OUTHIST,589)J,(RELATIVE_RESISTANCE(I,J),I=IUM,MX), &
                (RELATIVE_RESISTANCE(I,J),I=IABORTMIN,IABORTMAX),&
                (RELATIVE_RESISTANCE(I,J),I=1,ILU)
                589         FORMAT(I4,6(G12.4),/,3(G12.4))
            ENDDO
            WRITE(OUTHIST,506)
            WRITE(OUTHIST,597)
            597     FORMAT ('*****  REGOLITH ****')
            WRITE(OUTHIST,528)(I,I=IUM,MX), &
            (I,I=IABORTMIN,IABORTMAX),(I,I=1,ILU)
            528     FORMAT(' I= ',10I8,//)
            DO  J= JABORTMIN, JABORTMAX
                WRITE(OUTHIST,529)J,(REGOLITH(I,J),I=IUM,MX), &
                (REGOLITH(I,J),I=IABORTMIN,IABORTMAX),&
                (REGOLITH(I,J),I=1,ILU)
                529         FORMAT(I4,6(G12.4),/,3(G12.4))
            ENDDO
            WRITE(OUTHIST,506)
        ENDIF
        IF (MAXGRADIENT > 10000.0) THEN
            WRITE(OUTHIST,510)
            510     FORMAT('***ABORTED DUE TO HIGH GRADIENTS***')
            IF (FLUVIAL_AND_SLOPE_MODELING) CALL CHANNEL_PROPERTIES()
            CALL WRITE_FIRST_DATA_MATRIX()
            CALL WRITE_ALLUVIAL_LOCATIONS(OUTALLUV)
            CALL WRITE_BEDROCK_LOCATIONS(OUTROCK)
            CALL WRITE_LAKE_INFO(OUTLAKE)
            CALL WRITE_EROSION_DEPTH_INDEX(OUTEROSION_DEPTH_INDEX)
            CALL WRITE_DEFORMATION(OUTDEFORM)
            CALL WRITE_SEDIMENT_BASE()
            CALL WRITE_REGOLITH_THICKNESS()
            CALL WRITE_ROCK_RESISTANCE()
            CALL PRINT_SIMULATION_INFORMATION()
            CALL WRITE_ACCELERATED_EROSION_STATE(OUTHIGH)
            CALL WRITE_SHADED_RELIEF_IMAGE()
            CALL WRITE_IMAGE()
            IF (MODEL_GROUNDWATER) THEN
                CALL WRITE_GROUNDWATER_ELEVATION(IOTEMP1)
                CALL WRITE_GROUNDWATER_FLOW(IOTEMP1)
            ENDIF
            IF (DO_MORPHOMETRY) THEN
                CALL PRINT_MORPHOMETRY()
                WRITE_OUT_CHANNELS=.TRUE.
                CALL CALCULATE_TOPO_DIVERGENCE()
                CALL SUMMARIZE_CHANNELS()
                CALL WRITE_SECOND_DATA_MATRIX()
            ENDIF
            WRITE(OUTHIST,433)
            433     FORMAT(/,' ELEVATION')
            CALL SUMMARIZE_MATRIX_DATA(ELEVATION,THEAVG,THEMAX,THEMIN)
            WRITE(OUTHIST,434)
            434     FORMAT(/,' ERODECHANNEL')
            CALL SUMMARIZE_MATRIX_DATA(ERODE_CHANNEL,THEAVG,THEMAX,THEMIN)
            WRITE(OUTHIST,435)
            435     FORMAT(/,' SEDIMENT DIVERGENCE')
            CALL SUMMARIZE_MATRIX_DATA(CFW,THEAVG,THEMAX,THEMIN)
            WRITE(OUTHIST,436)
            436     FORMAT(/,' SEDIMENT YIELD')
            CALL SUMMARIZE_MATRIX_DATA(SEDIMENT_YIELD,THEAVG,THEMAX,THEMIN)
            WRITE(OUTHIST,437)
            437     FORMAT(/,' ERODE_SLOPE')
            CALL SUMMARIZE_MATRIX_DATA(ERODE_SLOPE,THEAVG,THEMAX,THEMIN)
            WRITE(OUTHIST,438)
            438     FORMAT(/,' GRADIENT')
            CALL SUMMARIZE_MATRIX_DATA(D8_GRADIENT,THEAVG,THEMAX,THEMIN)
            WRITE(OUTHIST,439)
            439     FORMAT(' SEDBASE')
            CALL SUMMARIZE_MATRIX_DATA(SEDIMENT_BASE,THEAVG,THEMAX,THEMIN)
            WRITE(OUTHIST,741)
            741     FORMAT(/,' BEDROCK')
            CALL SUMMARIZE_LOGICAL_MATRIX(IS_ROCK_SURFACE)
            WRITE(OUTHIST,742)
            742     FORMAT(/,' REGOLITH')
            CALL SUMMARIZE_MATRIX_DATA(REGOLITH,THEAVG,THEMAX,THEMIN)
            WRITE(OUTHIST,743)
            743     FORMAT(/,' RESISTANCE')
            CALL SUMMARIZE_MATRIX_DATA(RELATIVE_RESISTANCE,THEAVG,THEMAX,THEMIN)
            WRITE(OUTHIST,441)
            441     FORMAT(/,' IS_SEDIMENT_COVERED')
            CALL SUMMARIZE_LOGICAL_MATRIX(IS_SEDIMENT_COVERED)
            445   FORMAT (' TOTAL DELTA VOLUME CHANGE=',G12.5)
        ENDIF
        STOP
    END
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE FINALIZE_FLUVIAL_SLOPE_EROSION()
        USE ERODE_GLOBALS
        IMPLICIT NONE
        INTEGER I,J,IPY,IPX
        REAL (8) :: EEEMIN
        LOGICAL WRITEDETAIL
        !     ****************************************************************
        !      Close down shop
        !   CALLS: GRAD_DISCH_WRITE
        !     **********************************************************************
        WRITEDETAIL=.FALSE.
        WRITE(OUTHIST,511) ITERATION
        511 FORMAT(' FINAL ITERATION=',I6)
        IF (FLUVIAL_AND_SLOPE_MODELING.AND.WRITEDETAIL) THEN
            OPEN(77,FILE='STATE.DAT')
            IF (IS_X_PERIODIC) THEN
                IPX=1
            ELSE
                IPX=0
            ENDIF
            IF (IS_Y_PERIODIC) THEN
                IPY=1
            ELSE
                IPY=0
            ENDIF
            WRITE(77,998) MX,MY,IPX,IPY
            WRITE(77,996) CELL_SIZE
            DO I=1,MX
                DO J=1,MY
                    WRITE(77,997) FLOW_DIRECTION(I,J),ELEVATION(I,J),DRAINAGE_AREA(I,J), &
                    DISCHARGE(I,J),D8_GRADIENT(I,J),BASIN_NUMBER(I,J)
                ENDDO
            ENDDO
            CLOSE(77)
            998   FORMAT(I6,' ',I6,' ',I6,' ',I6)
            996   FORMAT(G13.6)
            997   FORMAT(I6,4(' ',G13.6),' ',I7)
            CALL GRAD_DISCH_WRITE()
        ENDIF
        CLOSE(OUTCHAN)
        CLOSE(OUTRECORD)
        CLOSE(OUTSUMMARY)
        CLOSE(OUTSOURCE)
        CLOSE(OUTREPORT)
        CLOSE(IOTEMP1)
        CLOSE(IOTEMP2)
        CLOSE(OUTCRATER)
        CLOSE(OUTHIST)
        999   EEEMIN=1.0E+36
        DO  J=1,MY
            DO  I=1,MX
                IF (ELEVATION(I,J) < EEEMIN) EEEMIN=ELEVATION(I,J)
            ENDDO
        ENDDO
        OPEN(OUTCONT,FILE=trim(INPUT_DIRECTORY_PATH) // 'CONTINUE.DAT',FORM='UNFORMATTED')
        WRITE(OUTCONT) ITERATION,PRESENT_TIME,EREFERENCE  &
        , ELEVATION,RELATIVE_RESISTANCE,SEDIMENT_BASE,INITIAL_ELEVATION,PREVIOUS_ELEVATION,REGOLITH,DEFORMATION &
        , SEDIMENT_YIELD  &
        ,CFNE,CFNW,CFW,CFN,IS_SEDIMENT_COVERED,ACCELERATED_EROSION,DO_ACCELERATED_EROSION,IS_ROCK_SURFACE,SEDIMENT_DEPOSITED &
        ,EROSION_DEPTH_INDEX,FLOW_DIRECTION,IDO,IFILE1,IFILE2 &
        ,ELAPSEDTIME,DEPOSITWORK,ERODEWORK,CRATERWORK,LAVAWORK,SLOPEWORK &
        ,CUMULATIVE_ELEVATION_CHANGE,CUMULATIVE_CRATERING_CHANGE,CUMULATIVE_EOLIAN_CHANGE,CUMULATIVE_LAVA_CHANGE
        CLOSE(OUTCONT)
        RETURN
    END
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

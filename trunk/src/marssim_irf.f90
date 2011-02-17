    !       ####################################################################################################
    !       this is the main source file for the marssim landform evolution model
    !       copyright (c) 2009 alan d. howard
    !       developer can be contacted by ah6p`virginia.edu and department of environmental sciences, p.o. box 400123,
    !                  university of virginia, charlottesville, va 22904-4123
    !       this program is free software; you can redistribute it and/or modify it under the terms of the gnu general public license
    !         as published by the free software foundation; either version 2 of the
    !         license, or (at your option) any later version.
    !       this program is distributed in the hope that it will be useful, but without any warranty;
    !          without even the implied warranty of merchantability or fitness for a particular purpose. see the gnu
    !          general public license for more details.
    !        you should have received a copy of the gnu general public license along with this program; if not, write to
    !          the free software foundation, inc., 51 franklin street, fifth floor, boston, ma 02110-1301 usa.
    !          a web link:   http://www.gnu.org/licenses/gpl-2.0.txt
    !       ####################################################################################################


    !        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    SUBROUTINE READ_INPUT_PARAMETERS()
		USE ERODE_GLOBALS
        USE CRATER_GLOBALS
        USE EOLIAN_GLOBALS
        USE LAVA_GLOBALS
        USE SEDDEBUG_GLOBALS
        USE SEDROUTE_GLOBALS
        USE AREA_GLOBALS
        USE LAKE_GLOBALS
        USE ACCRETION_GLOBALS
        IMPLICIT NONE
        INTEGER :: I,J,IXDEBUG,IWIDTH,JWIDTH,RANDDIRUSE,STICKYUSE
        INTEGER :: USEWET,STATISTDO
        INTEGER :: ISEDIMENT,ISEDROUTE,ISEDDIFFUSE,ISLOPEUSE,IREADALLUV
        INTEGER :: IVARRATEUSE,RANDTHRESHUSE,HIGHRATEUSE,WEATHER2USE
        INTEGER :: VARYIELDUSE,INEW_SIMULATION,DEFORMUSE
        INTEGER :: IRESINPUT,IREFLECT,IW,IE,IHORIZONTAL_LOWER_BOUNDARY,INON_ERODING_LOWER_BOUNDARY
        INTEGER :: BEDROCKHIGH,USEYPERIODIC,CRUSTUSE,USEXPERIODIC
        INTEGER :: READREGOLITH,IOTEMP3,BEDEXPLICIT,ISEDDEBUG
        INTEGER :: SMOOTHSEDUSE,NEWDIRECTIONUSE,DETACHUSE,IWATERLOWER
        INTEGER :: MX1,MY1,ISOFTCRATER,IQCONSTANT,IOCEANVAR
        INTEGER :: IDIVAVG,QQUSE,SEEPUSE,DOINRIVER
        INTEGER :: IDOCRATER,IDOLAVA,IDOEOLIAN,IDOERODE
        INTEGER :: ITAUAREA,IEVENTS,ISEEPAGEWEATHER,RADIATION_USE
        INTEGER :: IDOOCEAN,ROERINGUSE
        INTEGER :: IUSETOTALEXPOSE,IUSENORMAL,IDOLAKES
        INTEGER :: IEXPFLOW,IDOACCRETION,EXPOSECREEPUSE,IERROR
        INTEGER :: USEFLOWBOUND, INVERSEUSE, TOPEXPOSEUSE, EXPOSESMOOTH
        CHARACTER (80) :: TEXTLINE
        REAL (8) :: HIGH_CONVERGENCE_RUNOFF,LOW_CONVERGENCE_RUNOFF,EEEMIN
        REAL (8) :: READ_ERODIBILITY
        REAL (8) :: DARCIES,VISCOSITY,METRIC_PERMEABILITY,YEARLY_RECHARGE
        REAL (8) :: THEMIN,THEAVG,THEMAX
        REAL :: THEMINIMUM,DEFAULT_TIME_INCREMENT
        INTEGER (4) :: IFOLD
        INTEGER :: ISHOWCALC
        INTEGER :: IFLUXDEPEND,IMODEL_PELAGIC_DEPOSITION
        INTEGER :: ICHOSE
        REAL (8) :: TEMP1,TEMP2,TEMP3
        INTEGER, DIMENSION(:), ALLOCATABLE :: RANDSEED
        INTEGER :: NRANDOM
        EXTERNAL READ_ERODIBILITY
        !        *******************************************************************
        !   MODIFIES:  (MOST SIMULATION GLOBAL PARAMETERS AND MATRICES)
        !   CALLS:  READ_REGOLITH_THICKNESS, SETUP_DISTANCE_WEIGHTING, FIND_MODIFICATION_RANGE
        !           SUMMARIZE_LOGICAL_MATRIX, SUMMARIZE_MATRIX_DATA
        !
        !       **********************************************************************

        !        *******************************************************************
        !        *******************************************************************
        !          r e a d  i n  t h e  s i m u l a t i o n  p a r a m e t e r s
        !        *******************************************************************
        !        *******************************************************************
        !
        !        *******************************************************************
        !        ********************** boundary conditions ************************
        !        *******************************************************************

        IOTEMP3=43
        DEFAULT_TIME_INCREMENT=0.001
        READ(PARAMS,22) TEXTLINE
        22    FORMAT(A)
        WRITE(*,22) TEXTLINE
        READ(PARAMS,*) ISEED
        READ(PARAMS,*)INEW_SIMULATION,IDOERODE,IDOCRATER,IDOLAVA,IDOEOLIAN,    &
        IDOOCEAN, IDOACCRETION , IDOLAKES
        IF (INEW_SIMULATION > 0) THEN
            NEW_SIMULATION=.TRUE.
        ELSE
            NEW_SIMULATION=.FALSE.
        ENDIF
        IF (IDOCRATER > 0) THEN
            MODEL_IMPACT_CRATERING=.TRUE.
        ELSE
            MODEL_IMPACT_CRATERING=.FALSE.
        ENDIF
        IF (IDOLAVA > 0) THEN
            MODEL_LAVA_FLOWS=.TRUE.
        ELSE
            MODEL_LAVA_FLOWS=.FALSE.
        ENDIF
        IF (IDOEOLIAN > 0) THEN
            MODEL_EOLIAN_CHANGES=.TRUE.
        ELSE
            MODEL_EOLIAN_CHANGES=.FALSE.
        ENDIF
        IF (IDOERODE > 0) THEN
            FLUVIAL_AND_SLOPE_MODELING=.TRUE.
        ELSE
            FLUVIAL_AND_SLOPE_MODELING=.FALSE.
        ENDIF
        IF (IDOOCEAN > 0) THEN
            MODEL_OCEAN_LEVEL=.TRUE.
        ELSE
            MODEL_OCEAN_LEVEL=.FALSE.
        ENDIF
        IF (IDOACCRETION > 0) THEN
            MODEL_ACCRETION_AND_ABLATION=.TRUE.
        ELSE
            MODEL_ACCRETION_AND_ABLATION=.FALSE.
        ENDIF
        IF (IDOLAKES > 0) THEN
            MODEL_LAKE_EVAPORATION=.TRUE.
        ELSE
            MODEL_LAKE_EVAPORATION=.FALSE.
        ENDIF
        READ(PARAMS,*)MX,MY,MZ,INPUT_CELL_SIZE,VERTICAL_SCALING_FACTOR,CONVERT_TO_METERS
        IMMX=MX
        JMMX=MY
        LMMX=MX*MY
        RMMX=2*MAX(MX,MY)
        ALLOCATE(ELEVATION(MX,MY),INITIAL_ELEVATION(MX,MY),STAT=IERROR)
        ALLOCATE(CUMULATIVE_EOLIAN_CHANGE(MX,MY),STAT=IERROR)
        ALLOCATE(CUMULATIVE_ELEVATION_CHANGE(MX,MY),CUMULATIVE_LAVA_CHANGE(MX,MY),STAT=IERROR)
        ALLOCATE(CUMULATIVE_CRATERING_CHANGE(MX,MY),STAT=IERROR)
        IF (FLUVIAL_AND_SLOPE_MODELING.OR.MODEL_EOLIAN_CHANGES.OR.MODEL_IMPACT_CRATERING.OR.MODEL_ACCRETION_AND_ABLATION) THEN
            ALLOCATE(ERODE_SLOPE(MX,MY),STAT=IERROR)
        ENDIF
        IF (FLUVIAL_AND_SLOPE_MODELING.OR.MODEL_EOLIAN_CHANGES.OR.MODEL_ACCRETION_AND_ABLATION.OR.MODEL_IMPACT_CRATERING) THEN
            ALLOCATE(D8_GRADIENT(MX,MY),CFNW(MX,MY),STAT=IERROR)
            ALLOCATE(FLOW_DIRECTION(MX,MY),IDO(MX,MY),STAT=IERROR)
        ENDIF
        IF (FLUVIAL_AND_SLOPE_MODELING.OR.MODEL_ACCRETION_AND_ABLATION) THEN
            ALLOCATE(CFW(MX,MY),CFN(MX,MY),STAT=IERROR)
            ALLOCATE(CFNE(MX,MY),STAT=IERROR)
        ENDIF
        IF (FLUVIAL_AND_SLOPE_MODELING) THEN
            ALLOCATE(IS_SEDIMENT_COVERED(MX,MY),STAT=IERROR)
            ALLOCATE(REGOLITH(MX,MY),STAT=IERROR)
            ALLOCATE(IS_ROCK_SURFACE(MX,MY),STAT=IERROR)
            ALLOCATE(SEDIMENT_BASE(MX,MY),STAT=IERROR)
        ENDIF
        IF (FLUVIAL_AND_SLOPE_MODELING) THEN
            ALLOCATE(ERODE_CHANNEL(MX,MY),STAT=IERROR)
            ALLOCATE(RELATIVE_RESISTANCE(MX,MY),STAT=IERROR)
            ALLOCATE(CHANNEL_WIDTH(MX,MY),STAT=IERROR)
            ALLOCATE(SEDIMENT_YIELD(MX,MY),SEDIMENT_FLUX(MX,MY),STAT=IERROR)
            ALLOCATE(EQUILIBRIUM_GRADIENT(MX,MY),DRAINAGE_AREA(MX,MY),STAT=IERROR)
            ALLOCATE(MAXIMUM_SEDIMENT_YIELD(MX,MY),STAT=IERROR)
            ALLOCATE(PREVIOUS_ELEVATION(MX,MY),DISCHARGE(MX,MY),MAXIMUM_DISCHARGE(MX,MY),STAT=IERROR)
            ALLOCATE(DIVERGENCE(MX,MY),STAT=IERROR)
            ALLOCATE(ERODE_REGOLITH_CHANNEL(MX,MY),STAT=IERROR)
            ALLOCATE(DEFORMATION(MX,MY),I_OUTFLOW(LMMX),STAT=IERROR)
            ALLOCATE(J_OUTFLOW(LMMX),DOWNSTREAM_BASIN(LMMX),STAT=IERROR)
            ALLOCATE(BASIN_NUMBER(MX,MY),STAT=IERROR)
            ALLOCATE(EROSION_DEPTH_INDEX(MX,MY),ENCLOSED(LMMX),OVERFLOWS(LMMX),STAT=IERROR)
            ALLOCATE(SUBMERGED(MX,MY),STAT=IERROR)
            ALLOCATE(ACCELERATED_EROSION(MX,MY),DO_ACCELERATED_EROSION(MX,MY),STAT=IERROR)
            ALLOCATE(SEDIMENT_DEPOSITED(MX,MY),IS_INFLUENT_RIVER_LOCATION(MX,MY),STAT=IERROR)
            ALLOCATE(ERODING_LOWER_BOUNDARY(MX),STAT=IERROR)
            ALLOCATE(STARTING_ELEVATION(RMMX),ALLUVIAL_GRADIENT(RMMX),ACTUAL_GRADIENT(RMMX),STAT=IERROR)
            ALLOCATE(STEP_DISTANCE(RMMX),NEW_ELEVATION(RMMX),PROVISIONAL_ELEVATION(RMMX),STAT=IERROR)
            ALLOCATE(IS_BEDROCK_CHANNEL(RMMX),WATER_LEVEL(RMMX),STAT=IERROR)
            ALLOCATE(TOBEADDED(LMMX),LOCAL_BASIN_DISCHARGE(LMMX),QTOADD(LMMX),STAT=IERROR)
            ALLOCATE(NEEDTODO(LMMX),OUTER(LMMX),STAT=IERROR)
            ALLOCATE(BASIN_DRAINAGE_AREA(LMMX),LAKE_OUTLET_ELEVATION(LMMX),STAT=IERROR)
            ALLOCATE(LAKE_SURFACE_ELEVATION(LMMX),SORTING_VECTOR(LMMX),STAT=IERROR)
            ALLOCATE(NOMINAL_ERODE_SLOPE(MX,MY),STAT=IERROR)
            ALLOCATE(DONE_ONCE(MX,MY),STAT=IERROR)
            ALLOCATE(I_LOCATION(RMMX),J_LOCATION(RMMX),STAT=IERROR)
            ALLOCATE(IS_EXIT_BASIN(LMMX),STAT=IERROR)
            ALLOCATE(ROUTED_DISCHARGE(MX,MY),STAT=IERROR)
            ALLOCATE(LAKE_AREA(LMMX),STAT=IERROR)
            IF(MODEL_LAKE_EVAPORATION) THEN
                ALLOCATE(ROUTED_DISCHARGE(MX,MY),STAT=IERROR)
                ALLOCATE(LAKE_VOLUME(LMMX),STAT=IERROR)
                ALLOCATE(LOWEST_LAKE_ELEVATION(LMMX),NEW_BASIN_OUTFLUX(LMMX),STAT=IERROR)
                ALLOCATE(BASIN_OUTFLUX(LMMX),OLDOVERFLOWS(LMMX),STAT=IERROR)
                ALLOCATE(LAKEMIN(LMMX),BASINMIN(LMMX),STAT=IERROR)
                ALLOCATE(NEXTCYCLE(LMMX),NEWGEOMETRY(LMMX),STAT=IERROR)
                ALLOCATE(ISBORDER(LMMX),STAT=IERROR)
                ALLOCATE(BASIN_INFLUX(LMMX),STAT=IERROR)
            ENDIF
        ENDIF
        IF (MODEL_LAVA_FLOWS) THEN
            ALLOCATE(IS_LAVA_COVERED(MX,MY),ACTIVE_LAVA_FLOW(MX,MY),STAT=IERROR)
            ALLOCATE(LAVA_SOURCE_DIRECTION(MX,MY),ERUPTION_AGE(MX,MY),STAT=IERROR)
            ALLOCATE(ILOC(RMMX),JLOC(RMMX),LAVA_ELAPSED_TIME(RMMX),STAT=IERROR)
            ALLOCATE(LAVA_THICKNESS(RMMX),STAT=IERROR)
            ALLOCATE(LAVA_FLOW_PROBABILITY(9,RMMX),STAT=IERROR)
        ENDIF
        VERTICAL_SCALING=VERTICAL_SCALING_FACTOR*CONVERT_TO_METERS
        CELL_SIZE=INPUT_CELL_SIZE*CONVERT_TO_METERS
        READ(PARAMS,*)IVARRATEUSE, BOUNDARY_LOWERING_RATE
        IF (IVARRATEUSE > 0) THEN
            VARIABLE_EROSION_RATE = .TRUE.
        ELSE
            VARIABLE_EROSION_RATE = .FALSE.
        ENDIF
        READ(PARAMS,*)CHANNEL_TIMESTEP_SCALING,DEFAULT_CHANNEL_TIMESTEP, MAXIMUM_TIME_INCREMENT,MINUMUM_TIME_INCREMENT
        READ(PARAMS,*) SEDIMENT_TIMESTEP_SCALING, MASS_WASTING_TIMESTEP_SCALING
        READ(PARAMS,*) ICENT,JCENT,IWIDTH,JWIDTH
        IWINLOW=ICENT-IWIDTH
        IWINHIGH=ICENT+IWIDTH
        JWINLOW=JCENT-JWIDTH
        JWINHIGH=JCENT+JWIDTH
        IF (IWINLOW < 1) IWINLOW=1
        IF (IWINHIGH > MX) IWINHIGH=MX
        IF (JWINLOW < 1) JWINLOW=1
        IF (JWINHIGH > MY) JWINHIGH=MY
        READ(PARAMS,*) IXDEBUG, IEVENTS
        IF (IXDEBUG > 0) THEN
            DO_DEBUGGING=.TRUE.
        ELSE
            DO_DEBUGGING=.FALSE.
        ENDIF
        DO_EVENTS=.FALSE.
        IF (IEVENTS > 0) THEN
            EVENT_DONE=.TRUE.
            OPEN(72,FILE=trim(INPUT_DIRECTORY_PATH) // 'EVENTS.PRM')
            READ(72,*) NUMBER_OF_EVENTS,EVENT_TYPE
            IF (NUMBER_OF_EVENTS > 0) THEN
                DO I=1,NUMBER_OF_EVENTS
                    IF (EVENT_TYPE == 0) THEN
                        READ(72,*) EVENT_TIMES(I)
                    ELSE
                        READ(72,*) EVENT_ITERATIONS(I)
                    ENDIF
                    EVENT_DONE(I)=.FALSE.
                ENDDO
                EVENT_INDEX=1
                DO_EVENTS=.TRUE.
            ENDIF
        ENDIF
        READ(PARAMS,*) STATISTDO
        IF (STATISTDO > 0) THEN
            DO_MORPHOMETRY=.TRUE.
        ELSE
            DO_MORPHOMETRY=.FALSE.
        ENDIF
        !        ********************************************************************
        !        inelev.dat is initial elevations
        !        ********************************************************************
        OPEN(INDATA,FILE=trim(INPUT_DIRECTORY_PATH)//'INELEV.DAT')
        READ(PARAMS,*) IHORIZONTAL_LOWER_BOUNDARY,INON_ERODING_LOWER_BOUNDARY,USEYPERIODIC, &
        USEXPERIODIC,USEFLOWBOUND
        IF (INON_ERODING_LOWER_BOUNDARY > 0) THEN
            NON_ERODING_LOWER_BOUNDARY=.TRUE.
        ELSE
            NON_ERODING_LOWER_BOUNDARY=.FALSE.
        ENDIF
        IF (IHORIZONTAL_LOWER_BOUNDARY > 0) THEN
            HORIZONTAL_LOWER_BOUNDARY=.TRUE.
        ELSE
            HORIZONTAL_LOWER_BOUNDARY=.FALSE.
        ENDIF
        IF (USEYPERIODIC > 0) THEN
            IS_Y_PERIODIC=.TRUE.
            MYY=MY
        ELSE
            IS_Y_PERIODIC=.FALSE.
            MYY=MY-1
        ENDIF
        IF (USEXPERIODIC > 0) THEN
            IS_X_PERIODIC=.TRUE.
        ELSE
            IS_X_PERIODIC=.FALSE.
        ENDIF
        IF (USEFLOWBOUND>0) THEN
            DO_FLOW_BOUNDARIES=.TRUE.
            IS_X_PERIODIC=.FALSE.
            IS_Y_PERIODIC=.FALSE.
            HORIZONTAL_LOWER_BOUNDARY=.FALSE.
            NON_ERODING_LOWER_BOUNDARY=.FALSE.
        ELSE
            DO_FLOW_BOUNDARIES=.FALSE.
        ENDIF
        !        *******************************************************************
        !        ********************** output and recalculation timing *************
        !        *******************************************************************
        READ(PARAMS,22) TEXTLINE
        WRITE(*,22) TEXTLINE
        READ(PARAMS,*)STARTING_ITERATION, PRESENT_TIME
        TOTAL_ITERATIONS = STARTING_ITERATION
        READ(PARAMS,*)MAXIMUM_ITERATION,MAXIMUM_SIMULATION_TIME
        READ(PARAMS,*)ELEVATION_PRINT_INTERVAL
        READ(PARAMS,*)OUTPUT_PRINT_INTERVAL
        READ(PARAMS,*)RECALCULATE_GRADIENT_INTERVAL
        READ(PARAMS,*)WRITE_CHANGE_INTERVAL
        READ(PARAMS,*)KWRITE
        IF (KWRITE == 1)CRITICAL_SOURCE_DIVERGENCE=0.0
        IF (KWRITE == 2)CRITICAL_SOURCE_DIVERGENCE=0.1
        IF (KWRITE == 3)CRITICAL_SOURCE_DIVERGENCE=0.2
        IF (KWRITE == 4)CRITICAL_SOURCE_DIVERGENCE=0.4
        IF (KWRITE == 5)CRITICAL_SOURCE_DIVERGENCE=0.8
        IF (KWRITE == 6)CRITICAL_SOURCE_DIVERGENCE=1.6
        READ(PARAMS,550)(WRITETYPE(I),I=1,13)
        550 FORMAT(13I1)
        READ(PARAMS,*)ONEONLY
        IF (ONEONLY > 0) THEN
            WRITE_ABSOLUTE_ELEVATION = .TRUE.
        ELSE
            WRITE_ABSOLUTE_ELEVATION = .FALSE.
        ENDIF
        READ(PARAMS,*) IMAGE_OUTPUT_INTERVAL,DIVERGEINTERVAL
        !        *******************************************************************
        !        ********************* flow routing and drainage area ***************
        !        *******************************************************************
        READ(PARAMS,22) TEXTLINE
        WRITE(*,22) TEXTLINE
        READ(PARAMS,*) USEWET, IQCONSTANT
        IF (IQCONSTANT > 0) THEN
            RESCALE_DISCHARGES=.TRUE.
        ELSE
            RESCALE_DISCHARGES=.FALSE.
        ENDIF
        IF (USEWET > 0) THEN
            COMPLETE_RUNOFF=.TRUE.
        ELSE
            COMPLETE_RUNOFF=.FALSE.
        ENDIF
        READ(PARAMS,*) RANDDIRUSE
        IF (RANDDIRUSE > 0) THEN
            USE_RANDOM_FLOWROUTING=.TRUE.
        ELSE
            USE_RANDOM_FLOWROUTING=.FALSE.
        ENDIF
        READ(PARAMS,*) DISCHARGE_CONSTANT,DISCHARGE_EXPONENT
        DISCHARGE_SCALE_FACTOR=DISCHARGE_CONSTANT**(1.0/DISCHARGE_EXPONENT)
        TEMPQCONSTANT=DISCHARGE_CONSTANT
        LASTQCONSTANT=DISCHARGE_CONSTANT
        READ(PARAMS,*) AREAFACTOR, RAINDEPTH,RAINSTD
        READ(PARAMS,*) VARYIELDUSE,MEAN_CONVERGENCE_RUNOFF,HIGH_CONVERGENCE_RUNOFF,LOW_CONVERGENCE_RUNOFF
        READ(PARAMS,*) NCALCEVAP,NCHANGEEVAP,EVAPORATION_MEAN,EVAPORATION_STANDARD_DEVIATION
        READ(PARAMS,*) IMODEL_PELAGIC_DEPOSITION,WASHLOAD_FRACTION,PELAGICCREEP
        IF (IMODEL_PELAGIC_DEPOSITION > 0) THEN
            MODEL_PELAGIC_DEPOSITION=.TRUE.
        ELSE
            MODEL_PELAGIC_DEPOSITION=.FALSE.
        ENDIF
        IF(FLUVIAL_AND_SLOPE_MODELING.AND.MODEL_PELAGIC_DEPOSITION) THEN
            ALLOCATE(PELAGICCHANGE(MX,MY),PELAGICEDIFF(MX,MY),STAT=IERROR)
            ALLOCATE(PELAGIC_SEDIMENT_VOLUME(LMMX),STAT=IERROR)
        ENDIF
        !       *****************************************************************
        !        Read the location of an incoming river, its discharge and
        !           sediment load
        !       *****************************************************************
        READ(PARAMS,*) DOINRIVER
        IF (FLUVIAL_AND_SLOPE_MODELING) THEN
            IS_INFLUENT_RIVER_LOCATION=.FALSE.
        ENDIF
        IF (DOINRIVER > 0) THEN
            HAVE_INFLUENT_RIVERS=.TRUE.
            OPEN(76,FILE=trim(INPUT_DIRECTORY_PATH) // 'INRIVER.PRM')
            READ(76,*) NUMBER_OF_INFLUENT_RIVERS
            DO I=1,NUMBER_OF_INFLUENT_RIVERS
                READ(76,*) I_INFLUENT_RIVER(I),J_INFLUENT_RIVER(I),INFLUENT_RIVER_DISCHARGE(I),INFLUENT_RIVER_SEDIMENT_LOAD(I)
                IS_INFLUENT_RIVER_LOCATION(I_INFLUENT_RIVER(I),J_INFLUENT_RIVER(I))=.TRUE.
            ENDDO
            CLOSE(76)
        ELSE
            HAVE_INFLUENT_RIVERS=.FALSE.
        ENDIF

        !       *****************************************************************
        !        Calculate parameters for calculating specific runoff as a function
        !           of local divergence
        !       *****************************************************************
        IF (VARYIELDUSE > 0) THEN
            DIVERGENCE_DEPENDENT_RUNOFF=.TRUE.
            RUNOFF_CONSTANT_1=(1.0-LOW_CONVERGENCE_RUNOFF)/2.0
            RUNOFF_CONSTANT_2=(1.0+LOW_CONVERGENCE_RUNOFF)/2.0
            RUNOFF_CONSTANT_3=(1.0-RUNOFF_CONSTANT_2)/((HIGH_CONVERGENCE_RUNOFF-MEAN_CONVERGENCE_RUNOFF)* &
            DSQRT(RUNOFF_CONSTANT_1**2+(1.0-RUNOFF_CONSTANT_2)**2))
        ELSE
            DIVERGENCE_DEPENDENT_RUNOFF=.FALSE.
        ENDIF
        !        *******************************************************************
        !        ********************** bedrock channel parameters **************
        !        *******************************************************************
        READ(PARAMS,22) TEXTLINE
        WRITE(*,22) TEXTLINE
        READ(PARAMS,*)DETACHUSE, BEDEXPLICIT
        IF (BEDEXPLICIT > 0) THEN
            EXPLICIT_CHANNEL_BED_STATE=.TRUE.
        ELSE
            EXPLICIT_CHANNEL_BED_STATE=.FALSE.
        ENDIF
        IF (DETACHUSE > 0) THEN
            DO_FLUVIAL_DETACHMENT=.TRUE.
        ELSE
            DO_FLUVIAL_DETACHMENT=.FALSE.
        ENDIF
        READ(PARAMS,*) IFLUXDEPEND, ROCK_TENSILE_STRENGTH
        USE_SKLAR_BED_ABRASION=.FALSE.
        USE_WHIPPLE_BED_ABRASION=.FALSE.
        IF (IFLUXDEPEND == 1) THEN
            USE_SKLAR_BED_ABRASION=.TRUE.
            CHANNEL_TIMESTEP_SCALING=CHANNEL_TIMESTEP_SCALING*ROCK_TENSILE_STRENGTH**2
        ENDIF
        IF (IFLUXDEPEND > 1) THEN
            USE_WHIPPLE_BED_ABRASION=.TRUE.
        ENDIF
        READ(PARAMS,*)BEDROCK_DISCHARGE_EXPONENT
        READ(PARAMS,*)BEDROCK_GRADIENT_EXPONENT
        READ(PARAMS,*)BEDROCK_ERODIBILITY,MANNING,GRAVITY
        READ(PARAMS,*)DETACHMENT_CRITICAL_SHEAR
        PREVIOUS_CRITICAL_SHEAR=DETACHMENT_CRITICAL_SHEAR
        READ(PARAMS,*)REGOLITH_ERODIBILITY_FACTOR,REGOLITH_CRITICAL_SHEAR_FACTOR
        READ(PARAMS,*)CHANNEL_WIDTH_CONSTANT,CHANNEL_WIDTH_EXPONENT
        READ(PARAMS,*)ITAUAREA,VEGETATION_UPLAND_RESISTANCE,VEGETATION_CHANNEL_RESISTANCE, &
        VEGETATION_AREA_MINIMUM,VEGETATION_AREA_MAXIMUM
        IF (ITAUAREA > 0) THEN
            VARIABLE_VEGETATION_RESISTANCE=.TRUE.
            VEGETATION_FACTOR_SLOPE=(DLOG(VEGETATION_UPLAND_RESISTANCE)-DLOG(VEGETATION_CHANNEL_RESISTANCE))/ &
            (DLOG(VEGETATION_AREA_MINIMUM)-DLOG(VEGETATION_AREA_MAXIMUM))
            VEGETATION_FACTOR_INTERCEPT=DLOG(VEGETATION_UPLAND_RESISTANCE)- &
            VEGETATION_FACTOR_SLOPE*DLOG(VEGETATION_AREA_MINIMUM)
        ELSE
            VARIABLE_VEGETATION_RESISTANCE=.FALSE.
            VEGETATION_FACTOR_SLOPE=0.0
            VEGETATION_FACTOR_INTERCEPT=0.0
        ENDIF
        DISCHARGE_COEFFICIENT=WATER_DENSITY*(MANNING/CHANNEL_WIDTH_CONSTANT)**0.6*GRAVITY * &
        (9.8/GRAVITY)**0.3
        PREVIOUS_DISCHARGE_COEFFICIENT=DISCHARGE_COEFFICIENT
        !        *******************************************************************
        !        ********************** alluvial channel parameters *************
        !        *******************************************************************
        READ(PARAMS,22) TEXTLINE
        WRITE(*,22) TEXTLINE
        READ(PARAMS,*) ISEDDEBUG
        READ(PARAMS,*) ISEDIMENT,ISEDROUTE,ISEDDIFFUSE,IREFLECT
        IF (ISEDDEBUG > 0) THEN
            WRITE_SEDIMENT_DIAGNOSTICS=.TRUE.
        ELSE
            WRITE_SEDIMENT_DIAGNOSTICS=.FALSE.
        ENDIF
        IF (ISEDIMENT > 0) THEN
            DO_SEDIMENT_TRANSPORT=.TRUE.
        ELSE
            DO_SEDIMENT_TRANSPORT=.FALSE.
        ENDIF
        IF (ISEDROUTE > 0) THEN
            DO_SEDIMENT_ROUTING=.TRUE.
        ELSE
            DO_SEDIMENT_ROUTING=.FALSE.
        ENDIF
        IF (ISEDDIFFUSE > 0) THEN
            DO_SEDIMENT_DIFFUSION=.TRUE.
        ELSE
            DO_SEDIMENT_DIFFUSION=.FALSE.
        ENDIF
        IF (IREFLECT > 0) THEN
            NO_FLUX_LOWER_BOUNDARY=.TRUE.
        ELSE
            NO_FLUX_LOWER_BOUNDARY=.FALSE.
        ENDIF
        READ(PARAMS,*) SMOOTHSEDUSE,NEWDIRECTIONUSE,ALLUVIUM_SMOOTHING_FACTOR
        IF (SMOOTHSEDUSE > 0) THEN
            DO_ALLUVIAL_SMOOTHING=.TRUE.
        ELSE
            DO_ALLUVIAL_SMOOTHING=.FALSE.
        ENDIF
        IF (NEWDIRECTIONUSE > 0) THEN
            DO_ALLUVIAL_REROUTING=.TRUE.
        ELSE
            DO_ALLUVIAL_REROUTING=.FALSE.
        ENDIF
        READ(PARAMS,*) IREADALLUV
        READ(PARAMS,*) WATERUSE, OCEAN_ELEVATION, DELTA_FORESET_GRADIENT
        IF (WATERUSE > 0) THEN
            USE_AN_OCEAN = .TRUE.
        ELSE
            USE_AN_OCEAN = .FALSE.
        ENDIF
        READ(PARAMS,*) SEDIMENT_1_EXPONENT,SEDIMENT_2_EXPONENT,EFFECTIVE_DISCHARGE_RATIO
        READ(PARAMS,*) SEDIMENT_GRADIENT_EXPONENT,SEDIMENT_TRANSPORT_EXPONENT , TRANSPORTFACTOR,  &
        FLOW_FRACTION
        READ(PARAMS,*) SEDIMENT_POROSITY,SEDIMENT_SPECIFIC_GRAVITY, &
        GRAIN_SIZE, TRANSPORT_CRITICAL_DIM_SHEAR
        READ(PARAMS,*) STICKYUSE,STICKY_ROUTING_CRITICAL_VALUE
        IF (STICKYUSE > 0) THEN
            STICKY_SEDIMENT_ROUTING=.TRUE.
        ELSE
            STICKY_SEDIMENT_ROUTING=.FALSE.
        ENDIF
        READ(PARAMS,*) BEDLOAD_FRACTION
        SEDIMENT_DISCHARGE_FACTOR=((MANNING/CHANNEL_WIDTH_CONSTANT)**0.6)*EFFECTIVE_DISCHARGE_RATIO**SEDIMENT_2_EXPONENT &
        /(GRAIN_SIZE*(SEDIMENT_SPECIFIC_GRAVITY-1.0)*(GRAVITY/9.8)**0.3)
        SEDIMENT_CONSTANT=CHANNEL_WIDTH_CONSTANT  &
        *TRANSPORTFACTOR*FLOW_FRACTION*GRAIN_SIZE**1.5 &
        *SECONDS_PER_YEAR*(GRAVITY*(SEDIMENT_SPECIFIC_GRAVITY-1.0))**0.5 &
        / ((1.0-SEDIMENT_POROSITY)*BEDLOAD_FRACTION)
        SKLAR_MULT=(SKLAR_FACTOR*GRAVITY*(SEDIMENT_SPECIFIC_GRAVITY-1))/ROCK_TENSILE_STRENGTH**2
        SKLAR_MULT=SKLAR_MULT*(1.0-SEDIMENT_POROSITY)*1000.0*SEDIMENT_SPECIFIC_GRAVITY &
        *FLOW_FRACTION

        !        *******************************************************************
        !        ****************  regolith weathering parameters  ******************
        !        *******************************************************************
        READ(PARAMS,22) TEXTLINE
        WRITE(*,22) TEXTLINE
        READ(PARAMS,*)ROCK_WEATHERING_RATE,WEATHER_DECAY_RATE,INITIAL_REGOLITH_THICKNESS
        READ(PARAMS,*)WEATHER2USE,WEATHERING_TERM_2,WEATHERING_DECAY_2
        READ(PARAMS,*)ISEEPAGEWEATHER,SEEPAGE_WEATHERING_SCALING,SEEPAGE_WEATHERING_EXPONENT
        IF (ISEEPAGEWEATHER > 0) THEN
            SEEPAGE_WEATHERING=.TRUE.
        ELSE
            SEEPAGE_WEATHERING=.FALSE.
        ENDIF
        IF (WEATHER2USE > 0) THEN
            TWO_TERM_WEATHERING = .TRUE.
        ELSE
            TWO_TERM_WEATHERING = .FALSE.
        ENDIF
        READ(PARAMS,*)CRITICAL_BEDROCK_GRADIENT,WEATHER_MULT,WEATHER_DIVERGENCE
        READ(PARAMS,*) READREGOLITH
        !        *******************************************************************
        !        *****************  mass wasting parameters  ************************
        !        *******************************************************************
        READ(PARAMS,22) TEXTLINE
        WRITE(*,22) TEXTLINE
        READ(PARAMS,*)ISLOPEUSE
        IF (ISLOPEUSE > 0) THEN
            DO_MODEL_SLOPES=.TRUE.
        ELSE
            DO_MODEL_SLOPES=.FALSE.
        ENDIF
        READ(PARAMS,*)SLOPE_DIFFUSIVITY
        READ(PARAMS,*)ALVCREEPFAC
        READ(PARAMS,*)CRITICAL_GRADIENT_USE, ROERINGUSE
        IF (CRITICAL_GRADIENT_USE > 0) THEN
            USE_CRITICAL_SLOPE_CRADIENT =.TRUE.
        ELSE
            USE_CRITICAL_SLOPE_CRADIENT = .FALSE.
        ENDIF
        READ(PARAMS,*)CRITICAL_SLOPE_GRADIENT
        READ(PARAMS,*)SLOPE_GRADIENT_EXPONENT
        READ(PARAMS,*)MAXIMUM_DIFFUSIVITY_INCREASE
        READ(PARAMS,*)SLOPE_FAILURE_DIFFUSIVITY
        IF (ROERINGUSE > 0) THEN
            SLOPE_GRADIENT_EXPONENT=2.0
            SLOPE_FAILURE_DIFFUSIVITY=1.0
            USE_ROERING_MASS_WASTING=.TRUE.
        ELSE
            USE_ROERING_MASS_WASTING=.FALSE.
        ENDIF
        CRITICAL_GRADIENT_TERM = 1.0/CRITICAL_SLOPE_GRADIENT**SLOPE_GRADIENT_EXPONENT
        !        *******************************************************************
        !        *********************** bistable erosion  **************************
        !        *******************************************************************
        READ(PARAMS,22) TEXTLINE
        WRITE(*,22) TEXTLINE
        READ(PARAMS,*)HIGHRATEUSE,BEDROCKHIGH,LOW_EROSION_THRESHOLD,HIGH_EROSION_THRESHOLD
        READ(PARAMS,*)EROSION_RATE_CHANGE_LAG,BISTABLE_CRITICAL_SHEAR,BISTABLE_RUNOFF_FACTOR,BISTABLE_BEDROCK_ERODIBILITY
        IF (BEDROCKHIGH > 0) THEN
            USE_BISTABLE_BEDROCK=.TRUE.
        ELSE
            USE_BISTABLE_BEDROCK=.FALSE.
        ENDIF
        !        *******************************************************************
        !        ******** random variability of simulation parameters    ***************
        !        *******************************************************************
        READ(PARAMS,22) TEXTLINE
        WRITE(*,22) TEXTLINE
        READ(PARAMS,*) RANDTHRESHUSE, RANDDISCHUSE, &
        CRITICAL_SHEAR_VARIABILITY,DISCHARGE_COEFF_VARIATION,OMEGA_WEIGHT
        USE_RANDOM_DISCHARGE=.FALSE.
        IF (RANDTHRESHUSE > 0) THEN
            RANDOM_CRITICAL_SHEAR = .TRUE.
        ELSE
            RANDOM_CRITICAL_SHEAR = .FALSE.
            IF (RANDDISCHUSE > 0) THEN
                USE_RANDOM_DISCHARGE=.TRUE.
            ELSE
                USE_RANDOM_DISCHARGE=.FALSE.
            ENDIF
        ENDIF
        IF (OMEGA_WEIGHT < 1.0)    THEN
            OCORRECT = DSQRT((EXP(1.0)-1.0)/(DEXP(OMEGA_WEIGHT)-1.0))
        ELSE
            OCORRECT = 1.0
        ENDIF
        !        *******************************************************************
        !        ***************  3-d spatially variable bedrock resistance  ********
        !        *******************************************************************
        READ(PARAMS,22) TEXTLINE
        WRITE(*,22) TEXTLINE
        READ(PARAMS,*)SLOPEVARUSE
        READ(PARAMS,*)IRESINPUT,DIFFUSIVITY_VARIABILITY
        READ(PARAMS,*)CRUSTUSE,SURFACE_LAYER_THICKNESS,SURFACE_LAYER_RESISTANCE
        IF (HIGHRATEUSE > 0) THEN
            IF (CRUSTUSE > 0) THEN
                WRITE(OUTHIST,3567)
                3567          FORMAT(' CANNOT USE BOTH HIGHRATE THRESHOLD AND CRUST')
                STOP
            ENDIF
            BISTABLE_FLUVIAL_EROSION=.TRUE.
        ELSE
            BISTABLE_FLUVIAL_EROSION=.FALSE.
        ENDIF
        IF (CRUSTUSE > 0) THEN
            IF (HIGHRATEUSE > 0) THEN
                WRITE (OUTHIST,3567)
            ENDIF
            RESISTANT_SURFACE_LAYER=.TRUE.
        ELSE
            RESISTANT_SURFACE_LAYER=.FALSE.
        ENDIF
        !        *******************************************************************
        !        Calculate parameters for lognormal distribution for scaling input
        !           resistance values
        !        *******************************************************************
        IF (IRESINPUT > 0) THEN
            SCALE_3D_ROCK_ERODIBILITY=.TRUE.
            SIGMANORM=DLOG(DIFFUSIVITY_VARIABILITY)/1.29
            SIGMASQ=SIGMANORM**2/2.0
            RANDMULT=1.0/DEXP(-SIGMASQ)
        ELSE
            SCALE_3D_ROCK_ERODIBILITY=.FALSE.
        ENDIF
        READ(PARAMS,*)VERTICAL_RESISTANCE_SCALING
        !        *******************************************************************
        !        **************** rock and surface deformation *********************
        !        *******************************************************************
        READ(PARAMS,22) TEXTLINE
        WRITE(*,22) TEXTLINE
        READ(PARAMS,*) DEFORMUSE,DEFORMSCALE
        !        *******************************************************************
        !        *********************** groundwater flow **************************
        !        *******************************************************************
        READ(PARAMS,22) TEXTLINE
        WRITE(*,22) TEXTLINE
        READ (PARAMS,*) SEEPUSE, SEEPAGE_ITERATION_INTERVAL,IWATERLOWER,ISHOWCALC &
        ,IEXPFLOW
        IF (SEEPUSE > 0) THEN
            MODEL_GROUNDWATER=.TRUE.
        ELSE
            MODEL_GROUNDWATER=.FALSE.
        ENDIF
        IF (MODEL_GROUNDWATER) THEN
            ALLOCATE(TRANSMISSIVITY_TERM(MX,MY),WATER_ELEVATION(MX,MY),STAT=IERROR)
            ALLOCATE(FILTERED_GROUNDWATER_FLUX(MX,MY),STAT=IERROR)
            ALLOCATE(GROUNDWATER_FLUX(MX,MY),STAT=IERROR)
            ALLOCATE(FIXED_HEAD(MX,MY),STAT=IERROR)
        ENDIF
        IF (IEXPFLOW > 0) THEN
            EXPONENTIAL_PERMEABILITY_DECAY=.TRUE.
        ELSE
            EXPONENTIAL_PERMEABILITY_DECAY=.FALSE.
        ENDIF
        IF (IWATERLOWER > 0) THEN
            PERMEABILITY_RESCALING=.TRUE.
        ELSE
            PERMEABILITY_RESCALING=.FALSE.
        ENDIF
        IF (ISHOWCALC > 0) THEN
            SHOW_GROUNDWATER_CALCULATIONS=.TRUE.
        ELSE
            SHOW_GROUNDWATER_CALCULATIONS=.FALSE.
        ENDIF
        READ(PARAMS,*) YEARLY_RECHARGE,VISCOSITY,DARCIES,GROUNDWATER_DEPTH_SCALE &
        ,GROUNDWATER_FLOW_FRACTION,INITIAL_GROUNDWATER_DEPTH &
        ,EPOWER
        IF (EXPONENTIAL_PERMEABILITY_DECAY) THEN
            INITIAL_GROUNDWATER_DEPTH=GROUNDWATER_DEPTH_SCALE*INITIAL_GROUNDWATER_DEPTH
            GROUNDWATER_DEPTH_SCALE=-DLOG(0.5D0)/(GROUNDWATER_DEPTH_SCALE**EPOWER)
        ENDIF
        GROUNDWATER_RECHARGE_RATE=YEARLY_RECHARGE/SECONDS_PER_YEAR
        METRIC_PERMEABILITY=DARCIES/1.0E+12
        HYDRAULIC_CONDUCTIVITY=METRIC_PERMEABILITY*WATER_DENSITY*GRAVITY/VISCOSITY
        GROUNDWATER_SCALE_FACTOR=HYDRAULIC_CONDUCTIVITY/(CELL_SIZE*CELL_SIZE*2.0)
        IF (EXPONENTIAL_PERMEABILITY_DECAY) GROUNDWATER_SCALE_FACTOR=GROUNDWATER_SCALE_FACTOR/GROUNDWATER_DEPTH_SCALE
        TEMP1=GROUNDWATER_RECHARGE_RATE/GROUNDWATER_SCALE_FACTOR
        WRITE(*,4215) HYDRAULIC_CONDUCTIVITY,GROUNDWATER_SCALE_FACTOR,TEMP1
        WRITE(OUTHIST,4215) HYDRAULIC_CONDUCTIVITY,GROUNDWATER_SCALE_FACTOR,TEMP1
        4215  FORMAT(' HYDRAULIC CONDUCTIVITY=',G12.5,' GROUNDWATER_SCALE_FACTOR=',G12.5, &
        ' GROUNDWATER_RECHARGE_RATE/GROUNDWATER_SCALE_FACTOR=',G12.5)
        READ(PARAMS,*) MAXIMUM_GROUNDWATER_ITERATIONS,MAXIMUM_GOUNDWATER_ERROR,GROUNDWATER_RELAXATION_CONST
        IF (MODEL_GROUNDWATER.AND.(GROUNDWATER_FLOW_FRACTION > 0.5)) THEN
            SEDIMENT_DISCHARGE_FACTOR=SEDIMENT_DISCHARGE_FACTOR/EFFECTIVE_DISCHARGE_RATIO
            SEDIMENT_CONSTANT=SEDIMENT_CONSTANT/FLOW_FRACTION
            SKLAR_MULT=SKLAR_MULT/FLOW_FRACTION
        ENDIF
        WRITE(OUTHIST,32001) SEDIMENT_DISCHARGE_FACTOR,SEDIMENT_CONSTANT,SKLAR_MULT
        WRITE(*,32001) SEDIMENT_DISCHARGE_FACTOR,SEDIMENT_CONSTANT,SKLAR_MULT
        32001 FORMAT('CALCULATED SEDIMENT_DISCHARGE_FACTOR=',G12.5,' SEDIMENT_CONSTANT=',G12.5 &
        ,' SLARMULT=',G12.5)
        READ (PARAMS,*) QQUSE,IDIVAVG
        IF (QQUSE > 0) THEN
            USE_GROUNDWATER_FLOW=.TRUE.
        ELSE
            USE_GROUNDWATER_FLOW=.FALSE.
        ENDIF
        IF (IDIVAVG > 0) THEN
            SEEPAGE_AVERAGING=.TRUE.
        ELSE
            SEEPAGE_AVERAGING=.FALSE.
        ENDIF
        !       *********************************************************************
        !
        !       Part 2:  Read in parameters for eolian erosion-deposition
        !
        !       *********************************************************************
        READ(PARAMS,22) TEXTLINE
        WRITE(*,22) TEXTLINE
        READ(PARAMS,*) EOLIAN_EVENT_PROBABILITY
        IF (FLUVIAL_AND_SLOPE_MODELING) EOLIAN_EVENT_PROBABILITY=EOLIAN_EVENT_PROBABILITY/DEFAULT_TIME_INCREMENT
        READ(PARAMS,*) EOLIAN_TIME_INCREMENT, IUSETOTALEXPOSE,IUSENORMAL
        IF (IUSETOTALEXPOSE > 0) THEN
            USE_TOTAL_EXPOSURE=.TRUE.
        ELSE
            USE_TOTAL_EXPOSURE=.FALSE.
        ENDIF
        IF (IUSENORMAL > 0) THEN
            DEFAULT_EOLIAN_PROCESS=.TRUE.
        ELSE
            DEFAULT_EOLIAN_PROCESS=.FALSE.
        ENDIF
        READ(PARAMS,*) MINIMUM_EOLIAN_DEPOSIT_RATE,MAXIMUM_EOLIAN_DEPOSIT_RATE
        READ(PARAMS,*) ICHOSE
        EOLIAN_CONSTANT_1=(MAXIMUM_EOLIAN_DEPOSIT_RATE-MINIMUM_EOLIAN_DEPOSIT_RATE)/2.0
        EOLIAN_CONSTANT_2=(MAXIMUM_EOLIAN_DEPOSIT_RATE+MINIMUM_EOLIAN_DEPOSIT_RATE)/2.0
        IF (ICHOSE == 1) THEN
            READ(PARAMS,*) EXPOSURE_10_PERCENT,EXPOSURE_90_PERCENT
            EXPOSURE_50_PERCENT=(EXPOSURE_10_PERCENT+EXPOSURE_90_PERCENT)/2.0
            EOLIAN_CONSTANT_3=(MAXIMUM_EOLIAN_DEPOSIT_RATE-EOLIAN_CONSTANT_2)/((EXPOSURE_90_PERCENT-EXPOSURE_50_PERCENT)* &
            DSQRT(EOLIAN_CONSTANT_1**2+(MAXIMUM_EOLIAN_DEPOSIT_RATE-EOLIAN_CONSTANT_2)**2))
        ENDIF
        IF (ICHOSE == 2) THEN
            READ(PARAMS,*) EXPOSURE_50_PERCENT,EXPOSURE_90_PERCENT
            EOLIAN_CONSTANT_3=(MAXIMUM_EOLIAN_DEPOSIT_RATE-EOLIAN_CONSTANT_2)/((EXPOSURE_90_PERCENT-EXPOSURE_50_PERCENT)*  &
            DSQRT(EOLIAN_CONSTANT_1**2+(MAXIMUM_EOLIAN_DEPOSIT_RATE-EOLIAN_CONSTANT_2)**2))
        ENDIF
        IF (ICHOSE == 3) THEN
            IF (MINIMUM_EOLIAN_DEPOSIT_RATE < 0.0) THEN
                READ(PARAMS,*) ZERO_PERCENT_EXPOSURE,EXPOSURE_90_PERCENT
                TEMP2=EOLIAN_CONSTANT_2/DSQRT(EOLIAN_CONSTANT_1**2-EOLIAN_CONSTANT_2**2)
                TEMP3=(MAXIMUM_EOLIAN_DEPOSIT_RATE-EOLIAN_CONSTANT_2)/DSQRT(EOLIAN_CONSTANT_1**2+ &
                (MAXIMUM_EOLIAN_DEPOSIT_RATE-EOLIAN_CONSTANT_2)**2)
                EXPOSURE_50_PERCENT=(ZERO_PERCENT_EXPOSURE+EXPOSURE_90_PERCENT*TEMP2/TEMP3)/(1.0+TEMP2/TEMP3)
                EOLIAN_CONSTANT_3=TEMP3/(EXPOSURE_90_PERCENT-EXPOSURE_50_PERCENT)
            ELSE
                WRITE(*,350)
                350           FORMAT(' RL MUST BE LESS THAN ZERO')
                STOP
            ENDIF
        ENDIF
        IF (ICHOSE == 4) THEN
            IF (MINIMUM_EOLIAN_DEPOSIT_RATE < 0.0) THEN
                READ(PARAMS,*)RATE0,EXPOSURE_50_PERCENT
                IF (RATE0 >= ((MINIMUM_EOLIAN_DEPOSIT_RATE+MAXIMUM_EOLIAN_DEPOSIT_RATE)/2.0)) THEN
                    WRITE(*,361)
                    361               FORMAT('R0 MUST BE LESS THAT MEAN RATE')
                    STOP
                ENDIF
                IF (EXPOSURE_50_PERCENT <= 0) THEN
                    WRITE(*,362)
                    362               FORMAT(' E50 MUST BE GREATER THAN ZERO')
                    STOP
                ENDIF
                EOLIAN_CONSTANT_3=-(RATE0-EOLIAN_CONSTANT_2)/(EXPOSURE_50_PERCENT*  &
                DSQRT(EOLIAN_CONSTANT_1**2-(RATE0-EOLIAN_CONSTANT_2)**2))
            ELSE
                WRITE(*,350)
                STOP
            ENDIF
        ENDIF
        !       ********************************************************************
        !       Distance weighting factors used for either simulated eolian erosion/deposition or accretion/ablation
        !       ********************************************************************
        READ(PARAMS,*) DISTANCE_DECAY_FACTOR, WEIGHTING_CALCULATION_DISTANCE
        IF (DISTANCE_DECAY_FACTOR /= 0.0) THEN
            MAXIMUM_WEIGHT_DISTANCE=MIN(WEIGHTING_CALCULATION_DISTANCE,INT(-DLOG(0.05D0)/(DISTANCE_DECAY_FACTOR*CELL_SIZE)))
        ELSE
            MAXIMUM_WEIGHT_DISTANCE=WEIGHTING_CALCULATION_DISTANCE
        ENDIF
        WEIGHTING_DECAY_FACTOR=-DISTANCE_DECAY_FACTOR*CELL_SIZE
        WRITE(OUTHIST,3563) MAXIMUM_WEIGHT_DISTANCE,WEIGHTING_DECAY_FACTOR
        3563  FORMAT(' WEIGHTING RANGE=',I6,' WEIGHTING DECAY=',G12.5)
        CALL SETUP_DISTANCE_WEIGHTING()
        !       ********************************************************************
        !
        !       Part 3:  Read in lavaflow parameters
        !
        !       ********************************************************************
        !        *******************************************************************
        !         number_of_lava_sources is the number of lava sources in the simulation domain
        !           i_lava and j_lava hold the i and j coordinates of these sources
        !        *******************************************************************
        READ(PARAMS,22) TEXTLINE
        WRITE(*,22) TEXTLINE
        READ(PARAMS,*) NUMBER_OF_LAVA_SOURCES
        DO I=1,NUMBER_OF_LAVA_SOURCES
            READ(PARAMS,*) I_LAVA(I),J_LAVA(I)
        ENDDO
        !        *******************************************************************
        !         lava_event_probability is the probability of a lava flow event during a single
        !               iteration
        !        *******************************************************************
        READ(PARAMS,*) LAVA_EVENT_PROBABILITY
        IF (FLUVIAL_AND_SLOPE_MODELING) LAVA_EVENT_PROBABILITY=LAVA_EVENT_PROBABILITY/DEFAULT_TIME_INCREMENT
        !        *******************************************************************
        !         minimum_lava_flow_slope is the minimum gradient for lava flow at the flow source
        !        *******************************************************************
        READ(PARAMS,*) MINIMUM_LAVA_FLOW_SLOPE
        !        *******************************************************************
        !         lava_flow_thickness is the assumed thickness of individual lava flow deposits
        !          -- assumed to be spatially uniforn and temporally constant
        !        *******************************************************************
        READ(PARAMS,*) LAVA_FLOW_THICKNESS
        !        *******************************************************************
        !         minimum_flow_thickness is the minimum thickness of lava flow that can flow into
        !           an adjoining cell
        !        *******************************************************************
        READ(PARAMS,*) MINIMUM_FLOW_THICKNESS
        !        *******************************************************************
        !         new_segment_interval is the number if iterations between reassessment of
        !           which matrix points are potential sources for new flow segments
        !        *******************************************************************
        READ(PARAMS,*) NEW_SEGMENT_INTERVAL
        !        *******************************************************************
        !         source_change_interval is the number of iterations between changeover
        !           between different flow sources
        !        *******************************************************************
        READ(PARAMS,*) SOURCE_CHANGE_INTERVAL
        !        *******************************************************************
        !         eruption_stop_probability is the probability, per iteration, that the existing flow
        !           will solidify and stop being active_lava_flow -- presumably due to
        !           interruption of the flow source.  If the flow becomes inlactive
        !           a new flow starts from the source
        !        *******************************************************************
        READ(PARAMS,*) ERUPTION_STOP_PROBABILITY
        !        *******************************************************************
        !         this is the lower limit of probability for a cell to be a source for
        !           a new flow segment.  if the probability drops below this value the
        !           cell is no longer considered a possible flow source
        !        *******************************************************************
        READ(PARAMS,*) NO_FLOW_PROBABILITY
        !        *******************************************************************
        !         lava_gradient_weight is a parameter that determines how much the gradient
        !           between the edge of a flow and a neighboring point determines
        !           the probability of flow in that direction.
        !        *******************************************************************
        READ(PARAMS,*) LAVA_GRADIENT_WEIGHT
        !        *******************************************************************
        !         lava_duration_weight is a parameter that determines how rapidly a new cell
        !           diminishes in probability that it can be the source of a flow
        !           to a neighboring cell -- Presumably represents effects of cooling
        !           of emplaced lava
        !        *******************************************************************
        READ(PARAMS,*) LAVA_DURATION_WEIGHT
        !        *******************************************************************
        IBACK(2)=4
        IBACK(3)=5
        IBACK(4)=2
        IBACK(5)=3
        IBACK(6)=8
        IBACK(7)=9
        IBACK(8)=6
        IBACK(9)=7
        !       **********************************************************************
        !
        !       Part 4:  Read in cratering parameters
        !
        !       **********************************************************************
        READ(PARAMS,22) TEXTLINE
        WRITE(*,22) TEXTLINE
        READ(PARAMS,*) IMPACT_PROBABILITY, IFOLD, ISOFTCRATER
        IF (ISOFTCRATER > 0) THEN
            IS_REGOLITH_CRATER=.TRUE.
        ELSE
            IS_REGOLITH_CRATER=.FALSE.
        ENDIF
        READ(PARAMS,*) LARGE_CRATER_DEPTH_SCALE,LARGE_CRATER_DEPTH_EXPONENT, &
        LARGE_CRATER_RIM_SCALE,LARGE_CRATER_RIM_EXPONENT,TRANSITION_DIAMETER
        READ(PARAMS,*) SMALL_CRATER_DEPTH_SCALE,SMALL_CRATER_DEPTH_EXPONENT, &
        SMALL_CRATER_RIM_SCALE,SMALL_CRATER_RIM_EXPONENT
        READ(PARAMS,*) LARGE_CRATER_SHAPE_SCALE,LARGE_CRATER_SHAPE_EXPONENT,&
        SMALL_CRATER_SHAPE_SCALE,SMALL_CRATER_SHAPE_EXPONENT
        READ(PARAMS,*) CRATER_FREQUENCY_EXPONENT, FREQUENCY_CUTOFF_SCALING
        READ(PARAMS,*) SMALLEST_POSSIBLE_CRATER, SMALLEST_MODELED_CRATER,LARGEST_MODELED_CRATER
        READ(PARAMS,*) EJECTA_THICKNESS_VARIABILITY,NOISESD
        READ(PARAMS,*) INHERITANCE_PARAMETER, MAXIMUM_RIM_GRADIENT
        IF (SMALLEST_MODELED_CRATER < SMALLEST_POSSIBLE_CRATER) SMALLEST_MODELED_CRATER=SMALLEST_POSSIBLE_CRATER
        !       **********************************************************************
        !
        !       Part 5:  Read in ocean parameters
        !
        !       **********************************************************************
        READ(PARAMS,22) TEXTLINE
        WRITE(*,22) TEXTLINE
        READ(PARAMS,*) IOCEANVAR
        IF (MODEL_OCEAN_LEVEL.AND.(IOCEANVAR > 0)) THEN
            VARIABLE_OCEAN_ELEVATION=.TRUE.
            ALLOCATE(OCEAN_RECALCULATION_TIMES(2000),OCEAN_LEVELS(2000),STAT=IERROR)
        ELSE
            VARIABLE_OCEAN_ELEVATION=.FALSE.
        ENDIF
        !       **********************************************************************
        !
        !       Part 6:  Read in accretion parameters
        !
        !       **********************************************************************
        READ(PARAMS,22) TEXTLINE
        WRITE(*,22) TEXTLINE
        READ(PARAMS,*) ACCRETION_RATE
        READ(PARAMS,*) EXPOSECREEPUSE, RADIATION_USE, TOPEXPOSEUSE, INVERSEUSE, EXPOSESMOOTH
        IF (EXPOSECREEPUSE > 0) THEN
            EXPOSURE_DEPENDENT_CREEP=.TRUE.
            ALLOCATE(EXPOSE_RATE(MX,MY),SMOOTH_EXPOSE(MX,MY),STAT=IERROR)
        ELSE
            EXPOSURE_DEPENDENT_CREEP=.FALSE.
        ENDIF
        IF (RADIATION_USE > 0) THEN
            USE_SOLAR_EROSION=.TRUE.
        ELSE
            USE_SOLAR_EROSION=.FALSE.
        ENDIF
        IF (TOPEXPOSEUSE > 0) THEN
            USE_TOP_EXPOSURE=.TRUE.
        ELSE
            USE_TOP_EXPOSURE=.FALSE.
        ENDIF
        IF (INVERSEUSE > 0) THEN
            USE_INVERSE_EXPOSURE=.TRUE.
        ELSE
            USE_INVERSE_EXPOSURE=.FALSE.
        ENDIF
        IF (EXPOSESMOOTH > 0) THEN
            USE_EXPOSURE_SMOOTHING=.TRUE.
        ELSE
            USE_EXPOSURE_SMOOTHING=.FALSE.
        ENDIF
        READ(PARAMS,*) RAD_CONST,RAD_DUST_FACTOR,RAD_THRESH_CONVEXITY,RAD_DEPOSIT_RATE
        !        *******************************************************************
        !        *******************************************************************
        !          w r i t e  o u t  t h e  s i m u l a t i o n  p a r a m e t e r s
        !       !!! this needs work to write out all current simulation parameters
        !        *******************************************************************
        !        *******************************************************************
        !
        !        *******************************************************************
        !        ********************** boundary conditions ************************
        !        *******************************************************************
        WRITE(OUTHIST,500)
        500 FORMAT(' %%%%%%%%%%%%%%%MARS SIMULATION PROGRAM %%%%%%%%%%%%%%'/,&
        ' ******************  BOUNDARY CONDITIONS *******************')
        IF (NEW_SIMULATION) THEN
            WRITE(OUTHIST,487)
            487   FORMAT (' NEW START')
        ELSE
            WRITE(OUTHIST,486)
            486   FORMAT (' RESTART')
        ENDIF
        IF (FLUVIAL_AND_SLOPE_MODELING) THEN
            WRITE(OUTHIST,1487)
            1487 FORMAT(' MODELING FLUVIAL AND SLOPE PROCESSES')
        ENDIF
        IF (MODEL_IMPACT_CRATERING) THEN
            WRITE(OUTHIST,1486)
            1486 FORMAT(' MODELING IMPACT CRATERING')
        ENDIF
        IF (MODEL_LAVA_FLOWS) THEN
            WRITE(OUTHIST,1488)
            1488 FORMAT(' MODELING LAVA FLOWS')
        ENDIF
        IF (MODEL_EOLIAN_CHANGES) THEN
            WRITE(OUTHIST,1489)
            1489 FORMAT(' MODELING EOLIAN EROSION AND DEPOSITION')
        ENDIF
        IF (MODEL_OCEAN_LEVEL) THEN
            WRITE(OUTHIST,1490)
            1490 FORMAT(' MODELING VARIABLE OCEAN ELEVATION')
        ENDIF
        IF (MODEL_ACCRETION_AND_ABLATION) THEN
            WRITE(OUTHIST,1491)
            1491 FORMAT(' MODELING SLOPE-NORMAL ACCRETION AND ABLATION')
        ENDIF
        IF (MODEL_LAKE_EVAPORATION) THEN
            WRITE(*,8476)
            WRITE(OUTHIST,8476)
            8476 FORMAT(' MODELING LAKE EVAPORATION')
        ENDIF
        WRITE(OUTHIST,501) MX, MY, MZ,CELL_SIZE, VERTICAL_SCALING
        501 FORMAT(' MX=',I5,' MY=',I5, ' MZ=',I5,/,' HORIZONTAL CELL SIZE=',G12.5 &
            ,' VERTICAL SCALING SIZE=',G12.5)
        IF (.NOT.VARIABLE_EROSION_RATE) THEN
            WRITE(OUTHIST,503) BOUNDARY_LOWERING_RATE
            503 FORMAT(' CONSTANT SIMULATION PARAMETERS',/,' EROSION RATE=',G12.5)
        ELSE
            WRITE(OUTHIST,1492)
            1492 FORMAT(' USING TIME-VARYING PARAMETERS')
        ENDIF
        WRITE(OUTHIST,504) CHANNEL_TIMESTEP_SCALING,DEFAULT_CHANNEL_TIMESTEP &
           ,MAXIMUM_TIME_INCREMENT, MINUMUM_TIME_INCREMENT
        504 FORMAT(' FRACTION OF GRADIENT TIMESTEP IN MAX EROSION RATE='&
        ,G12.5,  &
        /,' MAXIMUM CHANNEL TIMESTEP=',G12.5,/, &
        ' MAXIMUM PERMITTED TIMESTEP=',G12.5,' MINIMUM TIMESTEP=',G12.5)
        WRITE(OUTHIST,934) SEDIMENT_TIMESTEP_SCALING,MASS_WASTING_TIMESTEP_SCALING
        934 FORMAT(' SEDIMENT AND SLOPE SCALING FACTORS FOR ',&
        'MAXIMUM EROSION RATE :' &
        ,/,2G12.5)
        IF (DO_DEBUGGING) THEN
            WRITE(OUTHIST,564) ICENT,JCENT
            564 FORMAT(' CENTER OF DEBUG WINDOW - I=',I5,' J=',I5)
        ENDIF
        IF (NON_ERODING_LOWER_BOUNDARY) THEN
            WRITE(OUTHIST,1495)
            1495 FORMAT(' THE LOWER BOUNDARY IS NOT ACTIVELY ERODED -' &
            ' IT IS CHANGED ONLY IN THE BOUNDARY CONDITIONS ROUTINE')
        ENDIF
        IF (HORIZONTAL_LOWER_BOUNDARY) THEN
            WRITE(OUTHIST,849)
            849   FORMAT(' LOWER BOUNDARY IS FORCED TO BE LEVEL')
        ELSE
            WRITE(OUTHIST,848)
            848   FORMAT(' AN UNEVEN LOWER BOUNDARY IS ALLOWED ')
        ENDIF
        IF (DO_FLOW_BOUNDARIES) THEN
            WRITE(OUTHIST,11490)
            11490 FORMAT(' THE OUTER BOUNDARIES ARE NON-EROSIONAL BUT FLOW AND SEDIMENT MAY EXIT')
        ENDIF
        IF (DO_EVENTS) THEN
            WRITE(OUTHIST,1493)
            1493 FORMAT(' DISCRETE EVENTS DURING SIMULATION ARE MODELED')
        ENDIF
        IF (DO_MORPHOMETRY) THEN
            WRITE(OUTHIST,1494)
            1494 FORMAT(' MORPHOMETRIC PARAMETERS ARE REPORTED')
        ENDIF
        IF (IS_X_PERIODIC) THEN
            WRITE(OUTHIST,1496)
            1496 FORMAT(' THE LATERAL BOUNDARIES ARE PERIODIC')
        ELSE
            WRITE(OUTHIST,1497)
            1497 FORMAT(' THE LATERAL BOUNDARIES ARE NON-FLUX')
        ENDIF
        IF (IS_Y_PERIODIC) THEN
            WRITE(OUTHIST,1498)
            1498 FORMAT(' THE TOP AND BOTTOM BOUNDARIES ARE PERIODIC')
        ELSE
            WRITE(OUTHIST,1499)
            1499 FORMAT(' THE TOP BOUNDARY IS NON-FLUX')
        ENDIF
        !        *******************************************************************
        !        ********************** output and recalculation timing *************
        !        *******************************************************************
        WRITE(OUTHIST,514) MAXIMUM_ITERATION, MAXIMUM_SIMULATION_TIME
        514 FORMAT(' ***************  OUTPUT AND RECALCULATION TIMING ******'&
        ,/,' TOTAL ITERATIONS=',I5,' MAXIMUM SIMULATION ELAPSED TIME=',G12.5)
        WRITE(OUTHIST,1500) STARTING_ITERATION,PRESENT_TIME
        1500 FORMAT(' STARTING ITERATION NUMBER=',I6,' STARTING TIME=',G12.5)
        WRITE(OUTHIST,515)ELEVATION_PRINT_INTERVAL,OUTPUT_PRINT_INTERVAL
        515 FORMAT(' ELEVATION OUTPUT INTERVAL=',I5, &
        ' PRINT INTERVAL=',I5)
        WRITE(OUTHIST,721)WRITE_CHANGE_INTERVAL
        721 FORMAT(' DIRECTION AND GRADIENT CHANGE PRINT INTERVAL=',I5)
        WRITE(OUTHIST,516)RECALCULATE_GRADIENT_INTERVAL
        516 FORMAT(' GRADIENT REEVALUATE INTERVAL=',I5)
        WRITE(OUTHIST,517)STARTING_ITERATION,PRESENT_TIME,CRITICAL_SOURCE_DIVERGENCE
        517 FORMAT(' STARTING ITERATION=',I6,' STARTING TIME=',G12.5,/, &
        ' MINIMUM DIVERGENCE FOR CHANNEL PLOTS=',G12.5)
        WRITE(OUTHIST,518)(WRITETYPE(I),I=1,10)
        518 FORMAT(' WRITE INDICES FOR OUTPUT =',10I2)
        WRITE(OUTHIST,519)
        519 FORMAT('      1 = ELEVATION')
        WRITE(OUTHIST,520)
        520 FORMAT('      2 = SEDIMENT THICKNESS')
        WRITE(OUTHIST,521)
        521 FORMAT('      3 = CHANNEL GRADIENT')
        WRITE(OUTHIST,523)
        523 FORMAT('      4 = SLOPE EROSION RATE')
        WRITE(OUTHIST,522)
        522 FORMAT('      5 = ADJUSTED CHANNEL EROSION RATE')
        WRITE(OUTHIST,551)
        551 FORMAT('      6 = OVERALL EROSION RATE')
        WRITE(OUTHIST,553)
        553 FORMAT('      7 = CONVERGENCE-DIVERGENCE')
        WRITE(OUTHIST,552)
        552 FORMAT('      8 = PROFILE CURVATURE')
        WRITE(OUTHIST,554)
        554 FORMAT('      9 = PLANFORM CURVATURE')
        WRITE(OUTHIST,558)
        558 FORMAT('     10 = LN(AREA)')
        WRITE(OUTHIST,559)
        559 FORMAT('     11 = LN(AREA/GRADIENT',/)
        WRITE(OUTHIST,784) IMAGE_OUTPUT_INTERVAL,DIVERGEINTERVAL
        784   FORMAT(' INTERVAL FOR PRINTING IMAGE FILE=',I5,/, &
        ' INTERVAL FOR CALCULATING DIVERGENCE=',I5,/, &
        ' **************** FLOW ROUTING AND DRAINAGE AREA *************')
        !        *******************************************************************
        !        ********************* flow routing and drainage area ***************
        !        *******************************************************************
        IF (.NOT.MODEL_LAKE_EVAPORATION) THEN
            IF (COMPLETE_RUNOFF) THEN
                WRITE(OUTHIST,547)
                547 FORMAT(' DEPRESSIONS ARE COMPLETELY FILLED AND SPILL')
            ELSE
                WRITE(OUTHIST,548)
                548 FORMAT(' DEPRESSIONS ARE DRY')
            ENDIF
        ENDIF
        IF (USE_RANDOM_FLOWROUTING) THEN
            WRITE(OUTHIST,495)
            495   FORMAT(' RANDOM DOWNSTREAM FLOW ROUTING')
        ELSE
            WRITE(OUTHIST,496)
            496   FORMAT(' DETERMINISTIC D8 DOWNSTREAM FLOW ROUTING')
        ENDIF
        WRITE(OUTHIST,494) AREAFACTOR, RAINDEPTH,RAINSTD
        494   FORMAT(' AREA CORRECTION FACTOR=',G12.5, 'RAINFALL MEAN DEPTH=',G12.5, &
              ' STANDARD DEVIATION OF RAINFALL=',G12.5)
        IF (DIVERGENCE_DEPENDENT_RUNOFF) THEN
            WRITE(OUTHIST,783) MEAN_CONVERGENCE_RUNOFF,HIGH_CONVERGENCE_RUNOFF,LOW_CONVERGENCE_RUNOFF,  &
            RUNOFF_CONSTANT_1,RUNOFF_CONSTANT_2,RUNOFF_CONSTANT_3
            783   FORMAT(' PROGRAM USES VARIABLE WATER YIELD BASED ON SLOPE' &
            ,' CONVERGENCE',/,' CONVERGENCE FOR MEAN WATER YIELD=',G12.5,&
            ' CONVERGENCE FOR 90% YIELD=',G12.5,' LOW WATER YIELD='&
            ,G12.5,/, &
            ' DERIVATIVE WATER YIELD CONSTANTS 1-3=',3(G12.5,' '))
        ENDIF
        IF (MODEL_PELAGIC_DEPOSITION) THEN
           WRITE(OUTHIST,1504) WASHLOAD_FRACTION,PELAGICCREEP
           1504 FORMAT(' PELAGIC DEPOSITION IS MODELED',/, &
           '     FRACTION OF SEDIMENT THAT IS DELIVERED OFFSHORE=',G12.5, &
           '     COEFFICIENT FOR PELAGIC SEDIMENT DIFFUSIVITY=',G12.5)
        ENDIF
        IF (HAVE_INFLUENT_RIVERS) THEN
           WRITE(OUTHIST,1505)
           1505 FORMAT(' MODELING INFLOW OF WATER AND SEDIMENT FROM ONE OR MORE RIVERS AT OUTER BOUNDARIES')
        ENDIF
        !        *******************************************************************
        !        ********************** bedrock channel parameters **************
        !        *******************************************************************
        IF (FLUVIAL_AND_SLOPE_MODELING.AND.DO_FLUVIAL_DETACHMENT) THEN
            WRITE(OUTHIST,505) BEDEXPLICIT
            505 FORMAT(' **********  BEDROCK CHANNEL PARAMETERS ************' &
            ,/,' INDEX FOR IMPLICIT (0) OR EXPLICIT (1) BEDROCK CHANNELS'&
            ,': ',I5 )
            IF (USE_SKLAR_BED_ABRASION.OR.USE_WHIPPLE_BED_ABRASION) THEN
            ELSE
               WRITE(OUTHIST,1506)
               1506 FORMAT(' USING SHEAR-STRESS DEPENDENT BEDROCK CHANNEL EROSION')
            ENDIF
            IF (USE_SKLAR_BED_ABRASION) THEN
               WRITE (OUTHIST,1508) ROCK_TENSILE_STRENGTH,SKLAR_MULT
               1508 FORMAT(' USING SKLAR BEDROCK EROSION, TENSILE STRENGTH=',G12.5, &
               ' RATE MULTIPLICATIVE FACTOR=',G12.5)
            ENDIF
            IF (USE_WHIPPLE_BED_ABRASION) THEN
               WRITE(OUTHIST,1509)
               1509 FORMAT(' USING WHIPPLE BEDROCK EROSION')
            ENDIF
            WRITE(OUTHIST,1507) BEDROCK_DISCHARGE_EXPONENT,BEDROCK_GRADIENT_EXPONENT
            1507 FORMAT(' BEDROCK_DISCHARGE_EXPONENT=',G12.5,' BEDROCK_GRADIENT_EXPONENT=',G12.5)
            IF (.NOT.VARIABLE_VEGETATION_RESISTANCE) THEN
                WRITE(OUTHIST,11508) DETACHMENT_CRITICAL_SHEAR
                11508 FORMAT(' CRITICAL SHEAR FOR CHANNEL EROSION=',G12.5)
            ENDIF
            WRITE(OUTHIST,506) BEDROCK_ERODIBILITY,DISCHARGE_COEFFICIENT
            506 FORMAT(' MATERIAL EROSION CONSTANT=',G12.5, &
               ' WATER EROSION CONSTANT=',G12.5)
            WRITE(OUTHIST,507) REGOLITH_ERODIBILITY_FACTOR,REGOLITH_CRITICAL_SHEAR_FACTOR
            507 FORMAT(' ERODIBILITY FACTOR FOR SLOPE DEBRIS', &
            ' IN CHANNEL ( >= 1.0): ',G12.5,/,   &
            ' CRITICAL SHEAR STRESS FACTOR FOR REGOLITH=',G12.5)
            WRITE(OUTHIST,731) CHANNEL_WIDTH_CONSTANT,CHANNEL_WIDTH_EXPONENT
            731 FORMAT(' CHANNEL WIDTH CONSTANT=',G12.5, &
            ' WIDTH EXPONENT=',G12.5)
            IF (VARIABLE_VEGETATION_RESISTANCE) THEN
                WRITE(OUTHIST,1510)
                1510 FORMAT(' AREA-DEPENDENT CRITICAL SHEAR STRESS DUE TO VEGETATION')
                WRITE(OUTHIST,1511) VEGETATION_UPLAND_RESISTANCE,VEGETATION_CHANNEL_RESISTANCE
                1511 FORMAT(' CRITICAL SHEAR STRESS ON UPLANDS=',G12.5, &
                ' CRITICAL SHEAR STRESS FOR LARGE AREA=',G12.5)
                WRITE(OUTHIST,1512) VEGETATION_AREA_MINIMUM,VEGETATION_AREA_MAXIMUM
                1512 FORMAT(' MAXIMUM AREA WITH FULL VEGETATION CRITICAL SHEAR=',G12.5 &
                ,/,' MINIMUM AREA FOR LARGE AREA CRITICAL SHEAR=',G12.5)
            ENDIF
            WRITE(OUTHIST,1520) MANNING
            1520 FORMAT(' MANNING RESISTANCE COEFFICENT=',G12.5)
        ELSE
            WRITE(OUTHIST,11506)
            11506 FORMAT('FLUVIAL DETACHMENT IS NOT MODELED')
        ENDIF
        WRITE(OUTHIST,1521) GRAVITY
        1521 FORMAT(' GRAVITATIONAL CONSTANT=',G12.5)
        !        *******************************************************************
        !        ********************** alluvial channel parameters *************
        !        *******************************************************************
        WRITE(OUTHIST,565)
        565   FORMAT(' ******************* ALLUVIAL CHANNEL PARAMETERS *******')
        IF (FLUVIAL_AND_SLOPE_MODELING.AND.DO_SEDIMENT_TRANSPORT) THEN
            IF (WRITE_SEDIMENT_DIAGNOSTICS) THEN
               WRITE(OUTHIST,4422)
                4422     FORMAT(' DOING SEDIMENT DEBUGGING')
            ELSE
                WRITE(OUTHIST,4423)
                4423     FORMAT(' NOT DOING SEDIMENT DEBUGGING')
            ENDIF
            IF (DO_SEDIMENT_ROUTING) THEN
                WRITE(OUTHIST,11507)
                11507 FORMAT(' MODELING DOWNSTREAM ROUTING OF SEDIMENT')
            ENDIF
            IF (DO_SEDIMENT_DIFFUSION) THEN
                WRITE(OUTHIST,21508)
                21508 FORMAT(' MODELING SEDIMENT TRANSPORT DIVERGENCE')
            ENDIF
            IF ((.NOT.IS_Y_PERIODIC).AND.NO_FLUX_LOWER_BOUNDARY) THEN
                WRITE(OUTHIST,11509)
                 11509 FORMAT(' THE LOWER BOUNDARY IS A NON-FLUX BARRIER TO SEDIMENT TRANSPORT')
            ELSE
                IF (.NOT.IS_Y_PERIODIC) THEN
                     WRITE(OUTHIST,11510)
                    11510 FORMAT(' SEDIMENT EXITS AT LOWER BOUNDARY')
                ENDIF
            ENDIF
            IF (DO_ALLUVIAL_SMOOTHING) THEN
                WRITE(OUTHIST,11511)
                11511 FORMAT(' SMOOTHING OF DEPOSITED SEDIMENT SURFACE')
                IF (DO_ALLUVIAL_REROUTING) THEN
                    WRITE(OUTHIST,11512)
                    11512 FORMAT(' FLOW DIRECTIONS ON SEDIMENT SURFACE CHANGES WITHIN AN ITEATION STEP')
                ENDIF
            ENDIF
            IF (IREADALLUV > 0) THEN
                WRITE(OUTHIST,567)
                567       FORMAT(' ALLUVIAL CHANNEL AREAS READ IN FROM FILE')
            ENDIF
            IF (USE_AN_OCEAN) THEN
                WRITE(OUTHIST,568) OCEAN_ELEVATION, DELTA_FORESET_GRADIENT
                568       FORMAT(' STANDING WATER LEVEL AT ELEVATION=',G12.5 &
                ,' WITH FORSET BED GRADIENT=',G12.5)
            ENDIF
            WRITE(OUTHIST,560) SEDIMENT_1_EXPONENT,SEDIMENT_2_EXPONENT
            560 FORMAT(' EXPONENTS FOR DRAINAGE AREA IN SEDIMENT TRANSPORT RATE=', &
            2G12.5)
            WRITE(OUTHIST,561) SEDIMENT_GRADIENT_EXPONENT,SEDIMENT_TRANSPORT_EXPONENT
            561 FORMAT(' EXPONENT FOR GRADIENT IN SEDIMENT TRANSPORT RATE=',&
            G12.5,/,' TRANSPORT EXPONENT IN SEDIMENT TRANSPORT RATE='&
            ,G12.5)
            WRITE(OUTHIST,1514) TRANSPORTFACTOR,SEDIMENT_POROSITY,SEDIMENT_SPECIFIC_GRAVITY,GRAIN_SIZE
            1514 FORMAT(' DIMENSIONLESS SHEAR STRESS MULTIPLIER=',G12.5,' SEDIMENT POROSITY=',G12.5 &
            ,/,' SEDIMENT SPECIFIC GRAVITY=',G12.5,' SEDIMENT GRAIN SIZE=',G12.5)
            WRITE(OUTHIST,562) SEDIMENT_CONSTANT,TRANSPORT_CRITICAL_DIM_SHEAR
            562 FORMAT(' MULTIPLICATIVE CONSTANT IN SEDIMENT TRANSPORT RATE=', &
            G12.5,/,' CRITICAL DIMENSIONLESS SHEAR STRESS',G12.5)
            WRITE(OUTHIST,1515) SEDIMENT_DISCHARGE_FACTOR
            1515 FORMAT(' DISCHARGE MULTIPLIER IN SEDIMENT TRANSPORT=',G12.5)
            WRITE(OUTHIST,1513) EFFECTIVE_DISCHARGE_RATIO,FLOW_FRACTION
            1513 FORMAT(' EFFECTIVE DISCHARGE FOR SEDIMENT TRANSPORT RELATIVE TO MEAN ANNUAL FLOOD=' &
            ,G12.5,/,' FRACTION OF YEAR WITH EFFECTIVE SEDIMENT TRANSPORT DISCHARGE=',G12.5)
            IF (STICKY_SEDIMENT_ROUTING) THEN
                WRITE(OUTHIST,497) STICKY_ROUTING_CRITICAL_VALUE
                497   FORMAT(' ALLUVIAL AREAS USE STICKY DIRECTIONS,',' WITH 0.5', &
                '    PROBABILITY OF CHANGE AT GRADIENT RATIO OF ',G12.5)
            ENDIF
        ELSE
            WRITE(OUTHIST,11513)
            11513 FORMAT(' FLUVIAL SEDIMENT TRANSPORT IS NOT MODELED')
        ENDIF
        !        *******************************************************************
        !        ****************  regolith weathering parameters  ******************
        !        *******************************************************************
        WRITE(OUTHIST,788) ROCK_WEATHERING_RATE,WEATHER_DECAY_RATE,INITIAL_REGOLITH_THICKNESS
        788   FORMAT(' ******************* REGOLITH WEATHERING PARAMETERS ****' &
        ,/,' MAXIMUM WEATHERING RATE=',G12.5,' DEPTH TO HALF RATE=', &
        G12.5,' INITIAL REGOLITH THICKNESS=',G12.5)
        IF (TWO_TERM_WEATHERING) THEN
            WRITE(OUTHIST,779) WEATHERING_TERM_2,WEATHERING_DECAY_2
            779    FORMAT(' TWO EXPONENTIAL WEATHERING RATE USED',/,&
            ' RATE CONSTANT 2=',G12.5,' DECAY CONSTANT 2=',G12.5)
        ENDIF
        WRITE(OUTHIST,787) CRITICAL_BEDROCK_GRADIENT,WEATHER_DIVERGENCE
        787   FORMAT(' GRADIENT VALUE GIVING DOUBLE WEATERING RATE=',G12.5,  &
        /,' DIVERGENCE VALUE GIVING DOUBLE WEATHERING RATE=', &
        G12.5)
        IF (WEATHER_DIVERGENCE > 0.0) THEN
            WEATHER_DIVERGENCE  = DLOG(2.0D0) / WEATHER_DIVERGENCE
        ENDIF
        WRITE(OUTHIST,786) CRITICAL_BEDROCK_GRADIENT,WEATHER_DIVERGENCE,READREGOLITH
        786   FORMAT(' SCALED WEATHER GRADIENT=',G12.5,' SCALED DIVERGENCE', &
        ' VALUE=',G12.5,/,' READ REGOLITH THICKNESS FROM FILE', &
        ' (0=NO, 1=YES): ',I5)
        IF (SEEPAGE_WEATHERING) THEN
            WRITE(OUTHIST,1516) SEEPAGE_WEATHERING_SCALING,SEEPAGE_WEATHERING_EXPONENT
            1516 FORMAT (' ROCK WEATHERING BY EMERGENT SEEPAGE IS MODELED' &
            ,/,' SEEPAGE WEATHERING RATE FACTOR=',G12.5,' SEEPAGE WEATHERING EXPONENT=',G12.5)
        ENDIF
        !        *******************************************************************
        !        *****************  mass wasting parameters  ************************
        !        *******************************************************************
        IF (DO_MODEL_SLOPES) THEN
            WRITE(OUTHIST,508) SLOPE_DIFFUSIVITY
            508 FORMAT(' *************** MASS WASTING PARAMETERS ***************'&
            ,/,' SLOPE EROSION IS DIRECTLY MODELED',/, &
            ' MASS MOVEMENT SLOPE SCALING=',G12.5)
            WRITE(OUTHIST,569) ALVCREEPFAC
            569   FORMAT(' RELATIVE RATE OF CREEP ON ALLUVIAL COVER =',G12.5)
            IF (USE_ROERING_MASS_WASTING) THEN
                WRITE(OUTHIST,11520)
                11520 FORMAT(' MASS WASTING USES ROERING PARAMETERIZATION')
            ELSE
                WRITE(OUTHIST,11521)
                11521 FORMAT(' HOWARD MASS WASTING PARAMETERIZATION USED')
            ENDIF
            IF (USE_CRITICAL_SLOPE_CRADIENT) WRITE(OUTHIST,513)CRITICAL_SLOPE_GRADIENT, &
            SLOPE_GRADIENT_EXPONENT,SLOPE_FAILURE_DIFFUSIVITY
            513 FORMAT(' CRITICAL SLOPE GRADIENT=',G12.5,&
            ' CRITICAL GRADIENT POWER=',G12.5,/, &
            ' CRITICAL GRADIENT MULTIPLIER - SLOPE_FAILURE_DIFFUSIVITY=',G12.5)
        ELSE
            WRITE(OUTHIST,438)
            438   FORMAT(' *******  SLOPE EROSION IS NOT DIRECTLY MODELED *******')
        ENDIF
        !        *******************************************************************
        !        *********************** bistable erosion  **************************
        !        *******************************************************************
        IF (BISTABLE_FLUVIAL_EROSION) THEN
            IF (USE_BISTABLE_BEDROCK) THEN
                WRITE(OUTHIST,1522)
                1522 FORMAT(' BISTABLE BEDROCK IS USED')
            ELSE
                WRITE(OUTHIST,1523)
                1523 FORMAT(' BISTABLE EROSION DOES NOT USE BISTABLE BEDROCK')
            ENDIF
            WRITE(OUTHIST,780) LOW_EROSION_THRESHOLD,HIGH_EROSION_THRESHOLD,EROSION_RATE_CHANGE_LAG, &
            BISTABLE_CRITICAL_SHEAR,BISTABLE_RUNOFF_FACTOR,BISTABLE_BEDROCK_ERODIBILITY
            780    FORMAT(' ************** BISTABLE EROSION *********************' &
            ,/,' USING EROSION-RATE DEPENDENT PARAMETERS',/, &
            ' LOW_EROSION_THRESHOLD=',G12.5,' HIGH_EROSION_THRESHOLD=',G12.5,' EROSION_RATE_CHANGE_LAG=', &
            G12.5,/,' BISTABLE_CRITICAL_SHEAR=',G12.5,' BISTABLE_RUNOFF_FACTOR=',G12.5,&
            ' BISTABLE_BEDROCK_ERODIBILITY=',G12.5)
        ENDIF
        !        *******************************************************************
        !        ******** random variability of simulation parameters    ***************
        !        *******************************************************************
        IF (RANDOM_CRITICAL_SHEAR) THEN
            WRITE(OUTHIST,498) CRITICAL_SHEAR_VARIABILITY,OMEGA_WEIGHT
            498   FORMAT('  *************** RANDOM VARIABILITY *******************'&
            ,/,' TEMPORALLY VARIABLE FLUVIAL EROSION THRESHOLD,' &
            ,' VARIANCE=',  &
            G12.5,' OMEGA=',G12.5)
        ENDIF
        IF (USE_RANDOM_DISCHARGE) THEN
            WRITE(OUTHIST,499) DISCHARGE_COEFF_VARIATION,OMEGA_WEIGHT
            499   FORMAT('  *************** RANDOM VARIABILITY *******************'&
            ,/,' TEMPORALLY VARIABLE FLUVIAL DISCHARGE, VARIANCE=',&
            G12.5,' OMEGA=',G12.5)
        ENDIF
        !        *******************************************************************
        !        ***************  3-d spatially variable bedrock resistance  ********
        !        *******************************************************************
        USE_3D_SLOPE_RESISTANCE = .FALSE.
        IF (SLOPEVARUSE > 0) THEN
            WRITE(OUTHIST,524)
            524 FORMAT(' *************** 3-D VARIABLE BEDROCK RESISTANCE ******'&
            ,/,' PROGRAM USES VARIABLE SLOPE RESISTANCE')
            WRITE(OUTHIST,525)DIFFUSIVITY_VARIABILITY
            525 FORMAT('   WITH SCALING VARIABLE=',G12.5,/)
            USE_3D_SLOPE_RESISTANCE = .TRUE.
        ENDIF
        IF (SCALE_3D_ROCK_ERODIBILITY) THEN
            WRITE(OUTHIST,805)
        ELSE
            WRITE(OUTHIST,806)
        ENDIF
        805   FORMAT(' INPUT RESISTANCE ARE RANDOM DEVIATES')
        806   FORMAT(' INPUT RESISTANCE ARE ABSOLUTE VALUES')
        WRITE(OUTHIST,785) VERTICAL_RESISTANCE_SCALING
        785   FORMAT(' RATIO OF VERTICAL TO HORIZONTAL SCALES FOR RESISTANCE=', &
        G12.5)
        IF (RESISTANT_SURFACE_LAYER) THEN
            WRITE(OUTHIST,1524) SURFACE_LAYER_THICKNESS, SURFACE_LAYER_RESISTANCE
            1524 FORMAT(' A RESISTANT SURFACE LAYER IS USED: LAYER THICKNESS=',G12.5, &
            ' RESISTANCE OF SURFACE LAYER AS COMPARED TO NORMAL BEDROCK=',G12.5)
        ENDIF
        !        *******************************************************************
        !        **************** rock and surface deformation *********************
        !        *******************************************************************
        IF (DEFORMUSE > 0) THEN
            DO_ROCK_DEFORMATION=.TRUE.
            WRITE(OUTHIST,776) DEFORMSCALE
            776   FORMAT(' ********* ROCK AND SURFACE DEFORMATION ****************' &
            ,/,' ROCK AND SURFACE DEFORMATION IS BEING USED WITH',/,  &
            ' DEFORMATION RATE SCALING OF ',G12.5,' WITH UPLIFT POSITIVE')
        ELSE
            DO_ROCK_DEFORMATION=.FALSE.
        ENDIF
        !        *******************************************************************
        !        **************** groundwater flow**********************************
        !        *******************************************************************
        IF (MODEL_GROUNDWATER) THEN
             WRITE(OUTHIST,1525)
             1525 FORMAT(' **************** GROUNDWATER FLOW PARAMETERS ******************')
             WRITE(OUTHIST,1526) SEEPAGE_ITERATION_INTERVAL
             1526 FORMAT(' GROUNDWATER FLOW CALCULATED EVERY ',I6,' ITERATIONS')
             IF (EXPONENTIAL_PERMEABILITY_DECAY) THEN
                 WRITE(OUTHIST,1527)
                 1527 FORMAT(' PERMEABILITY DECREASES EXPONENTIALLY WITH DEPTH BELOW SURFACE')
             ELSE
                 WRITE(OUTHIST,1528)
                 1528 FORMAT(' PERMEABILITY IS CONSTANT WITH DEPTH')
             ENDIF
             IF (PERMEABILITY_RESCALING) THEN
                 WRITE(OUTHIST,1529)
                 1529 FORMAT(' PERMEABILITY IS REFERENCED TO CURRENT LAND SURFACE ELEVATION')
             ELSE
                 WRITE(OUTHIST,1530)
                 1530 FORMAT(' PERMEABILITY IS REFERENCED TO STARTING LAND SURFACE ELEVATION')
             ENDIF
             IF (SHOW_GROUNDWATER_CALCULATIONS) THEN
                 WRITE(OUTHIST,1531)
                 1531 FORMAT(' DETAILS OF GROUNDWATER CALCULATIONS ARE PRINTED OUT')
             ENDIF
             WRITE(OUTHIST,1532) YEARLY_RECHARGE,VISCOSITY,DARCIES,GROUNDWATER_DEPTH_SCALE
             1532 FORMAT(' YEARLY RECHARGE=',G12.5,' GROUNDWATER VISCOSITY=',G12.5, &
             ' GROUNDWATER DEPTH SCALE=',G12.5)
             WRITE(OUTHIST,1533) GROUNDWATER_FLOW_FRACTION, INITIAL_GROUNDWATER_DEPTH, EPOWER
             1533 FORMAT(' FRACTION OF SURFACE WATER FLOW THAT IS GROUNDWATER=',G12.5,/, &
             ' INITIAL GROUNDWATER DEPTH BENEATH SURFACE=',G12.5,' RATE OF EXPOENTIAL DECAY=',G12.5)
             WRITE(*,4215) HYDRAULIC_CONDUCTIVITY,GROUNDWATER_SCALE_FACTOR,GROUNDWATER_RECHARGE_RATE
             WRITE(OUTHIST,14215) HYDRAULIC_CONDUCTIVITY,GROUNDWATER_SCALE_FACTOR,GROUNDWATER_RECHARGE_RATE
             14215  FORMAT(' HYDRAULIC CONDUCTIVITY=',G12.5,' GROUNDWATER_SCALE_FACTOR=',G12.5, &
             ' GROUNDWATER_RECHARGE_RATE=',G12.5)
             WRITE(OUTHIST,1534) MAXIMUM_GROUNDWATER_ITERATIONS, MAXIMUM_GOUNDWATER_ERROR,GROUNDWATER_RELAXATION_CONST
             1534 FORMAT(' MAXIMUM NUMBER OF ITERATIONS FOR FLOW CONVERGENCE=',I6,/, &
             ' MAXIMUM PERMITTED ERROR=',G12.5,' CALCULATION CONVERGENCE FACTOR=',G12.5)
             IF (GROUNDWATER_FLOW_FRACTION > 0.5) THEN
                 WRITE(OUTHIST,1535)
                 1535 FORMAT(' SEDIMENT TRANSPORT PARAMETERS RECALCULATED FOR GROUNDWATER FLOW')
                 WRITE(OUTHIST,32001) SEDIMENT_DISCHARGE_FACTOR,SEDIMENT_CONSTANT,SKLAR_MULT
                 WRITE(*,42001) SEDIMENT_DISCHARGE_FACTOR,SEDIMENT_CONSTANT,SKLAR_MULT
                 42001 FORMAT('CALCULATED SEDIMENT_DISCHARGE_FACTOR=',G12.5,' SEDIMENT_CONSTANT=',G12.5 &
                 ,' SLARMULT=',G12.5)
             ENDIF
             IF (USE_GROUNDWATER_FLOW) THEN
                 WRITE(OUTHIST,1536)
                 1536 FORMAT (' SEEPAGE MODELING USES GROUNDWATER FLUX (SUBSURFACE DISCHARGE)')
             ELSE
                 WRITE(OUTHIST,1537)
                 1537 FORMAT (' SEEPAGE MODELING USES GROUNDWATER FLUX DIVERGENCE (FLOW TO SURFACE)')
             ENDIF
             IF (SEEPAGE_AVERAGING) THEN
                 WRITE(OUTHIST,1538)
                 1538 FORMAT(' GROUNDWATER FLOW MEASUREMENTS ARE SPATIALLY AVERAGED')
             ENDIF
        ENDIF
        ! **************************************************************************
        !    Ocean parameters
        ! **************************************************************************
        IF (VARIABLE_OCEAN_ELEVATION) THEN
            WRITE(OUTHIST,1560)
            1560 FORMAT('******************MODELING TIME-VARYING OCEAN LEVELS*************')
        ENDIF
        !        *******************************************************************
        !         Set up initial rock resistance to unity and regolith thickness to
        !           initial_regolith_thickness
        !        *******************************************************************
        IF ( NEW_SIMULATION ) THEN
            DO  J = 1, MY
                DO  I = 1, MX
                    IF (FLUVIAL_AND_SLOPE_MODELING) THEN
                        IF (RESISTANT_SURFACE_LAYER) THEN
                            RELATIVE_RESISTANCE(I,J)=SURFACE_LAYER_RESISTANCE
                        ELSE
                            RELATIVE_RESISTANCE(I,J)=1.0
                        ENDIF
                    ENDIF
                    IF (FLUVIAL_AND_SLOPE_MODELING) REGOLITH(I,J)=INITIAL_REGOLITH_THICKNESS
                ENDDO
            ENDDO
        ENDIF
        THEMINIMUM=1.0E+15
        !        *******************************************************************
        !         Set initial sediment yield to zero, is_sediment_covered to false,
        !           is_rock_surface to false, accelerated_erosion to false, and last erosion rate to zero
        !        *******************************************************************
        IF (NEW_SIMULATION) THEN
                READ(INDATA,*)MX1,MY1
            IF ((MX1 /= MX).OR.(MY1 /= MY)) THEN
                WRITE(*,1121) MX, MX1,MY,MY1
                1121      FORMAT(' INCOMPATIBLE INPUT DIMENSIONS, MX=',I5,' MX1=',I5, &
                '  MY=',I5,' MY1=',I5)
                STOP
            ENDIF
            IF (FLUVIAL_AND_SLOPE_MODELING) THEN
                SEDIMENT_YIELD=0.0
                ACCELERATED_EROSION=.FALSE.
                DO_ACCELERATED_EROSION=.FALSE.
                PREVIOUS_ELEVATION=0.0
                DEFORMATION=0.0
            ENDIF
            IF (FLUVIAL_AND_SLOPE_MODELING) THEN
                IS_SEDIMENT_COVERED=.FALSE.
            ENDIF
            IF (FLUVIAL_AND_SLOPE_MODELING) THEN
                IS_ROCK_SURFACE=.FALSE.
            ENDIF
        ENDIF
        !        *********************************************************************
        !         Write out eolian parameters
        !        *********************************************************************
        IF (MODEL_EOLIAN_CHANGES) THEN
            WRITE(OUTHIST,2094)
            2094 FORMAT('*********************EOLIAN PARAMETERS***************************')
            WRITE(OUTHIST,1539) EOLIAN_EVENT_PROBABILITY,EOLIAN_TIME_INCREMENT
            1539 FORMAT(' EOLIAN EVENT PROBABILITY=',G12.5,' EOLIAN TIME INCREMENT=',G12.5)
            IF (USE_TOTAL_EXPOSURE) THEN
                WRITE(OUTHIST,1540)
                1540 FORMAT(' ALL CELLS WITHIN WINDOW USED FOR WEIGHTING')
            ELSE
                WRITE(OUTHIST,1541)
                1541 FORMAT(' ONLY VISIBLE CELLS WITHIN WINDOW USED FOR WEIGHTING')
            ENDIF
            IF (DEFAULT_EOLIAN_PROCESS) THEN
                WRITE(OUTHIST,1542)
                1542 FORMAT(' DEPOSITION AND EROSION ARE NORMAL TO THE LAND SURFACE')
            ELSE
                WRITE(OUTHIST,1543)
                1543 FORMAT(' DEPOSITION AND EROSION OCCUR IN VERTICAL DIRECTION')
            ENDIF
            WRITE(OUTHIST,2100) MINIMUM_EOLIAN_DEPOSIT_RATE,MAXIMUM_EOLIAN_DEPOSIT_RATE
            2100  FORMAT(' INPUT PARAMETERS: MINIMUM_EOLIAN_DEPOSIT_RATE=',G12.5, &
            ' MAXIMUM_EOLIAN_DEPOSIT_RATE=',G12.5)
            SELECT CASE (ICHOSE)
                CASE (1)
                    WRITE(OUTHIST,1544) EXPOSURE_10_PERCENT, EXPOSURE_90_PERCENT
                    1544 FORMAT(' 10TH AND 90TH PERCENT EXPOSURE INDEX VALUES: ', G12.5,' ',G12.5)
                CASE (2)
                    WRITE(OUTHIST,1545) EXPOSURE_50_PERCENT, EXPOSURE_90_PERCENT
                    1545 FORMAT(' 50TH AND 90TH PERCENT EXPOSURE INDEX VALUES: ', G12.5,' ',G12.5)
                CASE (3)
                    WRITE(OUTHIST,1546) ZERO_PERCENT_EXPOSURE,EXPOSURE_90_PERCENT
                    1546 FORMAT(' 0TH AND 90TH PERCENT EXPOSURE INDEX VALUES: ',G12.5,' ',G12.5)
                CASE(4)
                    WRITE(OUTHIST,1547) RATE0, EXPOSURE_50_PERCENT
                    1547 FORMAT(' RATE AT 0.0 EXPOSURE AND 90TH PERCENT EXPOSURE INDEX VALUES: ',G12.5,' ',G12.5)
                CASE DEFAULT
                    WRITE(OUTHIST,1548)
                    1548 FORMAT(' INVALID CHOICE OF INPUT PARAMETERS')
            END SELECT
        ENDIF
        IF (MODEL_EOLIAN_CHANGES.OR.MODEL_ACCRETION_AND_ABLATION.OR.USE_SOLAR_EROSION) THEN
            WRITE(OUTHIST,1549) DISTANCE_DECAY_FACTOR,WEIGHTING_CALCULATION_DISTANCE
            1549 FORMAT(' DECAY FACTOR FOR DISTANCE WEIGHTING=',G12.5,' MAXIMUM CALCULATION DISTANCE=',I5)
            WRITE(OUTHIST,2095) MAXIMUM_WEIGHT_DISTANCE,WEIGHTING_DECAY_FACTOR
            2095    FORMAT(' MAXIMUM RANGE=',I5, &
            ' WEIGHTING_DECAY_FACTOR=',G12.5)
            DO  I=1,MX
                IF (WEIGHTING_DECAY_FACTOR /= 0.0) THEN
                    WEIGHTX(I)=DEXP(WEIGHTING_DECAY_FACTOR*I)
                    WEIGHTD(I)=DEXP(WEIGHTING_DECAY_FACTOR*I*1.414)/(1.414)
                ELSE
                    WEIGHTX(I)=1.0
                    WEIGHTD(I)=1.0/(1.414)
                ENDIF
            ENDDO
        ENDIF
        !       ************************************************************************
        !        Write out lavaflow parameters
        !       ************************************************************************
        IF (MODEL_LAVA_FLOWS) THEN
            WRITE(OUTHIST,3411)
            3411 FORMAT('**********************LAVA FLOW PARAMETERS*********************')
            WRITE(OUTHIST,3412) NUMBER_OF_LAVA_SOURCES
            3412   FORMAT(' NUMBER_OF_LAVA_SOURCES=',I5)
            DO  I=1,NUMBER_OF_LAVA_SOURCES
                WRITE(OUTHIST,3414) I,I_LAVA(I),J_LAVA(I)
                3414   FORMAT(' I=',I5,' I_LAVA=',I5,' J_LAVA=',I5)
            ENDDO
            WRITE(OUTHIST,3430) MINIMUM_LAVA_FLOW_SLOPE
            3430   FORMAT(' MINIMUM_LAVA_FLOW_SLOPE=',G12.5)
            WRITE(OUTHIST,3440) LAVA_FLOW_THICKNESS
            3440    FORMAT(' LAVA_FLOW_THICKNESS=',G12.5)
            WRITE(OUTHIST,3445) MINIMUM_FLOW_THICKNESS
            3445      FORMAT(' MINIMUM_FLOW_THICKNESS=',G12.5)
            WRITE(OUTHIST,3455) NEW_SEGMENT_INTERVAL
            3455   FORMAT(' NEW_SEGMENT_INTERVAL=',I8)
            WRITE(OUTHIST,3456) SOURCE_CHANGE_INTERVAL
            3456   FORMAT(' SOURCE_CHANGE_INTERVAL=',I8)
            WRITE(OUTHIST,3480) ERUPTION_STOP_PROBABILITY
            3480    FORMAT(' PROGCLOG=',G12.5)
            WRITE(OUTHIST,3485) NO_FLOW_PROBABILITY
            3485    FORMAT(' NO_FLOW_PROBABILITY=',G12.5)
            WRITE(OUTHIST,3501) LAVA_GRADIENT_WEIGHT
            3501    FORMAT(' LAVA_GRADIENT_WEIGHT=',G12.5)
            WRITE(OUTHIST,3510) LAVA_DURATION_WEIGHT
            3510    FORMAT(' LAVA_DURATION_WEIGHT=',G12.5)
        ENDIF
        !       **********************************************************************
        !
        !        Write cratering parameters
        !
        !       **********************************************************************
        IF (MODEL_IMPACT_CRATERING) THEN
            WRITE(OUTHIST,1550) IMPACT_PROBABILITY
            1550 FORMAT('***********DOING IMPACT CRATERING*************',/, &
            ' PROBABILITY OF IMPACT EVENT=',G12.5)
            IF (FLUVIAL_AND_SLOPE_MODELING) THEN
                IF(IS_REGOLITH_CRATER) THEN
                    WRITE(OUTHIST,1551)
                    1551 FORMAT(' CRATERS AND EJECTA ARE MODELED AS REGOLITH')
                ELSE
                    WRITE(OUTHIST,1552)
                    1552 FORMAT(' CRATERS AND EJECTA ARE MODELED AS BEDROCK')
                ENDIF
            ENDIF
            IF ((IS_X_PERIODIC).AND.(IS_Y_PERIODIC).AND.(IFOLD > 0)) &
            THEN
                DO_EJECTA_WRAPAROUND=.TRUE.
                WRITE(OUTHIST,4958)
                4958       FORMAT(' EJECTA MODELING FOLDS OVER LATERAL BOUNDS ')
            ELSE
                DO_EJECTA_WRAPAROUND=.FALSE.
            ENDIF
            WRITE(OUTHIST,4841) LARGE_CRATER_DEPTH_SCALE,LARGE_CRATER_DEPTH_EXPONENT, &
            LARGE_CRATER_RIM_SCALE,LARGE_CRATER_RIM_EXPONENT,TRANSITION_DIAMETER
            4841   FORMAT(' LARGE_CRATER_DEPTH_SCALE=',G12.5,' LARGE_CRATER_DEPTH_EXPONENT=', &
            G12.5,/,' LARGE_CRATER_RIM_SCALE=',G12.5, &
            ' LARGE_CRATER_RIM_EXPONENT=',G12.5,' SIMPLE-COMPLEX TRANSITION DIAMETER=',&
            G12.5)
            WRITE(OUTHIST,4821)NIMPACTS,DIFFINTERVAL,NOISEINTERVAL
            4821   FORMAT(' NIMPACTS=',I8,/,&
            ' DIFFINTERVAL=',I8,' NOISEINTERVAL=',I8)
            WRITE(OUTHIST,4822)CRATER_FREQUENCY_EXPONENT,FREQUENCY_CUTOFF_SCALING, &
            SMALLEST_MODELED_CRATER,LARGEST_MODELED_CRATER
            4822   FORMAT(' CRATER_FREQUENCY_EXPONENT=',G12.5,' FREQUENCY_CUTOFF_SCALING=',G12.5,/, &
            ' SMALLEST_MODELED_CRATER=',G12.5, &
            ' LARGEST_MODELED_CRATER=',G12.5)
            WRITE(OUTHIST,4823)EJECTA_THICKNESS_VARIABILITY,NOISESD
            4823   FORMAT(' EJECTA_THICKNESS_VARIABILITY=',G12.5,' NOISESD=',G12.5)
            WRITE(OUTHIST,4824)INHERITANCE_PARAMETER
            4824   FORMAT(' INHERITANCE_PARAMETER=',G12.5)
            IF (EJECTA_THICKNESS_VARIABILITY > 0.0) THEN
                RANDOM_EJECTA_THICKNESS=.TRUE.
                !       **********************************************************************
                !
                !        ******* Renormalize ejecta_thickness_variability for a lognormal distribution of
                !                deposit depths with mean unity
                !
               !       **********************************************************************
                   EJECTA_THICKNESS_VARIABILITY=EJECTA_THICKNESS_VARIABILITY/DSQRT(DEXP(1.0D0)*(DEXP(1.0D0)-1.0D0))
            ELSE
                RANDOM_EJECTA_THICKNESS=.FALSE.
            ENDIF
            IF (NOISESD > 0.0) THEN
                MICRONOISE=.TRUE.
            ELSE
                MICRONOISE=.FALSE.
            ENDIF
            ITOTHITS=0
            ALINV=1.0/CRATER_FREQUENCY_EXPONENT
            !       **********************************************************************
            !
            !       ******* Determine the maximum range and values of xmin,xmax,ymin,ymax
            !               using the largest crater size
            !
            !       **********************************************************************
            CCON=1.0/(1.0/LARGEST_MODELED_CRATER**CRATER_FREQUENCY_EXPONENT-1.0/SMALLEST_MODELED_CRATER**CRATER_FREQUENCY_EXPONENT)
            BCON=-CCON/SMALLEST_MODELED_CRATER**CRATER_FREQUENCY_EXPONENT
            DIAMETER=LARGEST_MODELED_CRATER
            RADIUS=DIAMETER/2.0
            CRATER_DEPTH=LARGE_CRATER_DEPTH_SCALE*DIAMETER**LARGE_CRATER_DEPTH_EXPONENT
            RIM_HEIGHT=LARGE_CRATER_RIM_SCALE*DIAMETER**LARGE_CRATER_RIM_EXPONENT
            INTERIOR_SHAPE_EXPONENT=LARGE_CRATER_SHAPE_SCALE*DIAMETER**LARGE_CRATER_SHAPE_EXPONENT
            !            interior_shape_exponent=diameter/(2.0*crater_depth)
            EXTERIOR_SHAPE_EXPONENT=2.0-(RIM_HEIGHT/((RIM_HEIGHT-CRATER_DEPTH)/2.0+CRATER_DEPTH/(INTERIOR_SHAPE_EXPONENT+2)))
            MAXRANGE=DMIN1(MX*CELL_SIZE,RADIUS*(1.0/0.01)**(1.0/(EXTERIOR_SHAPE_EXPONENT-1.0)))
            !       **********************************************************************
            !
            !        ******* The total playing field size is scaled so that all significant
            !                depostion from a 50 km crater off the target field will be
            !                accounted for
            !
            !       **********************************************************************
            CALL FIND_MODIFICATION_RANGE()
        ENDIF
        CALL RANDOM_SEED(SIZE=NRANDOM)
        ALLOCATE(RANDSEED(NRANDOM))
        RANDSEED=ISEED
        CALL RANDOM_SEED(PUT=RANDSEED)
        DEALLOCATE (RANDSEED)
        !    *************************************************************************
        !         Write accretion-ablation parameters
        ! ****************************************************************************
        IF (MODEL_ACCRETION_AND_ABLATION) THEN
            WRITE(OUTHIST,1553) ACCRETION_RATE
            1553 FORMAT('****************MODELING ACCRETION AND ABLATION***********',/, &
            ' ACCRETION RATE=',G12.5)
        ENDIF
        IF (EXPOSURE_DEPENDENT_CREEP) THEN
            WRITE(OUTHIST,4987)
            4987 FORMAT ('**********MODELING EXPOSURE-DEPENDENT-CREEP**************',/, &
            ' USES EXPOSURE-DEPENDENT CREEP')
        ENDIF
        IF (USE_SOLAR_EROSION) THEN
            WRITE(OUTHIST,4986) RAD_CONST, RAD_DUST_FACTOR,RAD_THRESH_CONVEXITY,RAD_DEPOSIT_RATE
            4986 FORMAT ('**********************USES SOLAR RADIATION EROSION************',/, &
            'RADIATION_EROSION_SCALING=',G12.5, &
            ' DUST ENHANCEMENT FACTIOR=',G12.5,/,' CRITICAL CONVEXITY FOR ICE DEPOSITION=',G12.5, &
            ' ICE_DEPOSITION_RATE=',G12.5)
        ENDIF
        IF (MODEL_ACCRETION_AND_ABLATION.OR.EXPOSURE_DEPENDENT_CREEP.OR.USE_SOLAR_EROSION) THEN
            IF (USE_TOP_EXPOSURE) THEN
                WRITE(OUTHIST,4984)
                4984 FORMAT ('USES TOP EXPOSURE CALCULATION')
            ENDIF
            IF (USE_INVERSE_EXPOSURE) THEN
                WRITE(OUTHIST,4983)
                4983 FORMAT(' USES INVERSE MEASURE OF EXPOSURE')
            ENDIF
            IF (USE_EXPOSURE_SMOOTHING) THEN
                WRITE(OUTHIST,4982)
                4982 FORMAT(' USES SMOOTHING OF EXPOSURE INDEX')
            ENDIF
        ENDIF
         !       **********************************************************************
        !
        !       ***** More variable initialization
        !
        !       **********************************************************************
        NOISECALL=.TRUE.
        DIFFUSECALL=.TRUE.
        IF (.NOT.NEW_SIMULATION) THEN
            !       *****************************************************************
            !        continue.dat contains information used to restart a run, including
            !           the current iteration, the maximum desired iteration, the present
            !           time, the present time increment, and the reference elevation at the
            !           close of the last run (because relative elevations are output)
            !       *****************************************************************
            OPEN(OUTCONT,FILE=trim(INPUT_DIRECTORY_PATH) // 'CONTINUE.DAT',FORM='UNFORMATTED')
            READ(OUTCONT) ITERATION,PRESENT_TIME,EREFERENCE,EEEMIN &
            , ELEVATION,RELATIVE_RESISTANCE,SEDIMENT_BASE,INITIAL_ELEVATION,PREVIOUS_ELEVATION,REGOLITH,DEFORMATION &
            , SEDIMENT_YIELD &
            ,CFNE,CFNW,CFW,CFN,IS_SEDIMENT_COVERED,ACCELERATED_EROSION,DO_ACCELERATED_EROSION,IS_ROCK_SURFACE,SEDIMENT_DEPOSITED  &
            ,EROSION_DEPTH_INDEX,FLOW_DIRECTION,IDO,IFILE1,IFILE2  &
            ,ELAPSEDTIME,DEPOSITWORK,ERODEWORK,CRATERWORK,LAVAWORK,SLOPEWORK &
            ,CUMULATIVE_ELEVATION_CHANGE,CUMULATIVE_CRATERING_CHANGE,CUMULATIVE_EOLIAN_CHANGE,CUMULATIVE_LAVA_CHANGE
            CLOSE(OUTCONT)
        ENDIF
        !        *******************************************************************
        !         Read in the initial elevations and determine the minimum elevation
        !        *******************************************************************
        IF(NEW_SIMULATION) THEN
            EREFERENCE=-1.0E25
            DO  I = 1, MX
                DO  J = 1, MY
                        READ (INDATA,*) ELEVATION(I,J)
                        ELEVATION(I,J)=ELEVATION(I,J)*VERTICAL_SCALING
                    INITIAL_ELEVATION(I,J)=ELEVATION(I,J)
                    IF (FLUVIAL_AND_SLOPE_MODELING) THEN
                        SEDIMENT_BASE(I,J)=ELEVATION(I,J)
                    ENDIF
                    IF (FLUVIAL_AND_SLOPE_MODELING.AND.RESISTANT_SURFACE_LAYER) THEN
                        PREVIOUS_ELEVATION(I,J)=ELEVATION(I,J)-SURFACE_LAYER_THICKNESS
                    ENDIF
                    IF (ELEVATION(I,J) > EREFERENCE) EREFERENCE=ELEVATION(I,J)
                    IF (ELEVATION(I,J) < THEMINIMUM) THEMINIMUM=ELEVATION(I,J)
                    IF (FLUVIAL_AND_SLOPE_MODELING) THEN
                        EROSION_DEPTH_INDEX(I,J)=-1
                    ENDIF
                ENDDO
            ENDDO
            IF (READREGOLITH > 0) THEN
                OPEN(IOTEMP3,FILE=trim(INPUT_DIRECTORY_PATH) // 'INREG.DAT')
                CALL READ_REGOLITH_THICKNESS(IOTEMP3)
                CLOSE(IOTEMP3)
            ENDIF
        ENDIF
        WRITE(OUTHIST,221)
        221   FORMAT(' INITIAL ELEVATION')
        CALL SUMMARIZE_MATRIX_DATA(ELEVATION,THEAVG,THEMAX,THEMIN)
        IF (IREADALLUV == 1) THEN
            WRITE(OUTHIST,222)
            222       FORMAT(' INITIAL ALLUVIAL BASE')
            CALL SUMMARIZE_MATRIX_DATA(SEDIMENT_BASE,THEAVG,THEMAX,THEMIN)
            WRITE(OUTHIST,223)
            223       FORMAT(' INITIAL SEDIMENT COVER')
            CALL SUMMARIZE_LOGICAL_MATRIX(IS_SEDIMENT_COVERED)
        ENDIF
        !       ********************************************************************
        !        If we use temporally variable erosion rates, read in a table.
        !           The first line (number_of_parameter_changes) tells how many erosion rate intervals
        !           there are to be read.  For each rate interval the closing time for
        !           that erosion rate catagory and the erosion rate during that interval
        !           are read in, as well as the potentially time-varying simulation
        !           parameters (see the subroutine doeroderate for definition of the
        !           individual simulation_parameters values)
        !       ********************************************************************
        IF (VARIABLE_EROSION_RATE) THEN
            READ(INRATES,*) NUMBER_OF_PARAMETER_CHANGES
            WRITE(OUTHIST,738) NUMBER_OF_PARAMETER_CHANGES
            738   FORMAT(' NUMBER OF EROSION RATE CATEGORIES=',I5,/, &
            'TABLE OF LOWEST ELEVATIONS AND THEIR EROSION RATES',/)
            DO  I=1,NUMBER_OF_PARAMETER_CHANGES
                READ(INRATES,*) TIMES_FOR_PARAMETER_CHANGES(I),EROSION_RATE_VALUES(I)
                READ(INRATES,*)(SIMULATION_PARAMETERS(I,J),J=1,12)
                WRITE(OUTHIST,737) TIMES_FOR_PARAMETER_CHANGES(I),EROSION_RATE_VALUES(I)
                WRITE(OUTHIST,740) (SIMULATION_PARAMETERS(I,J),J=1,12)
                737   FORMAT (2G12.5)
                740   FORMAT(5G12.5,/,5G12.5,/,2G12.5)
            ENDDO
        ENDIF
        IF (.NOT.IS_Y_PERIODIC) THEN
            DO  I=1,MX
                IW=I-1
                IE=I+1
                IF (IS_X_PERIODIC) THEN
                    IF (IE > MX) IE=1
                    IF (IW < 1) IW=MX
                ELSE
                    IF (IE > MX) IE=MX
                    IF (IW < 1) IW=1
                ENDIF
                IF (FLUVIAL_AND_SLOPE_MODELING) THEN
                    IF (HORIZONTAL_LOWER_BOUNDARY) THEN
                        ERODING_LOWER_BOUNDARY(I)=.FALSE.
                    ELSE
                        IF (NON_ERODING_LOWER_BOUNDARY) THEN
                            ERODING_LOWER_BOUNDARY(I)=.FALSE.
                        ELSE
                            IF ((ELEVATION(I,MY) <= ELEVATION(IW,MY)).AND.(ELEVATION(I,MY) <= ELEVATION(IE,MY))) &
                            THEN
                                ERODING_LOWER_BOUNDARY(I)=.FALSE.
                            ELSE
                                ERODING_LOWER_BOUNDARY(I)=.TRUE.
                            ENDIF
                        ENDIF
                    ENDIF
                ENDIF
            ENDDO
        ENDIF
        RETURN
    END
    !        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE REPORT_MAX()
        USE ERODE_GLOBALS
        INTEGER SUBDET,SEDDET
        REAL (8) :: DIF1
        IF (SUBMERGED(IDET,JDET)) THEN
            SUBDET=1
        ELSE
            SUBDET=0
        ENDIF
        IF (IS_SEDIMENT_COVERED(IDET,JDET)) THEN
            SEDDET=1
        ELSE
            SEDDET=0
        ENDIF
        DIF1=ELEVATION(IDET,JDET)-PREVIOUS_ELEVATION(IDET,JDET)
        WRITE(OUTHIST,3992) IDET,JDET,SUBDET,SEDDET, FLOW_DIRECTION(IDET,JDET) &
        ,IDO(IDET,JDET),DISCHARGE(IDET,JDET),D8_GRADIENT(IDET,JDET),MAXIMUM_DISCHARGE(IDET,JDET) &
        ,MAXIMUM_ELEVATION_CHANGE,ELEVATION(IDET,JDET),ERODE_CHANNEL(IDET,JDET),ERODE_SLOPE(IDET,JDET) &
        ,CFN(IDET,JDET),CFNE(IDET,JDET),CFW(IDET,JDET),CFNW(IDET,JDET) &
        ,SEDIMENT_YIELD(IDET,JDET),MAXIMUM_SEDIMENT_YIELD(IDET,JDET) &
        ,PREVIOUS_ELEVATION(IDET,JDET),ERODE_REGOLITH_CHANNEL(IDET,JDET),SEDIMENT_FLUX(IDET,JDET) &
        ,EQUILIBRIUM_GRADIENT(IDET,JDET),DIF1
        3992  FORMAT('I=',I4,' J=',I4,' SUB=',I1,' SED=',I1,' ID=',I2 &
        ,' IDO=',I2,' QV=',G12.5,' S=',G12.5,'QM=',G12.5,/, &
        'MC=',G12.5,' ELEVATION=',G12.5,' ECH=',G12.5,' ESL=',G12.5, &
        ' CFN=',G12.5,' CFNE=',G12.5,' CFW=',G12.5,' CFNW=',G12.5 &
        ,/,'SY=',G12.5,' SYM=',G12.5,' PREVIOUS_ELEVATION=',G12.5,' ERG=',G12.5, &
        ' SQ=',G12.5,' SE=',G12.5,' EDIFF=',G12.5)
        RETURN
    END

    MODULE MARSSIM_IRF
        USE ERODE_GLOBALS
        USE LAKE_GLOBALS
        USE CRATER_GLOBALS
        USE EOLIAN_GLOBALS
        USE LAVA_GLOBALS
        USE ACCRETION_GLOBALS

        IMPLICIT NONE
        INTEGER I,J,K,KLOW
	    REAL (8) :: EEEMIN
	    LOGICAL :: DODETAILS

        SAVE
        CONTAINS

        SUBROUTINE MARSSIM_INITIALIZE(ARG_INPUT_DIRECTORY_PATH, INPUT_FILE)
	    IMPLICIT NONE
	    CHARACTER *(*):: ARG_INPUT_DIRECTORY_PATH
	    CHARACTER *(*):: INPUT_FILE
   
	    INPUT_DIRECTORY_PATH=ARG_INPUT_DIRECTORY_PATH
        !       *********************************************************************
        !        This is the master program for the marssim landform evolution model.
        !        The basic algorithms, assumptions and structure of the portions of the
        !        model dealing with weathering, slope processes, fluvial erosion, fluvial transport
        !        and sedimentation have been sumamrized in several publications:
        !          howard, a. d., 1994, water res. research, 30(7), 2261-85.
        !          howard, a. d., 1997, earth surf. proc. landforms, 22, 211-227.
        !          howard, a. d., 1999, in incised river channels, s. darby and a. simon, eds.,
        !                   wiley, 277-300.
        !          fagherazzi, s., howard, a. d., and wiberg, p. l.,2004, j. geophys res. 109,
        !              doi:10.1029/2003jf000091.
        !          howard, a. d., 2007, geomorphology, 91, 332-363.
        !          barnhart, c. j., howard, a. d., and moore, j. m., 2008, j. geophys. res., 114,
        !            e01003, doi:10.1029/2008je003122
        !
        !        References to publications related to the more specialized planetary processes are
        !        given in individual subroutines.
        !
        !        Portability issues:  The program is almost entirely written in standard fortran 90.
        !          Most data file output and input is from ascii files and should be identical on all architectures.
        !          However, the program
        !            also writes out binary image files (Photoshop 'raw' files), mostly shaded relief images.
        !          These routines, "imagewrite", "shadewrite" , "write_laava+ages",
        !            and "colorwrite" might have to be modified for other compilers or architectures.
        !          These use the non-standard
        !            output control character '$', which does not terminate a record.
        !          Why is FORTRAN so backward in terms of output format control???
        !          For standard compilers such as gfortran the next best solution is to change the format statement from:
        !                 130 format(a1,$)
        !          to:
        !                 130 format(a1)
        !          and the corresponding write statements from:
        !                 write(5,130) achar
        !          to:
        !                 write(5,130,advance='no')
        !          The disadvantage of this is that a record terminator (CR & LF) or (LF) is appended to the end
        !            of the file and you need to deal with this when importing into, say,
        !            Photoshop by telling it to ignore the last 1 or 2 bytes.
        !          All image-writing routines that use the $ formatting. are in "read_and_write_data_files.f90".
        !          A replacement routine using the above approach with advance='no' is provided as "alternate_read_and_write.f90".
        !          Probably the best solution is to link in a standard image library, such as TIFF and to rewrite
        !            the routines to directly write formatted image files.
        !          But I'm happy with my PGI and INTEL compilers and raw images.
        !        ********************************************************************
        !        define input - output file unit numbers
        !   MODIFIES:  (MOST INPUT-OUTPUT FILE NUMBERS), SUBMERGED, ITERATION, PREVIOUS_ELEVATION, ELAPSEDTIME,
        !              WRITE_OUT_CHANNELS, PRESENT_TIME, EVENT_DONE, EVENT_INDEX, DO_EVENTS
        !              SUBGRADCHANGE, ABSGRADCHANGE, NUMIDCHANGE
        !   CALLS:  READ_INPUT_PARAMETERS, SETUP_FLUVIAL_SLOPE_EROSION, WRITE_SHADED_RELIEF_IMAGE,
        !           WRITE_COLOR_SHADED_RELIEF_IMAGE, WRITE_IMAGE, DO_FLUVIAL_AND_SLOPE, DEPOSIT_ICE,
        !           REPORT_MAX, FINDOCEAN_ELEVATION, WRITE_DEBUG_INFO, CHANNEL_PROPERTIES,WRITE_FIRST_DATA_MATRIX
        !           PRINT_MORPHOMETRY, CALCULATE_TOPO_DIVERGENCE, SUMMARIZE_CHANNELS, MAKE_EVENT
        !           WRITE_SECOND_DATA_MATRIX, WRITE_SOURDED_DISCHARGE, WRITE_SUBMERGED_LOCATIONS, WRITE_REPORT
        !           PRINT_SIMULATION_INFORMATION, DO_LAVA_FLOWS, WRITE_LAVA_INFO, WRITE_LAVA_AGES, DO_EOLOIAN_CHANGE
        !           DO_IMPACT_CRATERING, DO_EXPOSURE_DEPENDENT_CREEP, DO_ACCRETION_ABLATION, FINALIZE_FLUVIAL_SLOPE_EROSION
        !           WRITE_ALLUVIAL_LOCATIONS, WRITE_BEDROCK_LOCATIONS, WRITE_EROSION_DEPTH_INDEX, WRITE_DEFORMATION
        !           WRITE_SEDIMENT_BASE, WRITE_LAKE_INFO, WRITE_REGOLITH_THICKNESS, WRITE_ROCK_RESISTANCE
        !           PRINT_SIMULATION_INFORMATION, WRITE_ACCELERATED_EROSION_STATE, WRITE_GROUNDWATER_ELEVATION
        !           WRITE_GROUNDWATER_FLOW
        !
        !        ********************************************************************

        PARAMETER_CHANGE_INDEX=1

        INDATA=15
        PARAMS=16
        INRESIST=13
        OUTDATA=19
        OUTHIST=17
        OUTCHAN=18
        OUTRECORD=11
        OUTSUMMARY=12
        OUTSOURCE=14
        OUTREPORT=30
        OUTCRATER=73
        OUTSEDDEBUG=74

        OUTLAKE=75
        INRATES=32
        OUTIMGDAT=34
        OUTROCK=35
        OUTRESIST=36
        OUTGRAD=26
        OUTALLUV=25
        OUTEROSION_DEPTH_INDEX=45
        OUTCONT=46
        OUTELEV=47
        IOTEMP1=21
        IOTEMP2=22
        OUTBASE=23
        OUTDEFORM=48
        OUTHIGH=49
        OUTREGOLITH=51
        EEEMIN=0.0
        DODETAILS=.FALSE.

        !        ********************************************************************
        !        MARSSIM.PRM is file of simulation parameters
        !        ********************************************************************

        !OPEN(PARAMS,FILE='MARSSIM.PRM')
        OPEN(PARAMS,FILE=INPUT_FILE)
        !      ```````````
        OPEN (82,FILE='DEBUG.PRN')
        !      ``````````
        !        ********************************************************************
        !        resist.in is file of random or systematic erodibility
        !        values for a 3-d cube of dimensions mx,my,mz.  It is a direct access
        !           file read by the subroutine readerodiblity in weather.f
        !           the erodibilty matrix must be created by a separate program.
        !           for example, matrix_3d.f creates a fractal cube and readmat.f
        !           converts it into a fortran direct access file format. The i dimesion
        !           is horizontal (left to right), the j dimension is horizontal (top
        !           to bottom, and the k dimension is vertical (increasing downward)
        !        ********************************************************************
        OPEN(INRESIST,FILE=trim(INPUT_DIRECTORY_PATH) // 'RESIST.IN',ACCESS='DIRECT', &
        FORM='UNFORMATTED',RECL=4)
        !        ********************************************************************
        !        basin.lst is descriptive file of simulation parameters and
        !        basin statistics during simulation
        !        ********************************************************************
        OPEN(OUTHIST,FILE='BASIN.LST')
        !        ********************************************************************
        !        channel.dat is list of drainage network loations and flow directions
        !        ********************************************************************
        OPEN(OUTCHAN,FILE='CHANNEL.DAT')
        !        ********************************************************************
        !        record.dat is historical record of progression of simulation
        !        towards steady state
        !        ********************************************************************
        OPEN(OUTRECORD,FILE='RECORD.DAT')
        !        ********************************************************************
        !        report.dat is historical record of relief & erosion rate
        !        ********************************************************************
        OPEN(OUTREPORT,FILE='REPORT.PRN')
        !        ********************************************************************
        !        summary.dat is like basin.lst without text
        !        ********************************************************************
        OPEN(OUTSUMMARY,FILE='SUMMARY.DAT')
        !        ********************************************************************
        !        source.dat is output file of channel source locations
        !        ********************************************************************
        OPEN(OUTSOURCE,FILE='SOURCE.DAT')
        !       **********************************************************************
        !        inrates.dat is file of time-dependent erosion rates and simulation
        !           parameters
        !       **********************************************************************
        OPEN(INRATES,FILE=trim(INPUT_DIRECTORY_PATH) // 'INRATES.DAT')
        OPEN(OUTCRATER,FILE='CRATER.DAT')
        OPEN(78,FILE='STATISTICS.PRN')
        !        ********************************************************************
        !        Read in initial parameters & elevations and summarize
        !        ********************************************************************
        CALL READ_INPUT_PARAMETERS(trim(INPUT_DIRECTORY_PATH))
        CALL SETUP_FLUVIAL_SLOPE_EROSION()
        SUBMERGED=.FALSE.
        !       **********************************************************************
        !        Write the raw image file for initial conditions
        !       **********************************************************************
        IF (NEW_SIMULATION) THEN
            CALL WRITE_IMAGE()
            IF (USE_SOLAR_EROSION) THEN
                CALL WRITE_COLOR_SHADED_RELIEF_IMAGE()
                CALL WRITE_SHADED_RELIEF_IMAGE()
            ELSE
                CALL WRITE_SHADED_RELIEF_IMAGE()
            ENDIF
        ENDIF
        ELOWEST=-1.0E+25
        IF ( NEW_SIMULATION ) THEN
            KLOW=1
        ELSE
            KLOW=ITERATION-1
        ENDIF
   END SUBROUTINE

   SUBROUTINE MARSSIM_RUN()
        IMPLICIT NONE
        INTEGER NCALL,SOURCEUSE
        REAL (8) :: EXPONENTIAL_DISTRIBUTION,NETDIFF
        REAL (8) :: LAVAPROBTOUSE,EOPROBTOUSE,CRATERPROBTOUSE
        REAL (4) :: PROBVALUE
        REAL (8) :: RRAND
        EXTERNAL RRAND

        L200: DO  K=KLOW,MAXIMUM_ITERATION
            ITERATION=ITERATION+1
            TOTAL_ITERATIONS=TOTAL_ITERATIONS+1
            CALL GRADIENT_AND_FLOW_DIRECTION()
            IF (FLUVIAL_AND_SLOPE_MODELING) THEN
                IF (BISTABLE_FLUVIAL_EROSION) THEN
                ELSE
                    IF(RESISTANT_SURFACE_LAYER) THEN
                    ELSE
                        PREVIOUS_ELEVATION=ELEVATION
                    ENDIF
                ENDIF
                CALL DO_FLUVIAL_AND_SLOPE()
                IF (USE_SOLAR_EROSION) CALL DEPOSIT_ICE()
                IF (DODETAILS) THEN
                    CALL REPORT_MAX()
                ENDIF
                CALL FINDOCEAN_ELEVATION(PRESENT_TIME)
                IF (PRESENT_TIME >= MAXIMUM_SIMULATION_TIME) EXIT L200
                IF (DO_DEBUGGING) THEN
                    IF ((MAXIMUM_ITERATION-ITERATION) < 4) CALL WRITE_DEBUG_INFO()
                ENDIF
            ENDIF
            IF (FLUVIAL_AND_SLOPE_MODELING) THEN
                ELAPSEDTIME=ELAPSEDTIME+TIME_INCREMENT
            ELSE
                ELAPSEDTIME=ELAPSEDTIME+1.0
            ENDIF
            IF (MOD(ITERATION,ELEVATION_PRINT_INTERVAL) == 0) THEN
                IF (FLUVIAL_AND_SLOPE_MODELING) CALL CHANNEL_PROPERTIES()
                CALL WRITE_FIRST_DATA_MATRIX()
            ENDIF
            IF (MOD(ITERATION,ELEVATION_PRINT_INTERVAL) == 0) THEN
                IF (DO_MORPHOMETRY) THEN
                    CALL PRINT_MORPHOMETRY()
                    CALL CALCULATE_TOPO_DIVERGENCE()
                    WRITE_OUT_CHANNELS=.TRUE.
                    CALL SUMMARIZE_CHANNELS()
                    CALL WRITE_SECOND_DATA_MATRIX()
                ENDIF
            ENDIF
            !       ********************************************************************
            !        If it is time to write out a raw image of elevations, do so
            !       ***********************************************************************
            IF (MOD(ITERATION,IMAGE_OUTPUT_INTERVAL) == 0) THEN
                CALL WRITE_IMAGE()
                IF (USE_SOLAR_EROSION) THEN
                    CALL WRITE_COLOR_SHADED_RELIEF_IMAGE()
                    CALL WRITE_SHADED_RELIEF_IMAGE()
                ELSE
                    CALL WRITE_SHADED_RELIEF_IMAGE()
                ENDIF
                IF (FLUVIAL_AND_SLOPE_MODELING) CALL WRITE_ROUTED_DISCHARGE()
                IF (FLUVIAL_AND_SLOPE_MODELING) CALL WRITE_SUBMERGED_LOCATIONS()
            ENDIF
            !        ********************************************************************
            !        If it is time to print erosion rate statistics do so
            !        ********************************************************************
            IF (FLUVIAL_AND_SLOPE_MODELING) THEN
                IF(MOD(ITERATION,WRITE_CHANGE_INTERVAL) == 0) THEN
                    CALL WRITE_REPORT()
                ENDIF
                IF (MOD(ITERATION,OUTPUT_PRINT_INTERVAL) == 0) THEN
                    IF (MOD(ITERATION,ELEVATION_PRINT_INTERVAL) /= 0) THEN
                        CALL PRINT_SIMULATION_INFORMATION()
                    ENDIF
                ENDIF
                IF (MOD(ITERATION,ELEVATION_PRINT_INTERVAL) == 0) THEN
                    CALL WRITE_ALLUVIAL_LOCATIONS(OUTALLUV)
                    CALL WRITE_BEDROCK_LOCATIONS(OUTROCK)
                    CALL WRITE_EROSION_DEPTH_INDEX(OUTEROSION_DEPTH_INDEX)
                    CALL WRITE_DEFORMATION(OUTDEFORM)
                    CALL WRITE_SEDIMENT_BASE()
                    CALL WRITE_LAKE_INFO(OUTLAKE)
                    CALL WRITE_REGOLITH_THICKNESS()
                    CALL WRITE_ROCK_RESISTANCE()
                    CALL PRINT_SIMULATION_INFORMATION()
                    CALL WRITE_ACCELERATED_EROSION_STATE(OUTHIGH)
                    IF (MODEL_GROUNDWATER) THEN
                        CALL WRITE_GROUNDWATER_ELEVATION(IOTEMP1)
                        CALL WRITE_GROUNDWATER_FLOW(IOTEMP1)
                    ENDIF
                    OPEN(OUTCONT,FILE=trim(INPUT_DIRECTORY_PATH) // 'CONTINUE.DAT',FORM='UNFORMATTED')
                    WRITE(OUTCONT) ITERATION,PRESENT_TIME,EREFERENCE,EEEMIN  &
                    , ELEVATION,RELATIVE_RESISTANCE,SEDIMENT_BASE,INITIAL_ELEVATION,PREVIOUS_ELEVATION,REGOLITH,DEFORMATION &
                    , SEDIMENT_YIELD &
                    ,CFNE,CFNW,CFW,CFN,IS_SEDIMENT_COVERED,ACCELERATED_EROSION,DO_ACCELERATED_EROSION &
                    ,IS_ROCK_SURFACE,SEDIMENT_DEPOSITED &
                    ,EROSION_DEPTH_INDEX,FLOW_DIRECTION,IDO,IFILE1,IFILE2 &
                    ,ELAPSEDTIME,DEPOSITWORK,ERODEWORK,CRATERWORK,LAVAWORK,SLOPEWORK &
                    ,CUMULATIVE_ELEVATION_CHANGE,CUMULATIVE_CRATERING_CHANGE,CUMULATIVE_EOLIAN_CHANGE,CUMULATIVE_LAVA_CHANGE
                    CLOSE(OUTCONT)
                ENDIF
            ENDIF
            IF (MODEL_LAVA_FLOWS) THEN
                IF (FLUVIAL_AND_SLOPE_MODELING) THEN
                    LAVAPROBTOUSE=LAVA_EVENT_PROBABILITY*TIME_INCREMENT
                ELSE
                    LAVAPROBTOUSE=LAVA_EVENT_PROBABILITY
                ENDIF
                IF (LAVAPROBTOUSE < 1.0) THEN
                    PROBVALUE=RRAND()
                    IF (PROBVALUE <= LAVAPROBTOUSE) THEN
                        PROBVALUE=RRAND()
                        SOURCEUSE=INT(PROBVALUE*NUMBER_OF_LAVA_SOURCES+1.0)
                        IF (SOURCEUSE < 1) SOURCEUSE=1
                        IF (SOURCEUSE > NUMBER_OF_LAVA_SOURCES) SOURCEUSE=NUMBER_OF_LAVA_SOURCES
                        WRITE(*,7011) SOURCEUSE
                        7011           FORMAT(' LAVAFLOW FROM SOURCE ',I3)
                        CALL DO_LAVA_FLOWS(SOURCEUSE)
                        IF (MOD(ITERATION,ELEVATION_PRINT_INTERVAL) == 0) THEN
                            CALL WRITE_LAVA_INFO()
                            CALL WRITE_LAVA_AGES()
                        ENDIF

                    ENDIF
                ELSE
                    NCALL=INT(EXPONENTIAL_DISTRIBUTION()*LAVAPROBTOUSE+0.5)
                    IF (NCALL >= 1) THEN
                        DO  I=1,NCALL
                            PROBVALUE=RRAND()
                            SOURCEUSE=INT(PROBVALUE*NUMBER_OF_LAVA_SOURCES+1.0)
                            IF (SOURCEUSE < 1) SOURCEUSE=1
                            IF (SOURCEUSE > NUMBER_OF_LAVA_SOURCES) SOURCEUSE=NUMBER_OF_LAVA_SOURCES
                            WRITE(*,7011)
                            CALL DO_LAVA_FLOWS(SOURCEUSE)
                        ENDDO
                    ENDIF
                ENDIF
                IF (MOD(ITERATION,ELEVATION_PRINT_INTERVAL) == 0) THEN
                    CALL WRITE_LAVA_INFO()
                    CALL WRITE_LAVA_AGES()
                ENDIF
            ENDIF
            IF (MODEL_EOLIAN_CHANGES) THEN
                IF (FLUVIAL_AND_SLOPE_MODELING) THEN
                    EOPROBTOUSE=EOLIAN_EVENT_PROBABILITY*TIME_INCREMENT
                ELSE
                    EOPROBTOUSE=EOLIAN_EVENT_PROBABILITY
                ENDIF
                IF (EOPROBTOUSE < 1.0) THEN
                    PROBVALUE=RRAND()
                    IF (PROBVALUE <= EOPROBTOUSE) THEN
                        WRITE(*,7012)
                        7012                FORMAT(' EOLIAN')
                        CALL DO_EOLIAN_CHANGE()
                        IF (.NOT.FLUVIAL_AND_SLOPE_MODELING) THEN
                            PRESENT_TIME=PRESENT_TIME+EOLIAN_TIME_INCREMENT
                        ENDIF
                    ENDIF
                ELSE
                    NCALL=INT(EXPONENTIAL_DISTRIBUTION()*EOPROBTOUSE+0.5)
                    IF (NCALL >= 1) THEN
                        DO  I=1,NCALL
                            WRITE(*,7012)
                            CALL DO_EOLIAN_CHANGE()
                            IF (.NOT.FLUVIAL_AND_SLOPE_MODELING) THEN
                                PRESENT_TIME=PRESENT_TIME+EOLIAN_TIME_INCREMENT
                            ENDIF
                        ENDDO
                    ENDIF
                ENDIF
            ENDIF
            IF (MODEL_IMPACT_CRATERING) THEN
                IF (FLUVIAL_AND_SLOPE_MODELING) THEN
                    CRATERPROBTOUSE=IMPACT_PROBABILITY*TIME_INCREMENT
                ELSE
                    CRATERPROBTOUSE=IMPACT_PROBABILITY
                ENDIF
                IF (CRATERPROBTOUSE < 1.0) THEN
                    PROBVALUE=RRAND()
                    IF (PROBVALUE <= CRATERPROBTOUSE) THEN
                        CALL DO_IMPACT_CRATERING()
                    ENDIF
                ELSE
                    NCALL=INT(EXPONENTIAL_DISTRIBUTION()*CRATERPROBTOUSE+0.5)
                    IF (NCALL >= 1) THEN
                        DO  I=1,NCALL
                            CALL DO_IMPACT_CRATERING()
                        ENDDO
                    ENDIF
                ENDIF
            ENDIF
            IF (MODEL_ACCRETION_AND_ABLATION) THEN
                WRITE(*,3221)
                3221     FORMAT(' CALL ACCRETION')
                IF (EXPOSURE_DEPENDENT_CREEP) THEN
                    CALL DO_EXPOSURE_DEPENDENT_CREEP()
                ELSE
                    CALL DO_ACCRETION_ABLATION()
                ENDIF
            ENDIF
            NETDIFF=0.0
            DO J=1,MY
                DO I=1,MX
                    NETDIFF=NETDIFF+ELEVATION(I,J)-INITIAL_ELEVATION(I,J)
                ENDDO
            ENDDO
            NETDIFF=NETDIFF/(MX*MY)
            WRITE(OUTHIST,798) ITERATION,PRESENT_TIME,MAXIMUM_ELEVATION_CHANGE,OCEAN_ELEVATION &
            ,NETDIFF
            WRITE(*,798) ITERATION,PRESENT_TIME,MAXIMUM_ELEVATION_CHANGE,OCEAN_ELEVATION  &
            ,NETDIFF
            798           FORMAT(' I=',I5,' T=',G12.5,' MC=',G12.5  &
            ,' OCEAN=',G12.5,' NETDIF=',F12.5)

            IF (DO_EVENTS) THEN
                IF (EVENT_TYPE == 1) THEN
                    IF(K >= EVENT_ITERATIONS(EVENT_INDEX)) THEN
                        IF (.NOT.EVENT_DONE(EVENT_INDEX)) THEN
                            CALL MAKE_EVENT(EVENT_INDEX)
                            EVENT_DONE(EVENT_INDEX)=.TRUE.
                            EVENT_INDEX=EVENT_INDEX+1
                            IF ((EVENT_INDEX > 10).OR. &
                            (EVENT_INDEX > NUMBER_OF_EVENTS))  &
                            DO_EVENTS=.FALSE.
                        ENDIF
                    ENDIF
                ELSE
                    IF(PRESENT_TIME >= EVENT_TIMES(EVENT_INDEX)) THEN
                        IF (.NOT.EVENT_DONE(EVENT_INDEX)) THEN
                            CALL MAKE_EVENT(EVENT_INDEX)
                            EVENT_DONE(EVENT_INDEX)=.TRUE.
                            EVENT_INDEX=EVENT_INDEX+1
                            IF ((EVENT_INDEX > 10).OR.(EVENT_INDEX > NUMBER_OF_EVENTS)) &
                            DO_EVENTS=.FALSE.
                        ENDIF
                    ENDIF
                ENDIF
            ENDIF
            IF (PRESENT_TIME >= MAXIMUM_SIMULATION_TIME) EXIT L200
        ENDDO L200
   END SUBROUTINE

   SUBROUTINE MARSSIM_FINALIZE()
        IMPLICIT NONE
   	REAL (8) :: TOTQ,BOUNDQ
        IF (FLUVIAL_AND_SLOPE_MODELING) THEN
            IF (HORIZONTAL_LOWER_BOUNDARY.AND.(.NOT.IS_Y_PERIODIC)) THEN
                BOUNDQ=0.0
                TOTQ=0.0
                DO I=1,MX
                    BOUNDQ=BOUNDQ+DISCHARGE(I,MY-1)
                    DO J=1,MY-1
                        TOTQ=TOTQ+CONVERT_DISCHARGE*D2X
                    ENDDO
                ENDDO
                WRITE(OUTHIST,4398) TOTQ,BOUNDQ,CONVERT_DISCHARGE,D2X
                WRITE(*,4398) TOTQ,BOUNDQ,CONVERT_DISCHARGE,D2X
                4398  FORMAT('TOTQ=',G12.5,' BOUNDQ=',G12.5, 'CONVERT_DISCHARGE=',G12.5,'D2X=',G12.5)
            ENDIF
            IF(MOD(ITERATION,WRITE_CHANGE_INTERVAL) /= 0) THEN
                IF (CHANGECOUNT > 0.0) THEN
                    SUMGRADCHANGE=SUMGRADCHANGE/  &
                    (CHANGECOUNT*(MX-1)*(MY-1))
                    ABSGRADCHANGE=ABSGRADCHANGE/  &
                    (CHANGECOUNT*(MX-1)*(MY-1))
                    NUMBIDCHANGE=NUMBIDCHANGE/(CHANGECOUNT* &
                    (MX-1)*(MY-1))
                    WRITE(OUTRECORD,580) PRESENT_TIME,SUMGRADCHANGE, &
                    ABSGRADCHANGE, &
                    NUMBIDCHANGE
                ENDIF
            ENDIF
            580             FORMAT(4(' ',E15.7))
            IF(MOD(ITERATION,WRITE_CHANGE_INTERVAL) /= 0) THEN
                CALL WRITE_REPORT()
            ENDIF
            IF (MOD(ITERATION,OUTPUT_PRINT_INTERVAL) /= 0) THEN
                IF (MOD(ITERATION,ELEVATION_PRINT_INTERVAL) /= 0) THEN
                    CALL CHANNEL_PROPERTIES()
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
                ENDIF
                IF (MOD(ITERATION,IMAGE_OUTPUT_INTERVAL) /= 0) THEN
                    CALL WRITE_IMAGE()
                    IF (USE_SOLAR_EROSION) THEN
                        CALL WRITE_COLOR_SHADED_RELIEF_IMAGE()
                        CALL WRITE_SHADED_RELIEF_IMAGE()
                    ELSE
                        CALL WRITE_SHADED_RELIEF_IMAGE()
                    ENDIF
                    IF (FLUVIAL_AND_SLOPE_MODELING) CALL WRITE_ROUTED_DISCHARGE()
                    IF (FLUVIAL_AND_SLOPE_MODELING) CALL WRITE_SUBMERGED_LOCATIONS()
                ENDIF
                IF (MODEL_LAVA_FLOWS) THEN
                    IF (MOD(ITERATION,ELEVATION_PRINT_INTERVAL) /= 0) THEN
                       CALL WRITE_LAVA_INFO()
                       CALL WRITE_LAVA_AGES()
                    ENDIF
                ENDIF
            ENDIF
        ENDIF
        CALL FINALIZE_FLUVIAL_SLOPE_EROSION()
        CLOSE(OUTHIST)
        STOP
       END SUBROUTINE
   END MODULE MARSSIM_IRF



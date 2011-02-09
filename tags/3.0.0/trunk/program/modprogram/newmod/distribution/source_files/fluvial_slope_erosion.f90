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
    SUBROUTINE FLUVIAL_DETACHMENT(MAXCHAN,MAXREG,NOMINAL_CRITICAL_SHEAR,LOCAL_BEDROCK_ERODIBILITY)
        USE ERODE_GLOBALS, FLUVIAL_TIMESTEP_FACTOR=>CFN
        IMPLICIT NONE
        INTEGER(4) :: I,J,NNN
        REAL(8) :: AFACTOR,GFACTOR,TEMP,XAREA,A1TERM,A2TERM,QSSS,TRANSPORT_STAGE,SEDSATFACTOR
        REAL(8) :: SATRATIO,TEMP1,TEMPKW, LOCAL_CRITICAL_SHEAR
        REAL(8) :: LOCAL_REGOLITH_ERODIBILITY,LOCAL_REGOLITH_CRITICAL_SHEAR
        REAL(8) :: LOCAL_INVERSE_RESISTANCE,REGOLITH_CRITICAL_SHEAR
        REAL(8) :: SHEAR_STRESS_FACTOR
        REAL(8) :: BEDLOAD_FLUX
        REAL(8) :: NORMVAR,LOGSD,NORMMEAN,NORMSD,RRRR
        REAL(8) :: LOGNORMAL_RANDOM_DEVIATE1
        EXTERNAL :: LOGNORMAL_RANDOM_DEVIATE1
        REAL(8), INTENT (OUT) :: MAXCHAN, MAXREG, NOMINAL_CRITICAL_SHEAR, LOCAL_BEDROCK_ERODIBILITY
        LOGICAL :: WHATSUP,DODEBUGWINDOW
        !     **********************************************************************
        !      The following lines of code generate values of the critical shear
        !         stress for each iteration if a random threshold is used.  The values
        !         are autoregressive if omega_weight is greater than zero.  Similarly
        !         if a random discharge is used, the discharge-dependent variable
        !         discharge_coefficient is calculated with autoregressive weighting
        !  MODIFIES: PEREVIOUS_CRITICAL_SHEAR, CFN, ERODE_CHANNEL, ERODE_REGOLITH_CHANNEL
        !  CALLS: LOGNORMAL_RANDOM_DEVIATE1, LOCAL_VALUES, SEDIMENT_TRANSPORT_FLUX
        !     **********************************************************************
        DODEBUGWINDOW=.FALSE.
        TEMPKW=DISCHARGE_COEFFICIENT
        NOMINAL_CRITICAL_SHEAR=DETACHMENT_CRITICAL_SHEAR
        LOCAL_BEDROCK_ERODIBILITY=BEDROCK_ERODIBILITY
        IF (RANDOM_CRITICAL_SHEAR) THEN
            LOGSD=DETACHMENT_CRITICAL_SHEAR*CRITICAL_SHEAR_VARIABILITY*OCORRECT
            NORMVAR=DLOG((LOGSD**2+DETACHMENT_CRITICAL_SHEAR**2)/DETACHMENT_CRITICAL_SHEAR**2)
            NORMMEAN=DLOG(DETACHMENT_CRITICAL_SHEAR**2/(DETACHMENT_CRITICAL_SHEAR**2+LOGSD**2))
            NORMSD=DSQRT(NORMVAR)
            RRRR=LOGNORMAL_RANDOM_DEVIATE1(NORMMEAN,NORMSD)
            NOMINAL_CRITICAL_SHEAR=PREVIOUS_CRITICAL_SHEAR*(1.0-OMEGA_WEIGHT)+OMEGA_WEIGHT/RRRR
            PREVIOUS_CRITICAL_SHEAR=NOMINAL_CRITICAL_SHEAR
        ELSE
        ENDIF
        !      ********************************************************************
        !       This is the loop for determining the rate of erosion in bedrock or
        !         regolith-flowred channels - The FLUVIAL_TIMESTEP_FACTOR matrix keeps track of info needed to
        !         determine the largest stable timestep
        !      ********************************************************************
        L110: DO  J=1,MYY
            L111: DO  I=1,MX
                FLUVIAL_TIMESTEP_FACTOR(I,J)=1.0E+25
                NNN=BASIN_NUMBER(I,J)
                !      ********************************************************************
                !       Skip this point if we are submerged  (depressions are water filled
                !          and we are below outlet level for this basin)
                !      ********************************************************************
                IF (SUBMERGED(I,J)) THEN
                    CYCLE
                ELSE
                    !      ********************************************************************
                    !       The water erosion rate depends upon local area and gradient raised
                    !        to their respective scaling exponents
                    !      ********************************************************************
                    IF (DISCHARGE(I,J) > 0.0) THEN
                        AFACTOR = DISCHARGE(I,J)**BEDROCK_DISCHARGE_EXPONENT
                    ELSE
                        AFACTOR=0.0
                    ENDIF
                    GFACTOR = D8_GRADIENT(I,J)**BEDROCK_GRADIENT_EXPONENT
                    TEMP = AFACTOR*GFACTOR
                    CALL LOCAL_VALUES(I,J,LOCAL_INVERSE_RESISTANCE,NOMINAL_CRITICAL_SHEAR, &
                    LOCAL_CRITICAL_SHEAR, LOCAL_BEDROCK_ERODIBILITY, LOCAL_REGOLITH_CRITICAL_SHEAR,  &
                    LOCAL_REGOLITH_ERODIBILITY)
                    !      ********************************************************************
                    !       Two alternative mechanisms of bed erosion can be specified:
                    !          abrasion by saltating bedload
                    !          and the "classic' stream power model
                    !      ********************************************************************
                    SHEAR_STRESS_FACTOR=TEMPKW*TEMP
                    WHATSUP=.FALSE.
                    IF (USE_SKLAR_BED_ABRASION) THEN
                        !      ********************************************************************
                        !       The formulation for bed erosion by saltating sediment is from
                        !          sklar and dietrich (2003) -- personal communicaiton --
                        !          Assumes transport is not near stage of suspension
                        !      ********************************************************************
                        XAREA=DISCHARGE(I,J)
                        A1TERM=(XAREA)**SEDIMENT_1_EXPONENT
                        A2TERM=(XAREA)**SEDIMENT_2_EXPONENT
                        CALL SEDIMENT_TRANSPORT_FLUX(QSSS,TRANSPORT_STAGE,D8_GRADIENT(I,J),A1TERM,A2TERM)
                        QSSS=QSSS/CHANNEL_WIDTH(I,J)
                        BEDLOAD_FLUX=DMAX1(0.0D0,SEDIMENT_YIELD(I,J)/CHANNEL_WIDTH(I,J))
                        IF (QSSS > 0.0) THEN
                            SEDSATFACTOR=DMAX1(0.0D0,(1.0-(BEDLOAD_FLUX/QSSS)))
                        ELSE
                            SEDSATFACTOR=0.0
                        ENDIF
                        BEDLOAD_FLUX=BEDLOAD_FLUX*BEDLOAD_FRACTION
                        IF (TRANSPORT_STAGE > 0.0) THEN
                            ERODE_CHANNEL(I,J)=-SKLAR_MULT*BEDLOAD_FLUX*SEDSATFACTOR/  &
                            SQRT(TRANSPORT_STAGE)
                        ELSE
                            ERODE_CHANNEL(I,J)=0.0
                        ENDIF
                        IF (ERODE_CHANNEL(I,J) < 0.0) THEN
                            FLUVIAL_TIMESTEP_FACTOR(I,J)=-1.0/ERODE_CHANNEL(I,J)
                        ELSE
                            FLUVIAL_TIMESTEP_FACTOR(I,J)=1.0E+25
                        ENDIF
                    ELSE
                        !      ********************************************************************
                        !       This is the bed erosion formulation for the stream power (shear
                        !          stress) erosion dependency -- no dependency on sediment transport
                        !       We only have water erosion if the local shear stress (temp) exceeds
                        !          the critical shear stress
                        !      ********************************************************************
                        IF (SHEAR_STRESS_FACTOR  >  LOCAL_CRITICAL_SHEAR) THEN
                            ERODE_CHANNEL(I,J) = -(SHEAR_STRESS_FACTOR-LOCAL_CRITICAL_SHEAR)* &
                            LOCAL_BEDROCK_ERODIBILITY * LOCAL_INVERSE_RESISTANCE
                            WHATSUP=.TRUE.
                            !      *********************************************************************
                            !       If the stream power rate is multiplied by the bedload sediment
                            !        humped relationship as suggested by whipple and tucker, then
                            !        do so
                            !      *********************************************************************
                            IF (USE_WHIPPLE_BED_ABRASION) THEN
                                XAREA=DISCHARGE(I,J)
                                A1TERM=(XAREA)**SEDIMENT_1_EXPONENT
                                A2TERM=(XAREA)**SEDIMENT_2_EXPONENT
                                CALL SEDIMENT_TRANSPORT_FLUX(QSSS,TRANSPORT_STAGE,D8_GRADIENT(I,J) &
                                ,A1TERM,A2TERM)
                                SATRATIO=DMAX1(0.0D0,SEDIMENT_YIELD(I,J)/QSSS)
                                IF (QSSS > 0.0) THEN
                                    SEDSATFACTOR=DMAX1(0.0D00,  &
                                    (4.0D00*SATRATIO*(1.0D00-SATRATIO)))
                                ELSE
                                    SEDSATFACTOR=0.0
                                ENDIF
                                ERODE_CHANNEL(I,J)=ERODE_CHANNEL(I,J)*SEDSATFACTOR
                            ENDIF
                            TEMP1=SHEAR_STRESS_FACTOR-LOCAL_CRITICAL_SHEAR
                            IF (TEMP1 > 0.0) THEN
                                FLUVIAL_TIMESTEP_FACTOR(I,J) =D8_GRADIENT(I,J)/   &
                                (LOCAL_BEDROCK_ERODIBILITY*LOCAL_INVERSE_RESISTANCE*BEDROCK_GRADIENT_EXPONENT*TEMP1)
                            ELSE
                                FLUVIAL_TIMESTEP_FACTOR(I,J)=1.0E+25
                            ENDIF
                        ENDIF
                    ENDIF
                    IF (ERODE_CHANNEL(I,J) < MAXCHAN) MAXCHAN=ERODE_CHANNEL(I,J)
                    !     **********************************************************************
                    !       Determine regolith erosion rate
                    !     *****************************************************************
                    REGOLITH_CRITICAL_SHEAR=LOCAL_REGOLITH_CRITICAL_SHEAR
                    IF (SHEAR_STRESS_FACTOR > REGOLITH_CRITICAL_SHEAR) THEN
                        ERODE_REGOLITH_CHANNEL(I,J)=-(SHEAR_STRESS_FACTOR-REGOLITH_CRITICAL_SHEAR)*  &
                        LOCAL_REGOLITH_ERODIBILITY* REGOLITH_ERODIBILITY_FACTOR
                    ENDIF
                    IF (ERODE_REGOLITH_CHANNEL(I,J) < MAXREG) MAXREG=ERODE_REGOLITH_CHANNEL(I,J)
                ENDIF
                IF (DODEBUGWINDOW) THEN
                    IF ((I<61).AND.(I>149).AND.(J<=MY).AND.(J>(MY-11))) THEN
                       WRITE(OUTHIST,792) I,J,ERODE_CHANNEL(I,J),ERODE_REGOLITH_CHANNEL(I,J),REGOLITH_CRITICAL_SHEAR, &
                           SHEAR_STRESS_FACTOR,LOCAL_REGOLITH_ERODIBILITY,REGOLITH_ERODIBILITY_FACTOR,LOCAL_CRITICAL_SHEAR, &
                           LOCAL_BEDROCK_ERODIBILITY,LOCAL_INVERSE_RESISTANCE,AFACTOR,GFACTOR,NOMINAL_CRITICAL_SHEAR
                       792 FORMAT(' I=',I5,' J=',I5,' EC=',G12.5,' ERC=',G12.5,' RCS=',G12.5,' SSF=',G12.5,' LRE=',G12.5,/, &
                           ' REF=',G12.5,' LCS=',G12.5,' LBE=',G12.5,' LIR=',G12.5,' AF=',G12.6,' GF=',G12.5,' NCS=',G12.5)
                   ENDIF
                ENDIF
            ENDDO L111
        ENDDO L110
        RETURN
    END
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE DO_THE_EROSION()
        USE ERODE_GLOBALS, CHANNEL_STATE_CHANGE=>IDO, CHANNEL_EROSION_RATE=>CFNE,  &
        NET_EROSION_RATE=>CFNW,SED_FLUX_DIVERGENCE=>CFW, SED_SURFACE_ELEVATION=>CFN
        USE SEDDEBUG_GLOBALS
        IMPLICIT NONE
        INTEGER(4) :: I,J,NNN,LOCAL_DIRECTION,INM,JNM,NYYY,IP,JP,II,JJ,III,JJJ
        REAL(8) :: SUM_CHANGE,ABSOLUTE_CHANGE,CHANNEL_GRADIENT
        REAL(8) :: ESL,ECC,CWW,RMULT,RFAC,CHANNEL_CHANGE,ECR
        REAL(8) :: TERODE,MAXSEDDIV,MAXSLPDIV
        REAL(8) :: ECOMP,QQQQ,SEDIMENT_FRACTION
        REAL(8) :: ENWX,TEMP
        REAL(8) :: TALERD,TNERODE,TCHANERODE,ROUTECHANGE
        REAL(8) :: SLOPESED,SED_FLUX_DIVERGENCECHANGE
        REAL(8) :: FGSUP,FGDEP,BGSUP,BGDEP,S1,S2,ESTART,EFIN
        REAL(8) :: NEWS2,NEWS1,EDIFFSUM,SEDSUM,NEDIFFSUM
        REAL(8) :: LOCAL_BEDROCK_ERODIBILITY
        REAL(8) :: NOMINAL_CRITICAL_SHEAR
        REAL(8) :: REGFAC,FFFF,OLDREG,REGRATIO
        REAL(8) :: TSED_FLUX_DIVERGENCE
        REAL(8) :: CHMAX,CHMAX1
        REAL(8) :: MAXREG,MAXCHAN,XTEMP
        REAL(8) :: MINCASE1,MAXCASE1,AVGCASE1,NCASE1
        REAL(8) :: MINCASE2,MAXCASE2,AVGCASE2,NCASE2
        REAL(8) :: MINREGRAT,MAXREGRAT,AVGREGRAT,NREGRAT
        REAL(8) :: NBEDROCK,NREGOLITH,NBEDTOREG,NREGTOBED
        REAL(8) :: MINREGFAC,MAXREGFAC,AVGREGFAC,NREGFAC
        REAL(8) :: MINREGOLITH,MAXREGOLITH,AVGREGOLITH,NNREGOLITH
        REAL(8) :: DIFFEL,TOTERODE,TOTDEPOSIT,NTOTDEPOSIT
        REAL(8) ::  TEMP1,TEMP2,BEDLOAD_FLUX
        REAL(8) :: TSEDCHANGE,TSEDSLOPE,ELCHANGE
        REAL(8) :: THEAVG,THEMIN,THEMAX,MAXDISCHARGE,MINECHANNEL,ERODEDEPTH
        INTEGER(4) :: NFG,NBG,NSDELT,NDALLUV,NDDELT
        INTEGER(4) :: NSED_SURFACE_ELEVATION,NROUTECH,NSED_FLUX_DIVERGENCECH,NROUTESED,NNONTOAL,NALTONON
        INTEGER(4) :: ITEMP ,IICOUNT
        LOGICAL(4) :: BACKGRADE,DO_ROUTING,DO_CONVERT
        LOGICAL(4) :: TO_BEDROCK,TO_REGOLITH
        INTEGER(4) :: NNCASE1,NNCASE2,NNCASE3
        INTEGER(4) :: NSCASE1,NSCASE2,NSCASE3
        REAL(8) :: CN1NE,CN1NW,CN2NE,CN2NW,CN3NE,CN3NW
        REAL(8) :: CS1NE,CS1NW,CS2NE,CS2NW,CS3NE,CS3NW
        REAL(8) :: SQSUM,SYSUM,SQYNUM,SGSUM,SWSUM
        REAL(8) :: SSUMCH,SGREQ,XGREQ
        REAL(8) :: XQSUM,XYSUM,XQYNUM,XGSUM,XWSUM,XSUMCH
        REAL(8) :: CMIN,CAVG,DEBUG_CHANGE
        INTEGER(4) :: NPOS,NZER
        LOGICAL(4) :: IS_SEDIMENT_SOURCE,WRITEDETAIL
        INTEGER(4) :: IBT
        REAL(8) :: AVG_SED_FLUX_DIVERGE,NUMBER_SED_FLUX_DIVERGE
        !      ********************************************************************
        !       This subroutine does the actual lanscape erosion as a combination
        !         of mass-wasting, fluvial erosion, and fluvial sediment transport
        !  MODIFIES:  CFN, CFNW, CFNE, CFW, IDO, IS_ROCK_SURFACE, IS_SEDIMENT_COVERED, ACCELERATED_EROSION
        !             PRESENT_TIME, SEDIMENT_BASE, REGOLITTH, MAXIMUM_DISCHARGE, MAXIMUM_SEDIMENT_YIELD
        !             PREVIOUS_ELEVATION, RELATIVE_RESISTANCE, EQUILIBRIUM_GRADIENT, SEDIMENT_FLUX
        !             SLOPEWORK, SLOPEGRAV, ERODEWORK, ERODEGRAV, LAST_TIME_INCREMENT,
        !             TIME_INREMENT, SEDIMENT_DEPOSITED, ELEVATION, CUMULATIVE_ELEVATION_CHANGE
        !  CALLS:  DO_MASS_WASTING, DO_WEATHERING, FLUVIAL_DETACHMENT, SEDMENT_FLUX_DIVERGENCE
        !          ROUTE_DEDMENT, SMOOTHSED, BOUNDARY_CONDITIONS
        !      ********************************************************************
        DEBUG_CHANGE=20.0
        WRITEDETAIL=.FALSE.
        NPOS=0
        NZER=0
        CMIN=1.0E+25
        CAVG=0.0
        SQSUM=0.0
        SGSUM=0.0
        SYSUM=0.0
        SWSUM=0.0
        SGREQ=0.0
        SSUMCH=0.0
        SQYNUM=0.0
        XQSUM=0.0
        XGSUM=0.0
        XYSUM=0.0
        XWSUM=0.0
        XGREQ=0.0
        XQYNUM=0.0
        XSUMCH=0.0
        NNCASE1=0
        NNCASE2=0
        NNCASE3=0
        CN1NE=0.0
        CN1NW=0.0
        CN2NE=0.0
        CN2NW=0.0
        CN3NE=0.0
        CN3NW=0.0
        NSCASE1=0
        NSCASE2=0
        NSCASE3=0
        CS1NE=0.0
        CS1NW=0.0
        CS2NE=0.0
        CS2NW=0.0
        CS3NE=0.0
        CS3NW=0.0
        NNREGOLITH=0.0
        MINREGOLITH=1.0E+25
        MAXREGOLITH=-MINREGOLITH
        AVGREGOLITH=0.0
        NREGFAC=0.0
        MINREGFAC=1.0E+25
        MAXREGFAC=-MINREGFAC
        AVGREGFAC=0.0
        NBEDROCK=0.0
        NREGOLITH=0.0
        NBEDTOREG=0.0
        NREGTOBED=0.0
        MINCASE1=1.0E+25
        MINCASE2=MINCASE1
        MAXCASE1=-MINCASE1
        MAXCASE2=-MINCASE1
        MINREGRAT=MINCASE1
        MAXREGRAT=-MINCASE1
        AVGREGRAT=0.0
        NREGRAT=0.0
        NCASE1=0.0
        NCASE2=0.0
        AVGCASE1=0.0
        AVGCASE2=0.0
        CHMAX=1.0E+25
        CHMAX1=CHMAX
        MAXREG=CHMAX
        MAXCHAN=CHMAX
        ABSOLUTE_CHANGE=0.0
        ESTART=0.0
        EFIN=0.0
        NSDELT=0
        ROUTECHANGE=0.0
        SED_FLUX_DIVERGENCECHANGE=0.0
        NDDELT=0
        NDALLUV=0
        TMULT = 1.0E+10
        MAXSEDDIV=0.0
        MAXSLPDIV=0.0
        SLOPESED=0.0
        FGSUP=0.0
        BGSUP=0.0
        FGDEP=0.0
        BGDEP=0.0
        NFG=0
        NBG=0
        NROUTESED=0
        AVG_SED_FLUX_DIVERGE=0.0
        NUMBER_SED_FLUX_DIVERGE=0.0
        ICASE=0
        !      ********************************************************************
        !       Initialize the erosion rate matrices
        !       matrix channel_state_change(i,j) is used to indicate whether a given matrix point is to be
        !         converted to a different type from its current type
        !         - a value of 0 means that there is no change
        !         - a value of 1 means that it is a bedrock channel converted to sediment
        !              covered (alluvial)
        !         - a value of 2 means that it is an alluvial channel converted to bedrock
        !         The conversions are not actually made until after the topography is
        !         modified
        !      ********************************************************************
        CHANNEL_STATE_CHANGE=0
        ERODE_SLOPE = 0.0
        NOMINAL_ERODE_SLOPE =0.0
        ERODE_REGOLITH_CHANNEL = 0.0
        ERODE_CHANNEL = 0.0
        !      *********************************************************************
        !        Call the subroutine to weather bedrock
        !      *********************************************************************
        CALL DO_WEATHERING()
        !      ********************************************************************
        !       If we are directly modeling mass wasting call the subroutine
        !      ********************************************************************
        IF (DO_MODEL_SLOPES) THEN
            CALL DO_MASS_WASTING()
            IF (ABORT_SIMULATION) RETURN
        ENDIF
        !   ***********************************************************************
        !    Calculate channel bed erosion
        !   ***********************************************************************
        IF (DO_FLUVIAL_DETACHMENT) THEN
            CALL FLUVIAL_DETACHMENT(MAXCHAN,MAXREG,NOMINAL_CRITICAL_SHEAR,LOCAL_BEDROCK_ERODIBILITY)
        ENDIF
        !      ********************************************************************
        !       This determines the sediment transport divergence on sediment-covered
        !          locations - we don't call if do_sediment_transport is false
        !
        !       sed_flux_divergence(i,j) is used to bookkeep transport divergence
        !
        !       First we zero the matrix
        !      ********************************************************************
        IF (DO_SEDIMENT_TRANSPORT) THEN
            SEDIMENT_FLUX=0.0
            EQUILIBRIUM_GRADIENT=0.0
            MAXIMUM_DISCHARGE=DISCHARGE
            MAXIMUM_SEDIMENT_YIELD=SEDIMENT_YIELD
            SED_FLUX_DIVERGENCE=0.0
            CALL SEDIMENT_FLUX_DIVERGENCE()
        ENDIF
        !      ********************************************************************
        !       This sets the default timescale scaling parameter (maximum_channel_timestep) based
        !         upon mass wasting
        !      ********************************************************************
        MAXIMUM_CHANNEL_TIMESTEP=DEFAULT_CHANNEL_TIMESTEP
        !      ********************************************************************
        !       If there is bedrock channel erosion then the maximum timestep for bedrock
        !         channel erosion equals the minimum value of sed_surface_elevation(i,j).  Also calculate
        !         minimum timestep values for slope mass-wasting and (now commented out) for
        !         alluvial channel sediment divergence.
        !       In the next few lines sed_surface_elevation is an alias for &
        !           FLUVIAL_TIMESTEP_FACTOR in Subroutine Fluvial_Detachment
        !       Later in the subroutine sed_surface_elevation means what it says.
        !      ********************************************************************
        IF (DO_FLUVIAL_DETACHMENT) THEN
            IF (LOCAL_BEDROCK_ERODIBILITY > 0.0) THEN
                L115: DO  J = 1, MYY
                    L116: DO  I = 1, MX
                        IF ((D8_GRADIENT(I,J) > GRADCUT).AND.(MAXIMUM_CHANNEL_TIMESTEP  >  SED_SURFACE_ELEVATION(I,J))) &
                        MAXIMUM_CHANNEL_TIMESTEP = SED_SURFACE_ELEVATION(I,J)
                    ENDDO L116
                ENDDO L115
                MAXIMUM_CHANNEL_TIMESTEP=MAXIMUM_CHANNEL_TIMESTEP*CHANNEL_TIMESTEP_SCALING
            ENDIF
        ENDIF
        IF (DO_MODEL_SLOPES) THEN
            L1115: DO  J=1,MYY
                L1116: DO  I=1,MX
                    IF (ABS(ERODE_SLOPE(I,J)) > MAXSLPDIV) &
                    MAXSLPDIV=ABS(ERODE_SLOPE(I,J))
                ENDDO L1116
            ENDDO L1115
            21116       FORMAT('MAXSLOPEDIV1=',G12.5)
            MAXSLPDIV=MASS_WASTING_TIMESTEP_SCALING/MAXSLPDIV
        ELSE
            MAXSLPDIV=1.0E+20
        ENDIF
        IF ((.NOT.DO_FLUVIAL_DETACHMENT).AND.(.NOT.DO_SEDIMENT_TRANSPORT)) THEN
            MAXSLPDIV=MASS_WASTING_TIMESTEP_SCALING/MAXSLPDIV
        ENDIF
        21118 FORMAT('MAXSLOPEDIV3=',G12.5)
        IF (DO_SEDIMENT_TRANSPORT.AND.(MAXSEDDIV > 0.0)) THEN
            MAXSEDDIV=SEDIMENT_TIMESTEP_SCALING/MAXSEDDIV
        ELSE
            MAXSEDDIV=1.0E+20
        ENDIF
        !      *********************************************************************
        !       The actual time increment is the minimum of the timescales for the erosional
        !         processes
        !      *********************************************************************
        LAST_TIME_INCREMENT=TIME_INCREMENT
        TIME_INCREMENT = MIN(MAXIMUM_CHANNEL_TIMESTEP,MAXSLPDIV,MAXSEDDIV,MAXIMUM_TIME_INCREMENT)
        IF (TIME_INCREMENT < MINUMUM_TIME_INCREMENT)TIME_INCREMENT=MINUMUM_TIME_INCREMENT
        IF (ITERATION < 2) THEN
            TIME_INCREMENT=1.0E-5
            LAST_TIME_INCREMENT=TIME_INCREMENT
        ENDIF
        IF (WRITEDETAIL) THEN
            WRITE(OUTHIST,696) MAXIMUM_CHANNEL_TIMESTEP,MAXSLPDIV,MAXSEDDIV,TIME_INCREMENT
            696 FORMAT(' CHAN,SLP,SED LIMITS, TIMESCALE: ',4G12.5)
        ENDIF
        TIMECHANGE=1.0
        !      *********************************************************************
        !        Find the equilibrium gradient for alluvial channels for the given sediment
        !         load and drainage area (discharge surrogate)
        !      *********************************************************************
        IF (DO_SEDIMENT_TRANSPORT) THEN
            N116: DO J=1,MY
                M117: DO I=1,MX
                    !      *********************************************************************
                    !        sed_surface_elevation holds the new alluvial surface elevations calculated by sediment routing
                    !          The value of 1.0e+25 indicates no sediment has so far been routed through that location
                    !          The equilibium alluvial channel gradient is calculated first
                    !      *********************************************************************
                    SED_SURFACE_ELEVATION(I,J)=1.0E+25
                    BEDLOAD_FLUX=SEDIMENT_YIELD(I,J)
                    CALL EQUILIBRIUM_SEDIMENT_GRADIENT(BEDLOAD_FLUX,QQQQ,MAXIMUM_DISCHARGE(I,J))
                    EQUILIBRIUM_GRADIENT(I,J)=QQQQ
                ENDDO M117
            ENDDO N116
            !      *********************************************************************
            !       This loop calls the sediment routing subroutine for locations where a
            !         non-alluvial channel disgorges sediment (measured by sediment_yield(i,j) onto an
            !         alluvial channel
            !      *********************************************************************
            NEWS2=0.0
            NEWS1=0.0
            EDIFFSUM=0.0
            SEDSUM=0.0
            NEDIFFSUM=0.0
            L117: DO  J=1,MYY
                L118: DO  I=1,MX
                    IF ((.NOT.IS_SEDIMENT_COVERED(I,J)).OR.IS_INFLUENT_RIVER_LOCATION(I,J)) THEN
                        IF (ELEVATION(I,J) < OCEAN_ELEVATION) CYCLE L118
                        IP=I
                        JP=J
                        LOCAL_DIRECTION=FLOW_DIRECTION(I,J)
                        IF (IS_INFLUENT_RIVER_LOCATION(I,J).AND.(LOCAL_DIRECTION < 2)) THEN
                            INM=IP
                            JNM=JP
                        ELSE
                            IF (LOCAL_DIRECTION < 2) THEN
                                CYCLE L118
                            ENDIF
                            JNM=JP+DOWNSTREAM(LOCAL_DIRECTION,2)
                            INM=IP+DOWNSTREAM(LOCAL_DIRECTION,1)
                            IF (IS_X_PERIODIC) THEN
                                IF (INM < 1) INM=MX
                                IF (INM > MX) INM=1
                            ENDIF
                            IF (IS_Y_PERIODIC) THEN
                                IF (JNM < 1) JNM=MY
                                IF (JNM > MY) JNM=1
                            ELSE
                                IF (JNM == MY) CYCLE L118
                            ENDIF
                        ENDIF
                        !      *********************************************************************
                        !       We have found a location where a non-alluvial channel impinges on an
                        !             alluvial channel
                        !      *********************************************************************
                        IF (IS_SEDIMENT_COVERED(INM,JNM).OR.&
                        (ELEVATION(INM,JNM) < OCEAN_ELEVATION) &
                        .OR.SUBMERGED(INM,JNM)    ) THEN
                            !      *********************************************************************
                            !       We call for sediment routing only if do_sediment_routing is true and if the amount of
                            !         sediment supplied is large enough
                            !      *********************************************************************
                            DO_ROUTING=.TRUE.
                            IF (DO_ROUTING) THEN
                                IF (DO_SEDIMENT_ROUTING) THEN
                                    NROUTESED=NROUTESED+1
                                    IS_SEDIMENT_SOURCE=.FALSE.
                                    IF (HAVE_INFLUENT_RIVERS.AND.IS_INFLUENT_RIVER_LOCATION(IP,JP))THEN
                                        IS_SEDIMENT_SOURCE=.TRUE.
                                    ENDIF
                                    IF (DO_FLOW_BOUNDARIES) THEN
                                        IF(FLOW_DIRECTION(I,J) == 1) THEN
                                        ELSE
                                            CALL ROUTE_SEDIMENT(I,J,TIME_INCREMENT,&
                                            BACKGRADE,IS_SEDIMENT_SOURCE,S1,S2)
                                        ENDIF
                                    ELSE
                                        CALL ROUTE_SEDIMENT(I,J,TIME_INCREMENT,&
                                        BACKGRADE,IS_SEDIMENT_SOURCE,S1,S2)
                                    ENDIF
                                    IF (IS_SEDIMENT_SOURCE) THEN
                                        IF (BACKGRADE) THEN
                                            ITEMP=1
                                        ELSE
                                            ITEMP=0
                                        ENDIF
                                        XTEMP=LAST_TIME_INCREMENT*SEDIMENT_YIELD(I,J)
                                    ENDIF
                                    !      *********************************************************************
                                    !       nontoa1, nbg,bgsup,bgdep,nfg,fgsup & fgdep are for bookkeeping the amount of
                                    !           deposited and eroded sediment - just for debugging purposes
                                    !      *********************************************************************
                                    NEWS2=NEWS2+S2*TIME_INCREMENT
                                    NEWS1=NEWS1+S1*TIME_INCREMENT
                                    IF (ISNAN(NEWS2).OR.ISNAN(NEWS1)) THEN
                                        WRITE(OUTHIST,8432) I,J,NEWS1,NEWS2,S1,S2 &
                                        ,TIME_INCREMENT,SEDIMENT_YIELD(I,J),SEDBIAS
                                        WRITE(*,8432) I,J,NEWS1,NEWS2,S1,S2 &
                                        ,TIME_INCREMENT,SEDIMENT_YIELD(I,J),SEDBIAS
                                        8432 FORMAT(' BAD SED ROUTING,I=',I5,' J=',I5,' NEWS1=',G12.5,' NEWS2=' &
                                        ,G12.5,/,' S1=',G12.5,' S2=',G12.5,' TINC=',G12.5,&
                                        ' SEDIMENT_YIELD=',G12.5,' SEDBIAS=',G12.5)
                                        WRITE(*,8431)
                                        8431 FORMAT(' IN THE ROUTING')
                                        ABORT_SIMULATION=.TRUE.
                                        RETURN
                                    ENDIF
                                ENDIF
                            ENDIF
                        ENDIF
                    ENDIF
                ENDDO L118
            ENDDO L117
            !      *********************************************************************
            !        Here we convert the changed elevations (sed_surface_elevation(i,j) due to sediment routing
            !         into equivalent sediment flux divergences (sed_flux_divergence(i,j))
            !      *********************************************************************
            IF (DO_SEDIMENT_ROUTING) THEN
                IF (DO_ALLUVIAL_SMOOTHING) THEN
                    L178: DO J=1,MYY
                        L179: DO I=1,MX
                            IF (SED_SURFACE_ELEVATION(I,J) < 1.0E+24) THEN
                                IF (SED_SURFACE_ELEVATION(I,J) > ELEVATION(I,J)) THEN
                                    SEDIMENT_DEPOSITED(I,J)=.TRUE.
                                    IF (SUBMERGED(I,J)) SEDIMENT_DEPOSITED(I,J)=.FALSE.
                                ELSE
                                    SEDIMENT_DEPOSITED(I,J)=.FALSE.
                                ENDIF
                            ELSE
                                SEDIMENT_DEPOSITED(I,J)=.FALSE.
                            ENDIF
                        ENDDO L179
                    ENDDO L178
                    CALL SMOOTHSED()
                ENDIF
                NEWS1=NEWS1
                NEWS2=NEWS2
                IF (WRITEDETAIL) THEN
                    WRITE(OUTHIST,1778) NEWS1,NEWS2
                    1778  FORMAT(' S1SUM=',G12.5,' S2SUM=',G12.5)
                ENDIF
            ENDIF
            NROUTECH=0
            NSED_FLUX_DIVERGENCECH=0
            NSED_SURFACE_ELEVATION=0
            L106: DO  J=1,MYY
                L107: DO  I=1,MX
                    !      *********************************************************************
                    !        In order to maintain stability when using
                    !             routing, we disable local divergence for positive values, letting
                    !             routing do the job - this is a kludge
                    !      *********************************************************************
                    IF (SED_FLUX_DIVERGENCE(I,J) > 0.0) SED_FLUX_DIVERGENCE(I,J)=0.0
                    !      *********************************************************************
                    !         If some sediment has been routed through this cell then we know that
                    !         because sed_surface_elevation(i,j) < 1.0e24
                    !      *********************************************************************
                    IF (SED_SURFACE_ELEVATION(I,J) < 1.0E24) THEN
                        SED_FLUX_DIVERGENCE(I,J)=(SED_SURFACE_ELEVATION(I,J)-ELEVATION(I,J))/TIME_INCREMENT
                        AVG_SED_FLUX_DIVERGE=AVG_SED_FLUX_DIVERGE+SED_FLUX_DIVERGENCE(I,J)
                        NUMBER_SED_FLUX_DIVERGE=NUMBER_SED_FLUX_DIVERGE+1.0
                    ENDIF
                ENDDO L107
            ENDDO L106
            IF (NUMBER_SED_FLUX_DIVERGE > 0.0) AVG_SED_FLUX_DIVERGE=AVG_SED_FLUX_DIVERGE/NUMBER_SED_FLUX_DIVERGE
        ENDIF
        !      ********************************************************************
        !       This determines baselevel lowering for simulations where the southern matrix
        !         edge at j=my is lowered at a constant rate
        !         Also the time increment is determined
        !         Also update the erosional and depositional work against gravity and distance
        !      ********************************************************************
        TIME_INCREMENT=TIME_INCREMENT*TIMECHANGE
        SLOPEWORK=SLOPEWORK+SLOPETOADD*TIME_INCREMENT
        ERODEWORK=ERODEWORK+ERODETOADD*TIME_INCREMENT
        SLOPEGRAV=SLOPEGRAV+SLGRAVTOADD*TIME_INCREMENT
        ERODEGRAV=ERODEGRAV+ERGRAVTOADD*TIME_INCREMENT
        !      *********************************************************************
        !         This zeros out matrices used to bookkeep the rate of channel erosion
        !         (channel_erosion_rate(i,j) and the overall erosion amount (net_erosion_rate(i,j))
        !      *********************************************************************
        CHANNEL_EROSION_RATE=0.0
        NET_EROSION_RATE=0.0
        !      ********************************************************************
        !       In this loop we finally determine the erosion rate based upon
        !        echannel, erode_slope, and alluvial sediment flux divergence
        !      ********************************************************************
        L120: DO  J = 1 , MYY
            M120: DO  I = 1 , MX
                !      *********************************************************************
                !        Define local variables and the overall (potential) bedrock channel erosion
                !             rate (channel_change) based on esl,ecc,regolith_erodibility_factor,and channel width
                !        The rate depends upon relative values of slope divergence and channel
                !             erosion potential
                !      *********************************************************************
                TO_BEDROCK=.FALSE.
                TO_REGOLITH=.TRUE.
                IF (EXPLICIT_CHANNEL_BED_STATE) THEN
                    REGFAC = -REGOLITH(I,J)
                    IF (.NOT.IS_ROCK_SURFACE(I,J)) THEN
                        REGFAC = REGOLITH(I,J)/TIME_INCREMENT
                        IF (IS_SEDIMENT_COVERED(I,J)) REGFAC=MAX(REGFAC, &
                        (ELEVATION(I,J)-SEDIMENT_BASE(I,J)))
                        NNREGOLITH=NNREGOLITH+1.0
                        AVGREGOLITH=AVGREGOLITH+REGOLITH(I,J)
                        IF (REGOLITH(I,J) > MAXREGOLITH) &
                        MAXREGOLITH=REGOLITH(I,J)
                        IF (REGOLITH(I,J) < MINREGOLITH) &
                        MINREGOLITH=REGOLITH(I,J)
                    ENDIF
                    ESL = ERODE_SLOPE(I,J)
                    ECC = ERODE_CHANNEL(I,J)
                    ECR = ERODE_REGOLITH_CHANNEL(I,J)
                    CWW = CHANNEL_WIDTH(I,J)
                    IF (ECR /= 0.0) THEN
                        REGRATIO=ECC/ECR
                    ELSE
                        REGRATIO=1.0
                    ENDIF
                    AVGREGRAT=AVGREGRAT+REGRATIO
                    IF (REGRATIO > MAXREGRAT) MAXREGRAT=REGRATIO
                    IF (REGRATIO < MINREGRAT) MINREGRAT=REGRATIO
                    NREGRAT=NREGRAT+1.0
                    FFFF=(1.0-REGRATIO)
                    RFAC = CWW/CELL_SIZE
                    RFAC = DMIN1(RFAC,1.0D0)
                    !      ********************************************************************
                    !           Case i
                    !      ********************************************************************
                    IF (.NOT.IS_ROCK_SURFACE(I,J)) THEN
                        CHANNEL_CHANGE = ESL+ECR*RFAC
                        CHANNEL_EROSION_RATE(I,J)=-CHANNEL_CHANGE+ESL
                        NSCASE1=NSCASE1+1
                        IF (CHANNEL_CHANGE < (-REGFAC)) THEN
                            TO_REGOLITH=.FALSE.
                            NSCASE2=NSCASE2+1
                            NSCASE1=NSCASE1-1
                            !      ********************************************************************
                            !           Case ii
                            !      ********************************************************************
                            CHANNEL_CHANGE= ECC*RFAC-REGFAC*FFFF+ESL*REGRATIO
                            TO_BEDROCK=.TRUE.
                            !      *******************************************************************
                            !        Case iii for regolith slopes
                            !      *******************************************************************
                            IF (ESL < (-REGFAC)) THEN
                                CHANNEL_CHANGE=ESL+ECC*RFAC
                            ENDIF
                            CHANNEL_EROSION_RATE(I,J)=-CHANNEL_CHANGE+ESL
                        ENDIF
                    ELSE
                        !      *******************************************************************
                        !        Case iii for bedrock slopes
                        !      *******************************************************************
                        CHANNEL_CHANGE=ESL+ECC*RFAC
                        CHANNEL_EROSION_RATE(I,J)=-ECC*RFAC
                        NSCASE3=NSCASE3+1
                    ENDIF
                ELSE
                    !      *********************************************************************
                    !       This is the code for implicit treatment of bedrock and regolith
                    !         channel erosion
                    !      *********************************************************************
                    REGFAC = 0.0
                    IF (IS_SEDIMENT_COVERED(I,J)) REGFAC=MIN(REGFAC,  &
                    (SEDIMENT_BASE(I,J)-ELEVATION(I,J)))/TIME_INCREMENT
                    ESL = ERODE_SLOPE(I,J)
                    ECC = ERODE_CHANNEL(I,J)
                    ECR = ERODE_REGOLITH_CHANNEL(I,J)
                    CWW = CHANNEL_WIDTH(I,J)
                    IF (ECR /= 0.0) THEN
                        REGRATIO=ECC/ECR
                    ELSE
                        REGRATIO=1.0
                    ENDIF
                    FFFF=(1.0-REGRATIO)
                    RFAC = CWW/CELL_SIZE
                    RFAC = DMIN1(RFAC,1.0D0)
                    RMULT = RFAC+REGRATIO*(1.0-RFAC)
                    !      ********************************************************************
                    !           Case i
                    !      ********************************************************************
                    CHANNEL_CHANGE = ESL+ECR*RFAC
                    NNCASE1=NNCASE1+1
                    IF (CHANNEL_CHANGE < REGFAC) THEN
                        !      ********************************************************************
                        !           Case ii
                        !      ********************************************************************
                        CHANNEL_CHANGE=(RFAC*(REGFAC+ECC)+REGRATIO*ESL)/RMULT
                        NNCASE1=NNCASE1-1
                        NNCASE2=NNCASE2+1
                        IF (ESL < REGFAC) THEN
                            !      ********************************************************************
                            !           Case iii
                            !      ********************************************************************
                            CHANNEL_CHANGE=ESL+ECC*RFAC/RMULT
                            NNCASE2=NNCASE2-1
                            NNCASE3=NNCASE3+1
                        ENDIF
                    ENDIF
                    CHANNEL_EROSION_RATE(I,J)=-CHANNEL_CHANGE+ESL
                ENDIF
                !      ********************************************************************
                !       Look first at non-sediment covered areas
                !      ********************************************************************
                IF (.NOT.IS_SEDIMENT_COVERED(I,J)) THEN
                    !      *********************************************************************
                    !        sum_change is the net rate of change of elevation from all processes
                    !      *********************************************************************
                    SUM_CHANGE=CHANNEL_CHANGE
                    IF (DO_SEDIMENT_TRANSPORT) THEN
                        !      ********************************************************************
                        !       Check for the case where the net sediment transporting capacity is
                        !         less than the sediment yield -- if so, convert to sediment cover by
                        !             setting channel_state_change(i,j) to 1
                        !      ********************************************************************
                        CHANNEL_GRADIENT=D8_GRADIENT(I,J)
                        QQQQ=SEDIMENT_FLUX(I,J)
                        NNN=BASIN_NUMBER(I,J)
                        IF (.NOT. SUBMERGED(I,J)) THEN
                            IF ((QQQQ) < SEDIMENT_YIELD(I,J)) THEN
                                CHANNEL_STATE_CHANGE(I,J)=1
                                ICASE(13)=ICASE(13)+1
                                IF (ITERATION == 8) THEN
                                    IF (MOD(ICASE(13),100) == 0) THEN
                                        IJDEBUG=IJDEBUG+1
                                        IDEBUG(IJDEBUG)=I
                                        JDEBUG(IJDEBUG)=J
                                    ENDIF
                                ENDIF
                                SQSUM=SQSUM+QQQQ
                                SGSUM=SGSUM+D8_GRADIENT(I,J)
                                SWSUM=SWSUM+DISCHARGE(I,J)
                                SYSUM=SYSUM+SEDIMENT_YIELD(I,J)
                                SSUMCH=SSUMCH+SUM_CHANGE
                                SGREQ=SGREQ+EQUILIBRIUM_GRADIENT(I,J)
                                SQYNUM=SQYNUM+1.0
                            ENDIF
                        ENDIF
                        IF ((DO_SEDIMENT_ROUTING).AND.(SED_SURFACE_ELEVATION(I,J) < 1.0E+24)) THEN
                            SUM_CHANGE=(SED_SURFACE_ELEVATION(I,J)-ELEVATION(I,J))*TIMECHANGE/ &
                            TIME_INCREMENT
                            CHANNEL_EROSION_RATE(I,J)=0.0
                            IF (SED_SURFACE_ELEVATION(I,J) > ELEVATION(I,J)) THEN
                                CHANNEL_STATE_CHANGE(I,J)=1
                                ICASE(14)=ICASE(14)+1
                            ENDIF
                        ENDIF
                    ENDIF
                ELSE
                    IF (DO_SEDIMENT_TRANSPORT) THEN
                        SUM_CHANGE=SED_FLUX_DIVERGENCE(I,J)
                        !      *********************************************************************
                        !        This adds in the contribution to deposition or erosion from slope erosion
                        !         slopesed is for debugging - keeping track of slope erosion or deposition
                        !         on alluvial channels
                        !      *********************************************************************
                        SUM_CHANGE=SUM_CHANGE+ESL
                        SLOPESED=SLOPESED+ESL
                        !      ********************************************************************
                        !       enwx is the new elevation at this point
                        !      ********************************************************************
                        ENWX=ELEVATION(I,J)+SUM_CHANGE*TIME_INCREMENT
                        !      ********************************************************************
                        !       This is the case where we remove all sediment and reconvert to
                        !           regolith by settin channel_state_change(i,j) to 2
                        !       We limit the overall erosion rate in this case to the rate of bedrock
                        !         channel erosion (channel_change)
                        !      ********************************************************************
                        IF ((ENWX < SEDIMENT_BASE(I,J)).AND.(CHANNEL_CHANGE < SUM_CHANGE)) THEN
                            IF (.NOT.SUBMERGED(I,J)) THEN
                                CHANNEL_STATE_CHANGE(I,J)=2
                                ICASE(15)=ICASE(15)-1
                                SUM_CHANGE=CHANNEL_CHANGE
                            ENDIF
                        ENDIF
                        !       *********************************************************************
                        !        This is a kludge to prevent sediment-covered areas below threshold from
                        !         remaining uneroded forever
                        !       *********************************************************************
                        IF (DO_SEDIMENT_ROUTING) THEN
                            IF ((SED_FLUX_DIVERGENCE(I,J) == 0.0).AND. &
                            (SEDIMENT_BASE(I,J) == ELEVATION(I,J))) THEN
                                SUM_CHANGE=CHANNEL_CHANGE
                                IF (.NOT.SUBMERGED(I,J)) THEN
                                    !                         channel_state_change(i,j)=2
                                    ICASE(16)=ICASE(16)-1
                                    XQSUM=XQSUM+SEDIMENT_FLUX(I,J)
                                    XGSUM=XGSUM+D8_GRADIENT(I,J)
                                    XWSUM=XWSUM+DISCHARGE(I,J)
                                    XYSUM=XYSUM+SEDIMENT_YIELD(I,J)
                                    XSUMCH=XSUMCH+SUM_CHANGE
                                    XGREQ=XGREQ+EQUILIBRIUM_GRADIENT(I,I)
                                    XQYNUM=XQYNUM+1.0
                                ENDIF
                            ENDIF
                        ENDIF
                        !      *********************************************************************
                        !        End of kludge and end of sediment covered locations
                        !      *********************************************************************
                    ENDIF
                ENDIF
                !      *********************************************************************
                !         Here we have done all of our erosion rate determinations and we are
                !         actually doing the erosion
                !         If we are converting to sediment-covered at this location, make the base
                !         of the sedimentary deposit equal to the present elevation prior to
                !         deposition
                !      *********************************************************************
                NET_EROSION_RATE(I,J)=SUM_CHANGE
                IF ((CHANNEL_STATE_CHANGE(I,J) == 1).AND.(.NOT.SUBMERGED(I,J))) THEN
                    SEDIMENT_BASE(I,J)=MIN(SEDIMENT_BASE(I,J),ELEVATION(I,J))
                ENDIF
                IF (ABS(SUM_CHANGE) > ABSOLUTE_CHANGE) THEN
                    ABSOLUTE_CHANGE=ABS(SUM_CHANGE)
                    IF (ABSOLUTE_CHANGE > DEBUG_CHANGE) THEN
                       CALL PRINT_AROUND(I,J)
                    ENDIF
                    MAXIMUM_ELEVATION_CHANGE=SUM_CHANGE
                    IDET=I
                    JDET=J
                ENDIF
                IF (ABS(SUM_CHANGE) > DEBUG_CHANGE) THEN
                    CALL PRINT_AROUND(I,J)
                ENDIF
                IF (DO_FLOW_BOUNDARIES) THEN
                    IF (FLOW_DIRECTION(I,J) /= 1) THEN
                        ELCHANGE=SUM_CHANGE*TIME_INCREMENT
                        ELEVATION(I,J) = ELEVATION(I,J) +ELCHANGE
                        CUMULATIVE_ELEVATION_CHANGE(I,J)=CUMULATIVE_ELEVATION_CHANGE(I,J)+ELCHANGE
                    ENDIF
                ELSE
                    ELCHANGE=SUM_CHANGE*TIME_INCREMENT
                    ELEVATION(I,J) = ELEVATION(I,J) +ELCHANGE
                    CUMULATIVE_ELEVATION_CHANGE(I,J)=CUMULATIVE_ELEVATION_CHANGE(I,J)+ELCHANGE
                ENDIF

                IF (SUBMERGED(I,J).OR.IS_SEDIMENT_COVERED(I,J)) THEN
                ELSE
                    IF (IS_ROCK_SURFACE(I,J)) THEN
                        NBEDROCK=NBEDROCK+1.0
                        !      ********************************************************************
                        !       If a bedrock slope element receives more input of regolith
                        !         from mass wasting plus local weathering than it
                        !         can potentially export, then it is converted to regolith
                        !       ******** NOW DISABLED *********
                        !      ********************************************************************
                        !  TEMP1=NOMINAL_ERODE_SLOPE(I,J)+ECR*RFAC
                        !  IF (SUM_CHANGE <= TEMP1) THEN
                        !                       is_rock_surface(i,j)=.false.
                        !                       regolith(i,j)=0.0
                        !                       write(outhist,5443) nominal_erode_slope,temp1,sum_change
                        !  5443  FORMAT(' TO REG 0, NOMINAL_ERODE_SLOPE=',G12.5,' T1=',G12.5,' SC=',G12.5)
                        !  ENDIF
                    ELSE
                        !      ********************************************************************
                        !       If a regolith-mantled slope erodes more rapidly than the potential
                        !         weathering rate, then that location is converted to bedrock
                        !      *******************************************************************
                        NREGOLITH=NREGOLITH+1.0
                        OLDREG=REGOLITH(I,J)
                        REGOLITH(I,J)=REGOLITH(I,J)+ELCHANGE
                        IF (REGOLITH(I,J) <= 0.0) THEN
                            REGOLITH(I,J)=0.0
                        ENDIF
                    ENDIF
                ENDIF
                !      *********************************************************************
                !         If this is a sediment covered area and we lower the elevation below the
                !         present sedimentary base, then we lower the sedimentary base - for
                !         non-sedimentary areas the sediment base is the present elevation
                !      *********************************************************************
                IF (.NOT.SUBMERGED(I,J)) THEN
                    IF (ELEVATION(I,J) < SEDIMENT_BASE(I,J))  THEN
                        SEDIMENT_BASE(I,J)=ELEVATION(I,J)
                    ENDIF
                ENDIF
                IF(MOD((ITERATION+1),OUTPUT_PRINT_INTERVAL) == 0) THEN
                    CALL WRITE_DEBUG(I,J,ESL,ECC,ECR,CHANNEL_CHANGE,SUM_CHANGE)
                ENDIF
            ENDDO M120
        ENDDO L120
        IF ((.NOT.IS_Y_PERIODIC).AND.(.NOT.NO_FLUX_LOWER_BOUNDARY)) THEN
            DO I=1,MX
                IF ((ELEVATION(I,MY-1)-ELEVATION(I,MY)) > (DELTA_FORESET_GRADIENT*CELL_SIZE)) THEN
                    ELEVATION(I,MY-1)=ELEVATION(I,MY)+DELTA_FORESET_GRADIENT*CELL_SIZE
                ENDIF
            ENDDO
        ENDIF
        IF (WRITEDETAIL) THEN
            WRITE(OUTHIST,51123) NSCASE1,NSCASE2,NSCASE3,NNCASE1, &
            NNCASE2,NNCASE3
            51123  FORMAT(6I8)
        ENDIF
        !      *********************************************************************
        !       Here we apply the boundary conditions
        !      *********************************************************************
        CALL BOUNDARY_CONDITIONS()
        !      *********************************************************************
        !       We increment the time and then do some bookkeeping
        !      *********************************************************************
        PRESENT_TIME = PRESENT_TIME +TIME_INCREMENT
        IF (DO_SEDIMENT_TRANSPORT) THEN
            TOTERODE=0.0
            TOTDEPOSIT=0.0
            NTOTDEPOSIT=0.0
            TERODE=0.0
            SEDIMENT_FRACTION=0.0
            TALERD=0.0
            TNERODE=0.0
            TCHANERODE=0.0
            TSED_FLUX_DIVERGENCE=0.0
            TSEDSLOPE=0.0
            TSEDCHANGE=0.0
            NYYY=0
            L735: DO  I=1,MX
                M735: DO  J=1,MYY
                    DIFFEL=ELEVATION(I,J)-INITIAL_ELEVATION(I,J)
                    IF (DIFFEL > 0.0) THEN
                        NTOTDEPOSIT=NTOTDEPOSIT+1
                        TOTDEPOSIT=TOTDEPOSIT+DIFFEL
                    ELSE
                        TOTERODE=TOTERODE+DIFFEL
                    ENDIF
                    TERODE=TERODE+NET_EROSION_RATE(I,J)
                    IF (IS_SEDIMENT_COVERED(I,J)) THEN
                        SEDIMENT_FRACTION=SEDIMENT_FRACTION+1.0
                        NYYY=NYYY+1
                        TALERD=TALERD+NET_EROSION_RATE(I,J)
                        TSED_FLUX_DIVERGENCE=TSED_FLUX_DIVERGENCE+SED_FLUX_DIVERGENCE(I,J)
                        TSEDSLOPE=TSEDSLOPE+ERODE_SLOPE(I,J)
                        TSEDCHANGE=TSEDCHANGE+NET_EROSION_RATE(I,J)-ERODE_SLOPE(I,J)
                    ELSE
                        TNERODE=TNERODE+NET_EROSION_RATE(I,J)
                        TCHANERODE=TCHANERODE+CHANNEL_EROSION_RATE(I,J)
                    ENDIF
                ENDDO M735
            ENDDO L735
            TERODE=TERODE*TIME_INCREMENT
            TALERD=TALERD*TIME_INCREMENT
            TCHANERODE=TCHANERODE*TIME_INCREMENT
            TNERODE=TNERODE*TIME_INCREMENT
            TSED_FLUX_DIVERGENCE=TSED_FLUX_DIVERGENCE*TIME_INCREMENT
            TSEDSLOPE=TSEDSLOPE*TIME_INCREMENT
            TSEDCHANGE=TSEDCHANGE*TIME_INCREMENT
            TERODE=TERODE*D2X
            TALERD=TALERD*D2X
            TSEDSLOPE=TSEDSLOPE*D2X
            TSEDCHANGE=TSEDCHANGE*D2X
            TNERODE=TNERODE*D2X
            TCHANERODE=TCHANERODE*D2X
            !     #####################################################
            !       experimental procedure to correct depositional bias
            !      ####################################################
            IF ((FORCE_SEDIMENT_CONSERVATION).AND.(NYYY > 0)) THEN
                IF (ITERATION > 1) THEN
                    SEDBIAS=BIAS_PARAMETER*NEWS1/TSEDCHANGE+  &
                    (1.0-BIAS_PARAMETER)*LASTSEDBIAS
                ELSE
                    SEDBIAS=NEWS1/TSEDCHANGE
                ENDIF
                IF (SEDBIAS < 0.9) SEDBIAS=0.9
                IF (SEDBIAS > 1.3) SEDBIAS=1.3
                LASTSEDBIAS=SEDBIAS
                WRITE(OUTHIST,2731) SEDBIAS
                2731  FORMAT(' SEDBIAS CORRECTION=',G12.5)
            ENDIF
            2732        CONTINUE
            !     #####################################################
            !       end of experimental procedure to correct depositional bias
            !      ####################################################
            IF(MOD((ITERATION+1),OUTPUT_PRINT_INTERVAL) == 0) THEN
                WRITE(OUTHIST,731) TERODE
                731         FORMAT(' TOTAL ELEVATION CHANGE=',G12.5)
                WRITE(OUTHIST,732) TALERD,NYYY
                732         FORMAT(' TOTAL ALLUVIAL SURFACE CHANNEL EROSION=',G12.5, &
                ' ON',I6,' POINTS')
                WRITE(OUTHIST,21732) TSEDSLOPE,TSEDCHANGE
                21732         FORMAT(' SLOPE EROSION ON SED. AREAS=',G12.5,/, &
                ' ALLUVIAL SURFACE CHANGE LESS SLOPE EROSION=',G12.5)
                WRITE(OUTHIST,738) TNERODE,TCHANERODE
                738           FORMAT(' TOTAL NON-SEDIMENTARY AREA EROSION=',G12.5,/, &
                ' TOTAL CHANNEL EROSION=',G12.5)
                SEDIMENT_FRACTION=SEDIMENT_FRACTION/(MX*MY)
                WRITE(OUTHIST,734)SEDIMENT_FRACTION
                734           FORMAT(' FRACTION OF SEDIMENT-COVERED POINTS=',G12.5)
            ENDIF
            IF (WRITEDETAIL) THEN
                WRITE(OUTHIST,731) TERODE
                WRITE(OUTHIST,732) TALERD,NYYY
                WRITE(OUTHIST,21732) TSEDSLOPE,TSEDCHANGE
                WRITE(OUTHIST,738) TNERODE,TCHANERODE
                SEDIMENT_FRACTION=SEDIMENT_FRACTION/(MX*MY)
                WRITE(OUTHIST,734)SEDIMENT_FRACTION
            ENDIF
            IF (ISNAN(NEWS2).OR.ISNAN(NEWS1)) THEN
                WRITE(OUTHIST,8332) I,J,NEWS1,NEWS2,S1,S2 &
                ,TIME_INCREMENT,SEDIMENT_YIELD(I,J),SEDBIAS
                WRITE(*,8332) I,J,NEWS1,NEWS2,S1,S2 &
                ,TIME_INCREMENT,SEDIMENT_YIELD(I,J),SEDBIAS
                8332 format(' bad sed routing,i=',i5,' j=',i5,' news1=',g12.5,' news2=' &
                ,g12.5,/,' tinc=',g12.5,' sediment_yield=',g12.5,' sedbias=',g12.5)
                STOP
            ENDIF
            IF (WRITEDETAIL) THEN
                WRITE(OUTHIST,27341) NEWS1,NEWS2
               27341 FORMAT(' SED INPUT TO ROUTE=',G12.5,' SED ROUTED=',G12.5)
            ENDIF
            TEMP1=-TNERODE
            TEMP2=-TOTERODE
            IF (WRITEDETAIL) THEN
                WRITE(OUTSEDDEBUG,7746) ITERATION,NEWS1,NEWS2, &
                TCHANERODE,TEMP1,TALERD,NYYY, &
                TOTDEPOSIT,NTOTDEPOSIT,TEMP2,TSED_FLUX_DIVERGENCE
            ENDIF
            7746  FORMAT(I13,5(G12.5,' '),I13,4(G12.5,' '))
            7744  FORMAT(' SDVOL=',G12.5,' NSED=',G12.5,' SEDDIF=',G12.5)
            7745  FORMAT(' TODEP=',G12.5,' NTOT=',G12.5,' TOTERODE=',G12.5)
            IF(MOD((ITERATION+1),OUTPUT_PRINT_INTERVAL) == 0) THEN
                WRITE(OUTHIST,433)
                433         FORMAT(/,' ELEVATION')
                CALL SUMMARIZE_MATRIX_DATA(ELEVATION,THEAVG,THEMAX,THEMIN)
                WRITE(OUTHIST,1433)
                1433          FORMAT(/,' AREA')
                CALL SUMMARIZE_MATRIX_DATA(DRAINAGE_AREA,THEAVG,THEMAX,THEMIN)
                WRITE(OUTHIST,1434)
                1434          FORMAT(/,' SURFACE DISCHARGE')
                CALL SUMMARIZE_MATRIX_DATA(DISCHARGE,THEAVG,MAXDISCHARGE,THEMIN)
                WRITE(OUTHIST,434)
                434         FORMAT(/,' ECHANNEL')
                CALL SUMMARIZE_MATRIX_DATA(ERODE_CHANNEL,THEAVG,THEMAX,MINECHANNEL)
                WRITE(OUTHIST,435)
                435         FORMAT(/,' SEDIMENT DIVERGENCE')
                CALL SUMMARIZE_MATRIX_DATA(SED_FLUX_DIVERGENCE,THEAVG,THEMAX,THEMIN)
                WRITE(OUTHIST,436)
                436         FORMAT(/,' SEDIMENT YIELD')
                CALL SUMMARIZE_MATRIX_DATA(SEDIMENT_YIELD,THEAVG,THEMAX,THEMIN)
                WRITE(OUTHIST,437)
                437         FORMAT(/,' ERODE_SLOPE')
                CALL SUMMARIZE_MATRIX_DATA(ERODE_SLOPE,THEAVG,THEMAX,THEMIN)
                WRITE(OUTHIST,438)
                438         FORMAT(/,' GRADIENT')
                CALL SUMMARIZE_MATRIX_DATA(D8_GRADIENT,THEAVG,THEMAX,THEMIN)
                WRITE(OUTHIST,439)
                439         FORMAT(' SEDBASE')
                CALL SUMMARIZE_MATRIX_DATA(SEDIMENT_BASE,THEAVG,THEMAX,THEMIN)
                WRITE(OUTHIST,440)
                440         FORMAT(' SEDIMENT_FLUX')
                CALL SUMMARIZE_MATRIX_DATA(SEDIMENT_FLUX,THEAVG,THEMAX,THEMIN)
                WRITE(OUTHIST,441)
                441         FORMAT(' EQUILIBRIUM_GRADIENT')
                CALL SUMMARIZE_MATRIX_DATA(EQUILIBRIUM_GRADIENT,THEAVG,THEMAX,THEMIN)
                WRITE(OUTHIST,442)
                442         FORMAT(' EROSION RATE')
                CALL SUMMARIZE_MATRIX_DATA(NET_EROSION_RATE,THEAVG,THEMAX,THEMIN)
                WRITE(OUTHIST,443)
                443         FORMAT(' CHANNEL EROSION')
                CALL SUMMARIZE_MATRIX_DATA(CHANNEL_EROSION_RATE,THEAVG,THEMAX,THEMIN)
                IF (MODEL_GROUNDWATER) THEN
                    WRITE(OUTHIST,1444)
                    1444          FORMAT(/,' GROUNDWATER DIVERGENCE')
                    CALL SUMMARIZE_MATRIX_DATA(GROUNDWATER_FLUX,THEAVG,THEMAX,THEMIN)
                    WRITE(OUTHIST,1442)
                    1442          FORMAT(/,' MAX DIVERGENCE')
                    CALL SUMMARIZE_MATRIX_DATA(FILTERED_GROUNDWATER_FLUX,THEAVG,THEMAX,THEMIN)
                ENDIF
                WRITE(OUTHIST,1439)
                1439          FORMAT(/,' REGOLITH')
                CALL SUMMARIZE_REGOLITH_DATA(REGOLITH)
                IICOUNT=0
                DO J=2,MY-1
                    DO I=2,MX-1
                        IICOUNT=IICOUNT+1
                        IF ((MOD(IICOUNT,10) == 0).OR.  &
                        (DISCHARGE(I,J) > (MAXDISCHARGE/10.0)) &
                        ) THEN
                            IF (.NOT.SUBMERGED(I,J)) THEN
                                IF (D8_GRADIENT(I,J) > 0.0) THEN
                                    TEMP=-ERODE_CHANNEL(I,J)
                                    ERODEDEPTH=INITIAL_ELEVATION(I,J)-ELEVATION(I,J)
                                    WRITE(78,1447) ITERATION,I,J,ERODEDEPTH,DISCHARGE(I,J), &
                                    D8_GRADIENT(I,J),TEMP,ERODE_SLOPE(I,J),NET_EROSION_RATE(I,J), &
                                    CHANNEL_EROSION_RATE(I,J),SEDIMENT_FLUX(I,J),SEDIMENT_YIELD(I,J),EQUILIBRIUM_GRADIENT(I,J)
                                ENDIF
                            ENDIF
                            1447 FORMAT(I6,',',2(I6,','),9(G13.6,','),G13.6)
                        ENDIF
                    ENDDO
                ENDDO
            ENDIF
        ENDIF
        !    *****************experimental***********************
        IF (DO_SEDIMENT_TRANSPORT) THEN
            L7120: DO  J=1,MYY
                M7120: DO  I=1,MX
                    IF (IS_SEDIMENT_COVERED(I,J)) THEN
                        IF (ELEVATION(I,J) <= SEDIMENT_BASE(I,J )) THEN
                            ICASE(17)=ICASE(17)-1
                        ENDIF
                        IF (ERODE_SLOPE(I,J) > 0.0 ) THEN
                            IF ( D8_GRADIENT(I,J) > (1.4*EQUILIBRIUM_GRADIENT(I,J)) ) THEN
                                DO_CONVERT=.FALSE.
                                L7122: DO  II=I-1,I+1
                                    III=II
                                    IF (IS_X_PERIODIC) THEN
                                        IF (III > MX) III=1
                                        IF (III < 1) III=MX
                                    ELSE
                                        IF (III > MX) CYCLE
                                        IF (III < 1) CYCLE
                                    ENDIF
                                    L7123: DO JJ=J-1,J+1
                                        JJJ=JJ
                                        IF (IS_Y_PERIODIC) THEN
                                            IF (JJJ > MY) JJJ=1
                                            IF (JJJ < 1) JJJ=MY
                                        ELSE
                                            IF (JJJ >= MY) CYCLE
                                            IF (JJJ < 1) CYCLE
                                        ENDIF
                                        IF (.NOT.IS_SEDIMENT_COVERED(III,JJJ)) THEN
                                            IF (ELEVATION(III,JJJ) > ELEVATION(I,J)) DO_CONVERT=.TRUE.
                                        ENDIF
                                    ENDDO L7123
                                ENDDO L7122
                                IF (DO_CONVERT) THEN
                                    CHANNEL_STATE_CHANGE(I,J)=2
                                    ICASE(18)=ICASE(18)-1
                                ENDIF
                            ENDIF
                        ENDIF
                    ENDIF
                ENDDO M7120
            ENDDO L7120
        ENDIF
        !      *********************************************************************
        !       This is where we actually change designations of locations to or from
        !         bedrock or alluvial
        !       matrix channel_state_change(i,j) is used to indicate whether a given matrix point is to be
        !         converted to a different type from its current type
        !         - a value of 0 means that there is no change
        !         - a value of 1 means that it is a bedrock channel converted to sediment
        !              covered (alluvial)
        !         - a value of 2 means that it is an alluvial channel converted to bedrock
        !      *********************************************************************
        NNONTOAL=0
        NALTONON=0
        L465: DO  J=1,MYY
            M465: DO  I=1,MX
                IF (CHANNEL_STATE_CHANGE(I,J) /= 0) THEN
                    IF (CHANNEL_STATE_CHANGE(I,J) == 1) THEN
                        NNONTOAL=NNONTOAL+1
                        IS_SEDIMENT_COVERED(I,J)=.TRUE.
                        IF (WRITEDETAIL) THEN
                            WRITE(OUTHIST,5440)
                            5440 FORMAT(' TO REG 3, CHANNEL_STATE_CHANGE=1')
                        ENDIF
                        IS_ROCK_SURFACE(I,J)=.FALSE.
                        REGOLITH(I,J)=DMAX1(REGOLITH(I,J),0.0D0)
                        ACCELERATED_EROSION(I,J)=.FALSE.
                        IF (SUBMERGED(I,J)) NSDELT=NSDELT+1
                    ELSE
                        IF (WRITEDETAIL) THEN
                            WRITE(OUTHIST,5439)
                            5439  FORMAT(' TO REG 4, SUBMERGED')
                        ENDIF
                        IS_SEDIMENT_COVERED(I,J)=.FALSE.
                        IS_ROCK_SURFACE(I,J)=.FALSE.
                        REGOLITH(I,J)=DMAX1(REGOLITH(I,J),0.0D0)
                        NALTONON=NALTONON+1
                    ENDIF
                ENDIF
                IF ((ELEVATION(I,J) < OCEAN_ELEVATION).AND. &
                (D8_GRADIENT(I,J) > DELTA_FORESET_GRADIENT)  &
                .AND.(ELEVATION(I,J) <= SEDIMENT_BASE(I,J))) THEN
                    IS_SEDIMENT_COVERED(I,J)=.FALSE.
                    IS_ROCK_SURFACE(I,J)=.TRUE.
                    REGOLITH(I,J)=0.0
                ENDIF
            ENDDO M465
        ENDDO L465
        IF (DO_SEDIMENT_TRANSPORT) THEN
            IF (WRITEDETAIL) THEN
                WRITE(OUTHIST,4734) NNONTOAL,NALTONON
                4734     FORMAT(' TO AL=',I13,' TO BED=',I13)
                L14734: DO  I=1,20
                    WRITE(OUTHIST,14735) I,ICASE(I)
                    14735       FORMAT(' I=',I3,' N=',I13)
                ENDDO L14734
                IF (SQYNUM > 0) THEN
                    SYSUM=SYSUM/SQYNUM
                    SQSUM=SQSUM/SQYNUM
                    SWSUM=SWSUM/SQYNUM
                    SGSUM=SGSUM/SQYNUM
                    SGREQ=SGREQ/SQYNUM
                    SSUMCH=SSUMCH/SQYNUM
                    WRITE(OUTHIST,14736) SQSUM,SYSUM,SWSUM,SGSUM,SGREQ,SSUMCH
                    14736       FORMAT(' CASE 13 AVERAGES,SEDIMENT_FLUX=',G12.5,' SEDIMENT_YIELD=',G12.5,&
                    ' Q=',G12.5,' GRAD=',G12.5, &
                    ' EQUILIBRIUM_GRADIENT=',G12.5,' DEL-E=',G12.5)
                ENDIF !    sqynum
                IF (XQYNUM > 0) THEN
                    XYSUM=XYSUM/XQYNUM
                    XQSUM=XQSUM/XQYNUM
                    XWSUM=XWSUM/XQYNUM
                    XGSUM=XGSUM/XQYNUM
                    XGREQ=XGREQ/XQYNUM
                    XSUMCH=XSUMCH/XQYNUM
                    WRITE(OUTHIST,14737) XQSUM,XYSUM,XWSUM,XGSUM,XGREQ,XSUMCH
                    14737       FORMAT(' CASE 16 AVERAGES, SEDIMENT_FLUX=',G12.5,' SEDIMENT_YIELD=',G12.5, &
                   ' Q=',G12.5,' GRAD=',G12.5, &
                    ' EQUILIBRIUM_GRADIENT=',G12.5,' DEL-E=',G12.5)
                ENDIF
            ENDIF
        ENDIF
        IF (BISTABLE_FLUVIAL_EROSION) THEN
            L466: DO  J=1,MYY
                M466: DO  I=1,MX
                    ECOMP=-((1.0-EROSION_RATE_CHANGE_LAG)*NET_EROSION_RATE(I,J)+EROSION_RATE_CHANGE_LAG*PREVIOUS_ELEVATION(I,J))
                    IF (ECOMP > HIGH_EROSION_THRESHOLD) THEN
                        IF (.NOT.IS_SEDIMENT_COVERED(I,J)) THEN
                            ACCELERATED_EROSION(I,J)=.TRUE.
                            DO_ACCELERATED_EROSION(I,J)=.TRUE.
                        ENDIF
                    ENDIF
                    IF (ECOMP < LOW_EROSION_THRESHOLD) ACCELERATED_EROSION(I,J)=.FALSE.
                    PREVIOUS_ELEVATION(I,J)=-ECOMP
                ENDDO M466
            ENDDO L466
        ENDIF
        IF (RESISTANT_SURFACE_LAYER) THEN
            L1465: DO  J=1,MYY
                M1465: DO  I=1,MX
                    IF (DO_ACCELERATED_EROSION(I,J)) THEN
                        RELATIVE_RESISTANCE(I,J)=1.0
                    ELSE
                        IF (.NOT.IS_SEDIMENT_COVERED(I,J)) THEN
                            IF ((ELEVATION(I,J) > PREVIOUS_ELEVATION(I,J)).AND.(ELEVATION(I,J) < INITIAL_ELEVATION(I,J))) THEN
                                RELATIVE_RESISTANCE(I,J)=SURFACE_LAYER_RESISTANCE
                            ELSE
                                RELATIVE_RESISTANCE(I,J)=1.0
                                DO_ACCELERATED_EROSION(I,J)=.TRUE.
                            ENDIF
                        ELSE
                            RELATIVE_RESISTANCE(I,J)=1.0
                        ENDIF
                    ENDIF
                ENDDO M1465
            ENDDO L1465
        ENDIF
        IP=64
        JP=243
        IF (DO_ACCELERATED_EROSION(IP,JP)) THEN
            IBT=1
        ELSE
            IBT=0
        ENDIF
        IF (WRITEDETAIL) THEN
            WRITE(OUTHIST,8223) IBT,RELATIVE_RESISTANCE(IP,JP),ELEVATION(IP,JP),PREVIOUS_ELEVATION(IP,JP), &
            ERODE_CHANNEL(IP,JP),ERODE_REGOLITH_CHANNEL(IP,JP),NET_EROSION_RATE(IP,JP), &
            CHANNEL_EROSION_RATE(IP,JP),D8_GRADIENT(IP,JP),ERODE_SLOPE(IP,JP),DISCHARGE(IP,JP)
            8223 FORMAT ('IBT=',I2,' RES=',G12.5,' ELEVATION=',G12.5,' PREVIOUS_ELEVATION=',G12.5,/ &
            ,' ECH=',G12.5,' EREG=',G12.5,' NET_EROSION_RATE=',G12.5,' CFBE=',G12.5,/ &
            ,' GRAD=',G12.5,' ERODE_SLOPE=',G12.5,' DISCHARGE=',G12.5)
        ENDIF
        RETURN
    END
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE WRITE_DEBUG(II,JJ,ESL,ECC,ECR,CHANCHANGE,SUMCHANGE)
        USE ERODE_GLOBALS
        IMPLICIT NONE
        INTEGER(4) :: I,J,II,JJ
        LOGICAL(4) :: DOXPROFILE,DOYPROFILE
        REAL(8) :: ESL,ECC,ECR,CHANCHANGE,SUMCHANGE
        REAL(8) :: ISUB,IROCK,ISED,ROUTEDIV,IIII
        DOXPROFILE=.FALSE.
        DOYPROFILE=.FALSE.
        IF (MX < 10) DOYPROFILE=.TRUE.
        IF (MY < 10) DOXPROFILE=.TRUE.
        IF (DOXPROFILE) THEN
            IF(JJ /= MY/2) RETURN
            I=II
            !                do 1100 i=1,mx
            IIII=DFLOAT(I)
            IF (SUBMERGED(I,JJ)) THEN
                ISUB=1.0
            ELSE
                ISUB=0.0
            ENDIF
            IF (IS_ROCK_SURFACE(I,JJ)) THEN
                IROCK=1.0
            ELSE
                IROCK=0.0
            ENDIF
            IF (CFN(I,JJ) < 1.0E+24) THEN
                ROUTEDIV=CFN(I,JJ)
            ELSE
                ROUTEDIV=0.0
            ENDIF
            IF (IS_SEDIMENT_COVERED(I,JJ)) THEN
                ISED=1.0
            ELSE
                ISED=0.0
            ENDIF
            WRITE(82,1130) IIII,ELEVATION(I,JJ),INITIAL_ELEVATION(I,JJ),SEDIMENT_BASE(I,JJ) &
            ,D8_GRADIENT(I,JJ),DISCHARGE(I,JJ), &
            ISUB,IROCK,ISED,   &
            ESL,ECC,ECR,CHANCHANGE,SUMCHANGE,CFNE(I,JJ),CFNW(I,JJ) &
            ,ROUTEDIV,CFW(I,JJ),REGOLITH(I,JJ) &
            ,EQUILIBRIUM_GRADIENT(I,JJ),SEDIMENT_FLUX(I,JJ)
            1130           FORMAT(21(G12.5,' '))
        ELSE
            IF (DOYPROFILE) THEN
                IF (II /= MX/2) RETURN
                J=JJ
                IIII=DFLOAT(J)
                IF (SUBMERGED(II,J)) THEN
                    ISUB=1.0
                ELSE
                    ISUB=0.0
                ENDIF
                IF (IS_ROCK_SURFACE(II,J)) THEN
                    IROCK=1.0
                ELSE
                    IROCK=0.0
                ENDIF
                IF (CFN(II,J) < 1.0E+24) THEN
                    ROUTEDIV=CFN(II,J)
                ELSE
                    ROUTEDIV=0.0
                ENDIF
                IF (IS_SEDIMENT_COVERED(II,J)) THEN
                    ISED=1.0
                ELSE
                    ISED=0.0
                ENDIF
                WRITE(82,1130) IIII,ELEVATION(II,J),INITIAL_ELEVATION(II,J),SEDIMENT_BASE(II,J) &
                ,D8_GRADIENT(II,J),DISCHARGE(II,J),  &
                ISUB,IROCK,ISED, &
                ESL,ECC,ECR,CHANCHANGE,SUMCHANGE,CFNE(II,J),CFNW(II,J) &
                ,ROUTEDIV,CFW(II,J),REGOLITH(II,J) &
                ,EQUILIBRIUM_GRADIENT(II,J),SEDIMENT_FLUX(II,J)
            ENDIF
        ENDIF
        RETURN
    END

    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE PRINT_AROUND(II,JJ)
        USE ERODE_GLOBALS, CHANNEL_STATE_CHANGE=>IDO, CHANNEL_EROSION_RATE=>CFNE,  &
        NET_EROSION_RATE=>CFNW,SED_FLUX_DIVERGENCE=>CFW, SED_SURFACE_ELEVATION=>CFN
        IMPLICIT NONE
        INTEGER(4) :: I,J,IW,IE,JN,JS,II,JJ
        REAL(8) :: TEMP1
        INTEGER(4) :: ITEMP1,ITEMP2,ITEMP3
        IW=MAX(1,II-1)
        IE=MIN(MX,II+1)
        JN=MAX(1,JJ-1)
        JS=MIN(MYY,JJ+1)
        WRITE(OUTHIST,100) II,JJ,TIME_INCREMENT
        100   FORMAT(' +SED AT I=',I5,' J=',I5,' TI=',G12.5)
        DO J=JN,JS
            DO I=IW,IE
                IF (IS_SEDIMENT_COVERED(I,J)) THEN
                    ITEMP1=1
                ELSE
                    ITEMP1=0
                ENDIF
                IF (IS_ROCK_SURFACE(I,J)) THEN
                    ITEMP2=1
                ELSE
                    ITEMP2=0
                ENDIF
                IF (SEDIMENT_DEPOSITED(I,J)) THEN
                    ITEMP3=1
                ELSE
                    ITEMP3=0
                ENDIF
                TEMP1=SED_SURFACE_ELEVATION(I,J)-ELEVATION(I,J)
                WRITE(OUTHIST,110) I,J,CHANNEL_STATE_CHANGE(I,J),ITEMP1,ITEMP2,ITEMP3, &
                SED_SURFACE_ELEVATION(I,J),ELEVATION(I,J),TEMP1,SED_FLUX_DIVERGENCE(I,J),ERODE_SLOPE(I,J), &
                ERODE_CHANNEL(I,J),ERODE_REGOLITH_CHANNEL(I,J),CHANNEL_EROSION_RATE(I,J),REGOLITH(I,J), &
                D8_GRADIENT(I,J),CHANNEL_WIDTH(I,J)
                110   FORMAT(    'I=',I5,' J=',I5,' CHANNEL_STATE_CHANGE=',I5,&
                ' SC=',I5,' RK=',I5,' SB=',I5,/, &
                '         SED_SURFACE_ELEVATION=',G12.5,' ELEVATION=',G12.5,' DELS=',G12.5, &
                ' SED_FLUX_DIVERGENCE=',G12.5  &
                ,/,'        ESL=',G12.5,' ECC=',G12.5,' ECR=',G12.5,' CHN ER RATE=',G12.5,' REG=',G12.5,&
                ' GRAD=',G12.5,' CWW=',G12.5)
            ENDDO
        ENDDO
        RETURN

    END
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE LOCAL_VALUES(I,J,LOCAL_INVERSE_RESISTANCE,NOMINAL_CRITICAL_SHEAR, &
        LOCAL_CRITICAL_SHEAR, LOCAL_BEDROCK_ERODIBILITY, LOCAL_REGOLITH_CRITICAL_SHEAR, &
        LOCAL_REGOLITH_ERODIBILITY)
        USE ERODE_GLOBALS
        IMPLICIT NONE
        REAL(8), INTENT(OUT)  :: LOCAL_INVERSE_RESISTANCE,LOCAL_CRITICAL_SHEAR
        REAL(8) , INTENT(OUT) :: LOCAL_REGOLITH_CRITICAL_SHEAR,LOCAL_REGOLITH_ERODIBILITY
        REAL(8), INTENT(IN) :: NOMINAL_CRITICAL_SHEAR
        REAL(8), INTENT(INOUT) :: LOCAL_BEDROCK_ERODIBILITY
        INTEGER(4), INTENT(IN) :: I,J
        LOCAL_INVERSE_RESISTANCE = 1.0/RELATIVE_RESISTANCE(I,J)
        LOCAL_CRITICAL_SHEAR = NOMINAL_CRITICAL_SHEAR &
        / LOCAL_INVERSE_RESISTANCE
        IF (VARIABLE_VEGETATION_RESISTANCE) THEN
            LOCAL_REGOLITH_ERODIBILITY=BEDROCK_ERODIBILITY
            IF(DRAINAGE_AREA(I,J) <= VEGETATION_AREA_MINIMUM) THEN
                LOCAL_REGOLITH_CRITICAL_SHEAR=VEGETATION_UPLAND_RESISTANCE
            ELSE
                IF (DRAINAGE_AREA(I,J) > VEGETATION_AREA_MAXIMUM) THEN
                    LOCAL_REGOLITH_CRITICAL_SHEAR=VEGETATION_CHANNEL_RESISTANCE
                ELSE
                    LOCAL_REGOLITH_CRITICAL_SHEAR= &
                    DEXP(DLOG(DRAINAGE_AREA(I,J))*VEGETATION_FACTOR_SLOPE+VEGETATION_FACTOR_INTERCEPT)
                ENDIF
            ENDIF
            LOCAL_CRITICAL_SHEAR=LOCAL_REGOLITH_CRITICAL_SHEAR*REGOLITH_CRITICAL_SHEAR_FACTOR
        ELSE
            IF (BISTABLE_FLUVIAL_EROSION) THEN
                IF (ACCELERATED_EROSION(I,J)) THEN
                    IF (USE_BISTABLE_BEDROCK) THEN
                        LOCAL_CRITICAL_SHEAR=BISTABLE_CRITICAL_SHEAR*LOCAL_CRITICAL_SHEAR/DETACHMENT_CRITICAL_SHEAR
                        LOCAL_BEDROCK_ERODIBILITY=BISTABLE_BEDROCK_ERODIBILITY
                        LOCAL_REGOLITH_CRITICAL_SHEAR=BISTABLE_CRITICAL_SHEAR/REGOLITH_CRITICAL_SHEAR_FACTOR
                        LOCAL_REGOLITH_ERODIBILITY=BISTABLE_BEDROCK_ERODIBILITY
                    ENDIF
                ELSE
                    LOCAL_REGOLITH_ERODIBILITY=BEDROCK_ERODIBILITY
                    LOCAL_REGOLITH_CRITICAL_SHEAR=DETACHMENT_CRITICAL_SHEAR/REGOLITH_CRITICAL_SHEAR_FACTOR
                ENDIF
            ELSE
                LOCAL_REGOLITH_ERODIBILITY=BEDROCK_ERODIBILITY
                LOCAL_REGOLITH_CRITICAL_SHEAR=DETACHMENT_CRITICAL_SHEAR/REGOLITH_CRITICAL_SHEAR_FACTOR
            ENDIF
        ENDIF
        RETURN
    END
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


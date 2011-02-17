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
    !     **********************************************************************
    !       these subroutines just read and write data and initialize variables
    !     **********************************************************************

    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE READ_ALLUVIAL_LOCATIONS(IUNIT)
        USE ERODE_GLOBALS
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: IUNIT
        INTEGER :: I,J,IIII(IMMX),MX1,MY1,ISTATUS
        !      *********************************************************************
        !       Reads matrix of integers indicating whether given matrix locations have an
        !         alluvial cover (1) or not (0) and sets up the logical is_sediment_covered matrix
        !      *********************************************************************
        L200: DO
            READ(IUNIT,130,IOSTAT=ISTATUS) MX1,MY1
            IF (ISTATUS < 0) RETURN
            130 FORMAT(2I5)
            IF ((MX1 /= MX).OR.(MY1 /= MY)) THEN
                WRITE(*,1121) MX, MX1,MY,MY1
                1121      FORMAT(' INCOMPATIBLE INPUT DIMENSIONS, MX=',I5,' MX1=',I5, &
                '  MY=',I5,' MY1=',I5)
                STOP
            ENDIF
            DO I=1,MX
                READ(IUNIT,110)(IIII(J),J=1,MY)
                110 FORMAT(256I1)
                DO J=1,MY
                    IF (IIII(J) == 1) THEN
                        IS_SEDIMENT_COVERED(I,J)=.TRUE.
                    ELSE
                        IS_SEDIMENT_COVERED(I,J)=.FALSE.
                    ENDIF
                ENDDO
            ENDDO
        ENDDO L200
        RETURN
    END
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE READ_BISTABLE_LOCATIONS(IUNIT)
        USE ERODE_GLOBALS
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: IUNIT
        INTEGER :: I,J,IIII(IMMX),ISTATUS
        !      *********************************************************************
        !       Reads matrix of integers indicating whether given matrix locations have an
        !         alluvial cover (1) or not (0) and sets up the logical is_sediment_covered matrix
        !      *********************************************************************
        L200: DO
            READ(IUNIT,130,IOSTAT=ISTATUS) MX,MY
            IF (ISTATUS < 0) RETURN
            130 FORMAT(2I5)
            DO  I=1,MX
                READ(IUNIT,110)(IIII(J),J=1,MY)
                110 FORMAT(256I1)
                DO  J=1,MY
                    IF (IIII(J) == 1) THEN
                        ACCELERATED_EROSION(I,J)=.TRUE.
                        DO_ACCELERATED_EROSION(I,J)=.TRUE.
                    ELSE
                        ACCELERATED_EROSION(I,J)=.FALSE.
                        DO_ACCELERATED_EROSION(I,J)=.FALSE.
                    ENDIF
                ENDDO
            ENDDO
        ENDDO L200
        RETURN
    END
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE READ_BEDROCK_LOCATIONS(IUNIT)
        USE ERODE_GLOBALS
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: IUNIT
        INTEGER :: I,J,IIII(IMMX),MX1,MY1,ISTATUS
        !      *********************************************************************
        !       Reads matrix of integers indicating whether given matrix locations
        !       are exposed bedrock (1) or regolith (0) covered
        !      *********************************************************************
        L200: DO
            READ(IUNIT,130,IOSTAT=ISTATUS) MX1,MY1
            IF (ISTATUS < 0) RETURN
            130 FORMAT(2I5)
            IF ((MX1 /= MX).OR.(MY1 /= MY)) THEN
                WRITE(*,1121) MX, MX1,MY,MY1
                1121      FORMAT(' INCOMPATIBLE INPUT DIMENSIONS, MX=',I5,' MX1=',I5, &
                '  MY=',I5,' MY1=',I5)
                STOP
            ENDIF
            DO  I=1,MX
                READ(IUNIT,110)(IIII(J),J=1,MY)
                110 FORMAT(256I1)
                DO  J=1,MY
                    IF (IIII(J) == 1) THEN
                        IS_ROCK_SURFACE(I,J)=.TRUE.
                    ELSE
                        IS_ROCK_SURFACE(I,J)=.FALSE.
                    ENDIF
                ENDDO
            ENDDO
        ENDDO L200
        RETURN
    END
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE READ_SEDIMENT_BASE(IUNIT,EEEMIN)
        !      *********************************************************************
        !        Reads matrix of elevations of the base of the sedimentary deposits -
        !         generally only used if restarting
        !      *********************************************************************
        USE ERODE_GLOBALS
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: IUNIT
        INTEGER :: I,J,ISTATUS
        REAL (8), INTENT(IN) :: EEEMIN
        L200: DO
            READ(IUNIT,*,IOSTAT=ISTATUS) MX,MY
            IF (ISTATUS < 0) RETURN
            DO  I = 1, MX
                DO  J = 1, MY
                    READ(IUNIT,*) SEDIMENT_BASE(I,J)
                    IF (.NOT.NEW_SIMULATION) THEN
                        SEDIMENT_BASE(I,J)=SEDIMENT_BASE(I,J)+EEEMIN
                    ENDIF
                ENDDO
            ENDDO
            READ(IUNIT,*) MX,MY
            DO  I = 1, MX
                DO  J = 1, MY
                    READ(IUNIT,*) REGOLITH(I,J)
                ENDDO
            ENDDO
        ENDDO L200
        RETURN
    END
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE READ_REGOLITH_THICKNESS(IUNIT)
        !      *********************************************************************
        !       Reads matrix of regolith thickness
        !      *********************************************************************
        USE ERODE_GLOBALS
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: IUNIT
        INTEGER :: I,J,ISTATUS
        L200: DO
            READ(IUNIT,*,IOSTAT=ISTATUS) MX,MY
            IF (ISTATUS < 0) RETURN
            DO  I = 1, MX
                DO  J = 1, MY
                    READ(IUNIT,*) REGOLITH(I,J)
                    IF (NEW_SIMULATION) THEN
                        IF (REGOLITH(I,J) <= 0.0) THEN
                            REGOLITH(I,J)=0.0
                            IS_ROCK_SURFACE(I,J)=.TRUE.
                        ELSE
                            IS_ROCK_SURFACE(I,J)=.FALSE.
                        ENDIF
                    ENDIF
                ENDDO
            ENDDO
        ENDDO L200
        RETURN
    END
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE READ_DEFORMATION(IUNIT)
        !      *********************************************************************
        !        Reads matrix of cumulative surface and rock deformation
        !         generally only used if restarting
        !      *********************************************************************
        USE ERODE_GLOBALS
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: IUNIT
        INTEGER :: I,J,MX1,MY1
        READ(IUNIT,*) MX1,MY1
        IF ((MX1 /= MX).OR.(MY1 /= MY)) THEN
            WRITE(*,1121) MX, MX1,MY,MY1
            1121      FORMAT(' INCOMPATIBLE INPUT DIMENSIONS, MX=',I5,' MX1=',I5, &
            '  MY=',I5,' MY1=',I5)
            STOP
        ENDIF
        DO  I = 1, MX1
            DO  J = 1, MY1
                READ(IUNIT,*) DEFORMATION(I,J)
            ENDDO
        ENDDO
        RETURN
    END
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE WRITE_DEBUG_INFO()
        USE ERODE_GLOBALS
        IMPLICIT NONE
        WRITE(OUTHIST,433)
        433         FORMAT(/,' ELEVATION')
        CALL PRINT_REAL_MATRIX_DATA(ELEVATION)
        WRITE(OUTHIST,1433)
        1433          FORMAT(/,' AREA')
        CALL PRINT_REAL_MATRIX_DATA(DRAINAGE_AREA)
        WRITE(OUTHIST,1434)
        1434          FORMAT(/,' SURFACE DISCHARGE')
        CALL PRINT_REAL_MATRIX_DATA(DISCHARGE)
        WRITE(OUTHIST,434)
        434         FORMAT(/,' ECHANNEL')
        CALL PRINT_REAL_MATRIX_DATA(ERODE_CHANNEL)
        WRITE(OUTHIST,435)
        435         FORMAT(/,' SEDIMENT DIVERGENCE')
        CALL PRINT_REAL_MATRIX_DATA(CFW)
        WRITE(OUTHIST,436)
        436         FORMAT(/,' SEDIMENT YIELD')
        CALL PRINT_REAL_MATRIX_DATA(SEDIMENT_YIELD)
        WRITE(OUTHIST,437)
        437         FORMAT(/,' ERODE_SLOPE')
        CALL PRINT_REAL_MATRIX_DATA(ERODE_SLOPE)
        WRITE(OUTHIST,438)
        438         FORMAT(/,' GRADIENT')
        CALL PRINT_REAL_MATRIX_DATA(D8_GRADIENT)
        WRITE(OUTHIST,439)
        439         FORMAT(' SEDBASE')
        CALL PRINT_REAL_MATRIX_DATA(SEDIMENT_BASE)
        WRITE(OUTHIST,440)
        440         FORMAT(' SEDIMENT_FLUX')
        CALL PRINT_REAL_MATRIX_DATA(SEDIMENT_FLUX)
        WRITE(OUTHIST,441)
        441         FORMAT(' EQUILIBRIUM_GRADIENT')
        CALL PRINT_REAL_MATRIX_DATA(EQUILIBRIUM_GRADIENT)
        WRITE(OUTHIST,442)
        442         FORMAT(' EROSION RATE')
        CALL PRINT_REAL_MATRIX_DATA(CFNW)
        WRITE(OUTHIST,443)
        443         FORMAT(' CHANNEL EROSION')
        CALL PRINT_REAL_MATRIX_DATA(CFNE)
        IF (MODEL_GROUNDWATER) THEN
            WRITE(OUTHIST,1444)
            1444          FORMAT(/,' GROUNDWATER DIVERGENCE')
            CALL PRINT_REAL_MATRIX_DATA(GROUNDWATER_FLUX)
            WRITE(OUTHIST,1442)
            1442          FORMAT(/,' MAX DIVERGENCE')
            CALL PRINT_REAL_MATRIX_DATA(FILTERED_GROUNDWATER_FLUX)
        ENDIF
        WRITE(OUTHIST,1439)
        1439          FORMAT(/,' REGOLITH')
        CALL PRINT_REAL_MATRIX_DATA(REGOLITH)
        WRITE(OUTHIST,2434)
        2434         FORMAT(/,' BEDROCK')
        CALL PRINT_LOGICAL_MATRIX_DATA(IS_ROCK_SURFACE)
        WRITE(OUTHIST,2435)
        2435         FORMAT(/,' SEDIMENT COVER')
        CALL PRINT_LOGICAL_MATRIX_DATA(IS_SEDIMENT_COVERED)
        RETURN
    END
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE WRITE_GRADIENT_INFO()
        USE ERODE_GLOBALS
        IMPLICIT NONE
        !      *********************************************************************
        !      Writes summary of alluvial cover, actual, and alluvial channel
        !      gradients
        !      *********************************************************************
        RETURN
    END
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE WRITE_ALLUVIAL_LOCATIONS(IUNIT)
        USE ERODE_GLOBALS
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: IUNIT
        INTEGER :: I,J,IIII(JMMX)
        INTEGER (4) :: IXXX
        !      *********************************************************************
        !       Writes matrix of integers indicating whether given matrix locations have an
        !         alluvial cover (1) or not (0)
        !      *********************************************************************
        !      *********************************************************************
        !      alluvial dat is integer output file of alluvial (1) and non-alluvial matrix
        !         locations
        !      *********************************************************************
        OPEN(IUNIT,FILE='ALLUVIAL.DAT',POSITION='APPEND')
        WRITE(IUNIT,130) MX,MY
        130 FORMAT(2I5)
        DO  I=1,MX
            DO  J=1,MY
                IF (IS_SEDIMENT_COVERED(I,J)) THEN
                    IIII(J)=1
                ELSE
                    IIII(J)=0
                ENDIF
            ENDDO
            WRITE(IUNIT,110)(IIII(J),J=1,MY)
            110     FORMAT(256I1)
        ENDDO
        IXXX=IUNIT
        CLOSE(IUNIT)
        RETURN
    END
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE WRITE_BEDROCK_LOCATIONS(IUNIT)
        USE ERODE_GLOBALS
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: IUNIT
        INTEGER :: I,J,IIII(JMMX)
        INTEGER (4) :: IXXX
        !      *********************************************************************
        !       Writes matrix of integers indicating whether given matrix locations
        !         expose bedrock (1) or are regolith covered (0)
        !      *********************************************************************
        !     **********************************************************************
        !      bedrock.dat is integer output file of bedrock (1) and regolith (0)
        !         slopes as a matrix
        !     **********************************************************************
        OPEN(IUNIT,FILE='BEDROCK.DAT',POSITION='APPEND')
        WRITE(IUNIT,130) MX,MY
        130 FORMAT(2I5)
        DO  I=1,MX
            DO  J=1,MY
                IF (IS_ROCK_SURFACE(I,J)) THEN
                    IIII(J)=1
                ELSE
                    IIII(J)=0
                ENDIF
            ENDDO
            WRITE(IUNIT,110)(IIII(J),J=1,MY)
            110     FORMAT(256I1)
        ENDDO
        IXXX=IUNIT
        CLOSE(IUNIT)
        RETURN
    END
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE WRITE_LAKE_INFO(IUNIT)
        USE ERODE_GLOBALS
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: IUNIT
        INTEGER :: I,J,IIII(JMMX)
        INTEGER (4) :: IXXX
        !      *********************************************************************
        !       Writes matrix of integers indicating whether given matrix locations
        !         expose bedrock (1) or are regolith covered (0)
        !      *********************************************************************
        !     **********************************************************************
        !      bedrock.dat is integer output file of bedrock (1) and regolith (0)
        !         slopes as a matrix
        !     **********************************************************************
        OPEN(IUNIT,FILE='SUBMERGED.DAT',POSITION='APPEND')
        WRITE(IUNIT,130) MX,MY
        130 FORMAT(2I5)
        DO  I=1,MX
            DO  J=1,MY
                IF (SUBMERGED(I,J)) THEN
                    IIII(J)=1
                ELSE
                    IIII(J)=0
                ENDIF
            ENDDO
            WRITE(IUNIT,110)(IIII(J),J=1,MY)
            110     FORMAT(256I1)
        ENDDO
        IXXX=IUNIT
        CLOSE(IUNIT)
        RETURN
    END
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE WRITE_EROSION_DEPTH_INDEX(IUNIT)
        USE ERODE_GLOBALS
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: IUNIT
        INTEGER :: I,J
        INTEGER (4) :: IXXX
        !      *********************************************************************
        !       Writes matrix of integers indicating the k (vertical) location of the
        !         bedrock interface in the 3-d resistance cube - k goes vertically
        !         downward
        !      *********************************************************************
        !     **********************************************************************
        !      erosion_depth_index.dat is the k (vertical) index of the present location of the
        !         bedrock interface relative to the 3-d cube of bedrock resistance
        !     **********************************************************************
        OPEN(IUNIT,FILE='EROSION_DEPTH_INDEX.DAT',POSITION='APPEND')
        WRITE(IUNIT,130) MX,MY
        130 FORMAT(2I5)
        DO  I=1,MX
            DO  J=1,MY
                WRITE(IUNIT,110) EROSION_DEPTH_INDEX(I,J)
                110     FORMAT(I4)
            ENDDO
        ENDDO
        IXXX=IUNIT
        CLOSE(IUNIT)
        RETURN
    END
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE WRITE_ACCELERATED_EROSION_STATE(IUNIT)
        USE ERODE_GLOBALS
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: IUNIT
        INTEGER :: I,J, IIII(JMMX)
        INTEGER (4) :: IXXX
        !      *********************************************************************
        !       Writes matrix of integers indicating whether given matrix locations
        !         are have high erosion rates (1) or are normal (0) in bistable landscapes
        !      *********************************************************************
        OPEN(IUNIT,FILE='BISTABLE.DAT',POSITION='APPEND')
        WRITE(IUNIT,130) MX,MY
        130 FORMAT(2I5)
        DO  I=1,MX
            DO  J=1,MY
                IF (ACCELERATED_EROSION(I,J)) THEN
                    IIII(J)=1
                ELSE
                    IIII(J)=0
                ENDIF
            ENDDO
            WRITE(IUNIT,110)(IIII(J),J=1,MY)
            110     FORMAT(256I1)
        ENDDO
        IXXX=IUNIT
        CLOSE(IUNIT)
        RETURN
    END
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE WRITE_SEDIMENT_BASE()
        USE ERODE_GLOBALS
        IMPLICIT NONE
        INTEGER :: I,J
        REAL (8) :: TEMP1
        !     **********************************************************************
        !      Writes output file of the elevation of the base of alluvial sediment
        !     **********************************************************************
        OPEN(OUTBASE,FILE='OUTBASE.DAT',POSITION='APPEND')
        WRITE(OUTBASE,500) MX,MY
        500 FORMAT(2I5)
        DO  I = 1, MX
            DO  J = 1, MY
                TEMP1=SEDIMENT_BASE(I,J)/VERTICAL_SCALING
                WRITE(OUTBASE,501) TEMP1
                501     FORMAT(G13.6)
            ENDDO
        ENDDO
        CLOSE(OUTBASE)
        RETURN
    END
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE WRITE_REGOLITH_THICKNESS()
        USE ERODE_GLOBALS
        IMPLICIT NONE
        INTEGER :: I,J
        REAL (8) :: TEMP1
        !     **********************************************************************
        !      Writes output file of the regolith thickness
        !     **********************************************************************
        OPEN(OUTREGOLITH,FILE='REGOLITH.DAT',POSITION='APPEND')
        WRITE(OUTREGOLITH,500) MX,MY
        500 FORMAT(2I5)
        DO  I = 1, MX
            DO  J = 1, MY
                TEMP1 = REGOLITH(I,J)/VERTICAL_SCALING
                WRITE(OUTREGOLITH,501) TEMP1
            ENDDO
        ENDDO
        501     FORMAT(G13.6)
        CLOSE(OUTREGOLITH)
        RETURN
    END
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE WRITE_DEFORMATION(IOUNIT)
        USE ERODE_GLOBALS
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: IOUNIT
        INTEGER :: I,J
        REAL (8) :: TEMP1
        !     **********************************************************************
        !      Writes output file of cumulative deformation of the surface and
        !         bedrock
        !     **********************************************************************
        !     ********************************************************************
        !      deform.dat is a file of cumulative surface (and rock) deformation
        !        for now, only open during and for final output
        !     ********************************************************************
        OPEN(IOUNIT,FILE='DEFORM.DAT',POSITION='APPEND')
        500 FORMAT(2I5)
        501     FORMAT(G13.6)
        WRITE(IOUNIT,500) MX,MY
        DO  I = 1, MX
            DO  J = 1, MY
                TEMP1 = DEFORMATION(I,J)/VERTICAL_SCALING
                WRITE(IOUNIT,501) TEMP1
            ENDDO
        ENDDO
        CLOSE(IOUNIT)
        RETURN
    END
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE WRITE_ROCK_RESISTANCE()
        USE ERODE_GLOBALS
        IMPLICIT NONE
        INTEGER :: I,J
        REAL (8) :: TEMP1
        !     **********************************************************************
        !      Writes the present value of rock resistance at the bedock interface
        !     **********************************************************************
        OPEN(OUTRESIST,FILE='RESIST.OUT',POSITION='APPEND')
        WRITE(OUTRESIST,500) MX,MY
        500 FORMAT(2I5)
        DO  I = 1, MX
            DO  J = 1, MY
                TEMP1 = RELATIVE_RESISTANCE(I,J)
                WRITE(OUTRESIST,501) TEMP1
                501     FORMAT(G13.6)
            ENDDO
        ENDDO
        CLOSE(OUTRESIST)
        RETURN
    END
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE WRITE_FIRST_DATA_MATRIX()
        USE ERODE_GLOBALS
        IMPLICIT NONE
        INTEGER :: I,J,K
        INTEGER (4) :: IXXX
        LOGICAL :: DOXPROFILE,DOYPROFILE
        REAL (8) :: TEMP1
        !     **********************************************************************
        !      Writes the output file of one or more matrices of info on drainage basin
        !         spatial variables
        !     **********************************************************************
        !      ********************************************************************
        !      outelev.dat is output file of final elevations
        !      ********************************************************************
        OPEN(OUTDATA,FILE='OUTELEV.DAT',POSITION='APPEND')
        DOXPROFILE=.FALSE.
        DOYPROFILE=.FALSE.
        IF (MX < 10) DOYPROFILE=.TRUE.
        IF (MY < 10) DOXPROFILE=.TRUE.
        L600: DO  K=1,6
            IF(WRITETYPE(K) > 0) THEN
                IF (DOXPROFILE.OR.DOYPROFILE) THEN
                ELSE
                    WRITE(OUTDATA,500) MX,MY
                ENDIF
                500         FORMAT(2I5)
                IF (DOXPROFILE) THEN
                    J=MY/2
                    DO  I=1,MX
                        IF (K == 1) THEN
                            IF (WRITE_ABSOLUTE_ELEVATION) THEN
                                TEMP1 = ELEVATION(I,J)/VERTICAL_SCALING
                            ELSE
                                TEMP1 = (ELEVATION(I,J) - ELEVATION(1,MY))/VERTICAL_SCALING
                            ENDIF
                        ENDIF
                        IF (K == 2) TEMP1 = (ELEVATION(I,J)-SEDIMENT_BASE(I,J))/VERTICAL_SCALING
                        IF (K == 3) TEMP1 = D8_GRADIENT(I,J)
                        IF (K == 4) TEMP1 = ERODE_SLOPE(I,J)
                        IF (K == 5) TEMP1 = CFNE(I,J)
                        IF (K == 6) TEMP1 = CFNW(I,J)
                        WRITE(OUTDATA,501) TEMP1
                    ENDDO
                    WRITE(OUTDATA,1101)
                    1101        FORMAT(' ')
                ELSE
                    IF (DOYPROFILE) THEN
                        I=MX/2
                        DO  J=1,MY
                            IF (K == 1) THEN
                                IF (WRITE_ABSOLUTE_ELEVATION) THEN
                                    TEMP1 = ELEVATION(I,J)/VERTICAL_SCALING
                                ELSE
                                    TEMP1 = (ELEVATION(I,J) - ELEVATION(1,MY))/VERTICAL_SCALING
                                ENDIF
                            ENDIF
                            IF (K == 2) TEMP1 = (ELEVATION(I,J)-SEDIMENT_BASE(I,J))/VERTICAL_SCALING
                            IF (K == 3) TEMP1 = D8_GRADIENT(I,J)
                            IF (K == 4) TEMP1 = ERODE_SLOPE(I,J)
                            IF (K == 5) TEMP1 = CFNE(I,J)
                            IF (K == 6) TEMP1 = CFNW(I,J)
                            WRITE(OUTDATA,501) TEMP1
                        ENDDO
                        WRITE(OUTDATA,1101)
                    ELSE
                        DO  I = 1, MX
                            DO  J = 1, MY
                                IF (K == 1) THEN
                                    IF (WRITE_ABSOLUTE_ELEVATION) THEN
                                        TEMP1 = ELEVATION(I,J)/VERTICAL_SCALING
                                    ELSE
                                        TEMP1 = (ELEVATION(I,J) - ELEVATION(1,MY))/VERTICAL_SCALING
                                    ENDIF
                                ENDIF
                                IF (K == 2) TEMP1 = (ELEVATION(I,J)-SEDIMENT_BASE(I,J))/VERTICAL_SCALING
                                IF (K == 3) TEMP1 = D8_GRADIENT(I,J)
                                IF (K == 4) TEMP1 = ERODE_SLOPE(I,J)
                                IF (K == 5) TEMP1 = CFNE(I,J)
                                IF (K == 6) TEMP1 = CFNW(I,J)
                                WRITE(OUTDATA,501) TEMP1
                                501             FORMAT(G13.6)
                            ENDDO
                        ENDDO
                    ENDIF
                ENDIF
            ENDIF
        ENDDO L600
        IXXX=OUTDATA
        CLOSE(OUTDATA)
        RETURN
    END
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE WRITE_SECOND_DATA_MATRIX()
        USE ERODE_GLOBALS
        IMPLICIT NONE
        INTEGER :: I,J,K
        INTEGER (4) :: IXXX
        REAL (8) :: TEMP1
        !     **********************************************************************
        !      Writes one or more matrices of output spatial valiables
        !     **********************************************************************
        !      ********************************************************************
        !      outelev.dat is output file of final elevations
        !      ********************************************************************
        OPEN(OUTDATA,FILE='OUTELEV.DAT',POSITION='APPEND')
        DO  K=7,11
            IF(WRITETYPE(K) > 0) THEN
                WRITE(OUTDATA,500) MX,MY
                500         FORMAT(2I5)
                DO  I = 1, MX
                    DO  J = 1, MY
                        IF (K == 7) TEMP1 = CFW(I,J)
                        IF (K == 8) TEMP1 = CFN(I,J)
                        IF (K == 9) TEMP1 = CFNE(I,J)
                        IF (K == 10)TEMP1=  DLOG(DRAINAGE_AREA(I,J))
                        IF (K == 11) THEN
                            IF ((D8_GRADIENT(I,J) > 0.0).AND.(DRAINAGE_AREA(I,J) > 0.0)) THEN
                                TEMP1 = DLOG(DRAINAGE_AREA(I,J)/D8_GRADIENT(I,J))
                            ELSE
                                TEMP1 = 0.0
                            ENDIF
                        ENDIF
                        IF (K == 12) TEMP1=REGOLITH(I,J)/VERTICAL_SCALING
                        IF (K == 13) TEMP1=RELATIVE_RESISTANCE(I,J)
                        WRITE(OUTDATA,501) TEMP1
                        501             FORMAT(G13.6)
                    ENDDO
                ENDDO
            ENDIF
        ENDDO
        IXXX=OUTDATA
        CLOSE(OUTDATA)
        RETURN
    END
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE WRITE_DATA_SAMPLE(I,J)
        USE ERODE_GLOBALS
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: I, J
        INTEGER :: ISUBM,ISEDC
        !      *********************************************************************
        !       Writes out debugging info for central debug point
        !      *********************************************************************
        IF (IS_SEDIMENT_COVERED(I,J)) THEN
            ISEDC=1
        ELSE
            ISEDC=0
        ENDIF
        IF (SUBMERGED(I,J)) THEN
            ISUBM=1
        ELSE
            ISUBM=0
        ENDIF
        WRITE(OUTHIST,100) I,J,DRAINAGE_AREA(I,J),FLOW_DIRECTION(I,J),IDO(I,J),ISEDC,ISUBM,ELEVATION(I, &
        J),SEDIMENT_BASE(I,J),D8_GRADIENT(I,J),ERODE_CHANNEL(I,J),ERODE_SLOPE(I,J), &
        CFW(I,J),CFN(I,J), &
        CFNE(I,J),CFNW(I,J),SEDIMENT_FLUX(I,J),EQUILIBRIUM_GRADIENT(I,J),SEDIMENT_YIELD(I,J), &
        MAXIMUM_ELEVATION_CHANGE
        100   FORMAT(' I=',I5,' J=',I5,' AR=',G10.3,' ID=',I5,'IDO=',I5, &
        ' ISEDC=',  &
        I5,' ISUBM=',I5,/,' ELEVATION=',G12.5,' SEDB=',G12.5,' GRAD=', &
        G12.5,' ECHAN=',G12.5,' ESL=',G12.5,/,' CFW=',G12.5, &
        ' CFN=',G12.5,' CFNE=',G12.5,' CFNW=',G12.5,/,' SQ=', &
        G12.5,' GREQ=',G12.5,' SY=',G12.5,' CHG=',G12.5)
        RETURN
    END
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE WRITE_REPORT()
        !      *********************************************************************
        !       This subroutine writes out periodic reports on erosion rate and
        !         relief:
        !              maximum relief
        !              average relative elevation
        !              average gradient
        !              average erosion rate
        !      *********************************************************************
        USE ERODE_GLOBALS
        IMPLICIT NONE
        INTEGER :: I, J
        REAL (8) :: EP, EM, EV,  AG,ER, SC
        REAL (8) :: BF,AF,HF,HTF
        INTEGER :: IXXX
        EM =-1.0E+20
        EP =1.0E+20
        EV = 0.0
        SC = 0.0
        AG = 0.0
        BF=0.0
        AF=0.0
        HF=0.0
        HTF=0.0
        DO  J=1,MYY
            DO  I=1,MX
                IF (ACCELERATED_EROSION(I,J)) HF=HF+1.0
                IF (DO_ACCELERATED_EROSION(I,J)) HTF=HTF+1.0
                IF (IS_ROCK_SURFACE(I,J)) BF=BF+1.0
                IF (IS_SEDIMENT_COVERED(I,J)) AF=AF+1.0
                IF (ELEVATION(I,J)  <  EP) EP = ELEVATION(I,J)
            ENDDO
        ENDDO
        DO  J=1,MYY
            DO  I=1,MX
                SC = SC + CFNW(I,J)
                ER = ELEVATION(I,J)-EP
                EV = EV + ER
                IF (ER  >  EM) EM = ER
                AG = AG + D8_GRADIENT(I,J)
            ENDDO
        ENDDO
        EV = EV/(MX*(MYY))
        AG = AG/((MX)*(MYY))
        SC = SC/(MX*(MYY))
        HF=HF/(MX*MYY)
        HTF=HTF/(MX*MYY)
        BF=BF/(MX*MYY)
        AF=AF/(MX*MYY)
        WRITE(OUTREPORT,200) ITERATION,PRESENT_TIME,EM,EV,AG,SC,HF,HTF,BF,AF
        200 FORMAT(' ',I5,9(' ',G12.5))
        IXXX=OUTREPORT
        RETURN
    END
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE WRITE_IMAGE()
        USE ERODE_GLOBALS
        IMPLICIT NONE
        INTEGER :: I,J
        INTEGER (1) :: BVAL
        CHARACTER :: ACHAR
        REAL (8) :: MAXELEV,MINELEV,AAAA,BBBB
        REAL (8) :: ZMAX,ZMIN
        INTEGER :: ITEMP,IFILE
        !     **********************************************************************
        !      Writes out a raw image file of elevations normalized to 256 byte range
        !         and an information file on the iteration and actual elevation range
        !     **********************************************************************
        !     **********************************************************************
        !      relele???.raw is a file of raw byte-normalized elevation data -- an x,y
        !         image of relative elevations normalized to 256 values, 0 is lowest
        !         and 255 is highest elevation
        !     **********************************************************************
        !     **********************************************************************
        !      topo.dat is a record of minimum elevation, time, and elevation range
        !         for each output record in image.raw
        !     **********************************************************************
        IFILE=57
        RFILENAME='RELELE'//RNUM3//RNUM2//RNUM1//RSUFFIX
        OPEN(IFILE,FILE=RFILENAME)
        OPEN(OUTIMGDAT,FILE='TOPO.DAT',POSITION='APPEND')
        MAXELEV=-1.0E+25
        MINELEV=-MAXELEV
        DO  I=1,MX
            DO  J=1,MY
                IF (ELEVATION(I,J) > MAXELEV) MAXELEV=ELEVATION(I,J)
                IF (ELEVATION(I,J) < MINELEV) MINELEV=ELEVATION(I,J)
            ENDDO
        ENDDO
        ZMAX=MAXELEV/VERTICAL_SCALING
        ZMIN=MINELEV/VERTICAL_SCALING
        AAAA=256.0/(MAXELEV-MINELEV)
        BBBB=-AAAA*MINELEV
        DO  J=1,MY
            DO  I=1,MX
                ITEMP=INT(AAAA*ELEVATION(I,J)+BBBB)
                IF (ITEMP < 0) ITEMP=0
                IF (ITEMP > 255) ITEMP=255
                BVAL=ITEMP
                ACHAR=CHAR(BVAL)
                WRITE(IFILE,130,ADVANCE='NO') ACHAR
                130   FORMAT(A1)
            ENDDO
        ENDDO
        WRITE(OUTIMGDAT,140) ITERATION,MX,MY,ZMAX,ZMIN,PRESENT_TIME
        140   FORMAT(3(' ',I6),3(' ',E15.7))
        CLOSE(OUTIMGDAT)
        CLOSE(57)
        RETURN
    END
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE WRITE_SHADED_RELIEF_IMAGE()
        USE ERODE_GLOBALS
        IMPLICIT NONE
        REAL (8) :: HALFPI,PI,RADIAN,NOMANGLE,DTXX
        REAL (8) :: ETEMP(IMMX,JMMX)
        REAL (8) :: YDIV,SLX,SLY,SLOPE,SLMAX,SLSUN,ANGLE,SLOPEANGLE,NEWVALUE
        INTEGER :: I,J,IW,IE,JS,JN,SLOPENUMB,IMX,JMX
        REAL (8) ::     NOMSLMAX,RRRR,RANGE,SLOPEAVG,XDIV,MAXVAL,MINVAL
        INTEGER :: II,JJ,IISTART,IIEND,JJSTART,JJEND,XREFLECT,YREFLECT
        REAL (8) :: XVAL1,XVAL2,XVAL3,XVAL4,XAVG
        INTEGER (1) :: BVAL
        INTEGER :: ITEMP
        CHARACTER :: ACHAR
        IMX=MX
        JMX=MY
        IF (IS_X_PERIODIC) THEN
            XREFLECT=1
        ELSE
            XREFLECT=0
        ENDIF
        IF (IS_Y_PERIODIC) THEN
            YREFLECT=1
        ELSE
            YREFLECT=0
        ENDIF
        !     **********************************************************************
        !     Makes a raw shaded relief raw image of the simulated terrain and
        !             writes it out to a file
        !     The file size depends on whether the lateral and/or top boundaries are
        !      periodic.  If both are false then the raw image size is (2*(MX-1),2*(MY-1))
        !     If the image is periodic in, say, the X dimension, then the horizontal image
        !       size increases to 2*(MX-1)+2*(MX/2+2) and similarly if the image is periodic
        !       in the Y dimension
        !     **********************************************************************
        HALFPI =1.5708
        PI=3.141593
        RADIAN=0.01745329
        NOMANGLE=25.0
        SLMAX=-1.0E+25
        SLOPEAVG=0.0
        SLOPENUMB= 0
        DO  J = 2,JMX-1
            DO  I = 2,IMX-1
                IW=I-1
                IE=I+1
                JS=J+1
                JN=J-1
                SLX = ( ELEVATION(IE,J)- ELEVATION(IW,J))/2.0
                SLY = ( ELEVATION(I,JN)- ELEVATION(I,JS))/2.0
                SLOPE =SQRT(SLX**2+SLY**2)
                IF (SLOPE > SLMAX) SLMAX = SLOPE
                SLOPENUMB = SLOPENUMB+1
                SLOPEAVG = SLOPEAVG+SLOPE
            ENDDO
        ENDDO
        NOMSLMAX = SIN(HALFPI-NOMANGLE*RADIAN)/  &
        COS(HALFPI-NOMANGLE*RADIAN)
        DTXX = SLMAX/NOMSLMAX
        SLMAX = NOMSLMAX
        SLOPEAVG = SLOPEAVG/(SLOPENUMB*DTXX)
        SLSUN=HALFPI-DATAN(SLMAX)
        DO  J = 1,JMX
            DO  I = 1,IMX
                IW=I-1
                IE=I+1
                JS=J+1
                JN=J-1
                YDIV=2.0*DTXX
                XDIV=2.0*DTXX
                IF (JS > JMX) THEN
                    YDIV=1.0*DTXX
                    JS=JMX
                ENDIF
                IF (JN < 1) THEN
                    YDIV=1.0*DTXX
                    JN=1
                ENDIF
                IF (IW < 1) THEN
                    XDIV=1.0*DTXX
                    IW=1
                ENDIF
                IF (IE > IMX) THEN
                    XDIV=1.0*DTXX
                    IE=IMX
                ENDIF
                SLX = ( ELEVATION(IE,J)- ELEVATION(IW,J))/XDIV
                SLY = ( ELEVATION(I,JS)- ELEVATION(I,JN))/YDIV
                SLOPE =SQRT(SLX**2+SLY**2)
                SLOPEANGLE= DATAN(SLOPE)
                IF ((SLY == 0.0).OR. (SLX == 0.0)) THEN
                    IF (SLY == 0.0) THEN
                        IF (SLX > 0.0) THEN
                            ANGLE =  HALFPI
                        ELSE
                            ANGLE = -HALFPI
                        ENDIF
                    ELSE
                        IF (SLY > 0.0) THEN
                            ANGLE = 0.0
                        ELSE
                            ANGLE = PI
                        ENDIF
                    ENDIF
                ELSE
                    ANGLE = DATAN ((SLX) / (SLY))
                    IF (((SLY)  <  0.0).AND.((SLX)  <  0.0)) THEN
                        ANGLE = ANGLE + PI
                    ELSE
                        IF ((SLY)  <  0.0) THEN
                            ANGLE = ANGLE + PI
                        ENDIF
                    ENDIF
                ENDIF
                ANGLE = ANGLE-HALFPI
                NEWVALUE =   &
                (COS(SLSUN)*COS(SLOPEANGLE)+SIN(SLSUN) &
                *SIN(SLOPEANGLE)*  &
                COS(ANGLE))
                IF (NEWVALUE < 0.0) NEWVALUE = 0.0
                ETEMP(I,J) = NEWVALUE
            ENDDO
        ENDDO
        MAXVAL = -1.0E+25
        MINVAL = -MAXVAL
        DO  J = 1,JMX
            DO  I = 1,IMX
                RRRR = ETEMP(I,J)
                IF (RRRR < MINVAL) MINVAL = RRRR
                IF (RRRR >  MAXVAL) MAXVAL = RRRR
            ENDDO
        ENDDO
        RANGE=MAXVAL-MINVAL
        RNUM1=CHAR(IFILE1)
        IF (IFILE2 > 57)  THEN
            RNUM2=CHAR(IFILE2+40)
        ELSE
            RNUM2=CHAR(IFILE2)
        ENDIF
        RNUM3=CHAR(IFILE3)
        IFILE1=IFILE1+1
        IF (IFILE1 > 57) THEN
            IFILE1=48
            IFILE2=IFILE2+1
            IF (IFILE2 > 57) THEN
                IFILE2=48
                IFILE3=IFILE3+1
            ENDIF
        ENDIF
        RFILENAME=RPREFIX//RNUM3//RNUM2//RNUM1//RSUFFIX
        OPEN(44,FILE=RFILENAME)
        IF (YREFLECT > 0) THEN
            JJSTART=JMX-JMX/4
            JJEND=JMX
            L1120: DO  J=JJSTART,JJEND
                L1121: DO  JJ=1,2
                    IF (XREFLECT > 0) THEN
                        IISTART=IMX-IMX/4
                        IIEND=IMX
                        L1123: DO  I=IISTART,IIEND
                            XVAL1=ETEMP(I,J)
                            IF (I == IIEND) THEN
                                XVAL2=ETEMP(1,J)
                                IF (J == JJEND) THEN
                                    XVAL4=ETEMP(1,1)
                                ELSE
                                    XVAL4=ETEMP(1,J+1)
                                ENDIF
                            ELSE
                                XVAL2=ETEMP(I+1,J)
                                IF (J == JJEND) THEN
                                    XVAL4=ETEMP(I+1,1)
                                ELSE
                                    XVAL4=ETEMP(I+1,J+1)
                                ENDIF
                            ENDIF
                            IF (J == JJEND) THEN
                                XVAL3=ETEMP(I,1)
                            ELSE
                                XVAL3=ETEMP(I,J+1)
                            ENDIF
                            L1124: DO  II=1,2
                                XAVG=XVAL1*WEIGHTS(II,JJ,1)+ &
                                XVAL2*WEIGHTS(II,JJ,2)+ &
                                XVAL3*WEIGHTS(II,JJ,3)+ &
                                XVAL4*WEIGHTS(II,JJ,4)
                                ITEMP=INT((XAVG-MINVAL)*253.0/RANGE)+1
                                IF (ITEMP < 0) ITEMP=0
                                IF (ITEMP > 255) ITEMP=255
                                BVAL=ITEMP
                                ACHAR=CHAR(BVAL)
                                WRITE(44,130,ADVANCE='NO') ACHAR
                            ENDDO L1124
                        ENDDO L1123
                    ENDIF
                    L1125: DO  I=1,IMX-1
                        XVAL1=ETEMP(I,J)
                        XVAL2=ETEMP(I+1,J)
                        IF (J == JJEND) THEN
                            XVAL4=ETEMP(I+1,1)
                            XVAL3=ETEMP(I,1)
                        ELSE
                            XVAL4=ETEMP(I+1,J+1)
                            XVAL3=ETEMP(I,J+1)
                        ENDIF
                        L1126: DO  II=1,2
                            XAVG=XVAL1*WEIGHTS(II,JJ,1)+  &
                            XVAL2*WEIGHTS(II,JJ,2)+ &
                            XVAL3*WEIGHTS(II,JJ,3)+ &
                            XVAL4*WEIGHTS(II,JJ,4)
                            ITEMP=INT((XAVG-MINVAL)*253.0/RANGE)+1
                            IF (ITEMP < 0) ITEMP=0
                            IF (ITEMP > 255) ITEMP=255
                            BVAL=ITEMP
                            ACHAR=CHAR(BVAL)
                            WRITE(44,130,ADVANCE='NO') ACHAR
                        ENDDO L1126
                    ENDDO L1125
                    IF (XREFLECT > 0) THEN
                        IISTART=0
                        IIEND=IMX/4
                        L1127: DO  I=IISTART,IIEND
                            XVAL2=ETEMP(I+1,J)
                            IF (I == IISTART) THEN
                                XVAL1=ETEMP(IMX,J)
                                IF (J == JJEND) THEN
                                    XVAL3=ETEMP(IMX,1)
                                ELSE
                                    XVAL3=ETEMP(IMX,J+1)
                                ENDIF
                            ELSE
                                XVAL1=ETEMP(I+1,J)
                                IF (J == JJEND) THEN
                                    XVAL3=ETEMP(I,1)
                                ELSE
                                    XVAL3=ETEMP(I,J+1)
                                ENDIF
                            ENDIF
                            IF (J == JJEND) THEN
                                XVAL4=ETEMP(I+1,1)
                            ELSE
                                XVAL4=ETEMP(I+1,J+1)
                            ENDIF
                            L1128: DO  II=1,2
                                XAVG=XVAL1*WEIGHTS(II,JJ,1)+ &
                                XVAL2*WEIGHTS(II,JJ,2)+  &
                                XVAL3*WEIGHTS(II,JJ,3)+ &
                                XVAL4*WEIGHTS(II,JJ,4)
                                ITEMP=INT((XAVG-MINVAL)*253.0/RANGE)+1
                                IF (ITEMP < 0) ITEMP=0
                                IF (ITEMP > 255) ITEMP=255
                                BVAL=ITEMP
                                ACHAR=CHAR(BVAL)
                                WRITE(44,130,ADVANCE='NO') ACHAR
                            ENDDO L1128
                        ENDDO L1127
                    ENDIF
                ENDDO L1121
            ENDDO L1120
        ENDIF
        L120: DO  J=1,JMX-1
            L121: DO  JJ=1,2
                IF (XREFLECT > 0) THEN
                    IISTART=IMX-IMX/4
                    IIEND=IMX
                    L123: DO  I=IISTART,IIEND
                        XVAL1=ETEMP(I,J)
                        IF (I == IIEND) THEN
                            XVAL2=ETEMP(1,J)
                            XVAL4=ETEMP(1,J+1)
                        ELSE
                            XVAL2=ETEMP(I+1,J)
                            XVAL4=ETEMP(I+1,J+1)
                        ENDIF
                        XVAL3=ETEMP(I,J+1)
                        L124: DO  II=1,2
                            XAVG=XVAL1*WEIGHTS(II,JJ,1)+  &
                            XVAL2*WEIGHTS(II,JJ,2)+ &
                            XVAL3*WEIGHTS(II,JJ,3)+ &
                            XVAL4*WEIGHTS(II,JJ,4)
                            ITEMP=INT((XAVG-MINVAL)*253.0/RANGE)+1
                            IF (ITEMP < 0) ITEMP=0
                            IF (ITEMP > 255) ITEMP=255
                            BVAL=ITEMP
                            ACHAR=CHAR(BVAL)
                            WRITE(44,130,ADVANCE='NO') ACHAR
                        ENDDO L124
                    ENDDO L123
                ENDIF
                L125: DO  I=1,IMX-1
                    XVAL1=ETEMP(I,J)
                    XVAL2=ETEMP(I+1,J)
                    XVAL4=ETEMP(I+1,J+1)
                    XVAL3=ETEMP(I,J+1)
                    L126: DO  II=1,2
                        XAVG=XVAL1*WEIGHTS(II,JJ,1)+ &
                        XVAL2*WEIGHTS(II,JJ,2)+ &
                        XVAL3*WEIGHTS(II,JJ,3)+ &
                        XVAL4*WEIGHTS(II,JJ,4)
                        ITEMP=INT((XAVG-MINVAL)*253.0/RANGE)+1
                        IF (ITEMP < 0) ITEMP=0
                        IF (ITEMP > 255) ITEMP=255
                        BVAL=ITEMP
                        ACHAR=CHAR(BVAL)
                        WRITE(44,130,ADVANCE='NO') ACHAR
                    ENDDO L126
                ENDDO L125
                IF (XREFLECT > 0) THEN
                    IISTART=0
                    IIEND=IMX/4
                    L127: DO  I=IISTART,IIEND
                        XVAL2=ETEMP(I+1,J)
                        IF (I == IISTART) THEN
                            XVAL1=ETEMP(IMX,J)
                            XVAL3=ETEMP(IMX,J+1)
                        ELSE
                            XVAL1=ETEMP(I+1,J)
                            XVAL3=ETEMP(I+1,J+1)
                        ENDIF
                        XVAL4=ETEMP(I+1,J+1)
                        L128: DO  II=1,2
                            XAVG=XVAL1*WEIGHTS(II,JJ,1)+ &
                            XVAL2*WEIGHTS(II,JJ,2)+ &
                            XVAL3*WEIGHTS(II,JJ,3)+ &
                            XVAL4*WEIGHTS(II,JJ,4)
                            ITEMP=INT((XAVG-MINVAL)*253.0/RANGE)+1
                            IF (ITEMP < 0) ITEMP=0
                            IF (ITEMP > 255) ITEMP=255
                            BVAL=ITEMP
                            ACHAR=CHAR(BVAL)
                            WRITE(44,130,ADVANCE='NO') ACHAR
                        ENDDO L128
                    ENDDO L127
                ENDIF
            ENDDO L121
        ENDDO L120
        IF (YREFLECT > 0) THEN
            JJSTART=0
            JJEND=JMX/4
            L2120: DO  J=JJSTART,JJEND
                L2121: DO  JJ=1,2
                    IF (XREFLECT > 0) THEN
                        IISTART=IMX-IMX/4
                        IIEND=IMX
                        L2123: DO  I=IISTART,IIEND
                            IF (J == JJSTART) THEN
                                XVAL1=ETEMP(I,JMX)
                            ELSE
                                XVAL1=ETEMP(I,J)
                            ENDIF
                            IF (I == IIEND) THEN
                                IF (J == JJSTART) THEN
                                    XVAL2=ETEMP(1,JMX)
                                ELSE
                                    XVAL2=ETEMP(1,J)
                                ENDIF
                                XVAL4=ETEMP(1,J+1)
                            ELSE
                                IF (J == JJSTART) THEN
                                    XVAL2=ETEMP(I+1,JMX)
                                ELSE
                                    XVAL2=ETEMP(I+1,J)
                                ENDIF
                                XVAL4=ETEMP(I+1,J+1)
                            ENDIF
                            XVAL3=ETEMP(I,J+1)
                            L2124: DO  II=1,2
                                XAVG=XVAL1*WEIGHTS(II,JJ,1)+  &
                                XVAL2*WEIGHTS(II,JJ,2)+ &
                                XVAL3*WEIGHTS(II,JJ,3)+ &
                                XVAL4*WEIGHTS(II,JJ,4)
                                ITEMP=INT((XAVG-MINVAL)*253.0/RANGE)+1
                                IF (ITEMP < 0) ITEMP=0
                                IF (ITEMP > 255) ITEMP=255
                                BVAL=ITEMP
                                ACHAR=CHAR(BVAL)
                                WRITE(44,130,ADVANCE='NO') ACHAR
                            ENDDO L2124
                        ENDDO L2123
                    ENDIF
                    L2125: DO  I=1,IMX-1
                        IF (J == JJSTART) THEN
                            XVAL1=ETEMP(I,JMX)
                            XVAL2=ETEMP(I+1,JMX)
                        ELSE
                            XVAL1=ETEMP(I,J)
                            XVAL2=ETEMP(I+1,J)
                        ENDIF
                        XVAL4=ETEMP(I+1,J+1)
                        XVAL3=ETEMP(I,J+1)
                        L2126: DO  II=1,2
                            XAVG=XVAL1*WEIGHTS(II,JJ,1)+ &
                            XVAL2*WEIGHTS(II,JJ,2)+ &
                            XVAL3*WEIGHTS(II,JJ,3)+ &
                            XVAL4*WEIGHTS(II,JJ,4)
                            ITEMP=INT((XAVG-MINVAL)*253.0/RANGE)+1
                            IF (ITEMP < 0) ITEMP=0
                            IF (ITEMP > 255) ITEMP=255
                            BVAL=ITEMP
                            ACHAR=CHAR(BVAL)
                            WRITE(44,130,ADVANCE='NO') ACHAR
                        ENDDO L2126
                    ENDDO L2125
                    IF (XREFLECT > 0) THEN
                        IISTART=0
                        IIEND=IMX/4
                        L2127: DO  I=IISTART,IIEND
                            IF (J == JJSTART) THEN
                                XVAL2=ETEMP(I+1,JMX)
                            ELSE
                                XVAL2=ETEMP(I+1,J)
                            ENDIF
                            IF (I == IISTART) THEN
                                IF (J == JJSTART) THEN
                                    XVAL1=ETEMP(IMX,JMX)
                                ELSE
                                    XVAL1=ETEMP(IMX,J)
                                ENDIF
                                XVAL3=ETEMP(IMX,J+1)
                            ELSE
                                IF (J == JJSTART) THEN
                                    XVAL1=ETEMP(I,JMX)
                                ELSE
                                    XVAL1=ETEMP(I,J)
                                ENDIF
                                XVAL3=ETEMP(I,J+1)
                            ENDIF
                            XVAL4=ETEMP(I+1,J+1)
                            L2128: DO  II=1,2
                                XAVG=XVAL1*WEIGHTS(II,JJ,1)+  &
                                XVAL2*WEIGHTS(II,JJ,2)+ &
                                XVAL3*WEIGHTS(II,JJ,3)+ &
                                XVAL4*WEIGHTS(II,JJ,4)
                                ITEMP=INT((XAVG-MINVAL)*253.0/RANGE)+1
                                IF (ITEMP < 0) ITEMP=0
                                IF (ITEMP > 255) ITEMP=255
                                BVAL=ITEMP
                                ACHAR=CHAR(BVAL)
                                WRITE(44,130,ADVANCE='NO') ACHAR
                            ENDDO L2128
                        ENDDO L2127
                    ENDIF
                ENDDO L2121
            ENDDO L2120
        ENDIF
        130   FORMAT(A1)
        CLOSE(44)
        RETURN
    END
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE WRITE_GROUNDWATER_FLOW(IUNIT)
        USE ERODE_GLOBALS
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: IUNIT
        INTEGER :: I,J
        REAL (8) :: TEMP
        OPEN(IUNIT,FILE='QQ.DAT',POSITION='APPEND')
        I=MX
        J=MY
        WRITE(IUNIT,440) I,J
        440   FORMAT(2I6)
        DO  I=1,MX
            DO  J=1,MY
                TEMP=FILTERED_GROUNDWATER_FLUX(I,J)
                WRITE(IUNIT,450) TEMP
                450   FORMAT(G12.5)
            ENDDO
        ENDDO
        CLOSE(IUNIT)
        RETURN
    END
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE WRITE_GROUNDWATER_ELEVATION(IUNIT)
        USE ERODE_GLOBALS
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: IUNIT
        INTEGER :: I,J
        REAL (8) :: TEMP
        OPEN(IUNIT,FILE='EWATER.DAT',POSITION='APPEND')
        I=MX
        J=MY
        WRITE(IUNIT,449) I,J
        449   FORMAT(2I6)
        DO  I=1,MX
            DO  J=1,MY
                TEMP=WATER_ELEVATION(I,J)/VERTICAL_SCALING
                WRITE(IUNIT,459) TEMP
                459   FORMAT(G12.5)
            ENDDO
        ENDDO
        CLOSE(IUNIT)
        RETURN
    END
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE FIND_GROUNDWATER_FLUX()
        USE ERODE_GLOBALS
        IMPLICIT NONE
        INTEGER :: I,J
        REAL (8) :: WMIN,WAVG,WMAX,QQQ,DISCHARGE_FROM_CELL
        WMIN=1.0E+25
        WMAX=-WMIN
        WAVG=0.0
        DO J=1,MY
            DO I=1,MX
                QQQ=DISCHARGE_FROM_CELL(I,J)
                IF (QQQ > WMAX) WMAX=QQQ
                IF (QQQ < WMIN) WMIN=QQQ
                WAVG=WAVG+QQQ
            ENDDO
        ENDDO
        WAVG=WAVG/(MX*MY)
        WRITE(OUTHIST,100) WMIN,WAVG,WMAX
        100   FORMAT(' LOCAL DISCHARGE, MIN=',G12.5,' AVG=',G12.5,' MAX=',G12.5)
        RETURN
    END
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE WRITE_ROUTED_DISCHARGE()
        USE ERODE_GLOBALS
        USE AREA_GLOBALS
        USE LAKE_GLOBALS
        IMPLICIT NONE
        INTEGER :: I,J, ITEMP
        LOGICAL :: CONNECTBASIN
        INTEGER(4) :: IFILE
        INTEGER (1) :: BVAL
        CHARACTER :: ACHAR
        REAL (8) :: MAXELEV,MINELEV,AAAA,BBBB
        REAL (8) :: ZMAX,ZMIN,QLOCAL
        IFILE=53
        CONNECTBASIN=.TRUE.
        RFILENAME='DSCHRG'//RNUM3//RNUM2//RNUM1//RSUFFIX
        OPEN(IFILE,FILE=RFILENAME)
        MAXELEV=-1.0E+25
        MINELEV=-MAXELEV
        DO J=1,MY
            DO I=1,MX
                IF (ROUTED_DISCHARGE(I,J) > 0.0) THEN
                    QLOCAL=LOG10(ROUTED_DISCHARGE(I,J))
                ELSE
                    QLOCAL=LOG10(0.01*D2X*CONVERT_DISCHARGE)
                ENDIF
                IF (QLOCAL > MAXELEV) MAXELEV=QLOCAL
                IF (QLOCAL < MINELEV) MINELEV=QLOCAL
            ENDDO
        ENDDO
        ZMAX=MAXELEV
        ZMIN=MINELEV
        WRITE(OUTHIST,200) ZMAX,ZMIN
        WRITE(*,200) ZMAX,ZMIN
        200   FORMAT('Q ZMAX=',G12.5,' ZMIN=',G12.5)
        AAAA=256.0/(MAXELEV-MINELEV)
        BBBB=-AAAA*MINELEV
        DO J=1,MY
            DO I=1,MX
                IF (ROUTED_DISCHARGE(I,J) > 0.0) THEN
                    QLOCAL=LOG10(ROUTED_DISCHARGE(I,J))
                ELSE
                    QLOCAL=LOG10(0.01*D2X*CONVERT_DISCHARGE)
                ENDIF
                ITEMP=INT(AAAA*QLOCAL+BBBB)
                IF (ITEMP < 0) ITEMP=0
                IF (ITEMP > 255) ITEMP=255
                BVAL=ITEMP
                ACHAR=CHAR(BVAL)
                WRITE(IFILE,130,ADVANCE='NO') ACHAR
                130   FORMAT(A1)
            ENDDO
        ENDDO
        CLOSE(53)
        RETURN
    END

    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE WRITE_SUBMERGED_LOCATIONS()
        USE ERODE_GLOBALS
        USE AREA_GLOBALS
        USE LAKE_GLOBALS
        IMPLICIT NONE
        INTEGER :: I,J
        INTEGER(4) :: IFILE
        INTEGER (1) :: BVAL
        CHARACTER :: ACHAR
        INTEGER(4) :: ITEMP
        IFILE=53
        RFILENAME='SUBMRG'//RNUM3//RNUM2//RNUM1//RSUFFIX
        OPEN(IFILE,FILE=RFILENAME)
        DO J=1,MY
            DO I=1,MX
                IF (SUBMERGED(I,J)) THEN
                    ITEMP=0
                ELSE
                    ITEMP=255
                ENDIF
                BVAL=ITEMP
                ACHAR=CHAR(BVAL)
                WRITE(IFILE,130,ADVANCE='NO') ACHAR
                130   FORMAT(A1)
            ENDDO
        ENDDO
        CLOSE(53)
        RETURN
    END
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE GRAD_DISCH_WRITE()
        USE ERODE_GLOBALS
        USE AREA_GLOBALS
        USE LAKE_GLOBALS
        IMPLICIT NONE
        REAL (8) :: REFGRAD,ELCOMP,EL_LOW,STEP_DISTANCE
        INTEGER :: I,J,ILOW,JLOW,ICOMP,JCOMP
        OPEN(87,FILE='GRAD_DISCH.DAT')
        WRITE (87,100) MX,MY
        100   FORMAT(I6,' ',I6)
        DO I=1,MX
            DO J=1,MY
                IF ((I > 1).AND.(I < MX).AND.(J > 1).AND.(J < MY)) THEN
                    ELCOMP=1.0E+25
                    IF (D8_GRADIENT(I,J) > 0.0) THEN
                        ILOW=I+DOWNSTREAM(FLOW_DIRECTION(I,J),1)
                        JLOW=J+DOWNSTREAM(FLOW_DIRECTION(I,J),2)
                        EL_LOW=ELEVATION(ILOW,JLOW)
                        IF (FLOW_DIRECTION(I-1,J-1) == 9) THEN
                            ELCOMP=ELEVATION(I-1,J-1)
                            ICOMP=I-1
                            JCOMP=J-1
                        ENDIF
                        IF (FLOW_DIRECTION(I-1,J) == 5) THEN
                            IF (ELEVATION(I-1,J) < ELCOMP) THEN
                                ELCOMP=ELEVATION(I-1,J)
                                ICOMP=I-1
                                JCOMP=J
                            ENDIF
                        ENDIF
                        IF (FLOW_DIRECTION(I-1,J+1) == 8) THEN
                            IF (ELEVATION(I-1,J+1) < ELCOMP) THEN
                                ELCOMP=ELEVATION(I-1,J+1)
                                ICOMP=I-1
                                JCOMP=J+1
                            ENDIF
                        ENDIF
                        IF (FLOW_DIRECTION(I,J-1) == 2) THEN
                            IF (ELEVATION(I,J-1) < ELCOMP) THEN
                                ELCOMP=ELEVATION(I,J-1)
                                ICOMP=I
                                JCOMP=J-1
                            ENDIF
                        ENDIF
                        IF (FLOW_DIRECTION(I,J+1) == 4) THEN
                            IF (ELEVATION(I,J+1) < ELCOMP) THEN
                                ELCOMP=ELEVATION(I,J+1)
                                ICOMP=I
                                JCOMP=J+1
                            ENDIF
                        ENDIF
                        IF (FLOW_DIRECTION(I+1,J-1) == 6) THEN
                            IF (ELEVATION(I+1,J-1) < ELCOMP) THEN
                                ELCOMP=ELEVATION(I+1,J-1)
                                ICOMP=I+1
                                JCOMP=J-1
                            ENDIF
                        ENDIF
                        IF (FLOW_DIRECTION(I+1,J) == 3) THEN
                            IF (ELEVATION(I+1,J) < ELCOMP) THEN
                                ELCOMP=ELEVATION(I+1,J)
                                ICOMP=I+1
                                JCOMP=J
                            ENDIF
                        ENDIF
                        IF (FLOW_DIRECTION(I+1,J+1) == 7) THEN
                            IF (ELEVATION(I+1,J+1) < ELCOMP) THEN
                                ELCOMP=ELEVATION(I+1,J+1)
                                ICOMP=I+1
                                JCOMP=J+1
                            ENDIF
                        ENDIF
                        IF (ELCOMP < 1.0E+24) THEN
                            STEP_DISTANCE=SQRT(FLOAT((ILOW-ICOMP)**2+(JLOW-JCOMP)**2))*CELL_SIZE
                            IF (STEP_DISTANCE > 0.0) THEN
                                REFGRAD=(ELCOMP-EL_LOW)/STEP_DISTANCE
                            ELSE
                                REFGRAD=0.0
                            ENDIF
                            IF (REFGRAD < 0.0) REFGRAD=0.0
                        ELSE
                            REFGRAD=D8_GRADIENT(I,J)
                        ENDIF
                    ENDIF
                    WRITE(87,200) D8_GRADIENT(I,J),REFGRAD,DISCHARGE(I,J),ROUTED_DISCHARGE(I,J)
                    200          FORMAT(4(G13.6,','))
                ELSE
                    WRITE(87,200) D8_GRADIENT(I,J),D8_GRADIENT(I,J),DISCHARGE(I,J),ROUTED_DISCHARGE(I,J)
                ENDIF
            ENDDO
        ENDDO
        CLOSE(87)
        RETURN
    END
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE WRITE_COLOR_SHADED_RELIEF_IMAGE()
        USE ERODE_GLOBALS
        IMPLICIT NONE
        REAL (8) :: HALFPI,PI,RADIAN,NOMANGLE,DTXX
        REAL (8) :: ETEMP(IMMX,JMMX)
        REAL (8) :: YDIV,SLX,SLY,SLOPE,SLMAX,SLSUN,ANGLE,SLOPEANGLE,NEWVALUE
        INTEGER :: I,J,IW,IE,JS,JN,SLOPENUMB,IMX,JMX,IX,JX
        REAL (8) ::     NOMSLMAX,RRRR,RANGE,SLOPEAVG,XDIV,MAXVAL,MINVAL
        INTEGER :: II,JJ,IISTART,IIEND,JJSTART,JJEND,XREFLECT,YREFLECT
        REAL (8) :: XVAL1,XVAL2,XVAL3,XVAL4,XAVG
        INTEGER (1) :: BVAL
        INTEGER :: ITEMP
        CHARACTER :: ACHAR,BCHAR,CCHAR
        !     **********************************************************************
        !      This writes out a raw shaded-relief image of the surface with color-coding dependent
        !        upon surface state (dust-covered, bedrock, ice)
        !     The file size depends on whether the lateral and/or top boundaries are
        !      periodic.  If both are false then the raw image size is (6*(MX-1),6*(MY-1))
        !     If the image is periodic in, say, the X dimension, then the horizontal image
        !       size increases to 6*(MX-1)+6*(MX/2+2) and similarly if the image is periodic
        !       in the Y dimension.  The image is in RGB bytes for each image pixel
        !     **********************************************************************

        IMX=MX
        JMX=MY
        IF (IS_X_PERIODIC) THEN
            XREFLECT=1
        ELSE
            XREFLECT=0
        ENDIF
        IF (IS_Y_PERIODIC) THEN
            YREFLECT=1
        ELSE
            YREFLECT=0
        ENDIF
        HALFPI =1.5708
        PI=3.141593
        RADIAN=0.01745329
        NOMANGLE=25.0
        SLMAX=-1.0E+25
        SLOPEAVG=0.0
        SLOPENUMB= 0
        DO  J = 2,JMX-1
            DO  I = 2,IMX-1
                IW=I-1
                IE=I+1
                JS=J+1
                JN=J-1
                SLX = ( ELEVATION(IE,J)- ELEVATION(IW,J))/2.0
                SLY = ( ELEVATION(I,JN)- ELEVATION(I,JS))/2.0
                SLOPE =SQRT(SLX**2+SLY**2)
                IF (SLOPE > SLMAX) SLMAX = SLOPE
                SLOPENUMB = SLOPENUMB+1
                SLOPEAVG = SLOPEAVG+SLOPE
            ENDDO
        ENDDO
        NOMSLMAX = SIN(HALFPI-NOMANGLE*RADIAN)/  &
        COS(HALFPI-NOMANGLE*RADIAN)
        DTXX = SLMAX/NOMSLMAX
        SLMAX = NOMSLMAX
        SLOPEAVG = SLOPEAVG/(SLOPENUMB*DTXX)
        SLSUN=HALFPI-DATAN(SLMAX)
        DO  J = 1,JMX
            DO  I = 1,IMX
                IW=I-1
                IE=I+1
                JS=J+1
                JN=J-1
                YDIV=2.0*DTXX
                XDIV=2.0*DTXX
                IF (JS > JMX) THEN
                    YDIV=1.0*DTXX
                    JS=JMX
                ENDIF
                IF (JN < 1) THEN
                    YDIV=1.0*DTXX
                    JN=1
                ENDIF
                IF (IW < 1) THEN
                    XDIV=1.0*DTXX
                    IW=1
                ENDIF
                IF (IE > IMX) THEN
                    XDIV=1.0*DTXX
                    IE=IMX
                ENDIF
                SLX = ( ELEVATION(IE,J)- ELEVATION(IW,J))/XDIV
                SLY = ( ELEVATION(I,JS)- ELEVATION(I,JN))/YDIV
                SLOPE =SQRT(SLX**2+SLY**2)
                SLOPEANGLE= DATAN(SLOPE)
                IF ((SLY == 0.0).OR. (SLX == 0.0)) THEN
                    IF (SLY == 0.0) THEN
                        IF (SLX > 0.0) THEN
                            ANGLE =  HALFPI
                        ELSE
                            ANGLE = -HALFPI
                        ENDIF
                    ELSE
                        IF (SLY > 0.0) THEN
                            ANGLE = 0.0
                        ELSE
                            ANGLE = PI
                        ENDIF
                    ENDIF
                ELSE
                    ANGLE = DATAN ((SLX) / (SLY))
                    IF (((SLY)  <  0.0).AND.((SLX)  <  0.0)) THEN
                        ANGLE = ANGLE + PI
                    ELSE
                        IF ((SLY)  <  0.0) THEN
                            ANGLE = ANGLE + PI
                        ENDIF
                    ENDIF
                ENDIF
                ANGLE = ANGLE-HALFPI
                NEWVALUE =   &
                (COS(SLSUN)*COS(SLOPEANGLE)+SIN(SLSUN) &
                *SIN(SLOPEANGLE)*  &
                COS(ANGLE))
                IF (NEWVALUE < 0.0) NEWVALUE = 0.0

                ETEMP(I,J) = NEWVALUE
            ENDDO
        ENDDO
        MAXVAL = -1.0E+25
        MINVAL = -MAXVAL
        DO  J = 1,JMX
            DO  I = 1,IMX
                RRRR = ETEMP(I,J)
                IF (RRRR < MINVAL) MINVAL = RRRR
                IF (RRRR >  MAXVAL) MAXVAL = RRRR
            ENDDO
        ENDDO
        RANGE=MAXVAL-MINVAL
        RNUM1=CHAR(IFILE1)
        IF (IFILE2 > 57)  THEN
            RNUM2=CHAR(IFILE2+40)
        ELSE
            RNUM2=CHAR(IFILE2)
        ENDIF
        RNUM3=CHAR(IFILE3)
!        IFILE1=IFILE1+1
!        IF (IFILE1 > 57) THEN
!            IFILE1=48
!            IFILE2=IFILE2+1
!        ENDIF
        RFILENAME=RPREFIX//'_C'//RNUM3//RNUM2//RNUM1//RSUFFIX
        OPEN(44,FILE=RFILENAME)
        IF (YREFLECT > 0) THEN
            JJSTART=JMX-JMX/4
            JJEND=JMX
            L1120: DO  J=JJSTART,JJEND
                JX=J
                L1121: DO  JJ=1,2
                    IF (XREFLECT > 0) THEN
                        IISTART=IMX-IMX/4
                        IIEND=IMX
                        L1123: DO  I=IISTART,IIEND
                            IX=I
                            XVAL1=ETEMP(I,J)
                            IF (I == IIEND) THEN
                                XVAL2=ETEMP(1,J)
                                IF (J == JJEND) THEN
                                    XVAL4=ETEMP(1,1)
                                ELSE
                                    XVAL4=ETEMP(1,J+1)
                                ENDIF
                            ELSE
                                XVAL2=ETEMP(I+1,J)
                                IF (J == JJEND) THEN
                                    XVAL4=ETEMP(I+1,1)
                                ELSE
                                    XVAL4=ETEMP(I+1,J+1)
                                ENDIF
                            ENDIF
                            IF (J == JJEND) THEN
                                XVAL3=ETEMP(I,1)
                            ELSE
                                XVAL3=ETEMP(I,J+1)
                            ENDIF
                            L1124: DO  II=1,2
                                XAVG=XVAL1*WEIGHTS(II,JJ,1)+ &
                                XVAL2*WEIGHTS(II,JJ,2)+ &
                                XVAL3*WEIGHTS(II,JJ,3)+ &
                                XVAL4*WEIGHTS(II,JJ,4)
                                ITEMP=INT((XAVG-MINVAL)*253.0/RANGE)+1
                                IF (ITEMP < 0) ITEMP=0
                                IF (ITEMP > 255) ITEMP=255
                                BVAL=ITEMP
                                ACHAR=CHAR(BVAL)
                                BCHAR=CHAR(BVAL)
                                CCHAR=CHAR(BVAL)
                                IF (IS_SEDIMENT_COVERED(IX,JX)) THEN
                                    CCHAR=CHAR(0)
                                ELSE
                                    IF (.NOT.IS_ROCK_SURFACE(IX,JX)) THEN
                                        BCHAR=CHAR(0)
                                    ENDIF
                                ENDIF
                                WRITE(44,130,ADVANCE='NO') ACHAR,BCHAR,CCHAR
                            ENDDO L1124
                        ENDDO L1123
                    ENDIF
                    L1125: DO  I=1,IMX-1
                        IX=I
                        XVAL1=ETEMP(I,J)
                        XVAL2=ETEMP(I+1,J)
                        IF (J == JJEND) THEN
                            XVAL4=ETEMP(I+1,1)
                            XVAL3=ETEMP(I,1)
                        ELSE
                            XVAL4=ETEMP(I+1,J+1)
                            XVAL3=ETEMP(I,J+1)
                        ENDIF
                        L1126: DO  II=1,2
                            XAVG=XVAL1*WEIGHTS(II,JJ,1)+  &
                            XVAL2*WEIGHTS(II,JJ,2)+ &
                            XVAL3*WEIGHTS(II,JJ,3)+ &
                            XVAL4*WEIGHTS(II,JJ,4)
                            ITEMP=INT((XAVG-MINVAL)*253.0/RANGE)+1
                            IF (ITEMP < 0) ITEMP=0
                            IF (ITEMP > 255) ITEMP=255
                            BVAL=ITEMP
                            ACHAR=CHAR(BVAL)
                            BCHAR=CHAR(BVAL)
                            CCHAR=CHAR(BVAL)
                            IF (IS_SEDIMENT_COVERED(IX,JX)) THEN
                                CCHAR=CHAR(0)
                            ELSE
                                IF (.NOT.IS_ROCK_SURFACE(IX,JX)) THEN
                                    BCHAR=CHAR(0)
                                ENDIF
                            ENDIF
                            WRITE(44,130,ADVANCE='NO') ACHAR,BCHAR,CCHAR
                        ENDDO L1126
                    ENDDO L1125
                    IF (XREFLECT > 0) THEN
                        IISTART=0
                        IIEND=IMX/4
                        L1127: DO  I=IISTART,IIEND
                            XVAL2=ETEMP(I+1,J)
                            IF (I == IISTART) THEN
                                IX=1
                                XVAL1=ETEMP(IMX,J)
                                IF (J == JJEND) THEN
                                    XVAL3=ETEMP(IMX,1)
                                ELSE
                                    XVAL3=ETEMP(IMX,J+1)
                                ENDIF
                            ELSE
                                IX=I
                                XVAL1=ETEMP(I+1,J)
                                IF (J == JJEND) THEN
                                    XVAL3=ETEMP(I,1)
                                ELSE
                                    XVAL3=ETEMP(I,J+1)
                                ENDIF
                            ENDIF
                            IF (J == JJEND) THEN
                                XVAL4=ETEMP(I+1,1)
                            ELSE
                                XVAL4=ETEMP(I+1,J+1)
                            ENDIF
                            L1128: DO  II=1,2
                                XAVG=XVAL1*WEIGHTS(II,JJ,1)+ &
                                XVAL2*WEIGHTS(II,JJ,2)+  &
                                XVAL3*WEIGHTS(II,JJ,3)+ &
                                XVAL4*WEIGHTS(II,JJ,4)
                                ITEMP=INT((XAVG-MINVAL)*253.0/RANGE)+1
                                IF (ITEMP < 0) ITEMP=0
                                IF (ITEMP > 255) ITEMP=255
                                BVAL=ITEMP
                                ACHAR=CHAR(BVAL)
                                BCHAR=CHAR(BVAL)
                                CCHAR=CHAR(BVAL)
                                IF (IS_SEDIMENT_COVERED(IX,JX)) THEN
                                    CCHAR=CHAR(0)
                                ELSE
                                    IF (.NOT.IS_ROCK_SURFACE(IX,JX)) THEN
                                        BCHAR=CHAR(0)
                                    ENDIF
                                ENDIF
                                WRITE(44,130,ADVANCE='NO') ACHAR,BCHAR,CCHAR
                            ENDDO L1128
                        ENDDO L1127
                    ENDIF
                ENDDO L1121
            ENDDO L1120
        ENDIF
        L120: DO  J=1,JMX-1
            JX=J
            L121: DO  JJ=1,2
                IF (XREFLECT > 0) THEN
                    IISTART=IMX-IMX/4
                    IIEND=IMX
                    L123: DO  I=IISTART,IIEND
                        IX=I
                        XVAL1=ETEMP(I,J)
                        IF (I == IIEND) THEN
                            XVAL2=ETEMP(1,J)
                            XVAL4=ETEMP(1,J+1)
                        ELSE
                            XVAL2=ETEMP(I+1,J)
                            XVAL4=ETEMP(I+1,J+1)
                        ENDIF
                        XVAL3=ETEMP(I,J+1)
                        L124: DO  II=1,2
                            XAVG=XVAL1*WEIGHTS(II,JJ,1)+  &
                            XVAL2*WEIGHTS(II,JJ,2)+ &
                            XVAL3*WEIGHTS(II,JJ,3)+ &
                            XVAL4*WEIGHTS(II,JJ,4)
                            ITEMP=INT((XAVG-MINVAL)*253.0/RANGE)+1
                            IF (ITEMP < 0) ITEMP=0
                            IF (ITEMP > 255) ITEMP=255
                            BVAL=ITEMP
                            ACHAR=CHAR(BVAL)
                            BCHAR=CHAR(BVAL)
                            CCHAR=CHAR(BVAL)
                            IF (IS_SEDIMENT_COVERED(IX,JX)) THEN
                                CCHAR=CHAR(0)
                            ELSE
                                IF (.NOT.IS_ROCK_SURFACE(IX,JX)) THEN
                                    BCHAR=CHAR(0)
                                ENDIF
                            ENDIF
                            WRITE(44,130,ADVANCE='NO') ACHAR,BCHAR,CCHAR
                        ENDDO L124
                    ENDDO L123
                ENDIF
                L125: DO  I=1,IMX-1
                    IX=I
                    XVAL1=ETEMP(I,J)
                    XVAL2=ETEMP(I+1,J)
                    XVAL4=ETEMP(I+1,J+1)
                    XVAL3=ETEMP(I,J+1)
                    L126: DO  II=1,2
                        XAVG=XVAL1*WEIGHTS(II,JJ,1)+ &
                        XVAL2*WEIGHTS(II,JJ,2)+ &
                        XVAL3*WEIGHTS(II,JJ,3)+ &
                        XVAL4*WEIGHTS(II,JJ,4)
                        ITEMP=INT((XAVG-MINVAL)*253.0/RANGE)+1
                        IF (ITEMP < 0) ITEMP=0
                        IF (ITEMP > 255) ITEMP=255
                        BVAL=ITEMP
                        ACHAR=CHAR(BVAL)
                        BCHAR=CHAR(BVAL)
                        CCHAR=CHAR(BVAL)
                        IF (IS_SEDIMENT_COVERED(IX,JX)) THEN
                            CCHAR=CHAR(0)
                        ELSE
                            IF (.NOT.IS_ROCK_SURFACE(IX,JX)) THEN
                                BCHAR=CHAR(0)
                            ENDIF
                        ENDIF
                        WRITE(44,130,ADVANCE='NO') ACHAR,BCHAR,CCHAR
                    ENDDO L126
                ENDDO L125
                IF (XREFLECT > 0) THEN
                    IISTART=0
                    IIEND=IMX/4
                    L127: DO  I=IISTART,IIEND
                        XVAL2=ETEMP(I+1,J)
                        IF (I == IISTART) THEN
                            IX=1
                            XVAL1=ETEMP(IMX,J)
                            XVAL3=ETEMP(IMX,J+1)
                        ELSE
                            IX=I
                            XVAL1=ETEMP(I+1,J)
                            XVAL3=ETEMP(I+1,J+1)
                        ENDIF
                        XVAL4=ETEMP(I+1,J+1)
                        L128: DO  II=1,2
                            XAVG=XVAL1*WEIGHTS(II,JJ,1)+ &
                            XVAL2*WEIGHTS(II,JJ,2)+ &
                            XVAL3*WEIGHTS(II,JJ,3)+ &
                            XVAL4*WEIGHTS(II,JJ,4)
                            ITEMP=INT((XAVG-MINVAL)*253.0/RANGE)+1
                            IF (ITEMP < 0) ITEMP=0
                            IF (ITEMP > 255) ITEMP=255
                            BVAL=ITEMP
                            ACHAR=CHAR(BVAL)
                            BCHAR=CHAR(BVAL)
                            CCHAR=CHAR(BVAL)
                            IF (IS_SEDIMENT_COVERED(IX,JX)) THEN
                                CCHAR=CHAR(0)
                            ELSE
                                IF (.NOT.IS_ROCK_SURFACE(IX,JX)) THEN
                                    BCHAR=CHAR(0)
                                ENDIF
                            ENDIF
                            WRITE(44,130,ADVANCE='NO') ACHAR,BCHAR,CCHAR
                        ENDDO L128
                    ENDDO L127
                ENDIF
            ENDDO L121
        ENDDO L120
        IF (YREFLECT > 0) THEN
            JJSTART=0
            JJEND=JMX/4
            L2120: DO  J=JJSTART,JJEND
                L2121: DO  JJ=1,2
                    IF (XREFLECT > 0) THEN
                        IISTART=IMX-IMX/4
                        IIEND=IMX
                        L2123: DO  I=IISTART,IIEND
                            IX=I
                            IF (J == JJSTART) THEN
                                JX=1
                                XVAL1=ETEMP(I,JMX)
                            ELSE
                                JX=J
                                XVAL1=ETEMP(I,J)
                            ENDIF
                            IF (I == IIEND) THEN
                                IF (J == JJSTART) THEN
                                    XVAL2=ETEMP(1,JMX)
                                ELSE
                                    XVAL2=ETEMP(1,J)
                                ENDIF
                                XVAL4=ETEMP(1,J+1)
                            ELSE
                                IF (J == JJSTART) THEN
                                    XVAL2=ETEMP(I+1,JMX)
                                ELSE
                                    XVAL2=ETEMP(I+1,J)
                                ENDIF
                                XVAL4=ETEMP(I+1,J+1)
                            ENDIF
                            XVAL3=ETEMP(I,J+1)
                            L2124: DO  II=1,2
                                XAVG=XVAL1*WEIGHTS(II,JJ,1)+  &
                                XVAL2*WEIGHTS(II,JJ,2)+ &
                                XVAL3*WEIGHTS(II,JJ,3)+ &
                                XVAL4*WEIGHTS(II,JJ,4)
                                ITEMP=INT((XAVG-MINVAL)*253.0/RANGE)+1
                                IF (ITEMP < 0) ITEMP=0
                                IF (ITEMP > 255) ITEMP=255
                                BVAL=ITEMP
                                ACHAR=CHAR(BVAL)
                                BCHAR=CHAR(BVAL)
                                CCHAR=CHAR(BVAL)
                                IF (IS_SEDIMENT_COVERED(IX,JX)) THEN
                                    CCHAR=CHAR(0)
                                ELSE
                                    IF (.NOT.IS_ROCK_SURFACE(IX,JX)) THEN
                                        BCHAR=CHAR(0)
                                    ENDIF
                                ENDIF
                                WRITE(44,130,ADVANCE='NO') ACHAR,BCHAR,CCHAR
                            ENDDO L2124
                        ENDDO L2123
                    ENDIF
                    L2125: DO  I=1,IMX-1
                        IX=I
                        IF (J == JJSTART) THEN
                            JX=1
                            XVAL1=ETEMP(I,JMX)
                            XVAL2=ETEMP(I+1,JMX)
                        ELSE
                            JX=J
                            XVAL1=ETEMP(I,J)
                            XVAL2=ETEMP(I+1,J)
                        ENDIF
                        XVAL4=ETEMP(I+1,J+1)
                        XVAL3=ETEMP(I,J+1)
                        L2126: DO  II=1,2
                            XAVG=XVAL1*WEIGHTS(II,JJ,1)+ &
                            XVAL2*WEIGHTS(II,JJ,2)+ &
                            XVAL3*WEIGHTS(II,JJ,3)+ &
                            XVAL4*WEIGHTS(II,JJ,4)
                            ITEMP=INT((XAVG-MINVAL)*253.0/RANGE)+1
                            IF (ITEMP < 0) ITEMP=0
                            IF (ITEMP > 255) ITEMP=255
                            BVAL=ITEMP
                            ACHAR=CHAR(BVAL)
                            BCHAR=CHAR(BVAL)
                            CCHAR=CHAR(BVAL)
                            IF (IS_SEDIMENT_COVERED(IX,JX)) THEN
                                CCHAR=CHAR(0)
                            ELSE
                                IF (.NOT.IS_ROCK_SURFACE(IX,JX)) THEN
                                    BCHAR=CHAR(0)
                                ENDIF
                            ENDIF
                            WRITE(44,130,ADVANCE='NO') ACHAR,BCHAR,CCHAR
                        ENDDO L2126
                    ENDDO L2125
                    IF (XREFLECT > 0) THEN
                        IISTART=0
                        IIEND=IMX/4
                        L2127: DO  I=IISTART,IIEND
                            IF (J == JJSTART) THEN
                                JX=1
                                XVAL2=ETEMP(I+1,JMX)
                            ELSE
                                JX=J
                                XVAL2=ETEMP(I+1,J)
                            ENDIF
                            IF (I == IISTART) THEN
                                IX=1
                                IF (J == JJSTART) THEN
                                    XVAL1=ETEMP(IMX,JMX)
                                ELSE
                                    XVAL1=ETEMP(IMX,J)
                                ENDIF
                                XVAL3=ETEMP(IMX,J+1)
                            ELSE
                                IX=I
                                IF (J == JJSTART) THEN
                                    XVAL1=ETEMP(I,JMX)
                                ELSE
                                    XVAL1=ETEMP(I,J)
                                ENDIF
                                XVAL3=ETEMP(I,J+1)
                            ENDIF
                            XVAL4=ETEMP(I+1,J+1)
                            L2128: DO  II=1,2
                                XAVG=XVAL1*WEIGHTS(II,JJ,1)+  &
                                XVAL2*WEIGHTS(II,JJ,2)+ &
                                XVAL3*WEIGHTS(II,JJ,3)+ &
                                XVAL4*WEIGHTS(II,JJ,4)
                                ITEMP=INT((XAVG-MINVAL)*253.0/RANGE)+1
                                IF (ITEMP < 0) ITEMP=0
                                IF (ITEMP > 255) ITEMP=255
                                BVAL=ITEMP
                                ACHAR=CHAR(BVAL)
                                BCHAR=CHAR(BVAL)
                                CCHAR=CHAR(BVAL)
                                IF (IS_SEDIMENT_COVERED(IX,JX)) THEN
                                    CCHAR=CHAR(0)
                                ELSE
                                    IF (.NOT.IS_ROCK_SURFACE(IX,JX)) THEN
                                        BCHAR=CHAR(0)
                                    ENDIF
                                ENDIF
                                WRITE(44,130,ADVANCE='NO') ACHAR,BCHAR,CCHAR
                            ENDDO L2128
                        ENDDO L2127
                    ENDIF
                ENDDO L2121
            ENDDO L2120
        ENDIF
        130   FORMAT(3A1)
        CLOSE(44)
        RETURN
    END
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        SUBROUTINE WRITE_LAVA_INFO()
        USE ERODE_GLOBALS
        USE LAVA_GLOBALS
        IMPLICIT NONE
        INTEGER :: I,J,IIII(IMMX),JJJJ(IMMX)
        !      *********************************************************************
        !       writes matrix of integers indicating whether given matrix locations
        !         expose lava (1) or are not lava covered (0)
        !      *********************************************************************
        OPEN(77,FILE='LAVA.DAT',POSITION='APPEND')
        OPEN(67,FILE='LACTIVE.DAT',POSITION='APPEND')
        WRITE(77,130) MX,MY
        WRITE(67,130) MX,MY
        130 FORMAT(2I5)
        DO  I=1,MX
            DO  J=1,MY
                IF (IS_LAVA_COVERED(I,J)) THEN
                    IIII(J)=1
                ELSE
                    IIII(J)=0
                ENDIF
                IF (ACTIVE_LAVA_FLOW(I,J)) THEN
                    JJJJ(J)=1
                ELSE
                    JJJJ(J)=0
                ENDIF
            ENDDO
            WRITE(77,110)(IIII(J),J=1,MY)
            WRITE(67,110)(JJJJ(J),J=1,MY)
            110     FORMAT(256I1)
        ENDDO
        CLOSE(77)
        CLOSE(67)
        RETURN
    END
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBROUTINE WRITE_LAVA_AGES()
        USE ERODE_GLOBALS
        USE LAVA_GLOBALS
        IMPLICIT NONE
        INTEGER :: I,J
        INTEGER (1) :: BVAL
        CHARACTER :: ACHAR
        REAL (8) :: MAXIAGE,MINIAGE,AAAA,BBBB,DAGE
        INTEGER :: ITEMP
        !     **********************************************************************
        !      writes out a raw image file of elevations normalized to 256 byte range
        !         and an information file on the iteration and actual elevation range
        !     **********************************************************************
        !     **********************************************************************
        !      image.raw is a file of raw byte-normalized elevation data -- an x,y
        !         image of relative elevations normalized to 256 values, 0 is lowest
        !         and 255 is highest elevation
        !     **********************************************************************
        OPEN(71,FILE='LAGE.DAT',POSITION='APPEND')
        OPEN(70,FILE='LAGE.RAW')
        !     **********************************************************************
        !      image.dat is a record of minimum elevation, time, and elevation range
        !         for each output record in image.raw
        !     **********************************************************************
        MAXIAGE=-1.0E+25
        MINIAGE=-MAXIAGE
        WRITE (71,200) MX,MY
        200   FORMAT(2I5)
        DO  I=1,MX
            DO  J=1,MY
                WRITE(71,210) ERUPTION_AGE(I,J)
                210   FORMAT (I8)
                IF (ERUPTION_AGE(I,J) < 1) ERUPTION_AGE(I,J)=1
                DAGE=DLOG(DFLOAT(ERUPTION_AGE(I,J)))
                IF (DAGE > MAXIAGE) MAXIAGE=DAGE
                IF (DAGE < MINIAGE) MINIAGE=DAGE
            ENDDO
        ENDDO
        AAAA=256.0/(MAXIAGE-MINIAGE)
        BBBB=-AAAA*MINIAGE
        DO  J=1,MY
            DO  I=1,MX
                DAGE=DLOG(DFLOAT(ERUPTION_AGE(I,J)))
                ITEMP=INT(AAAA*DAGE+BBBB)
                IF (ITEMP < 0) ITEMP=0
                IF (ITEMP > 255) ITEMP=255
                ITEMP=255-ITEMP
                BVAL=ITEMP
                ACHAR=CHAR(BVAL)
                WRITE(70,130,ADVANCE='NO') ACHAR
                130   FORMAT(A1)
            ENDDO
        ENDDO
        CLOSE(70)
        CLOSE(71)
        RETURN
    END


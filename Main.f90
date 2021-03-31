    INCLUDE 'Module_Constants_sub.f90'
    INCLUDE 'NumericalProbe_sub.f90'
    INCLUDE 'AbsorbingCoefficient_sub.f90'
    !INCLUDE 'ExchangeInterfaceData_sub.f90'
    INCLUDE 'ExchangeInterfaceDataNew_sub.f90'
    !INCLUDE 'ImportPerturbation_sub.f90'
    INCLUDE 'GetConservativeVariables_sub.f90'
    !INCLUDE 'SolidPMLSource_sub.f90'
    !INCLUDE 'GetMaxResidual_sub.f90'
    INCLUDE 'GetMeanFlow_sub.f90'
    INCLUDE 'GetMeshSize_sub.f90'
    INCLUDE 'GetOriginalVariables_sub.f90'
    INCLUDE 'InitializeAcousticField_sub.f90'
    INCLUDE 'LDDRK_sub.f90'
    INCLUDE 'Mesh_sub.f90'
    INCLUDE 'MeshSmoothing_sub.f90'
    INCLUDE 'Models_sub.f90'
    !INCLUDE 'BodyMoving_sub.f90'
    INCLUDE 'TransformCoordinate_sub.f90'
    !INCLUDE 'GetEffectMatrix_sub.f90'
    !INCLUDE 'GetDuDt_sub.f90'
    !INCLUDE 'GetNormalGradient_sub.f90'
    !INCLUDE 'ForcingSolver_sub.f90'
    !INCLUDE 'CorrectAcousticField_sub.f90'
    !INCLUDE 'Distribution_sub.f90'
!    INCLUDE 'lapack_LINUX.a'
!    INCLUDE 'blas_LINUX.a'
    !INCLUDE 'SOR_sub.f90'
    !INCLUDE 'DooliteSolve_sub.f90'
    INCLUDE 'OutputToTextFile_sub.f90'
    INCLUDE 'Filter_sub.f90'
    INCLUDE 'ReviseConserVariable_sub.f90'
    !INCLUDE 'CholeskyDecompsition_sub.f90'
    !INCLUDE 'SparseSolver_sub.f90'
    !INCLUDE 'MPI_CholeskyDecompsition_sub.f90'
    !INCLUDE 'PCGM.f90'
    !INCLUDE 'PCGSM.f90'
    INCLUDE 'PrintParameterSet_sub.f90'
    INCLUDE 'InitialModule_sub.f90'
    !INCLUDE 'OutputRestartFlowField_sub.f90'
    !INCLUDE 'RestartFlowField_sub.f90'
    ! 
    INCLUDE 'ObtainKineticEnergy_sub.f90'
    !INCLUDE 'GetEffectRange_sub.f90'
    !INCLUDE 'GetBoundaryForceDensity_sub.f90'
    !INCLUDE 'GetWallBCsError_sub.f90'
    !INCLUDE 'GetRHS_sub.f90'
    !INCLUDE 'SolidPenalizedSourceTerm_Sub.f90'
    ! 
    !INCLUDE 'LocalFilter_sub.f90'
    !INCLUDE 'AdaptiveDiscontinuityCapturingFilter_sub.f90'
PROGRAM AcousticScatteringParallelComputing
    !---------------------------------------------------
    !-----------------program introduction -----------------
    ! 
    ! 
    ! 
    !----------------------------------------------------
    !......
    !......
    !-----var------------------------------------------------------
    !-----iables------------------------------------------------------
    ! EffectRange:the influence width of Delta function for each lagrange
    !       point
    ! NormalVector:  the normal vector for each lagrange point
    ! Shapes:  the surface point for analysis model
    ! ShapesForP: shape points for pressure interpolation of solid surface
    ! Cell_S:  area of length for each cell of analysis model
    ! AV:  effect matrices for velocity
    ! AT:  effect matrices for temperature
    ! LagForce:  lagrange forcing
    ! SurfaceUVW: the fluid velocity for solid surface lagrange points
    ! Aprocessor:  (local)Effect matrices for each processor
    ! Atotal: temperal (global) effect matrices
    ! DuDt0: (local)velocity source for each processor	
    ! DuDt1: temperal (global) velocity source
    ! DuDt: (global) velocity source
    ! .....................
    ! MeshX, MeshY, MeshZ: Cartesian Mesh grid coordinates
    ! DeltaX, DeltaY, DeltaZ: Cartesian grid size for each direction
    ! U0, V0, W0, P0, ROU0: the background flow field
    ! U, V, W, P, ROU: acoustic field variables
    ! Q1~4, F1~4, G1~4, H1~4: conservative variables
    ! .....................
    ! Sigma~ : the absorbing coefficient for each side
    ! Computational zone 
    ! 	-----------------------------------------------------
    !   -       -               PML: SigmaY2               --
    !   -----------------------------------------------------
    !   -   PML	-                                   -   PML -
    !   -   S   -           Main                    -   S   -
    !   -   i   -               Zone                -   i   -
    !   -   g   -                                   -   g   -
    !   -   m   -                                   -   m   -	
    !   -   aX1	-                                   -   aX2 -
    !   -----------------------------------------------------
    !   -       -             PML: SigmaY2          -       -
    !   -----------------------------------------------------
    !.......................
    ! Au1~, Au2~, Au3~: auxiliary variables for each direction
    ! S1~4: source term for PML equations
    !-----intro------------------------------------------------------
    !-----duction------------------------------------------------------
    ! 
    USE MPI
    USE Constants_Model
    USE BackgroundFlowField
    USE VariableArtificialSelectiveDamping
    IMPLICIT NONE
    INTEGER:: I, J, k
    ! for lagrange points
    INTEGER:: PLN
    INTEGER, ALLOCATABLE:: EffectRange( :, : )
    REAL(KIND=8), ALLOCATABLE:: NormalVector( :, : ), Shapes( :, : ), BodyCentre( :, : )
    REAL(KIND=8), ALLOCATABLE:: BodyCentrePosition( :, : )
    REAL(KIND=8), ALLOCATABLE:: BodyCentreVeloci( :, : )
    REAL(KIND=8), ALLOCATABLE:: BodyCentreAcce( :, : )
    REAL(KIND=8), ALLOCATABLE:: TangentialVector( :, : )
    REAL(KIND=8), ALLOCATABLE:: InitialShapes( :, : )
    REAL(KIND=8), ALLOCATABLE:: ShapesForP( :, : ), Cell_S( : )
    REAL(KIND=8), ALLOCATABLE:: AV( :, : ), AT( :, : )
    REAL(KIND=8), ALLOCATABLE:: LagForce( :, : ), SurfUVW( :, : )
    REAL(KIND=8), ALLOCATABLE:: AProcessor( :, : ), ATotal( :, : )
    REAL(KIND=8), ALLOCATABLE:: DuDt0( : ), DuDt1( : ), DuDt( : )
    REAL(KIND=8), ALLOCATABLE:: DvDt0( : ), DvDt1( : ), DvDt( : )
    REAL(KIND=8), ALLOCATABLE:: DwDt0( : ), DwDt1( : ), DwDt( : )
    REAL(KIND=8), ALLOCATABLE:: DTDt0( : ), DTDt1( : ), DTDt( : )
    REAL(KIND=8), ALLOCATABLE:: DRouDt0( : ), DRouDt1( : ), DRouDt( : )
    INTEGER, ALLOCATABLE:: DuDtFlag( : ), DuDtFlagRoot( : )
    REAL(KIND=8), ALLOCATABLE:: X01( :, : ), X02( :, : ), X03( :, : ), X04( :, : )
    ! 
    ! for Cartesian grid information
    REAL(KIND=8), ALLOCATABLE:: MeshX( : ), MeshY( : ), MeshZ( : )
    REAL(KIND=8), ALLOCATABLE:: DeltaX( : ), DeltaY( : ), DeltaZ( : )
    ! 
    ! 
    ! acoustic field: nonconservative variables and conservative variables
    ! output to file for main processor
    ! 
    ! calculation
    REAL(KIND=8), ALLOCATABLE:: U( :, :, : ), V( :, :, : ), W( :, :, : )
    REAL(KIND=8), ALLOCATABLE:: Omiga( :, :, : ), Mach( :, :, : )
    REAL(KIND=8), ALLOCATABLE:: P( :, :, : ), ROU( :, :, : ), Te( :, :, : )
    REAL(KIND=8), ALLOCATABLE:: Q1( :, :, : ), Q2( :, :, : )
    REAL(KIND=8), ALLOCATABLE:: Q3( :, :, : ), Q4( :, :, : ), Q5( :, :, : )
    REAL(KIND=8), ALLOCATABLE:: F1( :, :, : ), F2( :, :, : )
    REAL(KIND=8), ALLOCATABLE:: F3( :, :, : ), F4( :, :, : ), F5( :, :, : )
    REAL(KIND=8), ALLOCATABLE:: G1( :, :, : ), G2( :, :, : )
    REAL(KIND=8), ALLOCATABLE:: G3( :, :, : ), G4( :, :, : ), G5( :, :, : )
    REAL(KIND=8), ALLOCATABLE:: H1( :, :, : ), H2( :, :, : )
    REAL(KIND=8), ALLOCATABLE:: H3( :, :, : ), H4( :, :, : ), H5( :, :, : )
    !REAL(KIND=8), ALLOCATABLE:: MaskSignal( :, :, : ), BodyForceDistriX( :, :, : ), &
    !&                                                  BodyForceDistriY( :, :, : )
    REAL(KIND=8), ALLOCATABLE:: BodyForceDistriX_Out( :, :, : ), BodyForceDistriY_Out( :, :, : )
    REAL(KIND=8), ALLOCATABLE:: CL_Temp(:), Cd_Temp(:)
    REAL(KIND=8), ALLOCATABLE:: CL_output(:), Cd_output(:)
    REAL(KIND=8), ALLOCATABLE:: Q2_Pre( :, :, : )
    REAL(KIND=8), ALLOCATABLE:: Q3_Pre( :, :, : )
    !REAL(KIND=8), ALLOCATABLE:: InertialForceX( : ), InertialForceY( : )
    ! 
    ! ------shoch capturing filter related. 
    REAL(KIND=8), ALLOCATABLE:: Shock_Sigma( :, :, : )
    ! 
    ! PML zone 
    REAL(KIND=8):: SigmaX1( -(NPML-1) : 0 )
    REAL(KIND=8):: SigmaX2( 1 : NPML )
    REAL(KIND=8):: SigmaY1( -(NPML-1) : 0 )
    REAL(KIND=8):: SigmaY2( 1 : NPML )
    REAL(KIND=8):: SigmaZ1( -(NPML-1) : 0 )
    REAL(KIND=8):: SigmaZ2( 1 : NPML )
    !!
    !MPI variables definition
    INTEGER:: IERR, NUMPROCS
    INTEGER:: MYID, MYROOT
    INTEGER:: MYLEFT, MYRIGHT, MYUPPER, MYLOWER, MYFORWARD, MYREAR
    INTEGER:: PX, PY, PZ
    INTEGER:: HTYPE, VTYPE
    INTEGER:: XTYPE, YTYPE, ZTYPE
    INTEGER:: YTYPE1
    INTEGER, ALLOCATABLE:: STATUS( :, : ), REQ(:)
    INTEGER:: COUNT1
    INTEGER,ALLOCATABLE:: BLOCKLENS( : ), INDICES( : )
    INTEGER:: SENDCNT
    INTEGER, ALLOCATABLE:: RECVCNT( : )
    INTEGER, ALLOCATABLE:: DISPLS( : )
    ! 
    ! grid number for each processor( MPI )
    ! NA: ghost point for MPI information exchange
    ! 		Here using DRP scheme, NA = 3.
    INTEGER:: XN1, YN2, ZN3
    ! 
    ! coordinate transformation
    REAL(KIND=8), ALLOCATABLE:: Jacobi( :, :, : )
    REAL(KIND=8), ALLOCATABLE:: KexiX( : ), EitaY( : ), TaoZ( : )
    ! 
    ! loop control variables
    INTEGER:: Tstep0, I0, J0, K0, MaxTimeStep, Loops, Flag
    REAL(KIND=8):: MaxResidual, MaxResidual0
    ! 
    ! output flag
    INTEGER:: IPX, IPY, IPZ
    INTEGER:: SizeN1, SizeN2, SizeN3
    INTEGER:: OUTPUTFlag, FLAG0, SIZE0, PmeanFlag, PrmsFlag
    INTEGER:: SizeNumber(3)
    INTEGER:: PmeanFlag1
    CHARACTER(LEN=80):: FilenameINPUT, FilenameOUTPUT, FilenameFile
    !  
    ! evaluate the parallel efficiency
    REAL(KIND=8):: STARTTIME, ENDTIME
    ! 
    !REAL(KIND=8):: FilterStartTime
    REAL(KIND=8), ALLOCATABLE:: Yposition( :, : )
    ! 
    ! restart parameter
    INTEGER:: Res_Loops
    ! 
    ! -------------added in Nov. 25, 2019,   for periodic BCs
    INTEGER:: NPML_Peri( 3, 2 )
    ! ------
    REAL(KIND=8):: Ek, Entropy
    REAL(KIND=8):: Ek_Sum, Entropy_Sum
    ! 
    ! -------------------------wall temperature boundary condition-----------
    ! INTEGER:: WallTempBCFlag
    ! WallTempBCFlag = 0: Dirichlet BCs;
    ! WallTempBCFlag = 1: Neumann BCs;
    ! 
    ! -------------------------evaluate the error of wall BCs----------------
    REAL(KIND=8):: L2_Error_Average
    ! 
    ! 
    REAL(KIND=8):: Temp(2), Temp1(2)
    !  
    INTEGER:: BodyForceOutFlag
    ! 
    INTEGER:: OutputloadingFlag(2)
    ! ------------------------------------CUT-OFF LINE-------------------------
    ! -------------------------------------------------------------------------
    !
    ! input the simulation parameter including the computational domain, MPI
    ! size ...... 
    ! CALL InitialModule(  )
    ! 
    ! 
    !!-------------------------------------STEP-1------------------------------
    ! -------------------------------------MPI PREPROCESSOR--------------------
    ! 
    !step- 1: initialize the parallel processors
    CALL MPI_INIT( IERR )
    CALL MPI_COMM_SIZE( MPI_COMM_WORLD, NUMPROCS, IERR )
    CALL MPI_COMM_RANK( MPI_COMM_WORLD, MYID, IERR )
    ! 
    ! 
    MYROOT = 0
    ! 
    !WRITE(*,*)MYID, 'Success!'
     ! input the simulation parameter including the computational domain, MPI
    ! size ......
    CALL InitialModule( MYID, MYROOT )
    !WRITE(*,*)MYID, 'Success!'
    !WRITE(*, *) MYID, NPX, NPY, NPZ
    ! 
    ! 
    ! for data exchanger
    ALLOCATE( RECVCNT( NUMPROCS ) )
    ALLOCATE( DISPLS( NUMPROCS ) )
    ! 
    ! 
    ! for series I/O: output results( MPI_ISEND, MPI_IRECV )
    IF ( MYID == MYROOT ) THEN
        SIZE0 = NUMPROCS - 1
        ALLOCATE( REQ( SIZE0 ) )
        ALLOCATE( STATUS( MPI_STATUS_SIZE, SIZE0 ) )
    ELSE
        SIZE0 = 1
        ALLOCATE( REQ( SIZE0:SIZE0 ) )
        ALLOCATE( STATUS( MPI_STATUS_SIZE, SIZE0 ) )
    END IF
    ! 
    ! 
    IF ( ( NPX * NPY * NPZ ) .NE. NUMPROCS ) THEN
        WRITE( *, * ) NPX, NPY, NPZ, 'The processor grid distribution doesn''t ', &
    &          ' conform with the processor summation! '
        STOP
        CALL MPI_FINALIZE( IERR )
    END IF
    ! 
    !!
    !get the current processor coordinates( PX, PY, PZ )
    PZ = MYID / ( NPX * NPY )
    PY = ( MYID - PZ * NPX * NPY ) / ( NPX )
    PX = ( MYID - PZ * NPX * NPY ) - PY * NPX
    ! 
    ! 
    !determine the topology relationship of all the processors
    !!
    MYLEFT = MYID - 1
    IF ( BCFLAG == 0 ) THEN
        IF ( MOD( MYID, NPX ) .EQ. 0 ) MYLEFT = MPI_PROC_NULL
    ELSEIF ( BCFLAG == 1 ) THEN
        IF ( MOD( MYID, NPX ) .EQ. 0 ) MYLEFT = NPX + MYLEFT
    END IF
    MYRIGHT = MYID + 1
    IF ( BCFLAG == 0 ) THEN
        IF ( MOD( MYRIGHT, NPX ) .EQ. 0 ) MYRIGHT = MPI_PROC_NULL
    ELSEIF ( BCFLAG == 1 ) THEN
        IF ( MOD( MYRIGHT, NPX ) .EQ. 0 ) MYRIGHT = MYRIGHT - NPX
    END IF
    MYREAR = MYID + NPX
    IF ( BCFlag == 0 ) THEN
        ! 
        ! free space for y direction: using PML boundary condition
        IF ( MYREAR .GE. ( PZ+1 )*NPX*NPY ) MYREAR = MPI_PROC_NULL
    ELSEIF( BCFlag == 1 ) THEN
        ! 
        ! periodic boundary condition for y direction: using periodic bc
        IF ( MYREAR .GE. ( PZ+1 )*NPX*NPY ) MYREAR = MYID - NPX * ( NPY - 1 )
    END IF
    MYFORWARD = MYID - NPX
    IF ( BCFlag == 0 ) THEN
        ! 
        ! free space for y -direction: using PML boundary condition
        IF ( MYFORWARD .LT. PZ*NPX*NPY ) MYFORWARD = MPI_PROC_NULL
    ELSEIF( BCFlag == 1 ) THEN
        ! 
        ! periodic boundary condition for y direction: using periodic bc
        IF ( MYFORWARD .LT. PZ*NPX*NPY ) MYFORWARD = MYID + NPX * ( NPY - 1 )
    END IF
    MYUPPER = MYID + NPX * NPY
    IF ( BCFLAG == 0 ) THEN
        IF ( MYUPPER .GE. NUMPROCS ) MYUPPER = MPI_PROC_NULL
    ELSEIF ( BCFLAG == 1 ) THEN
        IF ( MYUPPER .GE. NUMPROCS ) MYUPPER = MYUPPER - NUMPROCS
    END IF
    MYLOWER = MYID - NPX * NPY
    IF ( BCFLAG == 0 ) THEN
        IF( MYLOWER .LT. 0 ) MYLOWER = MPI_PROC_NULL 
    ELSEIF ( BCFLAG == 1 ) THEN
        IF( MYLOWER .LT. 0 ) MYLOWER = MYLOWER + NUMPROCS
    END IF
    ! ---------------------------------MPI PREPROCESSOR-------------------------
    ! ---------------------------------CUT-OFF LINE-----------------------------
    ! 
    !WRITE( *, FMT = "(7I3)" ) MYID, MYLEFT, MYRIGHT, MYFORWARD, MYREAR, MYLOWER, MYUPPER
    ! 
    !CALL MPI_BARRIER( MPI_COMM_WORLD, IERR )
    !!
    ! 
    ! ----------------------------------STEP-2----------------------------------
    ! ----------------------------------READ ORIGINAL MODEL DATA----------------
    OUTPUTFlag = 0
    PmeanFlag = 0
    PmeanFlag1 = 0
    PrmsFlag = 0
    FilenameINPUT = './INPUT/'
    FilenameOUTPUT = './OUTPUT/'
    FilenameFile = 'ShapesFile.dat'
    ! 
    ! 
    !step- 2: read the number for surface mesh
    IF ( MYID == MYROOT ) THEN
        OPEN( UNIT = 10, FILE = TRIM( TRIM( FilenameINPUT ) // TRIM( FilenameFile ) ) )
        READ( UNIT = 10, FMT = * ) NumBody
        CLOSE( UNIT = 10 )
    END IF
    ! 
    !!
    !bcast LN to every processor
    COUNT1 = 1
    CALL MPI_BCAST( NumBody, COUNT1, MPI_INTEGER, MYROOT, &
    &    MPI_COMM_WORLD, IERR )
    ALLOCATE( MultiLN( 1:NumBody ) )
    ! 
    ! 
    IF ( MYID == MYROOT ) THEN
        OPEN( UNIT = 10, FILE = TRIM( TRIM( FilenameINPUT ) // TRIM( FilenameFile ) ) )
        READ( UNIT=10, FMT = * ) NumBody
        READ( UNIT=10, FMT = * ) MultiLN( 1: NumBody )
        CLOSE( UNIT = 10 )
        LN = SUM( MultiLN )
    END IF
    COUNT1 = NumBody
    CALL MPI_BCAST( MultiLN, COUNT1, MPI_INTEGER, MYROOT, &
    &    MPI_COMM_WORLD, IERR )
    ! 
    ! 
    COUNT1 = 1
    CALL MPI_BCAST( LN, COUNT1, MPI_INTEGER, MYROOT, &
    &       MPI_COMM_WORLD, IERR )
    ! 
    !! 
    !WRITE(*,*) MYID, 'SUCCESS'
    ! 
    ! ---------------------------STEP-3-----------------------------------------------
    ! ---------------------------SOLID ARRAY SIZE DEFINITION--------------------------
    ! 
    !step- 3: allocate array variables for immersed body
    !   allocate variable related to the immersed boundary
    ALLOCATE( EffectRange( LN, 6 ), NormalVector( LN, 3 ) )
    ALLOCATE( Shapes( LN, 3 ), ShapesForP( LN, 3 ), Cell_S( LN ) )
    ALLOCATE( InitialShapes( LN, 3 ) )
    ALLOCATE( AV( LN, LN ), AT( LN, LN ) )
    ALLOCATE( LagForce( LN, 3 ), SurfUVW( LN, 3 ) )
    PLN = CEILING( 1.0d0 * LN / NUMPROCS )
    ALLOCATE( AProcessor( LN, PLN ) )
    ALLOCATE( ATotal( LN, NUMPROCS * PLN ) )
    ALLOCATE( DuDt0( LN ), DvDt0( LN ), DwDt0( LN ), DTdT0( LN ), DRouDt0( LN ) )
    ALLOCATE( DuDt( LN ), DvDt( LN ), DwDt( LN ), DTdT( LN ), DRouDt( LN )  )
    ALLOCATE( DuDt1( LN ), DvDt1( LN ), DwDt1( LN ), DTDt1( LN ), DRouDt1( LN ) )
    ALLOCATE( DuDtFlag( LN+1 ) )
    ALLOCATE( DuDtFlagRoot( LN ) )
    ALLOCATE( X01( LN, 1 ), X02( LN, 1 ), X03( LN, 1 ), X04( LN, 1 ) )
    ! 
    ! -----------------wall temperature boundary condition -----------------
    ALLOCATE( TangentialVector( LN, 3 ) )
    ALLOCATE( BodyCentre( NumBody, 3 ) )
    ALLOCATE( YPosition( NumBody, 3 ) )
    ALLOCATE( BodyCentrePosition( NumBody, 3 ) )
    ALLOCATE( BodyCentreVeloci( NumBody, 3 ) )
    ALLOCATE( BodyCentreAcce( NumBody, 3 ) )
    ! 
    ! 
    !!
    !read the original body coordinates for the main processor
    IF ( MYID == MYROOT ) THEN
        OPEN( UNIT = 10, FILE = TRIM( TRIM( FilenameINPUT ) // TRIM( FilenameFile ) ) )
        READ( UNIT = 10, FMT = * ) NumBody
        READ( UNIT = 10, FMT = * ) MultiLN( 1: NumBody )
        IF ( ModelFlag == 0 ) THEN
            DO I = 1, LN
                READ( UNIT = 10, FMT = * ) shapes( I, 1:3 )
            END DO
            Shapes( 1:LN , 3 ) = Centre( 3 )
        ELSE
            DO I = 1, LN
                READ( UNIT = 10, FMT = * ) shapes( I, 1:3 )
            END DO
        END IF
        CLOSE( UNIT = 10 )
    END IF
    ! 
    ! 
    ! bcast the shapes to other processors
    !bcast LN to every processor
    COUNT1 = LN*3
    CALL MPI_BCAST( Shapes(1,1), COUNT1, MPI_DOUBLE_PRECISION, MYROOT, &
    &       MPI_COMM_WORLD, IERR )
    ! ----------------------------------------------------------------------------------
    ! ----------------------------------------------------------------------------------
    ! ----------------------------STEP-4------------------------------------------------
    ! ----------------------------MODEL-------------------------------------------------
    !step- 4: input the analysis model
    ! 
    CALL Models( LRef, R_PML0, Shapes, BodyCentre, NormalVector, TangentialVector, Cell_S, LN, &
    &       NumBody, MultiLN, AOA, Centre, ModelFlag, 0 )
    InitialShapes = Shapes
    ! 
    ! 
    ! 
    !!----------------------------STEP-5-----------------------------------------------
    ! ---------------------------Mesh--------------------------------------------------
    ! 
    !step-5: get the size for Cartesian mesh
    CALL GetMeshSize( N1, N2, N3, StartEnd, CutPoints, DeltaXYZ, &
    &           R_PML, R_PML0, Shapes, Cell_S, LN, Centre, MaxX, MaxY, &
    &           MaxZ, SolidRelativePosi, q0, LRef, UniformGridRatio, NPML, &
    &           NPX, NPY, NPZ, ModelFlag, BCFlag, MYID, FilenameOUTPUT )
    ! 
    ALLOCATE( CL( NumBody ), Cd( NumBody ) )
    ! 
    !get the mesh number for each processor
    XN1 = ( 2 * NPML + N1 ) / NPX
    YN2 = ( 2 * NPML + N2 ) / NPY
    IF ( ModelFlag == 0 ) THEN
        ! 2-D Model
        ZN3 = 1
    ELSE
        ! 3-D Model
        ZN3 = ( 2 * NPML + N3 ) / NPZ
    END IF
    !!
    ! 
    ! ---------------------------------------STEP-6-------------------------------
    ! ---------------------------------Cartesian MESH ARRAY SIZE DEFINITION-------
    ! 
    !step-6: allocate array variables for Cartesian mesh of each processor
    ALLOCATE( MeshX( -(NPML-1) : (N1+NPML) ) )
    ALLOCATE( MeshY( -(NPML-1) : (N2+NPML) ) )
    ALLOCATE( MeshZ( -(NPML-1) : (N3+NPML) ) )
    ALLOCATE( DeltaX( -(NPML-1) : (N1+NPML-1) ) )
    ALLOCATE( DeltaY( -(NPML-1) : (N2+NPML-1) ) )
    ALLOCATE( DeltaZ( -(NPML-1) : (N3+NPML-1) ) )
    ALLOCATE( Jacobi( -(NPML-1):(N1+NPML), -(NPML-1):(N2+NPML), -(NPML-1):(N3+NPML) ) )
    ALLOCATE( KexiX( -(NPML-1) : ( N1 + NPML ) ) )
    ALLOCATE( EitaY( -(NPML-1) : ( N2 + NPML ) ) )
    ALLOCATE( TaoZ( -(NPML-1) : ( N3 + NPML ) ) )
    ALLOCATE( U( -(NA-1) : (XN1+NA), -(NA-1) : (YN2+NA), -(NA-1) : (ZN3+NA) ) )
    ALLOCATE( V( -(NA-1) : (XN1+NA), -(NA-1) : (YN2+NA), -(NA-1) : (ZN3+NA) ) )
    ALLOCATE( W( -(NA-1) : (XN1+NA), -(NA-1) : (YN2+NA), -(NA-1) : (ZN3+NA) ) )
    ALLOCATE( Omiga( -(NA-1) : (XN1+NA), -(NA-1) : (YN2+NA), -(NA-1) : (ZN3+NA) ) )
    ALLOCATE( P( -(NA-1) : (XN1+NA), -(NA-1) : (YN2+NA), -(NA-1) : (ZN3+NA) ) )
    ALLOCATE( ROU( -(NA-1) : (XN1+NA), -(NA-1) : (YN2+NA), -(NA-1) : (ZN3+NA) ) )
    ALLOCATE( Te( -(NA-1) : (XN1+NA), -(NA-1) : (YN2+NA), -(NA-1) : (ZN3+NA) ) )
    ALLOCATE( Mach( -(NA-1) : (XN1+NA), -(NA-1) : (YN2+NA), -(NA-1) : (ZN3+NA) ) )
    ALLOCATE( Q1( -(NA-1) : (XN1+NA), -(NA-1) : (YN2+NA), -(NA-1) : (ZN3+NA) ) )
    ALLOCATE( Q2( -(NA-1) : (XN1+NA), -(NA-1) : (YN2+NA), -(NA-1) : (ZN3+NA) ) )
    ALLOCATE( Q3( -(NA-1) : (XN1+NA), -(NA-1) : (YN2+NA), -(NA-1) : (ZN3+NA) ) )
    ALLOCATE( Q4( -(NA-1) : (XN1+NA), -(NA-1) : (YN2+NA), -(NA-1) : (ZN3+NA) ) )
    ALLOCATE( Q5( -(NA-1) : (XN1+NA), -(NA-1) : (YN2+NA), -(NA-1) : (ZN3+NA) ) )
    ALLOCATE( F1( -(NA-1) : (XN1+NA), -(NA-1) : (YN2+NA), -(NA-1) : (ZN3+NA) ) )
    ALLOCATE( F2( -(NA-1) : (XN1+NA), -(NA-1) : (YN2+NA), -(NA-1) : (ZN3+NA) ) )
    ALLOCATE( F3( -(NA-1) : (XN1+NA), -(NA-1) : (YN2+NA), -(NA-1) : (ZN3+NA) ) )
    ALLOCATE( F4( -(NA-1) : (XN1+NA), -(NA-1) : (YN2+NA), -(NA-1) : (ZN3+NA) ) )
    ALLOCATE( F5( -(NA-1) : (XN1+NA), -(NA-1) : (YN2+NA), -(NA-1) : (ZN3+NA) ) )
    ALLOCATE( G1( -(NA-1) : (XN1+NA), -(NA-1) : (YN2+NA), -(NA-1) : (ZN3+NA) ) )
    ALLOCATE( G2( -(NA-1) : (XN1+NA), -(NA-1) : (YN2+NA), -(NA-1) : (ZN3+NA) ) )
    ALLOCATE( G3( -(NA-1) : (XN1+NA), -(NA-1) : (YN2+NA), -(NA-1) : (ZN3+NA) ) )
    ALLOCATE( G4( -(NA-1) : (XN1+NA), -(NA-1) : (YN2+NA), -(NA-1) : (ZN3+NA) ) )
    ALLOCATE( G5( -(NA-1) : (XN1+NA), -(NA-1) : (YN2+NA), -(NA-1) : (ZN3+NA) ) )
    ALLOCATE( H1( -(NA-1) : (XN1+NA), -(NA-1) : (YN2+NA), -(NA-1) : (ZN3+NA) ) )
    ALLOCATE( H2( -(NA-1) : (XN1+NA), -(NA-1) : (YN2+NA), -(NA-1) : (ZN3+NA) ) )
    ALLOCATE( H3( -(NA-1) : (XN1+NA), -(NA-1) : (YN2+NA), -(NA-1) : (ZN3+NA) ) )
    ALLOCATE( H4( -(NA-1) : (XN1+NA), -(NA-1) : (YN2+NA), -(NA-1) : (ZN3+NA) ) )
    ALLOCATE( H5( -(NA-1) : (XN1+NA), -(NA-1) : (YN2+NA), -(NA-1) : (ZN3+NA) ) )
    !! 
    ALLOCATE( Shock_Sigma( -(NA-1) : (XN1+NA), -(NA-1) : (YN2+NA), -(NA-1) : (ZN3+NA) ) )
    ! 
    ! 
    ! total variables for outputting to file in main processor
    SizeN1 = N1 + 2 * NPML
    SizeN2 = N2 + 2 * NPML
    IF ( ModelFlag == 0 ) THEN
        ! 
        ! 2-D Model
        SizeN3 = 1
    ELSE
        ! 
        ! 3-D Model
        SizeN3 = N3 + 2 * NPML
    END IF
    ! ------
    ! ------
    ! create new datatype to output variables to file from each processor
    ! new datatype for sending process
    COUNT1 = YN2 * ZN3
    ALLOCATE( BLOCKLENS( COUNT1 ) )
    ALLOCATE( INDICES( COUNT1 ) )
    BLOCKLENS = XN1
    DO J = 1, ZN3
        DO I = 1, YN2
            INDICES( YN2*(J-1)+I ) = (J-1) * (XN1+2*NA) * (YN2+2*NA) + (I-1) * (XN1+2*NA)
        END DO
    END DO
    CALL MPI_TYPE_INDEXED( COUNT1, BLOCKLENS, INDICES, MPI_DOUBLE_PRECISION, HTYPE, IERR )
    CALL MPI_TYPE_COMMIT( HTYPE, IERR )
    ! 
    ! 
    ! new datattype for receving process
    DO J = 1, ZN3
        DO I = 1, YN2
            INDICES( YN2*(J-1)+I ) = (J-1) * SizeN1 * SizeN2 + (I-1) * SizeN1
        END DO
    END DO
    CALL MPI_TYPE_INDEXED( COUNT1, BLOCKLENS, INDICES, MPI_DOUBLE_PRECISION, VTYPE, IERR )
    CALL MPI_TYPE_COMMIT( VTYPE, IERR )
    ! --------
    ! --------
    !   
    ! 
    ! ------------
    ! ------------
    ! create new datatype for data exchange on the boundary for each processor
    ! 
    ! TYPE- 1: XTYPE ------for forward and rearward side data exchange
    DEALLOCATE( BLOCKLENS, INDICES )
    COUNT1 = NA * ZN3
    ALLOCATE( BLOCKLENS( COUNT1 ) )
    ALLOCATE( INDICES( COUNT1 ) )
    BLOCKLENS = XN1
    DO J = 1, ZN3
        DO I = 1, NA
            INDICES( NA*(J-1)+I ) = (J-1) * (XN1+2*NA) * (YN2+2*NA) + (I-1) * (XN1+2*NA)
        END DO 
    END DO
    CALL MPI_TYPE_INDEXED( COUNT1, BLOCKLENS, INDICES, MPI_DOUBLE_PRECISION, XTYPE, IERR )
    CALL MPI_TYPE_COMMIT( XTYPE, IERR )
    ! 
    ! TYPE- 2: YTYPE----------for left and right side data exchange
    DEALLOCATE( BLOCKLENS, INDICES )
    COUNT1 = YN2 * ZN3
    ALLOCATE( BLOCKLENS( COUNT1 ) )
    ALLOCATE( INDICES( COUNT1 ) )
    BLOCKLENS = NA
    DO J = 1, ZN3
        DO I = 1, YN2
            INDICES( YN2*(J-1)+I ) = (J-1) * (XN1+2*NA) * (YN2+2*NA) + (I-1) * (XN1+2*NA)
        END DO 
    END DO
    CALL MPI_TYPE_INDEXED( COUNT1, BLOCKLENS, INDICES, MPI_DOUBLE_PRECISION, YTYPE, IERR )
    CALL MPI_TYPE_COMMIT( YTYPE, IERR )
    ! 
    ! another type for Y
    DEALLOCATE( BLOCKLENS, INDICES )
    COUNT1 = ( YN2 + 2 * NA ) * ( ZN3 + 2 * NA )
    ALLOCATE( BLOCKLENS( COUNT1 ) )
    ALLOCATE( INDICES( COUNT1 ) )
    BLOCKLENS = NA
    DO J = 1, ZN3 + 2 * NA
        DO I = 1, YN2 + 2 * NA
            INDICES( (YN2+2*NA)*(J-1) + I ) = (J-1) * (XN1+2*NA) * (YN2+2*NA) + (I-1) * (XN1+2*NA)
        END DO
    END DO
    CALL MPI_TYPE_INDEXED( COUNT1, BLOCKLENS, INDICES, MPI_DOUBLE_PRECISION, YTYPE1, IERR )
    CALL MPI_TYPE_COMMIT( YTYPE1, IERR )
    ! 
    ! 
    ! TYPE- 3: ZTYPE ---------for upper and lower side data exchange
    DEALLOCATE( BLOCKLENS, INDICES )
    COUNT1 = YN2 * NA
    ALLOCATE( BLOCKLENS( COUNT1 ) )
    ALLOCATE( INDICES( COUNT1 ) )
    BLOCKLENS = XN1
    DO J = 1, NA
        DO I = 1, YN2
            INDICES( YN2*(J-1)+I ) = (J-1) * (XN1+2*NA) * (YN2+2*NA) + (I-1) * (XN1+2*NA)
        END DO
    END DO
    CALL MPI_TYPE_INDEXED( COUNT1, BLOCKLENS, INDICES, MPI_DOUBLE_PRECISION, ZTYPE, IERR )
    CALL MPI_TYPE_COMMIT( ZTYPE, IERR )
    ! -------------
    !step- 7: grid generation
    CALL Mesh( MeshX, MeshY, MeshZ, DeltaX, DeltaY, DeltaZ, N1, N2, N3, &
    &       NPML, StartEnd, CutPoints, DeltaXYZ, q0, ModelFlag, MYID, FilenameOUTPUT )
    ! 
    ! 
    ! mesh smoothing
    CALL MeshSmoothing( MeshX, MeshY, MeshZ, DeltaX, DeltaY, DeltaZ, N1, &
    &                     N2, N3, NPML, StartEnd, CutPoints, DeltaXYZ, q0, &
    &                     LRef, UniformGridRatio, ModelFlag, MYID, MYROOT, &
    &                     FilenameOUTPUT, 2 )
    ! 
    ! 
    !
    !step- 7.1: coordinate transformation/ get Jacobi matrix
    CALL TransformCoordinate( Jacobi, KexiX, EitaY, TaoZ, MeshX, MeshY, MeshZ, &
    &                   N1, N2, N3, NPML, DeltaXYZ, ModelFlag, CoorTransScheme )
    ! 
    ! 
    ! 
    ! ---------------------------------------------STEP-8-------------------------------
    ! ---------------------------------------------BACKGROUND FLOW----------------------
    ! 
    !step- 8: get the background flow field1
    CALL GetMeanFlow( U, V, W, P, ROU, NA, XN1, YN2, ZN3, MeshX, MeshY, MeshZ, &
    &       Centre, StartEnd, NPML, N1, N2, N3, PX, PY, PZ, NPX, NPY, NPZ, ModelFLag, &
    &       MYID, MYROOT, InitialFlag )
    !! 
    ! 
    ! 
    ! --------------------------------------STEP-9--------------------------------------
    ! -------------------------------------INITIALIZE FLOW FIELD------------------------
    ! 
    !step- 9: initialize acoutic field
    CALL InitializeAcousticField( U, V, W, P, ROU, NA, XN1, YN2, &
    &           ZN3, MeshX, MeshY, MeshZ, NPML, N1, N2, N3, PX, PY, PZ, &
    &                NPX, NPY, NPZ, ModelFlag, SourceFlag )
    ! 
    ! 
    !  
    ! 
    !read the original body coordinates for the main processor
    ! 
    K0 = 0
    DeltaT = DeltaXYZ * CFL / ( 1.0d0 + Ma )
    !Tstep0 = floor( 1d0 / DeltaT )
    !Tstep0 = 1
    Tstep0 = CEILING( TimeOutInterval / DeltaT )
    ! 
    MaxTimeStep = CEILING( TotalTime / DeltaT )
    ! 
    SurfUVW = 0.0d0
    ! FilterStartTime = 0.1d0*MaxX
    Yposition = 0.0d0
    ! 
    ! 
    ! -------------------------------------newly added in Nov. 30, 2019-----------
    NPML_Peri = NPML
    ! 
    ! revise the streamwise direction and vertical direction. 
    NPML_Peri( 1, 2 ) = N1 + NPML - N2
    NPML_Peri( 3, 2 ) = N3 + NPML - N2
    ! 
    IF ( MYID == MYROOT ) THEN
        DO I = 1, 3
            WRITE( *, FMT = "(2I4)" ) NPML_Peri( I, 1 : 2 )
        END DO
    END IF
    ! -----------------------------------------------------------------------------
    ! 
    Q1 = 0.0d0
    Q2 = 0.0d0
    Q3 = 0.0d0
    Q4 = 0.0d0
    Q5 = 0.0d0
    ! 
    Res_Loops = 1
    ! 
    ! 
    IF ( MYID == MYROOT ) THEN
        OPEN( UNIT = 15, FILE = TRIM( TRIM( FilenameOUTPUT ) // TRIM( 'TKE.dat' ) ) )
        WRITE( UNIT = 15, FMT = * ) "Time,   Kinetic Energy,   Enstrophy"
    END IF
    ! 
    ! output systematic parameters
    !CALL PrintParameterSet( MYID, MYROOT, NUMPROCS, DeltaXYZ, DeltaT, Tstep0 )
    ! 
    IF ( MYID==MYROOT ) THEN
        CALL CPU_TIME( STARTTIME )
    END IF
    ! 
    ! the device number 16~ can be used.
    DO Loops = Res_Loops, MaxTimeStep
        ! -------------------------------
        ! ----------------body moving-----------------
        ! 
        CALL GetConservativeVariables( Q1, Q2, Q3, Q4, Q5, F1, F2, F3, F4, F5, &
    &           G1, G2, G3, G4, G5, H1, H2, H3, H4, H5, &
    &           U, V, W, P, ROU, NA, XN1, YN2, ZN3, Jacobi, KexiX, EitaY, TaoZ, &
    &           MeshX, MeshY, MeshZ, Beita, SigmaX1, SigmaX2, SigmaY1, SigmaY2, &
    &           SigmaZ1, SigmaZ2, NPML, N1, N2, N3, MYID, PX, PY, PZ, NPX, NPY, NPZ, &
    &           Loops, DeltaT, Loops*DeltaT, DeltaXYZ, ModelFlag, BCFlag, SourceFlag, VSADFLAG, NSFLAG )
        ! 
        ! 
            ! add diffusion term to the governing equations: viscosity
            ! term------  N-S equation
        IF ( NSFLAG == 1 ) THEN
            ! Solve N-S equations
            CALL ViscosityStressTerm( F2, F3, F4, F5, G2, G3, G4, G5, H2, H3, H4, H5, U, V, W, &
    &                  P, ROU, NA, XN1, YN2, ZN3, Jacobi, KexiX, EitaY, TaoZ, NPML, NPML_Peri, N1, N2, N3, &
    &                       DeltaXYZ, MYID, PX, PY, PZ, MYLEFT, MYRIGHT, MYFORWARD, MYREAR, &
    &                           MYUPPER, MYLOWER, XTYPE, YTYPE1, ZTYPE, NPX, NPY, NPZ, &
    &                               ModelFlag, BCFLAG, DiffusiveScheme )
        END IF
        ! 
        !    CALL MPI_BARRIER( MPI_COMM_WORLD, IERR )
        !END IF
        Flag = MOD( Loops, 2 )
        !Flag = 0
        !
        ! 
        CALL LDDRK( Q1, Q2, Q3, Q4, Q5, F1, F2, F3, F4, F5, G1, G2, G3, G4, G5, &
    &           H1, H2, H3, H4, H5, U, V, W, P, ROU, Shock_Sigma, Jacobi, KexiX, EitaY, &
    &           TaoZ, MeshX, MeshY, MeshZ, SigmaX1, SigmaX2, SigmaY1, SigmaY2, SigmaZ1, &
    &           SigmaZ2, Beita, NA, XN1, YN2, ZN3, N1, N2, N3, NPML, NPML_Peri, DeltaT, DeltaXYZ, &
    &           Shapes, LN, PLN, Cell_S, NormalVector, InitialShapes, SurfUVW, StartEnd, &
    &           CutPoints, NumBody, MultiLN, PX, PY, PZ, NPX, NPY, NPZ, MYID, &
    &           MYROOT, NUMPROCS, MYLEFT, MYRIGHT, MYFORWARD, MYREAR, MYUPPER, MYLOWER, &
    &           XTYPE, YTYPE, YTYPE1, ZTYPE, Flag, Loops, VSADFLAG, ModelFlag, &
    &           BCFlag, SourceFlag, NSFLAG, CoupledFlag, FilterStartTime, &
    &           CoupledFilterFlag, ShockCapturingFilterFlag, ShockSensorFlag, &
    &           ShockFilterStencilFlag, BodyForceSolverFlag, DistriFlag, &
    &           ConvectiveScheme, DiffusiveScheme )
        !!
        ! 
        CALL ExchangeInterfaceDataNew( Q1, NA, XN1, YN2, ZN3, NPML, NPML_Peri, XTYPE, YTYPE1, &
    &               ZTYPE, MYLEFT, MYRIGHT, MYFORWARD, MYREAR, MYUPPER, MYLOWER, &
    &                   PX, PY, PZ, NPX, NPY, NPZ, ModelFlag, BCFlag, 1 )
        ! 
        CALL ExchangeInterfaceDataNew( Q2, NA, XN1, YN2, ZN3, NPML, NPML_Peri, XTYPE, YTYPE1, &
    &               ZTYPE, MYLEFT, MYRIGHT, MYFORWARD, MYREAR, MYUPPER, MYLOWER, &
    &                   PX, PY, PZ, NPX, NPY, NPZ, ModelFlag, BCFlag, 1 )
        ! 
        CALL ExchangeInterfaceDataNew( Q3, NA, XN1, YN2, ZN3, NPML, NPML_Peri, XTYPE, YTYPE1, &
    &               ZTYPE, MYLEFT, MYRIGHT, MYFORWARD, MYREAR, MYUPPER, MYLOWER, &
    &                   PX, PY, PZ, NPX, NPY, NPZ, ModelFlag, BCFlag, 1 )
        ! 
        IF ( ModelFlag == 1 ) THEN
            CALL ExchangeInterfaceDataNew( Q4, NA, XN1, YN2, ZN3, NPML, NPML_Peri, XTYPE, YTYPE1, &
    &                   ZTYPE, MYLEFT, MYRIGHT, MYFORWARD, MYREAR, MYUPPER, MYLOWER, &
    &                       PX, PY, PZ, NPX, NPY, NPZ, ModelFlag, BCFlag, 1 )
        END IF
        ! 
        CALL ExchangeInterfaceDataNew( Q5, NA, XN1, YN2, ZN3, NPML, NPML_Peri, XTYPE, YTYPE1, &
    &               ZTYPE, MYLEFT, MYRIGHT, MYFORWARD, MYREAR, MYUPPER, MYLOWER, &
    &                   PX, PY, PZ, NPX, NPY, NPZ, ModelFlag, BCFlag, 1 )
        !!
        !------------------------------------------------------------------------------------
        ! -----------------------------------------------------------------------------------
        ! 
        IF ( ( FilterFLAG == 1 ) .AND. MOD(Loops, 1) == 1 ) THEN
            ! get the revised conservative variables and then filter the
            ! conservative variables
            ! 
            ! x- direction filter
            ! excute the explicit filtering for the conservative variables.
            CALL Filter( Q1, Q2, Q3, Q4, Q5, NA, XN1, YN2, ZN3, ModelFlag, DeltaT*Loops, &
    &                                FilterStartTime, NPML, PX, PY, PZ, NPX, NPY, NPZ, 0 )
            ! 
            ! exchange the virtual mesh variables. 
            CALL ExchangeInterfaceDataNew( Q1, NA, XN1, YN2, ZN3, NPML, NPML_Peri, XTYPE, YTYPE1, &
    &                   ZTYPE, MYLEFT, MYRIGHT, MYFORWARD, MYREAR, MYUPPER, MYLOWER, &
    &                       PX, PY, PZ, NPX, NPY, NPZ, ModelFlag, BCFlag, 1 )
            !
            CALL ExchangeInterfaceDataNew( Q2, NA, XN1, YN2, ZN3, NPML, NPML_Peri, XTYPE, YTYPE1, &
    &                   ZTYPE, MYLEFT, MYRIGHT, MYFORWARD, MYREAR, MYUPPER, MYLOWER, &
    &                       PX, PY, PZ, NPX, NPY, NPZ, ModelFlag, BCFlag, 1 )
            CALL MPI_BARRIER( MPI_COMM_WORLD, IERR )
            ! 
            CALL ExchangeInterfaceDataNew( Q3, NA, XN1, YN2, ZN3, NPML, NPML_Peri, XTYPE, YTYPE1, &
    &                   ZTYPE, MYLEFT, MYRIGHT, MYFORWARD, MYREAR, MYUPPER, MYLOWER, &
    &                       PX, PY, PZ, NPX, NPY, NPZ, ModelFlag, BCFlag, 1 )
            !  
            IF ( ModelFlag == 1 ) THEN
                CALL ExchangeInterfaceDataNew( Q4, NA, XN1, YN2, ZN3, NPML, NPML_Peri, XTYPE, YTYPE1, &
    &                       ZTYPE, MYLEFT, MYRIGHT, MYFORWARD, MYREAR, MYUPPER, MYLOWER, &
    &                            PX, PY, PZ, NPX, NPY, NPZ, ModelFlag, BCFlag, 1 )
            END IF
            ! 
            CALL ExchangeInterfaceDataNew( Q5, NA, XN1, YN2, ZN3, NPML, NPML_Peri, XTYPE, YTYPE1, &
    &                   ZTYPE, MYLEFT, MYRIGHT, MYFORWARD, MYREAR, MYUPPER, MYLOWER, &
    &                       PX, PY, PZ, NPX, NPY, NPZ, ModelFlag, BCFlag, 1 )
            ! 
            ! y-direction filter
            ! excute the explicit filtering for the conservative variables.
            CALL Filter( Q1, Q2, Q3, Q4, Q5, NA, XN1, YN2, ZN3, ModelFlag, DeltaT*Loops, &
    &                                FilterStartTime, NPML, PX, PY, PZ, NPX, NPY, NPZ, 1 )
            ! 
            ! exchange the virtual mesh variables. 
            CALL ExchangeInterfaceDataNew( Q1, NA, XN1, YN2, ZN3, NPML, NPML_Peri, XTYPE, YTYPE1, &
    &                   ZTYPE, MYLEFT, MYRIGHT, MYFORWARD, MYREAR, MYUPPER, MYLOWER, &
    &                       PX, PY, PZ, NPX, NPY, NPZ, ModelFlag, BCFlag, 1 )
            !
            CALL ExchangeInterfaceDataNew( Q2, NA, XN1, YN2, ZN3, NPML, NPML_Peri, XTYPE, YTYPE1, &
    &                   ZTYPE, MYLEFT, MYRIGHT, MYFORWARD, MYREAR, MYUPPER, MYLOWER, &
    &                       PX, PY, PZ, NPX, NPY, NPZ, ModelFlag, BCFlag, 1 )
            ! 
            CALL ExchangeInterfaceDataNew( Q3, NA, XN1, YN2, ZN3, NPML, NPML_Peri, XTYPE, YTYPE1, &
    &                   ZTYPE, MYLEFT, MYRIGHT, MYFORWARD, MYREAR, MYUPPER, MYLOWER, &
    &                       PX, PY, PZ, NPX, NPY, NPZ, ModelFlag, BCFlag, 1 )
            ! 
            CALL MPI_BARRIER( MPI_COMM_WORLD, IERR )
            IF ( ModelFlag == 1 ) THEN
                CALL ExchangeInterfaceDataNew( Q4, NA, XN1, YN2, ZN3, NPML, NPML_Peri, XTYPE, YTYPE1, &
    &                       ZTYPE, MYLEFT, MYRIGHT, MYFORWARD, MYREAR, MYUPPER, MYLOWER, &
    &                           PX, PY, PZ, NPX, NPY, NPZ, ModelFlag, BCFlag, 1 )
            END IF
            ! 
            CALL ExchangeInterfaceDataNew( Q5, NA, XN1, YN2, ZN3, NPML, NPML_Peri, XTYPE, YTYPE1, &
    &                   ZTYPE, MYLEFT, MYRIGHT, MYFORWARD, MYREAR, MYUPPER, MYLOWER, &
    &                       PX, PY, PZ, NPX, NPY, NPZ, ModelFlag, BCFlag, 1 )
            ! 
            ! 
            IF ( ModelFlag == 1 ) THEN
            ! z direction
            ! added in 2 Dec. 2019
            CALL Filter( Q1, Q2, Q3, Q4, Q5, NA, XN1, YN2, ZN3, ModelFlag, DeltaT*Loops, &
    &                                FilterStartTime, NPML, PX, PY, PZ, NPX, NPY, NPZ, 2 )
            ! 
            ! exchange the virtual mesh variables. 
            CALL ExchangeInterfaceDataNew( Q1, NA, XN1, YN2, ZN3, NPML, NPML_Peri, XTYPE, YTYPE1, &
    &                   ZTYPE, MYLEFT, MYRIGHT, MYFORWARD, MYREAR, MYUPPER, MYLOWER, &
    &                       PX, PY, PZ, NPX, NPY, NPZ, ModelFlag, BCFlag, 1 )
            !
            CALL ExchangeInterfaceDataNew( Q2, NA, XN1, YN2, ZN3, NPML, NPML_Peri, XTYPE, YTYPE1, &
    &                   ZTYPE, MYLEFT, MYRIGHT, MYFORWARD, MYREAR, MYUPPER, MYLOWER, &
    &                       PX, PY, PZ, NPX, NPY, NPZ, ModelFlag, BCFlag, 1 )
            ! 
            CALL ExchangeInterfaceDataNew( Q3, NA, XN1, YN2, ZN3, NPML, NPML_Peri, XTYPE, YTYPE1, &
    &                   ZTYPE, MYLEFT, MYRIGHT, MYFORWARD, MYREAR, MYUPPER, MYLOWER, &
    &                       PX, PY, PZ, NPX, NPY, NPZ, ModelFlag, BCFlag, 1 )
            ! 
            CALL MPI_BARRIER( MPI_COMM_WORLD, IERR )
            IF ( ModelFlag == 1 ) THEN
                CALL ExchangeInterfaceDataNew( Q4, NA, XN1, YN2, ZN3, NPML, NPML_Peri, XTYPE, YTYPE1, &
    &                       ZTYPE, MYLEFT, MYRIGHT, MYFORWARD, MYREAR, MYUPPER, MYLOWER, &
    &                           PX, PY, PZ, NPX, NPY, NPZ, ModelFlag, BCFlag, 1 )
            END IF
            ! 
            CALL ExchangeInterfaceDataNew( Q5, NA, XN1, YN2, ZN3, NPML, NPML_Peri, XTYPE, YTYPE1, &
    &                   ZTYPE, MYLEFT, MYRIGHT, MYFORWARD, MYREAR, MYUPPER, MYLOWER, &
    &                       PX, PY, PZ, NPX, NPY, NPZ, ModelFlag, BCFlag, 1 )
            ! 
            END IF
            ! 
        END IF
        CALL GetOriginalVariables( U, V, W, P, ROU, Q1, Q2, Q3, Q4, Q5, NA, &
    &               XN1, YN2, ZN3, Jacobi, KexiX, EitaY, TaoZ, MeshX, MeshY, MeshZ, N1, N2, N3, &
!    &                   PX, PY, PZ, NPX, NPY, NPZ, NPML, DeltaT*(Loops+1.0), ModelFlag, 0 )
    &                   PX, PY, PZ, NPX, NPY, NPZ, NPML, DeltaT*Loops, ModelFlag, 0 )
            ! 
            ! 
        ! 
        CALL ObtainKineticEnergy( Ek, Entropy, U, V, W, ROU, XN1, YN2, ZN3, NA, NPML, NPML_Peri, &
    &                           MaxX, MaxY, MaxZ, DeltaXYZ, PX, PY, PZ, NPX, NPY, NPZ, ModelFlag )
        ! 
        ! 
        IF ( MOD( Loops, 5 ) == 0 ) THEN
            ! 
            CALL MPI_REDUCE(Ek, Ek_Sum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MYROOT, MPI_COMM_WORLD, IERR)
            ! 
            CALL MPI_REDUCE(Entropy, Entropy_Sum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MYROOT, MPI_COMM_WORLD, IERR)
            ! 
            IF ( MYID == MYROOT ) THEN
                WRITE( UNIT = 15, FMT = '(3F18.12)' ) Loops*DeltaT*Ma, Ek_Sum/Ma**2, Entropy_Sum/Ma**2
            END IF
        END IF
        ! 
        ! 
        ! 
        IF ( ( ( MOD( Loops, Tstep0 )==0 ) .AND. Loops <= floor(0.5d0*TotalTime/DeltaT) ) .OR. &
    &        ( ( MOD( Loops, Tstep0 )==0 ) .AND. Loops >= floor(0.9d0*TotalTime/DeltaT) ) ) THEN
           !IF ( 1 == 1 ) THEN
            OUTPUTFlag = OUTPUTFlag + 1
            ! output U
            SizeNumber( 1 ) = 1
            SizeNumber( 2 ) = 1
            SizeNumber( 3 ) = 1
            ! 
            ! 
            ! get vorticity omiga
            IF ( ModelFlag == 0 ) THEN
                ! 2-D 
                DO J = 1, YN2
                    DO I = 1, XN1
                        Omiga( I, J, 1 ) = KexiX(I+PX*XN1-NPML) * ( V(I+1, J, 1) - V(I-1, J, 1) ) &
    &                                      / (2.0d0 * DeltaXYZ) - EitaY(J+PY*YN2-NPML) * &
    &                                      ( U(I, J+1, 1)-U(I, J-1, 1) ) / ( 2.0d0 * DeltaXYZ )
                    END DO
                END DO
            ELSE
                ! 3-D
                DO K = 1, ZN3
                DO J = 1, YN2
                    DO I = 1, XN1
                        Omiga( I, J, K ) = ( V(I+1, J, K) - V(I-1, J, K) ) / (2.0d0 * DeltaXYZ) -  &
    &                                      ( U(I, J+1, K)-U(I, J-1, K) ) / ( 2.0d0 * DeltaXYZ )
                    END DO
                END DO
                END DO
            END IF
            ! 
    !        CALL MPI_BARRIER( MPI_COMM_WORLD, IERR )
        END IF
        ! 
    END DO
    ! 
    IF ( MYID==MYROOT ) THEN
        CALL CPU_TIME( ENDTIME )
        WRITE(*,FMT='(f24.14)' )ENDTIME-STARTTIME
    END IF
    ! 
    IF ( MYID == MYROOT ) THEN
        CLOSE( UNIT = 15 )
    END IF
    ! 
    ! -----------------------------------------------------------------------
    !--------------- outer loops end! ---------                             -
    !                                                                       -
    ! -----------------------------------------------------------------------
    ! 
    !!
    !step- 11: return and exit
    CALL MPI_FINALIZE( IERR )
END PROGRAM AcousticScatteringParallelComputing


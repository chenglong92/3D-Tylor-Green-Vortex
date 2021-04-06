!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    INSTRUCTION (High-order CAA Parallel Solver)
! -----------------------------------------------------------------------------
! Language: FORTRAN 90 + MPICH
! -----------------------------------------------------------------------------
! Parallel Mode: 3D domain decomposition. 
! Grid Type: Structured Cartesian grid. 
! Numerical Schemes: 
!       2N Storage 5/6 stages LDDRK (Hu); 
!       7-points cetral DRP(Tam); 
!       9/11-points Optimized Filtering(Bogey); 
!       Ghost Layer: >4; 
! -----------------------------------------------------------------------------
! Files: 
!    "INPUT": Geometry file -- determine the minimum grid size
!    "OUTPUT": output file. 
!    "ParameterSetting": input file, to define the computational domain, 
!                        processor distribution ....
!    "SRC_Coupled_NewSpringVersion": Source Code. 
! -----------------------------------------------------------------------------
! Developer: Long Cheng
! E-mail:  LC759@cam.ac.uk   &  Long.Cheng@buaa.edu.cn
! Tel:     (0044) 07938458738
! JDB 02-1, Fluid Dynamics Group, Engineering Dept., Cambridge University (CUED)
! ------------------------------------------------------------------------------


Methods to set the processor number: 
    
    1. Compile: (default: mpiifort)

          mpiifort -o NS.exe SRC_Coupled_NewSpringVersion/Main.f90   

                        or

          use Makefile: make

    2. Revise the processor grid using the file:

          "~/ParameterSetting/InitialConstantsModule.txt"
       
       Here total processor number: NUMPROCS = NPX*NPY*NPZ
              It is better to choose from the following numbers: 
                   NPX = NPY = NPZ = 2    (NUMPROCS = 8); 
                   NPX = NPY = NPZ = 3    (NUMPROCS = 27); 
                   NPX = NPY = NPZ = 4    (NUMPROCS = 64); 
                   NPX = NPY = NPZ = 5    (NUMPROCS = 125);
                   NPX = NPY = NPZ = 6    (NUMPROCS = 216); 
                   NPX = NPY = NPZ = 7    (NUMPROCS = 343);
                   NPX = NPY = NPZ = 8    (NUMPROCS = 512); 
                   NPX = NPY = NPZ = 9    (NUMPROCS = 729);
                   NPX = NPY = NPZ = 10   (NUMPROCS = 1000); 
                   NPX = NPY = NPZ = 11   (NUMPROCS = 1331);
                   NPX = NPY = NPZ = 12   (NUMPROCS = 1728); 
                   NPX = NPY = NPZ = 13   (NUMPROCS = 2197);
                   NPX = NPY = NPZ = 14   (NUMPROCS = 2744); 
                   NPX = NPY = NPZ = 15   (NUMPROCS = 3375);
                   NPX = NPY = NPZ = 16   (NUMPROCS = 4096); 
                   NPX = NPY = NPZ = 17   (NUMPROCS = 4913);
                   NPX = NPY = NPZ = 18   (NUMPROCS = 5832); 
                   NPX = NPY = NPZ = 19   (NUMPROCS = 6859);
                   NPX = NPY = NPZ = 20   (NUMPROCS = 8000); 
                   NPX = NPY = NPZ = 21   (NUMPROCS = 9261);
                   NPX = NPY = NPZ = 22   (NUMPROCS = 10648); 
                   NPX = NPY = NPZ = 23   (NUMPROCS = 12167);
                   NPX = NPY = NPZ = 24   (NUMPROCS = 13824); 
                   NPX = NPY = NPZ = 25   (NUMPROCS = 15625);


      --- !!!! due to the ghost layer, the NPX,NPY,NPZ <= 25 if grid point = 256^3. ---


!**********************************************************************************
       "InitialConstantsModule.txt": 
                                 ......

                                 ......

           ! -------------------------PROCESSOR GRID----------------
           ! MPI processor grid:                                   -
           !   NPX: the number for the x direction;                -
           !   NPY: the number for the y direction;                -
           !   NPZ: the number for the z direction;                -
           ! -------------------------------------------------------
           ! ---                NPX, NPY, NPZ                      -
           ! -------------------------------------------------------
                                 4    4    4                (--------tO be revised.)
!***********************************************************************************

                                                               5 Dec. 2019
                                                               Written by Long Cheng

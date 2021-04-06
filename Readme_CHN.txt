!!!!!!!!!!!!!!!!!!中文版!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    程序使用说明书 (High-order CAA Parallel Solver) 
!!!!! 该程序主要是利用高阶CFD/CAA数值格式进行三维Taylor-Green Vortex(TGV)可压缩自由衰减湍流模拟，
从而测试数值格式以及边界条件在精细刻画湍流方面的正确性及能力，同时测试程序的并行效率(强扩展性)。


! -----------------------------------------------------------------------------
 Language(语言): FORTRAN 90 + MPICH
 Platform(测试平台): Linux集群(如FAEL Cluster以及Tianhe-2)


 -----------------------------------------------------------------------------
-------------------------------! 算法说明!----------------------------- 
 Parallel Mode(并行策略): 3D domain decomposition(三维域分解，XYZ方向分别进行MPI计算网格分块). 
 Grid Type(网格类型): Structured Cartesian grid(结构化Cartesian全正交网格). 
 数值方法: Finite Difference Method(FDM, 有限差分) 
 Numerical Schemes(数值格式): 
       时间离散: 2N Storage 5/6 stages LDDRK (Hu); ------ Hu Fangqiang以及Stanescu发展的2N低存储5/6交替的LDDRK格式 
       空间离散: 7-points cetral DRP(Tam);  ----- Tam发展的7点中心DRP格式
       高阶滤波格式: 9/11-points Optimized Filtering(Bogey);   ---- 9点或者11点优化高阶滤波格式
       边界处理鬼点层：Ghost Layer: >4;  -----用于MPI分区后的每个时间推进步相邻MPI分区数据交换
 Boundary Conditions(边界条件): XYZ方向均是周期性边界条件


! -------------------------------------------------------------------------------------------------------------------------------------------------
! 程序文件组成说明: 
    "./INPUT": 该文件夹中主要是提供几何文件表面网格点坐标， 这里由于不包含几何，因此该文件夹中ShapesFile.dat文件决定最小网格尺寸； 
    "./OUTPUT": 计算结果输出文件夹 
    "./ParameterSetting": 该文件夹中用于计算参数输入，主要包含InitialConstantsModule.txt, InitialBackgroundFlowModule.txt,  InitialAcousticModule.txt, InitialBodyMoingModule.txt四个文件以及./InitializedFlowField和./RestartFlowField两个文件夹，其主要功能简述如下： 
———————————————————————————————————————————————————————————————————
                        文件                             |                                                                         功能
———————————————————————————————————————————————————————————————————
     InitialConstantsModule.txt              |   输入问题维数，求解方程类型，计算域， 时间步长，数值格式选择，滤波格式及滤波方法，MPI并行核数分配等
 InitialBackgroundFlowModule.txt       |   输入流动条件，如特征马赫数，雷诺数等
     InitialAcousticModule.txt                |   小扰动初始化参数(这里仅作程序测试用，具体计算时没有用，流场初始化在子程序GetMeanFlow_sub.f90中自动进行)
   InitialBodyMoingModule.txt             |   指定物体运动参数(这里没有物体扰流，因此该文件没有用)
        ./InitializedFlowField                   |    当选择从外部文件导入初始化流场时读取该文件夹中数据文件(此处未用)
           ./RestartFlowField                    |    当需要设置流场计算重启时从该文件夹中读取(此处为用)
———————————————————————————————————————————————————————————————————
    -
    "SRC_Coupled_NewSpringVersion": 源程序，主要包含如下子程序，相应说明如下： 
———————————————————————————————————————————————————————————————————
                        文件                             |                                                                         功能
———————————————————————————————————————————————————————————————————
                    Main.f90                          |   主程序，基本按照顺序流程执行，基本分为网格生成，初始化流场，时间推进求解，滤波，文件输出等步骤
         Module_Constants_sub.f90          |   该Module定义计算参数变量或者常量，方便InitialModule_sub.f90读取计算参数
          NumericalProbe_sub.f90            |    该程序设置数值探针位置，从而对流场进行采样得到某些散点初的流场信息
      AbsorbingCoefficient_sub.f90         |   该程序设置外边界PML区域内的吸收系数(该处由于在XYZ方向上均采用了周期性边界条件，因此无自由边界区域)
  ExchangeInterfaceDataNew_sub.f90   |   该程序交换每一个时间层或者亚层推进后相邻MPI网格重叠区域流场信息
   GetConservativeVariables_sub.f90      |   获得计算网格下的守恒变量以及flux信息
          GetMeanFlow_sub.f90                |   初始化流场
          GetMeshSize_sub.f90                 |    得到Cartesian计算网格尺寸
     GetOriginalVariables_sub.f90          |    由计算得到的守恒变量反推得到原始速度，压力等原始流场变量
    InitializeAcousticField_sub.f90         |     初始化扰动场信息(未用)
             LDDRK_sub.f90                      |     LDDRK时间推进得到下一时间步守恒变量
              Mesh_sub.f90                       |     生成Cartesian网格坐标点
      MeshSmoothing_sub.f90              |     光顺网格，仅用于当网格为非均匀情况
            Models_sub.f90                      |     读取几何网格点文件，生成计算模型
  TransformCoordinate_sub.f90          |      坐标变换，得到Jacobi矩阵等坐标变换信息
    OutputToTextFile_sub.f90              |      将流场文件写到硬盘中
             Filter_sub.f90                        |      高阶滤波
  ReviseConserVariable_sub.f90          |      修正守恒变量(仅用于当流场中存在扰流物体时)
    PrintParameterSet_sub.f90             |     打印计算参数信息
        InitialModule_sub.f90                |      该Module包含一系列参数读取子程序，主要用于读取./ParameterSetting文件夹中参数信息
———————————————————————————————————————————————————————————————————
 

! ------------------------------------------------------------------------------------------------
! Developer(开发者): Long Cheng(成龙)
! Github: https://github.com/Reed-in-the-wind
! E-mail:  hasen_chen@163.com   &  Long.Cheng@buaa.edu.cn
! Tel:     +86 13126525737
! 通讯地址：北京市昌平区北京航空航天大学沙河校区国实FAEL课题组(2020年及之前)
！               北京市海淀区北清路156号中关村环保园华为北京研究所编译器实验室(2021年4月开始)
! ------------------------------------------------------------------------------------------------


                                                  *********************
                                                  **编译和使用案例**
                                                  *********************
在天河二号上部署计算实例： 
    
    步骤1. Compile(编译): (default: mpiifort)

          方法一(手动命令):  
                    在命令行输入: mpiifort -o NS.exe SRC_Coupled_NewSpringVersion/Main.f90   

                        or

          方法二(use Makefile): 
                    在命令行输入: make

    步骤2. 设置计算参数
          ————————————————————————————————————
          》./ParameterSetting/InitialBackgroundFlowModule.txt文件中设置特征Ma和Re数
         ————————————————————————————————————
                                    .......................
                                    .......................
                                    .......................
  	  ! -----------------------Inflow boundary conditions------
  	  ! Ma: the inflow Mach number.                           -
  	  ! Re: the inflow Reynolds number.                       -
  	  ! -------------------------------------------------------
  	  ! ---                   Ma, Re                          -
  	  ! -------------------------------------------------------
	                          0.1   1600
                  *******************************************
                  ***将马赫数Ma = 0.1, Re = 1600********
                  *******************************************

          ————————————————————————————————————
          》./ParameterSetting/InitialBackgroundFlowModule.txt文件中设置一般计算参数
         ————————————————————————————————————
                                    .......................
                                    .......................
                                    .......................
  	      ! ! ---------------------ANALYSIS MODEL DIMENSION--------
   	      ! analysis model dimension: 2-D or 3-D                  -
	      ! ModelFlag = 0: 2-D                                    -
   	      ! ModelFlag = 1: 3-D                                    -
 	      ! -------------------------------------------------------
	      ! ---                ModelFlag                          -
	      ! -------------------------------------------------------
                  	                  1                                  
	                    .......................
                                    .......................
                                    .......................
    	      ! ------------------------ N-S Flag ---------------------
                      ! NSFLAG == 1: solve N-S equations.                     -
                      ! NSFLAG == 0: solve Euler equations.                   -
                      ! -------------------------------------------------------
                      ! ---                      NSFLAG                       -
                      ! -------------------------------------------------------
                                                       1
	                    .......................
                                    .......................
                                    .......................
   	      ! -------------------COMPUTATION ZONE PARAMETER----------
  	      ! Computational zone parameter: begin from 0.           -
  	      ! MaxX, MaxY, MaxZ: computation zone length for three   -
   	      !                                          directions.  -
   	      ! SolidRelativePosi(3): relative postions for the solid.-
    	      ! UniformGridLength: stands for the uniform grid width  -
	      ! near the solid zone The real width is (blade chord    -
	      ! length+2*UniformGridLength*blade chord length )       -
  	      ! -------------------------------------------------------
	      ! ---               MaxX, MaxY, MaxZ                    -
	      ! ---            SolidRelativePosi(1:3)                 -
  	      ! ---            UniformGridRatio(1:3)
  	      ! -------------------------------------------------------
  	            6.283185306     6.283185306      6.283185306
  	            3.141592653     3.141592653      3.141592653
                  	                    0.3  0.3 0.3
	                    .......................
                                    .......................
                                    .......................
  	      ! -------------------------PROCESSOR GRID----------------
	      ! MPI processor grid:                                   -
	      !   NPX: the number for the x direction;                -
	      !   NPY: the number for the y direction;                -
	      !   NPZ: the number for the z direction;                -
	      ! -------------------------------------------------------
	      ! ---                NPX, NPY, NPZ                      -
	      ! -------------------------------------------------------
                  	            10    10   10
	                    .......................
                                    .......................
                                    .......................
                  **************************************************************************************
                  ***例如以上ModelFlag = 1代表三维问题                                                   ********
                  ***例如以上NSFLAG = 1代表N-S方程                                                   ***********
	  ***例如以上MaxX, MaxY, MaxZ代表计算域长度	                              *********
	  ***例如以上NPX, NPY, NPZ = 10代表使用10*10*10=1000个CPU核心进行计算*****
                  **************************************************************************************
                                   等等

    步骤3. 提交任务
             》 天河二号上使用slurm自动排队系统提交任务，在本文件夹下直接提交slurm_submit.sh实现，其内容为： 
              ———————————————slurm_submit.sh—————————
                        	#!/bin/bash
		#SBATCH -N 42                     (42个节点)
		#SBATCH -n 1000                （######采用1000个核心）
		####SBATCH -c 12             
		#SBATCH -t 600                     (限制运行时间为600s)
		##yhrun -n 4 ./NS.exe
		yhrun -n 1000 NS.exe >log.txt
             ————————————————————————————————
            在命令行输入： slurm slurm_submit.sh即可完成提交


**************PS::::!!!!!!!!!!!!!!!!!!!
**************
**************
*******注意以下三维MPI核心组合已经做过测试，可以使用：
       
       Here total processor number: NUMPROCS = NPX*NPY*NPZ
              It is better to choose from the following numbers: 
                   NPX = NPY = NPZ = 2    (NUMPROCS = 8);     ------NUMPROCS代表总核心数目
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
     ----(由于鬼点层限制，如果采用的网格点数为256*256*256，则NPX, NPY, NPZ <= 25)

**************PS::::!!!!!!!!!!!!!!!!!!!
**************
**************
另外，由于本程序是从本人利用高阶CFD/CAA数值格式以及浸入式边界方法开发的针对运动边界发声问题修改而来，因此有些地方稍显奇怪，比如网格尺度是由几何决定等，请根据使用说明并适当修改源程序使用，有问题请按照以上通讯方式和我联系！！！！！
                                                               5 Dec. 2019(首次书写)
                                                               Written by Long Cheng
                                                               2021年4月6号 再次修改

!--------------------------------------------------------------------------------------------------
! 全局变量模块文件
!--------------------------------------------------------------------------------------------------
! 全局变量模块
    ! 全局变量模块
    module global_Var
        use const_Var
        use precision_EC
        use block_Type_Def
        implicit none

        ! 输入参量
        integer                                     :: iDim                     ! 求解问题维度
        integer                                     :: total_Step               ! 总迭代步数
        integer                                     :: kStep_Save               ! 保存间隔迭代数
        integer                                     :: nVar                     ! 变量数目
        real(dp)                                    :: gscale                   ! 网格单位
        real(dp)                                    :: rho                      ! 材料密度
        real(dp)                                    :: cp                       ! 材料热容系数
        real(dp)                                    :: kh                       ! 热传导系数 
        real(dp)                                    :: wallq_Iso_q              ! 等壁面热流值
        real(dp)                                    :: wallc_Iso_CH             ! 等壁面对流换热系数
        real(dp)                                    :: wallc_Iso_T              ! 等壁面对热换热温度
        real(dp)                                    :: initT_Iso_T              ! 等初始温度场值
        real(dp)                                    :: dt                       ! 推进时间步
        real(dp)                                    :: real_Time                ! 真实时间
        real(dp),dimension(RKn)                     :: Ralpha                   ! 龙格库塔系数1
        real(dp),dimension(RKn)                     :: Rgamma                   ! 龙格库塔系数2
        real(dp),dimension(RKn)                     :: Rbeta                    ! 龙格库塔系数3
        character(len=str_Len )                     :: project_Name             ! 项目名称
        character(len=file_Len)                     :: mesh_File                ! 网格数据文件
        character(len=file_Len)                     :: inp_File                 ! 连接关系文件
        character(len=file_Len)                     :: save_Dir                 ! 保存文件夹
        character(len=file_Len)                     :: wallq_qfile              ! 非等壁面热流输入文件
        character(len=file_Len)                     :: wallc_CHfile             ! 非等壁面对流换热系数输入文件
        character(len=file_Len)                     :: wallc_Tfile              ! 非等壁面对流换热温度输入文件
        character(len=file_Len)                     :: initT_Tfile              ! 非等初始温度场输入文件
        character(len=file_Len)                     :: continued_File           ! 续算初始化内部温度场输入文件
        logical                                     :: Iflag_init               ! 续算控制符
        logical                                     :: wallq_Iso                ! 等壁面热流控制符
        logical                                     :: wallc_Iso                ! 等对流换热控制符
        logical                                     :: initT_Iso                ! 等初始温度场控制符

        
        namelist /prj/ project_Name,iDim,gscale
        namelist /step/ total_Step,kStep_Save,Iflag_init,continued_File
        namelist /files/ mesh_File,inp_File,save_Dir
        namelist /material/ rho,cp,kh
        namelist /wallq/ wallq_Iso,wallq_Iso_q,wallq_qfile
        namelist /wallc/ wallc_Iso,wallc_Iso_CH,wallc_Iso_T,wallc_CHfile,wallc_Tfile
        namelist /initT/ initT_Iso,initT_Iso_T,initT_Tfile

        ! 网格变量
        type(block_Type),pointer    ,dimension(:)   :: mesh                     ! 网格（数值）
        real(dp)        ,allocatable,dimension(:)   :: res_rms                  ! 均方根残差

        ! 暂态变量
        integer                                     :: num_Block                ! 网格块数


        
        
    end module global_Var

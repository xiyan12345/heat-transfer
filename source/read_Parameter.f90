!--------------------------------------------------------------------------------------------------
! 读取控制参数文件
!--------------------------------------------------------------------------------------------------
! 读入控制参数模块
    module read_Parameter
        use precision_EC
        use const_Var
        use global_Var
        implicit none
        
    contains
        subroutine read_Parameter_heat_transfer

            character(len=file_Len) :: input_File       ! 输入参数文件
            logical                 :: fexist           ! 文件存在逻辑符

            input_File = "../resource/input/control.dat"
            nVar       = 1                              ! 针对固体热传导求解器，物理量就一个：T

           ! 默认设置
            ! prj
            project_Name = ""
            iDim         = 2
            gscale       = 0.001_dp
            ! step
            total_Step     = 1000
            kStep_Save     = 100
            Iflag_init     = .true.
            continued_File = ""
            ! files
            mesh_File    = ""
            inp_File     = ""
            save_Dir     = ""
            ! material
            rho          = 900_dp
            cp           = 20_dp
            kh           = 2_dp
            ! wallq
            wallq_Iso    = .true.
            wallq_Iso_q  = 0
            wallq_qfile  = ""
            ! wallc
            wallc_Iso    = .true.
            wallc_Iso_CH = 0
            wallc_Iso_T  = 273.15
            wallc_CHfile = ""
            wallc_Tfile  = ""
            ! initT
            initT_Iso   = .true.
            initT_Iso_T = 273.15
            initT_Tfile = ""

           ! 读入控制参数文件
            inquire(file = input_File,exist = fexist)
            if (fexist) then 
                open(unit = un_control, file = input_File, status = "OLD")
            else 
                print *,"control.dat file does not exist!"
                print *,"Press any button to end the program"
                read(*,*) 
                stop
            end if

            read(unit = un_control, nml = prj)
            read(unit = un_control, nml = step)
            read(unit = un_control, nml = files)
            read(unit = un_control, nml = material)
            read(unit = un_control, nml = wallq)
            read(unit = un_control, nml = wallc)
            read(unit = un_control, nml = initT)

            close(un_control)

           ! 打开日志文件
            open(unit = un_log, file = trim(save_Dir)//"log.dat", status = "REPLACE")
            write(un_log,*) "**********************************************************************"
            write(un_log,*) "                        heat transfer solver                          "
            write(un_log,*) "                            calculating...                            "
            write(un_log,*) "**********************************************************************"

            write(un_log,*) " Read Parameter ... "
            print        * ," Read Parameter ... "


            ! 龙格库塔系数设置
            select case(RKn)
            case (1)
            case (2)
            case (3)
                Ralpha(1) = 1._dp      ; Rgamma(1) = 0._dp      ; Rbeta(1) = 1._dp      ;
                Ralpha(2) = 3._dp/4._dp; Rgamma(2) = 1._dp/4._dp; Rbeta(2) = 1._dp/4._dp;
                Ralpha(3) = 1._dp/3._dp; Rgamma(3) = 2._dp/3._dp; Rbeta(3) = 2._dp/3._dp;
            end select


        end subroutine read_Parameter_heat_transfer
        
        
    end module read_Parameter
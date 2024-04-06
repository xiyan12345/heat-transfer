!--------------------------------------------------------------------------------------------------
! 常量模块文件
!--------------------------------------------------------------------------------------------------
! 常量模块
    module const_Var
        use precision_EC
        implicit none
        ! 2D计算原理：
            ! 2D网格必须画在xy平面上，读入二维网格时，将2D转换为3D;
            ! 2D -> 3D：z向（K向）拉伸，即K=3
        ! 自定义边界类型：
            ! 由于热传导求解器的边界条件类型与流体求解的边界条件类型不一致，由于Generic格式中提供了三种自定义边界类型：
            ! 因此将它们进行自定义,各含义如下
            ! generic#1(8):直接热流边界条件
            ! generic#2(9):对流换热边界条件
            ! generuc#3(10):直接壁温边界条件（暂时没开发）
        ! unit说明：
            ! 1、6：显示屏
            ! 2、22：control.dat
            ! 3、23：log.dat
            ! 4、24：mesh.dat
            ! 5、25：inp.dat
            ! 6、26：result.dat

        integer,parameter       :: RKn              = 3     ! 龙阶-库塔阶数
        integer,parameter       :: LAP              = 2     ! 虚网格层数
        integer,parameter       :: zStretch         = 3     ! z向拉伸层数         
        integer ,parameter      :: BC_Direct_Q      = 8     ! 直接热流边界条件
        integer,parameter       :: BC_Convection    = 9     ! 对流换热边界条件
        integer,parameter       :: BC_Direct_T      = 10    ! 直接温度边界条件              
        integer,parameter       :: BC_Symmetry      = 3     ! 对称边界条件
        integer,parameter       :: un_control       = 22    ! control.dat
        integer,parameter       :: un_log           = 23    ! log.dat
        integer,parameter       :: un_mesh          = 24    ! mesh.dat
        integer,parameter       :: un_inp           = 25    ! inp.dat
        integer,parameter       :: un_result        = 26    ! result.dat
        integer,parameter       :: un_print         = 6     ! 显示屏
        real   , parameter      :: CFL_Value        = 0.8   ! CFL数
    end module const_Var
!--------------------------------------------------------------------------------------------------
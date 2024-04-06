!--------------------------------------------------------------------------------------------------
! CFL条件模块文件
!--------------------------------------------------------------------------------------------------
! CFL模块
    module CFL
        use const_Var
        use global_Var
        implicit none
        
    contains
        subroutine calculate_dt()
            implicit none

            real(dp)                :: min_delt         ! 网格尺度最小值

           ! 确定 min_delt(粗略算)
            min_delt = 0
           ! 计算dt
            dt = CFL_Value*min_delt/1
        end subroutine
        
    end module CFL
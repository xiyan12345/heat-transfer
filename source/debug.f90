!--------------------------------------------------------------------------------------------------
! debug模块文件
!--------------------------------------------------------------------------------------------------
! debug模块：此模块是用于debug之用
    module debug
        use precision_EC
        implicit none
        
    contains
        subroutine show_Flux(F)
            implicit none

           ! 输入\输出变量
            real(dp),dimension(:,:,:,:)             :: F
        end subroutine show_Flux
        
    end module debug
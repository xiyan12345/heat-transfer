!--------------------------------------------------------------------------------------------------
! 初始化模块文件
!--------------------------------------------------------------------------------------------------
! 初始化边界条件和内部温度场模块(未完成)
    module init
        use const_Var
        use global_Var
        implicit none
        
    contains
        ! 设置初始内部温度场（未完成）
        subroutine init_Field()
            if (Iflag_init) then 
                if (initT_Iso) then
                    call initT_Iso_Field        ! 等值初始化温度场
                else 
                    call initT_None_Iso_Field   ! （未完成）
                end if
            else 
                call initT_Field_File           ! （未完成）
            end if
        end subroutine init_Field

        ! 等温初始化温度场函数
        subroutine initT_Iso_Field()
            implicit none

           ! 中间参量
            integer                         :: m        ! 块索引
            type(block_Type),pointer        :: B        ! 块指针

           ! 等温初始化温度场
            do m = 1, num_Block
                B => mesh(m)
                B%U(:,1:B%nx-1,1:B%ny-1,1:B%nz-1) = initT_Iso_T
                ! if (m ==2 ) then 
                !     B%U(:,1:B%nx-1,1:B%ny-1,1:B%nz-1) = 500._dp
                ! end if
            end do
        end subroutine initT_Iso_Field
        
        ! 利用文件非等温初始化温度场函数（未完成）
        subroutine initT_None_Iso_Field()
        end subroutine initT_None_Iso_Field

        ! 读入续算文件初始化流场（未完成）
        subroutine initT_Field_File()
        end subroutine initT_Field_File

    end module init
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
            end do
        end subroutine initT_Iso_Field
        
        ! 利用文件非等温初始化温度场函数（未完成）
        subroutine initT_None_Iso_Field()
        end subroutine initT_None_Iso_Field

        ! 读入续算文件初始化流场（未完成）
        subroutine initT_Field_File()
        end subroutine initT_Field_File

        ! 初始化边界条件
        subroutine init_Boundary()
            implicit none

           ! 中间参数
            integer                     :: m        ! 块索引
            integer                     :: ksub     ! 边界索引
            type(block_Type) ,pointer   :: B        ! 块指针
            type(BC_MSG_Type),pointer   :: BC       ! 块边界指针

            do m = 1, num_Block
                B => mesh(m)
                do ksub = 1, B%subface
                    BC => B%BC_MSG(ksub)
                    select case (BC%bc)
                    case (3)
                        call init_Symmetry_Boundary(B,BC)
                    case (8)
                        call init_DirectQ_Boundary(B,BC)
                    case (9)
                        call init_Convection_Boundary(B,BC)
                    case (10)
                        call init_DirectT_Boundary(B,BC)
                    case (-1)
                        call init_Inner_Boundary(B,BC)
                    case default
                        write(un_print,*) "There isn't such boundary"
                        write(un_log,*) "There isn't such boundary"
                        stop
                    end select
                end do
            end do
        end subroutine init_Boundary

        ! 初始化对称边界条件
        subroutine init_Symmetry_Boundary(B,BC)
            implicit none 

           ! 输入/输出变量
            type(block_Type) ,pointer           :: B        ! 块指针
            type(BC_MSG_Type),pointer           :: BC       ! 边界指针

           ! 中间变量
            integer                             :: k        ! 维度指示
            integer,dimension(3)                :: db,de    ! 子面范围

            db(1) = BC%ib; db(2) = BC%jb; db(3) = BC%kb;
            de(1) = BC%ie; de(2) = BC%je; de(3) = BC%ke;
            k = sign(1,BC%face-4)*(mod(BC%face-1,3)+1)      ! k的符号代表面的正负，k的数字代表索引维度

            if ( k>0 ) then
            else 
            end if

        end subroutine init_Symmetry_Boundary

        ! 初始化直接热流边界条件
        subroutine init_DirectQ_Boundary(B,BC)
            implicit none 

            ! 输入/输出变量
             type(block_Type) ,pointer           :: B    ! 块指针
             type(BC_MSG_Type),pointer           :: BC   ! 边界指针
        end subroutine init_DirectQ_Boundary

        ! 初始化对流换热边界条件
        subroutine init_Convection_Boundary(B,BC)
            implicit none 

            ! 输入/输出变量
             type(block_Type) ,pointer           :: B    ! 块指针
             type(BC_MSG_Type),pointer           :: BC   ! 边界指针
        end subroutine init_Convection_Boundary

        ! 初始化直接温度边界条件
        subroutine init_DirectT_Boundary(B,BC)
            implicit none 

            ! 输入/输出变量
             type(block_Type) ,pointer           :: B    ! 块指针
             type(BC_MSG_Type),pointer           :: BC   ! 边界指针
        end subroutine init_DirectT_Boundary

        ! 初始化内部边界条件
        subroutine init_Inner_Boundary(B,BC)
            implicit none 

            ! 输入/输出变量
             type(block_Type) ,pointer           :: B    ! 块指针
             type(BC_MSG_Type),pointer           :: BC   ! 边界指针
        end subroutine init_Inner_Boundary
    end module init
!--------------------------------------------------------------------------------------------------
! 边界条件更新模块文件
!--------------------------------------------------------------------------------------------------
! 边界条件更新模块
    module Boundary
        use global_Var
        use const_Var
        implicit none
        
    contains
        ! 块内边界信息交换
        subroutine exchange_Inner_Boundary()
            implicit none

            integer                     :: m            ! 块索引
            integer                     :: ksub         ! 子面索引
            integer                     :: dim          ! 维度索引
            integer                     :: i,j,k        ! 主面节点索引
            integer                     :: i1,j1,k1     ! 连接子面节点索引
            integer                     :: dim_Link     ! 连接面维度及连接面类型
            integer                     :: dim_Link1    ! 连接子面维度及连接面类型
            integer,dimension(3)        :: faceb        ! 主子面起始坐标
            integer,dimension(3)        :: facee        ! 主子面终止坐标
            integer,dimension(3)        :: faceb1       ! 连接子面起始坐标
            integer,dimension(3)        :: facee1       ! 连接子面终止坐标
            integer,dimension(3)        :: owner_O      ! 主面索引原点
            integer,dimension(3)        :: owner_D      ! 主面索引距离
            integer,dimension(3)        :: conne_O      ! 连接面索引原点
            integer,dimension(3)        :: conne_D      ! 连接面索引距离
            type(block_Type),pointer    :: B            ! 块指针
            type(block_Type),pointer    :: B_Conne      ! 连接块指针
            type(BC_MSG_Type),pointer   :: BC           ! 子面指针
            

            do m = 1, num_Block
                B => mesh(m)
                do ksub = 1, B%subface
                    BC => B%BC_MSG(ksub)
                    if ( BC%face>=0 ) cycle         ! 物理边界条件跳过
                   ! 内边界
                    B_Conne => mesh(BC%nb1)
                      ! 主子面虚网格范围
                    faceb(1) = BC%ib; faceb(2) = BC%jb; faceb(3) = BC%kb;
                    facee(1) = BC%ie; facee(2) = BC%je; facee(3) = BC%ke;
                    dim_Link = sign(1,BC%face-4)*(mod(BC%face-1,3)+1)
                    faceb(abs(dim_Link)) = faceb(abs(dim_Link))+sign(1,dim_Link)
                    facee(abs(dim_Link)) = faceb(abs(dim_Link))
                   ! 连接子面网格范围
                    faceb1(1) = BC%ib1; faceb1(2) = BC%jb1; faceb1(3) = BC%kb1; 
                    facee1(1) = BC%ie1; facee1(2) = BC%je1; facee1(3) = BC%ke1;
                    dim_Link1 = sign(1,BC%face1-4)*(mod(BC%face1-1,3)+1)
                    faceb1(abs(dim_Link1)) = faceb1(abs(dim_Link1))-sign(1,dim_Link1)
                    facee1(abs(dim_Link1)) = faceb1(abs(dim_Link1))

                    ! 主面原点定位
                    owner_O(:) = faceb(:)
                    ! 连接子面原点定位
                    do dim = 1, 3
                        if ( BC%L(dim)>0 ) then
                            conne_O(dim) = faceb1(dim)
                        else 
                            conne_O(dim) = facee1(dim)
                        end if
                    end do

                    ! 相互赋值
                    do k = faceb(3), facee(3)
                        do j = faceb(2), facee(2)
                            do i = faceb(1), facee(1)
                                ! 主面距离
                                owner_D(1) = i -faceb(1)
                                owner_D(2) = j -faceb(2)
                                owner_D(3) = k -faceb(3)
                                ! 连接面距离
                                do dim = 1, 3
                                    conne_D(abs(BC%L(dim))) = owner_D(dim)*sign(1,abs(BC%L(dim)))
                                end do

                                ! (i,j,k) <-> (i1,j1,k1)
                                i1 = conne_O(1) + conne_D(1)
                                j1 = conne_O(2) + conne_D(2)
                                k1 = conne_O(3) + conne_D(3)

                                B%U(:,i,j,k) = B_Conne%U(:,i1,j1,k1)
                            end do
                        end do 
                    end do
                end do
            end do
        end subroutine exchange_Inner_Boundary

        ! 更新边界条件
        subroutine upgrade_Boundary()
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
                        call upgrade_Symmetry_Boundary(B,BC)
                    case (8)
                        call upgrade_DirectQ_Boundary(B,BC)
                    case (9)
                        call upgrade_Convection_Boundary(B,BC)
                    case (10)
                        call upgrade_DirectT_Boundary(B,BC)
                    case (-1)
                        call upgrade_Inner_Boundary(B,BC)
                    case default
                        write(un_print,*) "There isn't such boundary"
                        write(un_log,*) "There isn't such boundary"
                        stop
                    end select
                end do
            end do
        end subroutine upgrade_Boundary

        ! 更新对称边界条件
        subroutine upgrade_Symmetry_Boundary(B,BC)
            implicit none 

           ! 输入/输出变量
            type(block_Type) ,pointer           :: B        ! 块指针
            type(BC_MSG_Type),pointer           :: BC       ! 边界指针

           ! 中间变量
            integer                             :: k        ! 维度指示
            integer,dimension(3)                :: db,de    ! 子面范围
            integer,dimension(3)                :: db1,de1  ! 虚网格范围

            db(1) = BC%ib; db(2) = BC%jb; db(3) = BC%kb;
            de(1) = BC%ie; de(2) = BC%je; de(3) = BC%ke;
            k = sign(1,BC%face-4)*(mod(BC%face-1,3)+1)      ! k的符号代表面的正负，k的数字代表索引维度

            if ( k>0 ) then
            else 
            end if
            
        end subroutine upgrade_Symmetry_Boundary

        ! 更新直接热流边界条件
        subroutine upgrade_DirectQ_Boundary(B,BC)
            implicit none 

            ! 输入/输出变量
            type(block_Type) ,pointer           :: B    ! 块指针
            type(BC_MSG_Type),pointer           :: BC   ! 边界指针
        end subroutine upgrade_DirectQ_Boundary

        ! 更新对流换热边界条件
        subroutine upgrade_Convection_Boundary(B,BC)
            implicit none 

           ! 输入/输出变量
            type(block_Type) ,pointer           :: B    ! 块指针
            type(BC_MSG_Type),pointer           :: BC   ! 边界指针
        end subroutine upgrade_Convection_Boundary

        ! 更新直接温度边界条件
        subroutine upgrade_DirectT_Boundary(B,BC)
            implicit none 

           ! 输入/输出变量
            type(block_Type) ,pointer           :: B    ! 块指针
            type(BC_MSG_Type),pointer           :: BC   ! 边界指针
        end subroutine upgrade_DirectT_Boundary

        ! 更新内部边界条件
        subroutine upgrade_Inner_Boundary(B,BC)
            implicit none 

            ! 输入/输出变量
            type(block_Type) ,pointer           :: B    ! 块指针
            type(BC_MSG_Type),pointer           :: BC   ! 边界指针
        end subroutine upgrade_Inner_Boundary

    end module Boundary
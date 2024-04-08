!--------------------------------------------------------------------------------------------------
! 边界条件更新模块文件
!--------------------------------------------------------------------------------------------------
! 边界条件更新模块
    module Boundary
        use global_Var
        use const_Var
        implicit none
        
    contains

       ! 更新边界虚网格
        ! 更新所有块边界的虚网格
        subroutine upgrade_Ghost_Cell()
            implicit none

           ! 中间参数
            integer                     :: m        ! 块索引
            integer                     :: ksub     ! 边界索引
            type(block_Type) ,pointer   :: B        ! 块指针
            type(BC_MSG_Type),pointer   :: BC       ! 块边界指针

           ! 更新非角部区域虚网格
            do m = 1, num_Block
                B => mesh(m)
                do ksub = 1, B%subface
                    BC => B%BC_MSG(ksub)
                    select case (BC%bc)
                    case (3)
                        call upgrade_Symmetry_Ghost_Cell(B,BC)
                    case (8)
                        call upgrade_DirectQ_Ghost_Cell(B,BC)
                    case (9)
                        call upgrade_Convection_Ghost_Cell(B,BC)
                    case (10)
                        call upgrade_DirectT_Ghost_Cell(B,BC)
                    case (-1)
                        call upgrade_Inner_Ghost_Cell(B,BC)
                    case default
                        write(un_print,*) "There isn't such boundary"
                        write(un_log,*) "There isn't such boundary"
                        stop
                    end select
                end do
            end do
           ! 更新角部棱区域虚网格
            do m = 1, num_Block
                B => mesh(m)
                call upgrade_Corner_Cell(B)
            end do

           ! 更新内边界处、角部棱区域的虚网格
            do m = 1, num_Block
                B => mesh(m)
                do ksub = 1, B%subface
                    BC => B%BC_MSG(ksub)
                    if ( BC%bc<0 ) then         ! 内边界特殊处理
                        call upgrade_Inner_Edge_Ghost_Cell(B,BC)
                    end if
                    
                end do
            end do

        end subroutine upgrade_Ghost_Cell

        ! 更新块对称边界的非角部区域处虚网格
        subroutine  upgrade_Symmetry_Ghost_Cell(B,BC)
            implicit none 

           ! 输入/输出变量
            type(block_Type) ,pointer           :: B        ! 块指针
            type(BC_MSG_Type),pointer           :: BC       ! 边界指针

           ! 中间变量
            integer                             :: k        ! 子面维度指示
            integer                             :: dim      ! 维度索引
            integer,dimension(3)                :: faceb    ! 内部子面起始坐标
            integer,dimension(3)                :: facee    ! 内部子面终止坐标
            integer,dimension(3)                :: faceb1   ! 虚网格子面起始坐标
            integer,dimension(3)                :: facee1   ! 虚网格子面终止坐标

           ! 主子面虚网格范围
            faceb(1) = BC%ib; faceb(2) = BC%jb; faceb(3) = BC%kb;
            facee(1) = BC%ie; facee(2) = BC%je; facee(3) = BC%ke;
            k = sign(1,BC%face-4)*(mod(BC%face-1,3)+1)
            do dim = 1, 3
                if ( dim == abs(k) ) then
                ! 子面指示维度特殊处理
                faceb(dim) = faceb(dim)-(sign(1,BC%face-4)+1)/2
                facee(dim) = faceb(dim)
                else 
                ! 非子面指示维度正常处理：end节点-1
                facee(dim) = facee(dim)-1
                end if
            end do

            faceb1(:) = faceb(:); facee1(:) = facee(:);
            faceb1(abs(k)) = faceb1(abs(k))+sign(1,k)
            facee1(abs(k)) = faceb1(abs(k))

            B%U(:,faceb1(1):facee1(1),faceb1(2):facee1(2),faceb1(3):facee1(3)) =  B%U(:,faceb(1):facee(1),faceb(2):facee(2),faceb(3):facee(3))

        end subroutine  upgrade_Symmetry_Ghost_Cell

        ! 更新直接热流边界的非角部区域处虚网格
        subroutine upgrade_DirectQ_Ghost_Cell(B,BC)
            implicit none 

            ! 输入/输出变量
            type(block_Type) ,pointer           :: B    ! 块指针
            type(BC_MSG_Type),pointer           :: BC   ! 边界指针

           ! 中间变量
            integer                             :: k        ! 子面维度指示
            integer                             :: dim      ! 维度索引
            integer,dimension(3)                :: faceb    ! 内部子面1起始坐标
            integer,dimension(3)                :: facee    ! 内部子面1终止坐标
            integer,dimension(3)                :: faceb2    ! 内部子面2起始坐标
            integer,dimension(3)                :: facee2    ! 内部子面2终止坐标
            integer,dimension(3)                :: faceb1   ! 虚网格子面起始坐标
            integer,dimension(3)                :: facee1   ! 虚网格子面终止坐标

           ! 主子面虚网格范围
            faceb(1) = BC%ib; faceb(2) = BC%jb; faceb(3) = BC%kb;
            facee(1) = BC%ie; facee(2) = BC%je; facee(3) = BC%ke;
            k = sign(1,BC%face-4)*(mod(BC%face-1,3)+1)
            do dim = 1, 3
                if ( dim == abs(k) ) then
                ! 子面指示维度特殊处理
                faceb(dim) = faceb(dim)-(sign(1,BC%face-4)+1)/2
                facee(dim) = faceb(dim)
                else 
                ! 非子面指示维度正常处理：end节点-1
                facee(dim) = facee(dim)-1
                end if
            end do

            faceb1(:) = faceb(:); facee1(:) = facee(:);
            faceb1(abs(k)) = faceb1(abs(k))+sign(1,k)
            facee1(abs(k)) = faceb1(abs(k))

            faceb2(:) = faceb(:); facee2(:) = facee(:);
            faceb2(abs(k)) = faceb2(abs(k))-sign(1,k)
            facee2(abs(k)) = faceb2(abs(k))

            B%U(:,faceb1(1):facee1(1),faceb1(2):facee1(2),faceb1(3):facee1(3)) =  2*B%U(:,faceb(1):facee(1),faceb(2):facee(2),faceb(3):facee(3))&
                                                                                 &-B%U(:,faceb2(1):facee2(1),faceb2(2):facee2(2),faceb2(3):facee2(3))

        end subroutine upgrade_DirectQ_Ghost_Cell

        ! 更新对流换热边界的非角部区域处虚网格
        subroutine upgrade_Convection_Ghost_Cell(B,BC)
            implicit none 

           ! 输入/输出变量
            type(block_Type) ,pointer           :: B    ! 块指针
            type(BC_MSG_Type),pointer           :: BC   ! 边界指针

           ! 中间变量
            integer                             :: k        ! 子面维度指示
            integer                             :: dim      ! 维度索引
            integer,dimension(3)                :: faceb    ! 内部子面1起始坐标
            integer,dimension(3)                :: facee    ! 内部子面1终止坐标
            integer,dimension(3)                :: faceb2    ! 内部子面2起始坐标
            integer,dimension(3)                :: facee2    ! 内部子面2终止坐标
            integer,dimension(3)                :: faceb1   ! 虚网格子面起始坐标
            integer,dimension(3)                :: facee1   ! 虚网格子面终止坐标

           ! 主子面虚网格范围
            faceb(1) = BC%ib; faceb(2) = BC%jb; faceb(3) = BC%kb;
            facee(1) = BC%ie; facee(2) = BC%je; facee(3) = BC%ke;
            k = sign(1,BC%face-4)*(mod(BC%face-1,3)+1)
            do dim = 1, 3
                if ( dim == abs(k) ) then
                ! 子面指示维度特殊处理
                faceb(dim) = faceb(dim)-(sign(1,BC%face-4)+1)/2
                facee(dim) = faceb(dim)
                else 
                ! 非子面指示维度正常处理：end节点-1
                facee(dim) = facee(dim)-1
                end if
            end do

            faceb1(:) = faceb(:); facee1(:) = facee(:);
            faceb1(abs(k)) = faceb1(abs(k))+sign(1,k)
            facee1(abs(k)) = faceb1(abs(k))

            faceb2(:) = faceb(:); facee2(:) = facee(:);
            faceb2(abs(k)) = faceb2(abs(k))-sign(1,k)
            facee2(abs(k)) = faceb2(abs(k))

            B%U(:,faceb1(1):facee1(1),faceb1(2):facee1(2),faceb1(3):facee1(3)) =  2*B%U(:,faceb(1):facee(1),faceb(2):facee(2),faceb(3):facee(3))&
                                                                                 &-B%U(:,faceb2(1):facee2(1),faceb2(2):facee2(2),faceb2(3):facee2(3))
        end subroutine upgrade_Convection_Ghost_Cell

        ! 更新直接温度边界的非角部区域处虚网格
        subroutine upgrade_DirectT_Ghost_Cell(B,BC)
            implicit none 

           ! 输入/输出变量
            type(block_Type) ,pointer           :: B    ! 块指针
            type(BC_MSG_Type),pointer           :: BC   ! 边界指针
        end subroutine upgrade_DirectT_Ghost_Cell

        ! 更新块内边界的非角部区域处虚网格
        subroutine upgrade_Inner_Ghost_Cell(B,BC)
            implicit none

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
            


           ! 单元中心型：与之前的节点型不一致

           ! 内边界
            B_Conne => mesh(BC%nb1)
           ! 主子面虚网格范围
            faceb(1) = BC%ib; faceb(2) = BC%jb; faceb(3) = BC%kb;
            facee(1) = BC%ie; facee(2) = BC%je; facee(3) = BC%ke;
            dim_Link = sign(1,BC%face-4)*(mod(BC%face-1,3)+1)
            do dim = 1, 3
                if ( dim == abs(dim_Link) ) then
                    ! 连接维特殊处理
                    faceb(dim) = faceb(dim)+(sign(1,BC%face-4)-1)/2
                    facee(dim) = faceb(dim)
                else 
                    ! 非连接维正常处理：end节点-1
                    facee(dim) = facee(dim)-1
                end if
            end do
            
           ! 连接子面网格范围
            faceb1(1) = BC%ib1; faceb1(2) = BC%jb1; faceb1(3) = BC%kb1; 
            facee1(1) = BC%ie1; facee1(2) = BC%je1; facee1(3) = BC%ke1;
            dim_Link1 = sign(1,BC%face1-4)*(mod(BC%face1-1,3)+1)
            do dim = 1, 3
                if ( dim == abs(dim_Link1) ) then
                    ! 连接维特殊处理
                    faceb1(dim) = faceb1(dim)-(sign(1,BC%face1-4)+1)/2
                    facee1(dim) = faceb1(dim)
                else 
                    ! 非连接维正常处理：end节点-1
                    facee1(dim) = facee1(dim)-1
                end if
            end do

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
                  
        end subroutine upgrade_Inner_Ghost_Cell

        ! 更新角部区域虚网格
        subroutine upgrade_Corner_Cell(B)
            implicit none

           ! 输入\输出变量
            type(block_Type),pointer                :: B            ! 块指针

           ! 中间变量
            integer,dimension(3)                    :: s            ! 单元起始索引
            integer,dimension(3)                    :: e            ! 单元终止索引

            s(:) = 1;
            e(1)   = B%nx-1; e(2)   = B%ny-1; e(3)   = B%nz-1; 
           ! 针对每一个块插12条棱
            ! x向4条棱插值
            ! (start,start);(start,end);(end,start);(end,end)
            B%U(:,s(1):e(1),s(2)-1,s(3)-1) = B%U(:,s(1):e(1),s(2)-1,s(3))+B%U(:,s(1):e(1),s(2),s(3)-1)-B%U(:,s(1):e(1),s(2),s(3))
            B%U(:,s(1):e(1),s(2)-1,e(3)+1) = B%U(:,s(1):e(1),s(2)-1,e(3))+B%U(:,s(1):e(1),s(2),e(3)+1)-B%U(:,s(1):e(1),s(2),e(3))
            B%U(:,s(1):e(1),e(2)+1,s(3)-1) = B%U(:,s(1):e(1),e(2)+1,s(3))+B%U(:,s(1):e(1),e(2),s(3)-1)-B%U(:,s(1):e(1),e(2),s(3))
            B%U(:,s(1):e(1),e(2)+1,e(3)+1) = B%U(:,s(1):e(1),e(2)+1,e(3))+B%U(:,s(1):e(1),e(2),e(3)+1)-B%U(:,s(1):e(1),e(2),e(3))
            ! y向4条棱插值
            B%U(:,s(1)-1,s(2):e(2),s(3)-1) = B%U(:,s(1)-1,s(2):e(2),s(3))+B%U(:,s(1),s(2):e(2),s(3)-1)-B%U(:,s(1),s(2):e(2),s(3))
            B%U(:,s(1)-1,s(2):e(2),e(3)+1) = B%U(:,s(1)-1,s(2):e(2),e(3))+B%U(:,s(1),s(2):e(2),e(3)+1)-B%U(:,s(1),s(2):e(2),e(3))
            B%U(:,e(1)+1,s(2):e(2),s(3)-1) = B%U(:,e(1)+1,s(2):e(2),s(3))+B%U(:,e(1),s(2):e(2),s(3)-1)-B%U(:,e(1),s(2):e(2),s(3))
            B%U(:,e(1)+1,s(2):e(2),e(3)+1) = B%U(:,e(1)+1,s(2):e(2),e(3))+B%U(:,e(1),s(2):e(2),e(3)+1)-B%U(:,e(1),s(2):e(2),e(3))
            ! z向4条棱插值
            B%U(:,s(1)-1,s(2)-1,s(3):e(3)) = B%U(:,s(1)-1,s(2),s(3):e(3))+B%U(:,s(1),s(2)-1,s(3):e(3))-B%U(:,s(1),s(2),s(3):e(3))
            B%U(:,s(1)-1,e(2)+1,s(3):e(3)) = B%U(:,s(1)-1,e(2),s(3):e(3))+B%U(:,s(1),e(2)+1,s(3):e(3))-B%U(:,s(1),e(2),s(3):e(3))
            B%U(:,e(1)+1,s(2)-1,s(3):e(3)) = B%U(:,e(1)+1,s(2),s(3):e(3))+B%U(:,e(1),s(2)-1,s(3):e(3))-B%U(:,e(1),s(2),s(3):e(3))
            B%U(:,e(1)+1,e(2)+1,s(3):e(3)) = B%U(:,e(1)+1,e(2),s(3):e(3))+B%U(:,e(1),e(2)+1,s(3):e(3))-B%U(:,e(1),e(2),s(3):e(3))
           ! 针对内边界处的棱需要特殊处理
        end subroutine upgrade_Corner_Cell

       ! 更新内边界处、角部棱区域的虚网格
        subroutine upgrade_Inner_Edge_Ghost_Cell(B,BC)
            implicit none 

           ! 输入\输出变量
            type(block_Type),pointer            :: B                ! 块指针 
            type(BC_MSG_Type),pointer           :: BC               ! 块子面指针

           ! 中间变量
            integer                             :: i,j,k        ! 单元索引
            integer                             :: i1,j1,k1     ! 连接子面节点索引
            integer                             :: dim          ! 维度索引
            integer                             :: dim_Link     ! 连接面维度及连接面类型
            integer                             :: dim_Link1    ! 连接子面维度及连接面类型
            integer,dimension(3)                :: faceb        ! 主子面起始坐标
            integer,dimension(3)                :: facee        ! 主子面终止坐标
            integer,dimension(3)                :: facebC       ! 考虑角棱区域的主子面起始坐标
            integer,dimension(3)                :: faceeC       ! 考虑角棱区域的主子面终止坐标
            integer,dimension(3)                :: faceb1       ! 连接子面起始坐标
            integer,dimension(3)                :: facee1       ! 连接子面终止坐标
            integer,dimension(3)                :: owner_O      ! 主面索引原点
            integer,dimension(3)                :: owner_D      ! 主面索引距离
            integer,dimension(3)                :: conne_O      ! 连接面索引原点
            integer,dimension(3)                :: conne_D      ! 连接面索引距离
            type(block_Type),pointer            :: B_Conne      ! 连接块指针

           ! 单元中心型：与之前的节点型不一致
           ! 内边界
            B_Conne => mesh(BC%nb1)
           ! 主子面虚网格范围
            faceb(1) = BC%ib; faceb(2) = BC%jb; faceb(3) = BC%kb;
            facee(1) = BC%ie; facee(2) = BC%je; facee(3) = BC%ke;
            dim_Link = sign(1,BC%face-4)*(mod(BC%face-1,3)+1)
            do dim = 1, 3
                if ( dim == abs(dim_Link) ) then
                    ! 连接维特殊处理
                    faceb(dim) = faceb(dim)+(sign(1,BC%face-4)-1)/2
                    facee(dim) = faceb(dim)
                else 
                    ! 非连接维正常处理：end节点-1
                    facee(dim) = facee(dim)-1
                end if
            end do

           ! 连接子面网格范围
            faceb1(1) = BC%ib1; faceb1(2) = BC%jb1; faceb1(3) = BC%kb1; 
            facee1(1) = BC%ie1; facee1(2) = BC%je1; facee1(3) = BC%ke1;
            dim_Link1 = sign(1,BC%face1-4)*(mod(BC%face1-1,3)+1)
            do dim = 1, 3
                if ( dim == abs(dim_Link1) ) then
                    ! 连接维特殊处理
                    faceb1(dim) = faceb1(dim)-(sign(1,BC%face1-4)+1)/2
                    facee1(dim) = faceb1(dim)
                else 
                    ! 非连接维正常处理：end节点-1
                    facee1(dim) = facee1(dim)-1
                end if
            end do

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

           ! 考虑角棱区域的主子面起始/终止坐标设置
            facebC(1) = faceb(1)-1; facebC(2) = faceb(2)-1; facebC(3) = faceb(3)-1; 
            faceeC(1) = facee(1)+1; faceeC(2) = facee(2)+1; faceeC(3) = facee(3)+1; 
            facebC(abs(dim_Link)) = facebC(abs(dim_Link))+1;
            faceeC(abs(dim_Link)) = faceeC(abs(dim_Link))-1;

           ! 相互赋值
            do k = facebC(3), faceeC(3)
                do j = facebC(2), faceeC(2)
                    do i = facebC(1), faceeC(1)
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
        end subroutine upgrade_Inner_Edge_Ghost_Cell

       ! 更新边界值
        ! 更新单一块的边界值
        subroutine upgrade_Boundary(B,Fi,Fj,Fk)
            implicit none

           ! 输入\输出变量
            real(dp),dimension(:,:,:,:)         :: Fi           ! 单块I向面通量空间
            real(dp),dimension(:,:,:,:)         :: Fj           ! 单块J向面通量空间
            real(dp),dimension(:,:,:,:)         :: Fk           ! 单块k向面通量空间
            type(block_Type) ,pointer           :: B            ! 块指针

           ! 中间变量

            integer                             :: ksub         ! 边界索引
            type(BC_MSG_Type),pointer           :: BC       ! 块边界指针

            do ksub = 1, B%subface
                BC => B%BC_MSG(ksub)
                select case (BC%bc)
                case (3)
                    call upgrade_Symmetry_Boundary(B,BC,Fi,Fj,Fk)
                case (8)
                    call upgrade_DirectQ_Boundary(B,BC,Fi,Fj,Fk)
                case (9)
                    call upgrade_Convection_Boundary(B,BC,Fi,Fj,Fk)
                ! case (10)
                !     call upgrade_DirectT_Boundary(B,BC)
                case (-1)
                    call upgrade_Inner_Boundary(B,BC,Fi,Fj,Fk)
                case default
                    write(un_print,*) "There isn't such boundary"
                    write(un_log,*) "There isn't such boundary"
                    stop
                end select
            end do
        end subroutine upgrade_Boundary

        ! 更新单一块上对称边界的通量
        subroutine upgrade_Symmetry_Boundary(B,BC,Fi,Fj,Fk)
            implicit none

           ! 输入\输出变量
            real(dp),dimension(:,:,:,:)         :: Fi           ! 单块I向面通量空间
            real(dp),dimension(:,:,:,:)         :: Fj           ! 单块J向面通量空间
            real(dp),dimension(:,:,:,:)         :: Fk           ! 单块k向面通量空间
            type(block_Type) ,pointer           :: B            ! 块指针
            type(BC_MSG_Type),pointer           :: BC           ! 块边界指针

           ! 中间变量
            integer                             :: k            ! 面维度指示
            integer,dimension(3)                :: faceb,facee  ! 边界范围

           ! 子面对称边界上热流直接赋值为0
            k = sign(1,BC%face-4)*(mod(BC%face-1,3)+1)
            select case (abs(k))
            case (1)    ! i向面
                ! 查找边界面范围
                faceb(1) = BC%ib  ; faceb(2) = BC%jb  ; faceb(3) = BC%kb; 
                facee(1) = BC%ie-1; facee(2) = BC%je-1; facee(3) = BC%ke-1; 
                facee(abs(k)) = facee(abs(k))+1
                Fi(1,faceb(1):facee(1),faceb(2):facee(2),faceb(3):facee(3)) = 0 
            case (2)    ! j向面
                ! 查找边界面范围
                faceb(1) = BC%ib  ; faceb(2) = BC%jb  ; faceb(3) = BC%kb; 
                facee(1) = BC%ie-1; facee(2) = BC%je-1; facee(3) = BC%ke-1; 
                facee(abs(k)) = facee(abs(k))+1
                Fj(1,faceb(1):facee(1),faceb(2):facee(2),faceb(3):facee(3)) = 0
            case (3)    ! k向面
                ! 查找边界面范围
                faceb(1) = BC%ib  ; faceb(2) = BC%jb  ; faceb(3) = BC%kb; 
                facee(1) = BC%ie-1; facee(2) = BC%je-1; facee(3) = BC%ke-1; 
                facee(abs(k)) = facee(abs(k))+1
                Fk(1,faceb(1):facee(1),faceb(2):facee(2),faceb(3):facee(3)) = 0 
            end select
           

        end subroutine upgrade_Symmetry_Boundary

        ! 更新单一块上直接热流边界的通量
        subroutine upgrade_DirectQ_Boundary(B,BC,Fi,Fj,Fk)
            implicit none

            ! 输入\输出变量
             real(dp),dimension(:,:,:,:)         :: Fi           ! 单块I向面通量空间
             real(dp),dimension(:,:,:,:)         :: Fj           ! 单块J向面通量空间
             real(dp),dimension(:,:,:,:)         :: Fk           ! 单块k向面通量空间
             type(block_Type) ,pointer           :: B            ! 块指针
             type(BC_MSG_Type),pointer           :: BC           ! 块边界指针

           ! 中间变量
            integer                             :: k            ! 面维度指示
            integer,dimension(3)                :: faceb,facee  ! 边界范围

           ! 边界上热流分为两部分：1、等边界热流；2、非等边界热流
            if (wallq_Iso) then
               ! 等边界热流
                k = sign(1,BC%face-4)*(mod(BC%face-1,3)+1)
                select case (abs(k))
                case (1)    ! i向面
                    ! 查找边界面范围
                    faceb(1) = BC%ib  ; faceb(2) = BC%jb  ; faceb(3) = BC%kb; 
                    facee(1) = BC%ie-1; facee(2) = BC%je-1; facee(3) = BC%ke-1; 
                    facee(abs(k)) = facee(abs(k))+1
                    Fi(1,faceb(1):facee(1),faceb(2):facee(2),faceb(3):facee(3)) = wallq_Iso_q 
                case (2)    ! j向面
                    ! 查找边界面范围
                    faceb(1) = BC%ib  ; faceb(2) = BC%jb  ; faceb(3) = BC%kb; 
                    facee(1) = BC%ie-1; facee(2) = BC%je-1; facee(3) = BC%ke-1; 
                    facee(abs(k)) = facee(abs(k))+1
                    Fj(1,faceb(1):facee(1),faceb(2):facee(2),faceb(3):facee(3)) = wallq_Iso_q 
                case (3)    ! k向面
                    ! 查找边界面范围
                    faceb(1) = BC%ib  ; faceb(2) = BC%jb  ; faceb(3) = BC%kb; 
                    facee(1) = BC%ie-1; facee(2) = BC%je-1; facee(3) = BC%ke-1; 
                    facee(abs(k)) = facee(abs(k))+1
                    Fk(1,faceb(1):facee(1),faceb(2):facee(2),faceb(3):facee(3)) = wallq_Iso_q  
                end select
            else 
               ! 非等边界热流（未完成）
            end if
 
        end subroutine upgrade_DirectQ_Boundary

        ! 更新单一块上对流边界的通量
        subroutine upgrade_Convection_Boundary(B,BC,Fi,Fj,Fk)
            implicit none

            ! 输入\输出变量
             real(dp),dimension(:,:,:,:)         :: Fi                  ! 单块I向面通量空间
             real(dp),dimension(:,:,:,:)         :: Fj                  ! 单块J向面通量空间
             real(dp),dimension(:,:,:,:)         :: Fk                  ! 单块k向面通量空间
             type(block_Type) ,pointer           :: B                   ! 块指针
             type(BC_MSG_Type),pointer           :: BC                  ! 块边界指针

             ! 中间变量
            integer                             :: dim_Sign             ! 面维度指示
            integer                             :: i,j,k                ! 面维度索引
            integer                             :: cellI,cellJ,cellK    ! 单元维度索引
            integer,dimension(3)                :: faceb,facee          ! 边界范围

           ! 对流换热边界上的对流换热分为两部分：1、等对流换热；2、非等对流换热
            if (wallq_Iso) then
               ! 等对流换热
                dim_Sign = sign(1,BC%face-4)*(mod(BC%face-1,3)+1)
                select case (abs(dim_Sign))
                case (1)    ! i向面
                    ! 查找边界面范围
                    faceb(1) = BC%ib  ; faceb(2) = BC%jb  ; faceb(3) = BC%kb; 
                    facee(1) = BC%ie-1; facee(2) = BC%je-1; facee(3) = BC%ke-1; 
                    facee(abs(dim_Sign)) = facee(abs(dim_Sign))+1
                    do k = faceb(3), facee(3)
                        do j = faceb(2), facee(2)
                            do i = faceb(1), facee(1)
                                cellI = i-(sign(1,dim_Sign)+1)/2
                                cellJ = j
                                cellK = k
                                Fi(1,i,j,k) = wallc_Iso_CH*(wallc_Iso_T-B%U(1,cellI,cellJ,cellK))
                            end do
                        end do
                    end do
                case (2)    ! j向面
                    ! 查找边界面范围
                    faceb(1) = BC%ib  ; faceb(2) = BC%jb  ; faceb(3) = BC%kb; 
                    facee(1) = BC%ie-1; facee(2) = BC%je-1; facee(3) = BC%ke-1; 
                    facee(abs(dim_Sign)) = facee(abs(dim_Sign))+1
                    do k = faceb(3), facee(3)
                        do j = faceb(2), facee(2)
                            do i = faceb(1), facee(1)
                                cellI = i
                                cellJ = j-(sign(1,dim_Sign)+1)/2
                                cellK = k
                                Fj(1,i,j,k) = wallc_Iso_CH*(wallc_Iso_T-B%U(1,cellI,cellJ,cellK))
                            end do
                        end do
                    end do
                case (3)    ! k向面
                    ! 查找边界面范围
                    faceb(1) = BC%ib  ; faceb(2) = BC%jb  ; faceb(3) = BC%kb; 
                    facee(1) = BC%ie-1; facee(2) = BC%je-1; facee(3) = BC%ke-1; 
                    facee(abs(dim_Sign)) = facee(abs(dim_Sign))+1
                    do k = faceb(3), facee(3)
                        do j = faceb(2), facee(2)
                            do i = faceb(1), facee(1)
                                cellI = i
                                cellJ = j
                                cellK = k-(sign(1,dim_Sign)+1)/2
                                Fk(1,i,j,k) = wallc_Iso_CH*(wallc_Iso_T-B%U(1,cellI,cellJ,cellK))
                            end do
                        end do
                    end do
                end select
            else 
               ! 非等对流换热（未完成）
            end if
             
 
        end subroutine upgrade_Convection_Boundary

        ! 更新单一块上对流边界的通量(未完成)
        subroutine upgrade_Inner_Boundary(B,BC,Fi,Fj,Fk)
            implicit none

            ! 输入\输出变量
             real(dp),dimension(:,:,:,:)         :: Fi           ! 单块I向面通量空间
             real(dp),dimension(:,:,:,:)         :: Fj           ! 单块J向面通量空间
             real(dp),dimension(:,:,:,:)         :: Fk           ! 单块k向面通量空间
             type(block_Type) ,pointer           :: B            ! 块指针
             type(BC_MSG_Type),pointer           :: BC           ! 块边界指针
 
        end subroutine upgrade_Inner_Boundary


    end module Boundary
!--------------------------------------------------------------------------------------------------
! 读取网格节点和连接参数文件
!--------------------------------------------------------------------------------------------------
! 读网格参数并计算几何量模块
    module read_Mesh
        use const_Var
        use global_Var
        implicit none
        
    contains
        subroutine read_Mesh_2D()
            integer                     :: m            ! 块索引变量
            integer                     :: ksub         ! 子面索引变量
            integer                     :: i,j,k        ! 块的索引 
            integer,dimension(3)        :: ib,ie,ib1,ie1  ! 子面大小
            logical                     :: fexist       ! 文件存在逻辑符
            real(dp)                    :: dz           ! z向拉伸时的间隔
            type(block_Type),pointer    :: B            ! 块指针
            type(BC_MSG_Type),pointer   :: BC           ! 块边界指针

           ! 读入网格节点文件
            inquire(file = mesh_File, exist = fexist)
            if (fexist) then
                open(unit = un_mesh, file = mesh_File, status = "OLD")
            else 
                print *,"mesh file does not exist!"
                print *,"Press any button to end the program"
                read(*,*) 
                stop
            end if

            print *, " Read mesh ... "
            write(un_log,*) " Read mesh ... "
            read(un_mesh,*) num_Block
            allocate(mesh(num_Block))

            do m = 1, num_Block
                B => mesh(m)
                read(un_mesh,*) B%nx,B%ny,B%nz
                B%nz = zStretch     ! z向拉伸
                dz   = 10           ! z向拉伸时的间隔

               ! 空间分配
                ! 节点型参量分配
                allocate(B%x(1-LAP:B%nx+LAP,1-LAP:B%ny+LAP,1-LAP:B%nz+LAP))
                allocate(B%y(1-LAP:B%nx+LAP,1-LAP:B%ny+LAP,1-LAP:B%nz+LAP))
                allocate(B%z(1-LAP:B%nx+LAP,1-LAP:B%ny+LAP,1-LAP:B%nz+LAP))

                ! i向面参量分配
                allocate(B%Si(1-LAP:B%nx+LAP,1-LAP:B%ny+LAP-1,1-LAP:B%nz+LAP-1))
                allocate(B%ni1(1-LAP:B%nx+LAP,1-LAP:B%ny+LAP-1,1-LAP:B%nz+LAP-1))
                allocate(B%ni2(1-LAP:B%nx+LAP,1-LAP:B%ny+LAP-1,1-LAP:B%nz+LAP-1))
                allocate(B%ni3(1-LAP:B%nx+LAP,1-LAP:B%ny+LAP-1,1-LAP:B%nz+LAP-1))

                ! j向面参量分配
                allocate(B%Sj(1-LAP:B%nx+LAP-1,1-LAP:B%ny+LAP,1-LAP:B%nz+LAP-1))
                allocate(B%nj1(1-LAP:B%nx+LAP-1,1-LAP:B%ny+LAP,1-LAP:B%nz+LAP-1))
                allocate(B%nj2(1-LAP:B%nx+LAP-1,1-LAP:B%ny+LAP,1-LAP:B%nz+LAP-1))
                allocate(B%nj3(1-LAP:B%nx+LAP-1,1-LAP:B%ny+LAP,1-LAP:B%nz+LAP-1))

                ! k向面参量分配
                allocate(B%Sk(1-LAP:B%nx+LAP-1,1-LAP:B%ny+LAP-1,1-LAP:B%nz+LAP))
                allocate(B%nk1(1-LAP:B%nx+LAP-1,1-LAP:B%ny+LAP-1,1-LAP:B%nz+LAP))
                allocate(B%nk2(1-LAP:B%nx+LAP-1,1-LAP:B%ny+LAP-1,1-LAP:B%nz+LAP))
                allocate(B%nk3(1-LAP:B%nx+LAP-1,1-LAP:B%ny+LAP-1,1-LAP:B%nz+LAP))

                ! 单元参量分配
                allocate(B%vol(1-LAP:B%nx+LAP-1,1-LAP:B%ny+LAP-1,1-LAP:B%nz+LAP-1))
                allocate(B%U(nVar,1-LAP:B%nx+LAP-1,1-LAP:B%ny+LAP-1,1-LAP:B%nz+LAP-1))
                allocate(B%Un(nVar,1-LAP:B%nx+LAP-1,1-LAP:B%ny+LAP-1,1-LAP:B%nz+LAP-1))
                allocate(B%res(nVar,1-LAP:B%nx+LAP-1,1-LAP:B%ny+LAP-1,1-LAP:B%nz+LAP-1))

               ! 读入网格节点坐标
                read(un_mesh,*) (((B%x(i,j,k),i=1,B%nx),j=1,B%ny),k=1,1)
                read(un_mesh,*) (((B%y(i,j,k),i=1,B%nx),j=1,B%ny),k=1,1)
                read(un_mesh,*) (((B%z(i,j,k),i=1,B%nx),j=1,B%ny),k=1,1)

               ! 2D -> 3D：z向拉伸
                do k = 2, zStretch
                    B%x(:,:,k) = B%x(:,:,1)
                    B%y(:,:,k) = B%y(:,:,1)
                    B%z(:,:,k) = B%z(:,:,1)+dz*(k-1)
                end do
            end do
            close(un_mesh)
           ! 读入连接关系文件
            inquire(file = inp_File, exist = fexist)
            if (fexist) then
                open(unit = un_inp, file = inp_File, status = "OLD")
            else 
                print *,"inp file does not exist!"
                print *,"Press any button to end the program"
                read(*,*) 
            end if

            read(un_inp,*)
            read(un_inp,*) 

            do m = 1, num_Block
                B => mesh(m)
                read(un_inp,*)
                read(un_inp,*)
                read(un_inp,*) B%subface
                B%subface = B%subface+2                         ! 附加两个对称面
                allocate(B%BC_MSG(B%subface))

                do ksub = 1, B%subface-2                        ! 处理文件中的子面
                    BC => B%BC_MSG(ksub)
                    read(un_inp,*) ib(1),ie(1),ib(2),ie(2),BC%bc
                    ib(3) = -1; ie(3) = -1*zStretch;
                    if (BC%bc < 0 ) then
                        read(un_inp,*) ib1(1),ie1(1),ib1(2),ie1(2),BC%nb1
                        ib1(3) = -1; ie1(3) = -1*zStretch;
                    end if
                    call convert_BC(BC,ib,ie,ib1,ie1)           ! 确定子面和连接子面的信息
                end do

                do ksub = B%subface-1, B%subface                ! 附加两个对称面的子面信息处理
                    BC => B%BC_MSG(ksub)
                    ib(1) = 1; ie(1) = B%nx; ib(2) = 1; ie(2) = B%ny;
                    if (ksub == B%subface-1) then
                        ib(3) = 1; ie(3) = 1
                    else 
                        ib(3) = zStretch; ie(3) = zStretch
                    end if
                    BC%bc = 3
                    call convert_BC(BC,ib,ie,ib1,ie1)           ! 确定子面和连接子面的信息 
                end do

            end do
            close(un_inp)
           ! 针对每一块，生成虚网格坐标
            do m = 1, num_Block
                B => mesh(m)
                call upgrade_ghost(B,LAP)
            end do
           ! 针对每一块，计算几何量
            do m = 1, num_Block
                B => mesh(m)
                call compute_Geometry(B)
            end do

        end subroutine read_Mesh_2D

       ! 确定子面和连接子面信息函数
        subroutine convert_BC(BC,ib,ie,ib1,ie1)
            implicit none
           ! 输入/输出参数
            integer,dimension(:)        :: ib,ie 
            integer,dimension(:)        :: ib1,ie1
            type(BC_MSG_Type)           :: BC

           ! 暂态参数
            integer                     :: dim          ! 维度索引
            integer                     :: linkDim      ! 连接面维度索引
            integer,dimension(3)        :: owner        ! 主面维度标识
            integer,dimension(3)        :: conne        ! 连接维度标识

           ! 设置主面参数
            BC%ib = min(abs(ib(1)),abs(ie(1))); BC%ie = max(abs(ib(1)),abs(ie(1)))
            BC%jb = min(abs(ib(2)),abs(ie(2))); BC%je = max(abs(ib(2)),abs(ie(2)))
            BC%kb = min(abs(ib(3)),abs(ie(3))); BC%ke = max(abs(ib(3)),abs(ie(3)))
            
            do dim = 1, 3

                if (ib(dim) == ie(dim)) then 
                    if(ib(dim) == 1) then
                        owner(dim) = -1
                        BC%face    = dim
                    else 
                        owner(dim) = 1
                        BC%face    = dim+3
                    end if
                    cycle
                else if (ib(dim) > 0 ) then
                    owner(dim) = sign(1,ie(dim) -ib(dim))*2
                else 
                    owner(dim) = sign(1,ie(dim) -ib(dim))*3
                end if

            end do

           ! 设置连接面参数
            if (BC%bc > 0) then
            ! 物理边界：连接面信息全部置0
                BC%ib1 = 0; BC%jb1 = 0; BC%kb1 = 0;
                BC%ie1 = 0; BC%je1 = 0; BC%ke1 = 0;
                BC%nb1 = 0; BC%face1 = 0
            else 
            ! 内边界：设置连接面信息
                BC%ib1 = min(abs(ib1(1)),abs(ie1(1))); BC%ie1 = max(abs(ib1(1)),abs(ie1(1)))
                BC%jb1 = min(abs(ib1(2)),abs(ie1(2))); BC%je1 = max(abs(ib1(2)),abs(ie1(2)))
                BC%kb1 = min(abs(ib1(3)),abs(ie1(3))); BC%ke1 = max(abs(ib1(3)),abs(ie1(3)))

                do dim = 1, 3

                    if (ib1(dim) == ie1(dim)) then 
                        if(ib1(dim) == 1) then
                            conne(dim) = -1
                            BC%face1   = dim
                        else 
                            conne(dim) = 1
                            BC%face1   = dim+3
                        end if
                        cycle
                    else if ( ib1(dim) > 0 ) then
                        conne(dim) = sign(1,ie1(dim) -ib1(dim))*2
                    else 
                        conne(dim) = sign(1,ie1(dim) -ib1(dim))*3
                    end if
                end do

                ! 设置顺序连接符
                do dim = 1, 3
                    do linkDim = 1, 3
                        if (abs(owner(dim)) == abs(conne(linkDim))) then
                            BC%L(dim) = linkDim*SIGN(1,owner(dim)*conne(linkDim))
                            exit
                        end if
                    end do
                end do
                
            end if
          



            
        end subroutine convert_BC

       ! 产生虚网格坐标函数
        subroutine upgrade_ghost(B,LAP)
            implicit none
           ! 输入/输出参数
            type(block_Type),pointer        :: B
            integer                         :: LAP

           ! 暂态参数
            integer                         :: n            ! 虚网格层数索引
            integer                         :: ksub         ! 子面索引
            integer                         :: dim_Link     ! 连接面维度及连接面类型
            integer                         :: dim_Link1    ! 连接子面维度及连接面类型
            integer                         :: i,j,k        ! 主面节点索引
            integer                         :: i1,j1,k1     ! 连接子面节点索引
            integer                         :: dim          ! 维度索引
            integer,dimension(3)            :: faceb        ! 主子面起始坐标
            integer,dimension(3)            :: facee        ! 主子面终止坐标
            integer,dimension(3)            :: faceb1       ! 连接子面起始坐标
            integer,dimension(3)            :: facee1       ! 连接子面终止坐标
            integer,dimension(3)            :: owner_O      ! 主面索引原点
            integer,dimension(3)            :: owner_D      ! 主面索引距离
            integer,dimension(3)            :: conne_O      ! 连接面索引原点
            integer,dimension(3)            :: conne_D      ! 连接面索引距离
            type(BC_MSG_Type),pointer       :: BC           ! 子面指针
            type(block_Type) ,pointer       :: B_Conne      ! 连接子块指针


           ! 循环产生虚网格
            do n = 1, LAP

               ! 产生非角点区域虚网格坐标
                do ksub = 1, B%subface
                    BC => B%BC_MSG(ksub)
                    if ( BC%bc >= 0) then
                    ! 物理边界
                        faceb(1) = BC%ib; faceb(2) = BC%jb; faceb(3) = BC%kb;
                        facee(1) = BC%ie; facee(2) = BC%je; facee(3) = BC%ke;
                        dim_Link = sign(1,BC%face-4)*(mod(BC%face-1,3)+1)
                        faceb(abs(dim_Link)) = faceb(abs(dim_Link))+sign(1,dim_Link)*n
                        facee(abs(dim_Link)) = faceb(abs(dim_Link))
                        select case (abs(dim_Link))
                        case (1)
                            B%x(faceb(1),faceb(2):facee(2),faceb(3):facee(3)) = 2_dp*B%x(faceb(1)-sign(1,dim_Link),faceb(2):facee(2),faceb(3):facee(3))-B%x(faceb(1)-sign(1,dim_Link)*2,faceb(2):facee(2),faceb(3):facee(3))
                            B%y(faceb(1),faceb(2):facee(2),faceb(3):facee(3)) = 2_dp*B%y(faceb(1)-sign(1,dim_Link),faceb(2):facee(2),faceb(3):facee(3))-B%y(faceb(1)-sign(1,dim_Link)*2,faceb(2):facee(2),faceb(3):facee(3))
                            B%z(faceb(1),faceb(2):facee(2),faceb(3):facee(3)) = 2_dp*B%z(faceb(1)-sign(1,dim_Link),faceb(2):facee(2),faceb(3):facee(3))-B%z(faceb(1)-sign(1,dim_Link)*2,faceb(2):facee(2),faceb(3):facee(3))
                        case (2)
                            B%x(faceb(1):facee(1),faceb(2),faceb(3):facee(3)) = 2_dp*B%x(faceb(1):facee(1),faceb(2)-sign(1,dim_Link),faceb(3):facee(3))-B%x(faceb(1):facee(1),faceb(2)-sign(1,dim_Link)*2,faceb(3):facee(3))
                            B%y(faceb(1):facee(1),faceb(2),faceb(3):facee(3)) = 2_dp*B%y(faceb(1):facee(1),faceb(2)-sign(1,dim_Link),faceb(3):facee(3))-B%y(faceb(1):facee(1),faceb(2)-sign(1,dim_Link)*2,faceb(3):facee(3))
                            B%z(faceb(1):facee(1),faceb(2),faceb(3):facee(3)) = 2_dp*B%z(faceb(1):facee(1),faceb(2)-sign(1,dim_Link),faceb(3):facee(3))-B%z(faceb(1):facee(1),faceb(2)-sign(1,dim_Link)*2,faceb(3):facee(3))
                        case (3)
                            B%x(faceb(1):facee(1),faceb(2):facee(2),faceb(3)) = 2_dp*B%x(faceb(1):facee(1),faceb(2):facee(2),faceb(3)-sign(1,dim_Link))-B%x(faceb(1):facee(1),faceb(2):facee(3),faceb(3)-sign(1,dim_Link)*2)
                            B%y(faceb(1):facee(1),faceb(2):facee(2),faceb(3)) = 2_dp*B%y(faceb(1):facee(1),faceb(2):facee(2),faceb(3)-sign(1,dim_Link))-B%y(faceb(1):facee(1),faceb(2):facee(3),faceb(3)-sign(1,dim_Link)*2)
                            B%z(faceb(1):facee(1),faceb(2):facee(2),faceb(3)) = 2_dp*B%z(faceb(1):facee(1),faceb(2):facee(2),faceb(3)-sign(1,dim_Link))-B%z(faceb(1):facee(1),faceb(2):facee(3),faceb(3)-sign(1,dim_Link)*2)
                        end select
                    else 
                    ! 内边界
                        B_Conne => mesh(BC%nb1)
                       ! 主子面虚网格范围
                        faceb(1) = BC%ib; faceb(2) = BC%jb; faceb(3) = BC%kb;
                        facee(1) = BC%ie; facee(2) = BC%je; facee(3) = BC%ke;
                        dim_Link = sign(1,BC%face-4)*(mod(BC%face-1,3)+1)
                        faceb(abs(dim_Link)) = faceb(abs(dim_Link))+sign(1,dim_Link)*n
                        facee(abs(dim_Link)) = faceb(abs(dim_Link))
                       ! 连接子面网格范围
                        faceb1(1) = BC%ib1; faceb1(2) = BC%jb1; faceb1(3) = BC%kb1; 
                        facee1(1) = BC%ie1; facee1(2) = BC%je1; facee1(3) = BC%ke1;
                        dim_Link1 = sign(1,BC%face1-4)*(mod(BC%face1-1,3)+1)
                        faceb1(abs(dim_Link1)) = faceb1(abs(dim_Link1))-sign(1,dim_Link1)*n
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

                                    B%x(i,j,k) = B_Conne%x(i1,j1,k1)
                                    B%y(i,j,k) = B_Conne%y(i1,j1,k1)
                                    B%z(i,j,k) = B_Conne%z(i1,j1,k1)
                                end do
                            end do
                        end do
                    end if  
                end do

               ! 产生角点区域虚网格坐标
               ! 12条棱插值
                ! x向4条棱插值
                call get_Corner_Value(B%x,1,1,1,B%nx,B%ny,B%nz,n)
                call get_Corner_Value(B%y,1,1,1,B%nx,B%ny,B%nz,n)
                call get_Corner_Value(B%z,1,1,1,B%nx,B%ny,B%nz,n)

                    
                
                
            end do
        end subroutine upgrade_ghost

        ! 角部插值函数
        subroutine get_Corner_Value(x,ib,jb,kb,ie,je,ke,n)
            implicit none
           ! 输入/输出参数
            integer                                                         :: n
            integer                                                         :: ib,jb,kb     ! 原始块起点
            integer                                                         :: ie,je,ke     ! 原始块终点 
            real(dp),dimension(ib-LAP:ie+LAP,jb-LAP:je+LAP,kb-LAP:ke+LAP)   :: x

            
           ! 中间参数
            integer                     :: i,j,k        ! 节点索引
            integer                     :: p            ! 遍历变量
            integer                     :: dim          ! 维度索引
            integer,dimension(3)        :: sign_Cor     ! 角部区域维度符号变量
            integer,dimension(3)        :: old_Cord     ! 角点区域旧角点坐标
            integer,dimension(3)        :: new_Cord     ! 角点区域新角点坐标
            integer,dimension(3)        :: b,e          ! 新角点坐标
            real(dp)                    :: x1,x2        ! 计算暂态变量
            


            b(1) = ib-n; b(2) = jb-n; b(3) = kb-n;
            e(1) = ie+n; e(2) = je+n; e(3) = ke+n;

            do p = 1, n
               ! 插x向4条棱
                ! start:代表相减；end：代表相加；（old,new)；(new,old)
                ! (start,start)
                x(ib:ie,jb-p,b(3)) = x(ib:ie,jb-p+1,b(3))+x(ib:ie,jb-p,b(3)+1)-x(ib:ie,jb-p+1,b(3)+1)
                x(ib:ie,b(2),kb-p) = x(ib:ie,b(2)+1,kb-p)+x(ib:ie,b(2),kb-p+1)-x(ib:ie,b(2)+1,kb-p+1)
                ! (start,end）
                x(ib:ie,jb-p,e(3)) = x(ib:ie,jb-p+1,e(3))+x(ib:ie,jb-p,e(3)-1)-x(ib:ie,jb-p+1,e(3)-1)
                x(ib:ie,b(2),ke+p) = x(ib:ie,b(2)+1,ke+p)+x(ib:ie,b(2),ke+p-1)-x(ib:ie,b(2)+1,ke+p-1)
                ! (end,start)
                x(ib:ie,je+p,b(3)) = x(ib:ie,je+p-1,b(3))+x(ib:ie,je+p,b(3)+1)-x(ib:ie,je+p-1,b(3)+1)
                x(ib:ie,e(2),kb-p) = x(ib:ie,e(2)-1,kb-p)+x(ib:ie,e(2),kb-p+1)-x(ib:ie,e(2)-1,kb-p+1)
                ! (end,end)
                x(ib:ie,je+p,e(3)) = x(ib:ie,je+p-1,e(3))+x(ib:ie,je+p,e(3)-1)-x(ib:ie,je+p-1,e(3)-1)
                x(ib:ie,e(2),ke+p) = x(ib:ie,e(2)-1,ke+p)+x(ib:ie,e(2),ke+p-1)-x(ib:ie,e(2)-1,ke+p-1)

               ! 插y向4条棱
                ! start:代表相减；end：代表相加；（old,new)；(new,old)
                ! (start,start)（problem)
                x(ib-p,jb:je,b(3)) = x(ib-p+1,jb:je,b(3))+x(ib-p,jb:je,b(3)+1)-x(ib-p+1,jb:je,b(3)+1)
                x(b(1),jb:je,kb-p) = x(b(1)+1,jb:je,kb-p)+x(b(1),jb:je,kb-p+1)-x(b(1)+1,jb:je,kb-p+1)
                ! (start,end）
                x(ib-p,jb:je,e(3)) = x(ib-p+1,jb:je,e(3))+x(ib-p,jb:je,e(3)-1)-x(ib-p+1,jb:je,e(3)-1)
                x(b(1),jb:je,ke+p) = x(b(1)+1,jb:je,ke+p)+x(b(1),jb:je,ke+p-1)-x(b(1)+1,jb:je,ke+p-1)
                ! (end,start)
                x(ie+p,jb:je,b(3)) = x(ie+p-1,jb:je,b(3))+x(ie+p,jb:je,b(3)+1)-x(ie+p-1,jb:je,b(3)+1)
                x(e(1),jb:je,kb-p) = x(e(1)-1,jb:je,kb-p)+x(e(1),jb:je,kb-p+1)-x(e(1)-1,jb:je,kb-p+1)
                ! (end,end)(problem)
                x(ie+p,jb:je,e(3)) = x(ie+p-1,jb:je,e(3))+x(ie+p,jb:je,e(3)-1)-x(ie+p-1,jb:je,e(3)-1)
                x(e(1),jb:je,ke+p) = x(e(1)-1,jb:je,ke+p)+x(e(1),jb:je,ke+p-1)-x(e(1)-1,jb:je,ke+p-1)
            
               ! 插z向4条棱
                ! start:代表相减；end：代表相加；（old,new)；(new,old)
                ! (start,start)
                x(ib-p,b(2),kb:ke) = x(ib-p+1,b(2),kb:ke)+x(ib-p,b(2)+1,kb:ke)-x(ib-p+1,b(2)+1,kb:ke)
                x(b(1),jb-p,kb:ke) = x(b(1)+1,jb-p,kb:ke)+x(b(1),jb-p+1,kb:ke)-x(b(1)+1,jb-p+1,kb:ke)
                ! (start,end)
                x(ib-p,e(2),kb:ke) = x(ib-p+1,e(2),kb:ke)+x(ib-p,e(2)-1,kb:ke)-x(ib-p+1,e(2)-1,kb:ke)
                x(b(1),je+p,kb:ke) = x(b(1)+1,je+p,kb:ke)+x(b(1),je+p-1,kb:ke)-x(b(1)+1,je+p-1,kb:ke)
                ! (end,start)
                x(ie+p,b(2),kb:ke) = x(ie+p-1,b(2),kb:ke)+x(ie+p,b(2)+1,kb:ke)-x(ie+p-1,b(2)+1,kb:ke)
                x(e(1),jb-p,kb:ke) = x(e(1)-1,jb-p,kb:ke)+x(e(1),jb-p+1,kb:ke)-x(e(1)-1,jb-p+1,kb:ke)
                ! (end,end)
                x(ie+p,e(2),kb:ke) = x(ie+p-1,e(2),kb:ke)+x(ie+p,e(2)-1,kb:ke)-x(ie+p-1,e(2)-1,kb:ke)
                x(e(1),je+p,kb:ke) = x(e(1)-1,je+p,kb:ke)+x(e(1),je+p-1,kb:ke)-x(e(1)-1,je+p-1,kb:ke)

            end do

           ! 插8个角点
            ! (start,start,start)
            old_Cord(1) = ib  ; old_Cord(2) = jb  ; old_Cord(3) = kb  ; 
            new_Cord(1) = b(1); new_Cord(2) = b(2); new_Cord(3) = b(3);
            do dim =1,3
                sign_Cor(dim) = sign(1,new_Cord(dim)-old_Cord(dim))
            end do

            do k = old_Cord(3)+sign_Cor(3), new_Cord(3),sign_Cor(3)
                do j = old_Cord(2)+sign_Cor(2), new_Cord(2),sign_Cor(2)
                    do i = old_Cord(1)+sign_Cor(1), new_Cord(1),sign_Cor(1)
                        x1 = (x(i-sign_Cor(1),j,k)+x(i,j-sign_Cor(2),k)+x(i,j,k-sign_Cor(3)))/3_dp
                        x2 = (x(i,j-sign_Cor(2),k-sign_Cor(3))+x(i-sign_Cor(1),j,k-sign_Cor(3))+x(i-sign_Cor(1),j-sign_Cor(2),k))/3_dp
                        x(i,j,k) = 2_dp*x1-x2
                    end do
                end do
            end do
            ! (start,end,start)
            old_Cord(1) = ib  ; old_Cord(2) = je  ; old_Cord(3) = kb  ; 
            new_Cord(1) = b(1); new_Cord(2) = e(2); new_Cord(3) = b(3);
            do dim =1,3
                sign_Cor(dim) = sign(1,new_Cord(dim)-old_Cord(dim))
            end do

            do k = old_Cord(3)+sign_Cor(3), new_Cord(3),sign_Cor(3)
                do j = old_Cord(2)+sign_Cor(2), new_Cord(2),sign_Cor(2)
                    do i = old_Cord(1)+sign_Cor(1), new_Cord(1),sign_Cor(1)
                        x1 = (x(i-sign_Cor(1),j,k)+x(i,j-sign_Cor(2),k)+x(i,j,k-sign_Cor(3)))/3_dp
                        x2 = (x(i,j-sign_Cor(2),k-sign_Cor(3))+x(i-sign_Cor(1),j,k-sign_Cor(3))+x(i-sign_Cor(1),j-sign_Cor(2),k))/3_dp
                        x(i,j,k) = 2_dp*x1-x2
                    end do
                end do
            end do
            ! (end,start,start)
            old_Cord(1) = ie  ; old_Cord(2) = jb  ; old_Cord(3) = kb  ; 
            new_Cord(1) = e(1); new_Cord(2) = b(2); new_Cord(3) = b(3);
            do dim =1,3
                sign_Cor(dim) = sign(1,new_Cord(dim)-old_Cord(dim))
            end do

            do k = old_Cord(3)+sign_Cor(3), new_Cord(3),sign_Cor(3)
                do j = old_Cord(2)+sign_Cor(2), new_Cord(2),sign_Cor(2)
                    do i = old_Cord(1)+sign_Cor(1), new_Cord(1),sign_Cor(1)
                        x1 = (x(i-sign_Cor(1),j,k)+x(i,j-sign_Cor(2),k)+x(i,j,k-sign_Cor(3)))/3_dp
                        x2 = (x(i,j-sign_Cor(2),k-sign_Cor(3))+x(i-sign_Cor(1),j,k-sign_Cor(3))+x(i-sign_Cor(1),j-sign_Cor(2),k))/3_dp
                        x(i,j,k) = 2_dp*x1-x2
                    end do
                end do
            end do
            ! (end,end,start)
            old_Cord(1) = ie  ; old_Cord(2) = je  ; old_Cord(3) = kb  ; 
            new_Cord(1) = e(1); new_Cord(2) = e(2); new_Cord(3) = b(3);
            do dim =1,3
                sign_Cor(dim) = sign(1,new_Cord(dim)-old_Cord(dim))
            end do

            do k = old_Cord(3)+sign_Cor(3), new_Cord(3),sign_Cor(3)
                do j = old_Cord(2)+sign_Cor(2), new_Cord(2),sign_Cor(2)
                    do i = old_Cord(1)+sign_Cor(1), new_Cord(1),sign_Cor(1)
                        x1 = (x(i-sign_Cor(1),j,k)+x(i,j-sign_Cor(2),k)+x(i,j,k-sign_Cor(3)))/3_dp
                        x2 = (x(i,j-sign_Cor(2),k-sign_Cor(3))+x(i-sign_Cor(1),j,k-sign_Cor(3))+x(i-sign_Cor(1),j-sign_Cor(2),k))/3_dp
                        x(i,j,k) = 2_dp*x1-x2
                    end do
                end do
            end do
            ! (start,start,end)
            old_Cord(1) = ib  ; old_Cord(2) = jb  ; old_Cord(3) = ke  ; 
            new_Cord(1) = b(1); new_Cord(2) = b(2); new_Cord(3) = e(3);
            do dim =1,3
                sign_Cor(dim) = sign(1,new_Cord(dim)-old_Cord(dim))
            end do

            do k = old_Cord(3)+sign_Cor(3), new_Cord(3),sign_Cor(3)
                do j = old_Cord(2)+sign_Cor(2), new_Cord(2),sign_Cor(2)
                    do i = old_Cord(1)+sign_Cor(1), new_Cord(1),sign_Cor(1)
                        x1 = (x(i-sign_Cor(1),j,k)+x(i,j-sign_Cor(2),k)+x(i,j,k-sign_Cor(3)))/3_dp
                        x2 = (x(i,j-sign_Cor(2),k-sign_Cor(3))+x(i-sign_Cor(1),j,k-sign_Cor(3))+x(i-sign_Cor(1),j-sign_Cor(2),k))/3_dp
                        x(i,j,k) = 2_dp*x1-x2
                    end do
                end do
            end do
            ! (start,end,end)
            old_Cord(1) = ib  ; old_Cord(2) = je  ; old_Cord(3) = ke  ; 
            new_Cord(1) = b(1); new_Cord(2) = e(2); new_Cord(3) = e(3);
            do dim =1,3
                sign_Cor(dim) = sign(1,new_Cord(dim)-old_Cord(dim))
            end do

            do k = old_Cord(3)+sign_Cor(3), new_Cord(3),sign_Cor(3)
                do j = old_Cord(2)+sign_Cor(2), new_Cord(2),sign_Cor(2)
                    do i = old_Cord(1)+sign_Cor(1), new_Cord(1),sign_Cor(1)
                        x1 = (x(i-sign_Cor(1),j,k)+x(i,j-sign_Cor(2),k)+x(i,j,k-sign_Cor(3)))/3_dp
                        x2 = (x(i,j-sign_Cor(2),k-sign_Cor(3))+x(i-sign_Cor(1),j,k-sign_Cor(3))+x(i-sign_Cor(1),j-sign_Cor(2),k))/3_dp
                        x(i,j,k) = 2_dp*x1-x2
                    end do
                end do
            end do
            ! (end,start,end)
            old_Cord(1) = ie  ; old_Cord(2) = jb  ; old_Cord(3) = ke  ; 
            new_Cord(1) = e(1); new_Cord(2) = b(2); new_Cord(3) = e(3);
            do dim =1,3
                sign_Cor(dim) = sign(1,new_Cord(dim)-old_Cord(dim))
            end do

            do k = old_Cord(3)+sign_Cor(3), new_Cord(3),sign_Cor(3)
                do j = old_Cord(2)+sign_Cor(2), new_Cord(2),sign_Cor(2)
                    do i = old_Cord(1)+sign_Cor(1), new_Cord(1),sign_Cor(1)
                        x1 = (x(i-sign_Cor(1),j,k)+x(i,j-sign_Cor(2),k)+x(i,j,k-sign_Cor(3)))/3_dp
                        x2 = (x(i,j-sign_Cor(2),k-sign_Cor(3))+x(i-sign_Cor(1),j,k-sign_Cor(3))+x(i-sign_Cor(1),j-sign_Cor(2),k))/3_dp
                        x(i,j,k) = 2_dp*x1-x2
                    end do
                end do
            end do
            ! (end,end,end)
            old_Cord(1) = ie  ; old_Cord(2) = je  ; old_Cord(3) = ke  ; 
            new_Cord(1) = e(1); new_Cord(2) = e(2); new_Cord(3) = e(3);
            do dim =1,3
                sign_Cor(dim) = sign(1,new_Cord(dim)-old_Cord(dim))
            end do

            do k = old_Cord(3)+sign_Cor(3), new_Cord(3),sign_Cor(3)
                do j = old_Cord(2)+sign_Cor(2), new_Cord(2),sign_Cor(2)
                    do i = old_Cord(1)+sign_Cor(1), new_Cord(1),sign_Cor(1)
                        x1 = (x(i-sign_Cor(1),j,k)+x(i,j-sign_Cor(2),k)+x(i,j,k-sign_Cor(3)))/3_dp
                        x2 = (x(i,j-sign_Cor(2),k-sign_Cor(3))+x(i-sign_Cor(1),j,k-sign_Cor(3))+x(i-sign_Cor(1),j-sign_Cor(2),k))/3_dp
                        x(i,j,k) = 2_dp*x1-x2
                    end do
                end do
            end do

        end subroutine get_Corner_Value

        ! 网格块计算几何量函数
        subroutine compute_Geometry(B)
            implicit none
           ! 输入/输出参数
            type(block_Type),pointer :: B

           ! 中间变量
            integer                  :: i,j,k       ! 索引变量
            real(dp),dimension(3)    :: p1,p2,p3,p4 ! 点的位矢量
            real(dp),dimension(3)    :: p5,p6,p7,p8 ! 点的位矢量
            real(dp),dimension(3)    :: S           ! 面矢量
            real(dp),dimension(3)    :: S1,S2,S3    ! 面矢量
            real(dp),dimension(3)    :: S4,S5,S6    ! 面矢量
           ! 计算i向面对相关参数
            do k = 1-LAP, B%nz-1+LAP
                do j = 1-LAP, B%ny-1+LAP
                    do i = 1-LAP, B%nx+LAP
                        ! 点指定：p1->(i,j,k);p2->(i,j+1,k);p3->(i,j,k+1);p4->(i,j+1,k+1)
                        p1(1) = B%x(i,j,k)    ; p1(2) = B%y(i,j,k)    ; p1(3) = B%z(i,j,k)  ;
                        p2(1) = B%x(i,j+1,k)  ; p2(2) = B%y(i,j+1,k)  ; p2(3) = B%z(i,j+1,k)  ;
                        p3(1) = B%x(i,j,k+1)  ; p3(2) = B%y(i,j,k+1)  ; p3(3) = B%z(i,j,k+1)  ;
                        p4(1) = B%x(i,j+1,k+1); p4(2) = B%y(i,j+1,k+1); p4(3) = B%z(i,j+1,k+1);
                        ! 计算面矢量
                        S = compute_S(p1,p2,p3,p4)
                        ! 面矢量指定
                        B%ni1(i,j,k) = S(1)
                        B%ni2(i,j,k) = S(2)
                        B%ni3(i,j,k) = S(3)
                        B%Si(i,j,k)  = sqrt(S(1)**2+S(2)**2+S(3)**2)
                    end do
                end do
            end do
           ! 计算j向面的相关参数
            do k = 1-LAP, B%nz-1+LAP
                do j = 1-LAP, B%ny+LAP
                    do i = 1-LAP, B%nx-1+LAP
                        ! 点指定：p1->(i,j,k);p2->(i,j,k+1);p3->(i+1,j,k);p4->(i+1,j,k+1)
                        p1(1) = B%x(i,j,k)    ; p1(2) = B%y(i,j,k)    ; p1(3) = B%z(i,j,k)  ;
                        p2(1) = B%x(i,j,k+1)  ; p2(2) = B%y(i,j,k+1)  ; p2(3) = B%z(i,j,k+1)  ;
                        p3(1) = B%x(i+1,j,k)  ; p3(2) = B%y(i+1,j,k)  ; p3(3) = B%z(i+1,j,k)  ;
                        p4(1) = B%x(i+1,j,k+1); p4(2) = B%y(i+1,j,k+1); p4(3) = B%z(i+1,j,k+1);
                        ! 计算面矢量
                        S = compute_S(p1,p2,p3,p4)
                        ! 面矢量指定
                        B%nj1(i,j,k) = S(1)
                        B%nj2(i,j,k) = S(2)
                        B%nj3(i,j,k) = S(3)
                        B%Sj(i,j,k)  = sqrt(S(1)**2+S(2)**2+S(3)**2)
                    end do
                end do
            end do
           ! 计算k向面的相关参数
            do k = 1-LAP, B%nz+LAP
                do j = 1-LAP, B%ny-1+LAP
                    do i = 1-LAP, B%nx-1+LAP
                        ! 点指定：p1->(i,j,k);p2->(i+1,j,k);p3->(i,j+1,k);p4->(i+1,j+1,k)
                        p1(1) = B%x(i,j,k)    ; p1(2) = B%y(i,j,k)    ; p1(3) = B%z(i,j,k)  ;
                        p2(1) = B%x(i+1,j,k)  ; p2(2) = B%y(i+1,j,k)  ; p2(3) = B%z(i+1,j,k)  ;
                        p3(1) = B%x(i,j+1,k)  ; p3(2) = B%y(i,j+1,k)  ; p3(3) = B%z(i,j+1,k)  ;
                        p4(1) = B%x(i+1,j+1,k); p4(2) = B%y(i+1,j+1,k); p4(3) = B%z(i+1,j+1,k);
                        ! 计算面矢量
                        S = compute_S(p1,p2,p3,p4)
                        ! 面矢量指定
                        B%nk1(i,j,k) = S(1)
                        B%nk2(i,j,k) = S(2)
                        B%nk3(i,j,k) = S(3)
                        B%Sk(i,j,k)  = sqrt(S(1)**2+S(2)**2+S(3)**2)
                    end do
                end do
            end do
           ! 计算单元体体积
            do k = 1-LAP,B%nz-1+LAP
                do j = 1-LAP,B%ny-1+LAP
                    do i = 1-LAP,B%nx-1+LAP
                        ! 指定点:p1->(i,j,k)  ;p2->(i+1,j,k)  ;p3->(i,j+1,k)   ;p4->(i+1,j+1,k)
                        !       p5->(i,j,k+1);p6->(i+1,j,k+1) ;p7->(i,j+1,k+1);p8->(i+1,j+1,k+1)
                        p1(1) = B%x(i,j,k)      ;p1(2) = B%y(i,j,k)      ;p1(3) = B%z(i,j,k)      ;
                        p2(1) = B%x(i+1,j,k)    ;p2(2) = B%y(i+1,j,k)    ;p2(3) = B%z(i+1,j,k)    ;
                        p3(1) = B%x(i,j+1,k)    ;p3(2) = B%y(i,j+1,k)    ;p3(3) = B%z(i,j+1,k)    ;
                        p4(1) = B%x(i+1,j+1,k)  ;p4(2) = B%y(i+1,j+1,k)  ;p4(3) = B%z(i+1,j+1,k)  ;
                        p5(1) = B%x(i,j,k+1)    ;p5(2) = B%y(i,j,k+1)    ;p5(3) = B%z(i,j,k+1)    ;
                        p6(1) = B%x(i+1,j,k+1)  ;p6(2) = B%y(i+1,j,k+1)  ;p6(3) = B%z(i+1,j,k+1)  ;
                        p7(1) = B%x(i,j+1,k+1)  ;p7(2) = B%y(i,j+1,k+1)  ;p7(3) = B%z(i,j+1,k+1)  ;
                        p8(1) = B%x(i+1,j+1,k+1);p8(2) = B%y(i+1,j+1,k+1);p8(3) = B%z(i+1,j+1,k+1);
                        ! 指点面:S1->Si(i,j,k);S2->Sj(i,j,k);S3->Sk(i,j,k);
                        !        S4 ->Si(i+1,j,k); S5 ->Sj(i,j+1,k); S6 ->Sk(i,j,k+1)
                        S1(1) = B%ni1(i,j,k)  ; S1(2) = B%ni2(i,j,k)   ; S1(3) = B%ni3(i,j,k)  ;
                        S2(1) = B%nj1(i,j,k)  ; S2(2) = B%nj2(i,j,k)   ; S2(3) = B%nj3(i,j,k)  ;
                        S3(1) = B%nk1(i,j,k)  ; S3(2) = B%nk2(i,j,k)   ; S3(3) = B%nk3(i,j,k)  ;
                        S4(1) = B%ni1(i+1,j,k); S4(2) = B%ni2(i+1,j,k) ; S4(3) = B%ni3(i+1,j,k);
                        S5(1) = B%nj1(i,j+1,k); S5(2) = B%nj2(i,j+1,k) ; S5(3) = B%nj3(i,j+1,k);
                        S6(1) = B%nk1(i,j,k+1); S6(2) = B%nk2(i,j,k+1) ; S6(3) = B%nk3(i,j,k+1) 
                        ! 计算体积
                        B%vol(i,j,k) = compute_Volume2(p1,p2,p3,p4,p5,p6,p7,p8,S1,S2,S3,S4,S5,S6)
                        write(*,*) compute_Volume2(p1,p2,p3,p4,p5,p6,p7,p8,S1,S2,S3,S4,S5,S6)
                        write(*,*) compute_Volume1(p1,p2,p3,p4,p5,p6,p7,p8)
                    end do
                end do
            end do
        end subroutine compute_Geometry


        ! 计算控制面法向量函数
        function compute_S(p1,p2,p3,p4)
           ! 注释说明：1->4向量；2->3向量；（1->4)*(2->3)
            implicit none
            ! 输入/输出参数
            real(dp),intent(in),dimension(3)    :: p1,p2,p3,p4      ! 面上4点
            real(dp),dimension(3)               :: compute_S        ! 面的法向量

            ! 中间参数
            real(dp),dimension(3)       :: dr14,dr23

            ! 计算 1->4 和 2->3 的矢量值
            dr14 = p4-p1; dr23 = p3-p2

            ! 计算n1,n2,n3
            compute_S(1) = (dr14(2)*dr23(3)-dr14(3)*dr23(2))/2
            compute_S(2) = (dr14(3)*dr23(1)-dr14(1)*dr23(3))/2
            compute_S(3) = (dr14(1)*dr23(2)-dr14(2)*dr23(1))/2
        end function compute_S

        ! 计算单元体体积方法1：散度定理
        ! 函数1：未知单元面矢量
        function compute_Volume1(p1,p2,p3,p4,p5,p6,p7,p8)
           ! 注释：1->(i,j,k);2->(i+1,j,k);3->(i,j+1,k);4->(i+1,j+1,k)
           !      5->(i,j,k+1);6->(i+1,j,k+1);7->(i,j+1,k+1);8->(i+1,j+1,k+1)
            implicit none
           ! 输入/输出参数
            real(dp),intent(in),dimension(3)         :: p1,p2,p3,p4,p5,p6,p7,p8
            real(dp)                                 :: compute_Volume1

           ! 中间参数
            real(dp),dimension(3)       :: S1,S2,S3,S4,S5,S6
            real(dp),dimension(3)       :: r1,r2,r3,r4,r5,r6

           ! 计算面同向矢量
            ! 计算i-面
            S1 = compute_S(p1,p3,p5,p7)
            ! 计算j-面
            S2 = compute_S(p1,p5,p2,p6)
            ! 计算k-面
            S3 = compute_S(p1,p2,p3,p4)
            ! 计算i+面
            S4 = compute_S(p2,p4,p6,p8)
            ! 计算j+面
            S5 = compute_S(p3,p7,p4,p8)
            ! 计算k+面
            S6 = compute_S(p5,p6,p7,p8)

           ! 计算面中心矢量
            r1 = (p1+p3+p5+p7)/4
            r2 = (p1+p5+p2+p6)/4
            r3 = (p1+p2+p3+p4)/4
            r4 = (p2+p4+p6+p8)/4
            r5 = (p3+p7+p4+p8)/4
            r6 = (p5+p6+p7+p8)/4

           ! 计算方法
            compute_Volume1 = (-(dot_product(S1,r1)+dot_product(S2,r2)+dot_product(S3,r3))&
                             &+(dot_product(S4,r4)+dot_product(S5,r5)+dot_product(S6,r6)))/3
            
        end function compute_Volume1
        ! 函数2：已知单元面矢量
        function compute_Volume2(p1,p2,p3,p4,p5,p6,p7,p8,S1,S2,S3,S4,S5,S6)
            implicit none
            ! 输入/输出参数
            real(dp),intent(in),dimension(3)         :: p1,p2,p3,p4,p5,p6,p7,p8
            real(dp),intent(in),dimension(3)         :: S1,S2,S3,S4,S5,S6
            real(dp)                                 :: compute_Volume2

           ! 中间变量
            real(dp),dimension(3)       :: r1,r2,r3,r4,r5,r6

           ! 计算面中心矢量
            r1 = (p1+p3+p5+p7)/4
            r2 = (p1+p5+p2+p6)/4
            r3 = (p1+p2+p3+p4)/4
            r4 = (p2+p4+p6+p8)/4
            r5 = (p3+p7+p4+p8)/4
            r6 = (p5+p6+p7+p8)/4

           ! 计算方法
            compute_Volume2 = ((dot_product(S1,r1)+dot_product(S2,r2)+dot_product(S3,r3))&
                             &-(dot_product(S4,r4)+dot_product(S5,r5)+dot_product(S6,r6)))/3
        end function compute_Volume2


    end module read_Mesh
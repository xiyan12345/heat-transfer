!--------------------------------------------------------------------------------------------------
! ��ȡ����ڵ�����Ӳ����ļ�
!--------------------------------------------------------------------------------------------------
! ��������������㼸����ģ��
    module read_Mesh
        use const_Var
        use global_Var
        implicit none
        
    contains
        subroutine read_Mesh_2D()
            integer                     :: m            ! ����������
            integer                     :: ksub         ! ������������
            integer                     :: i,j,k        ! ������� 
            integer,dimension(3)        :: ib,ie,ib1,ie1  ! �����С
            logical                     :: fexist       ! �ļ������߼���
            real(dp)                    :: dz           ! z������ʱ�ļ��
            type(block_Type),pointer    :: B            ! ��ָ��
            type(BC_MSG_Type),pointer   :: BC           ! ��߽�ָ��

           ! ��������ڵ��ļ�
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
                B%nz = zStretch     ! z������
                dz   = 10           ! z������ʱ�ļ��

               ! �ռ����
                ! �ڵ��Ͳ�������
                allocate(B%x(1-LAP:B%nx+LAP,1-LAP:B%ny+LAP,1-LAP:B%nz+LAP))
                allocate(B%y(1-LAP:B%nx+LAP,1-LAP:B%ny+LAP,1-LAP:B%nz+LAP))
                allocate(B%z(1-LAP:B%nx+LAP,1-LAP:B%ny+LAP,1-LAP:B%nz+LAP))

                ! i�����������
                allocate(B%Si(1-LAP:B%nx+LAP,1-LAP:B%ny+LAP-1,1-LAP:B%nz+LAP-1))
                allocate(B%ni1(1-LAP:B%nx+LAP,1-LAP:B%ny+LAP-1,1-LAP:B%nz+LAP-1))
                allocate(B%ni2(1-LAP:B%nx+LAP,1-LAP:B%ny+LAP-1,1-LAP:B%nz+LAP-1))
                allocate(B%ni3(1-LAP:B%nx+LAP,1-LAP:B%ny+LAP-1,1-LAP:B%nz+LAP-1))

                ! j�����������
                allocate(B%Sj(1-LAP:B%nx+LAP-1,1-LAP:B%ny+LAP,1-LAP:B%nz+LAP-1))
                allocate(B%nj1(1-LAP:B%nx+LAP-1,1-LAP:B%ny+LAP,1-LAP:B%nz+LAP-1))
                allocate(B%nj2(1-LAP:B%nx+LAP-1,1-LAP:B%ny+LAP,1-LAP:B%nz+LAP-1))
                allocate(B%nj3(1-LAP:B%nx+LAP-1,1-LAP:B%ny+LAP,1-LAP:B%nz+LAP-1))

                ! k�����������
                allocate(B%Sk(1-LAP:B%nx+LAP-1,1-LAP:B%ny+LAP-1,1-LAP:B%nz+LAP))
                allocate(B%nk1(1-LAP:B%nx+LAP-1,1-LAP:B%ny+LAP-1,1-LAP:B%nz+LAP))
                allocate(B%nk2(1-LAP:B%nx+LAP-1,1-LAP:B%ny+LAP-1,1-LAP:B%nz+LAP))
                allocate(B%nk3(1-LAP:B%nx+LAP-1,1-LAP:B%ny+LAP-1,1-LAP:B%nz+LAP))

                ! ��Ԫ��������
                allocate(B%vol(1-LAP:B%nx+LAP-1,1-LAP:B%ny+LAP-1,1-LAP:B%nz+LAP-1))
                allocate(B%U(nVar,1-LAP:B%nx+LAP-1,1-LAP:B%ny+LAP-1,1-LAP:B%nz+LAP-1))
                allocate(B%Un(nVar,1-LAP:B%nx+LAP-1,1-LAP:B%ny+LAP-1,1-LAP:B%nz+LAP-1))
                allocate(B%res(nVar,1-LAP:B%nx+LAP-1,1-LAP:B%ny+LAP-1,1-LAP:B%nz+LAP-1))

               ! ��������ڵ�����
                read(un_mesh,*) (((B%x(i,j,k),i=1,B%nx),j=1,B%ny),k=1,1)
                read(un_mesh,*) (((B%y(i,j,k),i=1,B%nx),j=1,B%ny),k=1,1)
                read(un_mesh,*) (((B%z(i,j,k),i=1,B%nx),j=1,B%ny),k=1,1)

               ! 2D -> 3D��z������
                do k = 2, zStretch
                    B%x(:,:,k) = B%x(:,:,1)
                    B%y(:,:,k) = B%y(:,:,1)
                    B%z(:,:,k) = B%z(:,:,1)+dz*(k-1)
                end do
            end do
            close(un_mesh)
           ! �������ӹ�ϵ�ļ�
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
                B%subface = B%subface+2                         ! ���������Գ���
                allocate(B%BC_MSG(B%subface))

                do ksub = 1, B%subface-2                        ! �����ļ��е�����
                    BC => B%BC_MSG(ksub)
                    read(un_inp,*) ib(1),ie(1),ib(2),ie(2),BC%bc
                    ib(3) = -1; ie(3) = -1*zStretch;
                    if (BC%bc < 0 ) then
                        read(un_inp,*) ib1(1),ie1(1),ib1(2),ie1(2),BC%nb1
                        ib1(3) = -1; ie1(3) = -1*zStretch;
                    end if
                    call convert_BC(BC,ib,ie,ib1,ie1)           ! ȷ������������������Ϣ
                end do

                do ksub = B%subface-1, B%subface                ! ���������Գ����������Ϣ����
                    BC => B%BC_MSG(ksub)
                    ib(1) = 1; ie(1) = B%nx; ib(2) = 1; ie(2) = B%ny;
                    if (ksub == B%subface-1) then
                        ib(3) = 1; ie(3) = 1
                    else 
                        ib(3) = zStretch; ie(3) = zStretch
                    end if
                    BC%bc = 3
                    call convert_BC(BC,ib,ie,ib1,ie1)           ! ȷ������������������Ϣ 
                end do

            end do
            close(un_inp)
           ! ���ÿһ�飬��������������
            do m = 1, num_Block
                B => mesh(m)
                call upgrade_ghost(B,LAP)
            end do
           ! ���ÿһ�飬���㼸����
            do m = 1, num_Block
                B => mesh(m)
                call compute_Geometry(B)
            end do

        end subroutine read_Mesh_2D

       ! ȷ�����������������Ϣ����
        subroutine convert_BC(BC,ib,ie,ib1,ie1)
            implicit none
           ! ����/�������
            integer,dimension(:)        :: ib,ie 
            integer,dimension(:)        :: ib1,ie1
            type(BC_MSG_Type)           :: BC

           ! ��̬����
            integer                     :: dim          ! ά������
            integer                     :: linkDim      ! ������ά������
            integer,dimension(3)        :: owner        ! ����ά�ȱ�ʶ
            integer,dimension(3)        :: conne        ! ����ά�ȱ�ʶ

           ! �����������
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

           ! �������������
            if (BC%bc > 0) then
            ! ����߽磺��������Ϣȫ����0
                BC%ib1 = 0; BC%jb1 = 0; BC%kb1 = 0;
                BC%ie1 = 0; BC%je1 = 0; BC%ke1 = 0;
                BC%nb1 = 0; BC%face1 = 0
            else 
            ! �ڱ߽磺������������Ϣ
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

                ! ����˳�����ӷ�
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

       ! �������������꺯��
        subroutine upgrade_ghost(B,LAP)
            implicit none
           ! ����/�������
            type(block_Type),pointer        :: B
            integer                         :: LAP

           ! ��̬����
            integer                         :: n            ! �������������
            integer                         :: ksub         ! ��������
            integer                         :: dim_Link     ! ������ά�ȼ�����������
            integer                         :: dim_Link1    ! ��������ά�ȼ�����������
            integer                         :: i,j,k        ! ����ڵ�����
            integer                         :: i1,j1,k1     ! ��������ڵ�����
            integer                         :: dim          ! ά������
            integer,dimension(3)            :: faceb        ! ��������ʼ����
            integer,dimension(3)            :: facee        ! ��������ֹ����
            integer,dimension(3)            :: faceb1       ! ����������ʼ����
            integer,dimension(3)            :: facee1       ! ����������ֹ����
            integer,dimension(3)            :: owner_O      ! ��������ԭ��
            integer,dimension(3)            :: owner_D      ! ������������
            integer,dimension(3)            :: conne_O      ! ����������ԭ��
            integer,dimension(3)            :: conne_D      ! ��������������
            type(BC_MSG_Type),pointer       :: BC           ! ����ָ��
            type(block_Type) ,pointer       :: B_Conne      ! �����ӿ�ָ��


           ! ѭ������������
            do n = 1, LAP

               ! �����ǽǵ���������������
                do ksub = 1, B%subface
                    BC => B%BC_MSG(ksub)
                    if ( BC%bc >= 0) then
                    ! ����߽�
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
                    ! �ڱ߽�
                        B_Conne => mesh(BC%nb1)
                       ! ������������Χ
                        faceb(1) = BC%ib; faceb(2) = BC%jb; faceb(3) = BC%kb;
                        facee(1) = BC%ie; facee(2) = BC%je; facee(3) = BC%ke;
                        dim_Link = sign(1,BC%face-4)*(mod(BC%face-1,3)+1)
                        faceb(abs(dim_Link)) = faceb(abs(dim_Link))+sign(1,dim_Link)*n
                        facee(abs(dim_Link)) = faceb(abs(dim_Link))
                       ! ������������Χ
                        faceb1(1) = BC%ib1; faceb1(2) = BC%jb1; faceb1(3) = BC%kb1; 
                        facee1(1) = BC%ie1; facee1(2) = BC%je1; facee1(3) = BC%ke1;
                        dim_Link1 = sign(1,BC%face1-4)*(mod(BC%face1-1,3)+1)
                        faceb1(abs(dim_Link1)) = faceb1(abs(dim_Link1))-sign(1,dim_Link1)*n
                        facee1(abs(dim_Link1)) = faceb1(abs(dim_Link1))

                        ! ����ԭ�㶨λ
                        owner_O(:) = faceb(:)
                        ! ��������ԭ�㶨λ
                        do dim = 1, 3
                            if ( BC%L(dim)>0 ) then
                                conne_O(dim) = faceb1(dim)
                            else 
                                conne_O(dim) = facee1(dim)
                            end if
                        end do

                        ! �໥��ֵ
                        do k = faceb(3), facee(3)
                            do j = faceb(2), facee(2)
                                do i = faceb(1), facee(1)
                                    ! �������
                                    owner_D(1) = i -faceb(1)
                                    owner_D(2) = j -faceb(2)
                                    owner_D(3) = k -faceb(3)
                                    ! ���������
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

               ! �����ǵ���������������
               ! 12�����ֵ
                ! x��4�����ֵ
                call get_Corner_Value(B%x,1,1,1,B%nx,B%ny,B%nz,n)
                call get_Corner_Value(B%y,1,1,1,B%nx,B%ny,B%nz,n)
                call get_Corner_Value(B%z,1,1,1,B%nx,B%ny,B%nz,n)

                    
                
                
            end do
        end subroutine upgrade_ghost

        ! �ǲ���ֵ����
        subroutine get_Corner_Value(x,ib,jb,kb,ie,je,ke,n)
            implicit none
           ! ����/�������
            integer                                                         :: n
            integer                                                         :: ib,jb,kb     ! ԭʼ�����
            integer                                                         :: ie,je,ke     ! ԭʼ���յ� 
            real(dp),dimension(ib-LAP:ie+LAP,jb-LAP:je+LAP,kb-LAP:ke+LAP)   :: x

            
           ! �м����
            integer                     :: i,j,k        ! �ڵ�����
            integer                     :: p            ! ��������
            integer                     :: dim          ! ά������
            integer,dimension(3)        :: sign_Cor     ! �ǲ�����ά�ȷ��ű���
            integer,dimension(3)        :: old_Cord     ! �ǵ�����ɽǵ�����
            integer,dimension(3)        :: new_Cord     ! �ǵ������½ǵ�����
            integer,dimension(3)        :: b,e          ! �½ǵ�����
            real(dp)                    :: x1,x2        ! ������̬����
            


            b(1) = ib-n; b(2) = jb-n; b(3) = kb-n;
            e(1) = ie+n; e(2) = je+n; e(3) = ke+n;

            do p = 1, n
               ! ��x��4����
                ! start:���������end��������ӣ���old,new)��(new,old)
                ! (start,start)
                x(ib:ie,jb-p,b(3)) = x(ib:ie,jb-p+1,b(3))+x(ib:ie,jb-p,b(3)+1)-x(ib:ie,jb-p+1,b(3)+1)
                x(ib:ie,b(2),kb-p) = x(ib:ie,b(2)+1,kb-p)+x(ib:ie,b(2),kb-p+1)-x(ib:ie,b(2)+1,kb-p+1)
                ! (start,end��
                x(ib:ie,jb-p,e(3)) = x(ib:ie,jb-p+1,e(3))+x(ib:ie,jb-p,e(3)-1)-x(ib:ie,jb-p+1,e(3)-1)
                x(ib:ie,b(2),ke+p) = x(ib:ie,b(2)+1,ke+p)+x(ib:ie,b(2),ke+p-1)-x(ib:ie,b(2)+1,ke+p-1)
                ! (end,start)
                x(ib:ie,je+p,b(3)) = x(ib:ie,je+p-1,b(3))+x(ib:ie,je+p,b(3)+1)-x(ib:ie,je+p-1,b(3)+1)
                x(ib:ie,e(2),kb-p) = x(ib:ie,e(2)-1,kb-p)+x(ib:ie,e(2),kb-p+1)-x(ib:ie,e(2)-1,kb-p+1)
                ! (end,end)
                x(ib:ie,je+p,e(3)) = x(ib:ie,je+p-1,e(3))+x(ib:ie,je+p,e(3)-1)-x(ib:ie,je+p-1,e(3)-1)
                x(ib:ie,e(2),ke+p) = x(ib:ie,e(2)-1,ke+p)+x(ib:ie,e(2),ke+p-1)-x(ib:ie,e(2)-1,ke+p-1)

               ! ��y��4����
                ! start:���������end��������ӣ���old,new)��(new,old)
                ! (start,start)��problem)
                x(ib-p,jb:je,b(3)) = x(ib-p+1,jb:je,b(3))+x(ib-p,jb:je,b(3)+1)-x(ib-p+1,jb:je,b(3)+1)
                x(b(1),jb:je,kb-p) = x(b(1)+1,jb:je,kb-p)+x(b(1),jb:je,kb-p+1)-x(b(1)+1,jb:je,kb-p+1)
                ! (start,end��
                x(ib-p,jb:je,e(3)) = x(ib-p+1,jb:je,e(3))+x(ib-p,jb:je,e(3)-1)-x(ib-p+1,jb:je,e(3)-1)
                x(b(1),jb:je,ke+p) = x(b(1)+1,jb:je,ke+p)+x(b(1),jb:je,ke+p-1)-x(b(1)+1,jb:je,ke+p-1)
                ! (end,start)
                x(ie+p,jb:je,b(3)) = x(ie+p-1,jb:je,b(3))+x(ie+p,jb:je,b(3)+1)-x(ie+p-1,jb:je,b(3)+1)
                x(e(1),jb:je,kb-p) = x(e(1)-1,jb:je,kb-p)+x(e(1),jb:je,kb-p+1)-x(e(1)-1,jb:je,kb-p+1)
                ! (end,end)(problem)
                x(ie+p,jb:je,e(3)) = x(ie+p-1,jb:je,e(3))+x(ie+p,jb:je,e(3)-1)-x(ie+p-1,jb:je,e(3)-1)
                x(e(1),jb:je,ke+p) = x(e(1)-1,jb:je,ke+p)+x(e(1),jb:je,ke+p-1)-x(e(1)-1,jb:je,ke+p-1)
            
               ! ��z��4����
                ! start:���������end��������ӣ���old,new)��(new,old)
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

           ! ��8���ǵ�
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

        ! �������㼸��������
        subroutine compute_Geometry(B)
            implicit none
           ! ����/�������
            type(block_Type),pointer :: B

           ! �м����
            integer                  :: i,j,k       ! ��������
            real(dp),dimension(3)    :: p1,p2,p3,p4 ! ���λʸ��
            real(dp),dimension(3)    :: p5,p6,p7,p8 ! ���λʸ��
            real(dp),dimension(3)    :: S           ! ��ʸ��
            real(dp),dimension(3)    :: S1,S2,S3    ! ��ʸ��
            real(dp),dimension(3)    :: S4,S5,S6    ! ��ʸ��
           ! ����i�������ز���
            do k = 1-LAP, B%nz-1+LAP
                do j = 1-LAP, B%ny-1+LAP
                    do i = 1-LAP, B%nx+LAP
                        ! ��ָ����p1->(i,j,k);p2->(i,j+1,k);p3->(i,j,k+1);p4->(i,j+1,k+1)
                        p1(1) = B%x(i,j,k)    ; p1(2) = B%y(i,j,k)    ; p1(3) = B%z(i,j,k)  ;
                        p2(1) = B%x(i,j+1,k)  ; p2(2) = B%y(i,j+1,k)  ; p2(3) = B%z(i,j+1,k)  ;
                        p3(1) = B%x(i,j,k+1)  ; p3(2) = B%y(i,j,k+1)  ; p3(3) = B%z(i,j,k+1)  ;
                        p4(1) = B%x(i,j+1,k+1); p4(2) = B%y(i,j+1,k+1); p4(3) = B%z(i,j+1,k+1);
                        ! ������ʸ��
                        S = compute_S(p1,p2,p3,p4)
                        ! ��ʸ��ָ��
                        B%ni1(i,j,k) = S(1)
                        B%ni2(i,j,k) = S(2)
                        B%ni3(i,j,k) = S(3)
                        B%Si(i,j,k)  = sqrt(S(1)**2+S(2)**2+S(3)**2)
                    end do
                end do
            end do
           ! ����j�������ز���
            do k = 1-LAP, B%nz-1+LAP
                do j = 1-LAP, B%ny+LAP
                    do i = 1-LAP, B%nx-1+LAP
                        ! ��ָ����p1->(i,j,k);p2->(i,j,k+1);p3->(i+1,j,k);p4->(i+1,j,k+1)
                        p1(1) = B%x(i,j,k)    ; p1(2) = B%y(i,j,k)    ; p1(3) = B%z(i,j,k)  ;
                        p2(1) = B%x(i,j,k+1)  ; p2(2) = B%y(i,j,k+1)  ; p2(3) = B%z(i,j,k+1)  ;
                        p3(1) = B%x(i+1,j,k)  ; p3(2) = B%y(i+1,j,k)  ; p3(3) = B%z(i+1,j,k)  ;
                        p4(1) = B%x(i+1,j,k+1); p4(2) = B%y(i+1,j,k+1); p4(3) = B%z(i+1,j,k+1);
                        ! ������ʸ��
                        S = compute_S(p1,p2,p3,p4)
                        ! ��ʸ��ָ��
                        B%nj1(i,j,k) = S(1)
                        B%nj2(i,j,k) = S(2)
                        B%nj3(i,j,k) = S(3)
                        B%Sj(i,j,k)  = sqrt(S(1)**2+S(2)**2+S(3)**2)
                    end do
                end do
            end do
           ! ����k�������ز���
            do k = 1-LAP, B%nz+LAP
                do j = 1-LAP, B%ny-1+LAP
                    do i = 1-LAP, B%nx-1+LAP
                        ! ��ָ����p1->(i,j,k);p2->(i+1,j,k);p3->(i,j+1,k);p4->(i+1,j+1,k)
                        p1(1) = B%x(i,j,k)    ; p1(2) = B%y(i,j,k)    ; p1(3) = B%z(i,j,k)  ;
                        p2(1) = B%x(i+1,j,k)  ; p2(2) = B%y(i+1,j,k)  ; p2(3) = B%z(i+1,j,k)  ;
                        p3(1) = B%x(i,j+1,k)  ; p3(2) = B%y(i,j+1,k)  ; p3(3) = B%z(i,j+1,k)  ;
                        p4(1) = B%x(i+1,j+1,k); p4(2) = B%y(i+1,j+1,k); p4(3) = B%z(i+1,j+1,k);
                        ! ������ʸ��
                        S = compute_S(p1,p2,p3,p4)
                        ! ��ʸ��ָ��
                        B%nk1(i,j,k) = S(1)
                        B%nk2(i,j,k) = S(2)
                        B%nk3(i,j,k) = S(3)
                        B%Sk(i,j,k)  = sqrt(S(1)**2+S(2)**2+S(3)**2)
                    end do
                end do
            end do
           ! ���㵥Ԫ�����
            do k = 1-LAP,B%nz-1+LAP
                do j = 1-LAP,B%ny-1+LAP
                    do i = 1-LAP,B%nx-1+LAP
                        ! ָ����:p1->(i,j,k)  ;p2->(i+1,j,k)  ;p3->(i,j+1,k)   ;p4->(i+1,j+1,k)
                        !       p5->(i,j,k+1);p6->(i+1,j,k+1) ;p7->(i,j+1,k+1);p8->(i+1,j+1,k+1)
                        p1(1) = B%x(i,j,k)      ;p1(2) = B%y(i,j,k)      ;p1(3) = B%z(i,j,k)      ;
                        p2(1) = B%x(i+1,j,k)    ;p2(2) = B%y(i+1,j,k)    ;p2(3) = B%z(i+1,j,k)    ;
                        p3(1) = B%x(i,j+1,k)    ;p3(2) = B%y(i,j+1,k)    ;p3(3) = B%z(i,j+1,k)    ;
                        p4(1) = B%x(i+1,j+1,k)  ;p4(2) = B%y(i+1,j+1,k)  ;p4(3) = B%z(i+1,j+1,k)  ;
                        p5(1) = B%x(i,j,k+1)    ;p5(2) = B%y(i,j,k+1)    ;p5(3) = B%z(i,j,k+1)    ;
                        p6(1) = B%x(i+1,j,k+1)  ;p6(2) = B%y(i+1,j,k+1)  ;p6(3) = B%z(i+1,j,k+1)  ;
                        p7(1) = B%x(i,j+1,k+1)  ;p7(2) = B%y(i,j+1,k+1)  ;p7(3) = B%z(i,j+1,k+1)  ;
                        p8(1) = B%x(i+1,j+1,k+1);p8(2) = B%y(i+1,j+1,k+1);p8(3) = B%z(i+1,j+1,k+1);
                        ! ָ����:S1->Si(i,j,k);S2->Sj(i,j,k);S3->Sk(i,j,k);
                        !        S4 ->Si(i+1,j,k); S5 ->Sj(i,j+1,k); S6 ->Sk(i,j,k+1)
                        S1(1) = B%ni1(i,j,k)  ; S1(2) = B%ni2(i,j,k)   ; S1(3) = B%ni3(i,j,k)  ;
                        S2(1) = B%nj1(i,j,k)  ; S2(2) = B%nj2(i,j,k)   ; S2(3) = B%nj3(i,j,k)  ;
                        S3(1) = B%nk1(i,j,k)  ; S3(2) = B%nk2(i,j,k)   ; S3(3) = B%nk3(i,j,k)  ;
                        S4(1) = B%ni1(i+1,j,k); S4(2) = B%ni2(i+1,j,k) ; S4(3) = B%ni3(i+1,j,k);
                        S5(1) = B%nj1(i,j+1,k); S5(2) = B%nj2(i,j+1,k) ; S5(3) = B%nj3(i,j+1,k);
                        S6(1) = B%nk1(i,j,k+1); S6(2) = B%nk2(i,j,k+1) ; S6(3) = B%nk3(i,j,k+1) 
                        ! �������
                        B%vol(i,j,k) = compute_Volume2(p1,p2,p3,p4,p5,p6,p7,p8,S1,S2,S3,S4,S5,S6)
                        write(*,*) compute_Volume2(p1,p2,p3,p4,p5,p6,p7,p8,S1,S2,S3,S4,S5,S6)
                        write(*,*) compute_Volume1(p1,p2,p3,p4,p5,p6,p7,p8)
                    end do
                end do
            end do
        end subroutine compute_Geometry


        ! ��������淨��������
        function compute_S(p1,p2,p3,p4)
           ! ע��˵����1->4������2->3��������1->4)*(2->3)
            implicit none
            ! ����/�������
            real(dp),intent(in),dimension(3)    :: p1,p2,p3,p4      ! ����4��
            real(dp),dimension(3)               :: compute_S        ! ��ķ�����

            ! �м����
            real(dp),dimension(3)       :: dr14,dr23

            ! ���� 1->4 �� 2->3 ��ʸ��ֵ
            dr14 = p4-p1; dr23 = p3-p2

            ! ����n1,n2,n3
            compute_S(1) = (dr14(2)*dr23(3)-dr14(3)*dr23(2))/2
            compute_S(2) = (dr14(3)*dr23(1)-dr14(1)*dr23(3))/2
            compute_S(3) = (dr14(1)*dr23(2)-dr14(2)*dr23(1))/2
        end function compute_S

        ! ���㵥Ԫ���������1��ɢ�ȶ���
        ! ����1��δ֪��Ԫ��ʸ��
        function compute_Volume1(p1,p2,p3,p4,p5,p6,p7,p8)
           ! ע�ͣ�1->(i,j,k);2->(i+1,j,k);3->(i,j+1,k);4->(i+1,j+1,k)
           !      5->(i,j,k+1);6->(i+1,j,k+1);7->(i,j+1,k+1);8->(i+1,j+1,k+1)
            implicit none
           ! ����/�������
            real(dp),intent(in),dimension(3)         :: p1,p2,p3,p4,p5,p6,p7,p8
            real(dp)                                 :: compute_Volume1

           ! �м����
            real(dp),dimension(3)       :: S1,S2,S3,S4,S5,S6
            real(dp),dimension(3)       :: r1,r2,r3,r4,r5,r6

           ! ������ͬ��ʸ��
            ! ����i-��
            S1 = compute_S(p1,p3,p5,p7)
            ! ����j-��
            S2 = compute_S(p1,p5,p2,p6)
            ! ����k-��
            S3 = compute_S(p1,p2,p3,p4)
            ! ����i+��
            S4 = compute_S(p2,p4,p6,p8)
            ! ����j+��
            S5 = compute_S(p3,p7,p4,p8)
            ! ����k+��
            S6 = compute_S(p5,p6,p7,p8)

           ! ����������ʸ��
            r1 = (p1+p3+p5+p7)/4
            r2 = (p1+p5+p2+p6)/4
            r3 = (p1+p2+p3+p4)/4
            r4 = (p2+p4+p6+p8)/4
            r5 = (p3+p7+p4+p8)/4
            r6 = (p5+p6+p7+p8)/4

           ! ���㷽��
            compute_Volume1 = (-(dot_product(S1,r1)+dot_product(S2,r2)+dot_product(S3,r3))&
                             &+(dot_product(S4,r4)+dot_product(S5,r5)+dot_product(S6,r6)))/3
            
        end function compute_Volume1
        ! ����2����֪��Ԫ��ʸ��
        function compute_Volume2(p1,p2,p3,p4,p5,p6,p7,p8,S1,S2,S3,S4,S5,S6)
            implicit none
            ! ����/�������
            real(dp),intent(in),dimension(3)         :: p1,p2,p3,p4,p5,p6,p7,p8
            real(dp),intent(in),dimension(3)         :: S1,S2,S3,S4,S5,S6
            real(dp)                                 :: compute_Volume2

           ! �м����
            real(dp),dimension(3)       :: r1,r2,r3,r4,r5,r6

           ! ����������ʸ��
            r1 = (p1+p3+p5+p7)/4
            r2 = (p1+p5+p2+p6)/4
            r3 = (p1+p2+p3+p4)/4
            r4 = (p2+p4+p6+p8)/4
            r5 = (p3+p7+p4+p8)/4
            r6 = (p5+p6+p7+p8)/4

           ! ���㷽��
            compute_Volume2 = ((dot_product(S1,r1)+dot_product(S2,r2)+dot_product(S3,r3))&
                             &-(dot_product(S4,r4)+dot_product(S5,r5)+dot_product(S6,r6)))/3
        end function compute_Volume2


    end module read_Mesh
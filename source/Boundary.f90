!--------------------------------------------------------------------------------------------------
! �߽���������ģ���ļ�
!--------------------------------------------------------------------------------------------------
! �߽���������ģ��
    module Boundary
        use global_Var
        use const_Var
        implicit none
        
    contains
        ! ���ڱ߽���Ϣ����
        subroutine exchange_Inner_Boundary()
            implicit none

            integer                     :: m            ! ������
            integer                     :: ksub         ! ��������
            integer                     :: dim          ! ά������
            integer                     :: i,j,k        ! ����ڵ�����
            integer                     :: i1,j1,k1     ! ��������ڵ�����
            integer                     :: dim_Link     ! ������ά�ȼ�����������
            integer                     :: dim_Link1    ! ��������ά�ȼ�����������
            integer,dimension(3)        :: faceb        ! ��������ʼ����
            integer,dimension(3)        :: facee        ! ��������ֹ����
            integer,dimension(3)        :: faceb1       ! ����������ʼ����
            integer,dimension(3)        :: facee1       ! ����������ֹ����
            integer,dimension(3)        :: owner_O      ! ��������ԭ��
            integer,dimension(3)        :: owner_D      ! ������������
            integer,dimension(3)        :: conne_O      ! ����������ԭ��
            integer,dimension(3)        :: conne_D      ! ��������������
            type(block_Type),pointer    :: B            ! ��ָ��
            type(block_Type),pointer    :: B_Conne      ! ���ӿ�ָ��
            type(BC_MSG_Type),pointer   :: BC           ! ����ָ��
            

            do m = 1, num_Block
                B => mesh(m)
                do ksub = 1, B%subface
                    BC => B%BC_MSG(ksub)
                    if ( BC%face>=0 ) cycle         ! ����߽���������
                   ! �ڱ߽�
                    B_Conne => mesh(BC%nb1)
                      ! ������������Χ
                    faceb(1) = BC%ib; faceb(2) = BC%jb; faceb(3) = BC%kb;
                    facee(1) = BC%ie; facee(2) = BC%je; facee(3) = BC%ke;
                    dim_Link = sign(1,BC%face-4)*(mod(BC%face-1,3)+1)
                    faceb(abs(dim_Link)) = faceb(abs(dim_Link))+sign(1,dim_Link)
                    facee(abs(dim_Link)) = faceb(abs(dim_Link))
                   ! ������������Χ
                    faceb1(1) = BC%ib1; faceb1(2) = BC%jb1; faceb1(3) = BC%kb1; 
                    facee1(1) = BC%ie1; facee1(2) = BC%je1; facee1(3) = BC%ke1;
                    dim_Link1 = sign(1,BC%face1-4)*(mod(BC%face1-1,3)+1)
                    faceb1(abs(dim_Link1)) = faceb1(abs(dim_Link1))-sign(1,dim_Link1)
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

                                B%U(:,i,j,k) = B_Conne%U(:,i1,j1,k1)
                            end do
                        end do 
                    end do
                end do
            end do
        end subroutine exchange_Inner_Boundary

        ! ���±߽�����
        subroutine upgrade_Boundary()
            implicit none

           ! �м����
            integer                     :: m        ! ������
            integer                     :: ksub     ! �߽�����
            type(block_Type) ,pointer   :: B        ! ��ָ��
            type(BC_MSG_Type),pointer   :: BC       ! ��߽�ָ��

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

        ! ���¶ԳƱ߽�����
        subroutine upgrade_Symmetry_Boundary(B,BC)
            implicit none 

           ! ����/�������
            type(block_Type) ,pointer           :: B        ! ��ָ��
            type(BC_MSG_Type),pointer           :: BC       ! �߽�ָ��

           ! �м����
            integer                             :: k        ! ά��ָʾ
            integer,dimension(3)                :: db,de    ! ���淶Χ
            integer,dimension(3)                :: db1,de1  ! ������Χ

            db(1) = BC%ib; db(2) = BC%jb; db(3) = BC%kb;
            de(1) = BC%ie; de(2) = BC%je; de(3) = BC%ke;
            k = sign(1,BC%face-4)*(mod(BC%face-1,3)+1)      ! k�ķ��Ŵ������������k�����ִ�������ά��

            if ( k>0 ) then
            else 
            end if
            
        end subroutine upgrade_Symmetry_Boundary

        ! ����ֱ�������߽�����
        subroutine upgrade_DirectQ_Boundary(B,BC)
            implicit none 

            ! ����/�������
            type(block_Type) ,pointer           :: B    ! ��ָ��
            type(BC_MSG_Type),pointer           :: BC   ! �߽�ָ��
        end subroutine upgrade_DirectQ_Boundary

        ! ���¶������ȱ߽�����
        subroutine upgrade_Convection_Boundary(B,BC)
            implicit none 

           ! ����/�������
            type(block_Type) ,pointer           :: B    ! ��ָ��
            type(BC_MSG_Type),pointer           :: BC   ! �߽�ָ��
        end subroutine upgrade_Convection_Boundary

        ! ����ֱ���¶ȱ߽�����
        subroutine upgrade_DirectT_Boundary(B,BC)
            implicit none 

           ! ����/�������
            type(block_Type) ,pointer           :: B    ! ��ָ��
            type(BC_MSG_Type),pointer           :: BC   ! �߽�ָ��
        end subroutine upgrade_DirectT_Boundary

        ! �����ڲ��߽�����
        subroutine upgrade_Inner_Boundary(B,BC)
            implicit none 

            ! ����/�������
            type(block_Type) ,pointer           :: B    ! ��ָ��
            type(BC_MSG_Type),pointer           :: BC   ! �߽�ָ��
        end subroutine upgrade_Inner_Boundary

    end module Boundary
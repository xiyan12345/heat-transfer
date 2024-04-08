!--------------------------------------------------------------------------------------------------
! �߽���������ģ���ļ�
!--------------------------------------------------------------------------------------------------
! �߽���������ģ��
    module Boundary
        use global_Var
        use const_Var
        implicit none
        
    contains

       ! ���±߽�������
        ! �������п�߽��������
        subroutine upgrade_Ghost_Cell()
            implicit none

           ! �м����
            integer                     :: m        ! ������
            integer                     :: ksub     ! �߽�����
            type(block_Type) ,pointer   :: B        ! ��ָ��
            type(BC_MSG_Type),pointer   :: BC       ! ��߽�ָ��

           ! ���·ǽǲ�����������
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
           ! ���½ǲ�������������
            do m = 1, num_Block
                B => mesh(m)
                call upgrade_Corner_Cell(B)
            end do

           ! �����ڱ߽紦���ǲ��������������
            do m = 1, num_Block
                B => mesh(m)
                do ksub = 1, B%subface
                    BC => B%BC_MSG(ksub)
                    if ( BC%bc<0 ) then         ! �ڱ߽����⴦��
                        call upgrade_Inner_Edge_Ghost_Cell(B,BC)
                    end if
                    
                end do
            end do

        end subroutine upgrade_Ghost_Cell

        ! ���¿�ԳƱ߽�ķǽǲ�����������
        subroutine  upgrade_Symmetry_Ghost_Cell(B,BC)
            implicit none 

           ! ����/�������
            type(block_Type) ,pointer           :: B        ! ��ָ��
            type(BC_MSG_Type),pointer           :: BC       ! �߽�ָ��

           ! �м����
            integer                             :: k        ! ����ά��ָʾ
            integer                             :: dim      ! ά������
            integer,dimension(3)                :: faceb    ! �ڲ�������ʼ����
            integer,dimension(3)                :: facee    ! �ڲ�������ֹ����
            integer,dimension(3)                :: faceb1   ! ������������ʼ����
            integer,dimension(3)                :: facee1   ! ������������ֹ����

           ! ������������Χ
            faceb(1) = BC%ib; faceb(2) = BC%jb; faceb(3) = BC%kb;
            facee(1) = BC%ie; facee(2) = BC%je; facee(3) = BC%ke;
            k = sign(1,BC%face-4)*(mod(BC%face-1,3)+1)
            do dim = 1, 3
                if ( dim == abs(k) ) then
                ! ����ָʾά�����⴦��
                faceb(dim) = faceb(dim)-(sign(1,BC%face-4)+1)/2
                facee(dim) = faceb(dim)
                else 
                ! ������ָʾά����������end�ڵ�-1
                facee(dim) = facee(dim)-1
                end if
            end do

            faceb1(:) = faceb(:); facee1(:) = facee(:);
            faceb1(abs(k)) = faceb1(abs(k))+sign(1,k)
            facee1(abs(k)) = faceb1(abs(k))

            B%U(:,faceb1(1):facee1(1),faceb1(2):facee1(2),faceb1(3):facee1(3)) =  B%U(:,faceb(1):facee(1),faceb(2):facee(2),faceb(3):facee(3))

        end subroutine  upgrade_Symmetry_Ghost_Cell

        ! ����ֱ�������߽�ķǽǲ�����������
        subroutine upgrade_DirectQ_Ghost_Cell(B,BC)
            implicit none 

            ! ����/�������
            type(block_Type) ,pointer           :: B    ! ��ָ��
            type(BC_MSG_Type),pointer           :: BC   ! �߽�ָ��

           ! �м����
            integer                             :: k        ! ����ά��ָʾ
            integer                             :: dim      ! ά������
            integer,dimension(3)                :: faceb    ! �ڲ�����1��ʼ����
            integer,dimension(3)                :: facee    ! �ڲ�����1��ֹ����
            integer,dimension(3)                :: faceb2    ! �ڲ�����2��ʼ����
            integer,dimension(3)                :: facee2    ! �ڲ�����2��ֹ����
            integer,dimension(3)                :: faceb1   ! ������������ʼ����
            integer,dimension(3)                :: facee1   ! ������������ֹ����

           ! ������������Χ
            faceb(1) = BC%ib; faceb(2) = BC%jb; faceb(3) = BC%kb;
            facee(1) = BC%ie; facee(2) = BC%je; facee(3) = BC%ke;
            k = sign(1,BC%face-4)*(mod(BC%face-1,3)+1)
            do dim = 1, 3
                if ( dim == abs(k) ) then
                ! ����ָʾά�����⴦��
                faceb(dim) = faceb(dim)-(sign(1,BC%face-4)+1)/2
                facee(dim) = faceb(dim)
                else 
                ! ������ָʾά����������end�ڵ�-1
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

        ! ���¶������ȱ߽�ķǽǲ�����������
        subroutine upgrade_Convection_Ghost_Cell(B,BC)
            implicit none 

           ! ����/�������
            type(block_Type) ,pointer           :: B    ! ��ָ��
            type(BC_MSG_Type),pointer           :: BC   ! �߽�ָ��

           ! �м����
            integer                             :: k        ! ����ά��ָʾ
            integer                             :: dim      ! ά������
            integer,dimension(3)                :: faceb    ! �ڲ�����1��ʼ����
            integer,dimension(3)                :: facee    ! �ڲ�����1��ֹ����
            integer,dimension(3)                :: faceb2    ! �ڲ�����2��ʼ����
            integer,dimension(3)                :: facee2    ! �ڲ�����2��ֹ����
            integer,dimension(3)                :: faceb1   ! ������������ʼ����
            integer,dimension(3)                :: facee1   ! ������������ֹ����

           ! ������������Χ
            faceb(1) = BC%ib; faceb(2) = BC%jb; faceb(3) = BC%kb;
            facee(1) = BC%ie; facee(2) = BC%je; facee(3) = BC%ke;
            k = sign(1,BC%face-4)*(mod(BC%face-1,3)+1)
            do dim = 1, 3
                if ( dim == abs(k) ) then
                ! ����ָʾά�����⴦��
                faceb(dim) = faceb(dim)-(sign(1,BC%face-4)+1)/2
                facee(dim) = faceb(dim)
                else 
                ! ������ָʾά����������end�ڵ�-1
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

        ! ����ֱ���¶ȱ߽�ķǽǲ�����������
        subroutine upgrade_DirectT_Ghost_Cell(B,BC)
            implicit none 

           ! ����/�������
            type(block_Type) ,pointer           :: B    ! ��ָ��
            type(BC_MSG_Type),pointer           :: BC   ! �߽�ָ��
        end subroutine upgrade_DirectT_Ghost_Cell

        ! ���¿��ڱ߽�ķǽǲ�����������
        subroutine upgrade_Inner_Ghost_Cell(B,BC)
            implicit none

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
            


           ! ��Ԫ�����ͣ���֮ǰ�Ľڵ��Ͳ�һ��

           ! �ڱ߽�
            B_Conne => mesh(BC%nb1)
           ! ������������Χ
            faceb(1) = BC%ib; faceb(2) = BC%jb; faceb(3) = BC%kb;
            facee(1) = BC%ie; facee(2) = BC%je; facee(3) = BC%ke;
            dim_Link = sign(1,BC%face-4)*(mod(BC%face-1,3)+1)
            do dim = 1, 3
                if ( dim == abs(dim_Link) ) then
                    ! ����ά���⴦��
                    faceb(dim) = faceb(dim)+(sign(1,BC%face-4)-1)/2
                    facee(dim) = faceb(dim)
                else 
                    ! ������ά��������end�ڵ�-1
                    facee(dim) = facee(dim)-1
                end if
            end do
            
           ! ������������Χ
            faceb1(1) = BC%ib1; faceb1(2) = BC%jb1; faceb1(3) = BC%kb1; 
            facee1(1) = BC%ie1; facee1(2) = BC%je1; facee1(3) = BC%ke1;
            dim_Link1 = sign(1,BC%face1-4)*(mod(BC%face1-1,3)+1)
            do dim = 1, 3
                if ( dim == abs(dim_Link1) ) then
                    ! ����ά���⴦��
                    faceb1(dim) = faceb1(dim)-(sign(1,BC%face1-4)+1)/2
                    facee1(dim) = faceb1(dim)
                else 
                    ! ������ά��������end�ڵ�-1
                    facee1(dim) = facee1(dim)-1
                end if
            end do

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
                  
        end subroutine upgrade_Inner_Ghost_Cell

        ! ���½ǲ�����������
        subroutine upgrade_Corner_Cell(B)
            implicit none

           ! ����\�������
            type(block_Type),pointer                :: B            ! ��ָ��

           ! �м����
            integer,dimension(3)                    :: s            ! ��Ԫ��ʼ����
            integer,dimension(3)                    :: e            ! ��Ԫ��ֹ����

            s(:) = 1;
            e(1)   = B%nx-1; e(2)   = B%ny-1; e(3)   = B%nz-1; 
           ! ���ÿһ�����12����
            ! x��4�����ֵ
            ! (start,start);(start,end);(end,start);(end,end)
            B%U(:,s(1):e(1),s(2)-1,s(3)-1) = B%U(:,s(1):e(1),s(2)-1,s(3))+B%U(:,s(1):e(1),s(2),s(3)-1)-B%U(:,s(1):e(1),s(2),s(3))
            B%U(:,s(1):e(1),s(2)-1,e(3)+1) = B%U(:,s(1):e(1),s(2)-1,e(3))+B%U(:,s(1):e(1),s(2),e(3)+1)-B%U(:,s(1):e(1),s(2),e(3))
            B%U(:,s(1):e(1),e(2)+1,s(3)-1) = B%U(:,s(1):e(1),e(2)+1,s(3))+B%U(:,s(1):e(1),e(2),s(3)-1)-B%U(:,s(1):e(1),e(2),s(3))
            B%U(:,s(1):e(1),e(2)+1,e(3)+1) = B%U(:,s(1):e(1),e(2)+1,e(3))+B%U(:,s(1):e(1),e(2),e(3)+1)-B%U(:,s(1):e(1),e(2),e(3))
            ! y��4�����ֵ
            B%U(:,s(1)-1,s(2):e(2),s(3)-1) = B%U(:,s(1)-1,s(2):e(2),s(3))+B%U(:,s(1),s(2):e(2),s(3)-1)-B%U(:,s(1),s(2):e(2),s(3))
            B%U(:,s(1)-1,s(2):e(2),e(3)+1) = B%U(:,s(1)-1,s(2):e(2),e(3))+B%U(:,s(1),s(2):e(2),e(3)+1)-B%U(:,s(1),s(2):e(2),e(3))
            B%U(:,e(1)+1,s(2):e(2),s(3)-1) = B%U(:,e(1)+1,s(2):e(2),s(3))+B%U(:,e(1),s(2):e(2),s(3)-1)-B%U(:,e(1),s(2):e(2),s(3))
            B%U(:,e(1)+1,s(2):e(2),e(3)+1) = B%U(:,e(1)+1,s(2):e(2),e(3))+B%U(:,e(1),s(2):e(2),e(3)+1)-B%U(:,e(1),s(2):e(2),e(3))
            ! z��4�����ֵ
            B%U(:,s(1)-1,s(2)-1,s(3):e(3)) = B%U(:,s(1)-1,s(2),s(3):e(3))+B%U(:,s(1),s(2)-1,s(3):e(3))-B%U(:,s(1),s(2),s(3):e(3))
            B%U(:,s(1)-1,e(2)+1,s(3):e(3)) = B%U(:,s(1)-1,e(2),s(3):e(3))+B%U(:,s(1),e(2)+1,s(3):e(3))-B%U(:,s(1),e(2),s(3):e(3))
            B%U(:,e(1)+1,s(2)-1,s(3):e(3)) = B%U(:,e(1)+1,s(2),s(3):e(3))+B%U(:,e(1),s(2)-1,s(3):e(3))-B%U(:,e(1),s(2),s(3):e(3))
            B%U(:,e(1)+1,e(2)+1,s(3):e(3)) = B%U(:,e(1)+1,e(2),s(3):e(3))+B%U(:,e(1),e(2)+1,s(3):e(3))-B%U(:,e(1),e(2),s(3):e(3))
           ! ����ڱ߽紦������Ҫ���⴦��
        end subroutine upgrade_Corner_Cell

       ! �����ڱ߽紦���ǲ��������������
        subroutine upgrade_Inner_Edge_Ghost_Cell(B,BC)
            implicit none 

           ! ����\�������
            type(block_Type),pointer            :: B                ! ��ָ�� 
            type(BC_MSG_Type),pointer           :: BC               ! ������ָ��

           ! �м����
            integer                             :: i,j,k        ! ��Ԫ����
            integer                             :: i1,j1,k1     ! ��������ڵ�����
            integer                             :: dim          ! ά������
            integer                             :: dim_Link     ! ������ά�ȼ�����������
            integer                             :: dim_Link1    ! ��������ά�ȼ�����������
            integer,dimension(3)                :: faceb        ! ��������ʼ����
            integer,dimension(3)                :: facee        ! ��������ֹ����
            integer,dimension(3)                :: facebC       ! ���ǽ����������������ʼ����
            integer,dimension(3)                :: faceeC       ! ���ǽ����������������ֹ����
            integer,dimension(3)                :: faceb1       ! ����������ʼ����
            integer,dimension(3)                :: facee1       ! ����������ֹ����
            integer,dimension(3)                :: owner_O      ! ��������ԭ��
            integer,dimension(3)                :: owner_D      ! ������������
            integer,dimension(3)                :: conne_O      ! ����������ԭ��
            integer,dimension(3)                :: conne_D      ! ��������������
            type(block_Type),pointer            :: B_Conne      ! ���ӿ�ָ��

           ! ��Ԫ�����ͣ���֮ǰ�Ľڵ��Ͳ�һ��
           ! �ڱ߽�
            B_Conne => mesh(BC%nb1)
           ! ������������Χ
            faceb(1) = BC%ib; faceb(2) = BC%jb; faceb(3) = BC%kb;
            facee(1) = BC%ie; facee(2) = BC%je; facee(3) = BC%ke;
            dim_Link = sign(1,BC%face-4)*(mod(BC%face-1,3)+1)
            do dim = 1, 3
                if ( dim == abs(dim_Link) ) then
                    ! ����ά���⴦��
                    faceb(dim) = faceb(dim)+(sign(1,BC%face-4)-1)/2
                    facee(dim) = faceb(dim)
                else 
                    ! ������ά��������end�ڵ�-1
                    facee(dim) = facee(dim)-1
                end if
            end do

           ! ������������Χ
            faceb1(1) = BC%ib1; faceb1(2) = BC%jb1; faceb1(3) = BC%kb1; 
            facee1(1) = BC%ie1; facee1(2) = BC%je1; facee1(3) = BC%ke1;
            dim_Link1 = sign(1,BC%face1-4)*(mod(BC%face1-1,3)+1)
            do dim = 1, 3
                if ( dim == abs(dim_Link1) ) then
                    ! ����ά���⴦��
                    faceb1(dim) = faceb1(dim)-(sign(1,BC%face1-4)+1)/2
                    facee1(dim) = faceb1(dim)
                else 
                    ! ������ά��������end�ڵ�-1
                    facee1(dim) = facee1(dim)-1
                end if
            end do

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

           ! ���ǽ����������������ʼ/��ֹ��������
            facebC(1) = faceb(1)-1; facebC(2) = faceb(2)-1; facebC(3) = faceb(3)-1; 
            faceeC(1) = facee(1)+1; faceeC(2) = facee(2)+1; faceeC(3) = facee(3)+1; 
            facebC(abs(dim_Link)) = facebC(abs(dim_Link))+1;
            faceeC(abs(dim_Link)) = faceeC(abs(dim_Link))-1;

           ! �໥��ֵ
            do k = facebC(3), faceeC(3)
                do j = facebC(2), faceeC(2)
                    do i = facebC(1), faceeC(1)
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
        end subroutine upgrade_Inner_Edge_Ghost_Cell

       ! ���±߽�ֵ
        ! ���µ�һ��ı߽�ֵ
        subroutine upgrade_Boundary(B,Fi,Fj,Fk)
            implicit none

           ! ����\�������
            real(dp),dimension(:,:,:,:)         :: Fi           ! ����I����ͨ���ռ�
            real(dp),dimension(:,:,:,:)         :: Fj           ! ����J����ͨ���ռ�
            real(dp),dimension(:,:,:,:)         :: Fk           ! ����k����ͨ���ռ�
            type(block_Type) ,pointer           :: B            ! ��ָ��

           ! �м����

            integer                             :: ksub         ! �߽�����
            type(BC_MSG_Type),pointer           :: BC       ! ��߽�ָ��

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

        ! ���µ�һ���϶ԳƱ߽��ͨ��
        subroutine upgrade_Symmetry_Boundary(B,BC,Fi,Fj,Fk)
            implicit none

           ! ����\�������
            real(dp),dimension(:,:,:,:)         :: Fi           ! ����I����ͨ���ռ�
            real(dp),dimension(:,:,:,:)         :: Fj           ! ����J����ͨ���ռ�
            real(dp),dimension(:,:,:,:)         :: Fk           ! ����k����ͨ���ռ�
            type(block_Type) ,pointer           :: B            ! ��ָ��
            type(BC_MSG_Type),pointer           :: BC           ! ��߽�ָ��

           ! �м����
            integer                             :: k            ! ��ά��ָʾ
            integer,dimension(3)                :: faceb,facee  ! �߽緶Χ

           ! ����ԳƱ߽�������ֱ�Ӹ�ֵΪ0
            k = sign(1,BC%face-4)*(mod(BC%face-1,3)+1)
            select case (abs(k))
            case (1)    ! i����
                ! ���ұ߽��淶Χ
                faceb(1) = BC%ib  ; faceb(2) = BC%jb  ; faceb(3) = BC%kb; 
                facee(1) = BC%ie-1; facee(2) = BC%je-1; facee(3) = BC%ke-1; 
                facee(abs(k)) = facee(abs(k))+1
                Fi(1,faceb(1):facee(1),faceb(2):facee(2),faceb(3):facee(3)) = 0 
            case (2)    ! j����
                ! ���ұ߽��淶Χ
                faceb(1) = BC%ib  ; faceb(2) = BC%jb  ; faceb(3) = BC%kb; 
                facee(1) = BC%ie-1; facee(2) = BC%je-1; facee(3) = BC%ke-1; 
                facee(abs(k)) = facee(abs(k))+1
                Fj(1,faceb(1):facee(1),faceb(2):facee(2),faceb(3):facee(3)) = 0
            case (3)    ! k����
                ! ���ұ߽��淶Χ
                faceb(1) = BC%ib  ; faceb(2) = BC%jb  ; faceb(3) = BC%kb; 
                facee(1) = BC%ie-1; facee(2) = BC%je-1; facee(3) = BC%ke-1; 
                facee(abs(k)) = facee(abs(k))+1
                Fk(1,faceb(1):facee(1),faceb(2):facee(2),faceb(3):facee(3)) = 0 
            end select
           

        end subroutine upgrade_Symmetry_Boundary

        ! ���µ�һ����ֱ�������߽��ͨ��
        subroutine upgrade_DirectQ_Boundary(B,BC,Fi,Fj,Fk)
            implicit none

            ! ����\�������
             real(dp),dimension(:,:,:,:)         :: Fi           ! ����I����ͨ���ռ�
             real(dp),dimension(:,:,:,:)         :: Fj           ! ����J����ͨ���ռ�
             real(dp),dimension(:,:,:,:)         :: Fk           ! ����k����ͨ���ռ�
             type(block_Type) ,pointer           :: B            ! ��ָ��
             type(BC_MSG_Type),pointer           :: BC           ! ��߽�ָ��

           ! �м����
            integer                             :: k            ! ��ά��ָʾ
            integer,dimension(3)                :: faceb,facee  ! �߽緶Χ

           ! �߽���������Ϊ�����֣�1���ȱ߽�������2���ǵȱ߽�����
            if (wallq_Iso) then
               ! �ȱ߽�����
                k = sign(1,BC%face-4)*(mod(BC%face-1,3)+1)
                select case (abs(k))
                case (1)    ! i����
                    ! ���ұ߽��淶Χ
                    faceb(1) = BC%ib  ; faceb(2) = BC%jb  ; faceb(3) = BC%kb; 
                    facee(1) = BC%ie-1; facee(2) = BC%je-1; facee(3) = BC%ke-1; 
                    facee(abs(k)) = facee(abs(k))+1
                    Fi(1,faceb(1):facee(1),faceb(2):facee(2),faceb(3):facee(3)) = wallq_Iso_q 
                case (2)    ! j����
                    ! ���ұ߽��淶Χ
                    faceb(1) = BC%ib  ; faceb(2) = BC%jb  ; faceb(3) = BC%kb; 
                    facee(1) = BC%ie-1; facee(2) = BC%je-1; facee(3) = BC%ke-1; 
                    facee(abs(k)) = facee(abs(k))+1
                    Fj(1,faceb(1):facee(1),faceb(2):facee(2),faceb(3):facee(3)) = wallq_Iso_q 
                case (3)    ! k����
                    ! ���ұ߽��淶Χ
                    faceb(1) = BC%ib  ; faceb(2) = BC%jb  ; faceb(3) = BC%kb; 
                    facee(1) = BC%ie-1; facee(2) = BC%je-1; facee(3) = BC%ke-1; 
                    facee(abs(k)) = facee(abs(k))+1
                    Fk(1,faceb(1):facee(1),faceb(2):facee(2),faceb(3):facee(3)) = wallq_Iso_q  
                end select
            else 
               ! �ǵȱ߽�������δ��ɣ�
            end if
 
        end subroutine upgrade_DirectQ_Boundary

        ! ���µ�һ���϶����߽��ͨ��
        subroutine upgrade_Convection_Boundary(B,BC,Fi,Fj,Fk)
            implicit none

            ! ����\�������
             real(dp),dimension(:,:,:,:)         :: Fi                  ! ����I����ͨ���ռ�
             real(dp),dimension(:,:,:,:)         :: Fj                  ! ����J����ͨ���ռ�
             real(dp),dimension(:,:,:,:)         :: Fk                  ! ����k����ͨ���ռ�
             type(block_Type) ,pointer           :: B                   ! ��ָ��
             type(BC_MSG_Type),pointer           :: BC                  ! ��߽�ָ��

             ! �м����
            integer                             :: dim_Sign             ! ��ά��ָʾ
            integer                             :: i,j,k                ! ��ά������
            integer                             :: cellI,cellJ,cellK    ! ��Ԫά������
            integer,dimension(3)                :: faceb,facee          ! �߽緶Χ

           ! �������ȱ߽��ϵĶ������ȷ�Ϊ�����֣�1���ȶ������ȣ�2���ǵȶ�������
            if (wallq_Iso) then
               ! �ȶ�������
                dim_Sign = sign(1,BC%face-4)*(mod(BC%face-1,3)+1)
                select case (abs(dim_Sign))
                case (1)    ! i����
                    ! ���ұ߽��淶Χ
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
                case (2)    ! j����
                    ! ���ұ߽��淶Χ
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
                case (3)    ! k����
                    ! ���ұ߽��淶Χ
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
               ! �ǵȶ������ȣ�δ��ɣ�
            end if
             
 
        end subroutine upgrade_Convection_Boundary

        ! ���µ�һ���϶����߽��ͨ��(δ���)
        subroutine upgrade_Inner_Boundary(B,BC,Fi,Fj,Fk)
            implicit none

            ! ����\�������
             real(dp),dimension(:,:,:,:)         :: Fi           ! ����I����ͨ���ռ�
             real(dp),dimension(:,:,:,:)         :: Fj           ! ����J����ͨ���ռ�
             real(dp),dimension(:,:,:,:)         :: Fk           ! ����k����ͨ���ռ�
             type(block_Type) ,pointer           :: B            ! ��ָ��
             type(BC_MSG_Type),pointer           :: BC           ! ��߽�ָ��
 
        end subroutine upgrade_Inner_Boundary


    end module Boundary
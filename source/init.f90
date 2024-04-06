!--------------------------------------------------------------------------------------------------
! ��ʼ��ģ���ļ�
!--------------------------------------------------------------------------------------------------
! ��ʼ���߽��������ڲ��¶ȳ�ģ��(δ���)
    module init
        use const_Var
        use global_Var
        implicit none
        
    contains
        ! ���ó�ʼ�ڲ��¶ȳ���δ��ɣ�
        subroutine init_Field()
            if (Iflag_init) then 
                if (initT_Iso) then
                    call initT_Iso_Field        ! ��ֵ��ʼ���¶ȳ�
                else 
                    call initT_None_Iso_Field   ! ��δ��ɣ�
                end if
            else 
                call initT_Field_File           ! ��δ��ɣ�
            end if
        end subroutine init_Field

        ! ���³�ʼ���¶ȳ�����
        subroutine initT_Iso_Field()
            implicit none

           ! �м����
            integer                         :: m        ! ������
            type(block_Type),pointer        :: B        ! ��ָ��

           ! ���³�ʼ���¶ȳ�
            do m = 1, num_Block
                B => mesh(m)
                B%U(:,1:B%nx-1,1:B%ny-1,1:B%nz-1) = initT_Iso_T
            end do
        end subroutine initT_Iso_Field
        
        ! �����ļ��ǵ��³�ʼ���¶ȳ�������δ��ɣ�
        subroutine initT_None_Iso_Field()
        end subroutine initT_None_Iso_Field

        ! ���������ļ���ʼ��������δ��ɣ�
        subroutine initT_Field_File()
        end subroutine initT_Field_File

        ! ��ʼ���߽�����
        subroutine init_Boundary()
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

        ! ��ʼ���ԳƱ߽�����
        subroutine init_Symmetry_Boundary(B,BC)
            implicit none 

           ! ����/�������
            type(block_Type) ,pointer           :: B        ! ��ָ��
            type(BC_MSG_Type),pointer           :: BC       ! �߽�ָ��

           ! �м����
            integer                             :: k        ! ά��ָʾ
            integer,dimension(3)                :: db,de    ! ���淶Χ

            db(1) = BC%ib; db(2) = BC%jb; db(3) = BC%kb;
            de(1) = BC%ie; de(2) = BC%je; de(3) = BC%ke;
            k = sign(1,BC%face-4)*(mod(BC%face-1,3)+1)      ! k�ķ��Ŵ������������k�����ִ�������ά��

            if ( k>0 ) then
            else 
            end if

        end subroutine init_Symmetry_Boundary

        ! ��ʼ��ֱ�������߽�����
        subroutine init_DirectQ_Boundary(B,BC)
            implicit none 

            ! ����/�������
             type(block_Type) ,pointer           :: B    ! ��ָ��
             type(BC_MSG_Type),pointer           :: BC   ! �߽�ָ��
        end subroutine init_DirectQ_Boundary

        ! ��ʼ���������ȱ߽�����
        subroutine init_Convection_Boundary(B,BC)
            implicit none 

            ! ����/�������
             type(block_Type) ,pointer           :: B    ! ��ָ��
             type(BC_MSG_Type),pointer           :: BC   ! �߽�ָ��
        end subroutine init_Convection_Boundary

        ! ��ʼ��ֱ���¶ȱ߽�����
        subroutine init_DirectT_Boundary(B,BC)
            implicit none 

            ! ����/�������
             type(block_Type) ,pointer           :: B    ! ��ָ��
             type(BC_MSG_Type),pointer           :: BC   ! �߽�ָ��
        end subroutine init_DirectT_Boundary

        ! ��ʼ���ڲ��߽�����
        subroutine init_Inner_Boundary(B,BC)
            implicit none 

            ! ����/�������
             type(block_Type) ,pointer           :: B    ! ��ָ��
             type(BC_MSG_Type),pointer           :: BC   ! �߽�ָ��
        end subroutine init_Inner_Boundary
    end module init
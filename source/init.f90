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
                ! if (m ==2 ) then 
                !     B%U(:,1:B%nx-1,1:B%ny-1,1:B%nz-1) = 500._dp
                ! end if
            end do
        end subroutine initT_Iso_Field
        
        ! �����ļ��ǵ��³�ʼ���¶ȳ�������δ��ɣ�
        subroutine initT_None_Iso_Field()
        end subroutine initT_None_Iso_Field

        ! ���������ļ���ʼ��������δ��ɣ�
        subroutine initT_Field_File()
        end subroutine initT_Field_File

    end module init
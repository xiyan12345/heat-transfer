!--------------------------------------------------------------------------------------------------
! debugģ���ļ�
!--------------------------------------------------------------------------------------------------
! debugģ�飺��ģ��������debug֮��
    module debug
        use precision_EC
        implicit none
        
    contains
        subroutine show_Flux(F)
            implicit none

           ! ����\�������
            real(dp),dimension(:,:,:,:)             :: F
        end subroutine show_Flux
        
    end module debug
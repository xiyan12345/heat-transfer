!--------------------------------------------------------------------------------------------------
! CFL����ģ���ļ�
!--------------------------------------------------------------------------------------------------
! CFLģ��
    module CFL
        use const_Var
        use global_Var
        implicit none
        
    contains
        subroutine calculate_dt()
            implicit none

            real(dp)                :: min_delt         ! ����߶���Сֵ

           ! ȷ�� min_delt(������)
            min_delt = 0
           ! ����dt
            dt = CFL_Value*min_delt/1
        end subroutine
        
    end module CFL
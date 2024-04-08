!--------------------------------------------------------------------------------------------------
! ͨ��ģ���ļ�
!--------------------------------------------------------------------------------------------------
! ͨ��ģ��
    module Flux
        use block_Type_Def
        use read_Mesh
        implicit none
        
    contains
       ! ����i�����ճ��ͨ����δ��ɣ�:�����⣬ԭ������ʸ��������
        function viscous_Flux_I(B,i,j,k)
            implicit none

           ! ����/�������
            integer                                 :: i,j,k            ! ����I���������
            type(block_Type),pointer                :: B                ! ��ָ��
            real(dp)                                :: viscous_Flux_i   ! �������  

           ! �м��
            real(dp)                 :: vol         ! �μ���Ԫ���
            real(dp)                 :: T1,T2,T3    ! i-��,j-��,k-���ϵ�������
            real(dp)                 :: T4,T5,T6    ! i+��,j+��,k+���ϵ�������
            real(dp),dimension(3)    :: S1,S2,S3    ! i-��,j-��,k-���ϵķ�����
            real(dp),dimension(3)    :: S4,S5,S6    ! i+��,j+��,k+���ϵķ�����
            real(dp),dimension(3)    :: p1,p2,p3,p4 ! �μ���Ԫ�ĵ��λʸ��
            real(dp),dimension(3)    :: p5,p6,p7,p8 ! �μ���Ԫ�ĵ��λʸ��
            real(dp),dimension(3)    :: delt_T      ! ����I�����ϵ��¶��ݶ�
            real(dp),dimension(3)    :: S           ! ����I�����ϵķ���ʸ��

           ! ����ʸ��ָ��
            S(1) = B%ni1(i,j,k); S(2) = B%ni2(i,j,k); S(3) = B%ni3(i,j,k);
           ! ��ָ��:
            ! I���棨I,J,K)->(i,j,k),(i,j+1,k),(i,j,k+1),(i,j+1,k+1)
            ! p1->(i-1/2,j,k)  ,p2->(i+1/2,j,k)    ,p3->(i-1/2,j+1,k)  ,p4->(i+1/2,j+1,k)
            ! p5->(i-1/2,j,k+1),p6->(i+1/2,j,k+1)  ,p7->(i-1/2,j+1,k+1),p8->(i+1/2,j+1,k+1)
            p1(1) = 0.5_dp*(B%x(i-1,j,k)+B%x(i,j,k))    ;p1(2) = 0.5_dp*(B%y(i-1,j,k)+B%y(i,j,k))    ;p1(3) = 0.5_dp*(B%z(i-1,j,k)+B%z(i,j,k))
            p2(1) = 0.5_dp*(B%x(i+1,j,k)+B%x(i,j,k))    ;p2(2) = 0.5_dp*(B%y(i+1,j,k)+B%y(i,j,k))    ;p2(3) = 0.5_dp*(B%z(i+1,j,k)+B%z(i,j,k))
            p3(1) = 0.5_dp*(B%x(i-1,j+1,k)+B%x(i,j+1,k));p3(2) = 0.5_dp*(B%y(i-1,j+1,k)+B%y(i,j+1,k));p3(3) = 0.5_dp*(B%z(i-1,j+1,k)+B%z(i,j+1,k))
            p4(1) = 0.5_dp*(B%x(i+1,j+1,k)+B%x(i,j+1,k));p4(2) = 0.5_dp*(B%y(i+1,j+1,k)+B%y(i,j+1,k));p4(3) = 0.5_dp*(B%z(i+1,j+1,k)+B%z(i,j+1,k))

            p5(1) = 0.5_dp*(B%x(i-1,j,k+1)+B%x(i,j,k+1))    ;p5(2) = 0.5_dp*(B%y(i-1,j,k+1)+B%y(i,j,k+1))    ;p5(3) = 0.5_dp*(B%z(i-1,j,k+1)+B%z(i,j,k+1))
            p6(1) = 0.5_dp*(B%x(i+1,j,k+1)+B%x(i,j,k+1))    ;p6(2) = 0.5_dp*(B%y(i+1,j,k+1)+B%y(i,j,k+1))    ;p6(3) = 0.5_dp*(B%z(i+1,j,k+1)+B%z(i,j,k+1))
            p7(1) = 0.5_dp*(B%x(i-1,j+1,k+1)+B%x(i,j+1,k+1));p7(2) = 0.5_dp*(B%y(i-1,j+1,k+1)+B%y(i,j+1,k+1));p7(3) = 0.5_dp*(B%z(i-1,j+1,k+1)+B%z(i,j+1,k+1))
            p8(1) = 0.5_dp*(B%x(i+1,j+1,k+1)+B%x(i,j+1,k+1));p8(2) = 0.5_dp*(B%y(i+1,j+1,k+1)+B%y(i,j+1,k+1));p8(3) = 0.5_dp*(B%z(i+1,j+1,k+1)+B%z(i,j+1,k+1))

           ! ���㵥Ԫ���ⷨ��ʸ�������
            ! ����i-��
            S1 = -1*compute_S(p1,p3,p5,p7)
            ! ����j-��
            S2 = -1*compute_S(p1,p5,p2,p6)
            ! ����k-��
            S3 = -1*compute_S(p1,p2,p3,p4)
            ! ����i+��
            S4 = compute_S(p2,p4,p6,p8)
            ! ����j+��
            S5 = compute_S(p3,p7,p4,p8)
            ! ����k+��
            S6 = compute_S(p5,p6,p7,p8)
            
            ! �������
            vol = compute_Volume2(p1,p2,p3,p4,p5,p6,p7,p8,-1*S1,-1*S2,-1*S3,S4,S5,S6)
           ! ��������������
            ! ����i-��
            T1 = B%U(1,i-1,j,k)
            ! ����j-��
            T2 = 0.25_dp*(B%U(1,i-1,j-1,k)+B%U(1,i,j-1,k)+B%U(1,i-1,j,k)+B%U(1,i,j,k))
            ! ����k-��
            T3 = 0.25_dp*(B%U(1,i-1,j,k-1)+B%U(1,i,j,k-1)+B%U(1,i-1,j,k)+B%U(1,i,j,k))
            ! ����i+��
            T4 = B%U(1,i,j,k)
            ! ����j+��
            T5 = 0.25_dp*(B%U(1,i-1,j,k)+B%U(1,i,j,k)+B%U(1,i-1,j+1,k)+B%U(1,i,j+1,k))
            ! ����k+��
            T6 = 0.25_dp*(B%U(1,i-1,j,k)+B%U(1,i,j,k)+B%U(1,i-1,j,k+1)+B%U(1,i,j,k+1))

           ! ��������I�����ϵ��¶��ݶ�
            delt_T(:) = (T1*S1+T2*S2+T3*S3+T4*S4+T5*S5+T6*S6)/vol
           ! ��������I�����ϵ�ճ��ͨ��
            viscous_Flux_i = dot_product(delt_T,S)




        end function viscous_Flux_I

       ! ����i�����ճ��ͨ����δ��ɣ�
        function viscous_Flux_J(B,i,j,k)
            implicit none

           ! ����/�������
            integer                                 :: i,j,k            ! ����J���������
            type(block_Type),pointer                :: B                ! ��ָ��
            real(dp)                                :: viscous_Flux_j   ! �������  
        end function viscous_Flux_J
       ! ����i�����ճ��ͨ����δ��ɣ�
        function viscous_Flux_K(B,i,j,k)
            implicit none

           ! ����/�������
            integer                                 :: i,j,k            ! ����J���������
            type(block_Type),pointer                :: B                ! ��ָ��
            real(dp)                                :: viscous_Flux_K   ! �������  
        end function viscous_Flux_K
        
    end module Flux
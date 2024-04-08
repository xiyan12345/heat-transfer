!--------------------------------------------------------------------------------------------------
! 通量模块文件
!--------------------------------------------------------------------------------------------------
! 通量模块
    module Flux
        use block_Type_Def
        use read_Mesh
        implicit none
        
    contains
       ! 计算i向面的粘性通量（未完成）:有问题，原因是面矢量的问题
        function viscous_Flux_I(B,i,j,k)
            implicit none

           ! 输入/输出变量
            integer                                 :: i,j,k            ! 所求I向面的索引
            type(block_Type),pointer                :: B                ! 块指针
            real(dp)                                :: viscous_Flux_i   ! 输出变量  

           ! 中间变
            real(dp)                 :: vol         ! 次级单元体积
            real(dp)                 :: T1,T2,T3    ! i-面,j-面,k-面上的物理量
            real(dp)                 :: T4,T5,T6    ! i+面,j+面,k+面上的物理量
            real(dp),dimension(3)    :: S1,S2,S3    ! i-面,j-面,k-面上的法向量
            real(dp),dimension(3)    :: S4,S5,S6    ! i+面,j+面,k+面上的法向量
            real(dp),dimension(3)    :: p1,p2,p3,p4 ! 次级单元的点的位矢量
            real(dp),dimension(3)    :: p5,p6,p7,p8 ! 次级单元的点的位矢量
            real(dp),dimension(3)    :: delt_T      ! 所求I向面上的温度梯度
            real(dp),dimension(3)    :: S           ! 所求I向面上的法向矢量

           ! 法向矢量指定
            S(1) = B%ni1(i,j,k); S(2) = B%ni2(i,j,k); S(3) = B%ni3(i,j,k);
           ! 点指定:
            ! I向面（I,J,K)->(i,j,k),(i,j+1,k),(i,j,k+1),(i,j+1,k+1)
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

           ! 计算单元面外法向矢量及体积
            ! 计算i-面
            S1 = -1*compute_S(p1,p3,p5,p7)
            ! 计算j-面
            S2 = -1*compute_S(p1,p5,p2,p6)
            ! 计算k-面
            S3 = -1*compute_S(p1,p2,p3,p4)
            ! 计算i+面
            S4 = compute_S(p2,p4,p6,p8)
            ! 计算j+面
            S5 = compute_S(p3,p7,p4,p8)
            ! 计算k+面
            S6 = compute_S(p5,p6,p7,p8)
            
            ! 计算体积
            vol = compute_Volume2(p1,p2,p3,p4,p5,p6,p7,p8,-1*S1,-1*S2,-1*S3,S4,S5,S6)
           ! 计算面上物理量
            ! 计算i-面
            T1 = B%U(1,i-1,j,k)
            ! 计算j-面
            T2 = 0.25_dp*(B%U(1,i-1,j-1,k)+B%U(1,i,j-1,k)+B%U(1,i-1,j,k)+B%U(1,i,j,k))
            ! 计算k-面
            T3 = 0.25_dp*(B%U(1,i-1,j,k-1)+B%U(1,i,j,k-1)+B%U(1,i-1,j,k)+B%U(1,i,j,k))
            ! 计算i+面
            T4 = B%U(1,i,j,k)
            ! 计算j+面
            T5 = 0.25_dp*(B%U(1,i-1,j,k)+B%U(1,i,j,k)+B%U(1,i-1,j+1,k)+B%U(1,i,j+1,k))
            ! 计算k+面
            T6 = 0.25_dp*(B%U(1,i-1,j,k)+B%U(1,i,j,k)+B%U(1,i-1,j,k+1)+B%U(1,i,j,k+1))

           ! 计算所求I向面上的温度梯度
            delt_T(:) = (T1*S1+T2*S2+T3*S3+T4*S4+T5*S5+T6*S6)/vol
           ! 计算所求I向面上的粘性通量
            viscous_Flux_i = dot_product(delt_T,S)




        end function viscous_Flux_I

       ! 计算i向面的粘性通量（未完成）
        function viscous_Flux_J(B,i,j,k)
            implicit none

           ! 输入/输出变量
            integer                                 :: i,j,k            ! 所求J向面的索引
            type(block_Type),pointer                :: B                ! 块指针
            real(dp)                                :: viscous_Flux_j   ! 输出变量  
        end function viscous_Flux_J
       ! 计算i向面的粘性通量（未完成）
        function viscous_Flux_K(B,i,j,k)
            implicit none

           ! 输入/输出变量
            integer                                 :: i,j,k            ! 所求J向面的索引
            type(block_Type),pointer                :: B                ! 块指针
            real(dp)                                :: viscous_Flux_K   ! 输出变量  
        end function viscous_Flux_K
        
    end module Flux
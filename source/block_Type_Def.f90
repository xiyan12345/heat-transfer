!--------------------------------------------------------------------------------------------------
! 类型定义文件
!--------------------------------------------------------------------------------------------------
! 块类型文件
    module block_Type_Def
        use precision_EC
        implicit none

        type BC_MSG_Type                                                ! 块边界类型
            integer :: ib,ie,jb,je,kb,ke,bc,face                        ! 主块边界连接信息
            integer :: ib1,ie1,jb1,je1,kb1,ke1,nb1,face1                ! 连接块边界连接信息
            integer :: L(3)                                             ! 连接顺序符
        end type BC_MSG_Type

        type block_Type                 ! 网格块类型
            integer     :: nx,ny,nz     ! 网格节点范围
            integer     :: subface      ! 网格子面数量

            ! 节点型参量
            real(dp),allocatable,dimension(:,:,:)   :: x,y,z            ! 节点坐标x,y,z
            ! 面型参量
            real(dp),allocatable,dimension(:,:,:)   :: Si,ni1,ni2,ni3   ! i向面的面积和沿x,y,z坐标的法向量
            real(dp),allocatable,dimension(:,:,:)   :: Sj,nj1,nj2,nj3   ! j向面的面积和沿x,y,z坐标的法向量
            real(dp),allocatable,dimension(:,:,:)   :: Sk,nk1,nk2,nk3   ! i向面的面积和沿x,y,z坐标的法向量
            ! 单元型变量
            real(dp),allocatable,dimension(:,:,:)   :: vol              ! 单元体积
            real(dp),allocatable,dimension(:,:,:,:) :: U                ! 单元物理参量矢量
            real(dp),allocatable,dimension(:,:,:,:) :: Un               ! 单元暂态物理参量矢量
            real(dp),allocatable,dimension(:,:,:,:) :: res              ! 单元物理参量残差量
            ! 边界型变量
            type(BC_MSG_Type),allocatable,dimension(:)  :: BC_MSG       ! 块边界
 
        end type block_Type
        
    end module block_Type_Def
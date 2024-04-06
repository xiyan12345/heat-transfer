!--------------------------------------------------------------------------------------------------
! 残差模块文件
!--------------------------------------------------------------------------------------------------
! 残差模块
    module residual
        use global_Var
        use Boundary
        implicit none
        
    contains
        subroutine compute_Residual()
            implicit none

           ! 中间变量
            integer                                     :: m            ! 块索引变量
            integer                                     :: I,J,K        ! 块单元I,J,K向索引
            real(dp),allocatable,dimension(:,:,:,:)     :: Fi           ! 单块I向面通量空间
            real(dp),allocatable,dimension(:,:,:,:)     :: Fj           ! 单块J向面通量空间
            real(dp),allocatable,dimension(:,:,:,:)     :: Fk           ! 单块k向面通量空间
            type(block_Type),pointer                    :: B            ! 块指针

            do m = 1, num_Block
                B => mesh(m)
               ! 对于每一块临时申请面通量空间
                allocate(Fi(nVar,B%nx,B%ny-1,B%nz-1))                       
                allocate(Fj(nVar,B%nx-1,B%ny,B%nz-1)) 
                allocate(Fk(nVar,B%nx-1,B%ny-1,B%nz)) 

               ! 针对每一块计算边界面通量(未完成)
                call upgrade_Boundary(B,Fi,Fj,Fk)   ! 这个是编译器有问题

               ! 计算每一块内部面通量(未完成)
               ! 计算每一块所有单元的残差量(未完成)
                do K = 1, B%nz-1
                    do J = 1, B%ny-1
                        do I = 1, B%nx-1
                            B%res(:,I,J,K) = Fi(:,I,J,K)-Fi(:,I+1,J,K)+Fj(:,I,J,K)-Fj(:,I,J+1,K)
                        end do
                    end do
                end do

               ! 释放之前申请的临时面通量空间
                deallocate(Fi)
                deallocate(Fj)
                deallocate(Fk)
            end do
        end subroutine compute_Residual
        
    end module residual
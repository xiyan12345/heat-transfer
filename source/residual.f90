!--------------------------------------------------------------------------------------------------
! �в�ģ���ļ�
!--------------------------------------------------------------------------------------------------
! �в�ģ��
    module residual
        use global_Var
        use Boundary
        use Flux
        implicit none
        
    contains
        subroutine compute_Residual()
            implicit none

           ! �м����
            integer                                     :: m            ! ����������
            integer                                     :: I,J,K        ! �鵥ԪI,J,K������
            real(dp),allocatable,dimension(:,:,:,:)     :: Fi           ! ����I����ͨ���ռ�
            real(dp),allocatable,dimension(:,:,:,:)     :: Fj           ! ����J����ͨ���ռ�
            real(dp),allocatable,dimension(:,:,:,:)     :: Fk           ! ����k����ͨ���ռ�
            type(block_Type),pointer                    :: B            ! ��ָ��

            do m = 1, num_Block
                B => mesh(m)
               ! ����ÿһ����ʱ������ͨ���ռ�
                allocate(Fi(nVar,B%nx,B%ny-1,B%nz-1))                       
                allocate(Fj(nVar,B%nx-1,B%ny,B%nz-1)) 
                allocate(Fk(nVar,B%nx-1,B%ny-1,B%nz)) 

               ! ���ÿһ�����߽���ͨ��(δ���)
                call upgrade_Boundary(B,Fi,Fj,Fk)   

               ! ����ÿһ���ڲ���ͨ��(δ���)
                ! ����I�����ͨ��
                do k = 1, B%nz-1
                    do j = 1, B%ny-1
                        do i = 2, B%nx-1
                            Fi(1,i,j,k) = viscous_Flux_I(B,i,j,k)
                        end do
                    end do
                end do
                ! ����J�����ͨ��
                do k = 1, B%nz-1
                    do j = 2, B%ny-1
                        do i = 1, B%nx-1
                            FJ(1,i,j,k) = viscous_Flux_J(B,i,j,k)
                        end do
                    end do
                end do
                ! ����K�����ͨ��
                do k = 2, B%nz-1
                    do j = 1, B%ny-1
                        do i = 1, B%nx-1
                            FK(1,i,j,k) = viscous_Flux_K(B,i,j,k)
                        end do
                    end do
                end do
                
               ! ����ÿһ�����е�Ԫ�Ĳв���(δ���)
                do K = 1, B%nz-1
                    do J = 1, B%ny-1
                        do I = 1, B%nx-1
                            B%res(:,I,J,K) = Fi(:,I,J,K)-Fi(:,I+1,J,K)+Fj(:,I,J,K)-Fj(:,I,J+1,K)
                        end do
                    end do
                end do

               ! �ͷ�֮ǰ�������ʱ��ͨ���ռ�
                deallocate(Fi)
                deallocate(Fj)
                deallocate(Fk)
            end do
        end subroutine compute_Residual
        
    end module residual
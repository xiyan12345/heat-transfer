!--------------------------------------------------------------------------------------------------
! ���Ͷ����ļ�
!--------------------------------------------------------------------------------------------------
! �������ļ�
    module block_Type_Def
        use precision_EC
        implicit none

        type BC_MSG_Type                                                ! ��߽�����
            integer :: ib,ie,jb,je,kb,ke,bc,face                        ! ����߽�������Ϣ
            integer :: ib1,ie1,jb1,je1,kb1,ke1,nb1,face1                ! ���ӿ�߽�������Ϣ
            integer :: L(3)                                             ! ����˳���
        end type BC_MSG_Type

        type block_Type                 ! ���������
            integer     :: nx,ny,nz     ! ����ڵ㷶Χ
            integer     :: subface      ! ������������

            ! �ڵ��Ͳ���
            real(dp),allocatable,dimension(:,:,:)   :: x,y,z            ! �ڵ�����x,y,z
            ! ���Ͳ���
            real(dp),allocatable,dimension(:,:,:)   :: Si,ni1,ni2,ni3   ! i������������x,y,z����ķ�����
            real(dp),allocatable,dimension(:,:,:)   :: Sj,nj1,nj2,nj3   ! j������������x,y,z����ķ�����
            real(dp),allocatable,dimension(:,:,:)   :: Sk,nk1,nk2,nk3   ! i������������x,y,z����ķ�����
            ! ��Ԫ�ͱ���
            real(dp),allocatable,dimension(:,:,:)   :: vol              ! ��Ԫ���
            real(dp),allocatable,dimension(:,:,:,:) :: U                ! ��Ԫ�������ʸ��
            real(dp),allocatable,dimension(:,:,:,:) :: Un               ! ��Ԫ��̬�������ʸ��
            real(dp),allocatable,dimension(:,:,:,:) :: res              ! ��Ԫ��������в���
            ! �߽��ͱ���
            type(BC_MSG_Type),allocatable,dimension(:)  :: BC_MSG       ! ��߽�
 
        end type block_Type
        
    end module block_Type_Def
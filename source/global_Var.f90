!--------------------------------------------------------------------------------------------------
! ȫ�ֱ���ģ���ļ�
!--------------------------------------------------------------------------------------------------
! ȫ�ֱ���ģ��
    ! ȫ�ֱ���ģ��
    module global_Var
        use const_Var
        use precision_EC
        use block_Type_Def
        implicit none

        ! �������
        integer                                     :: iDim                     ! �������ά��
        integer                                     :: total_Step               ! �ܵ�������
        integer                                     :: kStep_Save               ! ������������
        integer                                     :: nVar                     ! ������Ŀ
        real(dp)                                    :: gscale                   ! ����λ
        real(dp)                                    :: rho                      ! �����ܶ�
        real(dp)                                    :: cp                       ! ��������ϵ��
        real(dp)                                    :: kh                       ! �ȴ���ϵ�� 
        real(dp)                                    :: wallq_Iso_q              ! �ȱ�������ֵ
        real(dp)                                    :: wallc_Iso_CH             ! �ȱ����������ϵ��
        real(dp)                                    :: wallc_Iso_T              ! �ȱ�����Ȼ����¶�
        real(dp)                                    :: initT_Iso_T              ! �ȳ�ʼ�¶ȳ�ֵ
        real(dp)                                    :: dt                       ! �ƽ�ʱ�䲽
        real(dp)                                    :: real_Time                ! ��ʵʱ��
        real(dp),dimension(RKn)                     :: Ralpha                   ! �������ϵ��1
        real(dp),dimension(RKn)                     :: Rgamma                   ! �������ϵ��2
        real(dp),dimension(RKn)                     :: Rbeta                    ! �������ϵ��3
        character(len=str_Len )                     :: project_Name             ! ��Ŀ����
        character(len=file_Len)                     :: mesh_File                ! ���������ļ�
        character(len=file_Len)                     :: inp_File                 ! ���ӹ�ϵ�ļ�
        character(len=file_Len)                     :: save_Dir                 ! �����ļ���
        character(len=file_Len)                     :: wallq_qfile              ! �ǵȱ������������ļ�
        character(len=file_Len)                     :: wallc_CHfile             ! �ǵȱ����������ϵ�������ļ�
        character(len=file_Len)                     :: wallc_Tfile              ! �ǵȱ�����������¶������ļ�
        character(len=file_Len)                     :: initT_Tfile              ! �ǵȳ�ʼ�¶ȳ������ļ�
        character(len=file_Len)                     :: continued_File           ! �����ʼ���ڲ��¶ȳ������ļ�
        logical                                     :: Iflag_init               ! ������Ʒ�
        logical                                     :: wallq_Iso                ! �ȱ����������Ʒ�
        logical                                     :: wallc_Iso                ! �ȶ������ȿ��Ʒ�
        logical                                     :: initT_Iso                ! �ȳ�ʼ�¶ȳ����Ʒ�

        
        namelist /prj/ project_Name,iDim,gscale
        namelist /step/ total_Step,kStep_Save,Iflag_init,continued_File
        namelist /files/ mesh_File,inp_File,save_Dir
        namelist /material/ rho,cp,kh
        namelist /wallq/ wallq_Iso,wallq_Iso_q,wallq_qfile
        namelist /wallc/ wallc_Iso,wallc_Iso_CH,wallc_Iso_T,wallc_CHfile,wallc_Tfile
        namelist /initT/ initT_Iso,initT_Iso_T,initT_Tfile

        ! �������
        type(block_Type),pointer    ,dimension(:)   :: mesh                     ! ������ֵ��
        real(dp)        ,allocatable,dimension(:)   :: res_rms                  ! �������в�

        ! ��̬����
        integer                                     :: num_Block                ! �������


        
        
    end module global_Var

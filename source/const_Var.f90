!--------------------------------------------------------------------------------------------------
! ����ģ���ļ�
!--------------------------------------------------------------------------------------------------
! ����ģ��
    module const_Var
        use precision_EC
        implicit none
        ! 2D����ԭ��
            ! 2D������뻭��xyƽ���ϣ������ά����ʱ����2Dת��Ϊ3D;
            ! 2D -> 3D��z��K�����죬��K=3
        ! �Զ���߽����ͣ�
            ! �����ȴ���������ı߽������������������ı߽��������Ͳ�һ�£�����Generic��ʽ���ṩ�������Զ���߽����ͣ�
            ! ��˽����ǽ����Զ���,����������
            ! generic#1(8):ֱ�������߽�����
            ! generic#2(9):�������ȱ߽�����
            ! generuc#3(10):ֱ�ӱ��±߽���������ʱû������
        ! unit˵����
            ! 1��6����ʾ��
            ! 2��22��control.dat
            ! 3��23��log.dat
            ! 4��24��mesh.dat
            ! 5��25��inp.dat
            ! 6��26��result.dat

        integer,parameter       :: RKn              = 3     ! ����-��������
        integer,parameter       :: LAP              = 2     ! ���������
        integer,parameter       :: zStretch         = 3     ! z���������         
        integer ,parameter      :: BC_Direct_Q      = 8     ! ֱ�������߽�����
        integer,parameter       :: BC_Convection    = 9     ! �������ȱ߽�����
        integer,parameter       :: BC_Direct_T      = 10    ! ֱ���¶ȱ߽�����              
        integer,parameter       :: BC_Symmetry      = 3     ! �ԳƱ߽�����
        integer,parameter       :: un_control       = 22    ! control.dat
        integer,parameter       :: un_log           = 23    ! log.dat
        integer,parameter       :: un_mesh          = 24    ! mesh.dat
        integer,parameter       :: un_inp           = 25    ! inp.dat
        integer,parameter       :: un_result        = 26    ! result.dat
        integer,parameter       :: un_print         = 6     ! ��ʾ��
        real   , parameter      :: CFL_Value        = 0.8   ! CFL��
    end module const_Var
!--------------------------------------------------------------------------------------------------
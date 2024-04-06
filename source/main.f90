program main
    use const_Var
    use read_Parameter
    use read_Mesh
    use init
    use Boundary
    use residual
    use post_Process
    implicit none


    integer                     :: iStep         ! ����������
    integer                     :: subStep       ! RK�ӵ���
    integer                     :: m             ! ������
    type(block_Type),pointer    :: B             ! ��ָ��

    call read_Parameter_heat_transfer

    call read_Mesh_2D

    call init_Field         ! ��ʼ���¶ȳ�

   ! �������㿪ʼ

    ! �߽���Ϣ��������һ�����⣺û�������
    real_time = 0._dp
    call calculate_dt()                         ! ȷ��ʱ�䲽����δ��ɣ�
    do iStep = 1, total_Step
        real_Time = real_Time+dt
        call exchange_Inner_Boundary()          ! �ڱ߽���Ϣ����

        do m = 1, num_Block
            B => mesh(m)
            B%Un(:,:,:,:) = B%U(:,:,:,:)

            do subStep = 1, RKn

               ! �����ڲ��߽�����(δ���)
                call upgrade_Boundary()         
               ! ����вδ��ɣ�
                call compute_Residual()
               ! ���¿������δ��ɣ�
                B%U(:,:,:,:) = Ralpha(subStep)*B%Un(:,:,:,:)+Rgamma*B%U(:,:,:,:)+Rbeta*B%res(:,:,:,:)*dt
            end do
        end do
    end do

    
    
    ! call save_Tecplot1("../resource/output/result.plt")

    
    write(un_print,*) "Successful !!!"
    write(un_log  ,*) "Successful !!!"
    close(un_log)
    read(*,*) 
    
end program main
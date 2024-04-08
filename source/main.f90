program main
    use const_Var
    use read_Parameter
    use read_Mesh
    use init
    use Boundary
    use CFL
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
        ! ÿ��ʱ���ƽ������ݴ�ԭʼ������
        do m = 1, num_Block
            B => mesh(m)
            B%Un(:,:,:,:) = B%U(:,:,:,:)
        end do

        do subStep = 1, RKn                     ! ����-�����ƽ�
           ! �����п���¿���������(δ���)
            call upgrade_Ghost_Cell()           
           ! �����п����вδ��ɣ�
            call compute_Residual()             
        !    ! �������п���ڲ�������
        !     do m = 1, num_Block
        !         B => mesh(m)
        !         B%U(:,:,:,:) = Ralpha(subStep)*B%Un(:,:,:,:)+Rgamma(subStep)*B%U(:,:,:,:)+Rbeta(subStep)*B%res(:,:,:,:)*dt
        !     end do

        end do
 
    end do

    
    
    call save_Tecplot_Debug("../resource/output/result.plt")

    
    write(un_print,*) "Successful !!!"
    write(un_log  ,*) "Successful !!!"
    close(un_log)
    read(*,*) 
    
end program main



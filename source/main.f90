program main
    use const_Var
    use read_Parameter
    use read_Mesh
    use init
    use Boundary
    use residual
    use post_Process
    implicit none


    integer                     :: iStep         ! 迭代步索引
    integer                     :: subStep       ! RK子迭代
    integer                     :: m             ! 块索引
    type(block_Type),pointer    :: B             ! 块指针

    call read_Parameter_heat_transfer

    call read_Mesh_2D

    call init_Field         ! 初始化温度场

   ! 迭代计算开始

    ! 边界信息交换存在一个问题：没考虑清楚
    real_time = 0._dp
    call calculate_dt()                         ! 确定时间步长（未完成）
    do iStep = 1, total_Step
        real_Time = real_Time+dt
        call exchange_Inner_Boundary()          ! 内边界信息交换

        do m = 1, num_Block
            B => mesh(m)
            B%Un(:,:,:,:) = B%U(:,:,:,:)

            do subStep = 1, RKn

               ! 更新内部边界条件(未完成)
                call upgrade_Boundary()         
               ! 计算残差（未完成）
                call compute_Residual()
               ! 更新块变量（未完成）
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
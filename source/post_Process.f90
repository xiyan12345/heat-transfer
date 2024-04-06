!--------------------------------------------------------------------------------------------------
! 后处理文件
!--------------------------------------------------------------------------------------------------
    ! 后处理模块
    module post_Process
        use precision_EC
        use global_Var
        use const_Var
        implicit none
        
    contains
        subroutine save_Tecplot1(fileName)      ! teacher format -> tecplot
            implicit none
            character(len = file_Len)   :: fileName

           ! 中间变量
            integer                     :: m            ! 网格块索引
            integer                     :: i,j,k        ! 节点索引
            integer                     :: LAPs         ! 保存时的虚网格层数
            type(block_Type),pointer    :: B            ! 块指针

            LAPs = LAP
           ! tecplot格式：结构型数据（ORDERED)，数据组织形式（BLOCK）
           ! 储存数据：x（Nodal），y（Nodal），z（Nodal），T(CELLCENTERED)
            open( unit = un_result, file = fileName, status = "REPLACE")
            write(un_result,"(A24)")  "variables = ""x"" ""y"" ""z"" "
            do m = 1, num_Block
                B => mesh(m)
                write(un_result,100) B%nx+2*LAPs,B%ny+2*LAPs,B%nz+2*LAPs
                do k = 1-LAPs, B%nz+LAPs
                    do j = 1-LAPs, B%ny+LAPs
                        do i = 1-LAPs, B%nx+LAPs
                            write(un_result,"(3G24.15)") B%x(i,j,k),B%y(i,j,k),B%z(i,j,k)
                        end do
                    end do
                end do

            end do
            close(un_result) 


            ! 输入/输出格式
100          FORMAT("ZONE I = ", I3, ",J = ", I3, ",K = ", I3, ",DATAPACKING = POINT")
        end subroutine save_Tecplot1
        
    end module post_Process
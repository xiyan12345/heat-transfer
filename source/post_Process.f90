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

            LAPs = 1
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

        subroutine save_Tecplot(fileName)
            implicit none

           ! 输入\输出变量
            character(len=file_Len)     :: fileName

           ! 中间变量
            integer                     :: m            ! 块索引
            integer                     :: count        ! 计数器
            integer                     :: i,j,k        ! i,j,k维索引
            integer                     :: LAPs         ! 保存时的虚网格层数
            character(len=Lstr_Len)     :: zone_format  ! zone控制字符串
            type(block_Type),pointer    :: B            ! 块指针


            LAPs = 1
            open(unit = un_result, file = fileName, status = "REPLACE")
            write(un_result,"(A27)",advance="NO") "variables = ""x"",""y"",""z"",""T"""
            zone_format = "('zone I=',I??,',J=',I??,',K=',I??,',DATAPACKING=BLOCK,VARLOCATION=([4]=CELLCENTERED)')"
            do m = 1, num_Block


                B => mesh(m)
                
                write(zone_format(13:14),"(I2.2)") int(log10(real(B%nx+LAPs*2))+1)
                write(zone_format(23:24),"(I2.2)") int(log10(real(B%ny+LAPs*2))+1)
                write(zone_format(33:34),"(I2.2)") int(log10(real(B%nz+LAPs*2))+1)
                write(un_result,*)
                write(un_result,zone_format,advance="NO") B%nx+LAPs*2,B%ny+LAPs*2,B%nz+LAPs*2

               ! 写入x值
                count = 0
                do k = 1-LAPs, B%nz+LAPs
                    do j = 1-LAPs, B%ny+LAPs
                        do i = 1-LAPs, B%nx+LAPs
                            if (mod(count,3)==0) write(un_result,*)
                            count = count+1
                            write(un_result,"(G24.15)",advance="NO") B%x(i,j,k)
                        end do
                    end do
                end do

               ! 写入Y值
                count = 0
                do k = 1-LAPs, B%nz+LAPs
                    do j = 1-LAPs, B%ny+LAPs
                        do i = 1-LAPs, B%nx+LAPs
                            if (mod(count,3)==0) write(un_result,*)
                            count = count+1
                            write(un_result,"(G24.15)",advance="NO") B%y(i,j,k)
                        end do
                    end do
                end do

               ! 写入z值
                count = 0
                do k = 1-LAPs, B%nz+LAPs
                    do j = 1-LAPs, B%ny+LAPs
                        do i = 1-LAPs, B%nx+LAPs
                            if (mod(count,3)==0) write(un_result,*)
                            count = count+1
                            write(un_result,"(G24.15)",advance="NO") B%z(i,j,k)
                        end do
                    end do
                end do

               ! 写入T值
                count = 0
                do k = 1-LAPs, B%nz-1+LAPs
                    do j = 1-LAPs, B%ny-1+LAPs
                        do i = 1-LAPs, B%nx-1+LAPs
                            if (mod(count,3)==0) write(un_result,*)
                            count = count+1
                            write(un_result,"(G24.15)",advance="NO") B%U(1,i,j,k)
                        end do
                    end do
                end do

            end do
            close(un_result)


        end subroutine save_Tecplot
        
    end module post_Process
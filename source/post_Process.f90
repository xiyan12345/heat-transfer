!--------------------------------------------------------------------------------------------------
! �����ļ�
!--------------------------------------------------------------------------------------------------
    ! ����ģ��
    module post_Process
        use precision_EC
        use global_Var
        use const_Var
        implicit none
        
    contains
        subroutine save_Tecplot1(fileName)      ! teacher format -> tecplot
            implicit none
            character(len = file_Len)   :: fileName

           ! �м����
            integer                     :: m            ! ���������
            integer                     :: i,j,k        ! �ڵ�����
            integer                     :: LAPs         ! ����ʱ�����������
            type(block_Type),pointer    :: B            ! ��ָ��

            LAPs = 1
           ! tecplot��ʽ���ṹ�����ݣ�ORDERED)��������֯��ʽ��BLOCK��
           ! �������ݣ�x��Nodal����y��Nodal����z��Nodal����T(CELLCENTERED)
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


            ! ����/�����ʽ
100          FORMAT("ZONE I = ", I3, ",J = ", I3, ",K = ", I3, ",DATAPACKING = POINT")
        end subroutine save_Tecplot1

        subroutine save_Tecplot(fileName)
            implicit none

           ! ����\�������
            character(len=file_Len)     :: fileName

           ! �м����
            integer                     :: m            ! ������
            integer                     :: count        ! ������
            integer                     :: i,j,k        ! i,j,kά����
            integer                     :: LAPs         ! ����ʱ�����������
            character(len=Lstr_Len)     :: zone_format  ! zone�����ַ���
            type(block_Type),pointer    :: B            ! ��ָ��


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

               ! д��xֵ
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

               ! д��Yֵ
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

               ! д��zֵ
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

               ! д��Tֵ
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
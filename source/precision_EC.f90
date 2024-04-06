!--------------------------------------------------------------------------------------------------
! 精度控制文件
!--------------------------------------------------------------------------------------------------
! 精度控制模块
    module precision_EC
        implicit none

        integer,parameter       :: dp =8            ! 单精度
        integer,parameter       :: file_Len  = 64   ! 文件路径字符串长度
        integer,parameter       :: Sstr_Len  = 32   ! short型字符串长度
        integer,parameter       :: str_Len   = 64   ! 常规型字符串长度
        integer,parameter       :: Lstr_Len = 128   ! long型字符串长度
        
    end module precision_EC

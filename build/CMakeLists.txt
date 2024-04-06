cmake_minimum_required(VERSION 3.5)

#指定编译器
set(CMAKE_Fortran_COMPILER "D:/yanghao/code/C++/compiler/mingw64/bin/gfortran.exe")
#指定语言
enable_language(Fortran)
#指定编译器命令参数
set(compiler_args "-ffree-line-length-512")
set(CMAKE_Fortran_FLAGS ${compiler_args})

#项目入口
Project(Pr)

#指定项目目录
get_filename_component(project_path ${CMAKE_SOURCE_DIR} DIRECTORY)
#指定输出目录
set(EXECUTABLE_OUTPUT_PATH ${project_path}/bin)

#设置编译类型
set(CMAKE_BUILD_TYPE Debug)

#获取所有源文件
file(GLOB_RECURSE all_source ${project_path}/source/*.f90)

#生成可执行文件
add_executable(main ${all_source})
cmake_minimum_required(VERSION 2.6)
project(libexpm)
enable_language(C)


### Flags
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -O3 -fPIC")


### External dependencies
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

#find_package(OpenMP)
#if (OPENMP_FOUND)
#  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
#endif()


### Source files
set(LIBEXPM_FILES
  src/expm.c
  src/utils.c
  src/normest.c
)

file(GLOB src_h src/*.h)
execute_process(COMMAND echo ${src_h})


### Make and install
add_library(libexpm SHARED ${LIBEXPM_FILES})
add_library(expm_static STATIC ${LIBEXPM_FILES})


target_link_libraries (libexpm ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

set_target_properties(libexpm PROPERTIES PREFIX "")
if(NOT WIN32)
  set_target_properties(libexpm PROPERTIES SUFFIX ".so")
  install(TARGETS libexpm LIBRARY DESTINATION ${CMAKE_SOURCE_DIR}/build/lib)
else(NOT WIN32)
  install(TARGETS libexpm RUNTIME DESTINATION ${CMAKE_SOURCE_DIR}/build/lib)
endif(NOT WIN32)

install(FILES ${src_h} DESTINATION ${CMAKE_SOURCE_DIR}/build/include)

#INSTALL(FILES build/libglm.so DESTINATION lib)

cmake_minimum_required(VERSION 2.6)
project(expm)
enable_language(C)


### Flags
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wfatal-errors -O3 -fPIC")


### External dependencies
find_package(LAPACK REQUIRED)
find_package(BLAS REQUIRED)

#find_package(OpenMP)
#if (OPENMP_FOUND)
#  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
#endif()


### Source files
set(EXPM_FILES
  expm.c
  utils.c
)


### Make and install
add_library(expm SHARED ${EXPM_FILES})
add_library(expm_static STATIC ${EXPM_FILES})


target_link_libraries (expm ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

set_target_properties(expm PROPERTIES PREFIX "")
if(NOT WIN32)
  set_target_properties(expm PROPERTIES SUFFIX ".so")
  install(TARGETS expm LIBRARY DESTINATION ${CMAKE_SOURCE_DIR})
else(NOT WIN32)
  install(TARGETS expm RUNTIME DESTINATION ${CMAKE_SOURCE_DIR})
endif(NOT WIN32)


#INSTALL(FILES build/libglm.so DESTINATION lib)

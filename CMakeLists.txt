cmake_minimum_required(VERSION 2.8)

PROJECT(postprocess)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(postprocess MACOSX_BUNDLE postprocess)

if(VTK_LIBRARIES)
    target_link_libraries(postprocess ${VTK_LIBRARIES})
else()
    target_link_libraries(postprocess vtkHybrid vtkWidgets)
endif()

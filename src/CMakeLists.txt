cmake_minimum_required(VERSION 3.12)
project(numsim)

add_executable(${PROJECT_NAME}
        main.cpp
        Computation.cpp
        settings.cpp
        StaggeredGrid.h StaggeredGrid.cpp
        Discretization.cpp Discretization.h
        DonorCell.cpp DonorCell.h
        CentralDifferences.h CentralDifferences.cpp
        SOR.h SOR.cpp 
		GaussSeidel.h GaussSeidel.cpp
		SORRedBlach.h SORRedBlack.cpp
		GaussSeidelRedBlack.h GaussSeidelRedBlack.cpp
        PressureSolver.cpp PressureSolver.h
        output_writer/output_writer_paraview.cpp output_writer/output_writer_paraview.h
        output_writer/output_writer_text.cpp output_writer/output_writer_text.h
        output_writer/output_writer.cpp output_writer/output_writer.h
        Array2D.cpp Array2D.h
        FieldVariable.cpp FieldVariable.h)

target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_SOURCE_DIR})
message("Including ${PROJECT_SOURCE_DIR}")

# Search for the external package "VTK"
find_package(VTK)

# Output various cmake variables for demonstration purpose
message("If VTK was found on the system: VTK_FOUND: ${VTK_FOUND}")
message("The directory of VTK: VTK_DIR: ${VTK_DIR}")
message("The include directory of VTK: VTK_INCLUDE_DIRS: ${VTK_INCLUDE_DIRS}")

# If an installation of vtk was found on the system
if (VTK_FOUND)
    include_directories(${VTK_INCLUDE_DIRS})               # add the include directory where the header files are for the compiler
    target_link_libraries(${PROJECT_NAME} ${VTK_LIBRARIES}) # add the libraries for the linker
endif(VTK_FOUND)

install(TARGETS ${PROJECT_NAME} RUNTIME DESTINATION ${PROJECT_SOURCE_DIR}/../build)
# Add additonial compile options to enable more warnings
add_compile_options(-Wall -Wextra)

# Set the version of the C++ standard to use, we use C++14, published in 2014
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

message("CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")

cmake_minimum_required(VERSION 3.12)
project(numsim_parallel)

add_executable(${PROJECT_NAME}
        main.cpp
        Computation.cpp Computation.h
	ComputationParallel.cpp Computation.h
	Partitioning.cpp Partitioning.h
        settings.cpp settings.h
        StaggeredGrid.h StaggeredGrid.cpp
        Discretization.cpp Discretization.h
        DonorCell.cpp DonorCell.h
        CentralDifferences.h CentralDifferences.cpp
        SOR.h SOR.cpp
		GaussSeidel.h GaussSeidel.cpp
		SORRedBlack.h SORRedBlack.cpp
		GaussSeidelRedBlack.h GaussSeidelRedBlack.cpp
        PressureSolver.cpp PressureSolver.h
        output_writer/output_writer_paraview.cpp output_writer/output_writer_paraview.h
        output_writer/output_writer_text.cpp output_writer/output_writer_text.h
        PressureSolver.cpp PressureSolver.h
        output_writer/output_writer_paraview_parallel.cpp output_writer/output_writer_paraview_parallel.h
        output_writer/output_writer_text_parallel.cpp output_writer/output_writer_text_parallel.h
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

find_package(MPI REQUIRED)

include_directories(${MPI_INCLUDE_PATH})
target_link_libraries(${PROJECT_NAME} ${MPI_LIBRARIES})

if(MPI_COMPILE_FLAGS)
  set_target_properties(${PROJECT_NAME} PROPERTIES
    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(${PROJECT_NAME} PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()

install(TARGETS ${PROJECT_NAME} RUNTIME DESTINATION ${PROJECT_SOURCE_DIR}/../build)
# Add additonial compile options to enable more warnings
add_compile_options(-Wall -Wextra)

# Set the version of the C++ standard to use, we use C++14, published in 2014
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

message("CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")

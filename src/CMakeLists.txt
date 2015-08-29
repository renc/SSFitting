

set (PM_SRC 
DecimationModel.cpp
DecimationModel.h
DecimationModel_1_sim_para.cpp
DecimationModel_2_resample.cpp
DecimationModel_3_fitting.cpp
DecimationModel_4_fitting_error.cpp
DecimationModel_5_fitting_error_byQAS.cpp
DecimationModel_boundingbox.cpp
DecimationModel_errerRefine.cpp
GUIFrame.cpp
GUIFrame.h
main.cpp
MeshModel.cpp
MeshModel.h
myreadme.txt
OpenMeshAll.h
OpenMeshMeshType.h
QuadraticBezierTriangles.cpp
QuadraticBezierTriangles.h 
VHierarchy.h

EquationSolver_lgq/EquationSolver.cpp
EquationSolver_lgq/EquationSolver.h 
EquationSolver_lgq/xVec3.h 

lib_IDSSsrc/bnd_loop.cpp
lib_IDSSsrc/bnd_loop.h 
lib_IDSSsrc/loopdomain.h 
lib_IDSSsrc/looprange.h 

${COMMON_DIR}/CLoop.cpp 
${COMMON_DIR}/Hoppe94LoopEval.cpp
${COMMON_DIR}/Stam98ClassicLoopEval.cpp 


) 

add_definitions(-D_USE_MATH_DEFINES) #openmesh need this 

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /wd4018")  
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4018") # # warning C4018: '<' : signed/unsigned mismatch 
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /wd4244")  
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4244") # 
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /wd4267")  
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4267") 
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /wd4819")   ##warning C4819: The file contains a character that cannot be represented in the current code page (936). Save the file in Unicode format to prevent data loss 
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4819")

add_executable(pm ${PM_SRC})
target_link_libraries(pm
	${TAUSE_LIB}
	${OPENMESH_LIB_NAME}
	${FLTK_LIB}
	OpenGL32.lib 
	glu32.lib
)
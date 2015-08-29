#ifndef dgpstudio_common_openmesh_all_h
#define dgpstudio_common_openmesh_all_h

//这里包含了所有需要用到的openmesh library的头文件, 作为openmesh library的入口.
//目的是方便当openmesh update时候只要编辑这个文件就ok了.

//#ifdef _OPENMESH_1_1_0
// version 1.1.0
//#include <OpenMesh/Core/IO/MeshIO.hh> 
//#include <OpenMesh/Core/Math/VectorT.hh>
//#include <OpenMesh/Core/Mesh/Types/TriMesh_ArrayKernelT.hh>
//#include <OpenMesh/Core/Mesh/Types/PolyMesh_ArrayKernelT.hh>
//#include <OpenMesh/Core/Utils/Noncopyable.hh>
//
//
//#include <OpenMesh/Tools/Utils/Timer.hh>
//#include <OpenMesh/Tools/Subdivider/Uniform/SubdividerT.hh>
//#include <OpenMesh/Tools/Subdivider/Uniform/LoopT.hh>  
//
//
//#pragma comment (lib, "OpenMeshCoreMD.lib") // only release MD
//#pragma comment (lib, "OpenMeshToolsMD.lib")// only release MD

//#else //_OPENMESH_2_0
//// version 2.0
#include <OpenMesh/Core/IO/MeshIO.hh> 
#include <OpenMesh/Core/Geometry/VectorT.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Utils/Noncopyable.hh>


#include <OpenMesh/Tools/Utils/Timer.hh>
#include <OpenMesh/Tools/Subdivider/Uniform/SubdividerT.hh>
#include <OpenMesh/Tools/Subdivider/Uniform/LoopT.hh>  

//#pragma comment (lib, "libOpenMeshCore.lib") // only release MD
//#pragma comment (lib, "libOpenMeshTools.lib")// only release MD
//
//#endif // 


#ifdef _MSC_VER
//#   pragma comment(lib, "zlib")
//#   pragma comment(lib, "ws2_32")
//#   pragma comment(lib, "winmm")
//#   pragma comment(lib, "imagehlp")
//#   pragma comment(lib, "gdi32")
//#   pragma comment(lib, "user32")
//#   pragma comment(lib, "kernel32")
//#   pragma comment(lib, "version")
//#   pragma comment(lib, "advapi32")
//#   pragma comment(lib, "png")
//#   pragma comment(lib, "jpeg")
//#   pragma comment(lib, "zip")
//#   ifdef _DEBUG
//// Don't link against G3D when building G3D itself.
//#      ifndef G3D_BUILDING_LIBRARY_DLL
//#         pragma comment(lib, "G3Dd.lib")
//#      endif
//#   else
//// Don't link against G3D when building G3D itself.
//#      ifndef G3D_BUILDING_LIBRARY_DLL
//#         pragma comment(lib, "G3D.lib")
//#      endif
//#   endif


#endif // _MSC_VER











#endif // dgpstudio_common_openmesh_all_h
# SSFitting
Subdivision Surface Fitting

What?
This is an old project back to 2007 at school, almost 8 years later it is built fixed and updated for references. 

Steps: (almost listed at the Edit menu)
decimate the mesh, build a progressive mesh during simplification, the pm is used during the refinement step. 
resample and subd the simplified mesh, build the fitting equation, compute and get the initial control mesh, calc the fitting error. 

## Build
####Dependent libraries: 
Fltk-1.3.0, Cmake is used to create visual studio project file, using v11 (visual studio 2012) and 32 bits setting. Build the fltk and fltkgl projects, with setting win32 and MDd. The output fltk.lib and fltkgl.lib are located at Fltk-1.3.0/lib folder. 
The Fltk has opengl/glu/glut built in, so we do not need freeglut library. <Fl/gl.h, glu.h, glut.h> is used, rather than <GL/gl.h, glu.h, glut.h>. 

OpenMesh-2.0-RC3, Cmake is used too, openmesh-2.0-RC3/lib/OpenMeshCore.lib OpenMeshTools.lib. 

Taucs, there is a taucs.sln solution file can be used, link the libararies in taucs-win32, the output is taucs.ib

####SSFitting: 
Cmake to create project file and build.   

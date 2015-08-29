
> 2008.01.02
附加库目录中去除
..\..\..\..\libgfx-1.1.0\vc7\libgfx\Release
附加依赖项中去除
libgfx.lib jpeg.lib png.lib zlib.lib zip.lib

> 2007.12.20
利用MVC的structure来将程序结果弄好些。

Error是GUIFrame中可以含有MeshModel的指针, 但是想在MeshModel中也有指向GUIFrame的指针就出错了,不会实现.
这是GUIFrame.h和MeshModel.h头文件相互引用的错误, 其中一个class A(或两个)的头文件中不包含另一个class B的头文件, 
而只是用前置声明, 如A.h在定义class A前用extern class B; 但是在A.cpp中再包含B.h.

*View/Controller*												
GUIFrame	所有的gui都在此做, 和mesh有关的event发送给Model来处理.

*Model*
MeshModel	对mesh的操作都在这里完成.


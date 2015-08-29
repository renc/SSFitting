#ifndef DGP_BASEMESH_H
#define DGP_BASEMESH_H


namespace DGP {
	
	enum MeshTypes
	// just a interface
	class BaseMesh {
	public:
		BaseMesh() {}
		virtual ~BaseMesh() = 0;
	};
}
#endif //DGP_BASEMESH_H
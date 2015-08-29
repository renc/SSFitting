#ifndef dgp_geometry_object_h
#define dgp_geometry_object_h

namespace DGP {
	
	class BaseMesh; 

	class GeometryObject {

	public:
		GeometryObject() : base_mesh_(NULL) 
		{
			 
		}
		~GeometryObject() {
			if (base_mesh_ != NULL) {
				delete base_mesh_;
				base_mesh_ = NULL;
			}
		}
		//

	private:
		BaseMesh * base_mesh_;
	};
}
#endif //dgp_geometry_object_h
#include <iostream>
#include <geodesic/polygonal_addressing.hh>

using namespace Spectrum;

constexpr int verts_per_face = 3;
constexpr int len = 8;

int main(int argc, char *argv[])
{
// These always work
   std::cerr << "Edge count for len = " << len << " is " << PolygonalAddressing<verts_per_face>::EdgeCount(len) << "\n";
   std::cerr << "Face count for len = " << len << " is " << PolygonalAddressing<verts_per_face>::FaceCount(len) << "\n";
   std::cerr << "Vert count for len = " << len << " is " << PolygonalAddressing<verts_per_face>::VertCount(len) << "\n";

// Test if our static initialization worked
   std::cerr << "vert_vert[0][0] = " << PolygonalAddressing<verts_per_face>::GetVertVert(0, 0) << "\n";
   std::cerr << "vert_vert[1][0] = " << PolygonalAddressing<verts_per_face>::GetVertVert(1, 0) << "\n";
   std::cerr << "vert_vert[2][0] = " << PolygonalAddressing<verts_per_face>::GetVertVert(2, 0) << "\n";
};


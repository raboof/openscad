/*
 *  OpenSCAD (www.openscad.org)
 *  Copyright (C) 2009-2011 Clifford Wolf <clifford@clifford.at> and
 *                          Marius Kintel <marius@kintel.net>
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  As a special exception, you have permission to link this program
 *  with the CGAL library and distribute executables, as long as you
 *  follow the requirements of the GNU GPL in regard to all of the
 *  software in the executable aside from CGAL.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

#include "export.h"
#include "polyset.h"
#include "polyset-utils.h"
#include "GeometryUtils.h"
#include "dxfdata.h"

#ifdef ENABLE_CGAL
#include "CGAL_Nef_polyhedron.h"
#include "cgal.h"
#include "cgalutils.h"

#define QUOTE(x__) # x__
#define QUOTED(x__) QUOTE(x__)

///// Tools for ASCII Decimal based formats (OBJ, STL, AMF, OFF, &c)
// 3d coordinate represented in ASCII base-10 Decimal
// example = "5.2 12.1 13.0" x=5.2 y=12.1 z=13.0
typedef std::string ascii_vertex;
// sequence of ascii_vertex, can be more than 3.
// example: {"3 4 5","5.2 12.1 13.0","0 0 0"}
typedef std::vector<ascii_vertex> ascii_face;
// triangle: convenience struct for ascii_faces with three points.
typedef struct ascii_triangle {
    ascii_vertex vs1;
    ascii_vertex vs2;
    ascii_vertex vs3;
} ascii_triangle_t;

// backwards compatabile
typedef ascii_triangle_t triangle;

/* Convert PolySet to sequence of ASCII coordinate vertexes and faces.
  Can produce faces with more than three points.  */
void PolySet_to_ASCII_Faces(const PolySet &ps, std::vector<ascii_vertex> &vertices, std::vector<ascii_face> &faces) {
	vertices.clear();
	faces.clear();
	std::map<ascii_vertex,int> vertmap;
	for (size_t i = 0; i < ps.polygons.size(); i++) {
		const Polygon *poly = &ps.polygons[i];
		ascii_face face;
		std::map<ascii_vertex,int> dup_detect;
		for (size_t j = 0; j < poly->size(); j++) {
			Vector3d v = poly->at(j);
			std::stringstream stream;
			stream << v.transpose();
			ascii_vertex coord3d(stream.str());
			if (vertmap.count(coord3d)==0) {
				vertmap[coord3d] = vertices.size();
				vertices.push_back(coord3d);
			}
			dup_detect[coord3d] = 1;
			face.push_back(coord3d);
		}
		if (dup_detect.size() == poly->size()) {
			// polygon face doesn't contain duplicate points
			faces.push_back(face);
		}
	}
}

// given a list of ASCII decimal 3d coordinates, find the given coordinate index
// Example: find '0.1 3.4 1.2' in '{{0,2,1},{1,2,2},{0.1,3.4,1.2}}' returns 2.
size_t find_index(std::vector<ascii_vertex> &vertices, ascii_vertex tofind) {
		return std::distance(vertices.begin(), std::find(vertices.begin(), vertices.end(), tofind));
}

/*!
   Converts the given 3d CGAL Polyhedron (non-nef) to ASCII coordinate triangles.
   It is assumed the input Polyhedron is all triangle faces.
*/
void CGALPolyhedron_to_ASCII_Triangles(const CGAL_Polyhedron &P, std::vector<ascii_vertex> &vertices, std::vector<ascii_face> &faces)
{
	typedef CGAL_Polyhedron::Vertex Vertex;
	typedef CGAL_Polyhedron::Vertex_const_iterator VCI;
	typedef CGAL_Polyhedron::Facet_const_iterator FCI;
	typedef CGAL_Polyhedron::Halfedge_around_facet_const_circulator HFCC;

	for (FCI fi = P.facets_begin(); fi != P.facets_end(); ++fi) {
		HFCC hc = fi->facet_begin();
		HFCC hc_end = hc;
		Vertex v1, v2, v3;
		v1 = *VCI((hc++)->vertex());
		v3 = *VCI((hc++)->vertex());
		do {
			v2 = v3;
			v3 = *VCI((hc++)->vertex());
			double x1 = CGAL::to_double(v1.point().x());
			double y1 = CGAL::to_double(v1.point().y());
			double z1 = CGAL::to_double(v1.point().z());
			double x2 = CGAL::to_double(v2.point().x());
			double y2 = CGAL::to_double(v2.point().y());
			double z2 = CGAL::to_double(v2.point().z());
			double x3 = CGAL::to_double(v3.point().x());
			double y3 = CGAL::to_double(v3.point().y());
			double z3 = CGAL::to_double(v3.point().z());
			std::stringstream stream;
			stream << x1 << " " << y1 << " " << z1;
			ascii_vertex vs1 = stream.str();
			stream.str("");
			stream << x2 << " " << y2 << " " << z2;
			ascii_vertex vs2 = stream.str();
			stream.str("");
			stream << x3 << " " << y3 << " " << z3;
			ascii_vertex vs3 = stream.str();
			if (std::find(vertices.begin(), vertices.end(), vs1) == vertices.end())
				vertices.push_back(vs1);
			if (std::find(vertices.begin(), vertices.end(), vs2) == vertices.end())
				vertices.push_back(vs2);
			if (std::find(vertices.begin(), vertices.end(), vs3) == vertices.end())
				vertices.push_back(vs3);

			if (vs1 != vs2 && vs1 != vs3 && vs2 != vs3) {
				// The above condition ensures that there are 3 distinct vertices, but
				// they may be collinear. If they are, the unit normal is meaningless
				// so the default value of "1 0 0" can be used. If the vertices are not
				// collinear then the unit normal must be calculated from the
				// components.
				ascii_face tri;
				tri.push_back(vs1);
				tri.push_back(vs2);
				tri.push_back(vs3);
				faces.push_back(tri);
			}
		} while (hc != hc_end); // for each triangle in the Facet
	} // for each Polyhedron Facet
}

/* Convert CGAL Nef Polyhedron to ASCII decimal coordiante vertexes and
 triangles. Only produces faces with 3 points. As of writing, this works
 because it first converts to CGAL Polyhedron (non-nef), which makes all
 faces triangular. Then it converts to ASCI coordinates. */
void NefPoly_to_ASCII_Triangles(const CGAL_Nef_polyhedron &root_N, std::vector<ascii_vertex> &vertices, std::vector<ascii_face> &faces)
{
	// since conversion to Polyhedron can fail, we test here first.
	if (!root_N.p3->is_simple()) {
		PRINT("Object isn't a valid 2-manifold! Modify your design.");
		return;
	}

	CGAL::Failure_behaviour old_behaviour = CGAL::set_error_behaviour(CGAL::THROW_EXCEPTION);
	try {
		CGAL_Polyhedron P;
		root_N.p3->convert_to_Polyhedron(P);
		CGALPolyhedron_to_ASCII_Triangles(P, vertices, faces);
	} catch (CGAL::Assertion_exception e) {
		PRINTB("CGAL error in CGAL_Nef_polyhedron3::convert_to_Polyhedron(): %s", e.what());
	}
	CGAL::set_error_behaviour(old_behaviour);
}

void ASCII_Faces_to_obj(std::vector<ascii_vertex> &vertices, std::vector<ascii_face> &faces, std::ostream &output)
{
	setlocale(LC_NUMERIC, "C"); // Ensure radix is . (not ,) in output
	fprintf(stderr, "outputting\n");
	output << "# WaveFront *.obj file (generated by OpenSCAD "
		<< QUOTED(OPENSCAD_VERSION)
#ifdef OPENSCAD_COMMIT
		<< " (git " << QUOTED(OPENSCAD_COMMIT) << ")"
#endif
		<< ")\n\n"
		<< "g Object\n";
	for (size_t i = 0; i < vertices.size(); i++) {
		std::string x,y,z;
		std::istringstream stream(vertices[i]);
		if (!(stream >> x >> y >> z)) return;
		output << "v " << x << " " << y << " " << z << "\n";
	}
	output << "\n";
	for (size_t i = 0; i < faces.size(); i++) {
		ascii_face face = faces[i];
		output << "f ";
		for (size_t j = 0; j < face.size(); j++) {
			output << 1 + find_index(vertices, face[j]) << " ";
		}
		output << "\n";
	}
	output << "\n# end WaveFront *.obj file (generated by OpenSCAD "
		<< QUOTED(OPENSCAD_VERSION) << ")\n";
	setlocale(LC_NUMERIC, ""); // Set default locale
}

void append_obj(const shared_ptr<const Geometry> &geom, std::vector<ascii_vertex> &vertices, std::vector<ascii_face> &faces) 
{
	if (const auto geomlist = dynamic_pointer_cast<const GeometryList>(geom)) {
		for(const auto &item : geomlist->getChildren()) {
			append_obj(item.second, vertices, faces);
		}
	}
	if (const auto N = dynamic_pointer_cast<const CGAL_Nef_polyhedron>(geom)) {
		if (!N->isEmpty()) NefPoly_to_ASCII_Triangles(*N, vertices, faces);
	}
	else if (const auto ps = dynamic_pointer_cast<const PolySet>(geom)) {
		// FIXME: Implement this without creating a Nef polyhedron
		CGAL_Nef_polyhedron *N = CGALUtils::createNefPolyhedronFromGeometry(*ps);
		if (!N->isEmpty()) NefPoly_to_ASCII_Triangles(*N, vertices, faces);
		delete N;
	}
	else if (dynamic_pointer_cast<const Polygon2d>(geom)) {
		assert(false && "Unsupported file format");
	} else {
		assert(false && "Not implemented");
	}
}

void export_obj(const shared_ptr<const Geometry> &geom, std::ostream &output)
{
	std::vector<ascii_vertex> vertices;
	std::vector<ascii_face> triangles;

	append_obj(geom, vertices, triangles);

	ASCII_Faces_to_obj(vertices, triangles, output);
}

#endif // ENABLE_CGAL

#pragma once

#include <iostream>
#include <functional>

#include <boost/range/algorithm.hpp>
#include <boost/range/adaptor/map.hpp>

#include "Tree.h"
#include "Camera.h"
#include "memory.h"

class PolySet;

enum class FileFormat {
	STL,
	OFF,
	AMF,
	_3MF,
	OBJ,
	DXF,
	SVG,
	NEFDBG,
	NEF3,
	CSG,
	AST,
	TERM,
	ECHO,
	PNG
};

bool canPreview(const FileFormat format);

void exportFileByName(const shared_ptr<const class Geometry> &root_geom, FileFormat format,
											const char *name2open, const char *name2display);

void export_stl(const shared_ptr<const Geometry> &geom, std::ostream &output);
void export_3mf(const shared_ptr<const Geometry> &geom, std::ostream &output);
void export_off(const shared_ptr<const Geometry> &geom, std::ostream &output);
void export_amf(const shared_ptr<const Geometry> &geom, std::ostream &output);
void export_obj(const shared_ptr<const Geometry> &geom, std::ostream &output);
void export_dxf(const shared_ptr<const Geometry> &geom, std::ostream &output);
void export_svg(const shared_ptr<const Geometry> &geom, std::ostream &output);
void export_nefdbg(const shared_ptr<const Geometry> &geom, std::ostream &output);
void export_nef3(const shared_ptr<const Geometry> &geom, std::ostream &output);

// void exportFile(const class Geometry *root_geom, std::ostream &output, FileFormat format);

enum class Previewer { OPENCSG, THROWNTOGETHER };
enum class RenderType { GEOMETRY, CGAL, OPENCSG, THROWNTOGETHER };

struct ExportFileFormatOptions {
	const std::map<const std::string, FileFormat> exportFileFormats{
		{"stl", FileFormat::STL},
		{"off", FileFormat::OFF},
		{"amf", FileFormat::AMF},
		{"3mf", FileFormat::_3MF},
		{"obj", FileFormat::OBJ},
		{"dxf", FileFormat::DXF},
		{"svg", FileFormat::SVG},
		{"nefdbg", FileFormat::NEFDBG},
		{"nef3", FileFormat::NEF3},
		{"csg", FileFormat::CSG},
		{"ast", FileFormat::AST},
		{"term", FileFormat::TERM},
		{"echo", FileFormat::ECHO},
		{"png", FileFormat::PNG},
	};
};

struct ViewOption {
	const std::string name;
	bool& value;
};

struct ViewOptions {
	Previewer previewer{Previewer::OPENCSG};
	RenderType renderer{RenderType::OPENCSG};

	std::map<std::string, bool> flags{
		{"axes", false},
		{"scales", false},
		{"edges", false},
		{"wireframe", false},
		{"crosshairs", false},
	};

	const std::vector<std::string> names() {
		std::vector<std::string> names;
		boost::copy(flags | boost::adaptors::map_keys, std::back_inserter(names));
		return names;
	}

	bool &operator[](const std::string &name) {
		return flags.at(name);
	}

	bool operator[](const std::string &name) const {
		return flags.at(name);
	}
	
};

class OffscreenView;

std::unique_ptr<OffscreenView> prepare_preview(Tree &tree, const ViewOptions& options, Camera camera);
bool export_png(const shared_ptr<const class Geometry> &root_geom, const ViewOptions& options, Camera camera, std::ostream &output);
bool export_png(const OffscreenView &glview, std::ostream &output);

namespace Export {

struct Triangle {
	std::array<int, 3> key;
	Triangle(int p1, int p2, int p3)
	{
		// sort vertices with smallest value first without
		// changing winding order of the triangle.
		// See https://github.com/nophead/Mendel90/blob/master/c14n_stl.py

        if (p1 < p2) {
            if (p1 < p3) {
                key = {p1, p2, p3}; // v1 is the smallest
			} else {
                key = {p3, p1, p2}; // v3 is the smallest
			}
		} else {
            if (p2 < p3) {
                key = {p2, p3, p1}; // v2 is the smallest
            } else {
                key = {p3, p1, p2}; // v3 is the smallest
			}
		}
	}
};

class ExportMesh {
public:
	ExportMesh(const PolySet &ps);

	bool foreach_vertex(const std::function<bool(const std::array<double, 3>&)> callback) const;
	bool foreach_triangle(const std::function<bool(const std::array<int, 3>&)> callback) const;

private:
	std::map<std::array<double, 3>, int> vertexMap;
	std::vector<Triangle> triangles;
};

}

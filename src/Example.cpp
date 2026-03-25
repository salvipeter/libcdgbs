#include "libcdgbs/Example.hpp"
#include <iostream>

namespace libcdgbs {
  double bounding_box_diagonal(const Mesh& mesh) {
    double min_x = std::numeric_limits<double>::max();
    double max_x = std::numeric_limits<double>::lowest();
    double min_y = std::numeric_limits<double>::max();
    double max_y = std::numeric_limits<double>::lowest();
    double min_z = std::numeric_limits<double>::max();
    double max_z = std::numeric_limits<double>::lowest();

    for (const auto& v : mesh.vertices()) {
      const auto& pt = mesh.point(v);
      min_x = std::min(min_x, static_cast<double>(pt[0]));
      max_x = std::max(max_x, static_cast<double>(pt[0]));
      min_y = std::min(min_y, static_cast<double>(pt[1]));
      max_y = std::max(max_y, static_cast<double>(pt[1]));
      min_z = std::min(min_z, static_cast<double>(pt[2]));
      max_z = std::max(max_z, static_cast<double>(pt[2]));
    }

    return sqrt(pow(max_x - min_x, 2) + pow(max_y - min_y, 2) + pow(max_z - min_z, 2));
  }


  void Example::say_hello(const std::string& filename, const SurfGBS::InputParams& params) {
    gbs = SurfGBS();
    gbs.debug_outputs = true;
    gbs.readMGBS(filename + ".mgbs", params);
    // gbs.readGBS(filename + ".gbs", params);
    gbs.compute_domain_boundary();
    gbs.compute_domain_mesh();
    gbs.writeOBJ(gbs.meshDomain, filename + ".obj");
    gbs.compute_local_parameters();

    double diag = bounding_box_diagonal(gbs.meshDomain);
    for (size_t loop = 0; loop < gbs.num_loops; ++loop) {
      for (size_t side = 0; side < gbs.num_sides[loop]; ++side) {
        for (auto v : gbs.meshDomain.vertices()) {
          auto h = gbs.h_coords[v.idx()][loop][side];
          gbs.meshDomain.point(v)[2] = h*diag;
        }
        gbs.writeOBJ(gbs.meshDomain, std::string("h") + std::to_string(loop) + std::to_string(side) + std::string(".obj"));
        for (auto v : gbs.meshDomain.vertices()) {
          auto s = gbs.s_coords[v.idx()][loop][side];
          gbs.meshDomain.point(v)[2] = s*diag;
        }
        gbs.writeOBJ(gbs.meshDomain, std::string("s") + std::to_string(loop) + std::to_string(side) + std::string(".obj"));
      }
    }
    for (auto v : gbs.meshDomain.vertices()) {
      gbs.meshDomain.point(v)[2] = 0.0;
    }

    gbs.compute_blend_functions();
    gbs.evaluate_mesh(gbs.meshSurface, true, SurfGBS::GBS);
    gbs.writeOBJ(gbs.meshSurface, filename + "_surf.obj");
    gbs.writeMGBS(filename + "_out.mgbs");

    // std::cout << "Hello from libcdgbs!" << std::endl;
  }
}

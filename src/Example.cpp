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


  void Example::say_hello(const std::string& filename, const SurfGBS::InputParams& params_) {
    SurfGBS::InputParams params = params_;
    auto config = toml::parse(filename);
    params.target_length = toml::find<double>(config, "triangle_edge_length");
    gbs = SurfGBS();
    gbs.config = config;
    gbs.debug_outputs = true;
    auto fname = toml::find<std::string>(config, "filename");
    if (fname.ends_with(".mgbs"))
      gbs.readMGBS(fname, params);
    else if (fname.ends_with(".gbs"))
      gbs.readGBS(fname, params);
    else if (fname.ends_with(".cgb"))
      gbs.readCGB(fname, params);
    else {
      std::cerr << "Unknown file format!" << std::endl;
      return;
    }
    gbs.compute_domain_boundary();
    gbs.compute_domain_mesh();
    gbs.compute_local_parameters();

    auto export_domain = toml::find<bool>(config, "export_domain");
    auto export_h = toml::find<bool>(config, "export_h");
    auto export_s = toml::find<bool>(config, "export_s");
    auto export_patch = toml::find<bool>(config, "export_patch");

    if (export_domain)
      gbs.writeOBJ(gbs.meshDomain, filename + ".obj");

    for (size_t loop = 0; loop < gbs.num_loops; ++loop) {
      for (size_t side = 0; side < gbs.num_sides[loop]; ++side) {
        if (export_h) {
          for (auto v : gbs.meshDomain.vertices()) {
            auto h = gbs.h_coords[v.idx()][loop][side];
            gbs.meshDomain.point(v)[2] = h;
          }
          gbs.writeOBJ(gbs.meshDomain, std::string("h") + std::to_string(loop) + std::to_string(side) + std::string(".obj"));
        }
        if (export_s) {
          for (auto v : gbs.meshDomain.vertices()) {
            auto s = gbs.s_coords[v.idx()][loop][side];
            gbs.meshDomain.point(v)[2] = s;
          }
          gbs.writeOBJ(gbs.meshDomain, std::string("s") + std::to_string(loop) + std::to_string(side) + std::string(".obj"));
        }
      }
    }
    for (auto v : gbs.meshDomain.vertices()) {
      gbs.meshDomain.point(v)[2] = 0.0;
    }

    if (export_patch) {
      gbs.compute_blend_functions();
      gbs.evaluate_mesh(gbs.meshSurface, true, SurfGBS::GBS);
      gbs.writeOBJ(gbs.meshSurface, filename + "_surf.obj");
    }
  }
}

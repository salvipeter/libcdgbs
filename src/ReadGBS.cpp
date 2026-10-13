#include "libcdgbs/SurfGBS.hpp"
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <array>

using namespace libcdgbs;

bool SurfGBS::readGBS(const std::string& filename, const InputParams& params)
{
  std::ifstream in(filename);
  if (!in) {
    std::cerr << "Cannot open file\n";
    return false;
  }

  int cw_type, dom_type;
  double par_dil, dom_tol;
  in >> cw_type >> dom_type >> par_dil >> dom_tol;

  std::vector<std::vector<Ribbon> > ribbons;
  const size_t num_loops = 1;
  std::vector<size_t> num_sides;
  std::vector<std::vector<size_t>> num_rows;
  std::vector<std::vector<size_t>> num_cols;

  ribbons.clear();
  num_sides.clear();
  num_rows.clear();
  num_cols.clear();

  ribbons.resize(1);
  num_sides.resize(1);
  num_rows.resize(1);
  num_cols.resize(1);

  in >> num_sides[0];                // <# of ribbons>
  num_rows[0].resize(num_sides[0]);
  num_cols[0].resize(num_sides[0]);

  Eigen::Vector3d ccp;
  in >> ccp(0) >> ccp(1) >> ccp(2);   // <ccp_x> <ccp_y> <ccp_z>

  // The trailer (after all sides; see gbs_format.txt) can flag the
  // LAST layer of a side as a "Ribbon CP" row that is NOT part of the
  // blend surface, so the sides must be buffered before the ribbons
  // are constructed:
  //   #
  //   <RibCP weight> x n
  //   <RibCP flag (0/1)> x n
  //   <curve res> <mesh res> <harmonic levels>   (ignored)
  struct RawSide {
    int degU = 3;
    int layers = 1;
    Geometry::DoubleVector knotsU;
    std::vector<std::array<double, 3>> pts; // row-major (layer rows)
  };
  std::vector<RawSide> raw(num_sides[0]);

  for (size_t side = 0; side < num_sides[0]; ++side) {
    RawSide& r = raw[side];
    int degS, degH;
    // Spec order (gbs_format.txt): s-degree, h-degree, layer count.
    // The s-degree is the U (along-boundary) degree; the vertical
    // degree follows the layer count (degH is informational).
    in >> degS >> degH >> r.layers;
    r.degU = degS;

    // read full knot‐vector line for U
    in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    std::string knotLine;
    std::getline(in, knotLine);
    std::istringstream ksi(knotLine);
    r.knotsU = Geometry::DoubleVector{
        std::istream_iterator<double>(ksi),
        std::istream_iterator<double>()
    };

    const size_t ncols = r.knotsU.size() - r.degU - 1;
    r.pts.resize(ncols * r.layers);
    for (size_t i = 0; i < size_t(r.layers); ++i) {
      for (size_t j = 0; j < ncols; ++j) {
        double x, y, z;
        in >> x >> y >> z;
        r.pts[i * ncols + j] = { x, y, z };
      }
    }
  }

  // Optional trailer: '#', RibCP weights, RibCP flags. A flagged
  // side's last layer is dropped (Ribbon CPs ignored for now).
  std::vector<int> ribbon_cp(num_sides[0], 0);
  {
    std::string marker;
    if (in >> marker && marker == "#") {
      std::vector<double> ribcp_weights(num_sides[0], 1.0);
      bool ok = true;
      for (size_t i = 0; ok && i < num_sides[0]; ++i)
        ok = bool(in >> ribcp_weights[i]);
      for (size_t i = 0; ok && i < num_sides[0]; ++i)
        ok = bool(in >> ribbon_cp[i]);
      if (!ok)
        std::fill(ribbon_cp.begin(), ribbon_cp.end(), 0);
    }
  }

  for (size_t side = 0; side < num_sides[0]; ++side) {
    const RawSide& r = raw[side];
    const size_t ncols = r.knotsU.size() - r.degU - 1;
    const int rows = std::max(1, r.layers - (ribbon_cp[side] ? 1 : 0));
    const int degV = rows - 1;

    Geometry::DoubleVector knotsV;
    knotsV.insert(knotsV.end(), degV + 1, 0.0);
    knotsV.insert(knotsV.end(), degV + 1, 1.0);

    num_cols[0][side] = ncols;
    num_rows[0][side] = rows;

    Geometry::PointVector ctrl;
    ctrl.resize(ncols * size_t(rows));
    for (size_t i = 0; i < size_t(rows); ++i)
      for (size_t j = 0; j < ncols; ++j) {
        const auto& q = r.pts[i * ncols + j];
        ctrl[j * size_t(rows) + i] = { q[0], q[1], q[2] };
      }

    ribbons[0].emplace_back(r.degU, degV, r.knotsU, knotsV, ctrl);
  }

  load_ribbons(ribbons, params);

  return true;
}
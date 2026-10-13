#include "libcdgbs/SurfGBS.hpp"
#include "libcdgbs/MatrixUtil.hpp"

#include <algorithm>

using namespace libcdgbs;
using namespace MatrixUtil;

template <typename T>
void concatenateVectors(
  const std::vector<T>& v1,
  const std::vector<T>& v2,
  std::vector<T>& vc,
  bool                  common_joint = false
)
{
  const size_t s1 = v1.size();
  const size_t s2 = v2.size();
  const size_t shift = (common_joint ? 1 : 0);

  if (s1 == 0 && s2 == 0) {
    return;
  }

  vc.clear();
  vc.resize(s1 + s2 - shift);

  for (size_t ii = 0; ii < s1; ++ii) {
    vc[ii] = v1[ii];
  }
  for (size_t ii = 0; ii < s2; ++ii) {
    if (ii >= shift) {
      vc[s1 + ii - shift] = v2[ii];
    }
  }
}

bool SurfGBS::compute_harmonic_parameters()
{
  //Initialize
  s_coords.clear();
  s_coords.resize(meshDomain.n_vertices());
  h_coords.clear();
  h_coords.resize(meshDomain.n_vertices());
  for (const auto v : meshDomain.vertices()) {
    auto& s = s_coords[v.idx()];
    auto& h = h_coords[v.idx()];
    s.resize(num_loops);
    h.resize(num_loops);
    for (size_t loop = 0; loop < num_loops; ++loop) {
      s[loop].resize(num_sides[loop]);
      h[loop].resize(num_sides[loop]);
    }
  }

  auto& mesh = meshDomain;
  const size_t num_vert = mesh.n_vertices();

  SparseMatrix QQ(num_vert, num_vert);

  buildMatrixVertexLaplace(mesh, QQ);

  const size_t num_sides_all = std::accumulate(num_sides.begin(), num_sides.end(), 0);

  //std::cout << "Computing h-coordinates..." << std::endl;
  { // Compute h-coordinates
    std::vector<VertexHandle> sides_pts;
    for (size_t li = 0; li < num_loops; ++li) {
      for (size_t si = 0; si < num_sides[li]; ++si) {
        for (size_t i = 0; i < domain_boundary_vertices[li][si].size() - 1; ++i) {
          auto vtx = domain_boundary_vertices[li][si][i];
          sides_pts.push_back(vtx);
        }
      }
    }

    const size_t num_cons = sides_pts.size();
    SparseMatrix Ch(num_cons, num_vert), KKT;
    DenseMatrix pp = DenseMatrix::Zero(num_vert, num_sides_all);
    DenseMatrix dh = DenseMatrix::Ones(num_cons, num_sides_all);
    DenseMatrix rhs, x;

    addConstraint2Matrix(mesh, sides_pts, Ch);

    size_t idx = 0;
    for (size_t loop = 0; loop < num_loops; ++loop) {
      const bool three_sided = (num_sides[loop] == 3);
      const bool is_periodic = periodic[loop];
      if(!is_periodic){
        for (size_t side = 0; side < num_sides[loop]; ++side) {
          const size_t side_m1 = prev(loop, side);
          const size_t side_p1 = next(loop, side);
          std::vector<VertexHandle> side_pts;
          std::vector<VertexHandle> side_pts_m1;
          std::vector<VertexHandle> side_pts_p1;
          // Constraint values collected ALONGSIDE the kept vertices:
          // each value is computed from its own vertex's parameter, so
          // windowed (non-prefix) subsets can never be paired with the
          // head of the full parameter array, and the ramps stay linear
          // in the parameter under non-uniform boundary sampling (the
          // index-synthesized ramps assumed the old uniform sampling).
          std::vector<double> side_val_m1;
          std::vector<double> side_val_p1;
          for (size_t i = 0; i < domain_boundary_vertices[loop][side].size(); ++i) {
            auto vtx = domain_boundary_vertices[loop][side][i];
            side_pts.push_back(vtx);
          }
          for (size_t i = (three_sided ? 1 : 0); i < domain_boundary_vertices[loop][side_m1].size(); ++i) {
            auto vtx = domain_boundary_vertices[loop][side_m1][i];
            const double t = domain_boundary_params[loop][side_m1][i];
            const size_t nseg = num_segments[loop][side_m1];
            // Window width in the parameter: the last segment under
            // restrict_params, the h-width when enabled, else the
            // whole side.
            const bool seg_window = restrict_params && nseg > 1;
            const double w = seg_window     ? 1.0 / double(nseg)
                             : use_h_widths ? h_widths[loop][side][0]
                                            : 1.0;
            const bool keep =
              seg_window ? (t >= double(nseg - 1) / nseg)
                         : (!use_h_widths ||
                            t > 1.0 - h_widths[loop][side][0]);
            if (keep) {
              side_pts_m1.push_back(vtx);
              // h rises 0 -> 1 from the shared corner (t = 1) across
              // the window; the clamp meets the unconstrained default
              // value 1 continuously at the window edge.
              side_val_m1.push_back(
                std::min(1.0, (1.0 - t) / std::max(w, 1e-12)));
            }
          }
          for (size_t i = 0; i < domain_boundary_vertices[loop][side_p1].size() - (three_sided ? 1 : 0); ++i) {
            auto vtx = domain_boundary_vertices[loop][side_p1][i];
            const double t = domain_boundary_params[loop][side_p1][i];
            const size_t nseg = num_segments[loop][side_p1];
            const bool seg_window = restrict_params && nseg > 1;
            const double w = seg_window     ? 1.0 / double(nseg)
                             : use_h_widths ? h_widths[loop][side][1]
                                            : 1.0;
            const bool keep =
              seg_window ? (t <= 1.0 / double(nseg))
                         : (!use_h_widths ||
                            t <= h_widths[loop][side][1]);
            if (keep) {
              side_pts_p1.push_back(vtx);
              side_val_p1.push_back(
                std::min(1.0, t / std::max(w, 1e-12)));
            }
          }

          addConstraint2RHS(
            mesh,
            side_pts,
            dh,
            idx,
            false,
            0.0,
            0.0,
            0.0,
            true,
            true
          );

          // Values are final (window-normalized, subset-aligned), so
          // both neighbors go through the values overload unreversed;
          // the former per-branch ramp/values calls are unified.
          addConstraint2RHS(
            mesh,
            side_pts_m1,
            side_val_m1,
            dh,
            idx,
            false,
            true,
            true
          );

          addConstraint2RHS(
            mesh,
            side_pts_p1,
            side_val_p1,
            dh,
            idx,
            false,
            true,
            true
          );

          // addConstraint2RHS(
          //   mesh,
          //   side_pts,
          //   dh,
          //   idx,
          //   false,
          //   0.0,
          //   0.0,
          //   0.0,
          //   true,
          //   true
          // );

          // addConstraint2RHS(
          //   mesh,
          //   side_pts_m1,
          //   dh,
          //   idx,
          //   true,
          //   0.0,
          //   1.0,
          //   0.0,
          //   true,
          //   true
          // );

          // addConstraint2RHS(
          //   mesh,
          //   side_pts_p1,
          //   dh,
          //   idx,
          //   true,
          //   0.0,
          //   0.0,
          //   1.0,
          //   true,
          //   true
          // );

          ++idx;
        }
      }
      else {
        //Periodic loop
        for(size_t si = 0; si < num_sides[loop]; ++si) {
          for (size_t side = 0; side < num_sides[loop]; ++side) {
            std::vector<VertexHandle> side_pts;
            for (size_t i = 0; i < domain_boundary_vertices[loop][side].size() - 1; ++i) {
              auto vtx = domain_boundary_vertices[loop][side][i];
              side_pts.push_back(vtx);
            }

            addConstraint2RHS(
              mesh,
              side_pts,
              dh,
              idx,
              false,
              0.0,
              0.0,
              0.0,
              true,
              true
            );
          }
          ++idx;
        }
      }
    }

    buildMatrixKKTSystem(QQ, Ch, KKT);
    KKT.makeCompressed();
    Ch.resize(0, 0);
    //CC.data().squeeze();
    buildMatrixKKTRHS(pp, dh, rhs);
    dh.resize(0, 0);

    if (!solveLinearSystem(KKT, rhs, x)) {
      return false;
    }

    for (auto v : mesh.vertices()) {
      size_t idx = 0;
      for (size_t loop = 0; loop < num_loops; ++loop) {
        for (size_t side = 0; side < num_sides[loop]; ++side) {
          h_coords[v.idx()][loop][side] = x(v.idx(), idx);
          ++idx;
        }
      }
    }

  }

  //Compute s-coordinates
  //std::cout << "Computing s-coordinates..." << std::endl;
  for (size_t loop = 0; loop < num_loops; ++loop) {
    const bool three_sided = (num_sides[loop] == 3);
    const bool is_periodic = periodic[loop];
    if(!is_periodic){
      for (size_t side = 0; side < num_sides[loop]; ++side) {
        // std::cout << "Computing s-coordinates for loop " << loop << " side " << side << std::endl;
        const size_t side_m1 = prev(loop, side);
        const size_t side_p1 = next(loop, side);
        std::vector<VertexHandle> side_pts;
        std::vector<VertexHandle> side_pts_m1;
        std::vector<VertexHandle> side_pts_p1;
        for (size_t i = 0; i < domain_boundary_vertices[loop][side].size() - 1; ++i) {
          auto vtx = domain_boundary_vertices[loop][side][i];
          side_pts.push_back(vtx);
        }
        for (size_t i = (three_sided ? 1 : 0); i < domain_boundary_vertices[loop][side_m1].size()- 1; ++i) {
          auto vtx = domain_boundary_vertices[loop][side_m1][i];
          size_t nseg = num_segments[loop][side_m1];
          if (nseg == 1 || !restrict_params) {
            if(!use_h_widths || domain_boundary_params[loop][side_m1][i] > 1.0 - h_widths[loop][side][0]) {
              side_pts_m1.push_back(vtx);
            }
          }
          else {
            if (domain_boundary_params[loop][side_m1][i] >= double(nseg - 1) / nseg) {
              side_pts_m1.push_back(vtx);
            }
          }
        }
        for (size_t i = 0; i < domain_boundary_vertices[loop][side_p1].size() - (three_sided ? 1 : 0); ++i) {
          auto vtx = domain_boundary_vertices[loop][side_p1][i];
          size_t nseg = num_segments[loop][side_p1];
          if (nseg == 1 || !restrict_params) {
            if(!use_h_widths || domain_boundary_params[loop][side_p1][i] <= h_widths[loop][side][1]) {
              side_pts_p1.push_back(vtx);
            }
          }
          else {
            if (domain_boundary_params[loop][side_p1][i] <= double(1) / nseg) {
              side_pts_p1.push_back(vtx);
            }
          }
        }

        std::vector<VertexHandle> sides_pts;
        concatenateVectors(
          side_pts_m1,
          side_pts,
          sides_pts,
          false
        );
        concatenateVectors(
          std::vector<VertexHandle>(sides_pts),
          side_pts_p1,
          sides_pts,
          false
        );

        // std::cout << "sides_pts size: " << sides_pts.size() << std::endl;

        const size_t num_cons = sides_pts.size();
        SparseMatrix CC(num_cons, num_vert), KKT;
        DenseMatrix pp = DenseMatrix::Zero(num_vert, 1);
        DenseMatrix dd;
        DenseMatrix rhs, x;

        addConstraint2Matrix(mesh, sides_pts, CC);

        // std::cout << "CC size: " << CC.rows() << " " << CC.cols() << std::endl;
        // //Printing non-zeroes of CC
        // for (int k = 0; k < CC.outerSize(); ++k) {
        //   for (SparseMatrix::InnerIterator it(CC, k); it; ++it) {
        //     std::cout << "CC[" << it.row() << "][" << it.col() << "] = " << it.value() << std::endl;
        //   }
        // }

        addConstraint2RHS(
          mesh,
          side_pts_m1,
          dd,
          0,
          false,
          0.0,
          0.0,
          0.0,
          false,
          false
        );

        addConstraint2RHS(
          mesh,
          side_pts,
          domain_boundary_params[loop][side],
          dd,
          0,
          false,
          false,
          false
        );

        addConstraint2RHS(
          mesh,
          side_pts_p1,
          dd,
          0,
          false,
          1.0,
          0.0,
          1.0,
          false,
          true
        );

        // addConstraint2RHS(
        //   mesh,
        //   side_pts_m1,
        //   dd,
        //   0,
        //   false,
        //   0.0,
        //   0.0,
        //   0.0,
        //   false,
        //   false
        // );

        // addConstraint2RHS(
        //   mesh,
        //   side_pts,
        //   dd,
        //   0,
        //   true,
        //   0.0,
        //   0.0,
        //   1.0,
        //   false,
        //   false
        // );

        // addConstraint2RHS(
        //   mesh,
        //   side_pts_p1,
        //   dd,
        //   0,
        //   false,
        //   1.0,
        //   0.0,
        //   1.0,
        //   false,
        //   true
        // );

        // //printing elements of dd
        // for (int k = 0; k < dd.rows(); ++k) {
        //   std::cout << "dd[" << k << "] = " << dd(k, 0) << std::endl;
        //   // draw dashes after every 100th element
        //   if (k % 100 == 99) {
        //     std::cout << "------------------------" << std::endl;
        //   }
        // }

        buildMatrixKKTSystem(QQ, CC, KKT);
        KKT.makeCompressed();
        CC.resize(0, 0);
        //CC.data().squeeze();
        buildMatrixKKTRHS(pp, dd, rhs);
        dd.resize(0, 0);

        if (!solveLinearSystem(KKT, rhs, x)) {
          return false;
        }

        for (auto v : mesh.vertices()) {
          s_coords[v.idx()][loop][side] = x(v.idx(), 0);
        }

      }
    }
    else {
      //Periodic loop
      std::vector<VertexHandle> side_pts;
      for(size_t si = 0; si < num_sides[loop]; ++si) {
        for (size_t i = 0; i < domain_boundary_vertices[loop][si].size() - 1; ++i) {
          auto vtx = domain_boundary_vertices[loop][si][i];
          side_pts.push_back(vtx);
        }
      }
      const size_t num_cons = side_pts.size();
      SparseMatrix CC(num_cons, num_vert), KKT;
      DenseMatrix pp = DenseMatrix::Zero(num_vert, 1);
      DenseMatrix dd;
      DenseMatrix rhs, x;

      addConstraint2Matrix(mesh, side_pts, CC);

      for(size_t si = 0; si < num_sides[loop]; ++si) {
        std::vector<VertexHandle> pts;
        for (size_t i = 0; i < domain_boundary_vertices[loop][si].size() - 1; ++i) {
          auto vtx = domain_boundary_vertices[loop][si][i];
          pts.push_back(vtx);
        }
        

        auto pars = domain_boundary_params[loop][si];
        const auto u_min = double(si) / double(num_sides[loop]);
        const auto u_max = double(si + 1) / double(num_sides[loop]);

        for(size_t i = 0; i < pars.size(); ++i) {
          pars[i] = u_min + pars[i] * (u_max - u_min);
        }

        addConstraint2RHS(
          mesh,
          pts,
          pars,
          dd,
          0,
          false,
          false,
          false
        );
      }

      // Adding cut constraints to fix the periodicity
      std::vector<double> cut_sign(mesh.n_halfedges(), 0);
      std::vector<Mesh::Point> grad_h;
      mesh.computeFaceGradientofFunction(
        [&](Mesh::VertexHandle vh) {
          return h_coords[vh.idx()][loop][0];
        },
        grad_h
      );

      std::vector<Mesh::Point>         cut_pts;
      std::vector<Mesh::HalfedgeHandle> cut_hes;
      VertexHandle vv = side_pts.front();
      if (!mesh.traceVectorFieldonFacesfromVertex(
        grad_h,
        vv,
        cut_pts,
        cut_hes
      )) {
        // replacing every vector in grad_h with a constant vector perpendicular to the edge at vv and trying again
        auto edge_dir = mesh.calc_edge_vector(mesh.halfedge_handle(vv)).normalized();
        edge_dir = (edge_dir.cross(Mesh::Point(0.0, 0.0, 1.0)).normalized());
        for (auto& g : grad_h) {
          g = edge_dir;
        }
        if(!mesh.traceVectorFieldonFacesfromVertex(
          grad_h,
          vv,
          cut_pts,
          cut_hes
        )) {
          return false;
        }
      }
      

      if(debug_outputs) {
        // Writing cut curve for debugging
        std::vector<std::vector<std::vector<Eigen::Vector3d> > > cut_curve_(1);
        cut_curve_.front().resize(1);
        cut_curve_.front().back().reserve(cut_pts.size());
        for(size_t i = 0; i < cut_pts.size(); ++i) {
          Eigen::Vector3d p;
          p << cut_pts[i][0], cut_pts[i][1], cut_pts[i][2];
          cut_curve_.front().back()[i] = p;
        }
        writeLoops(cut_curve_, std::string("cut_curve_" + std::to_string(loop) + ".obj"));
      }


      for (auto he : mesh.halfedges()) {
        cut_sign[he.idx()] = 0;
      }
      for (auto he : cut_hes) {
        cut_sign[he.idx()] = -1;
        cut_sign[mesh.opposite_halfedge_handle(he).idx()] = 1;
      }

      for (auto vv : mesh.vertices()) {
        for (auto he : mesh.voh_range(vv)) {
          int sign = cut_sign[he.idx()];
          int v2_idx = mesh.to_vertex_handle(he).idx();
          double ww = -QQ.coeffRef(vv.idx(), v2_idx);
          if (sign != 0) {
            pp(vv.idx(), 0) += sign * ww*1.0;
          }
        }
      }


      buildMatrixKKTSystem(QQ, CC, KKT);
      KKT.makeCompressed();
      CC.resize(0, 0);
      //CC.data().squeeze();
      buildMatrixKKTRHS(pp, dd, rhs);
      dd.resize(0, 0);

      if (!solveLinearSystem(KKT, rhs, x)) {
        return false;
      }

      for (auto v : mesh.vertices()) {
        double s = x(v.idx(), 0);
        s = s < 0.0 ? s + 1.0 : s;
        s = s > 1.0 ? s - 1.0 : s;
        for(size_t side = 0; side < num_sides[loop]; ++side) {
          s_coords[v.idx()][loop][side] = s;
        }
      }
      
    }
  }


  return true;
}

bool SurfGBS::compute_deformed_parameters()
{
  auto& mesh = meshDomain;

  // Detecting sides to deform based on turning angles (integrated curvature)
  for (size_t li = 0; li < num_loops; ++li) {
    for (size_t si = 0; si < num_sides[li]; ++si) {
      // const auto & cc = domain_boundary_curves[li][si];
      // auto t0 = (cc[1] - cc[0]).normalized();
      // auto t1 = (cc[cc.size() - 1] - cc[cc.size() - 2]).normalized();
      // double angle = std::acos(std::min(std::max(t0.dot(t1), -1.0), 1.0));
      // std::cout << "Turning angle for side " << si << " of loop " << li << " is " << angle * 180.0 / M_PI << " degrees." << std::endl;
      // if (angle > M_PI / 4) {
      //   side_deformed[li][si] = true;
      //   deform_splines[li][si].controlPoint(1,1) = { 0.0, 0.0, 0.0 }; //pull down the middle control point to create a "valley"
      //   std::cout << "Deforming side " << si << " of loop " << li << " with turning angle " << angle * 180.0 / M_PI << " degrees." << std::endl;
      // } else {
      //   side_deformed[li][si] = false;
      // }

      // For now, just deform ribbons with more than deg + 1 control points in u direction
      if (ribbons[li][si].numControlPoints()[0] > ribbons[li][si].basisU().degree() + 1) {
        side_deformed[li][si] = true;
        deform_splines[li][si].controlPoint(1,1) = { deform_value, deform_value, deform_value }; //pull down the middle coefficient
        std::cout << "Deforming side " << si << " of loop " << li << std::endl;
      } 
      else {
        side_deformed[li][si] = false;
      }

    }
  }

  h_coords_deformed = h_coords;

  // h_coords_deformed.clear();
  // h_coords_deformed.resize(meshDomain.n_vertices());
  // for (const auto v : meshDomain.vertices()) {
  //   auto& h = h_coords_deformed[v.idx()];
  //   h.resize(num_loops);
  //   for (size_t loop = 0; loop < num_loops; ++loop) {
  //     h[loop].resize(num_sides[loop]);
  //   }
  // }


  for (auto v : mesh.vertices()) {
    for (size_t li = 0; li < num_loops; ++li) {
      for (size_t si = 0; si < num_sides[li]; ++si) {
        if (side_deformed[li][si]) {
          const double s = s_coords[v.idx()][li][si];
          const double h = h_coords[v.idx()][li][si];
          const double h_def = deform_splines[li][si].eval(s, h)[0];
          h_coords_deformed[v.idx()][li][si] = h_def;
        }
      }
    }
  }

  return true;
}
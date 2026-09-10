#include <fstream>

#include "libcdgbs/SurfGBS.hpp"
#include "libcdgbs/LoopFlattener.hpp"
#include "libcdgbs/TriangleWrapper.hpp"
#include "libcdgbs/GeomUtil.hpp"

using namespace libcdgbs;
using Vec3 = Eigen::Vector3d;
using Vec2 = Eigen::Vector2d;
using SubCurve3D = std::vector<Vec3>;
using Curve3D = std::vector<SubCurve3D>;
using Curves3D = std::vector<Curve3D>;

//calculate arclength of each ribbon boundary B-spline curve
double getLength(const Geometry::BSSurface& rib, double u_start = 0.0, double u_end = 1.0, size_t num_samples = 200) {
  double length = 0.0;
  for (size_t i = 0; i < num_samples; ++i) {
    double ua = (1.0 - (double(i) / double(num_samples))) * u_start + (double(i) / double(num_samples)) * u_end;
    double ub = (1.0 - (double(i + 1) / double(num_samples))) * u_start + (double(i + 1) / double(num_samples)) * u_end;;
    auto pa = rib.eval(ua, 0.0);
    auto pb = rib.eval(ub, 0.0);
    length += (pb - pa).norm();
  }
  return length;
}

double trapezoid(const std::function<double(double)> &f, double a, double b,
                 size_t n, double s = 0) {
  if (n == 1)
    return (f(a) + f(b)) / 2 * (b - a);
  size_t k = 1 << (n - 2);
  double h = (b - a) / k;
  double x = a + h / 2;
  double sum = 0;
  for (size_t i = 0; i < k; ++i)
    sum += f(x + h * i);
  return (s + h * sum) / 2;
}

double simpson(const std::function<double(double)> &f, double a, double b,
               size_t max_iterations, double tolerance) {
  double s = 0;
  double s_prev = -std::numeric_limits<double>::max(), st_prev = s_prev;
  for (size_t i = 1; i <= max_iterations; ++i) {
    double st = trapezoid(f, a, b, i, st_prev);
    s = (st * 4 - st_prev) / 3;
    if (i > 5 && // avoid exiting too early
        (std::abs(s - s_prev) < tolerance * std::abs(s_prev) ||
         (s == 0 && s_prev == 0)))
      break;
    s_prev = s;
    st_prev = st;
  }
  return s;
}

// Makes at most 2^(iterations-1)+1 evaluations (with 1 derivative),
// but exits early if already achieved the specified relative tolerance
double arclength(const Geometry::BSSurface &curve, double from, double to,
                 size_t iterations, double tolerance) {
  auto result =
    simpson([&](double u) {
      Geometry::VectorMatrix der;
      curve.eval(u, 0.0, 1, der);
      return der[1][0].norm();
    }, from, to, iterations, tolerance);
  return result;
}


SurfGBS::SurfGBS()
{

}

void SurfGBS::load_ribbons(const std::vector<std::vector<Ribbon> >& ribbon_surfs, const InputParams& params)
{
  ribbons = ribbon_surfs;

  SurfGBS::target_length = params.target_length;
  SurfGBS::deform_value = params.deform_value;
  SurfGBS::restrict_params = params.restrict_params;
  SurfGBS::c1_merge = params.c1_merge;
  SurfGBS::global_inner_loop_scale = params.global_inner_loop_scale;
  SurfGBS::use_h_widths = params.use_h_widths;
  SurfGBS::arclength_sampling = params.arclength_sampling;
  SurfGBS::merge_inner_c1 = params.merge_inner_c1;
  SurfGBS::ribbon_init_placement = params.ribbon_init_placement;

  num_segments.clear();
  num_segments.resize(ribbons.size());
  for (size_t loop = 0; loop < ribbons.size(); ++loop) {
    num_segments[loop].resize(ribbons[loop].size(), 1);
  }
  init_data();

  if (params.merge_smooth_corners) {
    merge_smooth_corners();
  }
}

bool SurfGBS::writeLoops(const std::vector<std::vector<std::vector<Eigen::Vector3d> > > &loops,
  const std::string filename, bool edges) {
  size_t index = 1;
  std::ofstream f(filename);
  if(!f.is_open()) {
    return false;
  }
  for (const auto& loop : loops) {
    for (const auto& curve : loop) {
      size_t start = index;
      for (const auto& p : curve) {
        f << "v " << p[0] << ' ' << p[1] << ' ' << p[2] << std::endl;
        index++;
      }
      if(curve.size() > 1 && edges) {
        f << 'l';
        for (size_t i = start; i < index; ++i)
          f << ' ' << i;
        f << std::endl;
      }
    }
  }
  return true;
}

void SurfGBS::load_ribbons_and_evaluate(const std::vector<std::vector<Ribbon> >& ribbon_surfs, Mesh& mesh, const InputParams& params, SurfaceType surface_type)
{
  load_ribbons(ribbon_surfs, params);
  compute_domain_boundary();
  // writeLoops(domain_boundary_curves, "/tmp/boundary.obj");
  compute_domain_mesh();
  compute_local_parameters();
  compute_blend_functions();
  evaluate_mesh(mesh, true, surface_type);
}


void SurfGBS::init_data()
{
  num_sides.clear();
  num_rows.clear();
  num_cols.clear();
  side_res.clear();
  side_segment_res.clear();
  deg_h.clear();
  h_widths.clear();
  deform_splines.clear();
  side_deformed.clear();
  periodic.clear();

  num_loops = ribbons.size();
  num_sides.resize(num_loops);
  num_rows.resize(num_loops);
  num_cols.resize(num_loops);
  side_res.resize(num_loops);
  side_segment_res.resize(num_loops);
  deg_h.resize(num_loops);
  h_widths.resize(num_loops);
  deform_splines.resize(num_loops);
  side_deformed.resize(num_loops);
  periodic.resize(num_loops, false);
  for (size_t loop = 0; loop < num_loops; ++loop) {
    const size_t ns = ribbons[loop].size();
    num_sides[loop] = ns;
    deg_h[loop].resize(ns, 3);
    h_widths[loop].resize(ns, std::array<double, 2>({1.0, 1.0}));
    num_rows[loop].resize(ns);
    num_cols[loop].resize(ns);
    side_res[loop].resize(ns);
    side_segment_res[loop].resize(ns);
    for (size_t side = 0; side < ns; ++side) {
      //Nmber of B-spline control points in U and V directions
      num_rows[loop][side] = ribbons[loop][side].basisV().knots().size() - ribbons[loop][side].basisV().degree() - 1;
      num_cols[loop][side] = ribbons[loop][side].basisU().knots().size() - ribbons[loop][side].basisU().degree() - 1;
      side_res[loop][side] = 100; // default resolution
      side_segment_res[loop][side].resize(num_segments[loop][side], 100); // default resolution per segment
    }
    if(merge_inner_c1 && loop > 0) {
      periodic[loop] = true;
    }
  }

  // Set number of samples based on target edge length
  for (size_t loop = 0; loop < num_loops; ++loop) {
    for (size_t side = 0; side < num_sides[loop]; ++side) {
      auto rib = ribbons[loop][side];
      if(num_segments[loop][side] == 1) {
        auto curve_length = getLength(rib);

        auto num_samples = std::max(1.0, curve_length / target_length);
        side_res[loop][side] = std::max(static_cast<size_t>(num_samples), size_t(5));
        //side_res[loop][side] = 200; //forcing 200 for now
        side_segment_res[loop][side][0] = side_res[loop][side];
      }
      else {
        side_res[loop][side] = 0;
        for(size_t seg = 0; seg < num_segments[loop][side]; ++seg) {
          double u_start = double(seg) / (num_segments[loop][side]);
          double u_end = double(seg + 1) / (num_segments[loop][side]);
          auto curve_length = getLength(rib, u_start, u_end);

          auto num_samples = std::max(1.0, curve_length / target_length);
          const size_t seg_res = std::max(static_cast<size_t>(num_samples), size_t(5));
          side_res[loop][side] += seg_res;
          side_segment_res[loop][side][seg]= seg_res;
          if(seg > 0) {
            side_res[loop][side] -= 1; //remove duplicate point at segment boundary
          }
        }
      }
    }
  }

  // Initialize deformation splines (bi-quadratic B-splines with 3x3 control points)
  for (size_t loop = 0; loop < num_loops; ++loop) {
    for (size_t side = 0; side < num_sides[loop]; ++side) {
      const size_t deg_u = 2;
      const size_t deg_v = 2;
      const Geometry::DoubleVector knots {0.0, 0.0, 0.0, 1.0, 1.0, 1.0};

      Geometry::VectorVector cps;
      cps.push_back({ 0.0, 0.0, 0.0 });
      cps.push_back({ 0.5, 0.5, 0.5 });
      cps.push_back({ 1.0, 1.0, 1.0 });
      cps.push_back({ 0.0, 0.0, 0.0 });
      cps.push_back({ 0.5, 0.5, 0.5 });
      cps.push_back({ 1.0, 1.0, 1.0 });
      cps.push_back({ 0.0, 0.0, 0.0 });
      cps.push_back({ 0.5, 0.5, 0.5 });
      cps.push_back({ 1.0, 1.0, 1.0 });
      
      deform_splines[loop].push_back(Geometry::BSSurface(deg_u, deg_v, knots, knots, cps));
      side_deformed[loop].push_back(false);
    }
  }


}

bool SurfGBS::compute_domain_boundary()
{
  developed_boundary_curves.clear();
  developed_boundary_curves_normalized.clear();
  domain_boundary_curves.clear();
  developed_boundary_curves.resize(num_loops);
  developed_boundary_curves_normalized.resize(num_loops);
  domain_boundary_curves.resize(num_loops);
  domain_boundary_params.clear();
  domain_boundary_params.resize(num_loops);

  boundary_points.clear();
  boundary_normals.clear();
  boundary_tangents.clear();
  boundary_tangents_unnormalized.clear();
  boundary_crossderivatives.clear();

  boundary_points.resize(num_loops);
  boundary_normals.resize(num_loops);
  boundary_tangents.resize(num_loops);
  boundary_tangents_unnormalized.resize(num_loops);
  boundary_crossderivatives.resize(num_loops);

  //const size_t arclength_res = 200;

  //std::cout << "Arclength sampling: " << (arclength_sampling ? "enabled" : "disabled") << std::endl;

  for (size_t loop = 0; loop < num_loops; ++loop) {
    boundary_points[loop].resize(num_sides[loop]);
    boundary_normals[loop].resize(num_sides[loop]);
    boundary_tangents[loop].resize(num_sides[loop]);
    boundary_tangents_unnormalized[loop].resize(num_sides[loop]);
    boundary_crossderivatives[loop].resize(num_sides[loop]);
    domain_boundary_params[loop].resize(num_sides[loop]);
    for (size_t side = 0; side < num_sides[loop]; ++side) {
      const auto& rib = ribbons[loop][side];
      auto rib_rev = rib; rib_rev.reverseU();
      Geometry::PointVector cpts_fwd, cpts_rev;
      const auto nu = rib.numControlPoints()[0];
      for (size_t i = 0; i < nu; ++i) {
        cpts_fwd.push_back(rib.controlPoint(i, 0));
        cpts_rev.push_back(rib_rev.controlPoint(i, 0));
      }
      Geometry::BSCurve curve_fwd(rib.basisU().degree(), rib.basisU().knots(), cpts_fwd);
      Geometry::BSCurve curve_rev(rib_rev.basisU().degree(), rib_rev.basisU().knots(), cpts_rev);

      // Determine fixed orientation by lexicographic order of endpoints
      const auto p_start = rib.eval(0.0, 0.0);
      const auto p_end = rib.eval(1.0, 0.0);
      const bool reversed = (std::tie(p_start[0], p_start[1], p_start[2]) < std::tie(p_end[0], p_end[1], p_end[2])); 
      //(rib.eval(0.0, 0.0).normSqr() > rib.eval(1.0, 0.0).normSqr());

      const size_t res = side_res[loop][side];
      boundary_points[loop][side].reserve(res);
      boundary_normals[loop][side].reserve(res);
      boundary_tangents[loop][side].reserve(res);
      boundary_tangents_unnormalized[loop][side].reserve(res);
      boundary_crossderivatives[loop][side].reserve(res);
      domain_boundary_params[loop][side].reserve(res);

      const auto ns = num_segments[loop][side];
      if (ns > 1) {
        for (size_t i = 0; i < ns; ++i) {
          const double u_start = double(i) / double(ns);
          const double u_end = double(i + 1) / double(ns);
          const double u_start_rev = double(ns - 1 - i) / double(ns);
          const double u_end_rev = double(ns - 2 - i) / double(ns);
          const size_t seg_res = side_segment_res[loop][side][i];
          const double seg_len = curve_fwd.arcLength(u_start, u_end);
          const auto p_seg_start = rib.eval(u_start, 0.0);
          const auto p_seg_end = rib.eval(u_end, 0.0);
          const bool seg_reversed = (std::tie(p_seg_start[0], p_seg_start[1], p_seg_start[2]) < std::tie(p_seg_end[0], p_seg_end[1], p_seg_end[2])); 
          //bool seg_reversed = (rib.eval(u_start, 0.0).normSqr() > rib.eval(u_end, 0.0).normSqr());
          for(size_t j = 0; j < seg_res; ++j) {
            if(i > 0 && j == 0) {
              continue; //skip duplicate point at segment boundary
            }
            double uj = (double(j) / double(seg_res - 1));
            if(arclength_sampling && j > 0 && j < seg_res - 1) { // Find point with arclength u * len by bisection
              double target_len = uj * seg_len;
              double low_u = 0.0;
              double high_u = 1.0;
              double mid_u = 0.5;
              const size_t max_iters = 20;
              for(size_t iter = 0; iter < max_iters; ++iter) {
                mid_u = (low_u + high_u) * 0.5;
                double mid_uj = (1.0 - mid_u) * u_start + mid_u * u_end;
                double mid_len = seg_reversed ? (seg_len - curve_rev.arcLength(1.0 - u_end, 1.0 - mid_uj)) :
                  curve_fwd.arcLength(u_start, mid_uj);
                if(mid_len < target_len) {
                  low_u = mid_u;
                }
                else {
                  high_u = mid_u;
                }
              }
              uj = (low_u + high_u) * 0.5;
            }
            double u = (1.0 - uj) * u_start + uj * u_end;            
            Geometry::VectorMatrix duv;
            auto pt = rib.eval(u, 0.0, 1, duv);

            boundary_points[loop][side].push_back({ pt[0], pt[1], pt[2] });
            auto du = duv[1][0];
            auto dv = duv[0][1];
            auto nn = (du ^ dv).normalized();
            auto tt = du.normalized();
            boundary_tangents[loop][side].push_back({ tt[0], tt[1], tt[2] });
            boundary_tangents_unnormalized[loop][side].push_back({ du[0], du[1], du[2] });
            boundary_normals[loop][side].push_back({ nn[0], nn[1], nn[2] });
            boundary_crossderivatives[loop][side].push_back({ dv[0], dv[1], dv[2] });
            domain_boundary_params[loop][side].push_back(u);
          }
        }
      }
      else {
        const double len = curve_fwd.arcLength(0.0, 1.0);
        // std::vector<double> params_forward;
        // std::vector<double> params_reverse;
        // std::vector<Vec3> points_forward;
        // std::vector<Vec3> points_reverse;
        //std::cout << "Ribbon Loop " << loop << " Side " << side << " Length difference:" << std::abs(len - len_rev) << std::endl;
        for (size_t i = 0; i < res; ++i) {
          double u = double(i) / double(res - 1);

          if(arclength_sampling && i > 0 && i < res - 1) { // Find point with arclength u * len by bisection
            double target_len = u * len;
            double low_u = 0.0;
            double high_u = 1.0;
            double mid_u = 0.5;
            const size_t max_iters = 20;
            for(size_t iter = 0; iter < max_iters; ++iter) {
              mid_u = (low_u + high_u) * 0.5;
              double mid_len = reversed ? (len - curve_rev.arcLength(0.0, 1.0 - mid_u)) :
                curve_fwd.arcLength(0.0, mid_u);
              if(mid_len < target_len) {
                low_u = mid_u;
              }
              else {
                high_u = mid_u;
              }
            }
            u = (low_u + high_u) * 0.5;
          }
          // params_forward.push_back(u);
          // auto pt_fw = rib.eval(u, 0.0);
          // points_forward.push_back({ pt_fw[0], pt_fw[1], pt_fw[2] });

          // u = double(i) / double(res - 1);
          // if(arclength_sampling && i > 0 && i < res - 1) { // Find point with arclength u * len by bisection
          //   double target_len = u * len_rev;
          //   double low_u = 0.0;
          //   double high_u = 1.0;
          //   double mid_u = 0.5;
          //   const size_t max_iters = 20;
          //   for(size_t iter = 0; iter < max_iters; ++iter) {
          //     mid_u = (low_u + high_u) * 0.5;
          //     double mid_len = true ? (len - getLength(rib, 0.0, 1.0 - mid_u, arclength_res)) :
          //       getLength(rib_rev, 0.0, mid_u, arclength_res);
          //     if(mid_len < target_len) {
          //       low_u = mid_u;
          //     }
          //     else {
          //       high_u = mid_u;
          //     }
          //   }
          //   u = (low_u + high_u) * 0.5;
          // }
          // params_reverse.push_back(1.0 - u);
          // auto pt_rev = rib.eval(1.0 - u, 0.0);
          // points_reverse.push_back({ pt_rev[0], pt_rev[1], pt_rev[2] });

          Geometry::VectorMatrix duv;
          auto pt = rib.eval(u, 0.0, 1, duv);

          boundary_points[loop][side].push_back({ pt[0], pt[1], pt[2] });
          auto du = duv[1][0];
          auto dv = duv[0][1];
          auto nn = (du ^ dv).normalized();
          auto tt = du.normalized();
          boundary_tangents[loop][side].push_back({ tt[0], tt[1], tt[2] });
          boundary_tangents_unnormalized[loop][side].push_back({ du[0], du[1], du[2] });
          boundary_normals[loop][side].push_back({ nn[0], nn[1], nn[2] });
          boundary_crossderivatives[loop][side].push_back({ dv[0], dv[1], dv[2] });
          domain_boundary_params[loop][side].push_back(u);
        }
        // //reverse params_reverse for comparison
        // std::reverse(params_reverse.begin(), params_reverse.end());
        // std::reverse(points_reverse.begin(), points_reverse.end());
        // for(size_t i = 0; i < params_forward.size(); ++i) {
        //   std::cout << "Parameterization mismatch at Loop " << loop << " Side " << side << " Index " << i 
        //     << ": " << std::abs(params_forward[i] - params_reverse[i]) << std::endl;
        //   std::cout << "Position mismatch at Loop " << loop << " Side " << side << " Index " << i 
        //     << ": " << (points_forward[i] - points_reverse[i]).norm() << std::endl;
        // }
      }
    }

    if (debug_outputs) {
      writeLoops(boundary_points, "boundary_3D.obj");
    }

    // Flip vector ?
    auto va = GeomUtil::vectorArea(boundary_points[loop]);
    Eigen::Vector3d avg_n(0.0, 0.0, 0.0);
    size_t num_pts = 0;
    for (size_t side = 0; side < boundary_normals[loop].size(); ++side) {
      for (size_t ii = 0; ii < boundary_normals[loop][side].size(); ++ii) {
        avg_n += boundary_normals[loop][side][ii];
        num_pts++;
      }
    }
    avg_n /= num_pts;
    bool flip = loop > 0 && (va.dot(avg_n)) < 0;
    if(flip) {
      std::cout << "Flipping loop " << loop << " with vector area " << va.transpose() << " and average normal " << avg_n.transpose() << std::endl;
    }

    developed_boundary_curves[loop] =
      LoopFlattener::developLoop(boundary_points[loop], boundary_normals[loop], flip);

    developed_boundary_curves_normalized[loop] =
      LoopFlattener::normalizeLoopAngles(boundary_points[loop], boundary_normals[loop], flip);

    domain_boundary_curves[loop] =
      LoopFlattener::closeLoop(developed_boundary_curves_normalized[loop]);

    if(loop > 0) {
      projectCurves2LSPlane(boundary_points[loop], domain_boundary_curves[loop]);
    }

    auto cog = GeomUtil::centroid(domain_boundary_curves[loop]);
    GeomUtil::shiftByVector(domain_boundary_curves[loop], -cog);
    // GeomUtil::shiftByVector(developed_boundary_curves_normalized[loop], -cog);
    // GeomUtil::shiftByVector(developed_boundary_curves[loop], -cog);

    auto aligned = domain_boundary_curves[loop];
    GeomUtil::alignCurveLoopPCA(domain_boundary_curves[loop], aligned, false);
    domain_boundary_curves[loop] = aligned;



  }

  if (debug_outputs) {
    writeLoops(domain_boundary_curves, "boundary_uv_0.obj");
  }

  if (num_loops > 1) { // Handling inner loops
    //Computing perimeter surface
    std::vector<std::vector<Ribbon> > perimeter_ribbons(1, ribbons.front());
    scale_perimeter_ribbons(perimeter_ribbons);
    SurfGBS perimeter_gbs;
    perimeter_gbs.load_ribbons(perimeter_ribbons, {target_length, true, 0.0, true, false, 1.0, true, arclength_sampling, false, false});
    perimeter_gbs.num_segments[0] = num_segments[0];
    perimeter_gbs.side_segment_res[0] = side_segment_res[0];
    perimeter_gbs.compute_domain_boundary();
    perimeter_gbs.compute_domain_mesh();
    perimeter_gbs.compute_local_parameters();
    perimeter_gbs.compute_blend_functions();
    perimeter_gbs.evaluate_mesh(perimeter_gbs.meshSurface, true);

    if (debug_outputs) {
      writeOBJ(perimeter_gbs.meshSurface, "perimeter.obj");
    }

    //Projecting inner loops
    for (size_t loop = 1; loop < num_loops; ++loop) {
      auto curv_proj = domain_boundary_curves[loop];
      const auto& curv = boundary_points[loop];
      perimeter_gbs.projectCurves2Domain(curv, curv_proj);
      if (debug_outputs) {
        writeLoops({ curv_proj }, std::string("boundary_proj_") + std::to_string(loop) + ".obj");
      }
      auto cog = GeomUtil::centroid(curv_proj);
      GeomUtil::shiftByVector(curv_proj, -cog);
      GeomUtil::alignPointSets(domain_boundary_curves[loop], curv_proj);
      // Scaling down loop based on distance of centroid from perimeter
      // auto cog_perimeter = GeomUtil::centroid(domain_boundary_curves[0]);
      // auto cog_loop = GeomUtil::centroid(domain_boundary_curves[loop]);
      // double dist_cog = (cog_perimeter - cog_loop).norm();
      // // double scale = 1.0 - calc_inner_loop_scale(domain_boundary_curves[loop], 
      // //   boundary_normals[loop], boundary_tangents[loop]);
      // //scale = 1.0; //std::max(scale, 0.3);
      // double scale = calc_inner_loop_scale(loop);
      // std::cout << "Loop " << loop << " scale factor: " << scale << std::endl;
      // //scale = 1.0; //disable for now
      // // Scaling loop
      // for (auto& sub : domain_boundary_curves[loop]) {
      //   for (auto& p : sub) {
      //     p *= scale;
      //   }
      // }
      GeomUtil::shiftByVector(domain_boundary_curves[loop], cog);
    }
    if(global_inner_loop_scale <= 0.0) {
      position_inner_loops();
    }
    else if(global_inner_loop_scale < 1.0) {
      // Scaling down inner loops
      for(size_t loop = 1; loop < num_loops; ++loop) {
        const double scale = global_inner_loop_scale;
        std::cout << "Scaling down loop " << loop << " by factor " << scale << std::endl;
        const auto cog = GeomUtil::centroid(domain_boundary_curves[loop]);
        for (auto& sub : domain_boundary_curves[loop]) {
          for (auto& p : sub) {
            p = cog + scale * (p - cog);
          }
        }
      }
    }
    else {
      //resolve_self_intersections();
    }

  }

  simple_domain.init(domain_boundary_curves);

  if (debug_outputs) {
    writeLoops(domain_boundary_curves, "boundary_uv_1.obj");
    writeLoops(developed_boundary_curves, "boundary_developed.obj");
    writeLoops(developed_boundary_curves_normalized, "boundary_normalized.obj");
    writeLoops(boundary_points, "boundary_xyz.obj");
  }

  compute_auto_h_widths(30.0, 2.0);

  return true;
}

bool SurfGBS::compute_domain_mesh()
{
  auto triangle_wrapper = TriangleWrapper();
  meshDomain = triangle_wrapper.triangulate_loop(domain_boundary_curves, target_length);

  domain_boundary_vertices.clear();
  domain_boundary_vertices.resize(num_loops);
  size_t idx_vtx = 0;
  for (size_t loop_idx = 0; loop_idx < domain_boundary_curves.size(); ++loop_idx) {
    domain_boundary_vertices[loop_idx].resize(num_sides[loop_idx]);
    auto const& loop = domain_boundary_curves[loop_idx];
    size_t first_idx = idx_vtx;
    for (size_t si = 0; si < loop.size(); ++si) {
      auto sub = loop[si];
      for (size_t i = 0; i < sub.size(); ++i) {
        auto vtx = meshDomain.vertex_handle((si == loop.size() - 1 && i == sub.size() - 1) ? first_idx : idx_vtx);
        domain_boundary_vertices[loop_idx][si].push_back(vtx);
        if (i < sub.size() - 1) {
          ++idx_vtx;
        }
      }
    }
  }

  return true;
}

bool SurfGBS::compute_local_parameters()
{
  if (true) {
    // Side-based quadratic parameterization
    s_coords.resize(meshDomain.n_vertices());
    h_coords.resize(meshDomain.n_vertices());
    for (const auto &v : meshDomain.vertices()) {
      auto q = meshDomain.point(v);
      Vec3 p(q[0], q[1], 0);
      auto &s = s_coords[v.idx()];
      auto &h = h_coords[v.idx()];
      s.resize(num_loops);
      h.resize(num_loops);
      simple_domain.computeParameters(p, h_widths[0], s, h);
    }
  } else if (!compute_harmonic_parameters()) {
    std::cout << "Error computing harmonic parameters" << std::endl;
    return false;
  }

  // auto s_tmp = s_coords; s_coords.clear();
  // h_coords.clear();
  // compute_harmonic_parameters();
  // s_coords = s_tmp;

  // auto h_tmp = h_coords; h_coords.clear();
  // s_coords.clear();
  // compute_harmonic_parameters();
  // h_coords = h_tmp;

  if(!compute_deformed_parameters()) {
    std::cout << "Error computing deformed parameters" << std::endl;
    return false;
  }

  return true;
}

bool SurfGBS::compute_blend_functions()
{
  blend_functions.clear();
  blend_functions.resize(meshDomain.n_vertices());
  for (const auto v : meshDomain.vertices()) {
    auto& Bf = blend_functions[v.idx()];
    Bf.resize(num_loops);
    for (size_t loop = 0; loop < num_loops; ++loop) {
      Bf[loop].resize(num_sides[loop]);
      const bool is_periodic = periodic[loop];
      for (size_t side = 0; side < num_sides[loop]; ++side) {
        Bf[loop][side].resize(
          num_rows[loop][side],
          std::vector<double>(num_cols[loop][side], 0.0)
        );

        const auto deg_s = ribbons[loop][side].basisU().degree();
        const auto deg_h = SurfGBS::deg_h[loop][side];

        const auto &hs = side_deformed[loop][side] ? h_coords_deformed : h_coords;

        const auto s = std::min(std::max(s_coords[v.idx()][loop][side], 0.0), 1.0);
        const auto h = std::min(std::max(hs[v.idx()][loop][side], 0.0), 1.0);

        const auto& Bu = ribbons[loop][side].basisU();
        const auto& Bv = Geometry::BSBasis(deg_h, [deg_h]() {
          Geometry::DoubleVector knots(deg_h + 1, 0.0);
          knots.insert(knots.end(), deg_h + 1, 1.0);
          return knots;
          }());

        const size_t span_u = Bu.findSpan(s), span_v = Bv.findSpan(h);
        Geometry::DoubleVector Bh, Bs;
        Bu.basisFunctions(span_u, s, Bs);
        Bv.basisFunctions(span_v, h, Bh);
        
        if(!is_periodic){
          for (size_t row = 0; row < num_rows[loop][side]; ++row) {
            for (size_t col = 0; col <= deg_s; ++col) {
              const size_t ri = row;
              const size_t ci = col + span_u - deg_s;
              const double mu = get_mu(v, loop, side, ri, ci);
              Bf[loop][side][ri][ci] =
                mu * Bs[col] * Bh[row];
            }
          }
        }
        else {
          auto BSp = periodicC1CubicBSplineBasis(s, num_sides[loop]);
          //Extracting only the relevant C1 periodic cubic basis functions for this side
          double B1 = BSp[side * 2];
          double B2 = BSp[side * 2 + 1];

          for (size_t row = 0; row < num_rows[loop][side]; ++row) {
            Bf[loop][side][row][1] = B1 * Bh[row];
            Bf[loop][side][row][2] = B2 * Bh[row];
          }
        }
      }
    }
  }

  return true;
}

bool SurfGBS::evaluate_mesh(Mesh& mesh, bool reset, SurfaceType surface_type)
{
  if(surface_type == BIHARMONIC) {
    return evaluate_mesh_biharmonic(mesh, reset);
  }
  else if (surface_type == GBS) {
    return evaluate_mesh_gbs(mesh, reset);
  }
  else {
    std::cerr << "Unknown surface type in evaluate_mesh" << std::endl;
    return false;
  }
}

bool SurfGBS::evaluate_mesh_gbs(Mesh& mesh, bool reset)
{
  if (reset) {
    mesh = Mesh(meshDomain);
  }
  for (const auto v : mesh.vertices()) {

    OpenMesh::Vec3d pt(0.0, 0.0, 0.0);
    double sum = 0.0;
    for (size_t loop = 0; loop < num_loops; ++loop) {
      for (size_t side = 0; side < num_sides[loop]; ++side) {
        for (size_t row = 0; row < num_rows[loop][side]; ++row) {
          for (size_t col = 0; col < num_cols[loop][side]; ++col) {
            const auto Bf = blend_functions[v.idx()][loop][side][row][col];
            const auto cp = OpenMesh::Vec3d(ribbons[loop][side].controlPoint(col, row).data());
            pt += Bf * cp;
            sum += Bf;
          }
        }
      }
    }
    pt /= sum;
    mesh.point(mesh.vertex_handle(v.idx())) = pt;
  }

  mesh.request_face_normals();
  mesh.update_face_normals();
  mesh.request_vertex_normals();
  mesh.update_vertex_normals();

  return true;
}


double SurfGBS::get_mu(const Mesh::VertexHandle& vtx, size_t loop, size_t side, size_t row, size_t col) const
{
  const auto eps = 1e-10;

  const auto side_m1 = prev(loop, side);
  const auto side_p1 = next(loop, side);

  const auto& hs = side_deformed[loop][side] ? h_coords_deformed : h_coords;
  const auto& hs_m1 = side_deformed[loop][side_m1] ? h_coords_deformed : h_coords;
  const auto& hs_p1 = side_deformed[loop][side_p1] ? h_coords_deformed : h_coords;

  const auto h = hs[vtx.idx()][loop][side];
  const auto hm1 = hs_m1[vtx.idx()][loop][side_m1];
  const auto hp1 = hs_p1[vtx.idx()][loop][side_p1];

  const auto p = num_rows[loop][side];
  const auto alpha = pow(hm1, p) / (pow(hm1, p) + pow(h, p));
  const auto beta = pow(hp1, p) / (pow(hp1, p) + pow(h, p));

  double mu = 1.0;
  if ((hm1 < eps && h < eps) || (h < eps && hp1 < eps)) {
    mu = 0.5;
  }
  else {
    if (col < num_rows[loop][side_m1]) {
      mu = alpha;
    }
    if (col >= num_cols[loop][side] - num_rows[loop][side_p1]) {
      mu = beta;
    }
  }

  return mu;

}

void SurfGBS::projectCurves2Domain(
  const std::vector<std::vector<Eigen::Vector3d>>& curves_xyz,
  std::vector<std::vector<Eigen::Vector3d>>& curves_uv
) const {
  if (meshDomain.n_faces() != meshSurface.n_faces() || meshDomain.n_faces() == 0) {
    return;
  }
  curves_uv.clear();
  curves_uv.resize(curves_xyz.size());
  for (size_t side = 0; side < curves_xyz.size(); ++side) {
    size_t num_pts = curves_xyz[side].size();
    curves_uv[side].resize(num_pts);
    for (size_t ii = 0; ii < num_pts; ++ii) {
      Eigen::Vector3d pt = curves_xyz[side][ii];
      auto closest_f = meshSurface.findClosestFace({ pt[0], pt[1], pt[2] });
      auto closest_uv = meshDomain.interpolateInFace(closest_f.face, closest_f.pt_bary);
      curves_uv[side][ii] = {closest_uv[0], closest_uv[1], 0.0}; // project2Triangle_uv(pt, closest_f);
    }
  }
}

Eigen::Vector3d SurfGBS::project2Triangle_uv(Eigen::Vector3d pt, Mesh::FaceHandle ff) const
{
  auto f = meshSurface.make_smart(ff);
  auto v1 = f.halfedge().to();
  auto v2 = f.halfedge().next().to();
  auto v3 = f.halfedge().next().next().to();
  auto p1_xyz = meshSurface.point(v1);
  auto p2_xyz = meshSurface.point(v2);
  auto p3_xyz = meshSurface.point(v3);
  auto p1_uv = meshDomain.point(v1);
  auto p2_uv = meshDomain.point(v2);
  auto p3_uv = meshDomain.point(v3);

  // Project point to the plane of the triangle
  auto pp = Mesh::Point(pt[0], pt[1], pt[2]) - (meshSurface.normal(ff) | (Mesh::Point(pt[0], pt[1], pt[2]) - p1_xyz)) * meshSurface.normal(ff);

  // Barycentric coordinates
  double u, v, w;
  Mesh::barycentricCoordinates(pp, p1_xyz, p2_xyz, p3_xyz, u, v, w);

  auto pt_uv = u * p1_uv + v * p2_uv + w * p3_uv;
  return { pt_uv[0], pt_uv[1], 0.0 };
}

void SurfGBS::findPCAPlane(
  const std::vector<std::vector<Eigen::Vector3d>>& curves_xyz,
  Eigen::Vector3d& origin,
  Eigen::Vector3d& normal,
  Eigen::Vector3d& axis_u,
  Eigen::Vector3d& axis_v)
{
  // Compute centroid
  origin = Eigen::Vector3d::Zero();
  size_t total_points = 0;
  for (const auto& curve : curves_xyz) {
    for (const auto& pt : curve) {
      origin += pt;
      total_points++;
    }
  }
  origin /= total_points;

  // Compute covariance matrix
  Eigen::Matrix3d cov = Eigen::Matrix3d::Zero();
  for (const auto& curve : curves_xyz) {
    for (const auto& pt : curve) {
      Eigen::Vector3d centered = pt - origin;
      cov += centered * centered.transpose();
    }
  }
  cov /= total_points;

  // Eigen decomposition
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(cov);
  if (solver.info() != Eigen::Success) {
    std::cerr << "Error: Eigen decomposition failed!" << std::endl;
    return;
  }

  // The normal is the eigenvector corresponding to the smallest eigenvalue
  normal = solver.eigenvectors().col(0).normalized();
  axis_u = solver.eigenvectors().col(2).normalized(); // Largest eigenvalue
  axis_v = solver.eigenvectors().col(1).normalized(); // Second largest eigenvalue

  // Ensure right-handed coordinate system
  if ((axis_u.cross(axis_v)).dot(normal) < 0) {
    axis_v = -axis_v;
  }
}

void SurfGBS::projectCurves2LSPlane(
  const std::vector<std::vector<Eigen::Vector3d>>& curves_xyz,
  std::vector<std::vector<Eigen::Vector3d>>& curves_uv
) {
  
  curves_uv.clear();
  curves_uv.resize(curves_xyz.size());
  Eigen::Vector3d plane_origin, plane_normal, plane_axis_u, plane_axis_v;
  findPCAPlane(curves_xyz, plane_origin, plane_normal, plane_axis_u, plane_axis_v);
  std::cout << "PCA Plane origin: " << plane_origin.transpose() << ", normal: " << plane_normal.transpose() << std::endl;
  for (size_t side = 0; side < curves_xyz.size(); ++side) {
    size_t num_pts = curves_xyz[side].size();
    curves_uv[side].resize(num_pts);
    for (size_t ii = 0; ii < num_pts; ++ii) {
      Eigen::Vector3d pt = curves_xyz[side][ii];
      Eigen::Vector3d to_pt = pt - plane_origin;
      Eigen::Vector3d uv = to_pt - (plane_normal.dot(to_pt)) * plane_normal;
      double u = uv.dot(plane_axis_u);
      double v = uv.dot(plane_axis_v);
      curves_uv[side][ii] = {u, v, 0.0};
    }
  }
}


void SurfGBS::scale_perimeter_ribbons(std::vector<std::vector<Ribbon> >& perimeter_ribbons) const
{
  const double scale = 1.0;
  for (size_t loop = 0; loop < perimeter_ribbons.size(); ++loop) {
    for (size_t side = 0; side < perimeter_ribbons[loop].size(); ++side) {
      auto& rib = perimeter_ribbons[loop][side];
      size_t nc = rib.numControlPoints()[0];
      for (size_t col = 0; col < nc; ++col){
        auto cp0 = rib.controlPoint(col, 0);
        auto cp1 = rib.controlPoint(col, 1);
        auto d = cp1 - cp0;
        rib.controlPoint(col, 1) = cp0 + d*scale;
      }
    }
  }
}

double SurfGBS::calc_inner_loop_scale(const std::vector<std::vector<Eigen::Vector3d>>& loop, const std::vector<std::vector<Eigen::Vector3d>>& normals, const std::vector<std::vector<Eigen::Vector3d>>& tangents) 
{
  double sk = 0.0;
  double reff = 0.0;
  double P = 0.0;
  size_t total_points = 0;
  // Computing the scale factor based on normals and tangents
  Eigen::Vector3d plane_origin, plane_normal, plane_axis_u, plane_axis_v;
  findPCAPlane(loop, plane_origin, plane_normal, plane_axis_u, plane_axis_v);
  for(size_t side = 0; side < loop.size(); ++side) {
    for(size_t i = 0; i < loop[side].size(); ++i) {
      const auto p = loop[side][i];
      const auto n = normals[side][i];
      const auto t = tangents[side][i];
      sk += std::abs(plane_normal.cross(t).dot(n)) / (1 + std::abs(plane_normal.dot(n)));
      total_points++;
    }
    for(size_t i = 0; i < loop[side].size() - 1; ++i) {
      const auto p = loop[side][i];
      const auto p1 = loop[side][i + 1];
      const auto c = (p + p1) * 0.5;
      const auto l = (p1 - p).norm();
      const auto r = (c - plane_origin).norm();
      reff += r * l;
      P += l;
    }
  }
  sk /= total_points;
  reff /= P;
  return sk / reff;
}

double SurfGBS::calc_inner_loop_scale(size_t loop) const
{
  size_t total_points = 0; 
  double sk = 0.0;
  std::vector<double> skv;
  const auto cog = GeomUtil::centroid(boundary_points[loop]);
  Eigen::Vector3d plane_origin, plane_normal, plane_axis_u, plane_axis_v;
  findPCAPlane(boundary_points[loop], plane_origin, plane_normal, plane_axis_u, plane_axis_v);

  for(size_t side = 0; side < ribbons[loop].size(); ++side) {
    const auto &rib = ribbons[loop][side];
    const auto &points = boundary_points[loop][side];
    const auto &normals = boundary_normals[loop][side];
    const auto &tangents = boundary_tangents[loop][side];
    const auto &crossderivatives = boundary_crossderivatives[loop][side];
    for(size_t i = 0; i < points.size(); ++i) {
      const auto p = points[i];
      const auto n = normals[i];
      const auto t = tangents[i];
      const auto cd = crossderivatives[i];
      const auto cp0 = p;
      const auto cp1 = (p + cd);
      size_t min_l = 0;
      size_t min_s = 0;
      size_t min_p = 0;
      double min_d = std::numeric_limits<double>::max();
      for (size_t l = 0; l < num_loops; ++l) {
        for (size_t s = 0; s < num_sides[l]; ++s) {
          if (l == loop) {
            continue;
          }
          for(size_t pi = 0; pi < boundary_points[l][s].size(); ++pi) {
            const auto pt = boundary_points[l][s][pi];
            const auto ni = boundary_normals[l][s][pi];
            const auto dv = (pt- cp0);
            const auto di = dv.norm();
            if(di < min_d /*&& dv.normalized().dot(n) > 0 */) {
              min_d = di;
              min_l = l;
              min_s = s;
              min_p = pi;
            }
          }
        }
      }
      const auto cp3 = boundary_points[min_l][min_s][min_p];
      const auto cdi = boundary_crossderivatives[min_l][min_s][min_p];
      const auto cp2 = (cp3 + cdi);
      const Geometry::BSCurve spline(
        3, 
        { 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0 }, 
        { 
          {cp0[0], cp0[1], cp0[2]}, 
          {cp1[0], cp1[1], cp1[2]}, 
          {cp2[0], cp2[1], cp2[2]}, 
          {cp3[0], cp3[1], cp3[2]} 
        }
      ); // Constructing a cubic Hermite

      const double al = spline.arcLength(0.0, 1.0);
      // SubCurve3D spline_pts;
      // for(size_t ii = 0; ii < 100; ++ii) {
      //   auto pp = spline.eval(double(ii) / 99.0);
      //   spline_pts.push_back({pp[0], pp[1], pp[2]});
      // }
      //writeLoops( { {spline_pts} } , "closest_" + std::to_string(loop) + "_" + std::to_string(side) + "_" + std::to_string(i) + ".obj");
        //skv.push_back(min_d / al);

      const auto c2 = (cp3 - cp0).squaredNorm();
      const auto b2 = std::pow((cp3 - cp0).dot(plane_normal), 2);
      const auto a2 = c2 - b2;
      const auto a = std::sqrt(a2);
      const double delta = std::abs(al - a);
      const double R = (cog - cp0).norm();
      const double alpha = 1.0 - delta / R;
      sk += alpha;
      std::cout << "Loop " << loop << " side " << side << " point " << i 
      << " closest to loop " << min_l << " side " << min_s << " point " << min_p 
        << " distance: " << min_d << " arc length: " << al << " delta: " << delta << 
        " R: " << R << " a: " << a <<  " ratio: " << alpha << std::endl;
      skv.push_back(alpha);
      //sk = std::min(min_d / al, sk);
      //sk += min_d / al;
      ++total_points;
    }
  }
  sk /= total_points;
  std::sort(skv.begin(), skv.end());
  // return sk;//std::max(0.5, sk);
  const auto s_min = skv.front();
  const auto s_median = skv[skv.size()/2];
  const auto s_max = skv.back();
  return std::max(s_min, 0.5) ; // returning median
}

void SurfGBS::resolve_self_intersections(size_t max_iter)
{
  // Compute outer loop bounding box
  Eigen::Vector3d bb_min = Eigen::Vector3d::Constant(std::numeric_limits<double>::max());
  Eigen::Vector3d bb_max = Eigen::Vector3d::Constant(std::numeric_limits<double>::lowest());
  for (const auto& side : domain_boundary_curves[0]) {
    for (const auto& pt : side) {
      bb_min = bb_min.cwiseMin(pt);
      bb_max = bb_max.cwiseMax(pt);
    }
  }
  double bb_diag = (bb_max - bb_min).norm();

  const double scales[4] = {0.8, 0.6, 0.5, 0.4};
  double scale = 1.0;
  auto domain_boundary_curves_copy = domain_boundary_curves;
  bool self_intersection = false;
  size_t num_iters = 0;
  do {
    self_intersection = false;
    for(size_t loop = 1; loop < num_loops; ++loop) { // Skip inner loops
      auto& pts = domain_boundary_curves[loop];
        
      //Checking all other loops and sides
        int intersection_idx = -1;
        int intersection_side = -1;
        for(size_t other_loop = 0; other_loop < num_loops; ++other_loop) {
          if (other_loop == loop) {
            continue;
          }
          const auto & other_pts = domain_boundary_curves[other_loop];
          for(size_t other_side = 0; other_side < num_sides[other_loop]; ++other_side) {
            for(size_t side = 0; side < num_sides[loop]; ++side){
              for(size_t i = 0; i < pts[side].size(); ++i) {
                for(size_t j = 0; j < other_pts[other_side].size(); ++j) {
                  if( (pts[side][i] - other_pts[other_side][j]).norm() < 0.1 * bb_diag /*16.0*target_length*/ ) {
                    intersection_idx = static_cast<int>(i);
                    intersection_side = static_cast<int>(side);
                  }
                }
              }
            }
          
        }

      }
      if (intersection_idx == -1) {
        continue;
      }
      else {
        std::cout << "Found self-intersection in loop " << loop << " side " << intersection_side << " at index " << intersection_idx << std::endl;
        self_intersection = true;
      }
      // Eigen::Vector3d dir = (pts[intersection_side][intersection_idx] - GeomUtil::centroid(domain_boundary_curves[loop])).normalized();
      // // Move pts in the direction of loop cog
      // for (size_t si = 0; si < pts.size(); ++si) {
      //   for (size_t pi = 0; pi < pts[si].size(); ++pi) {
      //     pts[si][pi] += dir * 15.0 * target_length;
      //   }
      // }
      // std::cout << "Resolved self-intersection in loop " << loop << " side " << intersection_side << std::endl;
    }


    if (self_intersection) {
      //std::cout << "Warning: Could not resolve self-intersections after 10 iterations!" << std::endl;
      domain_boundary_curves = domain_boundary_curves_copy;
      // Scaling down inner loops
      for(size_t loop = 1; loop < num_loops; ++loop) {
        const double scale = scales[std::min(num_iters, size_t(3))];
        std::cout << "Scaling down loop " << loop << " by factor " << scale << std::endl;
        const auto cog = GeomUtil::centroid(domain_boundary_curves[loop]);
        for (auto& sub : domain_boundary_curves[loop]) {
          for (auto& p : sub) {
            p = cog + (p - cog) * scale;
          }
        }
      }
    }
  } while (self_intersection && ++num_iters < max_iter);
}

void SurfGBS::position_inner_loops()
{
  if( num_loops <= 1) {
    return;
  }
  std::vector<SubCurve3D> loops_3D(num_loops);
  std::vector<SubCurve3D> crossderiv_3D(num_loops);
  std::vector<LoopFlattener::Loop2D> loops_2D(num_loops);

  std::vector<LoopFlattener::DistConstraint> constraints;

  // Collecting points for each loop
  for( size_t loop = 0; loop < num_loops; ++loop) {
    const auto& curve_loop_3D = boundary_points[loop];
    const auto& curve_crossderiv_3D = boundary_crossderivatives[loop];
    const auto& curve_loop_2D = domain_boundary_curves[loop];

    loops_2D[loop].fixed = (loop == 0); // Fixing outer loop
    for(size_t side = 0; side < num_sides[loop]; ++side) {
      const auto& curve_3D = curve_loop_3D[side];
      const auto& crossderiv = curve_crossderiv_3D[side];
      const auto& curve_2D = curve_loop_2D[side];
      for(size_t i = 0; i < curve_3D.size(); ++i) {
        loops_3D[loop].push_back(curve_3D[i]);
        crossderiv_3D[loop].push_back(crossderiv[i]);
        loops_2D[loop].ref.push_back({curve_2D[i].x(), curve_2D[i].y()});
      }
    }
  }

  std::map<std::array<size_t, 4>, bool> constraint_found;
  // Collecting distance constraints between loops
  for(size_t loop = 1; loop < num_loops; ++loop) {
    const auto& inner_loop_3D = loops_3D[loop];
    const auto& inner_crossderiv_3D = crossderiv_3D[loop];
    const auto& inner_loop_2D = loops_2D[loop];

    // Adding diameter constraints within the same loop
    const size_t num_inner_pts = inner_loop_2D.ref.size();
    for(size_t i = 0; i < num_inner_pts; i = i + 5) {
      size_t j = (i + num_inner_pts / 2) % num_inner_pts;
      if(i == j) {
        continue;
      }
      if(constraint_found.find({loop, i, loop, j}) != constraint_found.end()) {
        continue;
      }
      const auto& pt_i = inner_loop_3D[i];
      const auto& pt_j = inner_loop_3D[j];
      const double arclen = (pt_i - pt_j).norm();
      double dist = arclen;
      LoopFlattener::DistConstraint constraint;
      constraint.loopA = loop;
      constraint.loopB = loop;
      constraint.idxA = i;
      constraint.idxB = j;
      constraint.targetDist = dist; 
      constraints.push_back(constraint);
      constraint_found[{loop, i, loop, j}] = true;
    }

    for(size_t other_loop = 0; other_loop < num_loops; ++other_loop) {
      if(other_loop == loop) {
        continue;
      }
      const auto& outer_loop_3D = loops_3D[other_loop];
      const auto& outer_crossderiv_3D = crossderiv_3D[other_loop];
      const auto& outer_loop_2D = loops_2D[other_loop];
      for(size_t i = 0; i < inner_loop_2D.ref.size(); ++i) {
        // Finding closest point on outer loop
        double min_dist = std::numeric_limits<double>::max();
        Eigen::Vector3d closest_pt;
        size_t j_closest = 0;
        const auto& pt_inner = inner_loop_3D[i];
        for(size_t j = 0; j < outer_loop_2D.ref.size(); ++j) {
          if(constraint_found.find({loop, i, other_loop, j}) != constraint_found.end()) {
            continue;
          }
          const auto& pt_outer = outer_loop_3D[j];
          // Compute distance from inner to outer loop by cubic Hermite spline arclength
          const Eigen::Vector3d cd_inner = inner_crossderiv_3D[i];
          const Eigen::Vector3d cd_outer = outer_crossderiv_3D[j];
          auto spline = Geometry::BSCurve(
            3, 
            { 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0 }, 
            { 
              {pt_inner[0], pt_inner[1], pt_inner[2]}, 
              { (pt_inner + cd_inner)[0], (pt_inner + cd_inner)[1], (pt_inner + cd_inner)[2]}, 
              { (pt_outer + cd_outer)[0], (pt_outer + cd_outer)[1], (pt_outer + cd_outer)[2]}, 
              {pt_outer[0], pt_outer[1], pt_outer[2]} 
            }
          ); // Constructing a cubic Hermite
          const double arclen = spline.arcLength(0.0, 1.0);
          // Using Euclidean distance for now

          double dist = arclen; //(pt_inner - pt_outer).norm();
          if(dist < min_dist) {
            min_dist = dist;
            closest_pt = pt_outer;
            j_closest = j;
          }
        }
        LoopFlattener::DistConstraint constraint;
        constraint.loopA = loop;
        constraint.loopB = other_loop;
        constraint.idxA = i;
        constraint.idxB = j_closest;
        constraint.targetDist = min_dist; // Keeping some margin
        constraints.push_back(constraint);
        constraint_found[{loop, i, other_loop, j_closest}] = true;
        // // Exporting spline for visualization
        // const auto& pt_outer = outer_loop_3D[j_closest];
        // const auto cd_inner = inner_crossderiv_3D[i];
        // const auto cd_outer = outer_crossderiv_3D[j_closest];
        // auto spline = Geometry::BSCurve(
        //   3, 
        //   { 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0 }, 
        //   { 
        //     {pt_inner[0], pt_inner[1], pt_inner[2]}, 
        //     { (pt_inner + cd_inner)[0], (pt_inner + cd_inner)[1], (pt_inner + cd_inner)[2]}, 
        //     { (pt_outer + cd_outer)[0], (pt_outer + cd_outer)[1], (pt_outer + cd_outer)[2]}, 
        //     {pt_outer[0], pt_outer[1], pt_outer[2]} 
        //   }
        // ); // Constructing a cubic Hermite
        // writeLoops( { { {pt_inner, pt_inner + cd_inner, pt_outer + cd_outer, pt_outer} } }, "constraint_loop" + std::to_string(loop) + "_pt" + std::to_string(i) + 
        //   "_to_loop" + std::to_string(other_loop) + "_pt" + std::to_string(j_closest) + ".obj");
      }
    }
  }

  // // Printing contstraints
  // std::cout << "Distance constraints between loops:" << std::endl;
  // for( const auto& constraint : constraints) {
  //   std::cout << "  Loop " << constraint.loopA << " point " << constraint.idxA 
  //     << " to Loop " << constraint.loopB << " point " << constraint.idxB 
  //     << " target distance: " << constraint.targetDist << std::endl;
  // }

  LoopFlattener::solveHolePlacements(loops_2D, constraints, 100, 1E-3, 1E-9, 1E-13, true);
  std::cout << "Positioned inner loops using " << constraints.size() << " distance constraints." << std::endl;
  std::cout << "New inner loop positions:" << std::endl;
  for( size_t loop = 0; loop < num_loops; ++loop) {
    if( loop == 0 ) {
      continue;
    }
    std::cout << "  Loop " << loop << ": t = (" << loops_2D[loop].tx << ", " << loops_2D[loop].ty << "); theta = " << loops_2D[loop].theta << "; s = " << LoopFlattener::sigmoid(loops_2D[loop].u) << std::endl;
  }

  auto transformed_domain_boundary_curves = domain_boundary_curves;
  // Transforming domain boundary curves
  for( size_t loop = 1; loop < num_loops; ++loop) {
    auto& curve_loop_2D = transformed_domain_boundary_curves[loop];
    const auto& L = loops_2D[loop];
    const double scale = LoopFlattener::sigmoid(L.u);
    const double c = std::cos(L.theta), s = std::sin(L.theta);
    Eigen::Matrix2d R;
    R << c, -s,
         s,  c;
    const auto t = Eigen::Vector2d(L.tx, L.ty);
    for(size_t side = 0; side < num_sides[loop]; ++side) {
      auto& curve_2D = curve_loop_2D[side];
      for(size_t i = 0; i < curve_2D.size(); ++i) {
        const auto p = Eigen::Vector2d(curve_2D[i].x(), curve_2D[i].y());
        auto new_pt = scale*R*p + t;
        transformed_domain_boundary_curves[loop][side][i] = Eigen::Vector3d(new_pt[0], new_pt[1], 0.0);
      }
    }
  }
  //writeLoops(transformed_domain_boundary_curves, "transformed_domain_boundary_curves.obj");
  domain_boundary_curves = transformed_domain_boundary_curves;
}



void SurfGBS::merge_smooth_corners() {
  auto new_ribbons = ribbons;
  for(size_t loop = 0; loop < 1; ++loop) { // Skip inner loops
    new_ribbons[loop].clear();
    new_ribbons[loop].push_back(ribbons[loop].front());
    std::vector<size_t> num_merged_sides(1, 1);
    for(size_t side = 1; side < num_sides[loop]; ++side) {
      const auto side_m1 = prev(loop, side);
      // const auto side_p1 = next(loop, side);
      const auto rib1 = new_ribbons[loop].back();
      const auto rib2 = ribbons[loop][side];
      // Check angle between sides
      Geometry::VectorMatrix duv1, duv2;
      rib1.eval(1.0, 0.0, 1.0, duv1);
      rib2.eval(0.0, 0.0, 1.0, duv2);
      const double angle = std::acos(std::clamp(duv1[1][0].normalized() * duv2[1][0].normalized(), -1.0, 1.0)) * 180.0 / M_PI;
      // std::cout << "Loop " << loop << " sides " << side_m1 << " and " << side << " angle: "
      //  << angle << std::endl;
      if (angle < 1e-3) {
        // Create new ribbon
        const size_t deg_u1 = rib1.basisU().degree();
        const size_t deg_u2 = rib2.basisU().degree();
        const size_t deg_v = rib1.basisV().degree();
        if(deg_u1 != deg_u2) {
          std::cerr << "Error: Cannot merge ribbons with different degrees!" << std::endl;
          continue;
        }
        else {
          std::cout << "Merging ribbons at loop " << loop << " sides " << side_m1 << " and " << side << std::endl;
        }
        size_t num_segments = num_merged_sides.back() + 1;
        Geometry::DoubleVector knots_u, knots_v; // ToDo: Merging non-uniform knotvectors
        knots_u.insert(knots_u.end(), deg_u1 + 1, 0.0);
        for (size_t i = 1; i < num_segments; ++i) {
          knots_u.insert(knots_u.end(), deg_u1 - (c1_merge ? 1 : 0), i / double(num_segments));
        }
        knots_u.insert(knots_u.end(), deg_u1 + 1, 1.0);
        knots_v = rib1.basisV().knots();
        Geometry::PointVector cpts;
        for(size_t c = 0; c < rib1.numControlPoints().at(0) - (c1_merge ? 1 : 0); ++c) {
          for(size_t r = 0; r < rib1.numControlPoints().at(1); ++r) {
            cpts.push_back(rib1.controlPoint(c, r));
          }
        }
        for (size_t c = 1; c < rib2.numControlPoints().at(0); ++c) {
          for (size_t r = 0; r < rib2.numControlPoints().at(1); ++r) {
            cpts.push_back(rib2.controlPoint(c, r));
          }
        }
        Geometry::BSSurface new_rib(deg_u1, deg_v, knots_u, knots_v, cpts);
        new_ribbons[loop].back() = new_rib;
        ++num_merged_sides.back();
      }
      else {
        new_ribbons[loop].push_back(ribbons[loop][side]);
        num_merged_sides.push_back(1);
      }
    }

    {// Merge last and first ribbons if needed
      size_t side = 0;
      const auto side_m1 = prev(loop, side);
      const auto rib1 = new_ribbons[loop].back();
      const auto rib2 = new_ribbons[loop].front();
      // Check angle between sides
      Geometry::VectorMatrix duv1, duv2;
      rib1.eval(1.0, 0.0, 1.0, duv1);
      rib2.eval(0.0, 0.0, 1.0, duv2);
      const double angle = std::acos(std::clamp(duv1[1][0].normalized() * duv2[1][0].normalized(), -1.0, 1.0)) * 180.0 / M_PI;
      // std::cout << "Loop " << loop << " sides " << side_m1 << " and " << side << " angle: "
      //  << angle << std::endl;
      if (angle < 1e-3) {
        // Create new ribbon
        const size_t deg_u1 = rib1.basisU().degree();
        const size_t deg_u2 = rib2.basisU().degree();
        const size_t deg_v = rib1.basisV().degree();
        if (deg_u1 != deg_u2) {
          std::cerr << "Error: Cannot merge ribbons with different degrees!" << std::endl;
          continue;
        }
        else {
          std::cout << "Merging ribbons at loop " << loop << " sides " << side_m1 << " and " << side << std::endl;
        }
        size_t num_segments = num_merged_sides.back() + num_merged_sides.front();
        Geometry::DoubleVector knots_u, knots_v; // ToDo: Merging non-uniform knotvectors
        knots_u.insert(knots_u.end(), deg_u1 + 1, 0.0);
        for (size_t i = 1; i < num_segments; ++i) {
          knots_u.insert(knots_u.end(), deg_u1 - (c1_merge ? 1 : 0), i / double(num_segments));
        }
        knots_u.insert(knots_u.end(), deg_u1 + 1, 1.0);
        knots_v = rib1.basisV().knots();
        Geometry::PointVector cpts;
        for (size_t c = 0; c < rib1.numControlPoints().at(0) - (c1_merge ? 1 : 0); ++c) {
          for (size_t r = 0; r < rib1.numControlPoints().at(1); ++r) {
            cpts.push_back(rib1.controlPoint(c, r));
          }
        }
        for (size_t c = 1; c < rib2.numControlPoints().at(0); ++c) {
          for (size_t r = 0; r < rib2.numControlPoints().at(1); ++r) {
            cpts.push_back(rib2.controlPoint(c, r));
          }
        }
        Geometry::BSSurface new_rib(deg_u1, deg_v, knots_u, knots_v, cpts);
        new_ribbons[loop].front() = new_rib;
        new_ribbons[loop].pop_back();
        num_merged_sides.front() += num_merged_sides.back();
        num_merged_sides.pop_back();
      }
    }
    num_segments[loop] = num_merged_sides;
  }
  ribbons = new_ribbons;
  init_data();
}

void SurfGBS::compute_auto_h_widths(double turning_angle, double length_ratio)
{
  auto turning_angle_rad = turning_angle * M_PI / 180.0;
  for (size_t loop = 0; loop < num_loops; ++loop) {
    size_t num_sides = ribbons[loop].size();
    if (num_sides == 1 || loop >= 1) { continue; }
    for (size_t side = 0; side < num_sides; ++side) {
      size_t side_m1 = prev(loop, side);
      size_t side_p1 = next(loop, side);
      const auto& s = domain_boundary_curves[loop][side];
      auto s_m1 = domain_boundary_curves[loop][side_m1];
      //Reverse s_m1 to have same orientation
      s_m1 = SubCurve3D(s_m1.rbegin(), s_m1.rend());
      const auto& s_p1 = domain_boundary_curves[loop][side_p1];
      std::vector<double> lengths_m1, lengths_p1, lengths_i;
      GeomUtil::getEdgeLengths(s_m1, lengths_m1, false);
      GeomUtil::getEdgeLengths(s_p1, lengths_p1, false);
      GeomUtil::getEdgeLengths(s, lengths_i, false);
      double len_m1 = std::accumulate(lengths_m1.begin(), lengths_m1.end(), 0.0); // Summing up edge lengths
      double len_p1 = std::accumulate(lengths_p1.begin(), lengths_p1.end(), 0.0); // Summing up edge lengths
      double len_i = std::accumulate(lengths_i.begin(), lengths_i.end(), 0.0);
      // compute cumulative lenghts
      std::vector<double> cum_lengths_m1(1, 0.0);
      cum_lengths_m1.reserve(lengths_m1.size() + 1);
      std::vector<double> cum_lengths_p1(1, 0.0);
      cum_lengths_p1.reserve(lengths_p1.size() + 1);
      std::vector<double> cum_lengths_i(1, 0.0);
      cum_lengths_i.reserve(lengths_i.size() + 1);
      for (const auto& l : lengths_m1) { cum_lengths_m1.push_back(cum_lengths_m1.back() + l); }
      for (const auto& l : lengths_p1) { cum_lengths_p1.push_back(cum_lengths_p1.back() + l); }
      for (const auto& l : lengths_i) { cum_lengths_i.push_back(cum_lengths_i.back() + l); }

      // Now compute h_widths based on length ratio and turning angle
      for(size_t ii = 0; ii < h_widths[loop][side].size(); ++ii) {
        h_widths[loop][side][ii] = 1.0; // Reset to max
        
        if(ii == 0 && num_segments[loop][side_m1] <= 1) {
          // side - 1
          double max_ratio = 0.0;
          for (size_t jj = 1; jj < cum_lengths_m1.size(); ++jj) {
            double l_m1 = cum_lengths_m1[jj];
            
            double ratio = l_m1 / len_i;
            max_ratio = std::max(ratio, max_ratio);
            if (ratio >= length_ratio) {
              h_widths[loop][side][ii] = std::min(h_widths[loop][side][ii], l_m1 / len_m1);
              break;
            }
          }

          double max_angle = 0.0;
          auto t_0 = (s_m1[1] - s_m1[0]).normalized();
          for(size_t jj = 1; jj < s_m1.size(); ++jj) {
            auto t_i = (s_m1[jj] - s_m1[0]).normalized();
            double angle = abs(GeomUtil::getAngle(t_0, t_i));
            max_angle = std::max(angle, max_angle);
            if (angle >= turning_angle_rad) {
              h_widths[loop][side][ii] = std::min(h_widths[loop][side][ii], cum_lengths_m1[jj] / len_m1);
              break;
            }
          }
          // // Print max values
          // std::cout << "Loop " << loop << " Side " << side << " h-widths[" << ii << "] max length ratio: " << max_ratio 
          //   << ", max turning angle (deg): " << max_angle * 180.0 / M_PI << std::endl;
        }
        else if (ii == 1 && num_segments[loop][side_p1] <= 1) {
          // side + 1
          double max_ratio = 0.0;
          for (size_t jj = 1; jj < cum_lengths_p1.size(); ++jj) {
            double l_p1 = cum_lengths_p1[jj];
            double ratio = l_p1 / len_i;
            max_ratio = std::max(ratio, max_ratio);
            if (ratio >= length_ratio) {
              h_widths[loop][side][ii] = std::min(h_widths[loop][side][ii], cum_lengths_p1[jj] / len_p1);
              break;
            }
          }

          double max_angle = 0.0;
          auto t_0 = (s_p1[1] - s_p1[0]).normalized();
          for(size_t jj = 1; jj < s_p1.size(); ++jj) {
            auto t_i = (s_p1[jj] - s_p1[0]).normalized();
            double angle = abs(GeomUtil::getAngle(t_0, t_i));
            max_angle = std::max(angle, max_angle);
            if (angle >= turning_angle_rad) {
              h_widths[loop][side][ii] = std::min(h_widths[loop][side][ii], cum_lengths_p1[jj] / len_p1);
              break;
            }
          }
          // // Print max values
          // std::cout << "Loop " << loop << " Side " << side << " h-widths[" << ii << "] max length ratio: " << max_ratio 
          //   << ", max turning angle (deg): " << max_angle * 180.0 / M_PI << std::endl;
        }
      }
    }
  }

  // // Print h-widths
  // for (size_t loop = 0; loop < num_loops; ++loop) {
  //   for (size_t side = 0; side < num_sides[loop]; ++side) {
  //     std::cout << "Loop " << loop << " Side " << side << " h-widths: ";
  //     for (const auto& hw : h_widths[loop][side]) {
  //       std::cout << hw << " ";
  //     }
  //     std::cout << std::endl;
  //   }
  // }
}

void SurfGBS::find_init_placement_scales()
{
  init_placement_scales.clear();
  init_placement_scales.resize(num_loops);
  for(size_t loop = 0; loop < num_loops; ++loop) {
    init_placement_scales[loop].resize(num_sides[loop], 1.0);
  }

  // TODO: Implement constant initial scaling based on distance of h = 0.33 points from closest pts on boundary
  
}

inline size_t circular_index(size_t i, int offset, size_t n) {
  return (i + n + (offset % static_cast<int>(n))) % n;
}

size_t SurfGBS::prev(size_t loop, size_t side) const {
  return circular_index(side, -1, num_sides[loop]);
}

size_t SurfGBS::next(size_t loop, size_t side) const {
  return circular_index(side, +1, num_sides[loop]);
}

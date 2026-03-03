#pragma once
#include <geometry.hh>
#include "Mesh.hpp"
#include <Eigen/Eigen>
#include "SimpleDomain.hpp"

namespace libcdgbs {
  class SurfGBS {
  public:
    using Ribbon = Geometry::BSSurface;
    using Point2D = Geometry::Point2D;
    using Point3D = Geometry::Point3D;
    using VertexHandle = Mesh::VertexHandle;

    enum DomainType {
      CURVED,
      POLYGONAL,
      REGULAR
    } domain_type;

    typedef enum {
      GBS,
      BIHARMONIC
    } SurfaceType;

    class InputParams {
    public:
        double target_length;
        bool merge_smooth_corners;
        double deform_value;
        bool restrict_params;
        bool c1_merge;
        double global_inner_loop_scale;
        bool use_h_widths;
        bool arclength_sampling;
        bool merge_inner_c1;
        bool ribbon_init_placement;

        InputParams() : target_length(3.0),
            merge_smooth_corners(true),
            deform_value(0.0),
            restrict_params(true),
            c1_merge(false),
          global_inner_loop_scale(1.0),
          use_h_widths(true),
          arclength_sampling(true),
          merge_inner_c1(false),
          ribbon_init_placement(false){
        };

        InputParams(double target_length_,
            bool merge_smooth_corners_,
            double deform_value_,
            bool restrict_params_,
            bool c1_merge_,
            double global_inner_loop_scale_,
            bool use_h_widths_,
            bool arclength_sampling_,
            bool merge_inner_c1_,
            bool ribbon_init_placement_
        ) : target_length(target_length_),
            merge_smooth_corners(merge_smooth_corners_),
            deform_value(deform_value_),
            restrict_params(restrict_params_),
            c1_merge(c1_merge_),
          global_inner_loop_scale(global_inner_loop_scale_),
          use_h_widths(use_h_widths_),
          arclength_sampling(arclength_sampling_),
          merge_inner_c1(merge_inner_c1_),
          ribbon_init_placement(ribbon_init_placement_){
        };
    };

    std::vector<std::vector<Ribbon> > ribbons;
    std::vector<std::vector<size_t> > num_segments;
    std::vector<bool> periodic;

    Mesh meshDomain;
    Mesh meshSurface;

    std::vector<std::vector<std::vector<Eigen::Vector3d> > > boundary_points;
    std::vector<std::vector<std::vector<Eigen::Vector3d> > > boundary_normals;
    std::vector<std::vector<std::vector<Eigen::Vector3d> > > boundary_tangents;
    std::vector<std::vector<std::vector<Eigen::Vector3d> > > boundary_tangents_unnormalized;
    std::vector<std::vector<std::vector<Eigen::Vector3d> > > boundary_crossderivatives;

    std::vector<std::vector<std::vector<Eigen::Vector3d> > > developed_boundary_curves;
    std::vector<std::vector<std::vector<Eigen::Vector3d> > > developed_boundary_curves_normalized;
    std::vector<std::vector<std::vector<Eigen::Vector3d> > > domain_boundary_curves;
    std::vector<std::vector<std::vector<VertexHandle> > > domain_boundary_vertices;
    std::vector<std::vector<std::vector<double> > > domain_boundary_params;

    SimpleDomain simple_domain;

    size_t num_loops;
    std::vector<size_t> num_sides;
    std::vector<std::vector<size_t>> num_rows;
    std::vector<std::vector<size_t>> num_cols;
    std::vector<std::vector<size_t>> deg_h;

    double target_length = 3.0;
    std::vector<std::vector<size_t> > side_res;
    std::vector<std::vector<std::vector<size_t> > > side_segment_res; // indices of corners to be merged

    std::vector<std::vector<std::array<double, 2> > > h_widths;

    std::vector<std::vector<std::vector<double> > > s_coords;
    std::vector<std::vector<std::vector<double> > > h_coords;
    std::vector<std::vector<std::vector<double> > > h_coords_deformed;
    std::vector<std::vector<std::vector<std::vector<std::vector<double> > > > > blend_functions;

    bool restrict_params = true;

    std::vector<std::vector<Ribbon> > deform_splines;
    std::vector<std::vector<bool > > side_deformed;
    double deform_value = 0.0;

    double global_inner_loop_scale = 1.0;

    std::vector<std::vector<double> > init_placement_scales;

    bool use_h_widths = true;

    bool arclength_sampling = true;

    bool merge_inner_c1 = false;

    bool ribbon_init_placement = false;

    bool debug_outputs = false;

    SurfGBS();

    void load_ribbons(const std::vector<std::vector<Ribbon> >& ribbon_surfs, const InputParams& params);
    void load_ribbons_and_evaluate(const std::vector<std::vector<Ribbon> >& ribbon_surfs, Mesh& mesh, const InputParams& params, SurfaceType surface_type = GBS);

    void init_data();

    bool readMGBS(const std::string& filename, const InputParams& params = InputParams());
    bool readMLP(const std::string& filename, const InputParams& params = InputParams());
    bool readCGB(const std::string& filename, const InputParams& params = InputParams());
    bool readGBS(const std::string& filename, const InputParams& params = InputParams());
    bool readGBP(const std::string& filename, const InputParams& params = InputParams());
    bool readNGBS(const std::string& filename, const InputParams& params = InputParams());

    bool writeOBJ(const Mesh& mesh, const std::string& filename);
    bool writeMGBS(const std::string& filename) const;
    static bool writeLoops(const std::vector<std::vector<std::vector<Eigen::Vector3d> > >& loops, const std::string filename = "/tmp/boundary.obj", bool edges = true);

    bool compute_domain_boundary();
    bool compute_domain_mesh();
    bool compute_local_parameters();
    bool compute_blend_functions();
    bool evaluate_mesh(Mesh& mesh, bool reset = true, SurfaceType surface_type = GBS);
    bool evaluate_mesh_gbs(Mesh& mesh, bool reset = true);
    bool evaluate_mesh_biharmonic(Mesh& mesh, bool reset = true);

    bool compute_harmonic_parameters();
    bool compute_deformed_parameters();

    static std::vector<double> periodicC1CubicBSplineBasis(double u, size_t N);

    void projectCurves2Domain(
      const std::vector<std::vector<Eigen::Vector3d>>& curves_xyz, 
      std::vector<std::vector<Eigen::Vector3d>>& curves_uv
    ) const;
    Eigen::Vector3d project2Triangle_uv(Eigen::Vector3d pt, Mesh::FaceHandle ff) const;

    static void findPCAPlane(
      const std::vector<std::vector<Eigen::Vector3d>>& curves_xyz,
      Eigen::Vector3d& origin, 
      Eigen::Vector3d& normal, 
      Eigen::Vector3d& axis_u, 
      Eigen::Vector3d& axis_v 
    );
    static void projectCurves2LSPlane(
      const std::vector<std::vector<Eigen::Vector3d>>& curves_xyz,
      std::vector<std::vector<Eigen::Vector3d>>& curves_uv
    );
    void scale_perimeter_ribbons(std::vector<std::vector<Ribbon> >& perimeter_ribbons) const;
    static double calc_inner_loop_scale(const std::vector<std::vector<Eigen::Vector3d>>& loop, const std::vector<std::vector<Eigen::Vector3d>>& normals, const std::vector<std::vector<Eigen::Vector3d>>& tangents);
    double calc_inner_loop_scale(size_t loop) const;
    void resolve_self_intersections(size_t max_iter = 4);
    void position_inner_loops();

    void find_init_placement_scales();

    void merge_smooth_corners();
    bool c1_merge = false;

    void compute_auto_h_widths(double turning_angle = 30.0, double length_ratio = 2.0);

    size_t prev(size_t loop, size_t side) const;
    size_t next(size_t loop, size_t side) const;


    // Point3D evaluate(const Point2D& uv) const;
    // Point3D evaluate(const Mesh::VertexHandle& vtx) const;

    // double get_s_coord(const Point2D& uv, size_t loop, size_t side) const;
    // double get_s_coord(const Mesh::VertexHandle& vtx, size_t loop, size_t side) const;

    // double get_h_coord(const Point2D& uv, size_t loop, size_t side) const;
    // double get_h_coord(const Mesh::VertexHandle& vtx, size_t loop, size_t side) const;

    //double get_mu(const Point2D& uv, size_t loop, size_t side) const;
    double get_mu(const VertexHandle& vtx, size_t loop, size_t side, size_t row, size_t col) const;
  };
}

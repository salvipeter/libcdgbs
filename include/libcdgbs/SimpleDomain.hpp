#pragma once

#include <array>
#include <variant>
#include <vector>

#include <Eigen/Core>

class SimpleDomain {
public:
  using Vec3 = Eigen::Vector3d;
  using Edge = std::vector<Vec3>;
  using Loop = std::vector<Edge>;
  using Segment = std::array<Vec3, 2>;
  using Parabola = std::array<Vec3, 3>;
  using EdgeCurve = std::variant<Segment, Parabola>;
  struct Circle { Eigen::Vector3d center; double radius, start_angle; size_t n_pts; };

  // Create the domain
  // - loops[i][j] is a polyline representing side j in loop i
  // These polylines are modified(!):
  // - On the boundary, endpoints are preserved
  // - All other points are moved to uniformly sampled positions of the underlying geometry
  void init(std::vector<Loop> &loops);

  // Compute s/h local parameters
  // - p is the point of evaluation
  // - s[i][j] & h[i][j] are the parameters associated with side j in loop i
  void computeParameters(const Vec3 &p,
                         std::vector<std::vector<double>> &s,
                         std::vector<std::vector<double>> &h) const;

private:
  std::vector<size_t> num_sides;      // Number of sides in each loop
  std::vector<double> num_points;     // Number of points in each side of loop 0
  std::vector<EdgeCurve> boundaries;  // Segment/parabola of each side in loop 0
  std::vector<Circle> holes;          // Circle representing loop k (0th element not used!)
};

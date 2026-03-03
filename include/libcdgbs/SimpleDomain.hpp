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

  void init(std::vector<SimpleDomain::Loop> &loops);
  void computeParameters(const Vec3 &p, const std::vector<size_t> &num_sides,
                         std::vector<std::vector<double>> &s,
                         std::vector<std::vector<double>> &h) const;

private:
  std::vector<double> num_points;
  std::vector<EdgeCurve> boundaries;
  std::vector<Circle> holes;
};

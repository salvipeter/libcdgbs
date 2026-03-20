#include "libcdgbs/SimpleDomain.hpp"

#include <algorithm>
#include <exception>
#include <limits>

#include <Eigen/Dense>

using Vec2 = Eigen::Vector2d;
using Vec3 = SimpleDomain::Vec3;
using Segment = SimpleDomain::Segment;
using Parabola = SimpleDomain::Parabola;
using EdgeCurve = SimpleDomain::EdgeCurve;

namespace {

// Selects the real root of smallest absolute value
double selectSmallestRoot(const Eigen::VectorXcd &roots, double tol = 1e-9) {
  double best = 0;
  double min_abs = std::numeric_limits<double>::max();
  for (int i = 0; i < roots.size(); ++i)
    if (std::abs(roots[i].imag()) < tol) {
      double r = roots[i].real();
      if (std::abs(r) < min_abs) {
        min_abs = std::abs(r);
        best = r;
      }
    }
  return best;
}

// Selects the (assumed to be only) real root in [0,1]
double select01Root(const Eigen::VectorXcd &roots, double imagtol = 1e-9, double tol = 1e-6) {
  for (int i = 0; i < roots.size(); ++i)
    if (std::abs(roots[i].imag()) < imagtol &&
        roots[i].real() > -tol && roots[i].real() < 1 + tol)
      return std::clamp(roots[i].real(), 0.0, 1.0);
  throw std::runtime_error("no root in [0,1] found");
}

// Solves cN x^N + ... + c1 x + c0 = 0
template<size_t N>
auto generalSolver(const std::array<double, N+1> &coeffs) {
  if (std::abs(coeffs[N]) < 1e-10) {
    if constexpr (N == 1)
      throw std::runtime_error("constant polynomial");
    else {
      std::array<double, N> smaller;
      std::copy(coeffs.begin(), coeffs.begin() + N, smaller.begin());
      return generalSolver<N-1>(smaller);
    }
  }
  Eigen::MatrixXd companion = Eigen::MatrixXd::Zero(N, N);
  for (size_t i = 0; i < N; ++i) {
    if (i > 0)
      companion(i, i - 1) = 1;
    companion(i, N - 1) = -coeffs[i] / coeffs[N];
  }
  Eigen::EigenSolver<Eigen::MatrixXd> es(companion);
  return es.eigenvalues();
}

template<size_t N, size_t M>
double det2d(const Eigen::Vector<double, N> &a, const Eigen::Vector<double, M>  &b) {
  return a[0] * b[1] - a[1] * b[0];
};

[[maybe_unused]]
std::pair<double, Parabola> parabolaThroughPoint(const Parabola &p, const Vec3 &q) {
  auto getNormal = [](const Vec3 &v) { return Vec2(-v[1], v[0]).normalized(); };
  auto n0 = getNormal(p[1] - p[0]);
  auto n2 = getNormal(p[2] - p[1]);

  // Solve system: n0.dot(m) = 1, n2.dot(m) = 1
  Eigen::Matrix2d mat;
  mat << n0.transpose(), n2.transpose();
  Vec2 b(1.0, 1.0);
  Vec2 m = mat.colPivHouseholderQr().solve(b);

  // 3. Define coefficients for Area Polynomials: L(d) = c2*d^2 + c1*d + c0
  auto getAreaPoly = [&](const Vec3 &Ps, const Vec2 &vs, const Vec3 &Pe, const Vec2 &ve) {
    Vec3 W = Pe - Ps, U = q - Ps; Vec2 V = ve - vs;
    double c0 = det2d<3,3>(W, U);
    double c1 = det2d<2,3>(V, U) - det2d<3,2>(W, vs);
    double c2 = -det2d<2,2>(V, vs);
    return Vec3(c0, c1, c2);
  };
  auto L1 = getAreaPoly(p[0], n0, p[1], m);
  auto L2 = getAreaPoly(p[1], m,  p[2], n2);
  auto L3 = getAreaPoly(p[0], n0, p[2], n2);

  // 4. Expand 4*L1*L2 - L3^2 to get Quartic Coefficients [a0, a1, a2, a3, a4]
  // (A+Bx+Cx^2)(D+Ex+Fx^2) expansion
  std::array<double, 5> coeffs;
  coeffs[0] = 4 * (L1[0] * L2[0]) - (L3[0] * L3[0]);
  coeffs[1] = 4 * (L1[0] * L2[1] + L1[1] * L2[0]) - (2 * L3[0] * L3[1]);
  coeffs[2] = 4 * (L1[0] * L2[2] + L1[1] * L2[1] + L1[2] * L2[0]) -
    (L3[1] * L3[1] + 2 * L3[0] * L3[2]);
  coeffs[3] = 4 * (L1[1] * L2[2] + L1[2] * L2[1]) - (2 * L3[1] * L3[2]);
  coeffs[4] = 4 * (L1[2] * L2[2]) - (L3[2] * L3[2]);

  auto best_d = selectSmallestRoot(generalSolver<4>(coeffs));

  return {best_d, { p[0] + best_d * Vec3(n0[0], n0[1], 0),
                    p[1] + best_d * Vec3(m[0],  m[1],  0),
                    p[2] + best_d * Vec3(n2[0], n2[1], 0) }};
}

[[maybe_unused]]
Parabola parabolize(const EdgeCurve &c) {
  if (std::holds_alternative<Parabola>(c))
    return std::get<Parabola>(c);
  const auto &s = std::get<Segment>(c);
  return { s[0], (s[0] + s[1]) / 2, s[1] };
}

struct ImplicitConic {
  double A, B, C, D, E, F;
  double eval(double x, double y) const {
    return A * x * x + B * x * y + C * y * y + D * x + E * y + F;
  }
  Vec2 grad(double x, double y) const {
    return { 2 * A * x + B * y + D, 2 * C * y + B * x + E };
  }
  double normeval(double x, double y) const {
    return eval(x, y) / grad(x, y).norm();
  }
  ImplicitConic operator+(const ImplicitConic &c) const {
    return { A + c.A, B + c.B, C + c.C, D + c.D, E + c.E, F + c.F };
  }
  ImplicitConic operator*(double s) const {
    return { A * s, B * s, C * s, D * s, E * s, F * s };
  }
};

// Based on Sederberg'84
auto implicitize(const Parabola &p) {
  Vec3 q0 = p[0], q1 = 2 * (p[1] - p[0]), q2 = p[0] - 2 * p[1] + p[2];
  double a0 = q0[0], a1 = q1[0], a2 = q2[0];
  double b0 = q0[1], b1 = q1[1], b2 = q2[1];
  // A x^2 + B xy + C y^2 + D x + E y + F = 0
  auto A = b2 * b2, B = -2 * a2 * b2, C = a2 * a2;
  auto D = -2 * a0 * b2 * b2 + a1 * b1 * b2 - a2 * b1 * b1 + 2 * a2 * b0 * b2;
  auto E = -2 * b0 * a2 * a2 + b1 * a1 * a2 - b2 * a1 * a1 + 2 * b2 * a0 * a2;
  Eigen::Matrix4d m;
  m << a2, a1, a0, 0,
        0, a2, a1, a0,
       b2, b1, b0, 0,
        0, b2, b1, b0;
  auto F = m.determinant();
  return ImplicitConic{A, B, C, D, E, F};
}

auto implicitize(const Segment &p) {
  Vec3 d = p[1] - p[0], n(-d[1], d[0], 0);
  return ImplicitConic{0, 0, 0, n[0], n[1], -n.dot(p[0])};
}

[[maybe_unused]]
auto implicitize(const EdgeCurve &c) {
  if (std::holds_alternative<Segment>(c))
    return implicitize(std::get<Segment>(c));
  return implicitize(std::get<Parabola>(c));
}

// Returns the (assumed to be only) real root in [0,1]
[[maybe_unused]]
double parabolaConicintersection(const Parabola &p, const ImplicitConic &c) {
  // Quadratic form: x(t) = ax t^2 + bx t + cx
  auto ax = p[0][0] - 2.0*p[1][0] + p[2][0];
  auto ay = p[0][1] - 2.0*p[1][1] + p[2][1];
  auto bx = -2.0*p[0][0] + 2.0*p[1][0];
  auto by = -2.0*p[0][1] + 2.0*p[1][1];
  auto cx = p[0][0];
  auto cy = p[0][1];

  // Precompute products for x^2, xy, y^2 expansions
  auto x2_4 = ax*ax, x2_3 = 2.0*ax*bx, x2_2 = 2.0*ax*cx + bx*bx, x2_1 = 2.0*bx*cx, x2_0 = cx*cx;
  auto y2_4 = ay*ay, y2_3 = 2.0*ay*by, y2_2 = 2.0*ay*cy + by*by, y2_1 = 2.0*by*cy, y2_0 = cy*cy;
  auto xy_4 = ax*ay, xy_3 = ax*by + bx*ay, xy_2 = ax*cy + bx*by + cx*ay, xy_1 = bx*cy + cx*by,
    xy_0 = cx*cy;
  auto x_2 = ax, x_1 = bx, x_0 = cx, y_2 = ay, y_1 = by, y_0 = cy;

  // Assemble coefficients
  std::array<double, 5> coeffs;
  coeffs[4] = c.A*x2_4 + c.B*xy_4 + c.C*y2_4;
  coeffs[3] = c.A*x2_3 + c.B*xy_3 + c.C*y2_3;
  coeffs[2] = c.A*x2_2 + c.B*xy_2 + c.C*y2_2 + c.D*x_2 + c.E*y_2;
  coeffs[1] = c.A*x2_1 + c.B*xy_1 + c.C*y2_1 + c.D*x_1 + c.E*y_1;
  coeffs[0] = c.A*x2_0 + c.B*xy_0 + c.C*y2_0 + c.D*x_0 + c.E*y_0 + c.F;

  return select01Root(generalSolver<4>(coeffs));
}

EdgeCurve fitParabola(const SimpleDomain::Edge &e) {
  // Do a LSQ fit on the central control point.
  // If it is close to the p0 p2 line,
  // or if it curves in the wrong direction, return a straight line.
  auto p0 = e.front(), p2 = e.back();
  // Centripetal parameterization
  std::vector<double> t;
  t.push_back(0);
  for (size_t i = 1; i < e.size(); ++i) {
    auto d = std::sqrt((e[i] - e[i-1]).norm());
    t.push_back(t.back() + d);
  }
  for (size_t i = 1; i < e.size(); ++i)
    t[i] /= t.back();
  double AtA = 0, Atx = 0, Aty = 0;
  for (size_t i = 0; i < e.size(); ++i) {
    auto t1 = (1 - t[i]) * (1 - t[i]), t2 = t[i] * t[i];
    double A = 2 * t[i] * (1 - t[i]);
    AtA += A * A;
    Atx += A * (e[i][0] - t1 * p0[0] - t2 * p2[0]);
    Aty += A * (e[i][1] - t1 * p0[1] - t2 * p2[1]);
  }
  Vec3 p1(Atx / AtA, Aty / AtA, 0);
  Vec3 d = (p2 - p0).normalized(), n(-d[1], d[0], 0);
  double signed_distance = (p1 - p0).dot(n);
  if (signed_distance < (p2 - p0).norm() * 0.1)
    return Segment{p0, p2};
  return Parabola{p0, p1, p2};
}

SimpleDomain::Circle fitCircle(const SimpleDomain::Loop &loop) {
  Vec3 center(0, 0, 0);
  size_t k = 0;
  for (const auto &edge : loop) {
    k += edge.size() - 1;
    for (size_t i = 1; i < edge.size(); ++i)
      center += edge[i];
  }
  center /= k;
  double r = 0;
  for (const auto &edge : loop)
    for (size_t i = 1; i < edge.size(); ++i)
      r += (edge[i] - center).norm();
  r /= k;
  Vec3 dev = (loop[0][0] - center).normalized();
  double offset = std::atan2(dev[1], dev[0]);
  return { center, r, offset, k };
}

Vec3 evalEdgeCurve(const EdgeCurve &c, double u) {
  if (std::holds_alternative<Segment>(c)) {
    const auto &p = std::get<Segment>(c);
    return p[0] * (1 - u) + p[1] * u;
  }
  const auto &p = std::get<Parabola>(c);
  return p[0] * (1 - u) * (1 - u) + p[1] * 2 * (1 - u) * u + p[2] * u * u;
}

} // anonymous namespace

void SimpleDomain::init(std::vector<SimpleDomain::Loop> &loops) {
  auto num_loops = loops.size();
  num_sides.clear();
  num_points.clear();
  boundaries.clear();
  holes.resize(num_loops);
  num_sides.push_back(loops[0].size());
  for (auto &edge : loops[0]) {
    auto c = fitParabola(edge);
    size_t n = edge.size() - 1;
    for (size_t i = 1; i < n; i++) {
      double u = (double)i / n;
      edge[i] = evalEdgeCurve(c, u);
    }
    num_points.push_back(edge.size());
    boundaries.push_back(c);
  }
  for (size_t loop = 1; loop < num_loops; ++loop) {
    num_sides.push_back(loops[loop].size());
    auto c = fitCircle(loops[loop]);
    size_t index = 0;
    for (auto &edge : loops[loop]) {
      for (size_t i = 0; i < edge.size(); ++i) {
        double phi = c.start_angle - 2 * M_PI * index / c.n_pts;
        edge[i] = c.center + Vec3(std::cos(phi), std::sin(phi), 0) * c.radius;
        index++;
      }
      index--;
    }
    holes[loop] = c;
  }
}

namespace {

void computeDistances(const Vec3 &p,
                      const std::vector<EdgeCurve> &boundaries,
                      const std::vector<SimpleDomain::Circle> &circles,
                      std::vector<double> &d) {
  size_t n = boundaries.size();
  for (size_t i = 0; i < n; ++i)
    // d.push_back(parabolaThroughPoint(parabolize(boundaries[i]), p).first);
    d.push_back(implicitize(boundaries[i]).normeval(p[0], p[1]));
  for (size_t loop = 1; loop < circles.size(); ++loop) {
    const auto &c = circles[loop];
    d.push_back((c.center - p).norm() - c.radius);
  }
}

void computeBoundaryS(const Vec3 &p,
                      const std::vector<EdgeCurve> &boundaries,
                      std::vector<double> &s) {
  auto x = p[0], y = p[1];
  auto n = boundaries.size();
  for (size_t i = 0; i < n; ++i) {
    size_t i1 = (i + 1) % n, i_1 = (i + n - 1) % n;
    auto left_i = implicitize(boundaries[i_1]), right_i = implicitize(boundaries[i1]);
    auto d1 = left_i.eval(x, y), d2 = right_i.eval(x, y);
    auto s0 = d1 / (d1 + d2);
    auto initial_s = left_i + (left_i + right_i) * (-s0);
    auto final_s = parabolaConicintersection(parabolize(boundaries[i]), initial_s);
    s.push_back(final_s);
  }
}

void computeHoleS(const Vec3 &p, const SimpleDomain::Circle &c, size_t num_sides,
                  std::vector<double> &s) {
  Vec3 dev = (p - c.center).normalized();
  double angle = std::atan2(dev[1], dev[0]);
  double pp = 2 * M_PI;
  for (size_t j = 0; j < num_sides; ++j)
    s.push_back(std::fmod(c.start_angle - angle + 10 * pp, pp) / pp);
}

size_t pointOnSide(const Vec3 &p, const std::vector<double> &d, double tol) {
  size_t m = d.size();
  for (size_t i = 0; i < m; ++i)
    if (d[i] < tol)
      return i;
  return m;
}

void computeSideH(size_t on_side, const std::vector<double> &d,
                  const std::vector<size_t> &num_sides,
                  std::vector<std::vector<double>> &h) {
  auto n = num_sides[0];
  if (on_side < n) {
    for (size_t i = 0; i < n; ++i)
      h[0].push_back(i == on_side ? 0 : 1);
    for (size_t i = 1; i < h.size(); ++i)
      for (size_t j = 0; j < num_sides[i]; ++j)
        h[i].push_back(1);
    size_t i_1 = (on_side + n - 1) % n, i1 = (on_side + 1) % n;
    auto s_i = d[i_1] / (d[i_1] + d[i1]);
    h[0][i_1] = s_i;
    h[0][i1] = 1 - s_i;
  } else {
    for (size_t i = 0; i < n; ++i)
      h[0].push_back(1);
    for (size_t i = 1; i < h.size(); ++i)
      for (size_t j = 0; j < num_sides[i]; ++j)
        h[i].push_back(i == on_side - n + 1 ? 0 : 1);
  }
}

void fixS(size_t on_side, size_t n, size_t n_pts, std::vector<double> &s) {
  n_pts--;
  s[on_side] = std::round(s[on_side] * n_pts) / n_pts;
  if (s[on_side] == 0)
    s[(on_side+n-1)%n] = 1;
  else if (s[on_side] == 1)
    s[(on_side+1)%n] = 0;
}

void computeInteriorH(const std::vector<double> &d,
                      const std::vector<size_t> &num_sides,
                      std::vector<std::vector<double>> &h) {
  size_t m = d.size(), n = num_sides[0];
  std::vector<double> prods(m); // prods[i] = 1/(di-1 di) or 1/di^2
  double sum = 0;
  for (size_t i = 0; i < m; ++i) {
    if (i < n) {
      size_t i_1 = (i + n - 1) % n;
      prods[i] = 1 / (d[i_1] * d[i]);
    } else
      prods[i] = 1 / (d[i] * d[i]);
    sum += prods[i];
  }
  for (size_t i = 0; i < n; ++i)
    h[0].push_back(1 - (prods[i] + prods[(i+1)%n]) / sum);
  for (size_t i = n; i < m; ++i)
    for (size_t j = 0; j < num_sides[i-n+1]; ++j)
      h[i-n+1].push_back(std::sqrt(1 - prods[i] / sum));
}

} // anonymous namespace

void SimpleDomain::computeParameters(const Vec3 &p,
                                     std::vector<std::vector<double>> &s,
                                     std::vector<std::vector<double>> &h) const {
  // Compute Euclidean distances
  auto num_loops = num_sides.size();
  std::vector<double> d;
  computeDistances(p, boundaries, holes, d);
  computeBoundaryS(p, boundaries, s[0]);
  for (size_t i = 1; i < num_loops; ++i)
    computeHoleS(p, holes[i], num_sides[i], s[i]);
  size_t on_side = pointOnSide(p, d, 1e-4); // returns d.size() when not
  if (on_side < d.size()) {
    computeSideH(on_side, d, num_sides, h);
    if (on_side < num_sides[0])
      fixS(on_side, num_sides[0], num_points[on_side], s[0]);
  } else
    computeInteriorH(d, num_sides, h);
}

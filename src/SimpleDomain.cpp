#include "libcdgbs/SimpleDomain.hpp"

using Vec3 = SimpleDomain::Vec3;

static SimpleDomain::EdgeCurve fitParabola(const SimpleDomain::Edge &e) {
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
    return SimpleDomain::Segment{p0, p2};
  return SimpleDomain::Parabola{p0, p1, p2};
}

static SimpleDomain::Circle fitCircle(const SimpleDomain::Loop &loop) {
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

static Vec3 evalEdgeCurve(const SimpleDomain::EdgeCurve &c, double u) {
  if (std::holds_alternative<SimpleDomain::Segment>(c)) {
    const auto &p = std::get<SimpleDomain::Segment>(c);
    return p[0] * (1 - u) + p[1] * u;
  }
  const auto &p = std::get<SimpleDomain::Parabola>(c);
  return p[0] * (1 - u) * (1 - u) + p[1] * 2 * (1 - u) * u + p[2] * u * u;
}

void SimpleDomain::init(std::vector<SimpleDomain::Loop> &loops) {
  auto num_loops = loops.size();
  num_points.clear();
  boundaries.clear();
  holes.resize(num_loops);
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

static void computeDistances(const Vec3 &p,
                             const std::vector<SimpleDomain::EdgeCurve> &boundaries,
                             const std::vector<SimpleDomain::Circle> &circles,
                             std::vector<double> &d,
                             std::vector<std::pair<double,double>> &dmax) {
  size_t n = boundaries.size();
  for (size_t i = 0; i < n; ++i) {
    const auto &edge = std::get<SimpleDomain::Segment>(boundaries[i]);
    Vec3 q = edge[0], t = (edge[1] - q).normalized(), ndir = {-t[1], t[0], 0};
    const Vec3 &pl = std::get<SimpleDomain::Segment>(boundaries[(i+n-1)%n])[0];
    const Vec3 &pr = std::get<SimpleDomain::Segment>(boundaries[(i+1)%n])[1];
    d.push_back(std::abs(ndir.dot(p - q)));
    dmax.push_back({std::abs(ndir.dot(pl - q)), std::abs(ndir.dot(pr - q))});
  }
  for (size_t loop = 1; loop < circles.size(); ++loop) {
    const auto &c = circles[loop];
    d.push_back((c.center - p).norm() - c.radius);
  }
}

static void computeBoundaryS(const Vec3 &p, size_t n,
                             const std::vector<double> &d,
                             const std::vector<std::pair<double,double>> &dmax,
                             std::vector<double> &s) {
  for (size_t i = 0; i < n; ++i) {
    size_t i1 = (i + 1) % n, i_1 = (i + n - 1) % n;
    auto d_1 = d[i_1] / dmax[i_1].second, d1 = d[i1] / dmax[i1].first;
    s.push_back(d_1 / (d_1 + d1));
  }
}

static void computeHoleS(const Vec3 &p, const SimpleDomain::Circle &c, size_t num_sides,
                         std::vector<double> &s) {
  Vec3 dev = (p - c.center).normalized();
  double angle = std::atan2(dev[1], dev[0]);
  double pp = 2 * M_PI;
  for (size_t j = 0; j < num_sides; ++j)
    s.push_back(std::fmod(c.start_angle - angle + 10 * pp, pp) / pp);
}

static size_t pointOnSide(const Vec3 &p, const std::vector<double> &d, double tol) {
  size_t m = d.size();
  for (size_t i = 0; i < m; ++i)
    if (d[i] < tol)
      return i;
  return m;
}

static void computeSideH(size_t on_side, const std::vector<double> &d,
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

static void fixS(size_t on_side, size_t n, size_t n_pts, std::vector<double> &s) {
  n_pts--;
  s[on_side] = std::round(s[on_side] * n_pts) / n_pts;
  if (s[on_side] == 0)
    s[(on_side+n-1)%n] = 1;
  else if (s[on_side] == 1)
    s[(on_side+1)%n] = 0;
}

static void computeInteriorH(const std::vector<double> &d,
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

void SimpleDomain::computeParameters(const Vec3 &p, const std::vector<size_t> &num_sides,
                                     std::vector<std::vector<double>> &s,
                                     std::vector<std::vector<double>> &h) const {
  // Compute Euclidean distances
  auto num_loops = num_sides.size();
  std::vector<double> d;
  std::vector<std::pair<double,double>> dmax; // prev & next maximal distances
  computeDistances(p, boundaries, holes, d, dmax);
  // Compute s
  computeBoundaryS(p, num_sides[0], d, dmax, s[0]);
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

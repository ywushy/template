#include <cstdio>
#include <cmath>
#include <vector>
#include <algorithm>

template<typename T>
class Point {
public:
    T x, y;
    Point() : x(0), y(0) {}
    Point(const T &x, const T &y) : x(x), y(y) {}
    Point(const Point<T> &p) : x(p.x), y(p.y) {}

    bool operator<(const Point<T> &p) const {
        return x < p.x || x == p.x && y < p.y;
    }
};

template<typename T> Point<T> operator+(const Point<T> &p1, const Point<T> &p2)	{
    return Point<T>(p1.x + p2.x, p1.y + p2.y);
}
template<typename T> Point<T> operator-(const Point<T> &p1, const Point<T> &p2) {
    return Point<T>(p1.x - p2.x, p1.y - p2.y);
}
template<typename T> Point<T> operator*(const T &k, const Point<T> &p) {
    return Point<T>(k * p.x, k * p.y);
}
template<typename T> Point<T> operator*(const Point<T> &p, const T &k) {
    return k * p;
}
template<typename T> Point<T> operator/(const Point<T> &p, const T &k) {
    return Point<T>(p.x / k, p.y / k);
}
template<typename T> T dot(const Point<T> &p1, const Point<T> &p2) {
    return p1.x * p2.x + p1.y * p2.y;
}
template<typename T> T cross(const Point<T> &p1, const Point<T> &p2) {
    return p1.x * p2.y - p1.y * p2.x;
}
template<typename T> T operator*(const Point<T> &p1, const Point<T> &p2) {
    return dot(p1, p2);
}
template<typename T> T operator^(const Point<T> &p1, const Point<T> &p2)	{
    return cross(p1, p2);
}
template<typename T> T norm2(const Point<T> &p) {
    return p * p;
}
template<typename T> T norm(const Point<T> &p) {
    return std::sqrt(norm2(p));
}
template<typename T> T distance2(const Point<T> &p1, const Point<T> &p2) {
    return norm2(p1 - p2);
}
template<typename T> T distance(const Point<T> &p1, const Point<T> &p2) {
    return norm(p1 - p2);
}
template<typename T> Point<T> rotate90(const Point<T> &p) {
    return Point<T>(p.y, -p.x);
}
template<typename T> Point<T> rotate180(const Point<T> &p) {
    return Point<T>(-p.x, -p.y);
}
template<typename T> Point<T> rotate270(const Point<T> &p) {
    return Point<T>(-p.y, p.x);
}
template<typename T> Point<T> normalise(const Point<T> &p) {
    return p / norm(p);
}

const double eps = 1e-8;

template<typename T>
class Circle {
public:
    Point<T> p;
    T r;

    Circle() : p(), r(0) {}
    Circle(const Point<T> &p, const T &r) : p(p), r(r) {}
    Circle(const T &x, const T &y, const T &r) : p(x, y), r(r) {}

    Circle(const Circle<T> &c) : p(c.p), r(c.r) {}
};

template<typename T>
class Line {
public:
    Point<T> p1, p2;

    Line() : p1(), p2() {}
    Line(const Point<T> &p1, const Point<T> &p2) : p1(p1), p2(p2) {}
    Line(const T &x1, const T &y1, const T &x2, const T &y2) : p1(x1, y1), p2(x2, y2) {}

    Line(const Line<T> &l) : p1(l.p1), p2(l.p2) {}

    // convert to form A * x + B * y + C = 0 without normalisation
    void toABC(T &A, T &B, T &C) const {
        A = p1.y - p2.y;
        B = p2.x - p1.x;
        C = p1 ^ p2;
    }
};

template<typename T> T distance(const Point<T> &p, const Line<T> &l) {
    return std::abs(normalise(l.p2 - l.p1) ^ (p - l.p1));
}
template<typename T> T distance(const Line<T> &l, const Point<T> &p) {
    return distance(p, l);
}

template<typename T>
Line<T> equidistantLine(const Point<T> &p1, const Point<T> &p2) {
    Point<T> o = (p1 + p2) / (T)2;
    return Line<T>(o, o + rotate90(p2 - p1));
}

template<typename T>
std::vector<Point<T> > intersection(const Line<T> &l1, const Line<T> &l2) {
    T A1, B1, C1, A2, B2, C2;

    l1.toABC(A1, B1, C1);
    l2.toABC(A2, B2, C2);

    T det = A1 * B2 - A2 * B1;

    if (std::abs(det) <= eps) return std::vector<Point<T> >();

    return std::vector<Point<T> >(1, Point<T>(C2 * B1 - C1 * B2, A2 * C1 - A1 * C2) / det);
}

template<typename T>
std::vector<Point<T> > intersection(const Circle<T> &c1, const Circle<T> &c2) {
    T dist2 = distance2(c1.p, c2.p);
    if (dist2 <= eps) return std::vector<Point<T> >();

    // check touching outside

    T rsum = c1.r + c2.r;
    T rsum2 = rsum * rsum;

    if (std::abs(dist2 - rsum2) <= eps) return std::vector<Point<T> >(1, (c1.p * c1.r + c2.p * c2.r) / rsum);
    if (dist2 > rsum2) return std::vector<Point<T> >();

    // check touching inside

    T rdiff = c1.r - c2.r;
    T rdiff2 = rdiff * rdiff;

    if (std::abs(dist2 - rdiff2) <= eps) return std::vector<Point<T> >(1, c1.p + (c2.p - c1.p) * c1.r / (c1.r - c2.r));
    if (dist2 < rdiff2) return std::vector<Point<T> >();

    // two points

    T dist = std::sqrt(dist2);

    T cosa = (c1.r * c1.r + dist2 - c2.r * c2.r) / 2 / c1.r / dist;
    T sina = std::sqrt(1 - cosa * cosa);

    Point<T> vec = (c2.p - c1.p) / dist;
    Point<T> o = c1.p + vec * c1.r * cosa;

    std::vector<Point<T> > res;

    res.push_back(o + rotate90(vec) * sina * c1.r);
    res.push_back(o + rotate270(vec) * sina * c1.r);

    return res;
}

template<typename T>
std::vector<Point<T> > intersection(const Circle<T> &c, const Line<T> &l) {
    // check touching outside

    Point<T> vec = normalise(l.p2 - l.p1);
    Point<T> o = l.p1 + vec * dot(vec, c.p - l.p1);

    T dist2 = distance2(o, c.p);

    if (std::abs(dist2 - c.r * c.r) <= eps) return std::vector<Point<T> >(1, o);
    if (dist2 > c.r * c.r) return std::vector<Point<T> >();

    std::vector<Point<T> > res;
    T len = std::sqrt(c.r * c.r - dist2);

    res.push_back(o + vec * len);
    res.push_back(o + rotate180(vec) * len);

    return res;
}

template<typename T>
std::vector<Point<T> > intersection(const Line<T> &l, const Circle<T> &c) {
    return intersection(c, l);
}

const double PI = 3.1415926535897932384626433832795;

Circle<double> c[3];

template<typename T>
Circle<T> equiangleCircle(const Circle<T> &cc1, const Circle<T> &cc2) {
    const Circle<T> &c1 = (cc1.r < cc2.r) ? cc1 : cc2;
    const Circle<T> &c2 = (cc1.r > cc2.r) ? cc1 : cc2;

    T dist = norm(c2.p - c1.p);
    Point<T> vec = (c2.p - c1.p) / dist;

//	printf("vec: %.8lf %.8lf\n", vec.x, vec.y);

    T din = c1.r * dist / (c2.r + c1.r);
    T dout = c1.r * dist / (c2.r - c1.r);

    Circle<T> res;
    res.r = (din + dout) / (T)2;
    res.p = c1.p - vec * (res.r - din);

    return res;
}

int main() {
    return 0;
}

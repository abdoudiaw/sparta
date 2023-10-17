#include <vector>
#include <string> // for std::string
#include <utility> // for std::pair
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

const std::string POINTS_FILENAME = "WALL/points.txt";
const std::string LINES_FILENAME = "WALL/lines.txt";

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point_2;
typedef std::vector<Point_2> PointVector;

namespace GeometryTools {

struct PointWall {
    int id;
    double x;
    double y;
};

struct Line {
    int start_id;
    int end_id;
};

}

extern std::vector<GeometryTools::PointWall> global_points;
extern std::vector<GeometryTools::Line> global_lines;

double distanceFromPointToSegment(const GeometryTools::PointWall& P, const GeometryTools::PointWall& A, const GeometryTools::PointWall& B);
K::Vector_2 getNormal(const GeometryTools::PointWall& A, const GeometryTools::PointWall& B);
double dotProduct(const K::Vector_2& A, const K::Vector_2& B);
double magnitude(const K::Vector_2& A);
double angleBetweenVectors(const K::Vector_3& a, const K::Vector_3& b);
std::tuple<double, double, double, K::Vector_3> computeAnglesWithNormal(const GeometryTools::PointWall &P, const K::Vector_3 &velocity, const K::Vector_3 &magneticField);

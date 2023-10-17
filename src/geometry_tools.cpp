#include "geometry_tools.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <stdexcept>

using namespace GeometryTools;

using GeometryTools::PointWall;
using GeometryTools::Line;


// Function to read points from a file and store them in a vector
std::vector<PointWall>  cachePointsFromFile() {
    std::vector<PointWall>  cachedPoints;
    std::ifstream file(POINTS_FILENAME);
    if (!file) {
        std::cerr << "Failed to open the file: " << POINTS_FILENAME << std::endl;
        return cachedPoints; // Return an empty vector on failure
    }
    PointWall point;
    while (file >> point.id >> point.x >> point.y) {
        cachedPoints.push_back(point);
    }
    return cachedPoints;
}

std::vector<Line> cacheLinesFromFile() {
    std::vector<Line> cachedLines;
    // const std::string filename = "lines.txt";
    std::ifstream file(LINES_FILENAME);
    if (!file) {
        std::cerr << "Failed to open the file: " << LINES_FILENAME << std::endl;
        return cachedLines;
    }
    Line line;
    while (file >> line.start_id >> line.end_id) {
        cachedLines.push_back(line);
    }
    return cachedLines;
}

std::vector<GeometryTools::PointWall> global_points = cachePointsFromFile();
std::vector<GeometryTools::Line> global_lines = cacheLinesFromFile();


// Function to calculate the distance from a point P to a line segment AB
double distanceFromPointToSegment(const PointWall& P, const PointWall& A, const PointWall& B) {
    K::Vector_2 a_to_b(B.x - A.x, B.y - A.y);
    K::Vector_2 a_to_p(P.x - A.x, P.y - A.y);

    double projection = a_to_p * a_to_b / (a_to_b.squared_length());

    if (projection < 0) return CGAL::squared_distance(Point_2(P.x, P.y), Point_2(A.x, A.y));
    if (projection > 1) return CGAL::squared_distance(Point_2(P.x, P.y), Point_2(B.x, B.y));
    
    Point_2 closest_point = Point_2(A.x, A.y) + projection * a_to_b;
    return CGAL::squared_distance(Point_2(P.x, P.y), closest_point);
}

// Function to calculate a normal vector from a line segment AB
K::Vector_2 getNormal(const PointWall& A, const PointWall& B) {
    return K::Vector_2(A.y - B.y, B.x - A.x);
}
double dotProduct(const K::Vector_2& A, const K::Vector_2& B) {
    return A.x() * B.x() + A.y() * B.y();
}

// Calculate the magnitude of a 3D vector
double magnitude(const K::Vector_2& A) {
    return std::sqrt(A.squared_length());
}

double clamp(double val, double min_val, double max_val) {
    return std::max(min_val, std::min(val, max_val));
}

double angleBetweenVectors(const K::Vector_3& a, const K::Vector_3& b) {
    double cosTheta = (a * b) / (sqrt(a.squared_length()) * sqrt(b.squared_length()));
    return std::acos(clamp(cosTheta, -1.0, 1.0)) * (180.0 / M_PI);
}

std::tuple<double, double, double, K::Vector_3> computeAnglesWithNormal(const PointWall &P, const K::Vector_3 &velocity, const K::Vector_3 &magneticField) {
    if (global_points.empty() || global_lines.empty()) {
        throw std::runtime_error("Error reading from files or files are empty.");
    }

    double min_distance = std::numeric_limits<double>::max();
    PointWall closest_point_start, closest_point_end;

    for (const Line& line : global_lines) {
        const PointWall& start_point = global_points[line.start_id - 1];
        const PointWall& end_point = global_points[line.end_id - 1];

        double current_distance = distanceFromPointToSegment(P, start_point, end_point);

        if (current_distance < min_distance) {
            min_distance = current_distance;
            closest_point_start = start_point;
            closest_point_end = end_point;
        }
    }

    K::Vector_2 tempNormal = getNormal(closest_point_start, closest_point_end);
    K::Vector_3 normal(tempNormal.x(), tempNormal.y(), 0.0);

    double angle_velocity_normal = angleBetweenVectors(normal, velocity);
    double angle_magneticField_normal = angleBetweenVectors(normal, magneticField);

    return std::tuple(angle_velocity_normal, angle_magneticField_normal, min_distance, normal);
}


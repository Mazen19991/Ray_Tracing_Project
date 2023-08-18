#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <Windows.h>

#include <cmath>
#include <stdio.h>
#include <random>

#define _USE_MATH_DEFINES
#include <math.h>


struct Vector3 {

	double x;
	double y;
	double z;

	Vector3() = default;

	Vector3(double x, double y, double z) {
		this->x = x;
		this->y = y;
		this->z = z;
	}

	Vector3(const Vector3& copy) {
		this->x = copy.x;
		this->y = copy.y;
		this->z = copy.z;
	}

	Vector3& operator +=(const Vector3& r) {
		this->x += r.x;
		this->y += r.y;
		this->z += r.z;

		return *this;
	}

	Vector3& operator -=(const Vector3& r) {
		this->x -= r.x;
		this->y -= r.y;
		this->z -= r.z;

		return *this;
	}

	Vector3& operator *=(double r) {
		this->x *= r;
		this->y *= r;
		this->z *= r;

		return *this;
	}

	Vector3& operator /=(double r) {
		this->x /= r;
		this->y /= r;
		this->z /= r;

		return *this;
	}

	double getIndex(int index) const {
		switch (index) {
		case 0:
			return this->x;
		case 1:
			return this->y;
		case 2:
			return this->z;
		default:
			printf("getIndex ERROR! Attempted out-of-bounds Vector3 access with index: %d\n", index);
			return 0.0;
		}
	}
	static Vector3 Random() {
		double x = rand() / (double)RAND_MAX;
		double y = rand() / (double)RAND_MAX;
		double z = rand() / (double)RAND_MAX;
		return Vector3(x, y, z).Normalized();
	}

	Vector3 Cross(const Vector3& r) const {
		return Vector3(
			this->y * r.z - this->z * r.y,
			this->z * r.x - this->x * r.z,
			this->x * r.y - this->y * r.x
		);
	}

	Vector3& setIndex(int index, double value) {
		switch (index) {
		case 0:
			this->x = value;
			break;
		case 1:
			this->y = value;
			break;
		case 2:
			this->z = value;
			break;
		default:
			printf("setIndex ERROR! Attempted out-of-bounds Vector3 access with index: %d\n", index);
		}
		return *this;
	}

	Vector3 MatrixMultiply(double m[3][3]) const {
		Vector3 result = Vector3(0, 0, 0);
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				double newVal = result.getIndex(i) + this->getIndex(j) * m[i][j];
				result.setIndex(i, newVal);
			}
		}

		return result;
	}

	double Dot(const Vector3& r) const {
		return this->x * r.x + this->y * r.y + this->z * r.z;
	}

	double Magnitude() const {
		return sqrt(this->x * this->x + this->y * this->y + this->z * this->z);
	}

	// Could be improved by using reciprocal sqrt with scalar multiplication
	Vector3 Normalized() {
		double m = this->Magnitude();
		return Vector3(this->x / m, this->y / m, this->z / m);
	}
};

double RandomInRange(double min, double max)
{
	return min + (max - min) * ((double)rand() / (double)RAND_MAX);
}

Vector3 RandomPointOnUnitSphere()
{
	// Generate random point in unit cube centered at origin
	Vector3 point(
		RandomInRange(-1.0, 1.0),
		RandomInRange(-1.0, 1.0),
		RandomInRange(-1.0, 1.0)
	);

	// If point is outside unit sphere, discard and generate new point
	if (point.Magnitude() > 1.0) {
		return RandomPointOnUnitSphere();
	}
	else {
		return point.Normalized();
	}
}

struct Sphere {
	Sphere() = default;
	Sphere(Vector3 center, double radius, Vector3 color, double specularExp, double reflection) {
		this->center = center;
		this->radius = radius;
		this->color = color;
		this->specularExp = specularExp;
		this->reflection = reflection;
	}
	Vector3 center;
	double radius;
	Vector3 color;
	double specularExp;
	double reflection;
};

struct Plane {
	Plane() = default;
	Plane(Vector3 normal, double distance, Vector3 color, double specularExp, double reflection) {
		this->normal = normal;
		this->distance = distance;
		this->color = color;
		this->specularExp = specularExp;
		this->reflection = reflection;
	}
	Vector3 normal;
	double distance;
	Vector3 color;
	double specularExp;
	double reflection;
};

struct Cylinder {
public:
	Cylinder() = default;
	Cylinder(Vector3 center, double radius, double height, Vector3 color, double specularExp, double reflection) {
		this->center = center;
		this->radius = radius;
		this->height = height;
		this->color = color;
		this->specularExp = specularExp;
		this->reflection = reflection;
	}

	Vector3 center;
	double radius;
	double height;
	Vector3 color;
	double specularExp;
	double reflection;
};

struct Cone {
public:
	Cone() = default;
	Cone(Vector3 center, double radius, double height, Vector3 color, double specularExp, double reflection) {
		this->center = center;
		this->radius = radius;
		this->height = height;
		this->color = color;
		this->specularExp = specularExp;
		this->reflection = reflection;
	}

	Vector3 center;
	double radius;
	double height;
	Vector3 color;
	double specularExp;
	double reflection;
};

struct Triangle {
	Triangle() = default;
	Triangle(Vector3 v1, Vector3 v2, Vector3 v3, Vector3 color, double specularExp, double reflection) {
		this->v1 = v1;
		this->v2 = v2;
		this->v3 = v3;
		this->color = color;
		this->specularExp = specularExp;
		this->reflection = reflection;
	}
	Vector3 v1;
	Vector3 v2;
	Vector3 v3;
	Vector3 color;
	double specularExp;
	double reflection;
};

enum LightType { AMBIENT_L, POINT_L, DIRECTIONAL_L };
struct Light {
	Light() = default;
	Light(LightType type, double intensity, Vector3 direction) {
		this->type = type;
		this->intensity = intensity;
		this->direction = direction;
	}

	// Convenience constructor for ambient light
	Light(LightType type, double intensity) {
		this->type = type;
		this->intensity = intensity;
		this->direction = Vector3(0, 0, 0);
	}
	LightType type;
	double intensity;
	Vector3 direction;
};

struct QFSolutions {
	QFSolutions() = default;
	QFSolutions(double t1, double t2) {
		this->t1 = t1;
		this->t2 = t2;
	}
	double t1;
	double t2;
};

struct RaySphereIntersection {
	RaySphereIntersection() = default;
	RaySphereIntersection(Sphere sphere, double t) {
		this->sphere = sphere;
		this->t = t;
	}
	Sphere sphere;
	double t;
};

struct RayPlaneIntersection {
	RayPlaneIntersection() = default;
	RayPlaneIntersection(Plane plane, double t) {
		this->plane = plane;
		this->t = t;
	}
	Plane plane;
	double t;
};

struct RayCylinderIntersection {
	RayCylinderIntersection() = default;
	RayCylinderIntersection(Cylinder cylinder, double t1, double t2) {
		this->cylinder = cylinder;
		this->t1 = t1;
		this->t2 = t2;
	}
	Cylinder cylinder;
	double t1;
	double t2;
};

struct RayConeIntersection {
	RayConeIntersection() = default;
	RayConeIntersection(Cone cone, double t) {
		this->cone = cone;
		this->t = t;
	}
	Cone cone;
	double t;
};

struct RayTriangleIntersection {
	RayTriangleIntersection() = default;
	RayTriangleIntersection(Triangle triangle, double t, Vector3 barycentricCoords) {
		this->triangle = triangle;
		this->t = t;
		this->barycentricCoords = barycentricCoords;
	}
	Triangle triangle;
	double t;
	Vector3 barycentricCoords;
};

// Negation
inline Vector3 operator -(const Vector3& r) {
	return Vector3(-r.x, -r.y, -r.z);
}

inline Vector3 operator +(const Vector3& l, const Vector3& r) {
	return Vector3(l.x + r.x, l.y + r.y, l.z + r.z);
}

inline Vector3 operator -(const Vector3& l, const Vector3& r) {
	return Vector3(l.x - r.x, l.y - r.y, l.z - r.z);
}

inline Vector3 operator *(double l, const Vector3& r) {
	return Vector3(l * r.x, l * r.y, l * r.z);
}

inline Vector3 operator *(const Vector3& l, double r) {
	return Vector3(l.x * r, l.y * r, l.z * r);
}

inline Vector3 operator /(const Vector3& l, double r) {
	return Vector3(l.x / r, l.y / r, l.z / r);
}

const LPCWSTR WNDCLASSNAME = L"GRAPHICS!";
const int WINDOW_WIDTH = 1000;
const int WINDOW_HEIGHT = 1000;
const int MAX_TRACE_RECURSIONS = 5;
bool quit = false;
bool draw = true;

HDC hDc;
HWND hWnd;

Sphere spheres[] = {
	Sphere(Vector3(-0.2, -0.1, 8), 0.2, Vector3(255, 0, 0), 60.0, 0.15),		 // Red
	Sphere(Vector3(4, 1.8, 10.0), 0.5, Vector3(255, 0, 0), 60.0, 0.15),		 // Red
	//Sphere(Vector3(2, 0, 5.0), 1.0, Vector3(0, 0, 255), 700.0, 0.35),		 // Blue
	//Sphere(Vector3(-2, 0.5, 4.5), 1.0, Vector3(0, 255, 0), 400.0, 0.3),		 // Green
	//Sphere(Vector3(0, -5001, 0), 5000, Vector3(255, 255, 0), 1000.0, 0.20),  // Yellow (floor)
	Sphere(Vector3(0, 1, 14.0), 1, Vector3(255, 102, 0), 0, 0)  // orange
};

Plane planes[] = {
	Plane(Vector3(0, 0, 1), 20.0, Vector3(255, 255, 255), 100000, 0.55),  // white
};

Cylinder cylinders[] = {
  Cylinder(Vector3(4, 3.0, 10.0), 1.0, 2.9, Vector3(0, 255, 0), 60.0, 0.15),  // green
  //Cylinder(Vector3(2, 0.5, 5.0), 0.7, 1.5, Vector3(0, 0, 255), 700.0, 0.35),  // Blue
  //Cylinder(Vector3(-2, 1.5, 4.5), 0.8, 3.0, Vector3(0, 255, 0), 400.0, 0.3),  // Green
  //Cylinder(Vector3(0, -5001, 0), 2.5, 5000, Vector3(255, 255, 0), 1000.0, 0.20),  // Yellow (floor)
  //Cylinder(Vector3(0, 3.0, 14.0), 1.5, 2.5, Vector3(255, 255, 255), 100.0, 0.55)  // White
};

Cone cones[] = {
	//Cone(Vector3(4, 3.0, 10.0), 1.0, 2.9, Vector3(0, 255, 0), 60.0, 0.15),  // Red
	Cone(Vector3(0, 3.0, 14.0), 4.0, 4.0, Vector3(0, 0, 255), 700.0, 0.35),  // Blue
	//Cone(Vector3(-2, 1.5, 4.5), 0.8, 3.0, Vector3(0, 255, 0), 400.0, 0.3),  // Green
	//Cone(Vector3(0, -5001, 0), 2.5, 5000, Vector3(255, 255, 0), 1000.0, 0.20),  // Yellow (floor)
	//Cone(Vector3(0, 3.0, 14.0), 1.5, 2.5, Vector3(255, 255, 255), 100.0, 0.55)  // White
};

Triangle triangles[] = {
	Triangle(Vector3(-3, -1, 9), Vector3(-1, -1, 9), Vector3(0, 0, 9), Vector3(255, 255, 0), 1000, 0.00001),  // yellow
	//Triangle(Vector3(1, 1, 0), Vector3(-1, -1, 0), Vector3(0, 0, 1), Vector3(0, 255, 0), 1000, 0.5)   // green
};

Light lights[] = {
	Light(LightType::AMBIENT_L, 0.2),
	Light(LightType::POINT_L, 0.6, Vector3(2, 1, 0)),
	Light(LightType::DIRECTIONAL_L, 0.25, Vector3(1, 4, 4))
};

// Converts (0,0)-is-top-left, increasing-y-goes-down coord to standard graph
// Necessary to make the viewport perspective correct
Vector3 WindowToCanvasCoord(int x, int y) {
	return Vector3(x - WINDOW_WIDTH / 2, -y + WINDOW_HEIGHT / 2, 0);
}

Vector3 CanvasToWindowCoord(int x, int y) {
	return Vector3(x + WINDOW_WIDTH / 2, -(y - WINDOW_HEIGHT / 2), 0);
}

Vector3 CanvasToViewportCoord(int canvasX, int canvasY, int vpWidth, int vpHeight) {
	double wRatio = (double)vpWidth / WINDOW_WIDTH;
	double hRatio = (double)vpHeight / WINDOW_HEIGHT;
	return Vector3(canvasX * wRatio, canvasY * hRatio, 1);
}

// INFINITY indicates no solution. t2 = INFINITY means one solution
QFSolutions SolveQuadraticEquation(double a, double b, double c) {
	double discriminant = (b * b) - (4 * a * c);
	if (discriminant < 0) {
		return QFSolutions(INFINITY, INFINITY);
	}

	if (discriminant == 0.0) {
		double t = -b / (2 * a);
		return QFSolutions(t, INFINITY);
	}

	double t1 = (-b + sqrt(discriminant)) / (2 * a);

	double t2 = (-b - sqrt(discriminant)) / (2 * a);

	return QFSolutions(t1, t2);
}

RaySphereIntersection ClosestSphereIntersection(const Vector3& origin, const Vector3& dir, double tMin, double tMax) {
	int numSpheres = sizeof(spheres) / sizeof(spheres[0]);
	double closestT = INFINITY;
	Sphere closestSphere = spheres[0];

	for (int i = 0; i < numSpheres; i++) {
		Sphere currSphere = spheres[i];
		Vector3 sphereToOrigin = origin - currSphere.center;

		double a = dir.Dot(dir);
		double b = 2 * sphereToOrigin.Dot(dir);
		double c = sphereToOrigin.Dot(sphereToOrigin) - (currSphere.radius * currSphere.radius);

		QFSolutions tVals = SolveQuadraticEquation(a, b, c);

		// Solution doesn't exist if it equals INFINITY
		if (tVals.t1 != INFINITY && tVals.t1 > tMin && tVals.t1 < tMax) {
			if (tVals.t1 < closestT) {
				closestT = tVals.t1;
				closestSphere = spheres[i];
			}
		}

		if (tVals.t2 != INFINITY && tVals.t2 > tMin && tVals.t2 < tMax) {
			if (tVals.t2 < closestT) {
				closestT = tVals.t2;
				closestSphere = spheres[i];
			}
		}
	}

	return RaySphereIntersection(closestSphere, closestT);
}

RayPlaneIntersection ClosestPlaneIntersection(const Vector3& origin, const Vector3& dir, double tMin, double tMax) {
	int numPlanes = sizeof(planes) / sizeof(planes[0]);
	double closestT = INFINITY;
	Plane closestPlane = planes[0];

	for (int i = 0; i < numPlanes; i++) {
		Plane currPlane = planes[i];
		double denom = currPlane.normal.Dot(dir);
		if (denom > 1e-6) {  // Check for non-parallel ray and plane
			double t = (currPlane.distance - origin.Dot(currPlane.normal)) / denom;
			if (t > tMin && t < tMax && t < closestT) {
				closestT = t;
				closestPlane = planes[i];
			}
		}
	}

	return RayPlaneIntersection(closestPlane, closestT);
}

RayCylinderIntersection ClosestCylinderIntersection(const Vector3& origin, const Vector3& dir, double tMin, double tMax) {
	int numCylinders = sizeof(cylinders) / sizeof(cylinders[0]);
	double closestT1 = INFINITY;
	double closestT2 = INFINITY;
	Cylinder closestCylinder = cylinders[0];

	for (int i = 0; i < numCylinders; i++) {
		Cylinder currCylinder = cylinders[i];
		Vector3 cylinderToOrigin = origin - currCylinder.center;

		// Calculate the intersection points of the ray and the cylinder
		double a = (dir.x * dir.x) + (dir.z * dir.z);
		double b = 2 * ((dir.x * cylinderToOrigin.x) + (dir.z * cylinderToOrigin.z));
		double c = (cylinderToOrigin.x * cylinderToOrigin.x) + (cylinderToOrigin.z * cylinderToOrigin.z) - (currCylinder.radius * currCylinder.radius);
		QFSolutions tVals = SolveQuadraticEquation(a, b, c);

		// Check if the intersection points are within the height of the cylinder
		if (tVals.t1 != INFINITY && tVals.t1 > tMin && tVals.t1 < tMax) {
			Vector3 intersectPoint = origin + (tVals.t1 * dir);
			double intersectY = intersectPoint.y;
			if (intersectY >= currCylinder.center.y - currCylinder.height / 2 && intersectY <= currCylinder.center.y + currCylinder.height / 2) {
				if (tVals.t1 < closestT1) {
					closestT1 = tVals.t1;
					closestCylinder = cylinders[i];
				}

			}
		}
		if (tVals.t2 != INFINITY && tVals.t2 > tMin && tVals.t2 < tMax) {
			Vector3 intersectPoint = origin + (tVals.t2 * dir);
			double intersectY = intersectPoint.y;
			if (intersectY >= currCylinder.center.y - currCylinder.height / 2 && intersectY <= currCylinder.center.y + currCylinder.height / 2) {
				if (tVals.t2 < closestT2) {
					closestT2 = tVals.t2;
					closestCylinder = cylinders[i];
				}
			}
		}
	}

	// Return the closest intersection point and the cylinder
	if (closestT1 < closestT2) {
		return RayCylinderIntersection(closestCylinder, closestT1, closestT2);
	}
	else {
		return RayCylinderIntersection(closestCylinder, closestT2, closestT1);
	}
}

RayConeIntersection ClosestConeIntersection(const Vector3& origin, const Vector3& dir, double tMin, double tMax) {
	int numCones = sizeof(cones) / sizeof(cones[0]);
	double closestT = INFINITY;
	Cone closestCone = cones[0];

	for (int i = 0; i < numCones; i++) {
		Cone currCone = cones[i];
		Vector3 coneToOrigin = origin - currCone.center;

		//Calculate the intersection points of the ray and the cone
		double a = (dir.x * dir.x) + (dir.z * dir.z) - ((dir.y * dir.y) * (currCone.radius * currCone.radius) / (currCone.height * currCone.height));
		double b = 2 * ((dir.x * coneToOrigin.x) + (dir.z * coneToOrigin.z) - ((dir.y * coneToOrigin.y) * (currCone.radius * currCone.radius) / (currCone.height * currCone.height)));
		double c = (coneToOrigin.x * coneToOrigin.x) + (coneToOrigin.z * coneToOrigin.z) - ((coneToOrigin.y * coneToOrigin.y) * (currCone.radius * currCone.radius) / (currCone.height * currCone.height));
		QFSolutions tVals = SolveQuadraticEquation(a, b, c);

		// Check if the intersection points are within the height of the cone
		if (tVals.t1 != INFINITY && tVals.t1 > tMin && tVals.t1 < tMax) {
			Vector3 intersectPoint = origin + (tVals.t1 * dir);
			double intersectY = intersectPoint.y;
			if (intersectY >= closestCone.center.y - closestCone.height / 2 && intersectY <= closestCone.center.y + closestCone.height / 2) {
				closestT = tVals.t1;
			}
		}
		if (tVals.t2 != INFINITY && tVals.t2 > tMin && tVals.t2 < tMax && tVals.t2 < closestT) {
			Vector3 intersectPoint = origin + (tVals.t2 * dir);
			double intersectY = intersectPoint.y;
			if (intersectY >= closestCone.center.y - closestCone.height / 2 && intersectY <= closestCone.center.y + closestCone.height / 2) {
				closestT = tVals.t2;
			}
		}
	}

	// Return the closest intersection point and the cone
	return RayConeIntersection(closestCone, closestT);
}

RayTriangleIntersection ClosestTriangleIntersection(const Vector3& origin, const Vector3& dir, double tMin, double tMax) {
	int numTriangles = sizeof(triangles) / sizeof(triangles[0]);
	double closestT = INFINITY;
	Triangle closestTriangle = triangles[0];
	Vector3 closestBarycentricCoords = Vector3(0, 0, 0);

	for (int i = 0; i < numTriangles; i++) {
		Triangle currTriangle = triangles[i];
		Vector3 e1 = currTriangle.v2 - currTriangle.v1;
		Vector3 e2 = currTriangle.v3 - currTriangle.v1;
		Vector3 s = origin - currTriangle.v1;
		Vector3 s1 = dir.Cross(e2);
		Vector3 s2 = s.Cross(e1);
		double denom = s1.Dot(e1);
		if (fabs(denom) > 1e-6) {  // Check for non-parallel ray and triangle
			double t = s2.Dot(e2) / denom;
			if (t > tMin && t < tMax && t < closestT) {
				Vector3 barycentricCoords = Vector3(s1.Dot(s) / denom, s2.Dot(dir) / denom, 1 - (s1.Dot(s) + s2.Dot(dir)) / denom);
				if (barycentricCoords.x >= 0 && barycentricCoords.y >= 0 && barycentricCoords.z >= 0) {  // Check if intersection is inside triangle
					closestT = t;
					closestTriangle = triangles[i];
					closestBarycentricCoords = barycentricCoords;
				}
			}
		}
	}

	return RayTriangleIntersection(closestTriangle, closestT, closestBarycentricCoords);
}


Vector3 ReflectAcrossNormal(const Vector3& ray, const Vector3& normal) {
	return 2 * normal * normal.Dot(ray) - ray;
}

double CalculateLightAtPoint(
	const Vector3& point,
	const Vector3& normal,
	const Vector3& camera,
	double specularExp)
{
	double totalIntensity = 0.0;
	int numLights = sizeof(lights) / sizeof(lights[0]);
	for (int i = 0; i < numLights; i++) {
		Light currLight = lights[i];
		if (currLight.type == LightType::AMBIENT_L) {
			totalIntensity += currLight.intensity;
		}
		else {
			// Calculate the diffusion for this point
			Vector3 pointToLight;
			int tMax;
			if (currLight.type == LightType::POINT_L) {
				pointToLight = currLight.direction - point;
				tMax = 1.0;
			}
			else if (currLight.type == LightType::DIRECTIONAL_L) {
				pointToLight = currLight.direction;
				tMax = INFINITY;
			}
			else {
				printf("ERROR! Unknown light type: %d", currLight.type);
				return 1.0; // Don't change existing light in case of error
			}

			// Check if this light source is obstructed by a sphere, a cylinder, a plane, or a cone
			RaySphereIntersection blocker1 = ClosestSphereIntersection(point, pointToLight, 0.0001, tMax);
			RayCylinderIntersection blocker2 = ClosestCylinderIntersection(point, pointToLight, 0.0001, tMax);
			RayPlaneIntersection blocker3 = ClosestPlaneIntersection(point, pointToLight, 0.0001, tMax);
			RayConeIntersection blocker4 = ClosestConeIntersection(point, pointToLight, 0.0001, tMax);
			RayTriangleIntersection blocker5 = ClosestTriangleIntersection(point, pointToLight, 0.0001, tMax);


			if (blocker1.t != INFINITY && blocker2.t1 != INFINITY && blocker2.t2 != INFINITY && blocker3.t != INFINITY && blocker4.t != INFINITY && blocker5.t != INFINITY) {
				double minT = min(min(blocker1.t, blocker2.t1), min(blocker2.t2, blocker3.t), min(blocker3.t, blocker4.t), min(blocker4.t, blocker5.t));
				if (blocker1.t == minT || blocker2.t1 == minT || blocker2.t2 == minT || blocker3.t == minT || blocker4.t == minT || blocker5.t == minT) {
					break;
				}
			}
			else if (blocker1.t != INFINITY || blocker2.t1 != INFINITY || blocker2.t2 != INFINITY || blocker3.t != INFINITY || blocker4.t != INFINITY || blocker5.t != INFINITY) {
				break;
			}
			double normalDotDir = normal.Dot(pointToLight);
			if (normalDotDir > 0) {
				double diffuse = currLight.intensity * normalDotDir
					/ (normal.Magnitude() * pointToLight.Magnitude());
				totalIntensity += diffuse;
			}
			if (specularExp > 0.0) {
				Vector3 pointToCamera = camera - point;
				Vector3 reflection = ReflectAcrossNormal(pointToLight, normal);
				if (reflection.Dot(pointToCamera) > 0) {
					double specular = currLight.intensity * pow(
						reflection.Dot(pointToCamera)
						/ (reflection.Magnitude() * pointToCamera.Magnitude()), specularExp);
					totalIntensity += specular;
				}
			}
		}
	}
	// Clamp at 1.0 because Windows will overflow RGB values over 255 with mod 255.
	return min(totalIntensity, 1.0);
}

//double CalculateLightAtPoint(
//	const Vector3& point,
//	const Vector3& normal,
//	const Vector3& camera,
//	double specularExp)
//{
//	double totalIntensity = 0.0;
//	int numLights = sizeof(lights) / sizeof(lights[0]);
//	for (int i = 0; i < numLights; i++) {
//		Light currLight = lights[i];
//		if (currLight.type == LightType::AMBIENT_L) {
//			totalIntensity += currLight.intensity;
//		}
//		else {
//			// Calculate the diffusion for this point
//			Vector3 pointToLight;
//			int tMax;
//			if (currLight.type == LightType::POINT_L) {
//				pointToLight = currLight.direction - point;
//				tMax = 1.0;
//			}
//			else if (currLight.type == LightType::DIRECTIONAL_L) {
//				pointToLight = currLight.direction;
//				tMax = INFINITY;
//			}
//			else {
//				printf("ERROR! Unknown light type: %d", currLight.type);
//				return 1.0; // Don't change existing light in case of error
//			}
//
//			// Check if this light source is obstructed by a sphere, a cylinder, a plane, or a cone
//			RaySphereIntersection blocker1 = ClosestSphereIntersection(point, pointToLight, 0.0001, tMax);
//			RayCylinderIntersection blocker2 = ClosestCylinderIntersection(point, pointToLight, 0.0001, tMax);
//			RayPlaneIntersection blocker3 = ClosestPlaneIntersection(point, pointToLight, 0.0001, tMax);
//			RayConeIntersection blocker4 = ClosestConeIntersection(point, pointToLight, 0.0001, tMax);
//			RayTriangleIntersection blocker5 = ClosestTriangleIntersection(point, pointToLight, 0.0001, tMax);
//
//			if (blocker1.t != INFINITY && blocker2.t1 != INFINITY && blocker2.t2 != INFINITY && blocker3.t != INFINITY && blocker4.t != INFINITY && blocker5.t != INFINITY) {
//				double minT = min(min(blocker1.t, blocker2.t1), min(blocker2.t2, blocker3.t), min(blocker3.t, blocker4.t), min(blocker4.t, blocker5.t));
//				if (blocker1.t == minT || blocker2.t1 == minT || blocker2.t2 == minT || blocker3.t == minT || blocker4.t == minT || blocker5.t == minT) {
//					break;
//				}
//			}
//			else if (blocker1.t != INFINITY || blocker2.t1 != INFINITY || blocker2.t2 != INFINITY || blocker3.t != INFINITY || blocker4.t != INFINITY || blocker5.t != INFINITY) {
//				break;
//			}
//
//			// Calculate soft shadows
//			int numSamples = 50;
//			double lightSourceArea = 0.0;
//			if (currLight.type == LightType::POINT_L) {
//				lightSourceArea = 4.0 * M_PI * currLight.radius * currLight.radius;
//			}
//			else if (currLight.type == LightType::DIRECTIONAL_L) {
//				// Assume light source is a rectangle with dimensions 2x2 centered at origin
//				lightSourceArea = 4.0;
//			}
//			double shadowIntensity = 0.0;
//			for (int j = 0; j < numSamples; j++) {
//				Vector3 samplePoint;
//				if (currLight.type == LightType::POINT_L) {
//					// Generate random point on light source sphere
//					samplePoint = currLight.direction + RandomPointOnUnitSphere() * currLight.radius;
//				}
//				else if (currLight.type == LightType::DIRECTIONAL_L) {
//					// Generate random point on light source rectangle
//					samplePoint = Vector3(RandomInRange(-1.0, 1.0), RandomInRange(-1.0, 1.0), 0.0);
//				}
//				Vector3 pointToSample = samplePoint - point;
//
//				if (blocker1.t == INFINITY && blocker2.t1 == INFINITY && blocker2.t2 == INFINITY && blocker3.t == INFINITY && blocker4.t == INFINITY && blocker5.t == INFINITY) {
//					shadowIntensity += 1.0;
//				}
//			}
//			shadowIntensity /= numSamples;
//
//			double normalDotDir = normal.Dot(pointToLight);
//			if (normalDotDir > 0) {
//				double diffuse = currLight.intensity * normalDotDir
//					/ (normal.Magnitude() * pointToLight.Magnitude()) * shadowIntensity;
//				totalIntensity += diffuse;
//			}
//			if (specularExp > 0.0) {
//				Vector3 pointToCamera = camera - point;
//				Vector3 reflection = ReflectAcrossNormal(pointToLight, normal);
//				if (reflection.Dot(pointToCamera) > 0) {
//					double specular = currLight.intensity * pow(
//						reflection.Dot(pointToCamera)
//						/ (reflection.Magnitude() * pointToCamera.Magnitude()),
//						specularExp
//					) * shadowIntensity;
//					totalIntensity += specular;
//				}
//			}
//		}
//	}
//	return totalIntensity;
//}

Vector3 TraceRay(const Vector3& origin, const Vector3& to, double tMin, double tMax, int recursions) {
	Vector3 defaultColor = Vector3(0, 0, 0);
	RaySphereIntersection closestSphere = ClosestSphereIntersection(origin, to, tMin, tMax);
	RayCylinderIntersection closestCylinder = ClosestCylinderIntersection(origin, to, tMin, tMax);
	RayConeIntersection closestCone = ClosestConeIntersection(origin, to, tMin, tMax);
	RayTriangleIntersection closestTriangle = ClosestTriangleIntersection(origin, to, tMin, tMax);
	RayPlaneIntersection closestPlane = ClosestPlaneIntersection(origin, to, tMin, tMax);
	if (closestSphere.t == INFINITY && closestCylinder.t1 == INFINITY && closestCylinder.t2 == INFINITY && closestCone.t == INFINITY && closestTriangle.t == INFINITY && closestPlane.t == INFINITY) {
		// Ray did not intersect any sphere, cylinder, or plane in the scene
		return defaultColor;
	}

	// Determine whether the ray intersects a sphere, cylinder, or plane first
	bool hitSphere = false;
	bool hitCylinder = false;
	bool hitCone = false;
	bool hitTriangle = false;
	Vector3 point;
	Vector3 intersectNormal;
	if (closestSphere.t < closestCylinder.t1 && closestSphere.t < closestCylinder.t2 && closestSphere.t < closestCone.t && closestSphere.t < closestTriangle.t && closestSphere.t < closestPlane.t) {
		hitSphere = true;
		point = origin + closestSphere.t * to; // Get intersect point from t
		intersectNormal = (point - closestSphere.sphere.center).Normalized();
	}
	else if (closestCylinder.t1 < closestCylinder.t2 && closestCylinder.t1 < closestCone.t && closestCylinder.t1 < closestTriangle.t && closestCylinder.t1 < closestPlane.t) {
		hitCylinder = true;
		point = origin + closestCylinder.t1 * to; // Get intersect point from t1
		intersectNormal = (point - closestCylinder.cylinder.center).Normalized();
	}
	else if (closestCylinder.t2 < closestCone.t && closestCylinder.t2 < closestTriangle.t && closestCylinder.t2 < closestPlane.t) {
		hitCylinder = true;
		point = origin + closestCylinder.t2 * to; // Get intersect point from t2
		intersectNormal = (point - closestCylinder.cylinder.center).Normalized();
	}
	else if (closestCone.t < closestTriangle.t && closestCone.t < closestPlane.t) {
		hitCone = true;
		point = origin + closestCone.t * to; // Get intersect point from t1
		intersectNormal = (point - closestCone.cone.center).Normalized();
	}
	else if (closestTriangle.t < closestPlane.t) {
		hitTriangle = true;
		point = origin + closestTriangle.t * to; // Get intersect point from t
		intersectNormal = (point - closestTriangle.barycentricCoords).Normalized();
	}
	else {
		point = origin + closestPlane.t * to;
		intersectNormal = closestPlane.plane.normal;
	}

	Vector3 localColor;
	if (hitSphere) {
		localColor = closestSphere.sphere.color *
			CalculateLightAtPoint(point, intersectNormal, -1 * to, closestSphere.sphere.specularExp);
	}
	else if (hitCylinder) {
		localColor = closestCylinder.cylinder.color *
			CalculateLightAtPoint(point, intersectNormal, -1 * to, closestCylinder.cylinder.specularExp);
	}
	else if (hitCone) {
		localColor = closestCone.cone.color *
			CalculateLightAtPoint(point, intersectNormal, -1 * to, closestCone.cone.specularExp);
	}
	else if (hitTriangle) {
		localColor = closestTriangle.triangle.color *
			CalculateLightAtPoint(point, intersectNormal, -1 * to, closestTriangle.triangle.specularExp);
	}
	else {
		localColor = closestPlane.plane.color *
			CalculateLightAtPoint(point, intersectNormal, -1 * to, closestPlane.plane.specularExp);
	}

	if (hitSphere && closestSphere.sphere.reflection <= 0.0 && hitCylinder && closestCylinder.cylinder.reflection <= 0.0 && hitCone && closestCone.cone.reflection <= 0.0 && hitTriangle && closestTriangle.triangle.reflection <= 0.0 && !hitSphere && !hitCylinder && !hitCone && !hitTriangle && closestPlane.plane.reflection <= 0.0 || recursions >= MAX_TRACE_RECURSIONS) {
		return localColor;
	}

	Vector3 reflectedColor = TraceRay(point, ReflectAcrossNormal(-to, intersectNormal), 0.01, tMax, recursions + 1);
	if (hitSphere) {
		return (localColor * (1 - closestSphere.sphere.reflection) + reflectedColor * closestSphere.sphere.reflection);
	}
	else if (hitCylinder) {
		return (localColor * (1 - closestCylinder.cylinder.reflection) + reflectedColor * closestCylinder.cylinder.reflection);
	}
	else if (hitCone) {
		return (localColor * (1 - closestCone.cone.reflection) + reflectedColor * closestCone.cone.reflection);
	}
	else if (hitTriangle) {
		return (localColor * (1 - closestTriangle.triangle.reflection) + reflectedColor * closestTriangle.triangle.reflection);
	}
	else {
		return (localColor * (1 - closestPlane.plane.reflection) + reflectedColor * closestPlane.plane.reflection);
	}
}

//Vector3 TraceRay(const Vector3& origin, const Vector3& to, double tMin, double tMax, int recursions) {
//	Vector3 defaultColor = Vector3(0, 0, 0);
//
//	RaySphereIntersection closestSphere = ClosestSphereIntersection(origin, to, tMin, tMax);
//	RayCylinderIntersection closestCylinder = ClosestCylinderIntersection(origin, to, tMin, tMax);
//	RayConeIntersection closestCone = ClosestConeIntersection(origin, to, tMin, tMax);
//	RayTriangleIntersection closestTriangle = ClosestTriangleIntersection(origin, to, tMin, tMax);
//	RayPlaneIntersection closestPlane = ClosestPlaneIntersection(origin, to, tMin, tMax);
//	if (closestSphere.t == INFINITY && closestCylinder.t1 == INFINITY && closestCylinder.t2 == INFINITY && closestCone.t == INFINITY && closestTriangle.t == INFINITY && closestPlane.t == INFINITY) {
//		// Ray did not intersect any sphere, cylinder, or plane in the scene
//		return defaultColor;
//	}
//
//	// Determine whether the ray intersects a sphere, cylinder, or plane first
//	bool hitSphere = false;
//	bool hitCylinder = false;
//	bool hitCone = false;
//	bool hitTriangle = false;
//	Vector3 point;
//	Vector3 intersectNormal;
//	if (closestSphere.t < closestCylinder.t1 && closestSphere.t < closestCylinder.t2 && closestSphere.t < closestCone.t && closestSphere.t < closestTriangle.t && closestSphere.t < closestPlane.t) {
//		hitSphere = true;
//		point = origin + closestSphere.t * to; // Get intersect point from t
//		intersectNormal = (point - closestSphere.sphere.center).Normalized();
//	}
//	else if (closestCylinder.t1 < closestCylinder.t2 && closestCylinder.t1 < closestCone.t && closestCylinder.t1 < closestTriangle.t && closestCylinder.t1 < closestPlane.t) {
//		hitCylinder = true;
//		point = origin + closestCylinder.t1 * to; // Get intersect point from t1
//		intersectNormal = (point - closestCylinder.cylinder.center).Normalized();
//	}
//	else if (closestCylinder.t2 < closestCone.t && closestCylinder.t1 < closestTriangle.t && closestCylinder.t2 < closestPlane.t) {
//		hitCylinder = true;
//		point = origin + closestCylinder.t2 * to; // Get intersect point from t2
//		intersectNormal = (point - closestCylinder.cylinder.center).Normalized();
//	}
//	else if (closestCone.t < closestTriangle.t && closestCone.t < closestPlane.t) {
//		hitCone = true;
//		point = origin + closestCone.t * to; // Get intersect point from t1
//		intersectNormal = (point - closestCone.cone.center).Normalized();
//	}
//	else if (closestTriangle.t < closestPlane.t) {
//		hitTriangle = true;
//		point = origin + closestTriangle.t * to; // Get intersect point from t1
//		intersectNormal = (point - closestTriangle.barycentricCoords).Normalized();
//	}
//	else {
//		point = origin + closestPlane.t * to;
//		intersectNormal = closestPlane.plane.normal;
//	}
//
//	// Naive anti-aliasing: average color over 4 rays
//	Vector3 finalColor = Vector3(0, 0, 0);
//	finalColor += TraceRay(origin, (to + Vector3(0.001, 0, 0)).Normalized(), tMin, tMax, recursions + 1);
//	finalColor += TraceRay(origin, (to + Vector3(-0.001, 0, 0)).Normalized(), tMin, tMax, recursions + 1);
//	finalColor += TraceRay(origin, (to + Vector3(0, 0.001, 0)).Normalized(), tMin, tMax, recursions + 1);
//	finalColor += TraceRay(origin, (to + Vector3(0, -0.001, 0)).Normalized(), tMin, tMax, recursions + 1);
//	finalColor /= 4;
//
//	Vector3 localColor;
//	if (hitSphere) {
//		localColor = closestSphere.sphere.color *
//			CalculateLightAtPoint(point, intersectNormal, -1 * to, closestSphere.sphere.specularExp);
//	}
//	else if (hitCylinder) {
//		localColor = closestCylinder.cylinder.color *
//			CalculateLightAtPoint(point, intersectNormal, -1 * to, closestCylinder.cylinder.specularExp);
//	}
//	else if (hitCone) {
//		localColor = closestCone.cone.color *
//			CalculateLightAtPoint(point, intersectNormal, -1 * to, closestCone.cone.specularExp);
//	}
//	else if (hitTriangle) {
//		localColor = closestTriangle.triangle.color *
//			CalculateLightAtPoint(point, intersectNormal, -1 * to, closestTriangle.triangle.specularExp);
//	}
//	else {
//		localColor = closestPlane.plane.color *
//			CalculateLightAtPoint(point, intersectNormal, -1 * to, closestPlane.plane.specularExp);
//	}
//	return finalColor * (1 - closestPlane.plane.reflection) + localColor * closestPlane.plane.reflection;
//}

// Processes all messages sent to the window.
LRESULT CALLBACK WinProc(HWND hWnd, UINT Msg, WPARAM wParam, LPARAM lParam) {
	switch (Msg) {
		// Sent when window is closed
	case WM_CLOSE: {
		PostQuitMessage(0);
		break;
	}
	}
	// Send everything else to default processor
	return DefWindowProc(hWnd, Msg, wParam, lParam);
}

// This is mostly just setting window params and then the main loop
int WINAPI WinMain(
	HINSTANCE hInstance,
	HINSTANCE hPrevInstance,
	LPSTR lpCmdLine,
	int nShowCmd)
{
	MSG msg;
	WNDCLASSEX ex;
	ex.cbSize = sizeof(WNDCLASSEX);
	ex.style = CS_VREDRAW | CS_HREDRAW | CS_OWNDC; // combine styles w/ bitwise OR
	ex.lpfnWndProc = WinProc;
	ex.cbClsExtra = 0;
	ex.cbWndExtra = 0;
	ex.hInstance = hInstance;
	ex.hIcon = LoadIcon(NULL, IDI_APPLICATION);
	ex.hCursor = LoadCursor(NULL, IDC_ARROW);
	ex.hbrBackground = (HBRUSH)GetStockObject(BLACK_BRUSH);
	ex.lpszMenuName = NULL;
	ex.lpszClassName = (LPCWSTR)WNDCLASSNAME;
	ex.hIconSm = NULL;
	RegisterClassEx(&ex);

	// Finally making the window
	hWnd = CreateWindowEx(
		NULL,
		(LPCWSTR)WNDCLASSNAME,
		L"Game Design 1 Project",
		WS_OVERLAPPEDWINDOW ^ WS_MAXIMIZEBOX ^ WS_MINIMIZEBOX,
		100, 100, WINDOW_WIDTH, WINDOW_HEIGHT, NULL, NULL, hInstance, NULL);
	ShowWindow(hWnd, SW_SHOW);
	UpdateWindow(hWnd);

	hDc = GetDC(hWnd);
	Vector3 cameraCoord = Vector3(1, 0, 0);

	// Message loop
	while (!quit) {
		if (draw) {
			for (int x = 0; x < WINDOW_WIDTH; x++) {
				for (int y = 0; y < WINDOW_HEIGHT; y++) {
					Vector3 canvasCoord = WindowToCanvasCoord(x, y);
					Vector3 vpCoord = CanvasToViewportCoord(canvasCoord.x, canvasCoord.y, 1, 1);
					Vector3 color = TraceRay(cameraCoord, vpCoord, 1.0, INFINITY, 3);
					SetPixel(hDc, x, y, RGB(color.x, color.y, color.z));
				}
			}
			draw = false;
		}
		if (PeekMessage(&msg, NULL, NULL, NULL, PM_REMOVE)) {
			quit = (msg.message == WM_QUIT);
			TranslateMessage(&msg);
			DispatchMessage(&msg);
		}
	}
	ReleaseDC(hWnd, hDc);

	return msg.lParam;
}
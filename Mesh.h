#ifndef MESH_H
#define MESH_H
#include <iostream>
#include <vector>
#include <algorithm>

#include <Eigen/Dense>

using namespace std;

#define MachineEpsilon (std::numeric_limits<float>::epsilon() * 0.5)
inline float gamma(int n) {
	return (n * MachineEpsilon) / (1 - n * MachineEpsilon);
}

// Ray Declarations
struct Ray {
	// Ray Public Methods
	Ray():tMax(INFINITY){}
	Ray(const Eigen::Vector3f &o, const Eigen::Vector3f &d)
		: o(o), d(d), tMax(INFINITY){}
	Eigen::Vector3f operator()(float t) const { return o + d * t; }
	friend std::ostream &operator<<(std::ostream &os, const Ray &r) {
		os << "[o=" << r.o << ", d=" << r.d << "]";
		return os;
	}

	// Ray Public Data
	Eigen::Vector3f o;
	Eigen::Vector3f d;
	float tMax;
};

struct Bounds3f {
	Eigen::Vector3f pMin, pMax;

	Bounds3f() : pMax(Eigen::Vector3f(0,0,0)), pMin(Eigen::Vector3f(0, 0, 0)) {}
	Bounds3f(const Eigen::Vector3f p) : pMin(p), pMax(p) {}

	Bounds3f(const Eigen::Vector3f &pMin, const Eigen::Vector3f &pMax)
		: pMin(pMin)
		, pMax(pMax) {}

	Bounds3f unionBounds(const Bounds3f &b) {
		if (pMax == Eigen::Vector3f(0,0,0) && pMin == Eigen::Vector3f(0,0,0))
			return b;
		Eigen::Vector3f x(
			std::min(b.pMin[0], pMin[0]),
			std::min(b.pMin[1], pMin[1]),
			std::min(b.pMin[2], pMin[2]));
		Eigen::Vector3f y(
			std::max(b.pMax[0], pMax[0]),
			std::max(b.pMax[1], pMax[1]),
			std::max(b.pMax[2], pMax[2]));
		return Bounds3f(x,y);
	}

	int MaximumExtent() const {
		Eigen::Vector3f d = pMax - pMin;
		if (d[0] > d[1] && d[0] > d[2])
			return 0;
		else if (d[1] > d[2])
			return 1;
		else
			return 2;
	}

	Eigen::Vector3f Offset(const Eigen::Vector3f& p) const {
		Eigen::Vector3f o = p - pMin;
		if (pMax[0] > pMin[0]) o[0] /= pMax[0] - pMin[0];
		if (pMax[1] > pMin[1]) o[1] /= pMax[1] - pMin[1];
		if (pMax[2] > pMin[2]) o[2] /= pMax[2] - pMin[2];
		return o;
	}

	Eigen::Vector3f getCentroid() {
		return (pMin + pMax) * 0.5f;
	}

	float SurfaceArea() const {
		Eigen::Vector3f d = pMax - pMin;
		return 2 * (d[0] * d[1] + d[0] * d[2] + d[1] * d[2]);
	}

	bool IntersectP(const Ray &ray, const Eigen::Vector3f &invDir,
		const int dirIsNeg[3]) const {
		// Check for ray intersection against $x$ and $y$ slabs
		float tMin = ((dirIsNeg[0]? pMax : pMin)[0] - ray.o[0]) * invDir[0];
		float tMax = ((dirIsNeg[0] ? pMin : pMax)[0] - ray.o[0]) * invDir[0];
		float tyMin = ((dirIsNeg[1] ? pMin : pMax)[1] - ray.o[1]) * invDir[1];
		float tyMax = ((dirIsNeg[1] ? pMin : pMax)[1] - ray.o[1]) * invDir[1];

		// Update _tMax_ and _tyMax_ to ensure robust bounds intersection
		tMax *= 1 + 2 * gamma(3);
		tyMax *= 1 + 2 * gamma(3);
		if (tMin > tyMax || tyMin > tMax) return false;
		if (tyMin > tMin) tMin = tyMin;
		if (tyMax < tMax) tMax = tyMax;

		// Check for ray intersection against $z$ slab
		float tzMin = ((dirIsNeg[2] ? pMin : pMax)[2] - ray.o[2]) * invDir[2];
		float tzMax = ((dirIsNeg[2] ? pMin : pMax)[2] - ray.o[2]) * invDir[2];

		// Update _tzMax_ to ensure robust bounds intersection
		tzMax *= 1 + 2 * gamma(3);
		if (tMin > tzMax || tzMin > tMax) return false;
		if (tzMin > tMin) tMin = tzMin;
		if (tzMax < tMax) tMax = tzMax;
		return (tMin < ray.tMax) && (tMax > 0);
	}

	bool IntersectP(const Ray &ray, float *hitt0,
		float *hitt1) const {
		float t0 = 0, t1 = INFINITY;
		for (int i = 0; i < 3; ++i) {
			// Update interval for _i_th bounding box slab
			float invRayDir = 1 / ray.d[i];
			float tNear = (pMin[i] - ray.o[i]) * invRayDir;
			float tFar = (pMax[i] - ray.o[i]) * invRayDir;

			// Update parametric interval from slab intersection $t$ values
			if (tNear > tFar) std::swap(tNear, tFar);

			// Update _tFar_ to ensure robust ray--bounds intersection
			tFar *= 1 + 2 * gamma(3);
			t0 = tNear > t0 ? tNear : t0;
			t1 = tFar < t1 ? tFar : t1;
			if (t0 > t1) return false;
		}
		if (hitt0) *hitt0 = t0;
		if (hitt1) *hitt1 = t1;
		return true;
	}
};

struct SurfaceInteraction {
	Eigen::Vector3f pHit;
	float t;
	SurfaceInteraction(Eigen::Vector3f pHit, float t):pHit(pHit),t(t){}
};

int maxDimension(Eigen::Vector3f v) {
	return (v[0] > v[1]) ? ((v[0] > v[2]) ? 0 : 2) : ((v[1] > v[2]) ? 1 : 2);
}
int maxComponent(Eigen::Vector3f v) {
	return std::max(std::max(v[0],v[1]),v[2]);
}
Eigen::Vector3f abs(Eigen::Vector3f v) {
	return Eigen::Vector3f(std::abs(v[0]), std::abs(v[1]), std::abs(v[2]));
}

class Mesh {
public:
	/*  Mesh Data  */
	vector<Eigen::Vector3f> points;

	/*  Functions  */
	// constructor
	Mesh(const std::vector<Eigen::Vector3f> &vPoints) {
        Eigen::Vector3f pmin(INFINITY, INFINITY, INFINITY);
        Eigen::Vector3f pmax(-INFINITY, -INFINITY, -INFINITY);
        for (int i = 0; i < vPoints.size(); ++i) {
			pmin[0] = std::min(pmin[0],vPoints[i][0]);
			pmin[1] = std::min(pmin[1],vPoints[i][1]);
			pmin[2] = std::min(pmin[2],vPoints[i][2]);
			pmax[0] = std::max(pmax[0],vPoints[i][0]);
			pmax[1] = std::max(pmax[1],vPoints[i][1]);
			pmax[2] = std::max(pmax[2],vPoints[i][2]);
        }
        bounds=Bounds3f();
        bounds.pMin= pmin;
        bounds.pMax= pmax;

		points = vPoints;
    }

	Bounds3f getBounds() { return this->bounds; }
	Eigen::Vector3f getCentroid() { return this->bounds.getCentroid(); }

	bool Intersect(Ray &ray, SurfaceInteraction *isect) const {
		for(int i=0;i<points.size()/3;++i) {
			// Get triangle vertices in _p0_, _p1_, and _p2_
			const Eigen::Vector3f &p0 = points[i * 3 + 0];
			const Eigen::Vector3f &p1 = points[i * 3 + 1];
			const Eigen::Vector3f &p2 = points[i * 3 + 2];

			// Perform ray--triangle intersection test

			// Transform triangle vertices to ray coordinate space

			// Translate vertices based on ray origin
			Eigen::Vector3f p0t = p0 - Eigen::Vector3f(ray.o);
			Eigen::Vector3f p1t = p1 - Eigen::Vector3f(ray.o);
			Eigen::Vector3f p2t = p2 - Eigen::Vector3f(ray.o);

			// Permute components of triangle vertices and ray direction
			int kz = maxDimension(abs(ray.d));
			int kx = kz + 1;
			if (kx == 3) kx = 0;
			int ky = kx + 1;
			if (ky == 3) ky = 0;
			Eigen::Vector3f d = Eigen::Vector3f(ray.d[kx], ray.d[ky], ray.d[kz]);
			p0t = Eigen::Vector3f(p0t[kx], p0t[ky], p0t[kz]);
			p1t = Eigen::Vector3f(p1t[kx], p1t[ky], p1t[kz]);
			p2t = Eigen::Vector3f(p2t[kx], p2t[ky], p2t[kz]);

			// Apply shear transformation to translated vertex positions
			float Sx = -d[0] / d[2];
			float Sy = -d[1] / d[2];
			float Sz = 1.f / d[2];
			p0t[0] += Sx * p0t[2];
			p0t[1] += Sy * p0t[2];
			p1t[0] += Sx * p1t[2];
			p1t[1] += Sy * p1t[2];
			p2t[0] += Sx * p2t[2];
			p2t[1] += Sy * p2t[2];

			// Compute edge function coefficients _e0_, _e1_, and _e2_
			float e0 = p1t[0] * p2t[1] - p1t[1] * p2t[0];
			float e1 = p2t[0] * p0t[1] - p2t[1] * p0t[0];
			float e2 = p0t[0] * p1t[1] - p0t[1] * p1t[0];

			// Fall back to double precision test at triangle edges
			if (sizeof(float) == sizeof(float) &&
				(e0 == 0.0f || e1 == 0.0f || e2 == 0.0f)) {
				double p2txp1ty = (double)p2t[0] * (double)p1t[1];
				double p2typ1tx = (double)p2t[1] * (double)p1t[0];
				e0 = (float)(p2typ1tx - p2txp1ty);
				double p0txp2ty = (double)p0t[0] * (double)p2t[1];
				double p0typ2tx = (double)p0t[1] * (double)p2t[0];
				e1 = (float)(p0typ2tx - p0txp2ty);
				double p1txp0ty = (double)p1t[0] * (double)p0t[1];
				double p1typ0tx = (double)p1t[1] * (double)p0t[0];
				e2 = (float)(p1typ0tx - p1txp0ty);
			}

			// Perform triangle edge and determinant tests
			if ((e0 < 0 || e1 < 0 || e2 < 0) && (e0 > 0 || e1 > 0 || e2 > 0))
				return false;
			float det = e0 + e1 + e2;
			if (det == 0) return false;

			// Compute scaled hit distance to triangle and test against ray $t$ range
			p0t[2] *= Sz;
			p1t[2] *= Sz;
			p2t[2] *= Sz;
			float tScaled = e0 * p0t[2] + e1 * p1t[2] + e2 * p2t[2];
			if (det < 0 && (tScaled >= 0 || tScaled < ray.tMax * det))
				return false;
			else if (det > 0 && (tScaled <= 0 || tScaled > ray.tMax * det))
				return false;

			// Compute barycentric coordinates and $t$ value for triangle intersection
			float invDet = 1 / det;
			float b0 = e0 * invDet;
			float b1 = e1 * invDet;
			float b2 = e2 * invDet;
			float t = tScaled * invDet;

			// Ensure that computed triangle $t$ is conservatively greater than zero

			// Compute $\delta_z$ term for triangle $t$ error bounds
			float maxZt = maxComponent(abs(Eigen::Vector3f(p0t[2], p1t[2], p2t[2])));
			float deltaZ = gamma(3) * maxZt;

			// Compute $\delta_x$ and $\delta_y$ terms for triangle $t$ error bounds
			float maxXt = maxComponent(abs(Eigen::Vector3f(p0t[0], p1t[0], p2t[0])));
			float maxYt = maxComponent(abs(Eigen::Vector3f(p0t[1], p1t[1], p2t[1])));
			float deltaX = gamma(5) * (maxXt + maxZt);
			float deltaY = gamma(5) * (maxYt + maxZt);

			// Compute $\delta_e$ term for triangle $t$ error bounds
			float deltaE =
				2 * (gamma(2) * maxXt * maxYt + deltaY * maxXt + deltaX * maxYt);

			// Compute $\delta_t$ term for triangle $t$ error bounds and check _t_
			float maxE = maxComponent(abs(Eigen::Vector3f(e0, e1, e2)));
			float deltaT = 3 *
				(gamma(3) * maxE * maxZt + deltaE * maxZt + deltaZ * maxE) *
				std::abs(invDet);
			if (t <= deltaT) return false;

			// Interpolate $(u,v)$ parametric coordinates and hit point
			Eigen::Vector3f pHit = b0 * p0 + b1 * p1 + b2 * p2;

			// Fill in _SurfaceInteraction_ from triangle hit
			*isect = SurfaceInteraction(pHit, t);

			ray.tMax = t;
			
		}
		
		return true;
	}

	bool IntersectP(const Ray &ray) const {
		// Get triangle vertices in _p0_, _p1_, and _p2_
		for (int i = 0; i < points.size() / 3; ++i) {
			// Get triangle vertices in _p0_, _p1_, and _p2_
			const Eigen::Vector3f &p0 = points[i * 3 + 0];
			const Eigen::Vector3f &p1 = points[i * 3 + 1];
			const Eigen::Vector3f &p2 = points[i * 3 + 2];

			// Perform ray--triangle intersection test

			// Transform triangle vertices to ray coordinate space

			// Translate vertices based on ray origin
			Eigen::Vector3f p0t = p0 - Eigen::Vector3f(ray.o);
			Eigen::Vector3f p1t = p1 - Eigen::Vector3f(ray.o);
			Eigen::Vector3f p2t = p2 - Eigen::Vector3f(ray.o);

			// Permute components of triangle vertices and ray direction
			int kz = maxDimension(abs(ray.d));
			int kx = kz + 1;
			if (kx == 3) kx = 0;
			int ky = kx + 1;
			if (ky == 3) ky = 0;
			Eigen::Vector3f d = Eigen::Vector3f(ray.d[kx], ray.d[ky], ray.d[kz]);
			p0t = Eigen::Vector3f(p0t[kx], p0t[ky], p0t[kz]);
			p1t = Eigen::Vector3f(p1t[kx], p1t[ky], p1t[kz]);
			p2t = Eigen::Vector3f(p2t[kx], p2t[ky], p2t[kz]);

			// Apply shear transformation to translated vertex positions
			float Sx = -d[0] / d[2];
			float Sy = -d[1] / d[2];
			float Sz = 1.f / d[2];
			p0t[0] += Sx * p0t[2];
			p0t[1] += Sy * p0t[2];
			p1t[0] += Sx * p1t[2];
			p1t[1] += Sy * p1t[2];
			p2t[0] += Sx * p2t[2];
			p2t[1] += Sy * p2t[2];

			// Compute edge function coefficients _e0_, _e1_, and _e2_
			float e0 = p1t[0] * p2t[1] - p1t[1] * p2t[0];
			float e1 = p2t[0] * p0t[1] - p2t[1] * p0t[0];
			float e2 = p0t[0] * p1t[1] - p0t[1] * p1t[0];

			// Fall back to double precision test at triangle edges
			if (sizeof(float) == sizeof(float) &&
				(e0 == 0.0f || e1 == 0.0f || e2 == 0.0f)) {
				double p2txp1ty = (double)p2t[0] * (double)p1t[1];
				double p2typ1tx = (double)p2t[1] * (double)p1t[0];
				e0 = (float)(p2typ1tx - p2txp1ty);
				double p0txp2ty = (double)p0t[0] * (double)p2t[1];
				double p0typ2tx = (double)p0t[1] * (double)p2t[0];
				e1 = (float)(p0typ2tx - p0txp2ty);
				double p1txp0ty = (double)p1t[0] * (double)p0t[1];
				double p1typ0tx = (double)p1t[1] * (double)p0t[0];
				e2 = (float)(p1typ0tx - p1txp0ty);
			}

			// Perform triangle edge and determinant tests
			if ((e0 < 0 || e1 < 0 || e2 < 0) && (e0 > 0 || e1 > 0 || e2 > 0))
				return false;
			float det = e0 + e1 + e2;
			if (det == 0) return false;

			// Compute scaled hit distance to triangle and test against ray $t$ range
			p0t[2] *= Sz;
			p1t[2] *= Sz;
			p2t[2] *= Sz;
			float tScaled = e0 * p0t[2] + e1 * p1t[2] + e2 * p2t[2];
			if (det < 0 && (tScaled >= 0 || tScaled < ray.tMax * det))
				return false;
			else if (det > 0 && (tScaled <= 0 || tScaled > ray.tMax * det))
				return false;

			// Compute barycentric coordinates and $t$ value for triangle intersection
			float invDet = 1 / det;
			float b0 = e0 * invDet;
			float b1 = e1 * invDet;
			float b2 = e2 * invDet;
			float t = tScaled * invDet;

			// Ensure that computed triangle $t$ is conservatively greater than zero

			// Compute $\delta_z$ term for triangle $t$ error bounds
			float maxZt = maxComponent(abs(Eigen::Vector3f(p0t[2], p1t[2], p2t[2])));
			float deltaZ = gamma(3) * maxZt;

			// Compute $\delta_x$ and $\delta_y$ terms for triangle $t$ error bounds
			float maxXt = maxComponent(abs(Eigen::Vector3f(p0t[0], p1t[0], p2t[0])));
			float maxYt = maxComponent(abs(Eigen::Vector3f(p0t[1], p1t[1], p2t[1])));
			float deltaX = gamma(5) * (maxXt + maxZt);
			float deltaY = gamma(5) * (maxYt + maxZt);

			// Compute $\delta_e$ term for triangle $t$ error bounds
			float deltaE =
				2 * (gamma(2) * maxXt * maxYt + deltaY * maxXt + deltaX * maxYt);

			// Compute $\delta_t$ term for triangle $t$ error bounds and check _t_
			float maxE = maxComponent(abs(Eigen::Vector3f(e0, e1, e2)));
			float deltaT = 3 *
				(gamma(3) * maxE * maxZt + deltaE * maxZt + deltaZ * maxE) *
				std::abs(invDet);
			if (t <= deltaT) return false;

			return true;
		}
	}

private:
	/*  Render data  */
	Bounds3f bounds;
	/*  Functions    */

};

#endif

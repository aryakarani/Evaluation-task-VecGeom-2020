#ifndef _POINT_H
#define _POINT_H

#include "vector.h"

struct point : public vector
{
	point() = default;
	point(float x, float y, float z) : vector(x, y, z) {}
	explicit point(const vector& v) : point(v.x, v.y, v.z) {}
};

inline point operator+(const point& p, const vector& v)
{ return point(p.x + v.x, p.y + v.y, p.z + v.z); }

inline point operator-(const point& p, const vector& v)
{ return point(p.x - v.x, p.y - v.y, p.z - v.z); }

inline vector operator-(const point& p1, const point& p2)
{ return vector(p1.x - p2.x, p1.y - p2.y, p1.z - p2.z); }

inline float distance(const point& a, const point& b)
{ return norm(a - b); }

inline float distance2(const point& a, const point& b)
{ return norm2(a - b); }

#endif

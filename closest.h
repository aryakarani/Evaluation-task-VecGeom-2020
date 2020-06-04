#ifndef _CLOSEST_H_
#define _CLOSEST_H_

struct point;

/* calculate closest point in triangle ABC to the origin */
point closest_point_triangle(const point &A, const point &B, const point &C);

/* calculate closest point in tetrahedron ABCD to the origin */
point closest_point_tetrahedron(const point &A, const point &B, const point &C, const point &D);

#endif

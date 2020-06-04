#include "point.h"
#include "vector.h"
#include "closest.h"

#include <cstdio>

int main(int argc, char **argv)
{
	point triangle[3] = {
		{-0.5f, -0.5f, 0.0f },
		{ 0.0f,  0.5f, 0.0f },
		{ 0.5f, -0.5f, 0.0f }
	};

	point tetrahedron[4] = {
		{ 0.0f, 0.0f, 0.0f },
		{ 1.0f, 0.0f, 0.0f },
		{ 0.0f, 1.0f, 0.0f },
		{ 0.0f, 0.0f, 1.0f },
	};

	float coords[] = {-2.0f, -1.0f, -0.5f, 0.0f, +0.5f, +1.0f, +2.0f};

	for (auto z : coords) {
		for (auto y : coords) {
			for (auto x : coords) {
				point a, b, c, d;
				vector v(x, y, z);

				a = triangle[0] - v;
				b = triangle[1] - v;
				c = triangle[2] - v;

				point p_tri = closest_point_triangle(a, b, c) + v;

				a = tetrahedron[0] - v;
				b = tetrahedron[1] - v;
				c = tetrahedron[2] - v;
				d = tetrahedron[3] - v;

				point p_tetra = closest_point_tetrahedron(a, b, c, d) + v;

				printf("[% f, % f, % f], ", v.x, v.y, v.z);
				printf("tri[% f, % f, % f] (%f), ",
						p_tri.x, p_tri.y, p_tri.z, norm(p_tri - v));
				printf("tetra[% f, % f, % f] (%f)\n",
						p_tetra.x, p_tetra.y, p_tetra.z, norm(p_tetra - v));
			}
		}
	}

	return 0;
}

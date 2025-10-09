/*
 * Copyright (c) 2025, Jonathan Debove
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include <math.h>

#include "utils.h"

/*
 * 2D Laplace equation
 * d2p/dx2 + d2p/dy2 = b
 *
 * p(i,j) = ((p(i+1,j)+p(i-1,j))* dy2 + (p(i,j+1)+p(i,j-1))* dx2 - dx2*dy2*b(i,j)) / (2*(dx2+dy2))
 *
 * IC: p=0 with 0<=x<=2 and 0<=y<=1
 *
 * BC: p(0,y)=p(2,y)=p(x,0)=p(x,1)=0
 *
 * b(0.5,0.25)=100, b(1.5,0.75)=-100, b=0 elsewhere
 */

int main(void)
{
	int nx = 20;
	int ny = 10;
	int nt = 500;
	double dx = 2.0 / (nx - 1);
	double dy = 1.0 / (ny - 1);

	int n = 0;
	int i = 0;
	int j = 0;

	struct mat p0, p1, b;
	mat_create(&p0, nx, ny);
	mat_create(&p1, nx, ny);
	mat_create(&b, nx, ny);

	*mat_at(b, (1*nx)/4, (1*ny)/4) = +100.0;
	*mat_at(b, (3*nx)/4, (3*ny)/4) = -100.0;

	mat_print(p0, "p", 0);

	double dx2 = dx * dx;
	double dy2 = dy * dy;
	double den = (dx2 + dy2) * 2.0;
	for (n = 1; n <= nt; n++) {
		// Equation
		for (j = 1; j < ny - 1; j++) {
			for (i = 1; i < nx - 1; i++) {
				*mat_at(p1, i, j) =
					(dy2*(*mat_at(p0,i+1,j) +
					      *mat_at(p0,i-1,j)) +
					 dx2*(*mat_at(p0,i,j+1) +
					      *mat_at(p0,i,j-1)) -
					 dx2*dy2*(*mat_at(b,i,j))) / den;
			}
		}
		// BC
		for (i = 0; i < nx; i++) {
			*mat_at(p1, i, 0) = 0.0;
			*mat_at(p1, i, ny-1) = 0.0;
		}
		for (j = 0; j < ny; j++) {
			*mat_at(p1, 0, j) = 0.0;
			*mat_at(p1, nx-1, j) = 0.0;
		}

		if (mat_dist(p0, p1) < 1e-6) {
			mat_print(p1, "p", n);
			break;
		}
		if (n % 100 == 0) mat_print(p1, "p", n);

		SWAP(struct mat, p0, p1);
	}

	mat_destroy(&p0);
	mat_destroy(&p1);

	return 0;
}

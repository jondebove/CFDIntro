/*
 * Copyright (c) 2025, Jonathan Debove
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "utils.h"

/*
 * 2D Diffusion
 * du/dt = nu*lapla(u)
 *
 * FD in time, CD in space
 * u(n+1,i,j) = u(n,i,j) + nu*dt/dx*(u(n,i+1,j)-2*u(n,i,j)+u(n,i-1,j))
 *                       + nu*dt/dy*(u(n,i,j+1)-2*u(n,i,j)+u(n,i,j-1))
 *
 * IC: u(0,0.5<=x,y<=1) = 2, u(0,elsewhere) = 1
 *
 * BC: u(n,0)=u(n,2)=1
 */

int main(void)
{
	int nx = 20;
	int ny = 20;
	int nt = 50;
	double dx = 2.0 / (nx - 1);
	double dy = 2.0 / (ny - 1);
	double dt = 0.01;
	double nu = 0.1;

	int n = 0;
	int i = 0;
	int j = 0;

	struct mat u0, u1;
	mat_create(&u0, nx, ny);
	mat_create(&u1, nx, ny);

	// IC
	for (i = 0; i < nx; i++) {
		double x = dx * i;
		for (j = 0; j < ny; j++) {
			double y = dy * j;
			*mat_at(u0, i, j) = x >= 0.5 && x <= 1.0 &&
				y >= 0.5 && y <= 1.0 ?  2.0 : 1.0;
		}
	}
	mat_print(u0, "u", 0);

	for (n = 1; n <= nt; n++) {
		// BC
		for (i = 0; i < nx; i++) {
			*mat_at(u1, i, 0) = 1.0;
			*mat_at(u1, i, ny-1) = 1.0;
		}
		for (j = 0; j < ny; j++) {
			*mat_at(u1, 0, j) = 1.0;
			*mat_at(u1, nx-1, j) = 1.0;
		}

		// Equation
		for (i = 1; i < nx - 1; i++) {
			for (j = 1; j < ny - 1; j++) {
				*mat_at(u1, i, j) = *mat_at(u0, i, j)
					+ nu*dt/dx *
					(*mat_at(u0, i+1, j) -
					 *mat_at(u0, i, j)*2 +
					 *mat_at(u0, i-1, j))
					+ nu*dt/dy *
					(*mat_at(u0, i, j+1) -
					 *mat_at(u0, i, j)*2 +
					 *mat_at(u0, i, j-1));
			}
		}
		mat_print(u1, "u", n);

		SWAP(struct mat, u0, u1);
	}

	mat_destroy(&u0);
	mat_destroy(&u1);

	return 0;
}

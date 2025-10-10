/*
 * Copyright (c) 2025, Jonathan Debove
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "utils.h"

/*
 * Diffusion equation
 * du/dt = nu*d2u/dx2
 *
 * Solution: u(t, x) = A*exp(i*k*x)*exp(-nu*k²*t)
 *
 * FD in time, CD in space
 * u(n+1,i) = u(n,i) + nu*dt/dx2*(u(n,i+)-2*u(n,i)+u(n,i-1))
 *
 * IC: u(0,0.5<=x<=1) = 2, u(0,elsewhere) = 1
 *
 * BC: u(n,0)=u(n,2)=1
 */

int main(void)
{
	int nx = 20;
	int nt = 50;
	double nu = 0.1;
	double dx = 2.0 / (nx - 1);
	double dt = 0.01;

	double nudtdx2 = nu * dt / (dx * dx);

	int n = 0;
	int i = 0;

	struct mat u0, u1;
	mat_create(&u0, nx, 1);
	mat_create(&u1, nx, 1);

	// IC
	for (i = 0; i < nx; i++) {
		double x = dx * i;
		*mat_at(u0, i, 0) = x >= 0.5 && x <= 1.0 ? 2.0 : 1.0;
	}
	mat_print(u0, "u", 0);

	for (n = 1; n <= nt; n++) {
		// BC
		*mat_at(u1, 0, 0) = 1.0;
		*mat_at(u1, nx-1, 0) = 1.0;

		// Equation
		for (i = 1; i < nx - 1; i++) {
			*mat_at(u1, i, 0) = *mat_at(u0, i, 0) + nudtdx2 *
				(*mat_at(u0, i+1, 0) -
				 *mat_at(u0, i, 0)*2 +
				 *mat_at(u0, i-1, 0));
		}
		mat_print(u1, "u", n);

		SWAP(struct mat, u0, u1);
	}

	mat_destroy(&u0);
	mat_destroy(&u1);

	return 0;
}

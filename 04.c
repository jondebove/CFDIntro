/*
 * Copyright (c) 2025, Jonathan Debove
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include <math.h>

#ifndef M_PI
#define M_PI 3.141592653589793
#endif

#include "utils.h"

/*
 * Burgers'equation
 * du/dt = nu*d2u/dx2
 *
 * FD in time, BD and CD in space
 * u(n+1,i) = u(n,i) - dt/dx*u(n,i)*(u(n,i)-u(n,i-1)
 *                   + nu*dt/dx2*(u(n,i+)-2*u(n,i)+u(n,i-1))
 *
 * IC: u(0, x) = -2*nu*1/phi*dphi/dx + 4
 *     with phi(x) = exp(-x²/(4*nu))+exp((x-2*pi)²/(4*nu))
 *
 * BC: u(n,0) = u(n,2*pi)
 */

static
double sol(double t, double x, double nu)
{
	t *= 4.0;
	x -= t;
	t += 4.0;
	double x1 = x - M_PI * 2.0;
	double phi0 = exp(-x*x/(nu*t));
	double phi1 = exp(-x1*x1/(nu*t));
	double dphi0 = 4.0/t*x*phi0;
	double dphi1 = 4.0/t*x1*phi1;
	return 4.0 + (dphi0 + dphi1) / (phi0 + phi1);
}

int main(void)
{
	int nx = 20;
	int nt = 50;
	double nu = 0.1;
	double dx = 2.0 * M_PI / (nx - 1);
	double dt = 0.01;

	double nudtdx2 = nu * dt / (dx * dx);

	int err = EXIT_SUCCESS;
	int n = 0;
	int i = 0;

	struct mat u0, u1;
	mat_create(&u0, nx, 1);
	mat_create(&u1, nx, 1);

	// IC
	for (i = 0; i < nx; i++) {
		double x = dx * i;
		*mat_at(u0, i, 0) = sol(0, x, nu);
	}
	mat_print(u0, 0);

	for (n = 1; n <= nt; n++) {
		// Equation
		for (i = 0; i < nx; i++) {
			double udtdx = *mat_at(u0, i, 0) * dt / dx;
			if (udtdx > 1.0) {
				WARN("CFL=%f > 1\n", udtdx);
			}
			int im = i > 0 ? i - 1 : nx - 1;
			int ip = i < nx - 1 ? i + 1 : 0;
			*mat_at(u1, i, 0) = *mat_at(u0, i, 0) - udtdx *
					(*mat_at(u0, i, 0) -
					 *mat_at(u0, im, 0)) +
					nudtdx2 *
					(*mat_at(u0, ip, 0) -
					 *mat_at(u0, i, 0) * 2 +
					 *mat_at(u0, im, 0));
		}
		mat_print(u1, n);

		SWAP(struct mat, u0, u1);

		/* Print solution with n < 0 */
		for (i = 0; i < nx; i++) {
			*mat_at(u1, i, 0) = sol(n * dt, i * dx, nu);
		}
		mat_print(u1, -n);
	}

	mat_destroy(&u0);
	mat_destroy(&u1);

	return err;
}

/*
 * Copyright (c) 2025, Jonathan Debove
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include <assert.h>
#include <math.h>

#include "utils.h"

/*
 * 2D Laplace equation
 * d2p/dx2 + d2p/dy2 = 0
 *
 * p(i,j) = ((p(i+1,j)+p(i-1,j))* dy2 + (p(i,j+1)+p(i,j-1))* dx2) / (2*(dx2+dy2))
 *
 * IC: p = 0 with 0<=x<=2 and 0<=y<=1
 *
 * BC: p(0,y)=0, p(2,y)=y, dp/dy=0 @ y=0,1
 *
 * Solution: p(x,y) = x/4 - 4*sum_{n odd} sinh(n*pi*x)*cos(n*pi*y)/
 *                                        ((n*pi)²*sinh(2*pi*n))
 */

static
double sol(double x, double y, double tol, int maxit)
{
	if (x == 0.0) return 0;
	else if (x == 2.0) return y;

	maxit *= 2;
	int n;
	double ans = x / 4.0;
	for (n = 1; n < maxit; n += 2) {
		double npi = n * M_PI;
		double tmp = 4.0 * cos(npi*y)/(npi*npi)*exp(npi*(x-2.0))*
			expm1(-2.0*x*npi)/expm1(-4.0*npi);
		ans -= tmp;
		if (fabs(tmp) <= tol * fabs(ans)) {
			return ans;
		}
	}
	WARNF("maxit=%d reached\n", n / 2);
	return ans;
}

static
double poisson_iter(struct mat p, double dx2, double dy2, double omega)
{
	assert(p.nrow > 0);
	assert(p.ncol > 0);
	assert(dx2 > 0);
	assert(dy2 > 0);

	double den = (dx2 + dy2) * 2.0;
	dx2 /= den;
	dy2 /= den;

	int i, j, k;
	double epsmax = 0.0;
	/* Red-black sweep */
	for (k = 0; k <= 1; k++) {
		for (j = 1; j < p.ncol - 1; j++) {
			for (i = 2 - (j + k) % 2; i < p.nrow - 1; i += 2) {
				double eps = omega * (
						(*mat_at(p,i+1,j) +
						 *mat_at(p,i-1,j)) * dy2 -
						(*mat_at(p,i,j)) +
						(*mat_at(p,i,j+1) +
						 *mat_at(p,i,j-1)) * dx2);
				double p1 = *mat_at(p, i, j) + eps;

				eps /= MAX(ABS(p1), ABS(*mat_at(p, i, j)));
				epsmax = MAX(epsmax, eps);

				*mat_at(p, i, j) = p1;
			}
		}
	}
	return epsmax;
}

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

	struct mat p;
	mat_create(&p, nx, ny);

	// BC
	for (j = 0; j < ny; j++) {
		*mat_at(p, 0, j) = 0.0;
		*mat_at(p, nx-1, j) = j*dy;
	}
	mat_print(p, "p", 0);

	double dx2 = dx * dx;
	double dy2 = dy * dy;

	double omega = 1.0;
	double rho = 1.0 - 1.0 / (nx * ny);

	for (n = 1; n <= nt; n++) {
		// Equation
		double eps = poisson_iter(p, dx2, dy2, omega);
		// BC
		for (i = 1; i < nx-1; i++) {
			*mat_at(p, i, 0) = *mat_at(p, i, 1);
			*mat_at(p, i, ny-1) = *mat_at(p, i, ny-2);
		}

		if (eps < 1e-8) break;
		if (n % 10 == 0) mat_print(p, "p", n);
		omega = n == 1 ? 2.0 / (2.0 - rho) : 4.0 / (4.0 - rho*omega);
	}
	mat_print(p, "p", n);
	WARNF("n=%d\n", n);

	// Analytical solution
	for (j = 0; j < ny; j++) {
		for (i = 0; i < nx; i++) {
			*mat_at(p, i, j) = sol(i*dx, j*dy, 1e-8, 100);
		}
	}
	mat_print(p, "p", -1);

	mat_destroy(&p);

	return 0;
}

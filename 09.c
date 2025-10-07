/*
 * Copyright (c) 2025, Jonathan Debove
 * SPDX-License-Identifier: BSD-2-Clause
 */

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

int main(void)
{
	int nx = 20;
	int ny = 10;
	int nt = 5000;
	double dx = 2.0 / (nx - 1);
	double dy = 1.0 / (ny - 1);

	int n = 0;
	int i = 0;
	int j = 0;

	struct mat p0, p1;
	mat_create(&p0, nx, ny);
	mat_create(&p1, nx, ny);

	// BC
	for (j = 0; j < ny; j++) {
		*mat_at(p0, 0, j) = 0.0;
		*mat_at(p0, nx-1, j) = j*dy;
	}
	mat_print(p0, "p", 0);

	double dx2 = dx * dx;
	double dy2 = dy * dy;
	double den = (dx2 + dy2) * 2.0;
	for (n = 1; n <= nt; n++) {
		// Equation
		for (j = 1; j < ny - 1; j++) {
			for (i = 1; i < nx - 1; i++) {
				*mat_at(p1, i, j) = (dy2*(*mat_at(p0,i+1,j)+
							*mat_at(p0,i-1,j)) +
						dx2*(*mat_at(p0,i,j+1)+
							*mat_at(p0,i,j-1))) / den;
			}
		}
		// BC
		for (i = 0; i < nx; i++) {
			*mat_at(p1, i, 0) = *mat_at(p1, i, 1);
			*mat_at(p1, i, ny-1) = *mat_at(p1, i, ny-2);
		}
		for (j = 0; j < ny; j++) {
			*mat_at(p1, 0, j) = 0.0;
			*mat_at(p1, nx-1, j) = j*dy;
		}

		if (mat_dist(p0, p1) < 1e-6) {
			mat_print(p1, "p", n);
			break;
		}
		if (n % 100 == 0) mat_print(p1, "p", n);

		SWAP(struct mat, p0, p1);
	}

	// Analytical solution
	for (j = 0; j < ny; j++) {
		for (i = 0; i < nx; i++) {
			*mat_at(p0, i, j) = sol(i*dx, j*dy, 1e-8, 100);
		}
	}
	mat_print(p0, "p", -1);

	mat_destroy(&p0);
	mat_destroy(&p1);

	return 0;
}

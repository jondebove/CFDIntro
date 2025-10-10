/*
 * Copyright (c) 2025, Jonathan Debove
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include <assert.h>
#include <math.h>

#include "utils.h"

/*
 * 2D Poisson equation
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

static
double poisson_iter(struct mat p, struct mat b,
		double dx2, double dy2, double omega)
{
	assert(p.nrow > 0);
	assert(p.ncol > 0);
	assert(p.nrow == b.nrow);
	assert(p.ncol == b.ncol);
	assert(dx2 > 0);
	assert(dy2 > 0);

	double den = (dx2 + dy2) * 2.0;
	double dxdy2 = dx2 * dy2 / den;
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
						 *mat_at(p,i,j-1)) * dx2 -
						(*mat_at(b,i,j)) * dxdy2);
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

	double dx2 = 2.0 / (nx - 1);
	double dy2 = 1.0 / (ny - 1);
	dx2 *= dx2;
	dy2 *= dy2;

	struct mat p, b;
	mat_create(&p, nx, ny);
	mat_create(&b, nx, ny);

	*mat_at(b, (1*nx)/4, (1*ny)/4) = +100.0;
	*mat_at(b, (3*nx)/4, (3*ny)/4) = -100.0;

	mat_print(p, "p", 0);

	// SOR
	double omega = 1.0;
	double rho = (dy2 * cos(M_PI / p.nrow) +
			dx2 * cos(M_PI / p.ncol)) / (dx2 + dy2);
	rho *= rho;

	int n;
	for (n = 1; n <= nt; n++) {
		double eps = poisson_iter(p, b, dx2, dy2, omega);
		mat_print(p, "p", n);
		if (eps < 1e-8) break;
		omega = n == 1 ? 2.0 / (2.0 - rho) : 4.0 / (4.0 - rho*omega);
	}
	WARNF("n=%d\n", n);

	mat_destroy(&p);
	mat_destroy(&b);

	return 0;
}

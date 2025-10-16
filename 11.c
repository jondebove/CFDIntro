/*
 * Copyright (c) 2025, Jonathan Debove
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "utils.h"

/*
 * Navier-Stikes equation: cavity flow
 * dU/dt + (U.grad)U = -grad(p)/rho + nu*lapla(U) with U=(u,v)
 *
 * Chorin projection method:
 * u(n+1/2,i,j) = u(n,i,j) - u(n,i,j)*dt/dx*(u(n,i,j)-u(n,i-1,j))
 *                       - v(n,i,j)*dt/dy*(u(n,i,j)-u(n,i,j-1))
 *                       + nu*dt/dx2*(u(n,i+1,j)-2*u(n,i,j)+u(n,i-1,j))
 *                       + nu*dt/dy2*(u(n,i,j+1)-2*u(n,i,j)+u(n,i,j-1))
 * v(n+1/2,i,j) = v(n,i,j) - u(n,i,j)*dt/dx*(v(n,i,j)-v(n,i-1,j))
 *                       - v(n,i,j)*dt/dy*(v(n,i,j)-v(n,i,j-1))
 *                       + nu*dt/dx2*(v(n,i+1,j)-2*v(n,i,j)+v(n,i-1,j))
 *                       + nu*dt/dy2*(v(n,i,j+1)-2*v(n,i,j)+v(n,i,j-1))
 *
 * p(i,j) = ((p(i+1,j)+p(i-1,j))* dy2 + (p(i,j+1)+p(i,j-1))* dx2 -
 *           dx2*dy2*b(i,j)) / (2*(dx2+dy2))
 * with b(i,j) = rho*((u(i+1,j)-u(i-1,j))/2dx+(v(i,j+1)-v(i,j-1))/2dy)
 *
 * u(n+1,i,j) = u(n+1/2,i,j) - dt/(2 rho dx)*(p(i+1,j)-p(i-1,j))
 * v(n+1,i,j) = v(n+1/2,i,j) - dt/(2 rho dy)*(p(i,j+1)-p(i,j-1))
 *
 * IC: u=v=p=0
 *
 * BC: u(n,i,2)=1, u(n,i,0)=u(n,0,j)=u(n,2,j)=0
 *     v(n,i,2)=v(n,i,0)=v(n,0,j)=v(n,2,j)=0
 *     p(i,2)=0, dpdy(i,0)=0, dpdx(0,j)=0, dpdx(2,j)=0
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
				epsmax = MAX(eps, epsmax);

				*mat_at(p, i, j) = p1;
			}
		}
	}
	return epsmax;
}

int main(void)
{
	int nx = 20;
	int ny = 20;
	int nt = 50;
	double dx = 2.0 / (nx - 1);
	double dy = 2.0 / (ny - 1);
	double dt = 0.01;
	double nu = 0.1;
	double rho = 1.0;

	double nudtdx2 = nu * dt / (dx * dx);
	double nudtdy2 = nu * dt / (dy * dy);
	WARNF("betax=%f betay=%f\n", nudtdx2, nudtdy2);

	int i, j, m, n;

	struct mat u0, u1, v0, v1, p, b;
	mat_create(&u0, nx, ny);
	mat_create(&u1, nx, ny);
	mat_create(&v0, nx, ny);
	mat_create(&v1, nx, ny);
	mat_create(&p, nx, ny);
	mat_create(&b, nx, ny);

	// BC
	for (i = 1; i < nx-1; i++) {
		*mat_at(u0, i, ny-1) = 1.0;
		*mat_at(u1, i, ny-1) = 1.0;
	}
	mat_print(u0, "u", 0);
	mat_print(v0, "v", 0);
	mat_print(p, "p", 0);

	for (n = 1; n <= nt; n++) {
		// Equation
		for (j = 1; j < ny - 1; j++) {
			for (i = 1; i < nx - 1; i++) {
				double cdtdx = *mat_at(u0, i, j) * dt / dx;
				double cdtdy = *mat_at(v0, i, j) * dt / dy;
				if (cdtdx > 1.0 || cdtdy > 1.0) {
					WARNF("CFLx=%f CFLy=%f\n", cdtdx, cdtdy);
				}

				*mat_at(u1, i, j) = *mat_at(u0, i, j)
					- cdtdx *
					(*mat_at(u0, i, j) -
					 *mat_at(u0, i-1, j))
					- cdtdy *
					(*mat_at(u0, i, j) -
					 *mat_at(u0, i, j-1))
					+ nudtdx2 *
					(*mat_at(u0, i+1, j) -
					 *mat_at(u0, i, j)*2 +
					 *mat_at(u0, i-1, j))
					+ nudtdy2 *
					(*mat_at(u0, i, j+1) -
					 *mat_at(u0, i, j)*2 +
					 *mat_at(u0, i, j-1));
				*mat_at(v1, i, j) = *mat_at(v0, i, j)
					- cdtdx *
					(*mat_at(v0, i, j) -
					 *mat_at(v0, i-1, j))
					- cdtdy *
					(*mat_at(v0, i, j) -
					 *mat_at(v0, i, j-1))
					+ nudtdx2 *
					(*mat_at(v0, i+1, j) -
					 *mat_at(v0, i, j)*2 +
					 *mat_at(v0, i-1, j))
					+ nudtdy2 *
					(*mat_at(v0, i, j+1) -
					 *mat_at(v0, i, j)*2 +
					 *mat_at(v0, i, j-1));
			}
		}
		// Pressure
		for (j = 1; j < ny - 1; j++) {
			for (i = 1; i < nx - 1; i++) {
				*mat_at(b, i, j) = rho / (dt * 2.0) * (
						(*mat_at(u1, i+1, j) -
						 *mat_at(u1, i-1, j)) / dx +
						(*mat_at(v1, i, j+1) -
						 *mat_at(v1, i, j-1)) / dy);
			}
		}
		double omega = 1.0;
		double r2 = 0.9995;
		double eps = 0.0;
		for (m = 0; m < 500; m++) {
			eps = poisson_iter(p, b, dx*dx, dy*dy, omega);
			for (i = 0; i < nx; i++) {
				*mat_at(p, i, 0) = *mat_at(p, i, 1);
			}
			for (j = 0; j < ny; j++) {
				*mat_at(p, 0, j) = *mat_at(p, 1, j);
				*mat_at(p, nx-1, j) = *mat_at(p, nx-2, j);
			}
			if (eps < 1e-4) {
				m++;
				break;
			}
			omega = m == 0 ? 2.0 / (2.0 - r2) : 4.0 / (4.0 - r2*omega);
		}
		WARNF("m=%d eps=%e\n", m, eps);

		// Projection
		for (j = 1; j < ny - 1; j++) {
			for (i = 1; i < nx - 1; i++) {
				*mat_at(u1, i, j) = *mat_at(u1, i, j) -
					dt/(dx*rho*2) *
					(*mat_at(p, i+1, j) -
					 *mat_at(p, i-1, j));
				*mat_at(v1, i, j) = *mat_at(v1, i, j) -
					dt/(dy*rho*2) *
					(*mat_at(p, i, j+1) -
					 *mat_at(p, i, j-1));
			}
		}
		mat_print(u1, "u", n);
		mat_print(v1, "v", n);
		mat_print(p, "p", n);

		SWAP(struct mat, u0, u1);
		SWAP(struct mat, v0, v1);
	}

	mat_destroy(&u0);
	mat_destroy(&u1);
	mat_destroy(&v0);
	mat_destroy(&v1);
	mat_destroy(&p);
	mat_destroy(&b);

	return 0;
}

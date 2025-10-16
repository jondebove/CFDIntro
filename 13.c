/*
 * Copyright (c) 2025, Jonathan Debove
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "utils.h"

/*
 * 1D linear convection
 * du/dt + c*du/dx = 0
 *
 * IC: u(0,0) = 1, u(0,elsewhere) = 0
 * BC: u(n,0) = 1
 */

enum scheme {
	UPWIND,
	LEAPFROG,
	LAX_FRIEDRICHS,
	LAX_WENDROFF,
	SCHEME_COUNT
};

static char const *scheme_names[SCHEME_COUNT] = {
	"UPWIND",
	"LEAPFROG",
	"LAX_FRIEDRICHS",
	"LAX_WENDROFF",
};

int main(void)
{
	int nx = 20;
	int nt = 50;
	double dx = 2.0 / (nx - 1);
	double dt = dx * 0.8;

	double c = 1.0;
	double cdtdx = c * dt / dx;
	WARNF("CFL=%f\n", cdtdx);

	enum scheme scheme = UPWIND;
	char const *s = scheme_names[scheme];

	int n = 0;
	int i = 0;

	struct mat u0, u1, u2;
	mat_create(&u0, nx, 1);
	mat_create(&u1, nx, 1);
	mat_create(&u2, nx, 1);

	// IC
	*mat_at(u0, 0, 0) = 1.0;
	mat_print(u0, s, 0);

	for (n = 1; n <= nt; n++) {
		SWAP(struct mat, u1, u2);
		SWAP(struct mat, u0, u1);
		for (i = 1; i < nx - 1; i++) {
			if (scheme == UPWIND || (scheme == LEAPFROG && n == 1)) {
				*mat_at(u0, i, 0) = *mat_at(u1, i, 0) -
					cdtdx *
					(*mat_at(u1, i, 0) -
					 *mat_at(u1, i-1, 0));
			} else if (scheme == LEAPFROG && n > 1) {
				*mat_at(u0, i, 0) = *mat_at(u2, i, 0) -
					cdtdx *
					(*mat_at(u1, i+1, 0) -
					 *mat_at(u1, i-1, 0));
			} else if (scheme == LAX_FRIEDRICHS) {
				*mat_at(u0, i, 0) =
					(*mat_at(u1, i+1, 0) +
					 *mat_at(u1, i-1, 0) -
					 cdtdx *
					 (*mat_at(u1, i+1, 0) -
					  *mat_at(u1, i-1, 0))) / 2.0;
			} else if (scheme == LAX_WENDROFF) {
				*mat_at(u0, i, 0) = *mat_at(u1, i, 0) -
					cdtdx / 2.0 *(
					 (*mat_at(u1, i+1, 0) -
					  *mat_at(u1, i-1, 0)) -
					 cdtdx *
					 (*mat_at(u1, i+1, 0) -
					  *mat_at(u1, i, 0)*2 +
					  *mat_at(u1, i-1, 0)));
			}
		}
		*mat_at(u0, 0, 0) = 1.0;
		*mat_at(u0, nx-1, 0) = *mat_at(u1, nx-1, 0) - cdtdx *
			(*mat_at(u1, nx-1, 0) - *mat_at(u1, nx-2, 0));
		mat_print(u0, s, n);
	}

	mat_destroy(&u0);
	mat_destroy(&u1);
	mat_destroy(&u2);

	return 0;
}

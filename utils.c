/*
 * Copyright (c) 2025, Jonathan Debove
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include <stdio.h>
#include <stdlib.h>

#include "utils.h"

struct mat *mat_create(struct mat *m, int nrow, int ncol)
{
	if (m && nrow > 0 && ncol > 0) {
		m->nrow = nrow;
		m->ncol = ncol;
		if ((m->data = calloc(nrow * ncol, sizeof(*m->data)))) {
			return m;
		}
	}
	return NULL;
}

void mat_destroy(struct mat *m)
{
	if (m) {
		free(m->data);
		m->nrow = 0;
		m->ncol = 0;
		m->data = NULL;
	}
}

void mat_print(struct mat m, char *s, int n)
{
	int i, j;
	static int header = 1;
	if (header) {
		puts("s n i j v");
		header = 0;
	}
	if (!s) s = "none";
	for (i = 0; i < m.nrow; i++) {
		for (j = 0; j < m.ncol; j++) {
			printf("%s %d %d %d %e\n", s, n, i, j,
					*mat_at(m, i, j));
		}
	}
}

double mat_dist(struct mat a, struct mat b)
{
	int n = a.nrow * a.ncol;
	assert(b.nrow * b.ncol == n);
	if (n <= 0) return 0.0;

	double ans = 0.0;
	while (n--) {
		double tmp = ABS(a.data[n] - b.data[n]);
		if (tmp != 0.0) {
			tmp /= MAX(ABS(a.data[n]), ABS(b.data[n]));
			ans = MAX(ans, tmp);
		}
	}
	return ans;
}

/*
 * Copyright (c) 2025, Jonathan Debove
 * SPDX-License-Identifier: BSD-2-Clause
 */

#ifndef UTILS_H
#define UTILS_H

#include <assert.h>
#include <stdio.h>

#ifndef M_PI
#define M_PI 3.14159265358979323844
#endif

#define WARNF(...) fprintf(stderr, __VA_ARGS__)

#define ABS(x) ((x) >= 0 ? (x) : -(x))
#define MIN(x, y) ((x) < (y) ? (x) : (y))
#define MAX(x, y) ((x) > (y) ? (x) : (y))

#define SWAP(type, x, y) do { \
	type swap__tmp = (x); \
	(x) = (y); \
	(y) = swap__tmp; \
} while (0)

struct mat {
	int nrow;
	int ncol;
	double *data;
};

struct mat *mat_create(struct mat *m, int nrow, int ncol);

void mat_destroy(struct mat *m);

static inline
double *mat_at(struct mat const m, int i, int j)
{
	assert(i >= 0 && i < m.nrow);
	assert(j >= 0 && j < m.ncol);

	return &m.data[i + m.nrow * j];
}

void mat_print(struct mat m, char *s, int n);

double mat_dist(struct mat a, struct mat b);

#endif /* UTILS_H */

/*
 * Copyright (c) 2025, Jonathan Debove
 * SPDX-License-Identifier: BSD-2-Clause
 */

#ifndef UTILS_H
#define UTILS_H

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#define WARN(...) fprintf(stderr, __VA_ARGS__)

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

static inline
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

static inline
void mat_destroy(struct mat *m)
{
	if (m) {
		free(m->data);
		m->nrow = 0;
		m->ncol = 0;
		m->data = NULL;
	}
}

static inline
double *mat_at(struct mat const m, int i, int j)
{
	assert(i >= 0 && i < m.nrow);
	assert(j >= 0 && j < m.ncol);

	return &m.data[i + m.nrow * j];
}

static inline
void mat_print(struct mat m, int n)
{
	int i, j;
	if (n == 0) puts("n i j u");
	for (i = 0; i < m.nrow; i++) {
		for (j = 0; j < m.ncol; j++) {
			printf("%d %d %d %e\n", n, i, j, *mat_at(m, i, j));
		}
	}
}

#endif /* UTILS_H */

#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>

#define COUNT_OF(x) (sizeof(x) / sizeof(x[0]))
#define C(m, x, y) (m->data[(x) * (m->size) + (y)])

typedef struct {
    float* data;
    uint32_t size;
} matrix_t;

typedef struct {
    matrix_t** data;
    uint32_t size;
} mblock_t;

matrix_t* matrix_alloc(uint32_t size) {
    // calloc will memset 0
    float* data = (float*)calloc(size * size, sizeof(float));
    matrix_t* m = malloc(sizeof(matrix_t));

    m->data = data;
    m->size = size;

    return m;
}

void matrix_free(matrix_t* m) {
    free(m->data);
    free(m);
}

void matrix_cpy(matrix_t* dst, const matrix_t* src) {
    assert(dst->size == src->size);
    memcpy(dst->data, src->data, src->size * src->size * sizeof(float));
}

void matrix_zero(matrix_t* dst) {
    memset(dst->data, 0, dst->size * dst->size * sizeof(float));
}

void matrix_print(const matrix_t* m) {
    for (uint32_t i = 0; i < m->size; i++) {
        for (uint32_t j = 0; j < m->size; j++) {
            printf("%10.2f", C(m, i, j));
        }
        printf("\n");
    }
}

mblock_t* mblock_alloc(uint32_t size, uint32_t msize) {
    matrix_t** data = (matrix_t**)calloc(size * size, sizeof(matrix_t*));
    mblock_t* b = malloc(sizeof(mblock_t));

    for (uint32_t i = 0; i < size * size; i++) {
        data[i] = matrix_alloc(msize);
    }

    b->data = data;
    b->size = size;

    return b;
}

void mblock_free(mblock_t* b) {
    for (uint32_t i = 0; i < b->size * b->size; i++) {
        matrix_free(b->data[i]);
    }

    free(b->data);
    free(b);
}

void mblock_cpy(mblock_t* dst, const mblock_t* src) {
    for (uint32_t i = 0; i < src->size; i++) {
        for (uint32_t j = 0; j < src->size; j++) {
            matrix_cpy(C(dst, i, j), C(src, i, j));
        }
    }
}

void mblock_zero(mblock_t* dst) {
    for (uint32_t i = 0; i < dst->size; i++) {
        for (uint32_t j = 0; j < dst->size; j++) {
            matrix_zero(C(dst, i, j));
        }
    }
}

void mblock_print(const mblock_t* m) {
    for (uint32_t i = 0; i < m->size; i++) {
        for (uint32_t j = 0; j < m->size; j++) {
            printf("== x: %d, y: %d ==\n", i, j);
            matrix_print(C(m, i, j));
            printf("==================\n");
        }
        printf("\n");
    }
}

void mblock_fill(mblock_t* dst, const matrix_t* src) {
    assert(src->size % dst->size == 0);
    for (uint32_t i = 0; i < src->size; i++) {
        for (uint32_t j = 0; j < src->size; j++) {
            matrix_t* m = C(dst, i / (dst->data[0]->size), j / (dst->data[0]->size));
            C(m, i % (dst->data[0]->size), j % (dst->data[0]->size)) = C(src, i, j);
        }
    }
}

void matrix_reconstruct(matrix_t* dst, const mblock_t* src) {
    assert(dst->size % src->size == 0);

    for (uint32_t i = 0; i < dst->size; i++) {
        for (uint32_t j = 0; j < dst->size; j++) {
            matrix_t* submatrix = C(src, i / (src->data[0]->size), j / (src->data[0]->size));
            C(dst, i, j) = C(submatrix, i % (src->data[0]->size), j % (src->data[0]->size));
        }
    }
}

float det(const matrix_t* m) {
    if (m->size == 1)
        return m->data[0];
    if (m->size == 2)
        return m->data[0] * m->data[3] - m->data[1] * m->data[2];

    float d = 0;
    for (uint32_t i = 0; i < m->size; i++) {
        matrix_t* minor = matrix_alloc(m->size - 1);
        for (uint32_t j = 1; j < m->size; j++) {
            for (uint32_t k = 0; k < m->size; k++) {
                if (k < i)
                    C(minor, j - 1, k) = C(m, j, k);
                else if (k > i)
                    C(minor, j - 1, k - 1) = C(m, j, k);
            }
        }
        d += ((i % 2 == 0) ? 1 : -1) * m->data[i] * det(minor);
        matrix_free(minor);
    }

    return d;
}

float cof(const matrix_t* m, uint32_t row, uint32_t col) {
    matrix_t* minor = matrix_alloc(m->size - 1);
    uint32_t minorRow = 0, minorCol = 0;

    for (uint32_t i = 0; i < m->size; i++) {
        if (i == row)
            continue;
        minorCol = 0;
        for (uint32_t j = 0; j < m->size; j++) {
            if (j == col)
                continue;
            C(minor, minorRow, minorCol) = C(m, i, j);
            minorCol++;
        }
        minorRow++;
    }

    float c = ((row + col) % 2 == 0 ? 1 : -1) * det(minor);
    matrix_free(minor);

    return c;
}

void matrix_inv(const matrix_t* m, matrix_t* out) {
    assert(m->size == out->size);

    float d = det(m);
    if (m->size == 1) {
        out->data[0] = 1 / d;
        return ;
    }
    if (d == 0) {
        for (uint32_t i = 0; i < m->size * m->size; i++) {
            out->data[i] = 0;
        }
        return;
    }

    for (uint32_t i = 0; i < m->size; i++) {
        for (uint32_t j = 0; j < m->size; j++) {
            C(out, j, i) = cof(m, i, j) / d;
        }
    }
}

void matrix_add(const matrix_t* a, const matrix_t* b, matrix_t* out) {
    assert(a->size == b->size);
    assert(a->size == out->size);

    for (uint32_t i = 0; i < a->size * a->size; i++) {
        out->data[i] = a->data[i] + b->data[i];
    }
}

void matrix_sub(const matrix_t* a, const matrix_t* b, matrix_t* out) {
    assert(a->size == b->size);
    assert(a->size == out->size);

    for (uint32_t i = 0; i < a->size * a->size; i++) {
        out->data[i] = a->data[i] - b->data[i];
    }
}

void matrix_mul(const matrix_t* a, const matrix_t* b, matrix_t* out) {
    assert(a->size == b->size);
    assert(a->size == out->size);

    matrix_t* tmp = matrix_alloc(a->size);

    for (uint32_t i = 0; i < a->size; i++) {
        for (uint32_t j = 0; j < a->size; j++) {
            for (uint32_t k = 0; k < a->size; k++) {
                C(tmp, i, j) += C(a, i, k) * C(b, k, j);
            }
        }
    }

    matrix_cpy(out, tmp);
    matrix_free(tmp);
}

void matrix_div(const matrix_t* a, const matrix_t* b, matrix_t* out) {
    assert(a->size == b->size);
    assert(a->size == out->size);

    matrix_t* inverseB = matrix_alloc(a->size);
    matrix_inv(b, inverseB);
    matrix_mul(a, inverseB, out);
    matrix_free(inverseB);
}

void matrix_mul_imm(const matrix_t* a, float imm, matrix_t* out) {
    assert(a->size == out->size);

    for (uint32_t i = 0; i < a->size * a->size; i++) {
        out->data[i] = a->data[i] * imm;
    }
}

float matrix_err(const matrix_t* a, const matrix_t* b) {
    assert(a->size == b->size);

    float norm = 0.0f;
    for (uint32_t i = 0; i < a->size * a->size; i++) {
        float diff = a->data[i] - b->data[i];
        norm += diff * diff;
    }

    return sqrt(norm);
}

void matrix_i(matrix_t* output) {
    for (uint32_t i = 0; i < output->size; ++i) {
        for (uint32_t j = 0; j < output->size; ++j) {
            C(output, i, j) = (i == j) ? 1.0f : 0.0f;
        }
    }
}

uint32_t matrix_sqrt(const matrix_t* a, matrix_t* output, float err, uint64_t maxiter) {
    matrix_t* Y = matrix_alloc(a->size);
    matrix_t* Z = matrix_alloc(a->size);
    matrix_t* nextY = matrix_alloc(a->size);
    matrix_t* nextZ = matrix_alloc(a->size);
    matrix_t* tmp = matrix_alloc(a->size);

    matrix_cpy(Y, a);
    matrix_i(Z);

    uint32_t i = 0;
    uint8_t done = 0;

    if (a->size == 1) {
        C(output, 0, 0) = sqrtf(C(a, 0, 0));
        return 0;
    }
    while (!done) {
        matrix_inv(Z, tmp);
        matrix_add(Y, tmp, tmp);
        matrix_mul_imm(tmp, 0.5, nextY);

        matrix_inv(Y, tmp);
        matrix_add(Z, tmp, tmp);
        matrix_mul_imm(tmp, 0.5, nextZ);

        done = (matrix_err(Y, nextY) < err && matrix_err(Z, nextZ) < err);

        matrix_cpy(Y, nextY);
        matrix_cpy(Z, nextZ);
        i++;
        if(i > maxiter) {
            done = true;
        }
    }

    matrix_cpy(output, Y);

    matrix_free(Y);
    matrix_free(Z);
    matrix_free(nextY);
    matrix_free(nextZ);

    return i;
}

void matrix_transpose(const matrix_t* a, matrix_t* out) {
    for (uint32_t i = 0; i < a->size; i++) {
        for (uint32_t j = 0; j < a->size; j++) {
            C(out, i, j) = C(a, j, i);
        }
    }
}

void matrix_random_l(matrix_t* m, int64_t min_val, int64_t max_val) {
    for (uint32_t i = 0; i < m->size; ++i) {
        for (uint32_t j = 0; j <= i; ++j) {
            if (min_val == max_val) {
                C(m, i, j) = min_val;
            } else {
                C(m, i, j) = rand() % (max_val - min_val + 1) + min_val;
            }
        }
    }
}

void matrix_random_pds(matrix_t* m, int64_t min_val, int64_t max_val) {
    matrix_t* tmp = matrix_alloc(m->size);

    matrix_random_l(m, min_val, max_val);
    matrix_transpose(m, tmp);
    matrix_mul(m, tmp, m);

    matrix_free(tmp);
}


uint32_t matrix_cmp(const matrix_t* m1, const matrix_t* m2, float perr) {
    uint32_t c = 0;
    if (m1->size != m2->size) {
        return UINT32_MAX;
    }

    for (uint32_t i = 0; i < m1->size; ++i) {
        for (uint32_t j = 0; j < m1->size; ++j) {
            float diff = fabsf(C(m1, i, j) - C(m2, i, j));
            float allowed_err = perr * fabsf(C(m1, i, j));
            if (diff > allowed_err) {
                c++;
            }
        }
    }

    return c;
}

void decomp2(const matrix_t* matrix, matrix_t* l) {
    matrix_t* a = matrix_alloc(matrix->size);

    matrix_cpy(a, matrix);
    matrix_zero(l);

    for (uint32_t k = 0; k < matrix->size - 1; k++) {
        C(l, k, k) = sqrt(C(a, k, k));
        for (uint32_t i = k + 1; i < matrix->size; i++) {
            C(l, i, k) = C(a, i, k) / C(l, k, k);
        }
        for (uint32_t j = k + 1; j < matrix->size; j++) {
            for (uint32_t i = j; i < matrix->size; i++) {
                C(a, i, j) -= C(l, i, k) * C(l, j, k);
            }
        }
    }
    // last element
    l->data[l->size * l->size - 1] = sqrt(a->data[l->size * l->size - 1]);

    matrix_free(a);
}

void decomp2_block(const mblock_t* mblock, mblock_t* l) {
    mblock_t* a = mblock_alloc(mblock->size, mblock->data[0]->size);
    matrix_t* tmp1 = matrix_alloc(mblock->data[0]->size);
    matrix_t* tmp2 = matrix_alloc(mblock->data[0]->size);
    uint32_t n = mblock->size;

    mblock_cpy(a, mblock);
    mblock_zero(l);

    for (uint32_t k = 0; k < n - 1; k++) {
        matrix_sqrt(C(a, k, k), C(l, k, k), 0.0001, 20);
        for (uint32_t i = k + 1; i < n; i++) {
            matrix_div(C(a, i, k), C(l, k, k), C(l, i, k));
        }
        for (uint32_t j = k + 1; j < n; j++) {
            for (uint32_t i = j; i < n; i++) {
                matrix_transpose(C(l, j, k), tmp1);
                matrix_mul(C(l, i, k), tmp1, tmp2);
                matrix_sub(C(a, i, j), tmp2, C(a, i, j));
            }
        }
    }
    matrix_sqrt(C(a, n - 1, n - 1), C(l, n - 1, n - 1), 0.0001, 20);

    mblock_free(a);
}

void decomp1(const float* matrix, float* lower, uint32_t n) {
    memset(lower, 0, n * n * sizeof(float));

    for (uint32_t i = 0; i < n; i++) {
        for (uint32_t j = 0; j <= i; j++) {
            float sum = 0;

            if (j == i) {
                for (uint32_t k = 0; k < j; k++) {
                    sum += pow(lower[j * n + k], 2);
                }
                lower[j * n + j] = sqrt(matrix[j * n + j] - sum);
            } else {
                for (uint32_t k = 0; k < j; k++) {
                    sum += (lower[i * n + k] * lower[j * n + k]);
                }
                lower[i * n + j] = (matrix[i * n + j] - sum) / lower[j * n + j];
            }
        }
    }
}

int32_t main(int argc, char**argv) {
    srand(time(NULL));

    const uint32_t n = 64;
    const uint32_t block_size = 8;
    assert(128 % 8 == 0);
    const uint32_t blocks_count = n / block_size;

    matrix_t* matrix = matrix_alloc(n);
    matrix_t* m_out = matrix_alloc(n);

    mblock_t* mblock = mblock_alloc(block_size, blocks_count);
    mblock_t* mb_output = mblock_alloc(block_size, blocks_count);

    matrix_t* decomposed = matrix_alloc(n);
    matrix_t* transposed = matrix_alloc(n);
    matrix_t* reconstructed = matrix_alloc(n);


    matrix_random_pds(matrix, -50, 50);
    mblock_fill(mblock, matrix);

    puts("=======================");
    puts("Matrix to decompose");
    matrix_print(matrix);
    puts("=======================");

    puts("Normal decomp");
    decomp2(matrix, m_out);
    matrix_print(m_out);
    puts("=======================");

    puts("Block decomp");
    decomp2_block(mblock, mb_output);
    matrix_reconstruct(decomposed, mb_output);
    matrix_print(decomposed);
    puts("=======================");

    puts("L*L.T");
    matrix_reconstruct(decomposed, mb_output);
    matrix_transpose(decomposed, transposed);
    matrix_mul(decomposed, transposed, reconstructed);
    matrix_print(reconstructed);
    puts("=======================");

    matrix_print(matrix);
    uint32_t errors = matrix_cmp(matrix, reconstructed, 0.01);
    if(errors == 0) {
        printf("WORKED!");
    } else {
        printf("UNWORKED! Err: %u", errors);
    }

    matrix_free(matrix);
    mblock_free(mblock);
    matrix_free(decomposed);
    matrix_free(transposed);
    matrix_free(reconstructed);


    return 0;
}

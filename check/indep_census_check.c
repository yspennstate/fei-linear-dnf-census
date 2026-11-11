/* Independent check of a census of linear monotone DNFs with essential dimension n.
   Input: one family per line, terms as decimal bitmasks over n variables.
   For every line it checks that the family is a linear antichain using all n variables,
   computes a canonical form (least sorted relabelled term list over relabellings that respect a
   vertex invariant) and the automorphism count, the integer Walsh transform of f = 1 - 2F,
   H/I and min_i E_i/Inf_i. Reports duplicates, sum n!/|Aut|, and the two largest values of each ratio.
   Written separately from the census code; shares no source with it. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

#define MAXT 64
static int n, N;
static int fam[MAXT], nt;
static int used[8], mp[8];
static uint64_t vinv[8];
static int best[MAXT], cur[MAXT], have_best;
static long long aut_count;

static int cmpint(const void *a, const void *b) { int x = *(const int *)a, y = *(const int *)b; return (x > y) - (x < y); }
static int cmpu64(const void *a, const void *b) { uint64_t x = *(const uint64_t *)a, y = *(const uint64_t *)b; return (x > y) - (x < y); }

static void try_relabel(void) {
    for (int k = 0; k < nt; k++) { int y = 0; for (int j = 0; j < n; j++) if (fam[k] >> j & 1) y |= 1 << mp[j]; cur[k] = y; }
    qsort(cur, nt, sizeof(int), cmpint);
    if (!have_best) { memcpy(best, cur, nt * sizeof(int)); have_best = 1; aut_count = 1; return; }
    int c = 0;
    for (int k = 0; k < nt; k++) if (cur[k] != best[k]) { c = cur[k] < best[k] ? -1 : 1; break; }
    if (c < 0) { memcpy(best, cur, nt * sizeof(int)); aut_count = 1; }
    else if (c == 0) aut_count++;
}

static uint64_t pos_inv[8];
static void rec2(int pos) {
    if (pos == n) { try_relabel(); return; }
    for (int v = 0; v < n; v++) if (!used[v] && vinv[v] == pos_inv[pos]) {
        used[v] = 1; mp[v] = pos; rec2(pos + 1); used[v] = 0;
    }
}

static double phi(double p, double q) {
    double s = p + q, r = 0;
    if (p > 0) r += p * log2(s / p);
    if (q > 0) r += q * log2(s / q);
    return r;
}

typedef struct { double v; char line[400]; } Top;
static void push_top(Top t[2], double v, const char *line) {
    if (v > t[0].v) { t[1] = t[0]; t[0].v = v; strncpy(t[0].line, line, 399); }
    else if (v > t[1].v) { t[1].v = v; strncpy(t[1].line, line, 399); }
}

typedef struct { unsigned char len; unsigned char b[MAXT]; } Canon;
static int cmpcanon(const void *a, const void *b) {
    const Canon *x = a, *y = b;
    if (x->len != y->len) return x->len - y->len;
    return memcmp(x->b, y->b, x->len);
}

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: %s n file\n", argv[0]); return 2; }
    n = atoi(argv[1]); N = 1 << n;
    FILE *fp = fopen(argv[2], "r"); if (!fp) { perror("open"); return 2; }
    size_t cap = 1 << 20, cnt = 0; Canon *canon = malloc(cap * sizeof(Canon));
    double fact = 1; for (int k = 2; k <= n; k++) fact *= k;
    long double labeled = 0;
    long bad = 0;
    Top topH[2] = {{-1, ""}, {-1, ""}}, topE[2] = {{-1, ""}, {-1, ""}};
    char line[4096];
    int F[256]; long W[256];
    while (fgets(line, sizeof line, fp)) {
        char *s = line; nt = 0; int x;
        int off;
        while (sscanf(s, "%d%n", &x, &off) == 1) { fam[nt++] = x; s += off; }
        if (!nt) continue;
        line[strcspn(line, "\n")] = 0;
        /* validity: masks in range, union is everything, pairwise intersections <= 1, antichain, distinct */
        int uni = 0, ok = 1;
        for (int a = 0; a < nt; a++) { if (fam[a] <= 0 || fam[a] >= N) ok = 0; uni |= fam[a]; }
        for (int a = 0; a < nt && ok; a++) for (int b = a + 1; b < nt; b++) {
            int c = __builtin_popcount(fam[a] & fam[b]);
            if (fam[a] == fam[b] || c > 1) ok = 0;
            if ((fam[a] & fam[b]) == fam[a] || (fam[a] & fam[b]) == fam[b]) ok = 0;
        }
        if (uni != N - 1) ok = 0;
        if (!ok) { bad++; if (bad < 10) fprintf(stderr, "invalid: %s\n", line); continue; }
        /* vertex invariant: multiset of incident term sizes */
        for (int v = 0; v < n; v++) { uint64_t t = 0; for (int a = 0; a < nt; a++) if (fam[a] >> v & 1) t += (uint64_t)1 << (4 * (__builtin_popcount(fam[a]) - 1)); vinv[v] = t; }
        uint64_t sorted[8]; memcpy(sorted, vinv, sizeof sorted); qsort(sorted, n, sizeof(uint64_t), cmpu64);
        memcpy(pos_inv, sorted, sizeof sorted);
        memset(used, 0, sizeof used); have_best = 0; aut_count = 0;
        rec2(0);
        if (cnt == cap) { cap *= 2; canon = realloc(canon, cap * sizeof(Canon)); }
        canon[cnt].len = (unsigned char)nt; for (int k = 0; k < nt; k++) canon[cnt].b[k] = (unsigned char)best[k]; cnt++;
        labeled += (long double)fact / aut_count;
        /* truth table and integer Walsh transform of f = 1 - 2F */
        for (int y = 0; y < N; y++) { int v = 0; for (int a = 0; a < nt; a++) if ((y & fam[a]) == fam[a]) { v = 1; break; } F[y] = v; W[y] = 1 - 2 * v; }
        for (int h = 1; h < N; h <<= 1) for (int i = 0; i < N; i += h << 1) for (int j = i; j < i + h; j++) { long u = W[j], w = W[j + h]; W[j] = u + w; W[j + h] = u - w; }
        double den = (double)N * N, H = 0, I = 0, minr = 1e300;
        for (int S = 0; S < N; S++) { double p = W[S] * (double)W[S] / den; if (p > 0) H -= p * log2(p); }
        for (int i = 0; i < n; i++) {
            double inf = 0, e = 0;
            for (int S = 0; S < N; S++) if (!(S >> i & 1)) {
                double p = W[S] * (double)W[S] / den, q = W[S | 1 << i] * (double)W[S | 1 << i] / den;
                inf += q; e += phi(p, q);
            }
            I += inf; if (e / inf < minr) minr = e / inf;
        }
        push_top(topH, H / I, line); push_top(topE, minr, line);
    }
    qsort(canon, cnt, sizeof(Canon), cmpcanon);
    long dup = 0; for (size_t k = 1; k < cnt; k++) if (!cmpcanon(&canon[k], &canon[k - 1])) dup++;
    printf("{\"n\":%d,\"lines_valid\":%zu,\"invalid\":%ld,\"duplicate_classes\":%ld,\"sum_n_fact_over_aut\":%.0Lf,\n", n, cnt, bad, dup, labeled);
    printf(" \"max_H_over_I\":%.15f,\"attainer_H\":\"%s\",\"second_H_over_I\":%.15f,\n", topH[0].v, topH[0].line, topH[1].v);
    printf(" \"max_min_E_over_Inf\":%.15f,\"attainer_E\":\"%s\",\"second_min_E_over_Inf\":%.15f}\n", topE[0].v, topE[0].line, topE[1].v);
    return 0;
}

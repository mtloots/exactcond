/* exactcond.c -- exact conditional reference distributions, pure C back-end.
 *
 * THE ONE IDEA. A nuisance parameter is removed by conditioning on a statistic that is sufficient
 * for it; what survives is a distribution free of that nuisance, and a reference drawn from it is
 * exact at every sample size. The package supplies the two conditionings that recur in this
 * programme, and nothing else:
 *
 *   DISCRETE, both margins fixed. In an r x c table the row and column totals are sufficient for
 *   the two sets of marginal parameters, so conditioning on them leaves the Fisher-Yates law, which
 *   carries no nuisance at all. It is sampled here EXACTLY, and by an argument that needs no
 *   appeal to Patefield's algorithm: the conditional law given both margins is precisely the law of
 *   the table formed by pairing a fixed vector of row labels with a UNIFORMLY RANDOM PERMUTATION of
 *   the column labels. A Fisher-Yates shuffle therefore samples it directly, in O(N) per draw.
 *
 *   CONTINUOUS, direction only. In a linear model under a null in which the tested structure does
 *   not enter the data-generating process, a statistic invariant to the regression coefficients and
 *   to the error scale depends on the data only through the residual DIRECTION, which is uniform on
 *   the unit sphere of the residual space whatever those nuisance parameters are. Drawing that
 *   direction gives an exact reference, and it stays exact when the statistic is a maximum over
 *   candidate structures -- which is what makes it valid AFTER a graph has been selected by fit.
 *
 * REPRODUCIBILITY. Every replicate draws from its OWN splitmix64 stream, seeded from the base seed
 * and the replicate index. Results are therefore identical whatever the thread count, and identical
 * across the R and Python fronts, which is what the parity harness checks. The stream is the same
 * generator used by the other C back-ends in this programme.
 *
 * WHAT IS NOT HERE. The regime change-in-dependence test is exact by characteristic-function
 * inversion of a Jacobi ensemble, not by conditioning and simulation, so it has no business in this
 * file. Nothing here is asymptotic and nothing here is a bootstrap.
 */
#include <math.h>
#include <stdlib.h>
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ---------------------------------------------------------------- streams */
static unsigned long long sm64(unsigned long long *s){
  unsigned long long z = (*s += 0x9E3779B97F4A7C15ULL);
  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
  z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
  return z ^ (z >> 31);
}
static double ec_unif(unsigned long long *s){ return (sm64(s) >> 11) * 0x1.0p-53; }
static unsigned long long ec_seed(unsigned long long base, long rep){
  unsigned long long s = base ^ (0x9E3779B97F4A7C15ULL * (unsigned long long)(rep + 1));
  (void)sm64(&s); (void)sm64(&s);
  return s;
}
/* EXACTLY uniform on {0,...,n-1}. Taking floor(n * unif) would be uniform only to within 2^-53 n,
   and a package whose entire claim is exactness should not open with an approximation. */
static unsigned long long ec_below(unsigned long long *s, unsigned long long n){
  unsigned long long t = (0ULL - n) % n, x;      /* 2^64 mod n */
  do { x = sm64(s); } while (x < t);
  return x % n;
}
static double ec_norm(unsigned long long *s){
  double u1, u2;
  do { u1 = ec_unif(s); } while (u1 <= 0.0);
  u2 = ec_unif(s);
  return sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
}
static void ec_shuffle(int *a, int N, unsigned long long *s){
  int i, j, t;
  for (i = N - 1; i > 0; i--){
    j = (int) ec_below(s, (unsigned long long)(i + 1));
    t = a[i]; a[i] = a[j]; a[j] = t;
  }
}

/* ------------------------------------------------------------- statistics */
/* which = 0  linear-by-linear, T = sum_ij u_i v_j x_ij. Goodman's association parameter is the
 *            natural ordinal edge statistic and is sufficient for it given both margins.
 * which = 1  score dispersion, T = sum_i n_i (vbar_i - vbar)^2, for NOMINAL rows against an
 *            ordinal column. Its components are the per-row mean scores themselves, so the test
 *            is on the estimand rather than on a rank proxy for it.
 * which = 2  Pearson X^2.      which = 3  likelihood ratio G^2.
 */
static double ec_stat(const int *x, int r, int c, const double *u, const double *v, int which){
  int i, j; double N = 0.0, T = 0.0;
  for (i = 0; i < r * c; i++) N += x[i];
  if (N <= 0.0) return 0.0;
  if (which == 0){
    for (i = 0; i < r; i++) for (j = 0; j < c; j++) T += u[i] * v[j] * x[i * c + j];
    return T;
  }
  if (which == 1){
    double gm = 0.0;
    for (i = 0; i < r; i++) for (j = 0; j < c; j++) gm += v[j] * x[i * c + j];
    gm /= N;
    for (i = 0; i < r; i++){
      double n = 0.0, s = 0.0;
      for (j = 0; j < c; j++){ n += x[i * c + j]; s += v[j] * x[i * c + j]; }
      if (n > 0.0){ double d = s / n - gm; T += n * d * d; }
    }
    return T;
  }
  { /* which = 2 or 3 */
    double *rs = (double*) calloc((size_t) r, sizeof(double));
    double *cs = (double*) calloc((size_t) c, sizeof(double));
    if (!rs || !cs){ free(rs); free(cs); return 0.0; }
    for (i = 0; i < r; i++) for (j = 0; j < c; j++){ rs[i] += x[i*c+j]; cs[j] += x[i*c+j]; }
    for (i = 0; i < r; i++) for (j = 0; j < c; j++){
      double e = rs[i] * cs[j] / N, o = x[i * c + j];
      if (e <= 0.0) continue;
      if (which == 2){ double d = o - e; T += d * d / e; }
      else if (o > 0.0) T += 2.0 * o * log(o / e);
    }
    free(rs); free(cs);
    return T;
  }
}

/* Build the fixed label templates from the two margins. */
static int ec_labels(const int *m, int k, int *lab){
  int i, j, p = 0;
  for (i = 0; i < k; i++) for (j = 0; j < m[i]; j++) lab[p++] = i;
  return p;
}

/* ------------------------------------------- one table, both margins fixed */
void C_ec_rc_null(const int *rows, const int *cols, const int *r, const int *c,
                  const int *B, const int *which, const double *u, const double *v,
                  const int *seed, double *out){
  int R = *r, C = *c, nb = *B, w = *which, i, N = 0;
  int *rl, *ct;
  for (i = 0; i < R; i++) N += rows[i];
  rl = (int*) malloc((size_t) N * sizeof(int));
  ct = (int*) malloc((size_t) N * sizeof(int));
  if (!rl || !ct){ free(rl); free(ct); return; }
  ec_labels(rows, R, rl);
  ec_labels(cols, C, ct);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (i = 0; i < nb; i++){
    unsigned long long s = ec_seed((unsigned long long) *seed, i);
    int *cl = (int*) malloc((size_t) N * sizeof(int));
    int *tb = (int*) calloc((size_t)(R * C), sizeof(int));
    int k;
    if (!cl || !tb){ free(cl); free(tb); out[i] = 0.0; continue; }
    memcpy(cl, ct, (size_t) N * sizeof(int));
    ec_shuffle(cl, N, &s);
    for (k = 0; k < N; k++) tb[rl[k] * C + cl[k]]++;
    out[i] = ec_stat(tb, R, C, u, v, w);
    free(cl); free(tb);
  }
  free(rl); free(ct);
}

/* Tail probabilities from an observed value and a null sample. out: p, midp, lo, hi.
 * The Monte Carlo interval is the Wilson interval for the underlying tail probability, and it is
 * reported because a p-value estimated from B draws is an ESTIMATE; quoting it bare invites the
 * reader to believe a precision the number does not have. */
static void ec_tail(double obs, const double *nul, int B, double *out){
  int i, ge = 0, eq = 0; double ph, z = 1.959963984540054, den, cen, hw;
  for (i = 0; i < B; i++){
    if (nul[i] > obs + 1e-12) ge++;
    else if (fabs(nul[i] - obs) <= 1e-12) eq++;
  }
  out[0] = (1.0 + ge + eq) / (B + 1.0);
  out[1] = (1.0 + ge + 0.5 * eq) / (B + 1.0);
  ph = (double)(ge + eq) / B;
  den = 1.0 + z * z / B;
  cen = (ph + z * z / (2.0 * B)) / den;
  hw  = z * sqrt(ph * (1.0 - ph) / B + z * z / (4.0 * (double) B * B)) / den;
  out[2] = cen - hw < 0.0 ? 0.0 : cen - hw;
  out[3] = cen + hw > 1.0 ? 1.0 : cen + hw;
}

void C_ec_pvalue(const double *obs, const double *nul, const int *B, double *out){
  ec_tail(*obs, nul, *B, out);
}

void C_ec_rc_test(const int *x, const int *r, const int *c, const int *B, const int *which,
                  const double *u, const double *v, const int *seed, double *out){
  int R = *r, C = *c, nb = *B, i, j;
  int *rows = (int*) calloc((size_t) R, sizeof(int));
  int *cols = (int*) calloc((size_t) C, sizeof(int));
  double *nul = (double*) malloc((size_t) nb * sizeof(double));
  double t4[4];
  if (!rows || !cols || !nul){ free(rows); free(cols); free(nul); return; }
  for (i = 0; i < R; i++) for (j = 0; j < C; j++){ rows[i] += x[i*C+j]; cols[j] += x[i*C+j]; }
  out[0] = ec_stat(x, R, C, u, v, *which);
  C_ec_rc_null(rows, cols, r, c, B, which, u, v, seed, nul);
  ec_tail(out[0], nul, nb, t4);
  out[1] = t4[0]; out[2] = t4[1]; out[3] = t4[2]; out[4] = t4[3];
  out[5] = nul[0];
  for (i = 1; i < nb; i++) if (nul[i] > out[5]) out[5] = nul[i];   /* largest null draw */
  free(rows); free(cols); free(nul);
}

/* ------------------------------------------- many items, one group label, ONE test
 * THE MULTI-ITEM OMNIBUS. m categorical items are recorded on n respondents who fall into g
 * groups. The question "does ANY item differ across the groups" is answered by one statistic,
 * the SUM or the MAXIMUM over items of a per-item test statistic, referred to the law obtained
 * by permuting the group labels. That law conditions on every item's margin and on the group
 * sizes, so it is exact conditional in the same sense as the r x c tests above, and it yields
 * ONE p-value for all m items: there is no multiplicity to correct because there is no
 * multiplicity. The per-item statistics are referred to the SAME permutations, which gives the
 * unadjusted per-item p-values a closed procedure uses after the omnibus has rejected.
 *
 * Two exchangeabilities are offered. Respondents exchangeable (cl absent): the labels are
 * permuted across respondents, which is what Fisher's test assumes. Clusters exchangeable (cl
 * given, every cluster inside one group): the labels are permuted across CLUSTERS, each cluster
 * keeping its members, which is the exact answer to intra-cluster correlation -- an interviewer
 * or a fieldworker team -- and costs nothing but the number of clusters as the effective sample.
 *
 * Statistics per item: which = 0 Pearson X^2 on the g x k_j table; which = 1 the score dispersion
 * sum_a n_a (vbar_a - vbar)^2 with scores 0..k_j-1, the Kruskal-Wallis analogue built from the
 * group means themselves. Pearson standardises each item by its own variance, which is what one
 * wants when SUMMING across items of different prevalence; the dispersion keeps the items on the
 * score scale, which is what one wants when the question is about a common shift.
 *
 * x: n * m ints, item-major (column-major), values 0..k_j - 1.  grp: n ints in 0..g-1.
 * cl: n ints in 0..ncl-1, ignored when ncl == 0.  out: 12 doubles, obs_sum, obs_max, then for
 * the sum p, midp, lo, hi, then the same four for the max, then the largest null sum and the
 * largest null max.  item_stat: m observed per-item statistics.  item_p: m per-item p-values. */
static double ec_item_stat(const int *xj, const int *grp, int n, int g, int k,
                           const int *ng, int *cnt, int which){
  int a, l; double N = (double) n, T = 0.0;
  for (a = 0; a < g * k; a++) cnt[a] = 0;
  for (a = 0; a < n; a++) cnt[grp[a] * k + xj[a]]++;
  if (which == 0){
    for (l = 0; l < k; l++){
      double cs = 0.0; for (a = 0; a < g; a++) cs += cnt[a * k + l];
      if (cs <= 0.0) continue;
      for (a = 0; a < g; a++){
        double e = ng[a] * cs / N; if (e <= 0.0) continue;
        double d = cnt[a * k + l] - e; T += d * d / e; }
    }
    return T;
  }
  { double gm = 0.0;
    for (a = 0; a < g; a++) for (l = 0; l < k; l++) gm += l * (double) cnt[a * k + l];
    gm /= N;
    for (a = 0; a < g; a++){ double s = 0.0;
      if (ng[a] <= 0) continue;
      for (l = 0; l < k; l++) s += l * (double) cnt[a * k + l];
      s = s / ng[a] - gm; T += ng[a] * s * s; }
    return T; }
}
static void ec_multi_all(const int *x, const int *grp, int n, int m, int g, const int *k,
                         int kmax, const int *ng, int *cnt, int which, double *st){
  int j; for (j = 0; j < m; j++) st[j] = ec_item_stat(x + (size_t) j * n, grp, n, g, k[j], ng, cnt, which);
  (void) kmax;
}
void C_ec_multi(const int *x, const int *n_, const int *m_, const int *g_, const int *grp,
                const int *k, const int *cl, const int *ncl_, const int *which, const int *B,
                const int *seed, double *out, double *item_stat, double *item_p){
  int n = *n_, m = *m_, g = *g_, ncl = *ncl_, nb = *B, w = *which, i, j, kmax = 0;
  int *ng = (int*) calloc((size_t) g, sizeof(int));
  int *clg = NULL;                                   /* group of each cluster */
  double *nsum = (double*) malloc((size_t) nb * sizeof(double));
  double *nmax = (double*) malloc((size_t) nb * sizeof(double));
  int *ge_item = (int*) calloc((size_t) m, sizeof(int));
  double t4[4];
  if (!ng || !nsum || !nmax || !ge_item){ free(ng); free(nsum); free(nmax); free(ge_item); return; }
  for (j = 0; j < m; j++) if (k[j] > kmax) kmax = k[j];
  for (i = 0; i < n; i++) ng[grp[i]]++;
  { int *cnt = (int*) malloc((size_t)(g * kmax) * sizeof(int));
    if (!cnt){ free(ng); free(nsum); free(nmax); free(ge_item); return; }
    ec_multi_all(x, grp, n, m, g, k, kmax, ng, cnt, w, item_stat);
    free(cnt); }
  out[0] = 0.0; out[1] = 0.0;
  for (j = 0; j < m; j++){ out[0] += item_stat[j]; if (j == 0 || item_stat[j] > out[1]) out[1] = item_stat[j]; }
  if (ncl > 0){
    clg = (int*) malloc((size_t) ncl * sizeof(int));
    if (!clg){ free(ng); free(nsum); free(nmax); free(ge_item); return; }
    for (i = 0; i < ncl; i++) clg[i] = -1;
    for (i = 0; i < n; i++) clg[cl[i]] = grp[i];   /* the R front has checked nesting */
  }
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (i = 0; i < nb; i++){
    unsigned long long s = ec_seed((unsigned long long) *seed, i);
    int *pg = (int*) malloc((size_t) n * sizeof(int));
    int *cnt = (int*) malloc((size_t)(g * kmax) * sizeof(int));
    int *pn = (int*) calloc((size_t) g, sizeof(int));
    double *st = (double*) malloc((size_t) m * sizeof(double));
    int a; double sm = 0.0, mx = 0.0;
    if (!pg || !cnt || !pn || !st){ free(pg); free(cnt); free(pn); free(st); nsum[i] = 0.0; nmax[i] = 0.0; continue; }
    if (ncl > 0){
      int *pc = (int*) malloc((size_t) ncl * sizeof(int));
      if (!pc){ free(pg); free(cnt); free(pn); free(st); nsum[i] = 0.0; nmax[i] = 0.0; continue; }
      memcpy(pc, clg, (size_t) ncl * sizeof(int));
      ec_shuffle(pc, ncl, &s);                       /* permute labels across clusters */
      for (a = 0; a < n; a++) pg[a] = pc[cl[a]];
      free(pc);
    } else {
      memcpy(pg, grp, (size_t) n * sizeof(int));
      ec_shuffle(pg, n, &s);                         /* permute labels across respondents */
    }
    for (a = 0; a < n; a++) pn[pg[a]]++;
    ec_multi_all(x, pg, n, m, g, k, kmax, pn, cnt, w, st);
    for (a = 0; a < m; a++){ sm += st[a]; if (a == 0 || st[a] > mx) mx = st[a]; }
    nsum[i] = sm; nmax[i] = mx;
#ifdef _OPENMP
#pragma omp critical
#endif
    { for (a = 0; a < m; a++) if (st[a] >= item_stat[a] - 1e-12) ge_item[a]++; }
    free(pg); free(cnt); free(pn); free(st);
  }
  ec_tail(out[0], nsum, nb, t4); out[2] = t4[0]; out[3] = t4[1]; out[4] = t4[2]; out[5] = t4[3];
  ec_tail(out[1], nmax, nb, t4); out[6] = t4[0]; out[7] = t4[1]; out[8] = t4[2]; out[9] = t4[3];
  out[10] = nsum[0]; out[11] = nmax[0];
  for (i = 1; i < nb; i++){ if (nsum[i] > out[10]) out[10] = nsum[i]; if (nmax[i] > out[11]) out[11] = nmax[i]; }
  for (j = 0; j < m; j++) item_p[j] = (1.0 + ge_item[j]) / (nb + 1.0);
  free(ng); free(clg); free(nsum); free(nmax); free(ge_item);
}

/* --------------------------------------- a DIFFERENTIAL NETWORK without a Gaussian
 * THE GAP THIS FILLS. A differential network asks which pairs of variables are associated
 * differently in two groups. The established route estimates a precision matrix per group and
 * compares partial correlations through a Fisher-z difference, whose null is chi-square BECAUSE
 * the data are assumed Gaussian; for counts or ordinal survey items that assumption is false, and
 * the usual repair is to make the data Gaussian first (a latent field, a whitening step) rather
 * than to drop the assumption.
 *
 * Dropping it costs nothing here. Take any edge statistic that measures how far the groups'
 * associations differ -- this uses the heterogeneity of the Fisher-transformed score correlations,
 * sum_a (n_a - 3)(z_a - zbar)^2, which is the same quantity the Gaussian theory uses -- and refer
 * it not to a chi-square but to the law obtained by PERMUTING THE GROUP LABELS. That law conditions
 * on the observed data and on the group sizes, so it is exact whatever the marginal distributions
 * are: no Gaussianity, no large sample, no correction for ties. What the Gaussian assumption bought
 * was a closed-form null; a permutation supplies the same null and nothing else changes.
 *
 * MULTIPLICITY, which is the real obstacle in a network. There are p(p-1)/2 edges, so a per-edge
 * test needs a correction that costs most of the power. The same permutations give ONE statistic
 * over the whole network -- the SUM of the edge statistics (many edges differing a little) or the
 * MAXIMUM (one edge differing a lot) -- and so ONE p-value for "is the network different at all",
 * with the per-edge p-values from those same permutations available for a closed follow-up once it
 * has rejected. Nothing is corrected because nothing is repeated.
 *
 * CLUSTERS. If respondents sit in clusters, the labels are permuted across CLUSTERS, each keeping
 * its members, exactly as in C_ec_multi. That is the exact treatment of an interviewer effect.
 *
 * x: n * p doubles, variable-major: ranks for continuous variables (which makes the correlation
 * Spearman's and the whole procedure distribution-free), integer codes for ordinal ones, 0/1 for
 * binary. An affine rescoring changes nothing. out: 12 doubles as in C_ec_multi. edge_stat, edge_p: m each,
 * m = p(p-1)/2, in the order (0,1), (0,2), ..., (0,p-1), (1,2), ... */
static void ec_edge_all(const double *x, const int *grp, int n, int p, int g, int m,
                        double *su, double *suu, double *suv, double *ng, double *st){
  int i, j, l, a, e;
  for (a = 0; a < g; a++){ ng[a] = 0.0;
    for (j = 0; j < p; j++){ su[a * p + j] = 0.0; suu[a * p + j] = 0.0; }
    for (e = 0; e < m; e++) suv[a * m + e] = 0.0; }
  for (i = 0; i < n; i++){
    a = grp[i]; ng[a] += 1.0;
    { double *sua = su + a * p, *suua = suu + a * p, *suva = suv + a * m;
      e = 0;
      for (j = 0; j < p; j++){
        double xj = x[(size_t) j * n + i];
        sua[j] += xj; suua[j] += xj * xj;
        for (l = j + 1; l < p; l++) suva[e++] += xj * x[(size_t) l * n + i];
      } }
  }
  e = 0;
  for (j = 0; j < p; j++) for (l = j + 1; l < p; l++){
    double zb = 0.0, w = 0.0, T = 0.0;
    double zs[64]; double ws[64]; int a2, gg = g < 64 ? g : 64;
    for (a = 0; a < gg; a++){
      double na = ng[a], r, den, z;
      double s1 = su[a * p + j], s2 = su[a * p + l];
      double q1 = suu[a * p + j] - s1 * s1 / (na > 0 ? na : 1.0);
      double q2 = suu[a * p + l] - s2 * s2 / (na > 0 ? na : 1.0);
      double c12 = suv[a * m + e] - s1 * s2 / (na > 0 ? na : 1.0);
      den = sqrt(q1 * q2);
      r = (na > 3.0 && den > 1e-12) ? c12 / den : 0.0;
      if (r > 0.999999) r = 0.999999; if (r < -0.999999) r = -0.999999;
      z = 0.5 * log((1.0 + r) / (1.0 - r));
      zs[a] = z; ws[a] = na > 3.0 ? na - 3.0 : 0.0;
      zb += ws[a] * z; w += ws[a];
    }
    if (w > 0.0){ zb /= w;
      for (a2 = 0; a2 < gg; a2++){ double d = zs[a2] - zb; T += ws[a2] * d * d; } }
    st[e] = T; e++;
  }
}
void C_ec_diffnet(const double *x, const int *n_, const int *p_, const int *g_, const int *grp,
                  const int *cl, const int *ncl_, const int *B, const int *seed,
                  double *out, double *edge_stat, double *edge_p){
  int n = *n_, p = *p_, g = *g_, ncl = *ncl_, nb = *B, i, j, m = p * (p - 1) / 2;
  int *clg = NULL;
  double *nsum = (double*) malloc((size_t) nb * sizeof(double));
  double *nmax = (double*) malloc((size_t) nb * sizeof(double));
  int *ge_edge = (int*) calloc((size_t) m, sizeof(int));
  double t4[4];
  if (p < 2 || g < 2 || g > 64){ free(nsum); free(nmax); free(ge_edge); return; }
  if (!nsum || !nmax || !ge_edge){ free(nsum); free(nmax); free(ge_edge); return; }
  { double *su = (double*) malloc((size_t)(g * p) * sizeof(double));
    double *suu = (double*) malloc((size_t)(g * p) * sizeof(double));
    double *suv = (double*) malloc((size_t)(g * m) * sizeof(double));
    double *ng = (double*) malloc((size_t) g * sizeof(double));
    if (!su || !suu || !suv || !ng){ free(su); free(suu); free(suv); free(ng);
      free(nsum); free(nmax); free(ge_edge); return; }
    ec_edge_all(x, grp, n, p, g, m, su, suu, suv, ng, edge_stat);
    free(su); free(suu); free(suv); free(ng); }
  out[0] = 0.0; out[1] = 0.0;
  for (j = 0; j < m; j++){ out[0] += edge_stat[j]; if (j == 0 || edge_stat[j] > out[1]) out[1] = edge_stat[j]; }
  if (ncl > 0){
    clg = (int*) malloc((size_t) ncl * sizeof(int));
    if (!clg){ free(nsum); free(nmax); free(ge_edge); return; }
    for (i = 0; i < ncl; i++) clg[i] = -1;
    for (i = 0; i < n; i++) clg[cl[i]] = grp[i];
  }
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (i = 0; i < nb; i++){
    unsigned long long s = ec_seed((unsigned long long) *seed, i);
    int *pg = (int*) malloc((size_t) n * sizeof(int));
    double *su = (double*) malloc((size_t)(g * p) * sizeof(double));
    double *suu = (double*) malloc((size_t)(g * p) * sizeof(double));
    double *suv = (double*) malloc((size_t)(g * m) * sizeof(double));
    double *ng = (double*) malloc((size_t) g * sizeof(double));
    double *st = (double*) malloc((size_t) m * sizeof(double));
    int a; double sm = 0.0, mx = 0.0;
    if (!pg || !su || !suu || !suv || !ng || !st){
      free(pg); free(su); free(suu); free(suv); free(ng); free(st);
      nsum[i] = 0.0; nmax[i] = 0.0; continue; }
    if (ncl > 0){
      int *pc = (int*) malloc((size_t) ncl * sizeof(int));
      if (!pc){ free(pg); free(su); free(suu); free(suv); free(ng); free(st);
        nsum[i] = 0.0; nmax[i] = 0.0; continue; }
      memcpy(pc, clg, (size_t) ncl * sizeof(int));
      ec_shuffle(pc, ncl, &s);
      for (a = 0; a < n; a++) pg[a] = pc[cl[a]];
      free(pc);
    } else {
      memcpy(pg, grp, (size_t) n * sizeof(int));
      ec_shuffle(pg, n, &s);
    }
    ec_edge_all(x, pg, n, p, g, m, su, suu, suv, ng, st);
    for (a = 0; a < m; a++){ sm += st[a]; if (a == 0 || st[a] > mx) mx = st[a]; }
    nsum[i] = sm; nmax[i] = mx;
#ifdef _OPENMP
#pragma omp critical
#endif
    { for (a = 0; a < m; a++) if (st[a] >= edge_stat[a] - 1e-12) ge_edge[a]++; }
    free(pg); free(su); free(suu); free(suv); free(ng); free(st);
  }
  ec_tail(out[0], nsum, nb, t4); out[2] = t4[0]; out[3] = t4[1]; out[4] = t4[2]; out[5] = t4[3];
  ec_tail(out[1], nmax, nb, t4); out[6] = t4[0]; out[7] = t4[1]; out[8] = t4[2]; out[9] = t4[3];
  out[10] = nsum[0]; out[11] = nmax[0];
  for (i = 1; i < nb; i++){ if (nsum[i] > out[10]) out[10] = nsum[i]; if (nmax[i] > out[11]) out[11] = nmax[i]; }
  for (j = 0; j < m; j++) edge_p[j] = (1.0 + ge_edge[j]) / (nb + 1.0);
  free(clg); free(nsum); free(nmax); free(ge_edge);
}

/* ------------------------------- two tables, conditioning on the total score
 * THE ORDINAL EDGE LAW, and why no tilting or reweighting is needed. Under H0: phi1 = phi2 = phi
 * the joint density of the pair, conditional on both margins WITHIN each table, is
 *
 *     pi(x1, x2) proportional to  [prod 1/x1_ij!] [prod 1/x2_ij!] exp(phi * (T1 + T2)),
 *
 * so conditioning FURTHER on S = T1 + T2 makes the exponential factor exp(phi * S) a constant.
 * The common association parameter does not merely cancel in a ratio, it leaves the density
 * altogether: what remains is a distribution over pairs of tables that carries no unknown
 * parameter at all, exactly as the odds ratio drops out of Zelen's two by two construction.
 *
 * That law is sampled here directly. A move is a 2 x 2 swap in EACH table whose score changes are
 * equal and opposite, so both sets of margins and the total score S are all preserved; 2 x 2 swaps
 * generate the whole fibre of tables with given margins, so the chain reaches every attainable
 * pair. The chain starts AT the observed pair, which is what makes this robust where the obvious
 * construction is not: estimating each table's null pmf separately and convolving them requires
 * the observed score to be reachable under NO association, and for a strongly associated table it
 * simply is not -- measured here at five standard deviations out, where two million independent
 * draws never once reach the observed value. Worse, that failure occurs precisely when BOTH tables
 * are strongly associated, which is the case in which the honest answer is a LARGE p-value, so a
 * convolution estimator fails in the direction that would be read as evidence of a difference.
 *
 * The price is that successive draws are dependent, so the accuracy here is governed by the chain
 * rather than by independent sampling; the acceptance rate is returned so that it can be judged.
 * out: midp, attainable size at alpha, T1, S, acceptance rate. */
static double ec_ll_score(const double *u, const double *v, int i1, int i2, int j1, int j2){
  return (u[i1] - u[i2]) * (v[j2] - v[j1]);
}
void C_ec_rc_two(const int *x1, const int *x2, const int *r, const int *c, const int *B,
                 const int *burn, const int *thin, const double *u, const double *v,
                 const double *alpha, const int *seed, double *out){
  int R = *r, C = *c, nb = *B, nburn = *burn, nthin = *thin, i, j, k, it;
  int *t1 = NULL, *t2 = NULL, *hist = NULL, *ord = NULL;
  double *Tdraw = NULL, *prob = NULL, *mp = NULL;
  long long acc = 0, tries = 0;
  int T1, T2, S, lo, hi, nat, obs;
  double tot = 0.0, cw = 0.0;
  unsigned long long st = ec_seed((unsigned long long) *seed, 0);
  if (R < 2 || C < 2){ out[0] = -1.0; return; }
  t1 = (int*) malloc((size_t)(R * C) * sizeof(int));
  t2 = (int*) malloc((size_t)(R * C) * sizeof(int));
  Tdraw = (double*) malloc((size_t) nb * sizeof(double));
  if (!t1 || !t2 || !Tdraw){ free(t1); free(t2); free(Tdraw); out[0] = -1.0; return; }
  memcpy(t1, x1, (size_t)(R * C) * sizeof(int));
  memcpy(t2, x2, (size_t)(R * C) * sizeof(int));
  T1 = (int) llround(ec_stat(x1, R, C, u, v, 0));
  T2 = (int) llround(ec_stat(x2, R, C, u, v, 0));
  S = T1 + T2;
  { int cur = T1;
    for (it = 0; it < nburn + nb * nthin; it++){
      int a1, b1, p1, q1, a2, b2, p2, q2; double d1, d2, ratio;
      a1 = (int) ec_below(&st, (unsigned long long) R);
      do { b1 = (int) ec_below(&st, (unsigned long long) R); } while (b1 == a1);
      p1 = (int) ec_below(&st, (unsigned long long) C);
      do { q1 = (int) ec_below(&st, (unsigned long long) C); } while (q1 == p1);
      a2 = (int) ec_below(&st, (unsigned long long) R);
      do { b2 = (int) ec_below(&st, (unsigned long long) R); } while (b2 == a2);
      p2 = (int) ec_below(&st, (unsigned long long) C);
      do { q2 = (int) ec_below(&st, (unsigned long long) C); } while (q2 == p2);
      tries++;
      d1 = ec_ll_score(u, v, a1, b1, p1, q1);
      d2 = ec_ll_score(u, v, a2, b2, p2, q2);
      if (fabs(d1 + d2) > 1e-9) goto record;              /* S would move: not in the fibre */
      if (t1[a1*C+p1] < 1 || t1[b1*C+q1] < 1) goto record;
      if (t2[a2*C+p2] < 1 || t2[b2*C+q2] < 1) goto record;
      ratio = ((double) t1[a1*C+p1] * t1[b1*C+q1]) /
              (((double) t1[a1*C+q1] + 1.0) * ((double) t1[b1*C+p1] + 1.0));
      ratio *= ((double) t2[a2*C+p2] * t2[b2*C+q2]) /
               (((double) t2[a2*C+q2] + 1.0) * ((double) t2[b2*C+p2] + 1.0));
      if (ratio >= 1.0 || ec_unif(&st) < ratio){
        t1[a1*C+p1]--; t1[b1*C+q1]--; t1[a1*C+q1]++; t1[b1*C+p1]++;
        t2[a2*C+p2]--; t2[b2*C+q2]--; t2[a2*C+q2]++; t2[b2*C+p2]++;
        cur += (int) llround(d1);
        acc++;
      }
record:
      if (it >= nburn && ((it - nburn) % nthin) == 0){
        int idx = (it - nburn) / nthin;
        if (idx < nb) Tdraw[idx] = (double) cur;
      }
    }
  }
  lo = hi = (int) Tdraw[0];
  for (i = 0; i < nb; i++){ int t = (int) Tdraw[i]; if (t < lo) lo = t; if (t > hi) hi = t; }
  if (T1 < lo) lo = T1; if (T1 > hi) hi = T1;
  nat = hi - lo + 1;
  hist = (int*) calloc((size_t) nat, sizeof(int));
  prob = (double*) calloc((size_t) nat, sizeof(double));
  ord = (int*) malloc((size_t) nat * sizeof(int));
  mp = (double*) calloc((size_t) nat, sizeof(double));
  if (!hist || !prob || !ord || !mp){ free(hist); free(prob); free(ord); free(mp);
    free(t1); free(t2); free(Tdraw); out[0] = -1.0; return; }
  for (i = 0; i < nb; i++) hist[(int) Tdraw[i] - lo]++;
  for (k = 0; k < nat; k++){ prob[k] = (double) hist[k] / nb; tot += prob[k]; }
  for (k = 0; k < nat; k++) prob[k] /= tot;
  for (k = 0; k < nat; k++) ord[k] = k;
  for (i = 1; i < nat; i++){ int t = ord[i]; j = i - 1;             /* insertion sort, prob asc */
    while (j >= 0 && prob[ord[j]] > prob[t]){ ord[j+1] = ord[j]; j--; } ord[j+1] = t; }
  for (i = 0; i < nat; i++){ cw += prob[ord[i]]; mp[ord[i]] = cw - 0.5 * prob[ord[i]]; }
  obs = T1 - lo;
  out[0] = (obs >= 0 && obs < nat) ? mp[obs] : 0.0;
  out[1] = 0.0;
  for (k = 0; k < nat; k++) if (mp[k] <= *alpha) out[1] += prob[k];
  out[2] = T1; out[3] = S;
  out[4] = tries > 0 ? (double) acc / (double) tries : 0.0;
  free(hist); free(prob); free(ord); free(mp); free(t1); free(t2); free(Tdraw);
}

/* ------------------------- uniform direction on the sphere of a subspace
 * M is the n x n residual maker (symmetric idempotent); A holds m symmetric n x n matrices stacked
 * row-major. Each draw returns max_k u' A_k u with u uniform on the unit sphere of the column space
 * of M. The maximum is taken INSIDE the draw, which is the whole point: it is the distribution of
 * the statistic a practitioner actually reports after choosing the graph by fit, and it is exact
 * for ANY selection rule that is a function of u. */
void C_ec_sphere_maxquad(const double *M, const double *A, const int *n, const int *m,
                         const int *B, const int *seed, double *out){
  int N = *n, K = *m, nb = *B, i;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (i = 0; i < nb; i++){
    unsigned long long s = ec_seed((unsigned long long) *seed, i);
    double *z = (double*) malloc((size_t) N * sizeof(double));
    double *e = (double*) malloc((size_t) N * sizeof(double));
    int a, b, k; double ss = 0.0, best = 0.0;
    if (!z || !e){ free(z); free(e); out[i] = 0.0; continue; }
    for (a = 0; a < N; a++) z[a] = ec_norm(&s);
    for (a = 0; a < N; a++){ double t = 0.0;
      for (b = 0; b < N; b++) t += M[a * N + b] * z[b];
      e[a] = t; ss += t * t; }
    ss = sqrt(ss);
    if (ss > 0.0) for (a = 0; a < N; a++) e[a] /= ss;
    for (k = 0; k < K; k++){
      const double *Ak = A + (size_t) k * N * N; double q = 0.0;
      for (a = 0; a < N; a++){ double t = 0.0;
        for (b = 0; b < N; b++) t += Ak[a * N + b] * e[b];
        q += e[a] * t; }
      if (k == 0 || q > best) best = q;
    }
    out[i] = best;
    free(z); free(e);
  }
}

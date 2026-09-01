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

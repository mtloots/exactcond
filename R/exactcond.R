## exactcond: exact conditional reference distributions. R front to the shared C back-end.

.EC_STATS <- c(linlin = 0L, dispersion = 1L, pearson = 2L, lr = 3L)

#' Integer scores for an ordinal margin
#'
#' @param k number of categories.
#' @param centre if TRUE, centre the scores on zero.
#' @return A numeric vector of length \code{k}.
#' @details The dispersion and linear-by-linear statistics are unchanged in distribution by an
#'   affine transformation of the scores, so \code{0:(k-1)} and a centred version give identical
#'   p-values; \code{centre} exists only to make output easier to read.
#' @examples ec_scores(3)
#' @export
ec_scores <- function(k, centre = FALSE) {
  s <- seq_len(k) - 1
  if (centre) s - mean(s) else as.numeric(s)
}

ec_check_scores <- function(u, v, r, c) {
  if (is.null(u)) u <- ec_scores(r)
  if (is.null(v)) v <- ec_scores(c)
  if (length(u) != r) stop("u must have one score per row, so length ", r)
  if (length(v) != c) stop("v must have one score per column, so length ", c)
  list(u = as.double(u), v = as.double(v))
}

#' Exact conditional test in a two-way table
#'
#' Conditions on BOTH margins, which are sufficient for the marginal parameters, and so leaves a
#' null distribution carrying no nuisance parameter at all. That distribution is sampled exactly,
#' by random permutation of the column labels; nothing here is asymptotic.
#'
#' @param x an r by c table or integer matrix of counts.
#' @param stat one of \code{"dispersion"}, \code{"linlin"}, \code{"pearson"}, \code{"lr"}.
#' @param u,v scores for the rows and columns; defaults \code{0:(r-1)} and \code{0:(c-1)}.
#'   \code{"dispersion"} ignores \code{u}, the row variable being treated as nominal.
#' @param B number of draws from the conditional law.
#' @param seed integer seed. Each replicate uses its own stream, so the result does not depend on
#'   the number of threads.
#' @return A list with the observed statistic, the exact conditional p-value and its mid-p
#'   companion, a Wilson interval for the p-value reflecting the Monte Carlo error in \code{B}
#'   draws, the largest null draw seen, and the components of the statistic.
#' @details \code{"dispersion"} is for NOMINAL rows against an ORDINAL column and is
#'   \eqn{\sum_i n_i (\bar v_i - \bar v)^2}. Its components are the row mean scores themselves, so
#'   in an application where the column score is a signed difference the statistic is built out of
#'   the very contrasts being reported, rather than out of ranks. \code{"linlin"} is Goodman's
#'   \eqn{\sum_{ij} u_i v_j x_{ij}}, sufficient for the linear-by-linear association parameter given
#'   both margins.
#'
#'   When \code{p.value} is at \code{1/(B+1)} the observed statistic exceeded every draw, so the
#'   test has hit its resolution floor and the honest report is an upper bound together with
#'   \code{max_null}, not a p-value.
#' @examples
#' tab <- matrix(c(3, 228, 72, 10, 548, 67, 7, 276, 131), nrow = 3, byrow = TRUE)
#' ec_rc_test(tab, v = c(-1, 0, 1), B = 2000, seed = 1)
#' @export
ec_rc_test <- function(x, stat = c("dispersion", "linlin", "pearson", "lr"),
                       u = NULL, v = NULL, B = 20000L, seed = 1L) {
  stat <- match.arg(stat)
  x <- as.matrix(x); storage.mode(x) <- "integer"
  if (any(is.na(x)) || any(x < 0L)) stop("x must hold non-negative counts with no NA")
  r <- nrow(x); c <- ncol(x)
  s <- ec_check_scores(u, v, r, c)
  z <- .C("C_ec_rc_test", as.integer(t(x)), as.integer(r), as.integer(c), as.integer(B),
          as.integer(.EC_STATS[[stat]]), s$u, s$v, as.integer(seed), out = double(6))$out
  n <- rowSums(x)
  structure(list(statistic = z[1], p.value = z[2], midp = z[3],
                 mc.lower = z[4], mc.upper = z[5], max_null = z[6],
                 row_means = as.vector(x %*% s$v) / n, n = n, B = B, stat = stat, table = x),
            class = "ec_rc")
}

#' @export
print.ec_rc <- function(x, ...) {
  cat("Exact conditional test in a two-way table, both margins fixed\n")
  cat(sprintf("statistic (%s) = %.6g\n", x$stat, x$statistic))
  fl <- 1 / (x$B + 1)
  if (x$p.value <= fl + 1e-12)
    cat(sprintf("p < %.3g (resolution floor); largest of %s null draws = %.6g\n",
                fl, format(x$B, big.mark = ",", scientific = FALSE), x$max_null))
  else
    cat(sprintf("p = %.4g   mid-p = %.4g   Monte Carlo 95%% interval [%.4g, %.4g]\n",
                x$p.value, x$midp, x$mc.lower, x$mc.upper))
  if (x$stat == "dispersion") {
    cat("row mean scores: ",
        paste(sprintf("%.4g", x$row_means), collapse = ", "), "\n", sep = "")
  }
  invisible(x)
}

#' Null draws of a table statistic given both margins
#'
#' @param rows,cols the fixed row and column totals.
#' @inheritParams ec_rc_test
#' @return A numeric vector of \code{B} draws.
#' @examples ec_rc_null(c(30, 40), c(35, 35), B = 100, seed = 2)
#' @export
ec_rc_null <- function(rows, cols, stat = c("dispersion", "linlin", "pearson", "lr"),
                       u = NULL, v = NULL, B = 20000L, seed = 1L) {
  stat <- match.arg(stat)
  rows <- as.integer(rows); cols <- as.integer(cols)
  if (sum(rows) != sum(cols)) stop("the row and column totals must agree")
  s <- ec_check_scores(u, v, length(rows), length(cols))
  .C("C_ec_rc_null", rows, cols, as.integer(length(rows)), as.integer(length(cols)),
     as.integer(B), as.integer(.EC_STATS[[stat]]), s$u, s$v, as.integer(seed),
     out = double(B))$out
}

#' Exact conditional test that two tables share an association parameter
#'
#' Conditions on both margins WITHIN each table and then on the total score across them. Under the
#' null of a common linear-by-linear parameter the joint density carries that parameter only through
#' the total score, which the conditioning has fixed, so it leaves the density altogether: what
#' remains carries no unknown parameter at all. Nothing is exchanged between the groups, so the test
#' stays valid when the two tables have very different margins, which is exactly where permuting a
#' group label fails.
#'
#' @param t1,t2 the two tables, of the same dimension.
#' @param u,v integer scores; required to be integer valued so the total score lies on a lattice.
#' @param alpha level at which the attainable size is reported.
#' @param B number of recorded draws.
#' @param burn moves discarded before recording begins.
#' @param thin moves between recorded draws.
#' @param seed integer seed.
#' @return A list with the conditional mid-p value, the attainable size at \code{alpha}, the two
#'   conditioning statistics, and the acceptance rate of the sampler.
#' @details The law is sampled by moves that make a two by two swap in EACH table with equal and
#'   opposite score changes, so both sets of margins and the total score are preserved, and the
#'   chain starts at the observed pair.
#'
#'   That last point is what makes this usable. The obvious construction instead estimates each
#'   table's null distribution separately and convolves them, which requires the observed score to
#'   be reachable under NO association; for a strongly associated table it is not, and no increase
#'   in \code{B} repairs it. Worse, the failure arrives precisely when BOTH tables are strongly
#'   associated, which is the case where the correct answer is a large p-value, so the convolution
#'   fails in the direction a reader would take for evidence of a difference.
#'
#'   Successive draws are dependent, so accuracy is governed by the chain rather than by
#'   independent sampling. Check \code{accept} and raise \code{burn} or \code{thin} if it is low.
#' @examples
#' a <- matrix(c(20, 8, 3, 9, 22, 7, 4, 8, 19), 3, byrow = TRUE)
#' b <- matrix(c(18, 9, 4, 8, 20, 8, 5, 9, 18), 3, byrow = TRUE)
#' ec_rc_two(a, b, B = 2000, seed = 3)
#' @export
ec_rc_two <- function(t1, t2, u = NULL, v = NULL, B = 20000L, burn = 20000L, thin = 20L,
                      alpha = 0.05, seed = 1L) {
  t1 <- as.matrix(t1); t2 <- as.matrix(t2)
  storage.mode(t1) <- "integer"; storage.mode(t2) <- "integer"
  if (!identical(dim(t1), dim(t2))) stop("the two tables must have the same dimension")
  r <- nrow(t1); c <- ncol(t1)
  if (r < 2 || c < 2) stop("both tables must have at least two rows and two columns")
  s <- ec_check_scores(u, v, r, c)
  if (any(abs(s$u - round(s$u)) > 1e-9) || any(abs(s$v - round(s$v)) > 1e-9))
    stop("ec_rc_two needs integer scores, so that the total score lies on a lattice")
  z <- .C("C_ec_rc_two", as.integer(t(t1)), as.integer(t(t2)), as.integer(r), as.integer(c),
          as.integer(B), as.integer(burn), as.integer(thin), s$u, s$v, as.double(alpha),
          as.integer(seed), out = double(5))$out
  if (z[1] < 0) stop("the sampler could not be initialised on these tables")
  list(midp = z[1], size = z[2], T1 = z[3], S = z[4], accept = z[5],
       B = B, burn = burn, thin = thin, alpha = alpha)
}

#' Exact null for a statistic maximised over candidate structures
#'
#' Under a null in which the tested structure does not enter the data-generating process, a
#' statistic invariant to the regression coefficients and to the error scale depends on the data
#' only through the residual direction, which is uniform on the unit sphere of the residual space
#' whatever those nuisance parameters are. Maximising over candidates inside each draw gives the
#' exact null of the statistic a practitioner actually reports after selecting a structure by fit,
#' and it is valid for ANY selection rule that is a function of that direction.
#'
#' @param M the residual maker, a symmetric idempotent n by n matrix.
#' @param A a single n by n symmetric matrix, or a list of them, one per candidate.
#' @inheritParams ec_rc_test
#' @return A numeric vector of \code{B} draws of the maximised quadratic form.
#' @examples
#' n <- 12; X <- cbind(1, seq_len(n))
#' M <- diag(n) - X %*% solve(crossprod(X), t(X))
#' A <- diag(n)[, c(n, seq_len(n - 1))]; A <- (A + t(A)) / 2
#' quantile(ec_sphere_maxquad(M, A, B = 500, seed = 4), 0.95)
#' @export
ec_sphere_maxquad <- function(M, A, B = 20000L, seed = 1L) {
  M <- as.matrix(M); n <- nrow(M)
  if (ncol(M) != n) stop("M must be square")
  if (!is.list(A)) A <- list(A)
  A <- lapply(A, function(a) { a <- as.matrix(a)
    if (!identical(dim(a), c(n, n))) stop("every candidate must be n by n"); t(a) })
  .C("C_ec_sphere_maxquad", as.double(t(M)), as.double(unlist(A)), as.integer(n),
     as.integer(length(A)), as.integer(B), as.integer(seed), out = double(B))$out
}

#' One exact test for a battery of items across groups
#'
#' Several categorical items are recorded on the same respondents, who fall into groups, and the
#' question is whether ANY item differs across the groups. Rather than test each item and pay a
#' multiplicity correction, one statistic is formed, the sum or the maximum over items of a
#' per-item statistic, and referred to the law obtained by permuting the group labels. That law
#' conditions on every item's margin and on the group sizes, so it is exact conditional in the
#' sense of \code{\link{ec_rc_test}}, and it gives ONE p-value for the whole battery. The per-item
#' statistics are referred to the same permutations, which supplies the unadjusted per-item
#' p-values a closed procedure uses once the omnibus has rejected.
#'
#' @param x a matrix or data frame of items, one column per item, each a factor or integer
#'   coded from zero (or from one, which is shifted), with respondents in rows.
#' @param group a factor or vector giving each respondent's group.
#' @param cluster an optional factor or vector giving each respondent's cluster, every cluster
#'   lying inside one group. When supplied the labels are permuted across clusters, each cluster
#'   keeping its members, which is the exact treatment of intra-cluster correlation such as an
#'   interviewer or fieldworker effect.
#' @param stat \code{"pearson"}, the Pearson statistic of each item's group-by-level table, which
#'   standardises each item by its own variance and is the natural choice for the sum; or
#'   \code{"dispersion"}, the dispersion of the group mean scores with scores \code{0:(k-1)},
#'   the Kruskal-Wallis analogue on the score scale.
#' @param B number of permutations.
#' @param seed integer seed; each permutation uses its own stream, so the result does not depend
#'   on the thread count.
#' @return A list with \code{sum} and \code{max}, each holding the observed statistic, the exact
#'   p-value, its mid-p companion, a Wilson interval for the Monte Carlo error and the largest
#'   null draw; \code{items}, a data frame of the per-item statistics and unadjusted per-item
#'   p-values from the same permutations; and the sizes of the groups and, if given, the number
#'   of clusters.
#' @details The sum has power against many items moving a little; the maximum against one item
#'   moving a lot, and it is the one that matches or beats a Bonferroni-corrected battery while
#'   holding its size. With one binary item and no clusters, \code{"pearson"} is the exact test of
#'   the two-way table and agrees with \code{\link{ec_rc_test}} in distribution.
#' @examples
#' set.seed(1)
#' g <- rep(c("a", "b", "c"), c(30, 40, 30))
#' x <- cbind(rbinom(100, 1, 0.9), rbinom(100, 1, 0.95), rbinom(100, 1, ifelse(g == "c", 0.7, 0.9)))
#' ec_multi_test(x, g, B = 2000, seed = 1)
#' @export
ec_multi_test <- function(x, group, cluster = NULL, stat = c("pearson", "dispersion"),
                          B = 20000L, seed = 1L) {
  stat <- match.arg(stat)
  x <- as.data.frame(x); n <- nrow(x); m <- ncol(x)
  if (m < 1) stop("x must hold at least one item")
  if (length(group) != n) stop("group must have one entry per row of x")
  cols <- lapply(x, function(v) {
    if (is.factor(v)) v <- as.integer(v) - 1L else if (is.logical(v)) v <- as.integer(v)
    else { v <- as.integer(v); if (min(v, na.rm = TRUE) == 1L) v <- v - 1L }
    v })
  if (any(vapply(cols, function(v) any(is.na(v)), TRUE))) stop("x must have no missing values")
  k <- vapply(cols, function(v) max(v) + 1L, 1L)
  if (any(vapply(cols, min, 1L) < 0L)) stop("items must be coded from zero or one")
  group <- as.factor(group); grp <- as.integer(group) - 1L; g <- nlevels(group)
  if (g < 2) stop("group must have at least two levels")
  if (is.null(cluster)) { cl <- integer(n); ncl <- 0L } else {
    if (length(cluster) != n) stop("cluster must have one entry per row of x")
    cluster <- as.factor(cluster)
    if (any(tapply(grp, cluster, function(z) length(unique(z))) > 1L))
      stop("every cluster must lie inside one group")
    cl <- as.integer(cluster) - 1L; ncl <- nlevels(cluster) }
  z <- .C("C_ec_multi", as.integer(unlist(cols)), as.integer(n), as.integer(m), as.integer(g),
          grp, as.integer(k), as.integer(cl), as.integer(ncl),
          as.integer(c(pearson = 0L, dispersion = 1L)[[stat]]), as.integer(B), as.integer(seed),
          out = double(12), item_stat = double(m), item_p = double(m))
  o <- z$out   # obs_sum, obs_max, sum: p midp lo hi, max: p midp lo hi, largest null sum, largest null max
  part <- function(i, j, mx) list(statistic = o[i], p.value = o[j], midp = o[j + 1],
                                  mc.lower = o[j + 2], mc.upper = o[j + 3], max_null = mx)
  structure(list(sum = part(1, 3, o[11]), max = part(2, 7, o[12]),
                 items = data.frame(item = names(x), statistic = z$item_stat, p.value = z$item_p),
                 n = as.vector(table(group)), clusters = ncl, B = B, stat = stat,
                 exchangeable = if (ncl > 0) "clusters" else "respondents"),
            class = "ec_multi")
}

#' @export
print.ec_multi <- function(x, ...) {
  cat(sprintf("One exact test for %d items across %d groups, %s exchangeable, %s statistic\n",
              nrow(x$items), length(x$n), x$exchangeable, x$stat))
  fl <- 1 / (x$B + 1)
  for (nm in c("sum", "max")) { r <- x[[nm]]
    if (r$p.value <= fl + 1e-12)
      cat(sprintf("%-4s %.6g   p < %.3g (resolution floor); largest null draw %.6g\n",
                  nm, r$statistic, fl, r$max_null))
    else
      cat(sprintf("%-4s %.6g   p = %.4g   mid-p = %.4g   Monte Carlo 95%% interval [%.4g, %.4g]\n",
                  nm, r$statistic, r$p.value, r$midp, r$mc.lower, r$mc.upper)) }
  invisible(x)
}

#' A differential network without a Gaussian assumption
#'
#' Which pairs of variables are associated differently in two or more groups, tested without
#' assuming the data are Gaussian. The edge statistic is the heterogeneity of the Fisher-transformed
#' correlations across groups, \eqn{\sum_a (n_a - 3)(z_a - \bar z)^2}, the same quantity the Gaussian
#' theory uses; what changes is the reference. Instead of a chi-square that holds only under
#' normality, the statistic is referred to the law obtained by permuting the group labels, which
#' conditions on the observed data and the group sizes and is therefore exact whatever the marginals
#' are. With \code{rank = TRUE} the correlations are Spearman's, so the whole procedure is
#' distribution-free.
#'
#' @param x a matrix or data frame, variables in columns, respondents in rows.
#' @param group a factor or vector giving each respondent's group.
#' @param cluster an optional cluster per respondent, every cluster inside one group; when given the
#'   labels are permuted across clusters, each keeping its members.
#' @param rank replace each variable by its ranks, average ties, before correlating. The default,
#'   and what makes the test distribution-free; set \code{FALSE} to use the values as scores, which
#'   is what an ordinal coding already is.
#' @param B number of permutations.
#' @param seed integer seed; each permutation has its own stream.
#' @return A list with \code{sum} and \code{max} over the \eqn{p(p-1)/2} edges, each holding the
#'   observed statistic, the exact p-value, its mid-p companion, a Wilson interval and the largest
#'   null draw; and \code{edges}, a data frame of every edge with its statistic and the unadjusted
#'   per-edge p-value from the same permutations.
#' @details A network has \eqn{p(p-1)/2} edges, so testing them one at a time needs a correction
#'   that costs most of the power. The same permutations give ONE statistic for the whole network,
#'   the sum of the edge statistics when many edges differ a little and the maximum when one edge
#'   differs a lot, hence ONE p-value for "is the network different at all"; the per-edge p-values
#'   are then a closed follow-up, valid without further correction once the omnibus has rejected.
#'   Nothing is corrected because nothing is repeated.
#' @examples
#' set.seed(1)
#' n <- 120; p <- 5
#' g <- rep(c("a", "b"), each = n / 2)
#' z <- matrix(rnorm(n * p), n)
#' z[g == "b", 2] <- 0.8 * z[g == "b", 1] + 0.6 * z[g == "b", 2]   # one edge differs
#' ec_diffnet_test(z, g, B = 999, seed = 1)
#' @export
ec_diffnet_test <- function(x, group, cluster = NULL, rank = TRUE, B = 20000L, seed = 1L) {
  x <- as.matrix(x); n <- nrow(x); p <- ncol(x)
  if (p < 2) stop("x must hold at least two variables")
  if (length(group) != n) stop("group must have one entry per row of x")
  if (any(is.na(x))) stop("x must have no missing values")
  storage.mode(x) <- "double"
  if (isTRUE(rank)) x <- apply(x, 2, function(v) rank(v, ties.method = "average"))
  group <- as.factor(group); grp <- as.integer(group) - 1L; g <- nlevels(group)
  if (g < 2) stop("group must have at least two levels")
  if (min(table(group)) < 4) stop("every group needs at least four respondents")
  if (is.null(cluster)) { cl <- integer(n); ncl <- 0L } else {
    if (length(cluster) != n) stop("cluster must have one entry per row of x")
    cluster <- as.factor(cluster)
    if (any(tapply(grp, cluster, function(z) length(unique(z))) > 1L))
      stop("every cluster must lie inside one group")
    cl <- as.integer(cluster) - 1L; ncl <- nlevels(cluster) }
  m <- p * (p - 1L) / 2L
  z <- .C("C_ec_diffnet", as.double(x), as.integer(n), as.integer(p), as.integer(g), grp,
          as.integer(cl), as.integer(ncl), as.integer(B), as.integer(seed),
          out = double(12), edge_stat = double(m), edge_p = double(m))
  o <- z$out
  part <- function(i, j, mx) list(statistic = o[i], p.value = o[j], midp = o[j + 1],
                                  mc.lower = o[j + 2], mc.upper = o[j + 3], max_null = mx)
  nm <- colnames(x); if (is.null(nm)) nm <- paste0("V", seq_len(p))
  ij <- which(upper.tri(matrix(0, p, p)), arr.ind = TRUE)
  ij <- ij[order(ij[, "row"], ij[, "col"]), , drop = FALSE]
  structure(list(sum = part(1, 3, o[11]), max = part(2, 7, o[12]),
                 edges = data.frame(from = nm[ij[, "row"]], to = nm[ij[, "col"]],
                                    statistic = z$edge_stat, p.value = z$edge_p),
                 n = as.vector(table(group)), clusters = ncl, B = B, ranked = isTRUE(rank),
                 exchangeable = if (ncl > 0) "clusters" else "respondents"),
            class = "ec_diffnet")
}

#' @export
print.ec_diffnet <- function(x, ...) {
  cat(sprintf("Differential network over %d edges, %d groups, %s exchangeable%s\n",
              nrow(x$edges), length(x$n), x$exchangeable,
              if (x$ranked) ", ranks" else ", raw scores"))
  fl <- 1 / (x$B + 1)
  for (nm in c("sum", "max")) { r <- x[[nm]]
    if (r$p.value <= fl + 1e-12)
      cat(sprintf("%-4s %.6g   p < %.3g (resolution floor); largest null draw %.6g\n",
                  nm, r$statistic, fl, r$max_null))
    else
      cat(sprintf("%-4s %.6g   p = %.4g   mid-p = %.4g   Monte Carlo 95%% interval [%.4g, %.4g]\n",
                  nm, r$statistic, r$p.value, r$midp, r$mc.lower, r$mc.upper)) }
  e <- x$edges[order(x$edges$p.value), ][seq_len(min(5L, nrow(x$edges))), ]
  cat("strongest edges:\n")
  for (i in seq_len(nrow(e)))
    cat(sprintf("  %s -- %s   statistic %.4g   p = %.4g\n", e$from[i], e$to[i], e$statistic[i], e$p.value[i]))
  invisible(x)
}

#' Tail probabilities from an observed value and a null sample
#'
#' @param obs the observed statistic.
#' @param null a vector of draws from the null.
#' @return A list with the p-value, its mid-p companion, and a Wilson interval reflecting the Monte
#'   Carlo error in the number of draws supplied.
#' @examples ec_pvalue(2, rnorm(1000))
#' @export
ec_pvalue <- function(obs, null) {
  null <- as.double(null)
  z <- .C("C_ec_pvalue", as.double(obs), null, as.integer(length(null)), out = double(4))$out
  list(p.value = z[1], midp = z[2], mc.lower = z[3], mc.upper = z[4])
}

#' exactcond: exact conditional reference distributions
#'
#' Reference distributions exact at every sample size, obtained by conditioning on a statistic
#' sufficient for the nuisance parameter rather than by an asymptotic approximation. Two
#' conditionings are supplied: the Fisher-Yates law of a two-way table given both margins, and the
#' uniform law of the residual direction on the sphere of the residual space.
#'
#' @section What is here and why:
#' \describe{
#'   \item{\code{\link{ec_rc_test}}}{one table, both margins fixed. The dispersion statistic is for
#'     nominal rows against an ordinal column and is built from the row mean scores themselves, so
#'     where the column score is a signed difference the test is on the reported contrasts rather
#'     than on ranks of a variable that is almost all ties.}
#'   \item{\code{\link{ec_rc_two}}}{two tables, conditioning further on the total score, which the
#'     common association parameter multiplies and which the conditioning has fixed, so that
#'     parameter leaves the density altogether.}
#'   \item{\code{\link{ec_sphere_maxquad}}}{a statistic maximised over candidate structures, whose
#'     exact null follows from the residual direction being uniform on the sphere. This is what
#'     makes inference valid AFTER a spatial weight matrix has been chosen by fit, where the usual
#'     chi-squared reference rejects a true null far too often.}
#' }
#'
#' @section What is deliberately absent:
#' Nothing here is a bootstrap and nothing is asymptotic. A test that is exact by analytic means,
#' such as one whose null is obtained by inverting a characteristic function, needs no simulation
#' and has no place in this package.
#'
#' @name exactcond-package
#' @aliases exactcond
#' @keywords internal
"_PACKAGE"

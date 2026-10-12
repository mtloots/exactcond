## The claim these tests defend is EXACTNESS. A test that only checked the code ran would pass
## against a subtly biased sampler, so the central test compares the sampler against the
## conditional law obtained by COMPLETE ENUMERATION, atom by atom.

test_that("the sampler reproduces the enumerated conditional law", {
  rs <- c(3, 2, 2); cs <- c(3, 2, 2); N <- sum(rs); v <- c(0, 1, 2)
  w <- numeric(0); Tv <- numeric(0)
  for (a in 0:3) for (b in 0:2) { cc <- rs[1] - a - b
    if (cc < 0 || cc > 2) next
    for (d in 0:(3 - a)) for (e in 0:(2 - b)) { f <- rs[2] - d - e
      if (f < 0 || f > 2 - cc) next
      g <- cs[1] - a - d; h <- cs[2] - b - e; i <- cs[3] - cc - f
      if (g < 0 || h < 0 || i < 0 || g + h + i != rs[3]) next
      x <- matrix(c(a, b, cc, d, e, f, g, h, i), 3, byrow = TRUE)
      w <- c(w, exp(sum(lfactorial(rs)) + sum(lfactorial(cs)) -
                    lfactorial(N) - sum(lfactorial(x))))
      Tv <- c(Tv, sum(outer(v, v) * x)) } }
  w <- w / sum(w)
  dr <- ec_rc_null(rs, cs, stat = "linlin", B = 200000L, seed = 21)
  tv <- sort(unique(Tv))
  pe <- vapply(tv, function(t) sum(w[Tv == t]), 0)
  ps <- vapply(tv, function(t) mean(abs(dr - t) < 1e-9), 0)
  expect_lt(max(abs(pe - ps)), 0.005)
})

test_that("a 2x2 agrees with the closed-form exact test", {
  tt <- matrix(c(9, 3, 2, 10), 2, byrow = TRUE)
  fe <- stats::fisher.test(tt, alternative = "greater")$p.value
  ec <- ec_rc_test(tt, stat = "linlin", B = 200000L, seed = 11)
  expect_lt(abs(ec$p.value - fe), 0.005)
})

test_that("the p-value is invariant to affine rescoring", {
  tab <- matrix(c(3, 228, 72, 10, 548, 67, 7, 276, 131), 3, byrow = TRUE)
  a <- ec_rc_test(tab, v = c(-1, 0, 1), B = 20000L, seed = 41)
  b <- ec_rc_test(tab, v = c(4, 7, 10), B = 20000L, seed = 41)
  expect_equal(a$p.value, b$p.value)
  expect_equal(a$midp, b$midp)
})

test_that("the dispersion components are the row mean scores themselves", {
  tab <- matrix(c(3, 228, 72, 10, 548, 67, 7, 276, 131), 3, byrow = TRUE)
  r <- ec_rc_test(tab, v = c(-1, 0, 1), B = 1000L, seed = 1)
  expect_equal(r$row_means, as.vector(tab %*% c(-1, 0, 1)) / rowSums(tab))
  expect_equal(r$statistic, sum(r$n * (r$row_means - sum(tab %*% c(-1, 0, 1)) / sum(tab))^2))
})

test_that("the result does not depend on the thread count", {
  old <- Sys.getenv("OMP_NUM_THREADS", unset = NA)
  Sys.setenv(OMP_NUM_THREADS = "1"); a <- ec_rc_null(c(30, 40), c(35, 35), B = 2000L, seed = 61)
  Sys.setenv(OMP_NUM_THREADS = "4"); b <- ec_rc_null(c(30, 40), c(35, 35), B = 2000L, seed = 61)
  if (is.na(old)) Sys.unsetenv("OMP_NUM_THREADS") else Sys.setenv(OMP_NUM_THREADS = old)
  expect_identical(a, b)
})

test_that("the sphere draws have the law they must have", {
  n <- 7; A <- matrix(0, n, n); A[1, 1] <- 1
  d <- ec_sphere_maxquad(diag(n), A, B = 40000L, seed = 51)
  expect_gt(suppressWarnings(stats::ks.test(d, "pbeta", 0.5, (n - 1) / 2)$p.value), 0.001)
})

test_that("two identical tables give no evidence, even when strongly associated", {
  ## This is the case the obvious construction cannot do at all. Estimating each table's null
  ## distribution separately and convolving them needs the observed score to be reachable under NO
  ## association, and for this table it lies five standard deviations out, where two million draws
  ## never reach it. The failure would arrive here, where the correct answer is a LARGE p-value.
  a <- matrix(c(20, 8, 3, 9, 22, 7, 4, 8, 19), 3, byrow = TRUE)
  z <- ec_rc_two(a, a, B = 20000L, seed = 3)
  expect_gt(z$midp, 0.5)
  expect_equal(z$S, 2 * z$T1)
  expect_gt(z$accept, 0.01)
})

test_that("opposed association in the two tables is detected", {
  a <- matrix(c(20, 8, 3, 9, 22, 7, 4, 8, 19), 3, byrow = TRUE)
  z <- ec_rc_two(a, a[, 3:1], B = 20000L, seed = 3)
  expect_lt(z$midp, 0.001)
})

test_that("the two-table test holds its size with UNEQUAL margins", {
  ## The margins are deliberately opposed, which is where permuting a group label is invalid.
  skip_on_cran()
  u <- v <- c(0, 1, 2)
  make_P <- function(rm, cm, phi) { P <- exp(phi * outer(u, v))
    for (it in 1:200) { P <- P * (rm / rowSums(P)); P <- t(t(P) * (cm / colSums(P))) }
    P / sum(P) }
  drw <- function(n, P) { idx <- sample(9, n, TRUE, prob = as.vector(P))
    table(factor((idx - 1) %% 3, 0:2), factor((idx - 1) %/% 3, 0:2)) }
  set.seed(20260901)
  A <- make_P(c(.5, .3, .2), c(.45, .35, .2), 0.6)
  Bm <- make_P(c(.15, .35, .5), c(.2, .3, .5), 0.6)
  rej <- mean(vapply(seq_len(200), function(i)
    ec_rc_two(drw(150, A), drw(150, Bm), B = 3000L, burn = 4000L, thin = 10L,
              seed = 1000L + i)$midp < 0.05, TRUE))
  expect_lt(abs(rej - 0.05), 0.045)
})

test_that("misuse is refused rather than quietly accommodated", {
  expect_error(ec_rc_null(c(10, 10), c(10, 11)), "totals must agree")
  expect_error(ec_rc_test(matrix(c(1, -1, 2, 3), 2)), "non-negative")
  expect_error(ec_rc_test(matrix(1:9, 3), v = c(1, 2)), "one score per column")
  expect_error(ec_rc_two(matrix(1:9, 3), matrix(1:9, 3), v = c(0, 0.5, 1)), "integer scores")
})

## ---- the multi-item omnibus ------------------------------------------------------------------
test_that("with one binary item and respondents exchangeable, the omnibus is the exact test of the table", {
  set.seed(5); g <- rep(c("a", "b", "c"), c(40, 50, 45)); x <- rbinom(135, 1, ifelse(g == "c", 0.6, 0.85))
  om <- ec_multi_test(matrix(x, ncol = 1), g, stat = "pearson", B = 40000L, seed = 9)
  rc <- ec_rc_test(table(g, x), stat = "pearson", B = 40000L, seed = 9)
  expect_equal(om$sum$statistic, rc$statistic)
  expect_lt(abs(om$sum$p.value - rc$p.value), 0.01)
})

test_that("the sum and the maximum do not depend on the order of the items", {
  set.seed(6); g <- rep(1:3, c(30, 30, 30)); x <- matrix(rbinom(270, 1, 0.8), ncol = 3)
  a <- ec_multi_test(x, g, B = 3000L, seed = 2); b <- ec_multi_test(x[, 3:1], g, B = 3000L, seed = 2)
  expect_equal(a$sum$statistic, b$sum$statistic); expect_equal(a$max$statistic, b$max$statistic)
  expect_equal(a$sum$p.value, b$sum$p.value); expect_equal(a$max$p.value, b$max$p.value)
})

test_that("the cluster permutation keeps every cluster whole and rejects a cluster that straddles groups", {
  set.seed(7); g <- rep(c("a", "b"), each = 60); cl <- rep(1:12, each = 10)
  x <- matrix(rbinom(120 * 2, 1, 0.7), ncol = 2)
  r <- ec_multi_test(x, g, cluster = cl, B = 2000L, seed = 3)
  expect_equal(r$clusters, 12L); expect_equal(r$exchangeable, "clusters")
  expect_error(ec_multi_test(x, g, cluster = rep(1:12, 10), B = 200L), "inside one group")
})

test_that("a strong cluster effect inflates the respondent permutation and not the cluster permutation", {
  ## ten clusters per group, a cluster effect large against the item variance, NO group effect
  set.seed(8); reject <- c(resp = 0L, clus = 0L); R <- 60L
  for (i in seq_len(R)) {
    g <- rep(c("a", "b"), each = 100); cl <- rep(1:20, each = 10)
    u <- rnorm(20, 0, 1.5)[cl]
    x <- matrix(rbinom(200 * 3, 1, plogis(1 + u)), ncol = 3)
    pr <- ec_multi_test(x, g, B = 400L, seed = i)$max$p.value
    pc <- ec_multi_test(x, g, cluster = cl, B = 400L, seed = i)$max$p.value
    reject <- reject + c(pr < 0.05, pc < 0.05) }
  expect_gt(reject[["resp"]], reject[["clus"]])
  expect_lte(reject[["clus"]], 9L)                       # about 5 per cent of 60, with slack
})

test_that("the omnibus does not depend on the thread count", {
  old <- Sys.getenv("OMP_NUM_THREADS", unset = NA)
  set.seed(9); g <- rep(1:3, c(20, 20, 20)); x <- matrix(rbinom(180, 1, 0.75), ncol = 3)
  Sys.setenv(OMP_NUM_THREADS = "1"); a <- ec_multi_test(x, g, B = 1000L, seed = 13)
  Sys.setenv(OMP_NUM_THREADS = "4"); b <- ec_multi_test(x, g, B = 1000L, seed = 13)
  if (is.na(old)) Sys.unsetenv("OMP_NUM_THREADS") else Sys.setenv(OMP_NUM_THREADS = old)
  expect_equal(a$sum$p.value, b$sum$p.value); expect_equal(a$items$p.value, b$items$p.value)
})

## ---- the differential network ----------------------------------------------------------------
test_that("the differential network finds a planted edge and leaves the others alone", {
  set.seed(11); n <- 200; p <- 6; g <- rep(c("a", "b"), each = n / 2)
  z <- matrix(rnorm(n * p), n); z[g == "b", 2] <- 0.8 * z[g == "b", 1] + 0.6 * z[g == "b", 2]
  r <- ec_diffnet_test(z, g, B = 999L, seed = 1)
  expect_lt(r$max$p.value, 0.01)
  e <- r$edges[order(r$edges$p.value), ]
  expect_equal(paste(e$from[1], e$to[1]), "V1 V2")
  expect_gt(min(e$p.value[-1]), 0.05)
})

test_that("its p-value is uniform under the null, with Gaussian and with heavy-tailed data", {
  skip_on_cran()
  for (kind in c("gauss", "t3")) {
    set.seed(if (kind == "gauss") 21 else 22)
    pv <- vapply(seq_len(200), function(i) {
      y <- if (kind == "t3") matrix(rt(120 * 5, df = 3), 120) else matrix(rnorm(120 * 5), 120)
      ec_diffnet_test(y, rep(c("a", "b"), each = 60), B = 199L, seed = i)$max$p.value }, 0)
    expect_lt(abs(mean(pv < 0.10) - 0.10), 0.06)
  }
})

test_that("an affine rescoring of a variable leaves the network alone", {
  set.seed(12); y <- matrix(rnorm(150 * 4), 150); g <- rep(c("a", "b"), each = 75)
  a <- ec_diffnet_test(y, g, B = 499L, seed = 5)
  y2 <- y; y2[, 2] <- 3.5 * y2[, 2] - 7
  b <- ec_diffnet_test(y2, g, B = 499L, seed = 5)
  expect_equal(a$edges$statistic, b$edges$statistic)
  expect_equal(a$max$p.value, b$max$p.value)
})

test_that("the cluster permutation holds its size where the respondent permutation does not", {
  skip_on_cran()
  set.seed(13); rej <- c(resp = 0L, clus = 0L)
  for (i in seq_len(80)) {
    g <- rep(c("a", "b"), each = 90); cl <- paste(g, rep(rep(1:9, each = 10), 2))
    y <- matrix(rnorm(180 * 5), 180) + setNames(rnorm(18, 0, 0.9), unique(cl))[cl]
    rej <- rej + c(ec_diffnet_test(y, g, B = 199L, seed = i)$max$p.value < 0.05,
                   ec_diffnet_test(y, g, cluster = cl, B = 199L, seed = i)$max$p.value < 0.05) }
  expect_gt(rej[["resp"]], rej[["clus"]])
  expect_lte(rej[["clus"]], 12L)
})

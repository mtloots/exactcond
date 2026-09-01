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

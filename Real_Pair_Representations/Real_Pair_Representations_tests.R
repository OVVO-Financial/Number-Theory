# Reproducible base-R tests for the real-pair/quaternion implementation.
# Run in the directory containing both scripts:
# Rscript Real_Pair_Representations_tests.R
source("Real_Pair_Representations.R")

close_vec <- function(actual, expected, tol = 1e-10) {
  length(actual) == length(expected) &&
    all(is.finite(actual)) && all(is.finite(expected)) &&
    max(abs(actual - expected)) <= tol * (1 + max(abs(expected)))
}

native_qmul <- function(v, w) {
  a <- v[1]; b <- v[2]; cc <- v[3]; d <- v[4]
  e <- w[1]; f <- w[2]; g <- w[3]; h <- w[4]
  c(a*e - b*f - cc*g - d*h,
    a*f + b*e + cc*h - d*g,
    a*g - b*h + cc*e + d*f,
    a*h + b*g - cc*f + d*e)
}

one <- psi(c(1, 0, 0, 0))
i <- psi(c(0, 1, 0, 0))
j <- psi(c(0, 0, 1, 0))
k <- psi(c(0, 0, 0, 1))
zero <- psi(c(0, 0, 0, 0))

stopifnot(close_vec(phi_inv(phi(2 + 3i)), 2 + 3i))
stopifnot(close_vec(real_pair_multiply(phi(2 + 3i), phi(-4 + 2i)),
                     phi((2 + 3i) * (-4 + 2i))))
stopifnot(close_vec(real_pair_divide(phi(2 + 3i), phi(-4 + 2i)),
                     phi((2 + 3i) / (-4 + 2i))))
stopifnot(close_vec(real_pair_exp(phi(0.5 + 0.3i)),
                     phi(exp(0.5 + 0.3i))))
stopifnot(close_vec(real_pair_multiply(phi(3 - 2i),
                                      real_pair_inverse(phi(3 - 2i))),
                     phi(1 + 0i)))

stopifnot(close_vec(quat_multiply(i, j), k))
stopifnot(close_vec(quat_multiply(j, i), -k))
stopifnot(close_vec(quat_multiply(j, j), -one))
stopifnot(close_vec(quat_inverse(j), -j)) # Inverse of the quaternion basis element j
stopifnot(close_vec(quat_inverse(k), -k))
stopifnot(close_vec(quat_multiply(j, quat_inverse(j)), one))
stopifnot(close_vec(quat_divide(k, i), -j))
stopifnot(close_vec(quat_exp(zero), one))
stopifnot(close_vec(quat_exp(pi * j), -one))
stopifnot(close_vec(quat_exp(1e-10 * j),
                     psi(c(cos(1e-10), 0, sin(1e-10), 0))))

set.seed(20261104)
for (trial in seq_len(1000L)) {
  v <- rnorm(4)
  w <- rnorm(4)
  p <- psi(v)
  q <- psi(w)
  stopifnot(close_vec(psi_inv(p), v))
  stopifnot(close_vec(quat_multiply(p, q), psi(native_qmul(v, w))))
  stopifnot(close_vec(quat_multiply(p, quat_inverse(p)), one,
                       tol = 1e-9))
  stopifnot(close_vec(quat_multiply(quat_inverse(p), p), one,
                       tol = 1e-9))
  stopifnot(close_vec(quat_conjugate(p), psi(c(v[1], -v[2:4]))))
  stopifnot(close_vec(quat_norm_sq(p), sum(v * v)))
  stopifnot(close_vec(quat_norm_sq(quat_multiply(p, q)),
                       quat_norm_sq(p) * quat_norm_sq(q)))
}
stopifnot(inherits(try(quat_inverse(zero), silent = TRUE), "try-error"))
cat("All deterministic and 1000 randomized base-R tests passed.\n")

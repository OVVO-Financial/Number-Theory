# Real-pair representations of complex numbers and quaternions (R).
# Core arithmetic below uses numeric vectors; complex values are used only
# in phi() and phi_inv() for conversion and independent reference tests.

.check_vec <- function(v, n, label) {
  if (!is.numeric(v) || length(v) != n || any(!is.finite(v))) {
    stop(sprintf("%s must be a finite numeric vector of length %d", label, n))
  }
  invisible(TRUE)
}

phi <- function(z) {
  if (length(z) != 1L || !is.finite(Re(z)) || !is.finite(Im(z))) {
    stop("z must be one finite complex (or real) scalar")
  }
  c(Re(z) - Im(z), Re(z) + Im(z))
}

phi_inv <- function(p) {
  .check_vec(p, 2L, "p")
  complex(real = (p[1] + p[2]) / 2,
          imaginary = (p[2] - p[1]) / 2)
}

real_pair_add <- function(p, q) {
  .check_vec(p, 2L, "p"); .check_vec(q, 2L, "q")
  p + q
}

real_pair_subtract <- function(p, q) {
  .check_vec(p, 2L, "p"); .check_vec(q, 2L, "q")
  p - q
}

real_pair_multiply <- function(p, q) {
  .check_vec(p, 2L, "p"); .check_vec(q, 2L, "q")
  s <- (p[1] * q[2] + p[2] * q[1]) / 2
  t <- (p[2] * q[2] - p[1] * q[1]) / 2
  c(s - t, s + t)
}

real_pair_conjugate <- function(p) {
  .check_vec(p, 2L, "p")
  c(p[2], p[1])
}

real_pair_norm_sq <- function(p) {
  .check_vec(p, 2L, "p")
  sum(p * p) / 2
}

real_pair_inverse <- function(p) {
  .check_vec(p, 2L, "p")
  m <- max(abs(p))
  if (m == 0) stop("zero has no multiplicative inverse")
  u <- p / m
  ((2 / sum(u * u)) / m) * c(u[2], u[1])
}

real_pair_divide <- function(p, q) {
  real_pair_multiply(p, real_pair_inverse(q))
}

real_pair_exp <- function(p) {
  .check_vec(p, 2L, "p")
  a <- (p[1] + p[2]) / 2
  b <- (p[2] - p[1]) / 2
  c(exp(a) * (cos(b) - sin(b)),
    exp(a) * (cos(b) + sin(b)))
}

# Quaternion convention: q = z1 + z2*j,
# z1 = a + bi, z2 = c + di, and j*z = Conj(z)*j.
psi <- function(q) {
  .check_vec(q, 4L, "q")
  c(q[1] - q[2], q[1] + q[2],
    q[3] - q[4], q[3] + q[4])
}

psi_inv <- function(p) {
  .check_vec(p, 4L, "p")
  c((p[1] + p[2]) / 2, (p[2] - p[1]) / 2,
    (p[3] + p[4]) / 2, (p[4] - p[3]) / 2)
}

quat_add <- function(p, q) {
  .check_vec(p, 4L, "p"); .check_vec(q, 4L, "q")
  p + q
}

quat_subtract <- function(p, q) {
  .check_vec(p, 4L, "p"); .check_vec(q, 4L, "q")
  p - q
}

quat_multiply <- function(p, q) {
  .check_vec(p, 4L, "p"); .check_vec(q, 4L, "q")
  p1 <- p[1:2]; p2 <- p[3:4]
  q1 <- q[1:2]; q2 <- q[3:4]
  first <- real_pair_multiply(p1, q1) -
    real_pair_multiply(p2, real_pair_conjugate(q2))
  second <- real_pair_multiply(p1, q2) +
    real_pair_multiply(p2, real_pair_conjugate(q1))
  c(first, second)
}

quat_conjugate <- function(p) {
  .check_vec(p, 4L, "p")
  # Quaternion conjugation swaps the first pair and negates the second.
  c(p[2], p[1], -p[3], -p[4])
}

quat_norm_sq <- function(p) {
  .check_vec(p, 4L, "p")
  sum(p * p) / 2
}

quat_inverse <- function(p) {
  .check_vec(p, 4L, "p")
  m <- max(abs(p))
  if (m == 0) stop("zero quaternion has no multiplicative inverse")
  u <- p / m
  ((2 / sum(u * u)) / m) * quat_conjugate(u)
}

# Right division, p / q := p * q^(-1) (not q^(-1) * p).
quat_divide <- function(p, q) {
  quat_multiply(p, quat_inverse(q))
}

quat_exp <- function(p) {
  .check_vec(p, 4L, "p")
  q <- psi_inv(p)
  a <- q[1]
  v <- q[2:4]
  r <- sqrt(sum(v * v))
  # sinc(r) = sin(r) / r, including its continuous value at zero.
  if (r < 1e-4) {
    r2 <- r * r
    sinc <- 1 - r2 / 6 + (r2 * r2) / 120
  } else {
    sinc <- sin(r) / r
  }
  psi(c(exp(a) * cos(r), exp(a) * sinc * v))
}

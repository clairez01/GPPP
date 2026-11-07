# Univariate slice sampler (Neal 2003): stepping-out + shrinkage.
#
# x_cur:     current parameter value (on whatever scale log_post expects)
# lp_cur:    log_post(x_cur), passed to avoid a redundant evaluation
# log_post:  function(x) -> log unnormalised density
# w:         initial stepping-out interval width
# max_steps: maximum stepping-out expansions per side
#
# Returns list(x = new_value, lp = log_post(new_value)).
slice_sample_1d <- function(x_cur, lp_cur, log_post, w = 1.0, max_steps = 50L) {
  log_u <- lp_cur + log(runif(1L))

  L <- x_cur - w * runif(1L)
  R <- L + w

  J <- floor(max_steps * runif(1L))
  K <- (max_steps - 1L) - J
  while (J > 0L && log_post(L) > log_u) { L <- L - w; J <- J - 1L }
  while (K > 0L && log_post(R) > log_u) { R <- R + w; K <- K - 1L }

  repeat {
    x_new  <- runif(1L, L, R)
    lp_new <- log_post(x_new)
    if (lp_new > log_u) return(list(x = x_new, lp = lp_new))
    if (x_new < x_cur) L <- x_new else R <- x_new
  }
}

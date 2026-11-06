# =============================================================================
# sdy2354_table3.R -- reproduce Application II (Table 3) at the published design
#
#   MBMM_SDY2354_CSV=/path/to/sdy2354_igg_long.csv \
#     Rscript sdy2354_table3.R
#
# Obtain SDY2354 from ImmPort under its user agreement. The input must be a
# prepared CSV with columns subject, arm, target, day and value. Participant-
# level study data are intentionally not distributed with mbmm.
#
# Runtime is longer than the earlier in-sample approximation because the outer
# bank is refitted 56 times per question for exact leave-one-out scoring.
# `vignette("sdy2354")` runs the same construction on synthetic data at a
# smaller size.
#
# THE CONSTRUCTION IS TWO-LEVEL, and that is the whole point.
#
#   inner   For each shared point (m_s, g) the per-target nuisance coordinates
#           (A0, C1, r_pre, r_post) are INTEGRATED OUT over a pool of 1,200
#           draws -- pool_marginal().  They are not competing bank atoms.  If
#           they were, any candidate that happened to fit would represent its
#           shared point, a flat trajectory would score the same at g > 0 as at
#           g < 0, and the class coordinate would stop being identified.
#
#   offset  Each assay reads the shared gain with its own calibration offset,
#           g + w*delta, integrated over five delta values -- jitter_marginal().
#           A stalk assay reads w_stalk, the three full HAs read w_breadth.
#
#   outer   For each scored participant, the bank is fitted over the shared
#           points (g, m_s, w_stalk, w_breadth) using the other 55 participants
#           -- bank_loo() -- then restricted to {g > 0} and renormalised.
#
# The design totals 16 m_s x 36 g x 1,200 pool = 691,200 curves per target.
# That number is the inner integration grid, NOT the bank's support: the bank
# has 28 shared g x 16 m_s x 3 x 3 = 4,032 atoms.
#
# The output should be compared with Table 3 of the manuscript. Monte Carlo
# integration and platform-level optimiser tolerances can produce small numeric
# differences, but the construction here matches the stated leave-one-out
# analysis rather than an in-sample approximation.
# =============================================================================

library(mbmm)

input <- Sys.getenv("MBMM_SDY2354_CSV", unset = "")
if (!nzchar(input))
  stop("Set MBMM_SDY2354_CSV to a prepared SDY2354 CSV obtained from ImmPort.")
if (!file.exists(input)) stop("Input file does not exist: ", input)
sdy2354 <- utils::read.csv(input, stringsAsFactors = FALSE)
needed <- c("subject", "arm", "target", "day", "value")
if (!all(needed %in% names(sdy2354)))
  stop("Input must contain columns: ", paste(needed, collapse = ", "))
sdy2354 <- sdy2354[needed]
sdy2354$value <- as.numeric(sdy2354$value)
sdy2354$day <- as.integer(sdy2354$day)
if (any(!is.finite(sdy2354$value)) || any(sdy2354$value <= 0))
  stop("`value` must contain finite positive measurements.")

TARGETS <- c("SinaiIgG", "H18_IgG", "H2_IgG", "H9_IgG")
STALK   <- "SinaiIgG"            # chimeric construct, low-coating protocol
SIGMA   <- 0.17                  # profiled; see profile_sigma()
LAMBDA  <- 1
SEED    <- 20261021

## Fixed log10 design limits from the frozen manuscript analysis. They are
## constants here rather than functions of the complete participant table, so
## a held-out participant cannot alter the support used to score them.
BOUNDS <- list(
  SinaiIgG = list(A0 = c(0.537426497940624, 4.471111865180439),
                  C1 = c(0.037426497940624, 4.671111865180439)),
  H18_IgG  = list(A0 = c(2.522180839413665, 5.445030055327142),
                  C1 = c(2.022180839413665, 5.645030055327142)),
  H2_IgG   = list(A0 = c(2.517242084547646, 5.469858044800942),
                  C1 = c(2.017242084547646, 5.669858044800942)),
  H9_IgG   = list(A0 = c(2.744277120801843, 5.544356301181178),
                  C1 = c(2.244277120801843, 5.744356301181178)))

MS_N  <- 16L                     # m_s grid points, log-spaced in [0.1, 20]
MPOOL <- 1200L                   # nuisance pool per target
G     <- seq(-1.5, 2.0, by = 0.1)                 # outer grid
KEEP  <- which(G >= -1.1 - 1e-9 & G <= 1.6 + 1e-9)  # interior shared range
W     <- c(0, 0.25, 0.5)                          # jitter widths
DELTA <- c(-0.8, -0.4, 0, 0.4, 0.8)               # quintile midpoints of U(-1,1)
OFF   <- round(outer(W, DELTA) / 0.1)             # exact grid steps, by construction

MS <- 10^seq(log10(0.1), log10(20), length.out = MS_N)
# dense where the fast exponential dominates, coarse over the long decay
SCREEN <- (c(seq(1, 60, 1), seq(62, 150, 2), seq(154, 588, 8)) - 1) / 30

traj <- function(t, th)
  andraud_trajectory(t, th$A0, th$C1, th$C2, th$D1, th$D2, th$ms)
raw  <- function(t, th)
  andraud_trajectory(t, th$A0, th$C1, th$C2, th$D1, th$D2, th$ms, floor_at = -Inf)

lhs <- function(n, d, seed) {
  set.seed(seed)
  m <- matrix(NA_real_, n, d)
  for (j in seq_len(d)) m[, j] <- (sample.int(n, n) - 1L + runif(n)) / n
  m
}

if (!all(TARGETS %in% unique(sdy2354$target)))
  stop("Input is missing one or more required targets: ",
       paste(setdiff(TARGETS, unique(sdy2354$target)), collapse = ", "))
sdy2354$subject <- factor(sdy2354$subject)
sdy2354$arm <- factor(sdy2354$arm,
                      levels = c("Placebo", "LAIV8-IIV5",
                                 "LAIV8-IIV5_Adj", "IIV8_Adj"))
sdy2354$target <- factor(sdy2354$target, levels = TARGETS)
if (anyNA(sdy2354$arm)) stop("Input contains an unrecognised arm label.")
subjects <- levels(sdy2354$subject)
arm <- sdy2354$arm[match(subjects, sdy2354$subject)]
X   <- stats::model.matrix(~ arm)

# ---- one question ----------------------------------------------------------
# ref = "C1": g = log10(D1/C1), the boost gain over the post-prime plateau (Q1)
# ref = "A0": g = log10(D1/A0), the regimen gain over pre-vaccine baseline (Q2)
fit_question <- function(ref) {
  V <- list()
  for (k in seq_along(TARGETS)) {
    tg <- TARGETS[k]
    o  <- subset(sdy2354, target == tg)
    # seeded by the target's ALPHABETICAL position, as the original does via
    # levels(factor(target)) -- the pool draws depend on it
    u  <- lhs(MPOOL, 4L, SEED + match(tg, sort(TARGETS)))
    A0 <- 10^(BOUNDS[[tg]]$A0[1] + u[, 1] * diff(BOUNDS[[tg]]$A0))
    C1 <- 10^(BOUNDS[[tg]]$C1[1] + u[, 2] * diff(BOUNDS[[tg]]$C1))
    r_pre  <- -3 + u[, 3] * 6
    r_post <- -3 + u[, 4] * 6

    ob <- data.frame(id = factor(o$subject, levels = subjects),
                     time = (o$day - 1) / 30, value = o$value)
    ob <- droplevels(ob)

    # inner integral, one m_s slice at a time so memory stays flat
    E <- do.call(cbind, lapply(MS, function(ms) {
      cand <- do.call(rbind, lapply(G, function(g) {
        D1 <- (if (ref == "A0") A0 else C1) * 10^g
        data.frame(A0 = A0, C1 = C1, C2 = r_pre * C1,
                   D1 = D1, D2 = r_post * D1, ms = ms)
      }))
      LL <- trajectory_loglik(ob, cand, traj, sigma = SIGMA)
      # a candidate that leaves the positive orthant between observations is
      # not a trajectory; give it a residual no measurement can match
      LL[, !screen_positive(cand, raw, grid = SCREEN, floor = 0.01)$keep] <- -1e6
      pool_marginal(LL, rep(seq_along(G), each = MPOOL))
    }))

    # participants with no data for this target contribute nothing
    Ef <- matrix(0, length(subjects), ncol(E))
    Ef[match(levels(ob$id), subjects), ] <- E
    V[[tg]] <- lapply(seq_along(W), function(c)
      jitter_marginal(Ef, n_grid = length(G), keep = KEEP, offsets = OFF[c, ]))
    message(sprintf("  %-9s done", tg))
  }

  # outer: targets add, stalk reads w_stalk and the full HAs read w_breadth
  blocks <- list()
  for (c2 in seq_along(W)) for (c1 in seq_along(W)) {
    b <- V[[STALK]][[c1]]
    for (tg in setdiff(TARGETS, STALK)) b <- b + V[[tg]][[c2]]
    blocks[[length(blocks) + 1L]] <- b
  }
  LLb  <- do.call(cbind, blocks)
  gcol <- rep(rep(G[KEEP], times = MS_N), times = length(W)^2)

  loo <- bank_loo(LLb, lambda = LAMBDA)
  cls1 <- gcol > 0
  class_score <- function(region) vapply(seq_len(nrow(LLb)), function(i) {
    wr <- restrict_bank(loo$weights[i, ], region)$weights
    class_loglik(LLb[i, , drop = FALSE], wr)[1]
  }, numeric(1))
  l1 <- class_score(cls1)
  l0 <- class_score(!cls1)
  mass1 <- rowSums(loo$weights[, cls1, drop = FALSE])
  list(gate = fit_gate(l1, l0, X, firth = TRUE), l1 = l1, l0 = l0,
       atoms = ncol(LLb), mass = range(mass1))
}

message("Q1: boost gain, g = log10(D1/C1)")
q1 <- fit_question("C1")
message("Q2: regimen gain, g = log10(D1/A0)")
q2 <- fit_question("A0")

tab <- data.frame(
  arm = levels(arm),
  n   = as.vector(table(arm)),
  Q1  = round(100 * tapply(q1$gate$pi, arm, mean), 1),
  Q2  = round(100 * tapply(q2$gate$pi, arm, mean), 1))
print(tab, row.names = FALSE)

cat(sprintf(paste0("\nbank atoms: %d  |  leave-one-out class-1 mass ranges: ",
                   "Q1 %.3f-%.3f  Q2 %.3f-%.3f\n"),
            q1$atoms, q1$mass[1], q1$mass[2], q2$mass[1], q2$mass[2]))
for (nm in c("Q1", "Q2")) {
  q <- if (nm == "Q1") q1 else q2
  tt <- gate_test(q$l1, q$l0, X, firth = TRUE, B = 4000, seed = 1)
  cat(sprintf("%s  T_PO = %.2f on %d df   permutation p = %.5f\n",
              nm, tt$statistic, tt$df, tt$p_perm))
}

# Bayesian estimate of the appa genotype effect on pE11 and 4G8 immunoreactivity
# in the anterior rhombencephalon (AR) and hypothalamus superior lobe (HSL).
#
# Run from the repository root:  Rscript R/fit_genotype_effect.R
# Outputs go to results/.

library(rstan)
library(tidyverse)
library(ggdist)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# --- Build long-format model input from the raw per-region CSVs ---------------
# Each row of a raw CSV is one slide carrying one appa+/+ and one appa-/- section.
# slide_ids = "unique": slides are numbered consecutively across regions (pE11).
# slide_ids = "shared": slides are numbered 1..n within each region, so slide k
#   in AR and slide k in HSL share one slide effect (4G8, as in the published fit).
read_antibody <- function(antibody, slide_ids) {
  df <- map_dfr(c("AR", "HSL"), function(region) {
    raw <- read_csv(sprintf("data/%s_intensity_%s.csv", antibody, region),
                    show_col_types = FALSE)
    tibble(slide_in_region = seq_len(nrow(raw)),
           appa_wt = raw[[1]], appa_ko = raw[[2]], brain_region = region)
  })
  df$slide <- if (slide_ids == "unique") seq_len(nrow(df)) else df$slide_in_region

  df %>%
    pivot_longer(c(appa_wt, appa_ko), names_to = "genotype", values_to = "signal") %>%
    mutate(is_wt = as.integer(genotype == "appa_wt"),
           log_signal = log(signal),
           region_id = as.integer(factor(brain_region)),   # 1 = AR, 2 = HSL
           slide_id = as.integer(factor(slide))) %>%
    arrange(slide_id, region_id, is_wt)
}

model <- stan_model("stan/appa.stan")

fit_antibody <- function(antibody, slide_ids) {
  df <- read_antibody(antibody, slide_ids)
  stan_data <- list(N = nrow(df), R = max(df$region_id), J = max(df$slide_id),
                    y = df$log_signal, is_wt = df$is_wt,
                    slide = df$slide_id, region = df$region_id)

  fit <- sampling(model, data = stan_data, seed = 1234, chains = 4,
                  iter = 2000, warmup = 1000, control = list(adapt_delta = 0.95))
  saveRDS(fit, sprintf("results/fit_%s.rds", antibody))
  fit
}

# --- Posterior summary of the genotype effect (log WT/KO ratio) ---------------
summarise_delta <- function(fit, antibody) {
  as.data.frame(fit) %>%
    select(AR = `delta_region[1]`, HSL = `delta_region[2]`) %>%
    pivot_longer(everything(), names_to = "region", values_to = "delta") %>%
    group_by(region) %>%
    summarise(median = median(delta),
              cri90_low = quantile(delta, 0.05),
              cri90_high = quantile(delta, 0.95),
              p_wt_gt_ko = mean(delta > 0),
              fold_wt_over_ko = exp(median)) %>%
    mutate(antibody = antibody, .before = 1)
}

plot_delta <- function(fit, antibody) {
  draws <- as.data.frame(fit) %>%
    select(AR = `delta_region[1]`, HSL = `delta_region[2]`) %>%
    pivot_longer(everything(), names_to = "region", values_to = "delta")

  ggplot(draws, aes(x = region, y = delta)) +
    stat_halfeye(.width = 0.9, slab_alpha = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "gray40") +
    labs(x = "Brain region",
         y = sprintf("Genotype effect, log(WT/KO) %s immunoreactivity", antibody),
         caption = "Point: posterior median; bar: 90% credible interval") +
    theme_minimal(base_size = 14)
}

fits <- list(pE11 = fit_antibody("pE11", slide_ids = "unique"),
             `4G8` = fit_antibody("4G8", slide_ids = "shared"))

summary_tbl <- imap_dfr(fits, summarise_delta)
print(summary_tbl)
write_csv(summary_tbl, "results/genotype_effect_summary.csv")

iwalk(fits, function(fit, antibody) {
  ggsave(sprintf("results/genotype_effect_%s.pdf", antibody),
         plot_delta(fit, antibody), width = 5, height = 5)
})

# Convergence diagnostics
iwalk(fits, function(fit, antibody) {
  s <- summary(fit)$summary
  n_div <- sum(sapply(get_sampler_params(fit, inc_warmup = FALSE),
                      function(m) sum(m[, "divergent__"])))
  cat(sprintf("%s: max Rhat %.4f, min n_eff %.0f, divergent transitions %d\n",
              antibody, max(s[, "Rhat"]), min(s[, "n_eff"]), n_div))
})

writeLines(capture.output(sessionInfo()), "results/sessionInfo.txt")

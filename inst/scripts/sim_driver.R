# Runs replicates of simulation_study.R (full configuration), one .rds file
# per replicate, with a shared queue: each worker takes the next replicate
# that is neither done nor claimed by another worker (claims are directories,
# created atomically). Set before sourcing: STUDY ("A" or "B"), SCEN (scenario
# ids), REPS (replicate ids), BUDGET_MIN (stop starting jobs after it).
MODE <- "full"; RUN_A <- FALSE; RUN_B <- FALSE; OUT_DIR <- "sim_full"
source("simulation_study.R")
rep_dir <- file.path(OUT_DIR, "reps"); dir.create(rep_dir, showWarnings = FALSE, recursive = TRUE)
scen <- if (STUDY == "A") expand.grid(D = cfg$D, M = cfg$M, truth = cfg$truths, stringsAsFactors = FALSE) else cfg$scen_B
scen$id <- seq_len(nrow(scen))
grid <- expand.grid(s = SCEN, r = REPS)
grid <- grid[order(grid$r, grid$s), ]
t_start <- Sys.time()
for (k in seq_len(nrow(grid))) {
  s <- grid$s[k]; r <- grid$r[k]
  f <- file.path(rep_dir, sprintf("%s_s%02d_r%03d.rds", STUDY, s, r))
  if (file.exists(f)) next
  if (as.numeric(Sys.time() - t_start, units = "mins") > BUDGET_MIN) { cat("budget reached\n"); break }
  if (!dir.create(paste0(f, ".lock"), showWarnings = FALSE)) next
  job <- list(sc = scen[s, ], rep = r)
  t0 <- Sys.time()
  out <- if (STUDY == "A") run_rep_A(job) else run_rep_B(job)
  saveRDS(out, f)
  unlink(paste0(f, ".lock"), recursive = TRUE)
  err <- if (STUDY == "A") out$row$error else out$error
  cat(sprintf("%s %s scen %d (D=%d, M=%d, %s) rep %d: %.0f s%s\n", format(Sys.time(), "%H:%M"), STUDY, s,
              job$sc$D, job$sc$M, job$sc$truth, r, as.numeric(Sys.time() - t0, units = "secs"),
              if (!is.null(err) && !is.na(err)) paste(" ERROR", err) else ""))
}
cat("worker done\n")

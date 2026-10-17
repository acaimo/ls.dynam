# Scripts for the paper

These scripts reproduce the analyses of Caimo and Gollini, *Latent space
dynamic network actor models with position-rate link*, with the functions of
**ls.dynam**. Find them with `system.file("scripts", package = "ls.dynam")`
and copy them to a working directory.

| Script | Content |
|---|---|
| `radoslaw_application.R` | Application to the e-mail data with one latent dimension: fits, figures, evidence, sensitivity, held-out scores and calibration by simulation. Writes `radoslaw_results.rds`. |
| `radoslaw_dimension.R` | Two and three latent dimensions; reads `radoslaw_results.rds`. |
| `simulation_study.R` | Simulation study (`MODE <- "full"` for the settings of the paper, `"quick"` to check that it runs). |
| `sim_driver.R`, `sim_aggregate.R` | Run the replicates of `simulation_study.R` in several processes, one file per replicate, and summarise them. |

The application needs the data file `out.radoslaw_email_email` of the KONECT
collection (<https://konect.cc/networks/radoslaw_email>), saved as
`out_radoslaw.txt` in the working directory (or set `DATA_FILE` before
sourcing the script). Settings can be changed by assigning them before
`source()`, for example `QUICK <- TRUE` for short chains or
`N_STARTS <- 30` for the number of random starts used in the paper.

The headers of the scripts give more details on running times. For the
paper, the main fit of the application with 31 starting configurations took
about 40 minutes, the fit with three dimensions about 10 hours, and the
simulation study with 20 replicates per scenario about 8 hours of computing
(4 hours on two cores).

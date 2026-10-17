# ls.dynam 0.1.0

* First version, with the functions used for the paper "Latent space dynamic
  network actor models with position-rate link" (Caimo and Gollini).
* Compared with the analysis scripts of the paper: the default link of
  `lsdynam_fit()` and `lsdynam_rate_stage()` is the additive link in latent
  popularity and excess self-affinity, `c("pop", "rho")` (`"none"` when
  `D = 0`); components of the link can be given by name (for example
  `link_contrast(fit, "pop")`), and `component` is now the second argument
  of `link_contrast()`; new functions `read_konect_events()`,
  `prepare_events()`, `plot_activity_effects()`, `plot_traces()` and the
  methods `summary()` and `plot()` for fitted models.
* The scripts of the paper (application, higher dimensions, simulation
  study) are in `inst/scripts` and use the package.

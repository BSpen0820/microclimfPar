# ---------------------------------------------------------------------------
# Local smoke test for the calm-wind free-convection fix (fix/calm-wind-free-convection)
#
# Runs entirely on your machine using the package's bundled example data
# (climdata, vegp, soilc, dtmcaerth) -- no HPC, no real tile data needed.
# Takes well under a minute (the point model + a tiny grid run).
#
# What it checks: whether below-canopy air temperature (Tz) stays physically
# plausible relative to the input air temperature, both under the example
# data's normal wind and under an artificially-forced very-calm-wind version
# of the same data (the actual failure regime for this bug).
#
# NOTE: runpointmodel() has its own internal floor
# (weather$windspeed[weather$windspeed < 0.5] <- 0.5, in R/Cppwrappers.R)
# that this script does not bypass -- so the "calm wind" run below is
# effectively capped at 0.5 m/s, not truly zero. That's still a meaningful
# stress test (prior analysis found many real tile cells were already
# pathological well above 0.5 m/s), but it is not the most extreme case.
#
# Just source() this after devtools::load_all() or devtools::install() +
# library(microclimfPara) on the fix/calm-wind-free-convection branch.
# ---------------------------------------------------------------------------

# devtools::load_all()   # uncomment if not already loaded/installed

run_and_report <- function(climdata_in, label) {
  micropoint <- runpointmodel(climdata_in, 0.05, dtmcaerth, vegp, soilc)
  mout <- runmicro(micropoint, 0.05, vegp, soilc, dtmcaerth)
  Tz <- mout$Tz

  cat("\n===", label, "===\n")
  cat("Input air temp range:  ", paste(round(range(climdata_in$temp, na.rm = TRUE), 1), collapse = " to "), "C\n")
  cat("Input windspeed range: ", paste(round(range(climdata_in$windspeed, na.rm = TRUE), 3), collapse = " to "), "m/s\n")
  cat("Output Tz range:       ", paste(round(range(Tz, na.rm = TRUE), 1), collapse = " to "), "C\n")
  cat("Max Tz:                ", round(max(Tz, na.rm = TRUE), 1), "C\n")
  cat("# Tz values > 45C:     ", sum(Tz > 45, na.rm = TRUE), "\n")
  cat("# Tz values > (max input temp + 15C), i.e. implausible: ",
      sum(Tz > (max(climdata_in$temp, na.rm = TRUE) + 15), na.rm = TRUE), "\n")

  invisible(mout)
}

cat("Loading bundled example data (climdata, vegp, soilc, dtmcaerth)...\n")

# --- Baseline: package's normal bundled example weather ---
run_and_report(climdata, "BASELINE (bundled example weather)")

# --- Stress test: force very calm wind (still floored at 0.5 m/s inside
#     runpointmodel() -- see note above) while keeping everything else
#     (solar, humidity, pressure) the same, so heat flux is still real ---
climdata_calm <- climdata
climdata_calm$windspeed <- pmin(climdata_calm$windspeed, 0.05)
run_and_report(climdata_calm, "CALM-WIND STRESS TEST (windspeed capped at 0.05 m/s pre-floor)")

cat("\nDone. A healthy result: Tz stays within a few degrees of the input air\n")
cat("temp range in both runs, and the 'implausible' count is 0 in both.\n")
cat("A blown-up result (the old bug): Tz spikes far above input air temp,\n")
cat("especially in the calm-wind run.\n")

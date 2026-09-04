# Physics audit findings — microclimfPara C++ core

**Date:** 2026-09-04
**Branch this audit was run against:** `fix/calm-wind-free-convection`
**Context:** After fixing the original calm-wind Tz singularity (`rhcanopy()`/`TVbelow()`) and the analogous bug in the snow model's own ground temperature (`snowoneB()`), Bryan asked for a broader sanity check of the rest of the physics for similar overlooked bugs.

## Methodology

Two independent subagents, each briefed as a boundary-layer/land-surface-exchange physics expert with literature-search access, audited the same source (`microclimfCpp.cpp`, `microclimfParallel.cpp`, and the two headers):

- **Review 1** was told about the two already-fixed bugs (as context/pattern to recognize, not to re-flag) and asked to find more instances of the same class of problem plus any other numerical singularities.
- **Review 2** was a **blind** audit — given no leads at all about what had already been found, just the general task of auditing the codebase for physical implausibility and numerical robustness issues. It also diffed against upstream `ilyamaclean/microclimf` to separate fork-introduced bugs from inherited ones.

**Bryan's rule for this exercise:** an issue counts as confirmed/worth prioritizing only if *both* independent reviews flagged it on their own. Everything else is recorded here for later triage but hasn't been cross-validated.

A **third** subagent was then used to specifically verify a claim from Review 2 that the `snowoneB()` free-convection fix committed just before this audit (commit `2807859`) was itself flawed. That claim was confirmed, and the fix has already been corrected and committed (commit `301f656`) — see "Already resolved" below.

## Already resolved (this branch, before and during this audit)

1. **`rhcanopy()`/`TVbelow()` calm-wind Tz singularity** — commits `0d68988`, `0988f08`. Free-convection velocity scale (Beljaars 1995 / CLM5) + reinstated Fix-3 far-field clamp.
2. **`snowoneB()` ground/snow surface temp — same singularity, separate code path** — commit `2807859`. Same free-convection treatment applied to the snow model's own ground-to-reference-height conductance.
3. **`snowoneB()` free-convection gating was wrong** — commit `301f656`. The fix in #2 used absorbed radiation (`RabsG`) as its buoyancy proxy, which is always positive whenever there's downwelling longwave (i.e. essentially always, including clear calm winter nights), so the free-convection term was engaging unconditionally instead of only when the surface is actually unstable — contrary to the CLM5/Beljaars literature it was based on. Independent estimate: was warming modeled `Tg` by up to ~12°C on clear calm winter nights, suppressing real nocturnal snow-surface cooling. Corrected by converting absorbed radiation to an estimated *net* radiation (using `Ts≈Ta` as the zeroth-order surface-temperature guess for the emission term) before gating on it, plus a melt-season guard. Verified: daytime peak behavior essentially unchanged (13.7–13.9°C vs. 12.4°C previously); full nighttime-cooling recovery **not yet isolated in a dedicated test** — worth checking against real HPC data when this work is picked back up.

## Confirmed — flagged independently by BOTH reviews

These are the highest-confidence findings; both subagents found them with no coordination.

### A. Snow albedo integer-division staircase — `snowalbCpp()`

**Location:** `microclimfCpp.cpp` (`snowalbCpp`, the shared function) **and two inline copies** in `microclimfParallel.cpp` (`Snow2Worker` and `GridMicroSnow2Worker`) — all three need to change together or they'll drift apart again.

```cpp
IntegerVector hs(tsteps);   // hours since last snowfall, int
alb[i] = (-9.8740 * std::log(hs[i] / 24) + 78.3434) / 100.0;   // hs[i]/24 is INTEGER division
```

Because `hs[i]` and `24` are both `int`, albedo only takes the value corresponding to whole days — it's pinned at 0.95 for a full 24h after snowfall, then steps discontinuously to 0.783, then 0.715 at 48h, etc., instead of decaying smoothly. At ~800 W/m² clear-sky incident radiation that's up to a **~134 W/m² step change in absorbed shortwave** landing at whatever clock hour the prior snowfall happened to end — a real melt-timing distortion in a region with frequent small snowfall events resetting the clock repeatedly through the season. Also, for `hs=0..23` the expression is technically `log(0)` = `-inf`, only rescued by the downstream `if (alb > 0.95) alb = 0.95` clamp — works by accident, not by design.

**Fix:** `hs[i] / 24.0` in all three locations, with a small floor on the log argument (e.g. `max(hs/24.0, ~0.2)`) chosen so the unclamped value equals the fresh-snow ceiling, rather than relying on the post-hoc clamp.

### B. `canopycondCpp()` — NaN stomatal conductance at zero leaf area

**Location:** `microclimfCpp.cpp`, `canopycondCpp()`.

```cpp
P_sun = (1-exp(-k*PAI))/k;                       // fine
Rshade_abs = Rdif*((1-exp(-PAI))/PAI)*(1-om);    // PAI in the denominator
```

At `PAI = 0` (bare ground, rock, scree, burned terrain, or a seasonally-dormant meadow if the PAI raster goes to zero) this is `0/0 = NaN`. Consumed downstream (`TVaboveground()`) as `if (gS > 0.0) gV = ...` — `NaN > 0.0` is `false` in C++, so this silently falls through to `gV = 0`, meaning **all latent heat loss is dropped from the canopy energy balance** for that cell/hour with no warning; the failure looks like a normal (if slightly warmer) result rather than an obvious blow-up.

**Fix:** floor `PAI` the same way `zeroplanedisCpp()` and `canopysnowintCpp()` already do elsewhere in this codebase (e.g. `if (pai < 0.001) pai = 0.001`), and make the downstream check NaN-aware (`if (std::isfinite(gS) && gS > 0.0)`) so a future NaN doesn't get silently absorbed the same way.

## Flagged by Review 1 only (not independently re-checked)

Review 1 was primed with the calm-wind/forced-convection pattern, so take these with that in mind — they're plausible and Claude spot-checked the most severe one against the actual source, but they have not been through the same blind cross-check as the two items above.

1. **`windCpp()`'s `gHa` — forced-convection-only, `gmin=0.0001`, `psi_h` hardwired to 0** *(Severe if real)*. This is the conductance used for the ordinary, non-snow above/below-canopy grid output (`Tz`/`Tg`/`tleaf` — i.e. most of what actually gets used). Estimated to saturate at the `dTmx` clamp ceiling (not just spike occasionally) on any low-vegetation cell (bare rock, scree, alpine turf, burn scars) with real sun, producing a flat unrealistic plateau. Spot-checked: the code at that line does exist as described (`gturbCpp(out.uf, tiw.d, tiw.zm, zref, 43, 0, 0.0001)`). Suggested fix: same free-convection treatment as `rhcanopy()`/`snowoneB()`, using `zref` as the length scale.
   - **Caveat from Review 2's independent "checked and fine" pass:** Review 2 traced `gturbCpp`'s denominator and concluded there's no *mathematical* singularity (every caller clamps `psih` so the denominator can't hit zero) — but Review 2 was checking for a divide-by-zero, not evaluating whether the `gmin=0.0001` floor is too small to keep the *output temperature* physically realistic under calm wind + strong sun, which is Review 1's actual concern. These aren't necessarily in conflict, but they weren't independently validating the same claim — worth resolving with a dedicated test before fixing.
2. **Exact-zero wind speed → NaN in the snow model's `umu`** *(High if real)*. Three sites (`pointmodelsnow()`, `pointmprocess()`, `microclimatemodel_wrapper()`) compute `umu = uf/ufps` with no floor and no zero-guard; unlike snow state that's carried forward hour-to-hour, one calm hour with wind reported as exactly `0.0` (common — anemometers have a starting threshold) could poison the entire subsequent snowpack trajectory.
3. **`pointmprocess()` uses the unclamped `PenmanMonteithCpp` for `T0p`** *(Medium-High)*, whose result (`dtrp`) is claimed to drive a ground-heat-flux ratio (`dtR = dtr/dtrp`) that, if inflated, would suppress `G` and amplify #1.
4. **`BigLeafCpp`'s free-convection floor (`gfreeCpp(...)  * 2 * pai`) vanishes as `pai → 0`** *(Medium)* — same underlying issue as the already-fixed bugs, but scaled by leaf area so it disappears for sparse/absent vegetation, including for the *ground* surface which doesn't disappear when the leaves do.
5. **`snowalbCpp()` integer division** — see "Confirmed" section A above (this is where Review 1 found it too).

## Flagged by Review 2 only (blind review; not independently re-checked)

Review 2 also diffed against upstream `ilyamaclean/microclimf` and tagged findings `[FORK]` (introduced by this fork) vs `[UPSTREAM]` (inherited). Ranked by Review 2's own severity assessment:

1. **Snowmelt proportional to snowpack depth rather than available energy** *(`[UPSTREAM]`, potentially severe)* — `microclimfCpp.cpp`, `snowoneB()`'s melt calculation (`Fm = 583.3 * out.Tg * S`, where `S` is snow water equivalent). Algebraically this makes melt fraction per hour `≈ 0.0063 × T_s` regardless of available energy — estimated to imply "melt energy" of 2000–3500 W/m² for a 1.2–3m deep pack at a few degrees above freezing, several times the solar constant. If real, this would make deep packs ablate far too fast (days instead of weeks). **High priority to verify** — this is a pre-existing (not fork-introduced) issue, so may have been present in all prior runs.
2. **`snowdaymov()` returns 0°C for every positive `reqhgt`** *(`[FORK]`, potentially severe for sub-nivean work)* — the day-smoothing helper (`microclimfCpp.cpp`, new in this fork, replaces upstream's `snowdayan`) is claimed to use the raw `reqhgt` in a damping-depth formula that's only valid for negative (below-ground) query heights; for the standard positive `reqhgt` case it's claimed to silently produce an all-zero smoothed series via `maCpp`'s unguarded negative-window-width path, which would make the sub-snow temperature output a function of snow depth and mean annual temperature only, with no real weather signal. **Worth checking given this feeds sub-nivean thermal-refuge estimates.**
3. **`snowdaymov()` unvalidated day-index bounds** *(`[FORK]`)* — claimed possible out-of-bounds read/UB if `snowdays` is passed as a raw 0/1 indicator instead of calendar-day indices, or if `Tgref`'s span is shorter than expected (e.g. a subset run).
4. **`radoneB()`'s beam-recovery clip may be discarding real low-sun-angle radiation on slopes** *(`[FORK]`, relates to the area of the already-fixed SOLR/zenith-angle bug, but a different specific claim)* — Review 2 argues the `cosz < 0.065 → beam = 0` gate (from the earlier SOLR fix) zeroes real slope-normal direct radiation that a tilted surface can receive even at low sun elevation, and that `Rbeam > 1352 → 0` (vs. upstream's cap-not-zero) discards energy rather than conserving it. **This specifically wasn't independently checked by Review 1**, since Review 1 was told this area was already fixed and out of scope — so this is a fresh claim, not a re-confirmation, and deserves its own look before acting.
5. **`radoneB()` uses a hardcoded `0.5` where `cosz`/`si` belongs** *(`[UPSTREAM]`)* — claimed to under/overestimate under-canopy direct radiation on snow depending on sun angle, worst at low winter sun angles under forest canopy.
6. **Debug instrumentation (the `MCF_DEBUG_RHCANOPY`-gated `Rcout` traces, including ones Claude added) writes non-atomic `static` variables from RcppParallel worker threads** *(`[FORK]`)* — a real data race if the debug env var is ever set with `parallel=TRUE`; also, `Rcpp::Rcout` from a non-main thread is itself unsafe per RcppParallel's own documentation. **Should be dropped before any production run with the debug env var active, and definitely before merging** — same conclusion as the "drop debug commits before merge" note already in `calm-wind-tz-bug.md`.
7. **RH computed from unclamped `Tz`, then `Tz` is clamped afterward** (`TVaboveground()`, `snowabovepoint()`) — when the clamp actually bites (i.e. exactly when the calm-wind solver has gone furthest off the rails), the reported `(Tz, rh)` pair is claimed to be thermodynamically inconsistent with each other, which would corrupt any downstream VPD/wet-bulb calculation.
8. **A grab-bag of NaN-domain violations that propagate silently rather than erroring**, of varying realism/severity — see Review 2's "Finding 10" table for the full list (includes: two-stream radiation edge cases at extreme leaf reflectance/transmittance combos, a `dtR = dtr/dtrp` division that's claimed unreachable in practice, and a few others). Lowest priority of the list; recorded for completeness.

Review 2 also has a substantial "checked and concluded fine" section (serial/parallel consistency broadly clean except item 6/A above; `gturbCpp`'s own denominator; tall-canopy `zref` guarding; the already-fixed `rhcanopy()` change judged sound; stomatal/soil-water conductance bounding) — worth reading directly if resuming this work, since it saves re-deriving those checks.

## Next steps (Bryan's note, 2026-09-04)

Bryan wants to come back to this later and create a **new branch off `fix/calm-wind-free-convection`** (not off `main`) to work through these — since the free-convection fixes on that branch are a prerequisite/related context for a lot of this list, and it hasn't been merged yet.

Suggested order when this is picked back up, given what's confirmed vs. not:
1. Fix the two cross-confirmed items (A, B above) — low-risk, well-understood, one-line-ish changes each.
2. Get a third/blind opinion on the highest-stakes single-review claims before acting on them — especially the snowmelt-scales-with-depth claim (#1 under Review-2-only) and the `windCpp()` conductance floor claim (#1 under Review-1-only), since both could be shaping a lot of existing output if real.
3. Drop the debug instrumentation (`MCF_DEBUG_RHCANOPY`-gated traces) before any parallel production run or merge, regardless of the rest — that one's a legitimate correctness/safety issue on its own, not dependent on further triage.

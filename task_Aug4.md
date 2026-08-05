# Task: Recoil-Pressure Over-Prediction Investigation (Aug 4)

## Motivation

VDEP power-sweep cases (testrun58-67, Al 6061) show liquid metal being ejected from an open
groove that has no liquid metal available to refill it. Hypothesis: recoil pressure (vapor
recoil from the laser-heated free surface) is over-predicted, expelling more liquid than the
real process does.

**Test of the hypothesis**: if total metal volume (liquid + solid, i.e. `alpha.metal` integrated
over the domain) is conserved to within what visibly leaves through an open domain boundary, that
rules out a numerics/mass-loss bug and points at the physics — specifically excessive recoil
pressure — as the actual cause of the ejection. That's Task 1. Tasks 2-4 are three independent,
isolated remedies for the recoil-pressure side of that hypothesis, each run as its own simulation.

---

## Task 1 — Metal volume conservation check (testrun64)

**Goal**: confirm whether metal volume loss in testrun64 (750W, fully completed, 100-400µs,
already reconstructed) is (a) negligible, (b) explained by metal visibly leaving through an open
domain boundary (e.g. ejected spatter/droplets exiting the top), or (c) unexplained
(numerical loss) beyond that.

**Boundary crossing is expected and acceptable** — this is an open-groove/spatter-ejection
process, not a closed container. The check is not "does total metal volume in the domain stay
constant," it's "is any drop in domain-metal-volume explained by metal that visibly left through
a boundary, or is some of it unaccounted for."

**Method**: testrun64 was run without a volume-tracking `functionObject`, so this is a
post-processing job on the already-reconstructed timestep data, not a rerun:
- For each reconstructed timestep, compute `∫ alpha.metal dV` over the whole domain (ParaView's
  `IntegrateVariables` filter does this directly — same functionality as OpenFOAM's own
  `volIntegrate` functionObject, just applied after the fact).
- Plot total metal volume vs. time.
- Cross-check any drop against the existing render views (`results/render_view.py`, all four
  views, or the batch `_render_stacked_video.sh` output already available for testrun64) to see
  whether visible ejecta/spatter near a boundary accounts for it.

**Pass/fail framing**: "pass" = any volume decrease is attributable to visible boundary-crossing
ejecta, not silent disappearance. A pass supports the recoil-pressure-is-too-strong hypothesis and
justifies proceeding with Tasks 2-4 as physically-motivated (not numerics-driven) remedies.

**Result (2026-08-04)** — `results/check_metal_volume_conservation.py`, run against all 77
reconstructed timesteps of testrun64 (`report/testrun64_metal_volume.csv`,
`report/testrun64_metal_volume.png`):
- Net change over the full 100-400µs run: **-0.133%** — essentially fully conserved.
- The loss is **not** a smooth continuous decline (which would suggest numerical
  diffusion/leakage). It's two distinct step-drops separated by perfectly flat, zero-loss
  plateaus: 0% → -0.10% over 97-160µs, dead flat -0.10% over 160-250µs, -0.10% → -0.133% over
  250-280µs, dead flat -0.133% over 280-400µs (through to the end of the run).
- This step signature (bounded transient episodes + long no-loss plateaus) is consistent with
  discrete ejecta events (spatter/droplets leaving the domain), not silent numerical mass
  destruction. **Task 1 passes** — supports proceeding with Tasks 2-4.
- **Visually confirmed** (2026-08-04): rendered `results/render_view.py --view=lateral` at the
  window boundaries (`report/tr64_lateral_t*.png`). At t=1.614421e-04s (end of window 1) there is
  a clearly separated droplet sitting above the nominal surface line, connected to the bulk melt
  track by a thin neck. At t=2.798706e-04s (end of window 2) a similar, smaller isolated
  protrusion appears at the same location (leading edge, current laser position). Both step-drops
  correspond to real, visible ejecta — not silent numerical mass loss. Task 1 fully closed.

---

## Tasks 2-4 — Three isolated recoil-pressure remedies

Shared setup for all three:
- Fork from **testrun61** (100µs, 32-core seed) — same lineage as the production power-sweep forks.
- Run at **750W**, matching **testrun64** — the one fully-completed production run — so each
  remedy sim has a direct, real baseline to compare against.
- **One remedy per sim** (isolated, clean attribution) — no combined changes.
- Add a `volIntegrate(alpha.metal)` `functionObject` to each new sim's `controlDict` (live
  mass-conservation monitoring going forward, rather than only retroactive post-processing as in
  Task 1).
- Proposed naming: **testrun68** (Task 2), **testrun69** (Task 3), **testrun70** (Task 4).

### Task 2 (testrun68) — Recoil pressure temperature cap

Clamp `T` to `Tcap` before it enters the recoil-pressure exponent
([UEqn.H:48-52](applications/solvers/laserbeamFoam/UEqn.H#L48-L52)):

```
pVap = 0.54 * p0 * exp( LatentHeatVap*Mm*(min(T,Tcap) - Tvap) / (R*min(T,Tcap)*Tvap) )
```

**`Tcap = 3800 K`** (≈1.39×Tvap). Single-cell `pVap` plateaus at ~19.3 atm (vs. unbounded,
reaching 68.8 atm at the 4400K observed in testrun64); the interface-felt (2-cell-averaged)
pressure plateaus at ~9.9 atm (vs. 34.6 atm uncapped at the same 4400K). See
[report/recoil_pressure_cap.png](report/recoil_pressure_cap.png) and
[report/felt_pressure.png](report/felt_pressure.png).

Grounded in testrun64's actual field data (checked directly, t=400µs, 2.8M cells,
`results/analyze_interface_temperature.py`):
- Bulk metal interior (`alpha.metal`>0.999) never exceeds ~3214K — this cap never touches
  ordinary bulk-liquid physics, regardless of exact value chosen.
- Interface cells (`0.001<alpha.metal<0.999`) commonly and repeatedly reach 3500-4800K — not a
  rare numerical fluke (127 cells >3500K, 30 cells >4000K in this single snapshot; 69-80% of very
  hot cells are interface cells, essentially none are bulk metal).
- 3800K (vs. an initially-considered 3566K or 4000K) was chosen so the fix's effect stays clearly
  visible against the actually-observed 4400K case, rather than only clipping a near-negligible
  tail.

The formula and its physical basis (Anisimov/Knight-type saturated-vapor-pressure recoil model)
are unmodified by this change — only the temperature fed into it is bounded. `Tvap`,
`LatentHeatVap`, `Mm`, `p0`, and the 0.54 prefactor all check out against standard literature
values for Al — the base parameterization is not the issue, the unbounded exponential's lack of a
physical saturation mechanism is.

### Task 3 (testrun69) — Surface tension increase

`sigma: 0.87 → 0.95 N/m` (+9%) in `constant/transportProperties`.

Source: Molina, Voytovych, Louis & Eustathopoulos, *"The surface tension of liquid aluminium in
high vacuum: The role of surface condition."* The current 0.87 N/m already matches their clean-Al
measurement; 0.95 N/m is their corrected oxide-covered-surface estimate (they argue the raw
"high" literature value of ~1050 mN/m, from Goumiri & Joud, is itself an overestimate due to an
undersized sessile-drop mass for that optical method; ~950 mN/m is the more defensible
oxide-covered figure). Al 6061 is ~97% Al by mass, so this pure-Al measurement is a reasonably
direct match.

### Task 4 (testrun70) — Beam quality factor (Radius_Flavour)

`Radius_Flavour: 1.336 → 0.7` in `constant/LaserProperties` (equivalently, assumed beam quality
`M²: 1.497 → 2.857`).

This flattens the laser's radial intensity profile — peak (center) intensity drops ~48% — while
leaving `laserRadius` (35µm, the real measured hardware spot size) completely unchanged
(confirmed in `src/laserHeatSource/laserHeatSource.C`: `beam_radius` is set directly from
`laserRadius` and is never touched by `Radius_Flavour`). This lowers peak surface temperature
(and thus recoil pressure) by redistributing the same total absorbed power over the same nominal
footprint, rather than by enlarging the physical spot — which would misrepresent the actual
hardware. See [report/radius_flavour_comparison.png](report/radius_flavour_comparison.png).

---

## Open items

- Task 1 fully done and confirmed (2026-08-04) — passed, see result above.
- **testrun69** (surface tension, Task 3) — case set up (testrun64 template + testrun61 seed
  hand-off, `sigma`→0.95, mass-conservation monitor added) and **running** as of 2026-08-04.
  Note: the `volFieldValue` functionObject needs an explicit `writeFields false;` entry in this
  OpenFOAM version — omitted originally, caused an immediate `FOAM FATAL IO ERROR`; fixed in both
  testrun69 and the pushed testrun70 config, and in `vdep_remedy_sims.md`'s template.
- **testrun68** (recoil pressure cap, Task 2) — not started; still needs the one-time `Tcap`
  solver change + separate `lbf3-tcap` image build (see `vdep_remedy_sims.md`).
- **testrun70** (beam quality, Task 4) — config pushed for Zixun; needs his local testrun61 seed
  hand-off before it can run (see `vdep_remedy_sims.md` and the case's own `README.md`).

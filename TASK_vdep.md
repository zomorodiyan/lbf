# Task: Replicate VDEP Experimental Results (Al 6061, laserbeamFoam)

## Ultimate Goal

Reproduce experimentally measured vapor depression (keyhole) depth and morphology
for Al 6061 single-track laser scans at varying power levels, using the `laserbeamFoam`
VOF solver. Validation target: keyhole depth vs. laser power curve from X-ray
synchrotron or ex-situ characterization data.

---

## Physical Setup

- **Material:** Al 6061 (T_sol=855 K, T_liq=925 K, k=167 W/m·K, ρ=2700 kg/m³)
- **Laser:** 1064 nm, 35 µm radius, Drude/Fresnel absorptivity (~11% at normal incidence)
- **Scan speed:** 6000 mm/s (6 m/s)
- **Powers of interest:** 650, 700, 750, 800, 850 W
- **Domain:** 0.64 × 0.7 × 3.2 mm (x lateral, y depth, z scan)
  - Gas headspace: 0.3 mm above surface
  - Substrate: 0.4 mm below surface
  - Track: z = 0.5 mm → 2.9 mm (2.4 mm), 0.5 mm lead-in, 0.3 mm run-out
- **Total track duration:** 400 µs
- **Mesh:** 16×18×80 base cells (~39 µm), 3-level AMR → ~5 µm at interface
- **Decomposition:** 16 cores, hierarchical n=(2 1 8); 32 cores n=(2 2 8) for heavy cases

---

## Strategy

### Original power sweep (constant kappa) — testruns 43–48

Warmup in tr43 (1000W→500W, 0–100 µs), then forked into tr44–47 (650/700/750/800 W).
tr48 = 850W continuation of tr43 state, constant kappa 167 W/m·K.

**Observation:** tr44–47 produce a keyhole (VDEP confirmed) but melt pool is too short —
keyhole walls are in near-constant contact with the solidification front on both sides,
forming unphysical vapor tails. Keyhole depth also too shallow (<200 µm).

### Improved series (temperature-dependent kappa) — testruns 50–52

To elongate the melt pool and deepen the keyhole:
- **Thermal conductivity:** changed from constant 167 W/m·K to linear
  `poly_kappa (232.0 -0.076 0 0 0 0 0 0)` → 167 W/m·K at solidus (855 K),
  ~80 W/m·K at 2000 K, positive past vaporization temperature.
- **Solver fix required:** `UEqn.H` line 48 used raw `T` in the recoil pressure
  `exp()` — caused SIGFPE when gas-phase cells had near-zero T with non-constant kappa.
  Fixed by applying `Tsafe = max(T, 300 K)` clamp (matching the guard already in TEqn.H).
  Docker image rebuilt after fix.

#### tr50 — warmup with new kappa
- Copy of tr43, linear kappa, 16 cores
- Run 0–100 µs (same warmup: 1000W → 500W)
- Provides clean initial state for tr51/tr52

#### tr51 — low-power stabilization leg (650 W, 0–25% of track)
- Fork from tr50 at t=100 µs
- Power: 650 W constant
- Duration: 100–200 µs (first 25% of the scan track)
- Purpose: allow melt pool to reach quasi-steady state at moderate power
  before ramping to target power

#### tr52 — high-power production run (850 W, 25–100% of track)
- Fork from tr51 at t=200 µs
- Power: 850 W constant
- Duration: 200–400 µs (remaining 75% of scan track)
- Purpose: capture quasi-steady keyhole morphology at 850 W with elongated melt pool

---

## Progress

- [x] Material properties set for Al 6061 (tr37 onward)
- [x] Domain geometry finalized: 0.64×0.7×3.2 mm, 16×18×80 base cells (tr43)
- [x] Warmup strategy designed: 1000W→500W, avoids ejection instability
- [x] tr43 complete (0–100 µs, constant kappa)
- [x] tr44 (650W), tr45 (700W), tr46 (750W), tr47 (800W) run to ~337–400 µs
- [x] tr47 fully reconstructed (75 snapshots, VTK ready)
- [x] tr46 reconstructed to 337 µs (62 snapshots, VTK ready)
- [x] tr48 (850W, constant kappa) running — continuation of tr43 state
- [x] Identified melt pool too short issue in tr44–47 (constant kappa)
- [x] Fixed SIGFPE bug in UEqn.H (Tsafe clamp for pVap exp())
- [x] Docker image rebuilt with fix
- [ ] Run tr50 (warmup, linear kappa) — wait for successful completion to t=100 µs
- [ ] Fork tr50 → tr51 (650 W, 100–200 µs)
- [ ] Fork tr51 → tr52 (850 W, 200–400 µs)
- [ ] Reconstruct and post-process tr51/tr52
- [ ] Extract keyhole depth vs. power and compare to experiment

---

## Key Lessons Learned

- Al 6061 has higher thermal conductivity (167 W/m·K vs ~110 for AlSi10Mg) —
  keyhole does not form at target powers from a cold start; 1000W spike needed.
- Reducing gas headspace below 0.3 mm causes severe timestep collapse during
  1000W startup due to ejected droplets reaching the atmosphere boundary.
- The violent startup phase (first ~10 µs at 1000W) is the most expensive part;
  timestep recovers to ~5–10 ns once keyhole stabilizes around t=10–20 µs.
- tr43 stage 1 (0.7mm deep domain) completed in 18 min vs tr39–41 (0.5mm domain)
  still crawling after 4h — domain depth matters more than base cell count.
- Running 4 cases simultaneously causes ~12x slowdown vs solo due to memory
  bandwidth contention. Run sequentially or max 2 at a time.
- Non-constant poly_kappa crashes with SIGFPE in UEqn.H because pVap exp() uses
  raw T without a lower bound — fixed by adding Tsafe = max(T, 300K) clamp.
- processor dirs are root-owned — always delete via Docker, not sudo rm.

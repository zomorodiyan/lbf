# Cu / Low-Carbon Steel Dissimilar Weld — Progress Log

Goal: Numerically replicate the macrosegregation patterns (peninsula and beach features) from Soysal & Kou (2016) for Cu/low-carbon steel GTAW butt welding.

---

## 2026-06-18

**Case setup — `compressibleLaserbeamFoam/Cu_LCS_interface`**
Assembled the full case: 5×2×6 mm domain (Cu at x=0–4 mm, steel at x=4–5 mm, gas layers top and bottom), laser scanning along z at 1 mm/s with a 0.5 mm beam radius, 8-core decomposition. Resolved several configuration bugs including inside-out block mesh vertex ordering, wrong ray direction (`V_incident (0 1 0)` → `(0 -1 0)`), and zero laser deposition in parallel (fixed by switching to `radialPolarHeatSource yes`).

**200 W test run — stable but slow**
Ran a low-power geometry-check at 200 W nominal (~25 W absorbed, ~12.6% Cu absorptivity). Solver was stable with no crashes, but the acoustic CFL constraint in `compressibleLaserbeamFoam` limited the timestep to ~6×10⁻⁷ s — making even a 50 ms test run take over 6 hours. Identified this as a fundamental limitation of the compressible solver architecture that cannot be overcome through input-file parameters.

---

## 2026-06-19

**New solver — `multicomponentLaserbeamFoam`**
Created a new solver by removing the per-phase compressible density-update and acoustic terms from `pEqn.H`, replacing them with a standard incompressible Poisson solve. Also removed the kinetic energy work terms from `TEqn.H` and `createFields.H`. The solver retains the full multi-component VOF framework, phase-change physics (melting/solidification, recoil pressure), and laser ray tracing.

**Successful compilation**
The new solver compiled cleanly inside Docker with no errors. Binary confirmed at `/root/OpenFOAM/user-v2506/platforms/linux64GccDPInt32Opt/bin/multicomponentLaserbeamFoam`.

---

## 2026-06-27

**Case conversion — `multiComponentlaserbeamFoam/Cu_LCS_interface`**
Converted the `Cu_LCS_interface` case from the compressible 5-phase format (Cu + Cu_vapour + steel + steel_vapour + air) to the incompressible 3-phase format (Cu + steel + air), writing `phaseProperties`, `physicalProperties.Cu/steel/air`, and updating `setFieldsDict`, `fvSchemes`, and `fvSolution`. Mesh initializes correctly and `setFields` completes without errors.

**Docker image rebuild**
Triggered a full image rebuild (`CACHE_BUST=$(date +%s)`) to bake the new solver binary into the image. Pending completion before first test run.

---

**First successful run — acoustic CFL resolved**
After resolving missing input files (`p`, `turbulenceProperties`, `div(rhoPhi,T)` scheme), the solver ran stably with no errors. The timestep grew from 1.2×10⁻⁶ s to 1.5×10⁻⁵ s within the first 0.1 ms and is still climbing — 20× larger than the compressible solver's fixed 6.7×10⁻⁷ s limit. Energy deposition is ramping correctly with the laser power schedule.

---

## Pending

- Let timestep settle to confirm it reaches maxDeltaT=1e-4 s (target: ~150× speedup over compressible)
- Run parallel (8-core) version to verify decomposition works
- Calibrate power to match GTAW conduction-mode pool dimensions from Soysal & Kou (2016)
- Full 3-second production run to capture macrosegregation

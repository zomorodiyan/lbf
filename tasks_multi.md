# Task: Simulate Cu/Steel Dissimilar Weld Macrosegregation

## Physical Goal

Replicate the **peninsula and beach macrosegregation features** from:
> Soysal & Kou, *Macrosegregation in dissimilar-metal fusion welding*, Acta Materialia 2016

**Experiment:** GTAW butt weld, pure Cu (Grade 110) to low-carbon steel, 1.2 mm thick plates.
60 A × 9.5 V = 570 W, 1 mm/s travel speed, electrode offset **1.5 mm into Cu side**.

**Key physics to capture:**
- Peninsula: steel partially solid while Cu side is fully liquid (T_LW < T_L,steel = 1798 K)
- Beach: Cu-rich zone at fusion boundary pushed by Marangoni flow
- Asymmetric melt pool — Cu melts earlier (T_liq = 1358 K vs 1798 K for steel)
- Cu/Fe immiscible in liquid state → sharp interface, no Fickian mixing

**Domain (Case: `tutorials/multiComponentlaserbeamFoam/Cu_LCS_interface/`):**
- x: 0–4 mm = Cu, 4–5 mm = steel, interface at x=4 mm
- y: 0–0.4 mm gas, 0.4–1.6 mm metal (1.2 mm thick), 1.6–2 mm gas
- z: 0–6 mm scan direction, laser track z=1.5–4.5 mm at 1 mm/s
- Beam center: x=2.5 mm (1.5 mm offset into Cu), y=1.7 mm (100 µm above metal surface)

---

## Solver Goal

Create `multicomponentLaserbeamFoam` by stripping compressible-gas behaviour out of
`compressibleLaserbeamFoam` — eliminating the acoustic CFL constraint that limited
timestep to ~6×10⁻⁷ s — while keeping the N-phase VOF, laser ray tracing, Marangoni,
and phase-change physics.

**Why not use `laserbeamFoam`?**
`laserbeamFoam` is hardwired to exactly two phases via
`immiscibleIncompressibleTwoPhaseMixture`. Adding a third metal phase is
structurally impossible without replacing the entire mixture model.

**Why not keep `compressibleLaserbeamFoam` as-is?**
The compressible gas path (ideal-gas EOS for air/vapour phases) introduces
acoustic time-step constraints (CFL ~ cell/sound speed ≈ 9 ns for 40 µm cells
in steel), requires a fully coupled energy–density–pressure loop, and crashes
with "Negative initial temperature" when any vapour phase alpha → 0 in a cell
while `hConstThermo` tries to invert a near-zero internal energy field. None
of this is needed for single-track concentration-field experiments.

---

## What changes, and why

### A. New solver directory (copy, do not modify the original)

Copy `applications/solvers/compressibleLaserbeamFoam/` to
`applications/solvers/multicomponentLaserbeamFoam/`.

Rename the main `.C` file and update the `Make/files` EXE target to
`$(FOAM_USER_APPBIN)/multicomponentLaserbeamFoam`.

### B. `pEqn.H` — replace compressible pressure equation with incompressible Poisson

**Current (compressible):** builds one `p_rghEqnComps[i]` per phase that
contains `psi*ddt(p_rgh)` (compressibility), then assembles a combined
compressible+incompressible system. After solving it calls `correctRho()` to
update every phase density from the pressure change.

**Target (incompressible-like):** solve the pure Poisson equation that
`laserbeamFoam` uses:

```
∇·(rAUf ∇p_rgh) = ∇·phiHbyA − vDot
```

`vDot` (volumetric source from phase-change / vaporisation) is kept because
the vaporisation model in `multiphaseMixtureThermo::solve()` still computes it
and it is needed for mass conservation across the liquid–vapour interface.

Remove: `PtrList<fvScalarMatrix> p_rghEqnComps`, the compressibility loop,
`p_rgh_0`, `mixture.correctRho()`, `phase.dgdt()` update, and `K` update.

Keep: `vDot` subtracted as RHS source (same as before), `ghf*fvc::snGrad(rho)`
buoyancy, surface tension, `constrainPressure`, non-orthogonal corrector loop.

### C. `TEqn.H` — remove compressible work terms

**Current:** includes kinetic energy and pressure-work terms scaled by `rCv()`:

```cpp
+ (fvc::div(fvc::absolute(phi, U), p)
 + fvc::ddt(rho, K) + fvc::div(mixture.rhoPhi(), K)) * mixture.rCv()
```

These terms are present in a compressible energy equation (internal energy
form) and are negligible at melt-pool velocities (< 10 m/s, Ma << 0.01).
They add stiffness and couple back to the pressure field in a way that depends
on the compressible EOS.

Remove these three lines from the `fvScalarMatrix TEqn(...)` block.
The rest of the energy equation (laser deposition, latent heat TRHS,
`mass_dot`, thermal diffusion via `alphaEff`) is unchanged.

### D. Thermophysical properties for gas/vapour phases — switch to `rhoConst` EOS

**Current:** air uses `equationOfState perfectGas` → ψ = 1/(RT), density
depends on p and T, `correctRho()` is active.

**Target:** all gas-state phases (air, SS316L_vapour, Ti64_vapour, and any
future vapour) use `equationOfState rhoConst` → constant density, ψ = 0,
`correctRho()` is a no-op for these phases.

For `rhoConst` the `mixture` subdictionary must include an explicit `rho`
value. Use ambient-temperature densities:
- air: 1.2 kg/m³
- metal vapour phases: ~0.5 kg/m³ (approximate; does not affect flow inside
  the melt pool, only the inertia of the near-vacuum vapour plume)

This change also eliminates the thermo-inversion crash: with `rhoConst`, the
`hConstThermo` temperature inversion is only called for cells where that phase
is present at finite alpha, and ψ = 0 means no pressure-driven density update
can push the internal energy negative.

Files to edit (for the SS316L_Ti64 tutorial case):
- `constant/thermophysicalProperties.air`
- `constant/thermophysicalProperties.SS316L_vapour`
- `constant/thermophysicalProperties.Ti64_vapour`

And the same three files for `Cu_LCS_interface/constant/`.

### E. `createFields.H` — remove `K` field and `createK.H`

The kinetic energy field `K` is only used in the compressible work terms of
`TEqn.H`. Once those terms are removed, `K` is unused. Remove the
`#include "createK.H"` line and any explicit `K = 0.5*magSqr(U)` update in
`pEqn.H`.

---

## File-by-file change summary

| File | Action |
|------|--------|
| `applications/solvers/multicomponentLaserbeamFoam/` | New directory (copy of compressibleLaserbeamFoam) |
| `multicomponentLaserbeamFoam.C` | Rename, update app name string |
| `Make/files` | Update EXE name |
| `pEqn.H` | Replace compressible system with incompressible Poisson |
| `TEqn.H` | Remove compressible work terms (3 lines) |
| `createFields.H` | Remove `createK.H` include |
| `thermophysicalProperties.air` (both cases) | Switch to `rhoConst` |
| `thermophysicalProperties.*_vapour` (both cases) | Switch to `rhoConst` |

---

## What is NOT changed

- `multiphaseMixtureThermo` library — unchanged; the N-phase VOF, alpha
  equations, vaporisation model (`solve()`), Marangoni, surface tension,
  interface compression, and latent heat machinery all remain.
- `UEqn.H` — momentum equation is identical.
- `update.H` — phase property averaging loop unchanged.
- `compressibleLaserbeamFoam` — original solver untouched.
- Tutorial case geometry, mesh, boundary conditions, `transportProperties`,
  all metal-phase thermophysical properties — unchanged.

---

## Runtime speedup settings (current run: `cu_lcs_500w_amr1`)

### Gas layer thickness — 400 µm (confirmed, no change needed)
Top and bottom gas layers in `constant/blockMeshDict` are each 400 µm.
This is already at the minimum practical thickness to avoid boundary effects on the
Marangoni flow at the free surface. No change was made.

### `maxAlphaCo` raised from 0.2 → 0.5 (`system/controlDict`)
```
maxAlphaCo      0.5;  // raised from 0.2 — allows ~2.5x larger timestep
```
**Effect:** the VOF Courant number at the metal/gas interface is the binding timestep
constraint (not the bulk Courant number `maxCo`). Relaxing this from 0.2 to 0.5 allows
timesteps roughly 2.5× larger, cutting total wall time by ~60%.

**Risk:** coarser time resolution of the interface may allow Cu/air and steel/air
boundaries to smear slightly. At conduction-mode velocities (< 0.5 m/s) and 100 µm
AMR cells the smearing is expected to be sub-cell and acceptable for a calibration run.
**Revert to 0.2 for final production run** if interface sharpness matters.

### `nAlphaSubCycles` reduced from 4 → 2 (`system/fvSolution`)
```
nAlphaSubCycles 2;  // reduced from 4 — halves alpha solve cost per timestep
```
**Effect:** each PIMPLE outer iteration solves the alpha (VOF) transport equation twice
instead of four times. Alpha is the most expensive per-cell solve in a 3-phase system.
Halving sub-cycles cuts the alpha solve cost by ~50%, equivalent to a further 1.5–2×
speedup on top of the `maxAlphaCo` change.

**Risk:** fewer sub-cycles reduces the accuracy of the interface compression step.
Smearing and non-monotone alpha profiles are more likely at large `deltaT`. Monitor
alpha fields in the first snapshot (t=0.05s) for unphysical over/undershoot.
**Revert to 4 for final production run** if smearing is visible.

### AMR level: 1 (200 µm base → 100 µm at metal/gas interface)
2-level AMR (200 → 50 µm) was tested but requires ~48 hours for t=0.5s (too slow).
1-level AMR targets ~12-hour wall time at the current settings.

---

## Progress

- [x] A. Copy solver directory and rename
- [x] B. Rewrite `pEqn.H` — incompressible Poisson, no compressibility terms
- [x] C. Edit `TEqn.H` — remove compressible work terms
- [x] D. Edit `createFields.H` + `UEqn.H` — remove K field, add pRefCell
- [x] E. Switch all phases to `rhoConst` EOS (Cu, steel, air)
- [x] F. Build solver binary in Docker image (confirmed present)
- [x] G. Set up `Cu_LCS_interface` case for new solver (3-phase: Cu + steel + air)
- [x] H. Serial test run — stable, no errors, deltaT grows to maxDeltaT=1e-4 s (**180× faster than compressible**)
- [x] I. Parallel 8-core test run — stable, confirmed working with `--allow-run-as-root`
- [x] I2. Speedup settings applied: maxAlphaCo=0.5, nAlphaSubCycles=2, 1-level AMR (temporary, see above)
- [x] I3. Global Tvap clip added to TEqn.H — caps all cells at 2835 K (Cu vaporisation point); prevents unphysical Marangoni runaway; see solver changes doc for full rationale
- [ ] J. Calibrate laser power to match GTAW melt pool dimensions from Soysal & Kou
- [ ] K. Full 3-second production run — capture peninsula and beach features
- [ ] L. Reconstruct and visualise in ParaView, compare to weld.png cross-section

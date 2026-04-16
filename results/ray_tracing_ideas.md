# Ray Tracing Physics — Discussion Notes

*Context: laserbeamFoam, ray tracing in laser powder bed fusion / keyhole welding*

---

## The Core Problem

With unlimited bounces, total absorptivity → 1 regardless of single-bounce Fresnel absorptivity.
After N bounces with single-bounce absorptivity A: total absorbed = `1 - (1-A)^N`.
At N=10, A=0.3 → 97% absorbed. At A=0.1 → 65%.
Material sensitivity (resistivity) only shows up in the first few bounces and is then washed out.
This explains why changing resistivity in the testruns gave nearly identical high absorptivity trends.

---

## Option 1 — Larger `rayPowerRelTol` (absolute power cutoff)

**What it does:** Kill a ray when its remaining power is less than `rayPowerRelTol × maxRayPower`.
Current default is `1e-6` (one millionth of peak ray power), effectively zero.
Raising to e.g. `0.01` kills a ray once it has deposited 99% of its energy.

**Location:** `LaserProperties` dict → `rayPowerRelTol`

- Pro: simple, no code change
- Con: not physically motivated; remaining power silently discarded (energy non-conservative);
  does not fix the resistivity sensitivity problem for deep keyholes where power depletes fast

---

## Option 2 — `maxBounces` (recommended short-term fix)

**What it does:** Hard cap on the number of reflections per ray.
In real keyholes, rays bounce O(3–8) times before escaping or depositing nearly all energy.

**Implementation needed:**
- Read `maxBounces` from `LaserProperties` dict (default 5)
- Increment `bounceCount_` on each specular reflection (currently tracked but never incremented)
- Kill ray when `bounceCount_ >= maxBounces`
- Deposit remaining power in current cell (energy conservative) rather than discarding

**Location of change:** `src/laserHeatSource/laserHeatSource.C` in the interface-hit block (~line 1100)

- Pro: directly fixes the resistivity/material sensitivity problem; physically motivated;
  computationally safe (prevents infinite loops in trapped-ray geometry)
- Con: requires choosing a value; when killing the ray, remaining power must go somewhere

**Suggested workflow:** start with `maxBounces 5`, sweep from 3–10 to check sensitivity.

The two options can be combined: kill a ray when EITHER `bounceCount >= maxBounces`
OR `power < rayPowerAbsTol`.

---

## Option 3 — Energy-conservative ray death

**What it does:** Whenever a ray is killed (by any mechanism — power, bounces, bounding box exit),
deposit the remaining power in the last occupied cell rather than discarding it.
This ensures `TotalQ` equals laser input power exactly.

Currently rays that exit the bounding box lose their remaining energy silently.

- Pro: physically correct; easy to verify (check TotalQ vs. input power)
- Con: may concentrate energy in cells at domain boundary if many rays exit there

---

## Option 4 — Diffuse/Lambertian scattering component

**What it does:** Real keyhole walls are not perfect mirrors — surface roughness adds a diffuse
(randomly scattered into hemisphere) component to each reflection.
Mix specular and diffuse with a roughness parameter `s` (0 = mirror, 1 = fully diffuse).

- Pro: physically richer; naturally reduces effective bounce depth; breaks the geometric trap
  where a ray bounces perfectly between parallel walls
- Con: more complex implementation; requires defining a roughness parameter per material

---

## Option 5 — Grazing-angle Fresnel behaviour

At near-grazing incidence the Fresnel equations give R → 1 (total reflection, zero absorption).
A ray entering a deep narrow keyhole at shallow angle may bounce many times with near-zero
absorption per bounce, making it effectively transparent to the material.
Worth checking whether keyhole geometry creates this regime.

No code change needed; diagnostic: log `theta_in` and `absorptivity` per bounce at debug level.

---

## Option 6 — Vapor plume / plasma absorption in keyhole cavity

**Physics:** Vaporized metal in the keyhole forms a plasma that absorbs laser radiation
via inverse Bremsstrahlung (free-electron absorption). At higher power densities this
can absorb 5–50% of the laser power before it reaches the keyhole walls.

**Current state:** Cells with `alpha < dep_cutoff` (gas/vapor) are fully transparent.
The `e_num_density` parameter is used only for the Fresnel equations on the metal surface,
not for plasma absorption in the cavity.

**To implement:**
- Identify cavity cells: `alpha < dep_cutoff` AND high temperature (or vapor indicator field)
- Apply Beer-Lambert attenuation along the ray path through cavity cells:
  `power *= exp(-mu_plasma * path_length)` where `mu_plasma` depends on plasma electron density
  and laser wavelength
- Plasma electron density in the cavity is much lower than bulk metal and is temperature-dependent

- Pro: physically complete for high-power keyhole regimes; makes the model qualitatively
  different and more realistic
- Con: requires plasma density / temperature model; substantially more implementation effort;
  may not be significant at lower power densities

**Recommendation:** assess whether current power density puts you in the plasma-dominated
regime before investing in this. Plasma absorption becomes significant when keyhole is deep
and laser power density is high (> ~10^6 W/cm²).

---

## Priority Summary

| Option | Effort | Impact | Recommendation |
|--------|--------|--------|----------------|
| `maxBounces` | Low | High — directly fixes resistivity sensitivity | Do first |
| Energy-conservative death | Low | Medium — physical correctness | Do alongside maxBounces |
| Larger `rayPowerRelTol` | None | Low — complementary cutoff | Tune after maxBounces |
| Diffuse scattering | Medium | Medium — better wall physics | After validation |
| Vapor plume absorption | High | High if high-power regime | Longer term |
| Grazing angle diagnostic | None | Diagnostic only | Run now to check |

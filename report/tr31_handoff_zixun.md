# testrun31_vdep_2_Al — Setup Handoff

## Context

testrun30_vdep_1_Al was the first AlSi10Mg VDEP simulation (1000 W nominal, 1500 mm/s,
35 µm beam radius). The keyhole drilled through the full 400 µm metal plate, indicating
the domain was too shallow and the effective absorbed power too high. testrun31 corrects
both issues.

---

## Changes from testrun30

### 1. Absorbed power: 700 W → 750 W

Trapp et al. (2017) measured in-situ absorptivity of aluminium powder during LPBF using
direct calorimetry at 1070 nm and reported:

- Conduction mode: ~15%
- Keyhole mode: ~75%

Reference: J. Trapp, A.M. Rubenchik, G. Guss, M.J. Matthews,
"In situ absorptivity measurements of metallic powders during laser powder-bed fusion
additive manufacturing", *Applied Materials Today*, 9, 341–349 (2017).
https://doi.org/10.1016/j.apmt.2017.08.006

Based on this, the correct input power for a 1000 W nominal laser in keyhole mode Al is
**750 W** (75%). This is higher than testrun30's 700 W — the wrong direction for
reducing keyhole depth — but it is the physically supported value and is used here.

### 2. Tvap: 2743 K → 2900 K

To compensate for the higher absorbed power AND to calibrate keyhole depth to the VDEP
experimental target of ~350 µm, Tvap is increased aggressively from the pure-Al boiling
point (2743 K at 1 atm) to **2900 K**.

Rationale:
- Pure Al Tvap is a lower bound; AlSi10Mg alloy effective boiling onset is higher
  (Raoult's law gives +22 K minimum; dynamic keyhole conditions push it further)
- Tvap appears in the Clausius-Clapeyron exponent — increasing it by 157 K reduces
  the evaporative recoil flux by ~40%, which is sufficient to compensate for the extra
  power and bring the keyhole depth down toward 350 µm
- This is a calibration target: if testrun31 depth is still > 350 µm, increase to
  3000 K in testrun32; if depth < 350 µm, step back toward 2850 K

In `constant/transportProperties`:
```
Tvap    2900;
```

### 3. Domain depth: 500 µm → 900 µm total

testrun30 had 100 µm gas + 400 µm metal. The keyhole hit the bottom wall, which
distorts late-time flow and pressure. testrun31 uses:

| Layer | testrun30 | testrun31 |
|-------|-----------|-----------|
| Gas above surface | 100 µm | 300 µm |
| Metal (solid) | 400 µm | 600 µm |
| Total y domain | 500 µm | 900 µm |

A 350 µm target keyhole depth sits comfortably within 600 µm of metal with 250 µm
of margin at the bottom.

### 4. Mesh: update blockMeshDict and setFieldsDict

Keep ~40 µm base cells in y → ceil(900/40) = **23 cells in y**.

In `system/blockMeshDict`: change y block to 0.9 mm with 23 cells.
```
vertices
(
    (-0.16  0     0)   // x: ±160 µm, y: 0→0.9 mm, z: unchanged
    ( 0.16  0     0)
    ( 0.16  0.9   0)
    (-0.16  0.9   0)
    (-0.16  0     2.5)
    ( 0.16  0     2.5)
    ( 0.16  0.9   2.5)
    (-0.16  0.9   2.5)
);
blocks
(
    hex (0 1 2 3 4 5 6 7) (8 23 64) simpleGrading (1 1 1)
);
```

In `system/setFieldsDict`: metal starts at y = 0.3 mm (300 µm gas layer):
```
boxToCell { box (-1 0.0003 -1) (1 0.0009 1); }  // metal y: 300–900 µm
```

### 5. Everything else: unchanged from testrun30

- Scan speed: 1500 mm/s
- Beam radius: 35 µm (Gaussian, Radius_Flavour 1.336)
- All AlSi10Mg material properties (rho, nu, cp, kappa, sigma, LatentHeat, etc.)
- Mesh decomposition: 32 cores, hierarchical (2 2 8)
- AMR: 2 levels, 40 µm → 10 µm at interface
- endTime: 0.0015 s

---

## Run command

Same as testrun30 — see TESTRUNS.md for the 32-core Docker launch command,
substituting `testrun31_vdep_2_Al` for the case name.

---

## Success criterion

Keyhole depth in the quasi-steady region (laser position z ≈ 1.0–1.5 mm) should
be **~350 µm** to match the VDEP experiment. Measure from the undisturbed metal
surface (y = 300 µm) to the deepest alpha.metal = 0.5 isosurface.

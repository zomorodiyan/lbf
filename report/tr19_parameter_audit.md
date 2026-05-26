# testrun19_316L — Parameter Audit
*Written 2026-05-21. Cross-reference: testrun30_vdep_1_Al was created from testrun19 as the AlSi10Mg baseline.*

---

## 1. Values in testrun19 that diverge from 316L literature

### 1a. `elec_resistivity = 1.6384e-7 Ω·m`  *(tuned, not physical)*
**File:** `constant/transportProperties`  
Measured 316L liquid resistivity is **~1.2–1.4×10⁻⁶ Ω·m** — roughly 8–9× higher.  
The comment says "A0=10% at normal incidence (Drude-Fresnel, lambda=1.064um)": this value was
back-calculated to reproduce a target absorptivity, not taken from bulk measurements.  
**Risk:** If you change laser wavelength or want to study absorptivity scaling, this is the
parameter to revisit, not a physical material property.

### 1b. `beta = 5.0e-6 1/K`  *(3× too small)*
**File:** `constant/transportProperties`  
Liquid 316L volumetric thermal expansion is typically **~1.5×10⁻⁵ 1/K**.  
The low value suppresses buoyancy-driven (Rayleigh-Bénard) flow relative to physical 316L.
For LPBF melt pools the Marangoni number dominates over the Grashof number by a large margin,
so the effect on melt pool shape is expected to be small — but not zero.

### 1c. `poly_kappa = (11.5 0.00307 ...)` — comment is wrong
**File:** `constant/transportProperties`  
Polynomial evaluates to:  
- at 300 K: 11.5 + 0.00307×300 = **12.4 W/(m·K)**  
- at 1800 K: 11.5 + 0.00307×1800 = **17.0 W/(m·K)**  

The inline comment says "~16.4 W/(m·K) at 300K, ~25 W/(m·K) at 1800K" — those numbers are
**wrong for these coefficients**. They describe a steeper polynomial that was not committed.  
Literature 316L: ~14–16 W/(m·K) at 300 K, rising to ~25 W/(m·K) near 1800 K.  
The polynomial underestimates conductivity by ~15–30%; heat accumulation in the melt pool
is consequently slightly over-predicted.

### 1d. `Marangoni_Constant = -2.28e-4 N/(m·K)` — correct for the model, weak vs. clean 316L
**File:** `constant/transportProperties`  
This is the *exact value from Morohoshi (2013) for pure iron, η=0* — internally consistent
with the surface tension study (see SURFACE_TENSION_STUDY.md).  
However, most literature measurements for clean 316L give dσ/dT ≈ **−4.3×10⁻⁴ to −5×10⁻⁴ N/(m·K)**.
Using the Morohoshi η=0 value thus underestimates outward Marangoni flow by ~2×.
This is a deliberate modeling choice tied to the Morohoshi series, not an error.

### 1e. `poly_cp = (500 0.0 ...)` — constant, ignores high-T rise
**File:** `constant/transportProperties`  
316L cp is ~500 J/(kg·K) at 300 K but rises to ~600–650 near 1400 K.
The constant value underestimates energy storage at high temperatures; melting depth is
likely slightly over-predicted as a result.

---

## 2. Additional flags found during testrun30_vdep_1_Al setup

### 2a. `initial/T` wall gradient — scales with conductivity
**File:** `initial/T` (all wall patches except atmosphere)  
```
fixedGradient -50.0  [K/m]
```
Heat flux out of the wall = kappa × 50 W/m².  
- 316L (kappa ≈ 12.4 W/m·K): ~620 W/m² loss  
- AlSi10Mg (kappa ≈ 165 W/m·K): **~8,250 W/m² loss** — 13× higher  

The gradient was tuned (or left as-is) for 316L. For AlSi10Mg it represents a much
stronger cooling boundary condition at frontAndBack, rightWall, and leftWall. Given the
small wall areas in LPBF (320 µm wide × 500 µm deep × 2.5 mm long), the total power
removed is still small relative to the 147 W laser, but if thermal gradients near the edges
matter for your result, this should be re-examined.

### 2b. `maxDeltaT = 2e-6 s` — may need reduction for AlSi10Mg
**File:** `system/controlDict`  
AlSi10Mg thermal diffusivity: α = kappa/(ρ·cp) ≈ 165/(2670×900) ≈ **6.9×10⁻⁵ m²/s**  
316L: α ≈ 12.4/(8000×500) ≈ **3.1×10⁻⁶ m²/s**  
→ AlSi10Mg diffuses heat ~22× faster per unit length.  

Fourier stability estimate over a 10 µm cell: dt_F = dx²/(4α) ≈ (10⁻⁵)²/(4×6.9×10⁻⁵) ≈ **3.6×10⁻⁷ s**.  
The current maxDeltaT (2×10⁻⁶ s) is about 5.5× above this. The temperature solver is
implicit so this is not a crash risk, but accuracy of within-timestep heat diffusion is
reduced. In practice, the free-surface CFL constraint (`maxCo=0.1`) usually controls
to dt ~ 10⁻⁸–10⁻⁷ s, keeping this academic — but worth checking in the early timesteps log.

### 2c. `timeVsLaserPower = 147.3684 W` — inherited from 316L, not reviewed
**File:** `constant/timeVsLaserPower`  
This power level was chosen for 316L keyhole/melt-pool conditions. Typical commercial
AlSi10Mg LPBF uses P=200–370 W, v=1000–1800 mm/s (much faster scan). The current
case uses v=200 mm/s (slow) and P=147 W — intentionally low-energy / slow-scan for
numerical stability, or carry-over from 316L? Confirm this is the intended operating point.

### 2d. Naming conflict with SURFACE_TENSION_STUDY.md
**File:** `SURFACE_TENSION_STUDY.md`  
That document lists `testrun30` as a planned 316L case (Morohoshi σ(T) for η=0.01).
The new aluminum case is `testrun30_vdep_1_Al` — a different directory, no file conflict.
However, the planned 316L testrun30 will need a distinct suffix (e.g. `testrun30_316L_sigma`)
to avoid confusion. Update SURFACE_TENSION_STUDY.md when that case is created.

### 2e. `laserRadius = 3.5e-5 m` (35 µm) — not changed, not material-specific but study-relevant
**File:** `constant/LaserProperties`  
The combination P=147 W / v=200 mm/s / r=35 µm gives a very different energy density than
commercial AlSi10Mg recipes. For AlSi10Mg the melt pool depth will be very different from
316L due to the lower melting point (880 K vs 1723 K) and much higher conductivity.
The 500 µm domain depth (y) copied from testrun19 is likely more than sufficient but the
melt pool may be relatively shallow — the domain depth is not a concern.

---

## 3. AlSi10Mg values still uncertain in testrun30

| Parameter | Value used | Uncertainty |
|---|---|---|
| `elec_resistivity` | 3.0×10⁻⁷ Ω·m | Liquid AlSi10Mg not well tabulated; affects Drude absorptivity |
| `Marangoni_Constant` | −1.5×10⁻⁴ N/(m·K) | Strongly alloy/impurity dependent; range −1e-4 to −3e-4 |
| `emS` | 0.09 | Polished Al; oxidised LPBF powder could be 0.10–0.30 |
| `poly_kappa` | (199 −0.112 …) | Linear fit; no data below solidus/above liquidus phase change |
| `beta` | 1.1×10⁻⁴ 1/K | Literature range for liquid Al alloys 0.9–1.2×10⁻⁴ |

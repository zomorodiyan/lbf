# Surface Tension Study

**Lineage:** tr19 → tr19_hires (1 µm write) → **tr26** (2-level refinement, 2.5 µm, restart t=1.111 ms) → tr27–30 (0.1 ms runs from tr26 snapshot)

## Cases (all restart from testrun26_316L)

| Case | Change vs. tr26 | Recompile |
|---|---|---|
| testrun27 | Full Morohoshi σ(T) for pure 316L (η=0); CSF sigma +25% → 2.3310 N/m | yes |
| testrun28 | Full Morohoshi σ(T) for 316L+PLC η=0.01 (carbon/oxygen); T*=1963 K, competing vortices | yes |
| testrun29 | `elec_resistivity` = 7.04×10⁻⁸ → A₀ = 5% (halved from 10%) | no |
| testrun30 | Full Morohoshi σ(T) poly for η=0.01; dσ/dT changes sign at T*=1963 K | **yes** |

## Material reference (Morohoshi 2013, `report/surface_tension_report.tex`)

| η | σ at 1873K (N/m) | dσ/dT (N/m/K) | Marangoni |
|---|---|---|---|
| 0 (pure) | 1.8648 | −2.28×10⁻⁴ | outward |
| 0.01 | 1.7675 | +3.41×10⁻⁴ | inward |
| 0.1  | 1.4209 | +1.244×10⁻³ | inward |

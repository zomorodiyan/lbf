# Multi-Component Laser Simulations — Objectives

## Purpose

This directory contains multi-material simulations using `compressibleLaserbeamFoam`.
The overarching goal is to use LaserbeamFoam as a tool to study **macrosegregation phenomena
in dissimilar-metal laser processing** — specifically to reproduce and explain experimentally
observed features that arise when two immiscible or thermophysically mismatched metals are
melted together.

---

## Target Article

**Macrosegregation in dissimilar-metal fusion welding**
T. Soysal, S. Kou
*Acta Materialia*, 2016
DOI: 10.1016/j.actamat.2016.(...) *(confirm from paper)*

### Experimental materials (from Section 2)

- **Cu:** commercially pure Grade 110 (Cu ≥ 99%, O = 0.04%, trace Ag)
- **Steel:** low-carbon steel
- Both: 1.2 mm thick, 60 mm wide, 150 mm long
- Process: gas-tungsten arc welding (GTAW), butt joint, no filler metal

### Key findings from the article

The authors butt-welded pure Cu to low-carbon steel and identified two macrosegregation
mechanisms governed by the relationship between the bulk weld metal liquidus temperature
T_LW and the individual material liquidus temperatures:

| Condition | Feature | Description |
|---|---|---|
| T_LW < T_L1 (steel) | **Peninsula** | An unmixed island of steel survives inside the weld pool — cooler liquid weld metal quickly freezes the intruding liquid steel before it can dissolve |
| T_LW > T_L2 (Cu) | **Beach** | An irregular zone of unmixed Cu forms at the fusion boundary — cooler liquid Cu layer quickly freezes as the liquid weld metal intrudes |

### What the weld cross-section shows (image `weld.png`)

The cross-sectional micrograph shows two panels (a) and (b) at 500 µm scale:

**(a) Top panel — macrosegregation features:**
- **Irregular beach of unmixed Cu** (orange, left): large irregular Cu-rich zone along the
  fusion boundary. Cooler liquid Cu layer freezes as the bulk weld metal intrudes.
- **Partially mixed Cu** zone transitioning into bulk weld metal (grey centre)
- **Peninsula of unmixed steel** (dark, right): a solid steel protrusion extending into
  the weld pool. Cooler liquid weld metal quickly freezes intruding liquid steel.
- **Flow during welding** arrows visible — convective flow patterns indicated in bulk

**(b) Bottom panel — layered structure:**
- **Layered weld metal**: bulk weld pool undercools into the metastable miscibility gap.
  Cu-rich liquid (matrix) plus Fe-rich liquid (spheres) form alternating layers, possibly
  due to rotational flow in the pool.
- **Regions containing spheres**: Fe-rich spheres, lighter than surrounding Cu-rich liquid,
  float to tops of layers — buoyancy-driven stratification.

### Primary simulation target

**Replicate the peninsula feature** — this is the most tractable for the solver because it
is driven by the liquidus temperature contrast and melt pool fluid dynamics, without
requiring CALPHAD. The peninsula forms on the steel side because T_LW < T_L1(steel):
liquid weld metal (bulk Cu-Fe mix, lower liquidus) solidifies against the steel face before
the steel can fully dissolve into the pool.

---

## Simulation Plan

### Case 1 — SS316L / Ti-6Al-4V (baseline, currently running)

**Path:** `SS316L_Ti64_interface/`

A first 3D bare-plate dissimilar-material run to establish the workflow, verify solver
behaviour at a material interface, and observe melt pool asymmetry driven by the ~230 K
solidus temperature contrast between 316L (1648 K) and Ti64 (1878 K).

**What we expect to see:** asymmetric melt pool (316L side melts deeper and wider),
Marangoni-driven cross-interface flow, different solidification front shapes on each side.

**Limitation:** 316L/Ti64 form brittle FeTi/Fe₂Ti intermetallics in reality — the solver
treats the interface as linearly mixed, which is physically approximate.

---

### Case 2 — Pure Cu / Low-Carbon Steel (peninsula & beach replication)

**Path:** `Cu_LCS_interface/`

**Goal:** replicate the peninsula and beach macrosegregation features reported by Soysal & Kou.

#### Experimental process parameters (Section 2 of article)

| Parameter | Value |
|---|---|
| Process | GTAW butt weld, no filler metal |
| Polarity | DC electrode negative (straight polarity) |
| Current | 60 A |
| Voltage | 9.5 V |
| Heat input | **570 W** (60 A × 9.5 V) |
| Travel speed | **1 mm/s** |
| Electrode offset | **1.5 mm into Cu side** from joint interface |
| Shielding gas | Ar, 11 L/min |
| Back-purge | Ar purge box |
| Welding mode | Manual, forward angle 60–70°, no torch oscillation |
| Plate dimensions | 1.2 mm thick × 60 mm wide × 150 mm long |

#### Key geometric insight: joint orientation vs. scan direction

In the experiment the arc **travels along the joint length** (z-direction in our domain).
The Cu/steel interface is a **longitudinal plane** — parallel to the x-z plane, separating
Cu (x < interface) from steel (x > interface). The arc is offset 1.5 mm into the Cu side
in the **x-direction** (across the joint width), not along the scan.

In our domain:
- x = width (5 mm total): Cu at x=0–4.0 mm, steel at x=4.0–5.0 mm, interface at **x=4.0 mm**
- Laser y position: **y=1.7 mm** (inside top gas layer, 100 µm above metal surface at y=1.6 mm)
- y = build/depth (2 mm total): gas below (0–0.4 mm), metal (0.4–1.6 mm), gas above (1.6–2 mm)
- z = scan direction (5 mm): laser travels along the joint, **parallel to the interface**

The 1.5 mm electrode offset into Cu is a physical compensation for Cu's high thermal
conductivity (~400 W/m·K vs ~50 W/m·K for steel): centering the arc on the interface would
funnel almost all heat into Cu, leaving the steel side unmelted. Offsetting into Cu balances
the melt pool across both materials. Beam center is at x=2.5 mm (1.5 mm from the interface at x=4.0 mm).

#### Domain schematics

**Top view (x-z plane, looking down in -y direction)**

```
                             ^ z=4.5mm (beam end)
                             |
                             | scan direction (+z)
                             |
z=6mm  +---------------------*--------+---------------------+
       |                     |        |                     |
       |                     |        |                     |
       |                     * z=3.0mm|                     |
       |                     |        |                     |
       |   Cu (Grade 110)    |        |  low-carbon steel   |
       |                     |        |                     |
       |                     * z=1.5mm|                     |
       |                     | (beam  |                     |
       |                     |  start)|                     |
z=0mm  +---------------------+--------+---------------------+
      x=0                 x=2.5mm  x=4.0mm              x=5.0mm
                          (beam    (interface)
                           track)
```

**Cross-section view (x-y plane, looking along z)**

```
y=2.0mm  +--------------------------------------------------------------+
         |         gas / atmosphere (laser entry from topWall)          |
y=1.6mm  +--------------------------------------------------------------+
         |                              |                                |
         |   Cu (Grade 110)             | low-carbon steel               |
         |   TSol=1356K  TLiq=1358K     | TSol=1768K  TLiq=1798K        |
         |                              |                                |
         |      o  <-- beam center      |                                |
         |         x=2.5mm             |                                |
         |                     interface at x=4.0mm                     |
y=0.4mm  +--------------------------------------------------------------+
         |              gas / atmosphere (bottomWall open)               |
y=0      +--------------------------------------------------------------+
        x=0           x=2.5mm        x=4.0mm                        x=5.0mm
        (Cu far wall) (beam center)  (interface)                 (steel far wall)
```


**Material contrast:**

| Property | Pure Cu (Grade 110) | Low-Carbon Steel |
|---|---|---|
| TSolidus | ~1356 K | ~1768 K |
| TLiquidus | ~1358 K | ~1798 K |
| ρ (liquid) | ~8000 kg/m³ | ~7000 kg/m³ |
| μ (liquid) | ~4×10⁻³ Pa·s | ~6×10⁻³ Pa·s |
| dσ/dT | ~−0.13×10⁻³ N/(m·K) | ~−0.43×10⁻³ N/(m·K) |
| Miscibility | Immiscible in liquid state (Cu-Fe phase diagram) | — |

**Key modelling choices:**
- `interfaceCompression (Cu steel) 1` — immiscible pair, maintain sharp liquid-liquid interface
- `interfaceDiffusion (Cu steel) 0` — no Fickian mixing
- `sigma (Cu steel) 0` — negligible liquid-liquid interfacial energy between immiscible metals
- Large liquidus contrast (~410 K) is the primary driver of peninsula formation

**What the solver can capture:**
- Peninsula: steel remaining partially solid while Cu side is fully liquid ✓
- Beach: Cu-rich zones pushed to fusion boundary by Marangoni-driven flow ✓
- Asymmetric melt pool geometry — deeper on Cu side (lower melting point) ✓
- Convective flow patterns consistent with the arrows shown in weld.png ✓

**What the solver cannot capture:**
- Layered structure from metastable liquid-liquid miscibility gap — requires CALPHAD ✗
- Fe-rich spheres floating due to buoyancy within Cu matrix — requires resolved
  liquid-liquid phase separation, not VOF mixing ✗

---

## Notes on Solver Limitations

Material properties at the dissimilar interface are computed by **linear volume-fraction
weighting** (rule of mixtures) of per-phase values — confirmed in `update.H` lines 68–70:

```cpp
TSolidus  += alpha1 * phaseTSolidus;
TLiquidus += alpha1 * phaseTLiquidus;
LatentHeat += alpha1 * phaseLH;
```

This is physically approximate for immiscible systems. The correct solidus/liquidus for a
Cu-Fe mixture is strongly nonlinear (eutectic, miscibility gap) and would require a
CALPHAD-coupled thermodynamic database (e.g. OpenCalphad with a validated Fe-Cu TDB file)
to compute accurately. This remains an open research problem for dissimilar-material LPBF.

# Paste into ParaView -> View -> Python Shell
# ─────────────────────────────────────────────────────────────────
# BEFORE running:
#   1. Load rays_laser0.vtk.series
#   2. Apply Filters -> Connectivity
#   3. Select the Connectivity filter in the Pipeline Browser
#   4. Paste this whole script into the Python Shell
# ─────────────────────────────────────────────────────────────────

from paraview.simple import *
import random

# ── Parameters ───────────────────────────────────────────────────
MODE        = 'random'    # 'random' or 'specific'
N_SELECT    = 10          # rays to show (MODE='random')
SEED        = 42          # change for a different draw
SPECIFIC    = [0, 45, 90, 135, 180, 225, 270, 315]
TOTAL_RAYS  = 360
TUBE_RADIUS = 1e-5        # adjust to your domain size
# ─────────────────────────────────────────────────────────────────

if MODE == 'random':
    random.seed(SEED)
    ids = sorted(random.sample(range(TOTAL_RAYS), N_SELECT))
else:
    ids = sorted(SPECIFIC)

print(f"Showing RegionIds: {ids}")

connectivity = GetActiveSource()

# Auto-detect whether RegionId is stored on points or cells
info = connectivity.GetDataInformation()
pt_arrays = [info.GetPointDataInformation().GetArrayInformation(i).GetName()
             for i in range(info.GetPointDataInformation().GetNumberOfArrays())]
cl_arrays = [info.GetCellDataInformation().GetArrayInformation(i).GetName()
             for i in range(info.GetCellDataInformation().GetNumberOfArrays())]

print(f"Point arrays: {pt_arrays}")
print(f"Cell  arrays: {cl_arrays}")

if 'RegionId' in pt_arrays:
    assoc = 'POINTS'
elif 'RegionId' in cl_arrays:
    assoc = 'CELLS'
else:
    raise RuntimeError("RegionId not found — is the Connectivity filter selected?")

print(f"RegionId is on {assoc}")

# Also detect where 'power' lives for coloring
color_assoc = 'POINTS' if 'power' in pt_arrays else ('CELLS' if 'power' in cl_arrays else None)

# One Threshold per selected ray
thresholds = []
for ray_id in ids:
    t = Threshold(Input=connectivity)
    t.Scalars        = [assoc, 'RegionId']
    t.LowerThreshold = ray_id
    t.UpperThreshold = ray_id
    t.ThresholdMethod = 'Between'
    thresholds.append(t)

combined = AppendDatasets(Input=thresholds) if len(thresholds) > 1 else thresholds[0]

# Threshold -> UnstructuredGrid; convert back to PolyData for Tube
surface = ExtractSurface(Input=combined)

tube = Tube(Input=surface)
tube.Radius       = TUBE_RADIUS
tube.NumberofSides = 8

display = Show(tube)

if color_assoc:
    ColorBy(display, (color_assoc, 'power'))
    GetColorTransferFunction('power')
    RescaleTransferFunctionToDataRange()
else:
    print("'power' field not found — color manually in ParaView")

Hide(connectivity)
RenderAllViews()

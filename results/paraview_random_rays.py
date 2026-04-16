# ParaView Programmable Filter — Random ray selector
# =====================================================
# BEFORE clicking Apply, in the filter properties panel:
#   - Set "Output Data Set Type" to "vtkPolyData"
#
# Steps:
#   1. Open rays_laser0.vtk.series in ParaView
#   2. Filters -> Alphabetical -> Programmable Filter
#   3. Set "Output Data Set Type" dropdown to "vtkPolyData"
#   4. Paste this script into the Script box
#   5. Adjust N_SELECT and SEED, then Apply
#   6. Apply Filters -> Tube on the result to make lines thick

import numpy as np
import vtk

# ── Parameters ────────────────────────────────────────
N_SELECT = 10   # number of rays to show
SEED     = 42   # change for a different random pick
# ──────────────────────────────────────────────────────

pdi = self.GetInputDataObject(0, 0)
pdo = self.GetOutputDataObject(0)

ray_index_arr = pdi.GetPointData().GetArray('rayIndex')

if ray_index_arr is None:
    print("ERROR: rayIndex field not found in this dataset.")
    print("This field requires a rebuilt solver. For now use:")
    print("  Filters -> Threshold on RegionId (after Connectivity filter)")
    pdo.ShallowCopy(pdi)
else:
    n_points = pdi.GetNumberOfPoints()
    ray_indices_np = np.array(
        [int(ray_index_arr.GetValue(i)) for i in range(n_points)]
    )

    unique_rays = np.unique(ray_indices_np)
    rng = np.random.default_rng(SEED)
    n = min(N_SELECT, len(unique_rays))
    selected = set(
        int(x) for x in rng.choice(unique_rays, n, replace=False)
    )
    print(f"Selected ray indices: {sorted(selected)}")

    power_arr = pdi.GetPointData().GetArray('power')

    new_points  = vtk.vtkPoints()
    new_lines   = vtk.vtkCellArray()
    new_power   = vtk.vtkFloatArray(); new_power.SetName('power')
    new_ray_idx = vtk.vtkIntArray();   new_ray_idx.SetName('rayIndex')

    pt_map = {}
    for i in range(n_points):
        if int(ray_index_arr.GetValue(i)) not in selected:
            continue
        new_id = new_points.InsertNextPoint(pdi.GetPoint(i))
        pt_map[i] = new_id
        new_power.InsertNextValue(power_arr.GetValue(i) if power_arr else 0.0)
        new_ray_idx.InsertNextValue(int(ray_index_arr.GetValue(i)))

    for c in range(pdi.GetNumberOfCells()):
        cell = pdi.GetCell(c)
        if cell.GetNumberOfPoints() != 2:
            continue
        p0, p1 = cell.GetPointId(0), cell.GetPointId(1)
        if p0 in pt_map and p1 in pt_map:
            new_lines.InsertNextCell(2)
            new_lines.InsertCellPoint(pt_map[p0])
            new_lines.InsertCellPoint(pt_map[p1])

    pdo.SetPoints(new_points)
    pdo.SetLines(new_lines)
    pdo.GetPointData().AddArray(new_power)
    pdo.GetPointData().AddArray(new_ray_idx)
    pdo.GetPointData().SetActiveScalars('power')

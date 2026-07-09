# Disk cleanup — 2026-07-09

Notes on freeing disk space across the WSL Linux filesystem and the Windows C:\ drive.
Kept for reference in case similar cleanup is needed again.

## WSL filesystem (`tutorials/laserbeamFoam/vdep/`)

Deleted `processor*/` directories for completed testruns whose results were already
reconstructed into `VTKs/` (per [TESTRUNS.md](TESTRUNS.md), processor dirs are root-owned
and must be removed via Docker, not `sudo rm`):

```bash
docker run --rm -v $(pwd):/workspace lbf3 bash -c "
  rm -rf /workspace/tutorials/laserbeamFoam/vdep/CASE/processor*"
```

Applied to testrun36–38, 43–48. Verified each case's `VTKs/` internal max time matched
or exceeded the latest `processor9` timestep before deleting, confirming reconstruction
was current.

**Result:** `vdep/` went from 323 GB → 149 GB (**freed 174 GB**).

## `~/.vscode-server` (WSL side)

- Deleted a 3.0 GB GitHub Copilot Chat local semantic-search index
  (`GitHub.copilot-chat/local-index.1.db`) — rebuilds automatically.
- Cleared `CachedExtensionVSIXs` (541 MB) — extension installer cache, re-downloads as needed.

**Result:** `~/.vscode-server` went from 4.9 GB → 1.4 GB (**freed 3.5 GB**).

## Windows C:\ drive

C:\ was at 1.9 TB used / 13 GB free. The largest item is the WSL Ubuntu virtual disk
(`ext4.vhdx`, ~1004 GB) — **left untouched by choice**, since vhdx files grow but don't
auto-shrink when files are deleted inside WSL, and reclaiming that space wasn't wanted here.

Other cleanup performed directly:

| Item | Freed |
|---|---|
| ParaView crash dumps (`AppData/Local/CrashDumps/*.dmp`) | 35 GB |
| Docker: unused `go-melt` image, dangling `<none>` image, build cache prune | ~66 GB |
| VSCode Roaming `WebStorage` cache (rebuilds automatically) | ~25 GB |
| `AppData/Local/Temp` (Windows temp files) | ~3.9 GB |
| `vscode-remote-wsl/stable` (stale VSCode Server version caches; the actually-running server lives in `~/.vscode-server/bin` on the WSL side) | 3.6 GB |

**Result so far:** C:\ free space went from 13 GB → 77 GB.

### Still outstanding: Docker's internal vhdx

`AppData/Local/Docker/wsl/disk/docker_data.vhdx` is **590 GB on disk** but `docker system df`
shows only ~148 GB of actual image data — same non-shrinking vhdx issue as WSL, just in
Docker Desktop's own hidden distro. Compacting it (via `wsl --shutdown` + `diskpart` from a
plain Windows PowerShell, run outside any WSL session) should reclaim roughly **430–440 GB**.
Steps documented in the conversation; not run automatically since it requires quitting Docker
Desktop and shutting down WSL, which would disconnect any active WSL-backed session.

## Key lesson

`.vhdx` virtual disks (WSL distros, Docker Desktop's backing store) never shrink on their own —
deleting files inside them frees space *inside* the filesystem but not on the Windows host disk.
Periodic compaction via `diskpart` (`compact vdisk`) is needed to actually reclaim that space on C:\.

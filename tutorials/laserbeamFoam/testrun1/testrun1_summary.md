# testrun1 Summary

## Overview

| Item                      | Value                      |
| **Total simulation time** | **0.6 ms** (6e-4 s)       |
| **Wall-clock runtime**    | **~9.5 hours** (34,178 s) |
| **Time steps taken**      | 49,725                     |
| **Written snapshots**     | 20                         |
| **Cores used**            | 54 (hierarchical 3×3×6)   |

## Mesh

| Item             | Value   |
| Base cells       | 90,000  |
| Base cell size   | 10 um   |
| Most refined cell size | 5 um |
| Peak cells (AMR) | 172,789 |
| Final cells      | 172,712 |
| Max allowed      | 200,000 |
| AMR growth factor | ~1.92× |

## Time Stepping

| Item       | Value                            |
| Average Δt | 12.1 ns                          |
| Smallest Δt | 1.9 ns                          |
| Largest Δt | 557 ns                           |
| Throughput (sim-time/wall-time) | ~17.6 ns sim-time per wall-second |

## Laser / Ray Tracing

| Item                  | Value                            |
| Rays per step         | 360 (10 radial × 36 angular)    |
| Deposited power       | avg 173 W / min 84 W / max 194 W |
| Absorption efficiency | ~86% of 200 W incident          |

## Solver Effort (Linear Iterations)

| Equation | Solver invocations | Avg iters. to converge | Total linear iters. | Invocations/step | % of linear iter. cost |
| T        | 1,042,056          | 1.1                    | 1,127,154           | 21.0             | 25.2                   |
| p_rgh    | 149,178            | 17.6                   | 2,624,080           | 3.0              | 58.7                   |
| pcorr    | 9,942              | 72.6                   | 721,939             | 0.2              | 16.1                   |

## Diagnostics

| Item                          | Status   |
| Continuity error (cumulative) | −3.24e-5 |

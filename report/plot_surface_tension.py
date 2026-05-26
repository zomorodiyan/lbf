#!/usr/bin/env python3
"""
Plot σ(T) and dσ/dT(T) from Morohoshi (2013) Fe-C-O model, 316L + PLC.
Marks the temperature at which dσ/dT = 0 (σ maximum) for each η curve.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# ── Bisection root finder ──────────────────────────────────────────────────────
def bisect(f, a, b, tol=0.2):
    for _ in range(60):
        c = (a + b) / 2
        if abs(b - a) < tol:
            return c
        if f(a) * f(c) < 0:
            b = c
        else:
            a = c
    return (a + b) / 2

# ── Dissolution model ──────────────────────────────────────────────────────────
rho_PLC = 1.145
rho_s   = 8.0
h       = 50e-4   # cm
D       = 500e-4  # cm
f_C, f_O = 0.632, 0.281
wC0, wO0  = 0.030, 0.001  # wt%

def compositions(eta):
    x = (eta * rho_PLC * h) / (eta * rho_PLC * h + rho_s * D)
    return x * f_C * 100 + (1-x) * wC0,  x * f_O * 100 + (1-x) * wO0

# ── Morohoshi (2013) ───────────────────────────────────────────────────────────
def K_moro(T, wC, wO):
    return 10**(-0.17*wO - 0.427*wC) * wO * np.exp(29300/T - 10.9)

def sigma_moro(T, wC, wO):
    Kv = K_moro(T, wC, wO)
    return 1925 - 0.455*(T-1808) - 0.155*T*np.log(1+Kv)

def dsigma_moro(T, wC, wO):
    Kv = K_moro(T, wC, wO)
    return -0.455 - 0.155*(np.log(1+Kv) - Kv/(1+Kv)*29300/T)

# ── Cases ──────────────────────────────────────────────────────────────────────
etas   = [0.001, 0.01]
colors = ['#2266CC', '#CC4400']

T = np.linspace(1658, 3068, 1000)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.2), sharex=True)
fig.subplots_adjust(wspace=0.38)

# ── Keyhole wall temperature band ─────────────────────────────────────────────
T_kh_lo, T_kh_hi = 2000, 3068
for ax in (ax1, ax2):
    ax.axvspan(T_kh_lo, T_kh_hi, color='#ffcc88', alpha=0.22, zorder=0)

# ── Plot curves and collect zero crossings ────────────────────────────────────
crossings = []   # (T_star, sigma_star, color, eta)
for eta, color in zip(etas, colors):
    wC, wO = compositions(eta)
    ax1.plot(T, sigma_moro(T, wC, wO),  color=color, lw=2.2, ls='-', label=f"η = {eta*100:g}%")
    ax2.plot(T, dsigma_moro(T, wC, wO), color=color, lw=2.2, ls='-', label=f"η = {eta*100:g}%")

    if dsigma_moro(1650, wC, wO) * dsigma_moro(3100, wC, wO) < 0:
        T_star = bisect(lambda t: dsigma_moro(t, wC, wO), 1650, 3100)
        crossings.append((T_star, sigma_moro(T_star, wC, wO), color, eta))

# ── Green lines + markers at zero-crossing temperatures ───────────────────────
dy1_map = {
    0.001: -65,
    0.01:  +60,
}

for T_star, sig_star, color, eta in crossings:
    pass

# ── Reference temperature lines ────────────────────────────────────────────────
ref_lines = [
    (1658, 'Solidus',  '#222222'),
    (1723, 'Liquidus', '#888888'),
    (3068, 'Boiling',  '#444444'),
]
for ax in (ax1, ax2):
    for Tv, lbl, col in ref_lines:
        ax.axvline(Tv, color=col, lw=1.6, ls='-', alpha=0.85, zorder=4)
    ax.set_xlim(1558, 3168)
    ax.grid(True, alpha=0.25)
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())

ax2.axhline(0, color='black', lw=1.4, ls='-')

# Reference line labels on σ panel only — horizontal, low for solidus/liquidus, high for boiling
ylo, yhi = ax1.get_ylim()
y_solidus  = ylo + (yhi - ylo) * 0.05
y_liquidus = ylo + (yhi - ylo) * 0.14
y_high     = ylo + (yhi - ylo) * 0.42
label_heights   = [y_solidus, y_liquidus, y_high]
label_offsets   = [+8, +8, -8]
label_aligns    = ['left', 'left', 'right']
for (Tv, lbl, col), ypos, xoff, ha in zip(ref_lines, label_heights, label_offsets, label_aligns):
    ax1.text(Tv + xoff, ypos, lbl, fontsize=18, color=col,
             va='bottom', ha=ha, fontweight='bold')

# Keyhole wall band label — placed on dσ/dT panel where vertical space is ample
T_kh_mid = (T_kh_lo + T_kh_hi) / 2
ylo2, yhi2 = ax2.get_ylim()
ax2.text(T_kh_mid, ylo2 + (yhi2 - ylo2) * 0.58, "keyhole wall\ntemperature",
         fontsize=20, color='#aa5500', ha='center', va='top', style='italic')

# ── Legends ────────────────────────────────────────────────────────────────────
handles, labels = ax1.get_legend_handles_labels()
for ax in (ax1, ax2):
    ax.legend(handles, labels, fontsize=20, loc='upper right', ncol=1,
              framealpha=0.92, edgecolor='gray',
              handlelength=1.8, labelspacing=0.4)

ax1.set_ylabel("σ (mN/m)", fontsize=22)
ax2.set_ylabel("dσ/dT (mN/m·K)", fontsize=22)
ax1.set_xlabel("T (K)", fontsize=22)
ax2.set_xlabel("T (K)", fontsize=22)
ax1.tick_params(labelsize=20)
ax2.tick_params(labelsize=20)

out = "/home/mzomoro1/bin/lbf3/report/surface_tension_curves.png"
fig.savefig(out, dpi=150, bbox_inches='tight')
print(f"Saved: {out}")
for T_star, _, _, eta in crossings:
    print(f"  η={eta}: T* = {T_star:.1f} K")
plt.close()

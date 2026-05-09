#!/usr/bin/env python3
"""3D surface mesh + 2D contour map of the co-resistance phi landscape.

Left panel  — 3D surface mesh: hillshaded surface with wireframe overlay
              shows the topology of co-resistance.
Right panel — 2D filled-contour map: top-down view with labelled isolines
              and clear block annotations.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd
from matplotlib.colors import LightSource
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 — registers 3d projection
from scipy.interpolate import RectBivariateSpline

ROOT = Path(__file__).resolve().parents[1]
EDGES_CSV = (
    ROOT / "outputs" / "analysis_outputs"
    / "cross_resistance_network" / "cross_resistance_edges.csv"
)
OUT = ROOT / "outputs" / "final_framework_outputs"

DRUG_ORDER = [
    "Ciprofloxacin",
    "Norfloxacin",
    "Amoxicillin-Clavulanic acid",
    "Ceftriaxone",
    "Ceftazidime",
    "Cefepime",
]
DRUG_LABELS = ["Cipro", "Norflox", "Amox-Clav", "CRO", "CAZ", "FEP"]

matplotlib.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "DejaVu Sans", "Liberation Sans"],
    "font.size": 8,
    "savefig.dpi": 300,
    "savefig.facecolor": "white",
})


def _build_phi_matrix(edges: pd.DataFrame, site: str = "ALL") -> np.ndarray:
    data = edges[edges["site"] == site]
    n = len(DRUG_ORDER)
    mat = np.eye(n, dtype=float)
    idx = {d: i for i, d in enumerate(DRUG_ORDER)}
    for _, row in data.iterrows():
        a = idx.get(row["drug_a"])
        b = idx.get(row["drug_b"])
        if a is not None and b is not None:
            mat[a, b] = float(row["phi"])
            mat[b, a] = float(row["phi"])
    return mat


def _make_smooth_grid(mat: np.ndarray, res: int = 90):
    n = mat.shape[0]
    x_raw = np.arange(n, dtype=float)
    spline = RectBivariateSpline(x_raw, x_raw, mat, kx=3, ky=3)
    xi = np.linspace(0, n - 1, res)
    Z = np.clip(spline(xi, xi), 0.0, 1.0)
    X, Y = np.meshgrid(xi, xi)
    return xi, X, Y, Z


def figure_phi_surface_mesh() -> Path:
    edges = pd.read_csv(EDGES_CSV)
    mat = _build_phi_matrix(edges, site="ALL")
    n = len(DRUG_ORDER)
    xi, X, Y, Z = _make_smooth_grid(mat)

    cmap = plt.cm.Blues

    # ── figure: 2 panels, 3D left, 2D right ──────────────────────────────
    fig = plt.figure(figsize=(13.5, 5.6), facecolor="white")
    gs = gridspec.GridSpec(
        1, 2,
        width_ratios=[1.35, 1],
        wspace=0.08,
        left=0.02, right=0.97,
        top=0.88, bottom=0.08,
    )
    ax3d = fig.add_subplot(gs[0], projection="3d", computed_zorder=False)
    ax2d = fig.add_subplot(gs[1])

    # ══════════════════════════════════════════════════════════════════════
    # LEFT — 3D surface mesh
    # ══════════════════════════════════════════════════════════════════════

    ls = LightSource(azdeg=300, altdeg=40)
    rgb = ls.shade(Z, cmap=cmap, vmin=0, vmax=1, blend_mode="soft")
    ax3d.plot_surface(X, Y, Z, facecolors=rgb, linewidth=0,
                      antialiased=True, shade=False, zorder=3)

    # Invisible surface for colorbar scalar mapping
    _csurf = ax3d.plot_surface(X, Y, Z, cmap=cmap, vmin=0, vmax=1, alpha=0)

    # Wireframe overlay
    step = 9
    ax3d.plot_wireframe(
        X[::step, ::step], Y[::step, ::step], Z[::step, ::step],
        color="white", linewidth=0.28, alpha=0.50, zorder=4,
    )

    # Floor contour shadow
    ax3d.contourf(X, Y, Z, zdir="z", offset=0.0, levels=12,
                  cmap="Blues", alpha=0.28, vmin=0, vmax=1, zorder=1)

    # Axes
    ax3d.set_xticks(range(n))
    ax3d.set_xticklabels(DRUG_LABELS, fontsize=6.5, ha="right")
    ax3d.set_yticks(range(n))
    ax3d.set_yticklabels(DRUG_LABELS, fontsize=6.5, ha="left")
    ax3d.set_zlim(0.0, 1.10)
    ax3d.set_zlabel("Phi (φ)", fontsize=8, labelpad=8)

    for pane in (ax3d.xaxis.pane, ax3d.yaxis.pane, ax3d.zaxis.pane):
        pane.fill = False
        pane.set_edgecolor("#cccccc")
    ax3d.grid(True, linewidth=0.3, color="#dddddd", alpha=0.5)

    # Rotate so both peaks visible from the side (perpendicular to diagonal)
    ax3d.view_init(elev=30, azim=145)
    ax3d.set_box_aspect([1.0, 1.0, 0.50])

    ax3d.set_title("A   3D phi landscape", fontsize=9, fontweight="bold",
                   pad=10, loc="left")

    # ══════════════════════════════════════════════════════════════════════
    # RIGHT — 2D filled-contour map (top-down view)
    # ══════════════════════════════════════════════════════════════════════

    cf = ax2d.contourf(xi, xi, Z, levels=22, cmap="Blues", vmin=0, vmax=1)

    # Labelled isolines at key phi thresholds
    iso_levels = [0.30, 0.50, 0.70, 0.90]
    cs = ax2d.contour(xi, xi, Z, levels=iso_levels,
                      colors="white", linewidths=0.75, alpha=0.85)
    ax2d.clabel(cs, fmt={v: f"φ={v:.2f}" for v in iso_levels},
                fontsize=6.5, inline=True, inline_spacing=4,
                colors="white")

    # Drug axis ticks
    ax2d.set_xticks(range(n))
    ax2d.set_xticklabels(DRUG_LABELS, rotation=35, ha="right", fontsize=7.5)
    ax2d.set_yticks(range(n))
    ax2d.set_yticklabels(DRUG_LABELS, fontsize=7.5)
    ax2d.tick_params(length=0)
    for sp in ax2d.spines.values():
        sp.set_linewidth(0.5)
        sp.set_edgecolor("#aaaaaa")

    # Block annotations — clear white-boxed labels
    _ann_kw = dict(ha="center", va="center", fontsize=8, fontweight="bold",
                   color="white",
                   bbox=dict(boxstyle="round,pad=0.35", facecolor="#1a3a6e",
                             alpha=0.80, edgecolor="none"))
    ax2d.text(0.5, 0.5, "Fluoroquinolone\nblock", **_ann_kw)
    ax2d.text(4.0, 4.0, "Cephalosporin /\nESBL block", **_ann_kw)
    ax2d.text(2.0, 2.0, "Amox-Clav\nintermediate",
              ha="center", va="center", fontsize=7, color="white",
              bbox=dict(boxstyle="round,pad=0.25", facecolor="#2a4a6e",
                        alpha=0.65, edgecolor="none"))

    # Shared colorbar anchored to right panel
    cb = fig.colorbar(cf, ax=ax2d, shrink=0.88, pad=0.03, aspect=22)
    cb.set_label("Phi (φ) correlation", fontsize=7.5)
    cb.ax.tick_params(labelsize=7, length=2)
    cb.outline.set_linewidth(0.5)

    ax2d.set_title("B   Top-down contour map", fontsize=9, fontweight="bold",
                   pad=10, loc="left")

    # ── overall title ─────────────────────────────────────────────────────
    fig.suptitle(
        "Co-resistance phi landscape — all sites\n"
        "Two resistance ecology peaks (fluoroquinolone and cephalosporin/ESBL) "
        "separated by a mid-phi valley",
        fontsize=9.5, fontweight="bold", y=0.98,
    )

    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / "figure_phi_surface_mesh.png"
    fig.savefig(path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return path


if __name__ == "__main__":
    path = figure_phi_surface_mesh()
    print(f"Saved: {path}")

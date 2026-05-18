#!/usr/bin/env python3
"""Build Nature Communications-style manuscript figures and LaTeX tables.

The figures are vector PDFs for Overleaf.  They intentionally avoid decorative
styling: one sans-serif font family, colorblind-aware colors, direct labels, and
source-data-driven panels.
"""

from __future__ import annotations

import math
import re
from pathlib import Path

import pandas as pd
from reportlab.lib import colors
from reportlab.lib.pagesizes import landscape, letter
from reportlab.lib.units import inch
from reportlab.pdfgen import canvas


ROOT = Path(__file__).resolve().parents[1]
FINAL = ROOT / "outputs" / "final_framework_outputs"
ANALYSIS = ROOT / "outputs" / "analysis_outputs"
FIG_DIR = ROOT / "manuscript" / "figures"
TABLE_DIR = ROOT / "manuscript" / "tables"
SOURCE_DIR = ROOT / "manuscript" / "source_data"


BLUE = colors.HexColor("#2B6CB0")
LIGHT_BLUE = colors.HexColor("#A9C8E8")
ORANGE = colors.HexColor("#C56B45")
RED = colors.HexColor("#B64B4B")
GRAY = colors.HexColor("#6B7280")
LIGHT_GRAY = colors.HexColor("#E5E7EB")
MID_GRAY = colors.HexColor("#9CA3AF")
DARK = colors.HexColor("#111827")
PALE = colors.HexColor("#F8FAFC")
GREEN = colors.HexColor("#2F855A")


DRUG_SHORT = {
    "Ciprofloxacin": "Cipro",
    "Norfloxacin": "Norflox",
    "Amoxicillin-Clavulanic acid": "Amox-Clav",
    "Ceftriaxone": "CRO",
    "Ceftazidime": "CAZ",
    "Cefepime": "FEP",
    "Oxacillin": "Oxa",
    "Penicillin": "Pen",
    "Erythromycin": "Ery",
    "E. coli / Cipro": "E. coli / Cipro",
    "E. coli / Amox-Clav": "E. coli / Amox-Clav",
    "S. aureus / Oxacillin": "S. aureus / Oxacillin",
}

DRUG_ORDER = ["Cipro", "Norflox", "Amox-Clav", "CRO", "CAZ", "FEP"]
SAUREUS_DRUG_ORDER = ["Oxa", "Pen", "Cipro", "Ery"]


def ensure_dirs() -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    SOURCE_DIR.mkdir(parents=True, exist_ok=True)


def short_drug(name: str) -> str:
    return DRUG_SHORT.get(str(name), str(name))


def parse_auc_ci(text: str) -> tuple[float, float, float]:
    if text is None or (isinstance(text, float) and math.isnan(text)):
        return math.nan, math.nan, math.nan
    match = re.match(r"\s*([0-9.]+)\s+\(([0-9.]+)-([0-9.]+)\)", str(text))
    if not match:
        try:
            value = float(text)
            return value, math.nan, math.nan
        except ValueError:
            return math.nan, math.nan, math.nan
    return tuple(float(match.group(i)) for i in range(1, 4))


def xmap(value: float, x0: float, x1: float, xmin: float, xmax: float) -> float:
    value = max(xmin, min(xmax, float(value)))
    return x0 + (value - xmin) / (xmax - xmin) * (x1 - x0)


def text(c: canvas.Canvas, x: float, y: float, s: str, size: int = 9, bold: bool = False, fill=DARK) -> None:
    c.setFillColor(fill)
    c.setFont("Helvetica-Bold" if bold else "Helvetica", size)
    c.drawString(x, y, s)


def centered_text(c: canvas.Canvas, x: float, y: float, s: str, size: int = 9, bold: bool = False, fill=DARK) -> None:
    c.setFillColor(fill)
    c.setFont("Helvetica-Bold" if bold else "Helvetica", size)
    c.drawCentredString(x, y, s)


def right_text(c: canvas.Canvas, x: float, y: float, s: str, size: int = 9, bold: bool = False, fill=DARK) -> None:
    c.setFillColor(fill)
    c.setFont("Helvetica-Bold" if bold else "Helvetica", size)
    c.drawRightString(x, y, s)


def panel_label(c: canvas.Canvas, x: float, y: float, label: str) -> None:
    text(c, x, y, label, size=13, bold=True, fill=DARK)


def white_page(c: canvas.Canvas, w: float, h: float) -> None:
    c.setFillColor(colors.white)
    c.rect(0, 0, w, h, stroke=0, fill=1)


def axis_auc(c: canvas.Canvas, x0: float, x1: float, y: float, xmin: float = 0.45, xmax: float = 0.90) -> None:
    c.setStrokeColor(LIGHT_GRAY)
    c.setLineWidth(0.6)
    for tick in [0.5, 0.6, 0.7, 0.8, 0.9]:
        x = xmap(tick, x0, x1, xmin, xmax)
        c.line(x, y, x, y + 205)
        centered_text(c, x, y - 13, f"{tick:.1f}", size=7, fill=GRAY)
    c.setStrokeColor(MID_GRAY)
    c.setLineWidth(1.0)
    c.line(xmap(0.5, x0, x1, xmin, xmax), y, xmap(0.5, x0, x1, xmin, xmax), y + 205)
    text(c, x0, y - 28, "AUC", size=8, fill=GRAY)


def auc_fill(value: float):
    """Color scale for centered AUC in the model-class matrix."""
    if pd.isna(value):
        return colors.white
    value = max(0.45, min(0.90, float(value)))
    t = (value - 0.45) / 0.45
    if value < 0.58:
        r0, g0, b0 = 250, 245, 241
        r1, g1, b1 = 197, 107, 69
    elif value < 0.68:
        r0, g0, b0 = 239, 246, 255
        r1, g1, b1 = 169, 200, 232
        t = (value - 0.58) / 0.10
    else:
        r0, g0, b0 = 219, 235, 248
        r1, g1, b1 = 43, 108, 176
        t = (value - 0.68) / 0.22
    t = max(0.0, min(1.0, t))
    return colors.Color(
        (r0 + (r1 - r0) * t) / 255,
        (g0 + (g1 - g0) * t) / 255,
        (b0 + (b1 - b0) * t) / 255,
    )


def marker(c: canvas.Canvas, x: float, y: float, color, shape: str, hollow: bool = False, size: float = 5) -> None:
    c.setStrokeColor(color)
    c.setFillColor(colors.white if hollow else color)
    c.setLineWidth(1.2)
    if shape == "circle":
        c.circle(x, y, size, stroke=1, fill=1)
    elif shape == "square":
        c.rect(x - size, y - size, 2 * size, 2 * size, stroke=1, fill=1)
    elif shape == "diamond":
        p = c.beginPath()
        p.moveTo(x,          y + size)
        p.lineTo(x + size,   y)
        p.lineTo(x,          y - size)
        p.lineTo(x - size,   y)
        p.close()
        c.drawPath(p, stroke=1, fill=1)
    else:
        c.circle(x, y, size, stroke=1, fill=1)


def figure_1_framework() -> Path:
    path = FIG_DIR / "figure_1_framework.pdf"
    w, h = 7.2 * inch, 3.05 * inch
    c = canvas.Canvas(str(path), pagesize=(w, h))
    white_page(c, w, h)
    text(c, 0.18 * inch, h - 0.32 * inch, "Fig. 1 | Background-matched MALDI-AMR audit", 12, True)
    text(c, 0.18 * inch, h - 0.52 * inch, "A model-agnostic test of whether focal-drug prediction survives co-resistance background control.", 7.8, False, GRAY)

    y = h - 1.45 * inch
    boxes = [
        ("Model predictions", "isolate ID, site/year,\norganism, drug, label,\nprobability score"),
        ("Raw transfer", "External AUC/AUPR\nbefore background control"),
        ("Background match", "Group isolates by labels\nof the other antibiotics"),
        ("Centered audit", "Subtract each stratum's\nmean model score"),
        ("Interpretation", "retained focal signal,\npartial retention, or\nbackground-driven collapse"),
    ]
    x0 = 0.18 * inch
    box_w = 1.16 * inch
    gap = 0.25 * inch
    box_h = 0.86 * inch
    for i, (title, body) in enumerate(boxes):
        x = x0 + i * (box_w + gap)
        c.setFillColor(PALE)
        c.setStrokeColor(LIGHT_GRAY)
        c.roundRect(x, y, box_w, box_h, 8, stroke=1, fill=1)
        text(c, x + 0.08 * inch, y + box_h - 0.20 * inch, title, 6.6, True, DARK)
        for j, line in enumerate(body.split("\n")):
            text(c, x + 0.08 * inch, y + box_h - 0.38 * inch - j * 0.14 * inch, line, 5.8, False, GRAY)
        if i < len(boxes) - 1:
            c.setStrokeColor(MID_GRAY)
            c.setLineWidth(1.0)
            ax = x + box_w + 0.05 * inch
            ay = y + box_h / 2
            c.line(ax, ay, ax + gap - 0.10 * inch, ay)
            c.line(ax + gap - 0.15 * inch, ay + 4, ax + gap - 0.10 * inch, ay)
            c.line(ax + gap - 0.15 * inch, ay - 4, ax + gap - 0.10 * inch, ay)

    panel_label(c, 0.22 * inch, 1.03 * inch, "Audit outputs")
    outputs = [
        ("raw AUC", "apparent external performance"),
        ("matched AUC", "performance in strata with both R and S isolates"),
        ("background-centered AUC", "ranking after stratum-level score shifts are removed"),
        ("matched retention", "how much of the target comparison remains interpretable"),
        ("cross-resistance network", "the label ecology the model could exploit"),
    ]
    for i, (name, desc) in enumerate(outputs):
        yy = 0.77 * inch - i * 0.16 * inch
        c.setFillColor(BLUE if i == 2 else LIGHT_GRAY)
        c.circle(0.34 * inch, yy + 2, 2.2, stroke=0, fill=1)
        text(c, 0.46 * inch, yy, name, 6.5, True, DARK)
        text(c, 2.02 * inch, yy, desc, 6.5, False, GRAY)

    c.showPage()
    c.save()
    return path


def figure_2_primary_audit(primary: pd.DataFrame) -> Path:
    path = FIG_DIR / "figure_2_primary_background_audit.pdf"
    c = canvas.Canvas(str(path), pagesize=landscape(letter))
    w, h = landscape(letter)
    white_page(c, w, h)
    text(c, 0.55 * inch, h - 0.45 * inch, "Fig. 2 | Focal-drug signal after background matching", 14, True)
    text(c, 0.55 * inch, h - 0.70 * inch, "Raw external AUC is compared with stratum-centered AUC. Centering removes score shifts shared by co-resistance strata.", 9, False, GRAY)

    data = _prep_auc_ci(primary)
    panels = [
        ("E. coli / Cipro", BLUE, 0.55 * inch,
         "Ciprofloxacin retains residual within-background ranking."),
        ("E. coli / Amox-Clav", ORANGE, 5.65 * inch,
         "Amox-Clav weakens or collapses after background control."),
    ]
    for pair, col, px, subtitle in panels:
        _pair_panel(c, pair, subtitle, col, px, h, data)
    c.showPage()
    c.save()
    return path


def model_class_column(row: pd.Series) -> str:
    if row["model_class"] == "CNN/Mega":
        return "CNN/Mega"
    if row["model_class"] == "LightGBM" and row["model_variant"] == "multi-task":
        return "LGBM multi"
    if row["model_class"] == "LightGBM" and row["model_variant"] == "single-task":
        return "LGBM single"
    if row["model_class"] == "Weis LR":
        return "Weis LR"
    return f"{row['model_class']} {row['model_variant']}".strip()


def figure_3_matrix_rows(matrix_df: pd.DataFrame) -> pd.DataFrame:
    matrix = matrix_df[matrix_df["status"].eq("complete")].copy()
    targets = [
        ("E. coli", "Ciprofloxacin", "A-2018"),
        ("E. coli", "Ciprofloxacin", "DRIAMS-C"),
        ("E. coli", "Ciprofloxacin", "DRIAMS-D"),
        ("E. coli", "Amoxicillin-Clavulanic acid", "A-2018"),
        ("E. coli", "Amoxicillin-Clavulanic acid", "DRIAMS-C"),
        ("E. coli", "Amoxicillin-Clavulanic acid", "DRIAMS-D"),
        ("S. aureus", "Oxacillin", "A-2018"),
        ("S. aureus", "Oxacillin", "DRIAMS-B"),
        ("S. aureus", "Oxacillin", "DRIAMS-C"),
    ]
    target_index = {target: index for index, target in enumerate(targets)}

    def organism_short(organism: str) -> str:
        if organism == "Escherichia coli":
            return "E. coli"
        if organism == "Staphylococcus aureus":
            return "S. aureus"
        return organism

    matrix["organism_short"] = matrix["organism"].map(organism_short)
    matrix["target_key"] = list(zip(matrix["organism_short"], matrix["drug"], matrix["site"]))
    matrix = matrix[matrix["target_key"].isin(target_index)].copy()
    matrix["target_order"] = matrix["target_key"].map(target_index)
    matrix["display_target"] = matrix.apply(
        lambda row: f"{row['organism_short']} / {short_drug(row['drug'])} / {row['site']}",
        axis=1,
    )
    matrix["model_column"] = matrix.apply(model_class_column, axis=1)
    matrix = matrix[matrix["model_column"].isin(["CNN/Mega", "LGBM multi", "LGBM single"])].copy()
    for col in ["raw_auc", "centered_auc", "matched_retention"]:
        matrix[col] = pd.to_numeric(matrix[col], errors="coerce")
    return matrix.sort_values(["target_order", "model_column"])


def figure_3_model_replication(model_df: pd.DataFrame, matrix_df: pd.DataFrame | None = None) -> Path:
    path = FIG_DIR / "figure_3_model_family_replication.pdf"
    if matrix_df is not None and not matrix_df.empty:
        focus = figure_3_matrix_rows(matrix_df)
        c = canvas.Canvas(str(path), pagesize=landscape(letter))
        w, h = landscape(letter)
        white_page(c, w, h)
        text(c, 0.55 * inch, h - 0.45 * inch, "Fig. 3 | Completed model-class background-audit matrix", 14, True)
        text(
            c,
            0.55 * inch,
            h - 0.70 * inch,
            "CNN/Mega, LightGBM multi-task and LightGBM single-task models are audited with the same co-resistance matching framework.",
            9,
            False,
            GRAY,
        )

        columns = ["CNN/Mega", "LGBM multi", "LGBM single"]
        x_label = 0.70 * inch
        x0 = 3.20 * inch
        cell_w = 1.80 * inch
        cell_h = 0.34 * inch
        y_top = h - 1.30 * inch
        centered_text(c, x_label + 1.00 * inch, y_top + 0.18 * inch, "target / site", 8.2, True, DARK)
        for i, col in enumerate(columns):
            centered_text(c, x0 + i * cell_w + cell_w / 2, y_top + 0.18 * inch, col, 8.2, True, DARK)

        row_targets = focus[["target_order", "display_target"]].drop_duplicates().sort_values("target_order")
        for row_i, target in enumerate(row_targets.itertuples(index=False)):
            y = y_top - (row_i + 1) * cell_h
            if row_i in {3, 6}:
                c.setStrokeColor(LIGHT_GRAY)
                c.setLineWidth(0.7)
                c.line(0.58 * inch, y + cell_h + 2, w - 0.60 * inch, y + cell_h + 2)
            text(c, x_label, y + 0.10 * inch, target.display_target, 7.6, False, DARK)
            row_subset = focus[focus["target_order"].eq(target.target_order)].set_index("model_column")
            for col_i, col in enumerate(columns):
                x = x0 + col_i * cell_w
                if col in row_subset.index:
                    record = row_subset.loc[col]
                    raw = float(record["raw_auc"])
                    centered = float(record["centered_auc"])
                    caution = "caution" in str(record["adequacy_label"])
                    c.setFillColor(auc_fill(centered))
                    c.setStrokeColor(MID_GRAY if caution else LIGHT_GRAY)
                    c.rect(x, y, cell_w - 0.04 * inch, cell_h - 0.04 * inch, stroke=1, fill=1)
                    label = f"{raw:.2f}->{centered:.2f}"
                    centered_text(c, x + (cell_w - 0.04 * inch) / 2, y + 0.10 * inch, label, 7.1, False, DARK)
                    if caution:
                        centered_text(c, x + cell_w - 0.14 * inch, y + 0.10 * inch, "*", 8.0, True, RED)
                else:
                    c.setFillColor(PALE)
                    c.setStrokeColor(LIGHT_GRAY)
                    c.rect(x, y, cell_w - 0.04 * inch, cell_h - 0.04 * inch, stroke=1, fill=1)
                    centered_text(c, x + (cell_w - 0.04 * inch) / 2, y + 0.10 * inch, "not run", 7.0, False, GRAY)

        text(
            c,
            0.70 * inch,
            0.86 * inch,
            "Cells show raw AUC -> background-centered AUC; color encodes centered AUC. * denotes low matched support.",
            7.2,
            False,
            GRAY,
        )
        c.setFillColor(auc_fill(0.50))
        c.rect(7.45 * inch, 0.80 * inch, 0.25 * inch, 0.13 * inch, stroke=0, fill=1)
        text(c, 7.76 * inch, 0.80 * inch, "chance/weak", 6.8, False, GRAY)
        c.setFillColor(auc_fill(0.65))
        c.rect(8.55 * inch, 0.80 * inch, 0.25 * inch, 0.13 * inch, stroke=0, fill=1)
        text(c, 8.86 * inch, 0.80 * inch, "moderate", 6.8, False, GRAY)
        c.setFillColor(auc_fill(0.78))
        c.rect(9.42 * inch, 0.80 * inch, 0.25 * inch, 0.13 * inch, stroke=0, fill=1)
        text(c, 9.73 * inch, 0.80 * inch, "retained", 6.8, False, GRAY)
        c.showPage()
        c.save()
        return path

    c = canvas.Canvas(str(path), pagesize=landscape(letter))
    w, h = landscape(letter)
    white_page(c, w, h)
    text(c, 0.55 * inch, h - 0.45 * inch, "Fig. 3 | Background sensitivity across model families", 14, True)
    text(c, 0.55 * inch, h - 0.70 * inch, "CNN and multi-task LGBM outputs are audited using the same co-resistance strata.", 9, False, GRAY)

    keep = model_df[model_df["drug"].isin(["Cipro", "Amox-Clav"]) & model_df["site"].isin(["A-2018", "DRIAMS-C", "DRIAMS-D"])].copy()
    keep["drug_order"] = keep["drug"].map({"Cipro": 0, "Amox-Clav": 1})
    keep["site_order"] = keep["site"].map({"A-2018": 0, "DRIAMS-C": 1, "DRIAMS-D": 2})
    keep = keep.sort_values(["drug_order", "site_order"])

    panels = [("CNN", "cnn_raw_auc", "cnn_centered_auc", 0.75 * inch), ("LGBM multi", "lgbm_raw_auc", "lgbm_centered_auc", 5.70 * inch)]
    xmin, xmax = 0.45, 0.90
    for model, raw_col, cen_col, px in panels:
        panel_label(c, px, h - 1.15 * inch, model)
        x0, x1 = px + 1.05 * inch, px + 4.05 * inch
        y0 = 1.55 * inch
        axis_auc(c, x0, x1, y0, xmin, xmax)
        for i, row in enumerate(keep.itertuples(index=False)):
            y = h - 1.75 * inch - i * 0.37 * inch
            color = BLUE if row.drug == "Cipro" else ORANGE
            right_text(c, x0 - 0.18 * inch, y - 3, f"{row.drug} {row.site}", 7.2, False, DARK)
            raw = float(getattr(row, raw_col))
            cen = float(getattr(row, cen_col))
            c.setStrokeColor(color)
            c.setLineWidth(1.1)
            c.line(xmap(raw, x0, x1, xmin, xmax), y, xmap(cen, x0, x1, xmin, xmax), y)
            marker(c, xmap(raw, x0, x1, xmin, xmax), y, color, "circle", size=4.4)
            marker(c, xmap(cen, x0, x1, xmin, xmax), y, color, "square", size=4.4)
    c.showPage()
    c.save()
    return path


def color_for_phi(phi: float):
    if math.isnan(phi):
        return colors.white
    t = max(0.0, min(1.0, phi))
    r0, g0, b0 = 247, 250, 252
    r1, g1, b1 = 43, 108, 176
    return colors.Color((r0 + (r1 - r0) * t) / 255, (g0 + (g1 - g0) * t) / 255, (b0 + (b1 - b0) * t) / 255)


def figure_4_cross_resistance() -> Path:
    path = FIG_DIR / "figure_4_cross_resistance_network.pdf"
    edges = pd.read_csv(ANALYSIS / "cross_resistance_network" / "cross_resistance_edges.csv")
    all_edges = edges[edges["site"].eq("ALL")]
    matrix = pd.DataFrame(index=DRUG_ORDER, columns=DRUG_ORDER, data=float("nan"))
    for drug in DRUG_ORDER:
        matrix.loc[drug, drug] = 1.0
    for _, row in all_edges.iterrows():
        a, b = short_drug(row["drug_a"]), short_drug(row["drug_b"])
        if a in matrix.index and b in matrix.columns:
            matrix.loc[a, b] = row["phi"]
            matrix.loc[b, a] = row["phi"]

    c = canvas.Canvas(str(path), pagesize=landscape(letter))
    w, h = landscape(letter)
    white_page(c, w, h)
    text(c, 0.55 * inch, h - 0.45 * inch, "Fig. 4 | Co-resistance blocks define exploitable background", 14, True)
    text(c, 0.55 * inch, h - 0.70 * inch, "Phi correlations across E. coli resistance labels show that focal drug labels are embedded in structured resistance ecology.", 9, False, GRAY)

    cell = 0.55 * inch
    x0 = 1.35 * inch
    y0 = h - 1.55 * inch
    for j, drug in enumerate(DRUG_ORDER):
        centered_text(c, x0 + j * cell + cell / 2, y0 + 0.15 * inch, drug, 8, True, DARK)
        right_text(c, x0 - 0.12 * inch, y0 - j * cell - cell / 2 - 3, drug, 8, True, DARK)
    for i, drug_i in enumerate(DRUG_ORDER):
        for j, drug_j in enumerate(DRUG_ORDER):
            val = float(matrix.loc[drug_i, drug_j])
            x = x0 + j * cell
            y = y0 - (i + 1) * cell
            c.setFillColor(color_for_phi(val))
            c.setStrokeColor(colors.white)
            c.rect(x, y, cell, cell, stroke=1, fill=1)
            if not math.isnan(val):
                centered_text(c, x + cell / 2, y + cell / 2 - 3, f"{val:.2f}", 7, False, DARK if val < 0.65 else colors.white)

    text(c, 5.35 * inch, h - 1.45 * inch, "Strongest resistant blocks", 10, True, DARK)
    top = all_edges.sort_values("phi", ascending=False).head(5)
    for i, row in enumerate(top.itertuples(index=False)):
        yy = h - 1.80 * inch - i * 0.35 * inch
        text(c, 5.35 * inch, yy, f"{short_drug(row.drug_a)} / {short_drug(row.drug_b)}", 8.5, True, DARK)
        text(c, 6.70 * inch, yy, f"phi={row.phi:.3f}", 8.5, False, GRAY)
        text(c, 7.50 * inch, yy, f"RR lift={row.rr_lift:.2f}", 8.5, False, GRAY)

    text(c, 5.35 * inch, 1.50 * inch, "Interpretation", 10, True, DARK)
    text(c, 5.35 * inch, 1.25 * inch, "Models can exploit spectra associated with these", 8.5, False, GRAY)
    text(c, 5.35 * inch, 1.05 * inch, "resistant blocks without learning a focal drug in isolation.", 8.5, False, GRAY)
    c.showPage()
    c.save()
    return path


def _pair_panel(c: canvas.Canvas, pair_label: str, subtitle: str, col,
                px: float, h: float, data: pd.DataFrame,
                xmin: float = 0.45, xmax: float = 0.90) -> None:
    """Draw a single pair panel (raw vs centered AUC per site) at position px."""
    panel_label(c, px, h - 1.15 * inch, pair_label)
    text(c, px, h - 1.38 * inch, subtitle, 8, False, GRAY)
    x0, x1 = px + 1.15 * inch, px + 4.55 * inch
    y0 = 1.55 * inch
    axis_auc(c, x0, x1, y0, xmin, xmax)
    sub = data[data["pair"].eq(pair_label)].copy()
    sub["site_order"] = sub["site"].map({"A-2018": 0, "DRIAMS-B": 1, "DRIAMS-C": 2, "DRIAMS-D": 3})
    sub = sub.sort_values("site_order")
    for i, row in enumerate(sub.itertuples(index=False)):
        y = h - 1.88 * inch - i * 0.47 * inch
        caution = "caution" in str(row.adequacy)
        right_text(c, x0 - 0.18 * inch, y - 3, row.site, 8, False, GRAY if caution else DARK)
        raw = getattr(row, "raw_auc_95ci_mid")
        cen = getattr(row, "stratum_centered_auc_95ci_mid")
        raw_x = xmap(raw, x0, x1, xmin, xmax)
        cen_x = xmap(cen, x0, x1, xmin, xmax)
        bar_col = col if not caution else MID_GRAY
        c.setStrokeColor(bar_col)
        c.setLineWidth(1.2)
        c.line(raw_x, y, cen_x, y)
        marker(c, raw_x, y, bar_col, "circle", hollow=False)
        marker(c, cen_x, y, bar_col, "square", hollow=caution)
        c.setStrokeColor(bar_col)
        c.setLineWidth(0.45)
        for lo_attr, hi_attr in [
            ("raw_auc_95ci_low", "raw_auc_95ci_high"),
            ("stratum_centered_auc_95ci_low", "stratum_centered_auc_95ci_high"),
        ]:
            lo_val = getattr(row, lo_attr)
            hi_val = getattr(row, hi_attr)
            if not math.isnan(lo_val) and not math.isnan(hi_val):
                lx = xmap(lo_val, x0, x1, xmin, xmax)
                hx = xmap(hi_val, x0, x1, xmin, xmax)
                c.line(lx, y, hx, y)
                for cx in [lx, hx]:
                    c.line(cx, y - 2.5, cx, y + 2.5)
        delta = getattr(row, "raw_to_centered_delta")
        text(c, x1 + 0.12 * inch, y - 3,
             f"Δ={delta:+.3f}; ret={row.matched_retention_pct:.1f}%", 7.2, False, GRAY)
    marker(c, px + 0.10 * inch, y0 - 0.18 * inch, DARK, "circle")
    text(c, px + 0.25 * inch, y0 - 0.22 * inch, "raw", 7.5, False, GRAY)
    marker(c, px + 0.75 * inch, y0 - 0.18 * inch, DARK, "square")
    text(c, px + 0.90 * inch, y0 - 0.22 * inch, "background-centered", 7.5, False, GRAY)


def _prep_auc_ci(data: pd.DataFrame) -> pd.DataFrame:
    data = data.copy()
    for col in ["raw_auc_95ci", "stratum_centered_auc_95ci"]:
        values = data[col].map(parse_auc_ci)
        data[col + "_mid"]  = values.map(lambda t: t[0])
        data[col + "_low"]  = values.map(lambda t: t[1])
        data[col + "_high"] = values.map(lambda t: t[2])
    return data


def _build_three_way_data(
    primary_df: pd.DataFrame,
    sa_summary_df: pd.DataFrame,
    ecoli_bg_df: pd.DataFrame,
    saureus_bg_df: pd.DataFrame,
) -> pd.DataFrame:
    """Combine raw MALDI, co-resistance-only, and background-centred AUC for figure 6."""
    data = _prep_auc_ci(primary_df)
    rows = []

    for pair_label, drug_long in [
        ("E. coli / Cipro",     "Ciprofloxacin"),
        ("E. coli / Amox-Clav", "Amoxicillin-Clavulanic acid"),
    ]:
        sub    = data[data["pair"].eq(pair_label)].copy()
        bg_sub = ecoli_bg_df[ecoli_bg_df["drug"].eq(drug_long)].copy()
        for _, row in sub.iterrows():
            site   = row["site"]
            bg_row = bg_sub[bg_sub["site"].eq(site)]
            bg_auc = float(bg_row["exact_background_auc"].iloc[0]) if len(bg_row) else float("nan")
            rows.append({
                "pair":        pair_label,
                "site":        site,
                "raw_auc":     row["raw_auc_95ci_mid"],
                "raw_lo":      row["raw_auc_95ci_low"],
                "raw_hi":      row["raw_auc_95ci_high"],
                "bg_only_auc": bg_auc,
                "centered_auc": row["stratum_centered_auc_95ci_mid"],
                "centered_lo": row["stratum_centered_auc_95ci_low"],
                "centered_hi": row["stratum_centered_auc_95ci_high"],
                "caution":     "caution" in str(row["adequacy"]),
            })

    sa_oxa  = sa_summary_df[sa_summary_df["drug"].eq("Oxacillin")].copy()
    bg_sa   = saureus_bg_df[saureus_bg_df["drug"].eq("Oxacillin")].copy()
    for _, row in sa_oxa.iterrows():
        site    = row["site"]
        bg_row  = bg_sa[bg_sa["site"].eq(site)]
        bg_auc  = float(bg_row["exact_background_auc"].iloc[0]) if len(bg_row) else float("nan")
        has_str = int(row["n_valid_strata"]) > 0
        rows.append({
            "pair":        "S. aureus / Oxacillin",
            "site":        site,
            "raw_auc":     float(row["raw_auc"]),
            "raw_lo":      float(row["raw_auc_ci_low"]),
            "raw_hi":      float(row["raw_auc_ci_high"]),
            "bg_only_auc": bg_auc,
            "centered_auc": float(row["stratum_centered_auc"])        if has_str else float("nan"),
            "centered_lo":  float(row["stratum_centered_auc_ci_low"]) if has_str else float("nan"),
            "centered_hi":  float(row["stratum_centered_auc_ci_high"]) if has_str else float("nan"),
            "caution":     "caution" in str(row["adequacy_label"]),
        })

    return pd.DataFrame(rows)


def figure_6_three_way_decomposition(
    primary_df: pd.DataFrame,
    sa_summary_df: pd.DataFrame,
    ecoli_bg_df: pd.DataFrame,
    saureus_bg_df: pd.DataFrame,
) -> Path:
    """Three-way AUC decomposition: raw MALDI, co-resistance-only, background-centred.

    Filled circle  = raw MALDI AUC
    Hollow square  = co-resistance-only AUC (no spectra, AST background only)
    Filled diamond = background-centred MALDI AUC
    """
    path = FIG_DIR / "figure_6_three_way_decomposition.pdf"
    c = canvas.Canvas(str(path), pagesize=landscape(letter))
    w, h = landscape(letter)
    white_page(c, w, h)
    text(c, 0.50 * inch, h - 0.44 * inch,
         "Fig. 6 | Three-way AUC decomposition: raw MALDI, co-resistance-only, background-centred", 13, True)
    text(c, 0.50 * inch, h - 0.67 * inch,
         "Filled circle = raw MALDI AUC;  hollow square = co-resistance-only AUC (no spectra);  filled diamond = background-centred MALDI AUC.  * = sparse matched support (cautionary).",
         8.5, False, GRAY)

    GREEN_SA = colors.HexColor("#2F855A")
    xmin, xmax = 0.44, 1.02
    site_order_map = {"A-2018": 0, "DRIAMS-B": 1, "DRIAMS-C": 2, "DRIAMS-D": 3}
    df = _build_three_way_data(primary_df, sa_summary_df, ecoli_bg_df, saureus_bg_df)

    panels_cfg = [
        ("E. coli / Cipro",       BLUE,     0.50 * inch,
         "MALDI exceeds co-resistance-only at DRIAMS-D;\ncentred AUC above chance at all interpretable sites."),
        ("E. coli / Amox-Clav",   ORANGE,   3.85 * inch,
         "Co-resistance-only matches or beats MALDI at external\nsites; centred AUC collapses to chance."),
        ("S. aureus / Oxacillin", GREEN_SA, 7.20 * inch,
         "MALDI exceeds co-resistance-only at A-2018 (+0.076)\nand DRIAMS-C (+0.098); centred AUC retained at both."),
    ]

    for pair_label, col, px, subtitle in panels_cfg:
        panel_label(c, px, h - 1.10 * inch, pair_label)
        for j, line in enumerate(subtitle.split("\n")):
            text(c, px, h - 1.32 * inch - j * 0.155 * inch, line, 7.5, False, GRAY)

        x0, x1 = px + 0.92 * inch, px + 3.10 * inch
        y0 = 1.55 * inch
        chart_h = 210

        c.setStrokeColor(LIGHT_GRAY)
        c.setLineWidth(0.6)
        for tick in [0.5, 0.6, 0.7, 0.8, 0.9, 1.0]:
            xt = xmap(tick, x0, x1, xmin, xmax)
            c.line(xt, y0, xt, y0 + chart_h)
            centered_text(c, xt, y0 - 13, f"{tick:.1f}", size=7, fill=GRAY)
        c.setStrokeColor(MID_GRAY)
        c.setLineWidth(1.0)
        xc = xmap(0.5, x0, x1, xmin, xmax)
        c.line(xc, y0, xc, y0 + chart_h)
        text(c, x0 - 0.05 * inch, y0 - 27, "AUC", size=8, fill=GRAY)

        sub = df[df["pair"].eq(pair_label)].copy()
        sub["_ord"] = sub["site"].map(site_order_map)
        sub = sub.sort_values("_ord")

        for i, row in enumerate(sub.itertuples(index=False)):
            y      = h - 1.78 * inch - i * 0.52 * inch
            caution = bool(row.caution)
            bar_col = col if not caution else MID_GRAY
            slabel  = row.site + ("*" if caution else "")
            right_text(c, x0 - 0.08 * inch, y - 3, slabel, 8, False,
                       GRAY if caution else DARK)

            raw = float(row.raw_auc)
            bg  = float(row.bg_only_auc)
            cen = float(row.centered_auc)

            # Span line across all three values
            vals = [v for v in [raw, bg, cen] if not math.isnan(v)]
            if len(vals) >= 2:
                lx = xmap(min(vals), x0, x1, xmin, xmax)
                rx = xmap(max(vals), x0, x1, xmin, xmax)
                c.setStrokeColor(LIGHT_GRAY if caution else bar_col)
                c.setLineWidth(0.9)
                c.line(lx, y, rx, y)

            if not math.isnan(raw):
                marker(c, xmap(raw, x0, x1, xmin, xmax), y, bar_col, "circle")
            if not math.isnan(bg):
                marker(c, xmap(bg, x0, x1, xmin, xmax), y, bar_col, "square", hollow=True)
            if not math.isnan(cen):
                marker(c, xmap(cen, x0, x1, xmin, xmax), y, bar_col, "diamond")

    # Shared legend
    lx, ly = 0.50 * inch, 0.92 * inch
    marker(c, lx + 0.10 * inch, ly, DARK, "circle")
    text(c, lx + 0.22 * inch, ly - 4, "Raw MALDI AUC", 7.5, False, GRAY)
    marker(c, lx + 1.60 * inch, ly, DARK, "square", hollow=True)
    text(c, lx + 1.72 * inch, ly - 4, "Co-resistance-only AUC", 7.5, False, GRAY)
    marker(c, lx + 3.40 * inch, ly, DARK, "diamond")
    text(c, lx + 3.52 * inch, ly - 4, "Background-centred AUC", 7.5, False, GRAY)
    text(c, lx, 0.55 * inch,
         "Co-resistance-only baseline: predicts focal drug from co-resistance AST signature only, without MALDI spectra (exact-background smoothed prevalence, leave-one-out).",
         7, False, GRAY)

    c.showPage()
    c.save()
    return path


COMPARISON_ROWS = [
    # (dimension, prior_maldi_amr, tripod_ai, this_framework, note)
    ("Reports external-site AUC",               "Yes",  "Required",    "Yes",  ""),
    ("Controls for co-resistance background",   "No",   "Not required","Yes",  "Core audit metric"),
    ("Quantifies signal collapse vs. retention","No",   "No",          "Yes",  "raw − centered delta"),
    ("Site-specific evaluation",                "Yes",  "Recommended", "Yes",  ""),
    ("Multi-drug panel (> 2 pairs)",            "Some", "N/A",         "Yes",  "saureus_panel, ecoli_mechanism6"),
    ("Second organism validation",              "No",   "N/A",         "Yes",  "E. coli + S. aureus"),
    ("Sensitivity analysis (stratum threshold)","No",   "No",          "Yes",  "2, 3, 5, 10 isolates per stratum"),
    ("Model-agnostic (accepts any CSV)",        "No",   "N/A",         "Yes",  "No DRIAMS/PyTorch dependency"),
    ("Open code + reproducible example data",   "Partial","N/A",       "Yes",  "github + SCHEMA.md"),
]


def figure_7_framework_comparison() -> Path:
    """Comparison table of this framework vs prior MALDI-AMR work and TRIPOD+AI."""
    path = FIG_DIR / "figure_7_framework_comparison.pdf"
    w, h = 7.8 * inch, 4.20 * inch
    c = canvas.Canvas(str(path), pagesize=(w, h))
    white_page(c, w, h)
    text(c, 0.18 * inch, h - 0.32 * inch,
         "Fig. 7 | Framework comparison", 12, True)
    text(c, 0.18 * inch, h - 0.52 * inch,
         "What the background-matched audit adds relative to prior MALDI-AMR publications and TRIPOD+AI.",
         7.8, False, GRAY)

    col_xs = [0.18 * inch, 3.10 * inch, 4.65 * inch, 5.88 * inch]
    col_ws = [2.88 * inch, 1.52 * inch, 1.20 * inch, 1.80 * inch]
    headers = ["Evaluation dimension", "Prior MALDI-AMR", "TRIPOD+AI", "This framework"]

    y_header = h - 0.82 * inch
    for i, (hdr, cx) in enumerate(zip(headers, col_xs)):
        text(c, cx, y_header, hdr, 7.5, True, DARK)

    c.setStrokeColor(MID_GRAY)
    c.setLineWidth(0.5)
    c.line(0.18 * inch, y_header - 5, w - 0.18 * inch, y_header - 5)

    for row_idx, (dim, prior, tripod, ours, note) in enumerate(COMPARISON_ROWS):
        y = y_header - 0.34 * inch - row_idx * 0.34 * inch
        bg = PALE if row_idx % 2 == 0 else colors.white
        c.setFillColor(bg)
        c.rect(0.12 * inch, y - 0.09 * inch, w - 0.24 * inch, 0.30 * inch, stroke=0, fill=1)

        text(c, col_xs[0], y, dim, 7.2, False, DARK)
        for val, cx in [(prior, col_xs[1]), (tripod, col_xs[2]), (ours, col_xs[3])]:
            is_yes = val.lower().startswith("yes")
            is_no  = val.lower().startswith("no")
            fill_col = GREEN if is_yes else (RED if is_no else GRAY)
            text(c, cx, y, val, 7.2, is_yes, fill_col)

    c.setStrokeColor(LIGHT_GRAY)
    c.setLineWidth(0.4)
    c.line(0.18 * inch, h - 0.82 * inch - len(COMPARISON_ROWS) * 0.34 * inch - 0.08 * inch,
           w - 0.18 * inch,
           h - 0.82 * inch - len(COMPARISON_ROWS) * 0.34 * inch - 0.08 * inch)
    text(c, 0.18 * inch, 0.25 * inch,
         "Prior MALDI-AMR: Weis et al. 2022 (Nat. Med.) and comparable publications. "
         "TRIPOD+AI: transparent reporting guideline (2024).",
         6.5, False, GRAY)
    c.showPage()
    c.save()
    return path


def figure_5_public_support(wgs: pd.DataFrame, enrichment: pd.DataFrame) -> Path:
    path = FIG_DIR / "figure_5_public_wgs_proteomic_support.pdf"
    c = canvas.Canvas(str(path), pagesize=landscape(letter))
    w, h = landscape(letter)
    white_page(c, w, h)
    text(c, 0.55 * inch, h - 0.45 * inch, "Fig. 5 | Public WGS-linked MALDI data support lineage encoding", 14, True)
    text(c, 0.55 * inch, h - 0.70 * inch, "A public Basel UPEC dataset links Bruker MALDI spectra, WGS-derived lineage labels, and susceptibility phenotypes.", 9, False, GRAY)

    panel_label(c, 0.65 * inch, h - 1.15 * inch, "a")
    text(c, 0.90 * inch, h - 1.15 * inch, "MALDI peak-feature AUC", 10, True, DARK)
    x0, x1 = 2.45 * inch, 4.85 * inch
    y_start = h - 1.65 * inch
    for tick in [0.5, 0.6, 0.7, 0.8, 0.9, 1.0]:
        x = xmap(tick, x0, x1, 0.5, 1.0)
        c.setStrokeColor(LIGHT_GRAY)
        c.line(x, y_start - 0.2 * inch, x, y_start - 1.85 * inch)
        centered_text(c, x, y_start - 2.02 * inch, f"{tick:.1f}", 7, False, GRAY)
    for i, row in enumerate(wgs.itertuples(index=False)):
        y = y_start - i * 0.48 * inch
        label = str(row.target).replace("_", " ")
        right_text(c, x0 - 0.15 * inch, y - 4, label, 8.5, False, DARK)
        color = BLUE if row.target == "ST131" else ORANGE
        c.setFillColor(color)
        c.rect(x0, y - 0.07 * inch, xmap(row.auc, x0, x1, 0.5, 1.0) - x0, 0.14 * inch, stroke=0, fill=1)
        text(c, xmap(row.auc, x0, x1, 0.5, 1.0) + 0.05 * inch, y - 4, f"{row.auc:.3f}", 8, True, DARK)

    panel_label(c, 5.25 * inch, h - 1.15 * inch, "b")
    text(c, 5.50 * inch, h - 1.15 * inch, "Overlap with published ST131 biomarkers", 10, True, DARK)
    x0, x1 = 7.45 * inch, 10.15 * inch
    y_start = h - 1.65 * inch
    max_fold = 3.4
    for tick in [1, 2, 3]:
        x = xmap(tick, x0, x1, 0, max_fold)
        c.setStrokeColor(LIGHT_GRAY)
        c.line(x, y_start - 0.2 * inch, x, y_start - 1.85 * inch)
        centered_text(c, x, y_start - 2.02 * inch, f"{tick}x", 7, False, GRAY)
    show = enrichment[enrichment["target"].isin(["ST131", "Ciprofloxacin_R", "Ceftriaxone_R"])].copy()
    for i, row in enumerate(show.itertuples(index=False)):
        y = y_start - i * 0.48 * inch
        label = str(row.target).replace("_", " ")
        right_text(c, x0 - 0.15 * inch, y - 4, label, 8.5, False, DARK)
        color = BLUE if row.target == "ST131" else ORANGE
        c.setFillColor(color)
        c.rect(x0, y - 0.07 * inch, xmap(row.fold_enrichment, x0, x1, 0, max_fold) - x0, 0.14 * inch, stroke=0, fill=1)
        text(
            c,
            xmap(row.fold_enrichment, x0, x1, 0, max_fold) + 0.05 * inch,
            y - 4,
            f"{row.fold_enrichment:.2f}x; p={row.empirical_p_ge_observed:.4f}",
            7.5,
            False,
            DARK,
        )
    c.showPage()
    c.save()
    return path


def latex_escape(value: object) -> str:
    text_value = "" if value is None or (isinstance(value, float) and math.isnan(value)) else str(value)
    replacements = {
        "\\": r"\textbackslash{}",
        "&": r"\&",
        "%": r"\%",
        "$": r"\$",
        "#": r"\#",
        "_": r"\_",
        "{": r"\{",
        "}": r"\}",
        "~": r"\textasciitilde{}",
        "^": r"\textasciicircum{}",
    }
    for old, new in replacements.items():
        text_value = text_value.replace(old, new)
    return text_value


def write_latex_table(path: Path, caption: str, label: str, df: pd.DataFrame, notes: str = "") -> None:
    cols = list(df.columns)
    align = "l" + "c" * (len(cols) - 1)
    lines = [
        r"\begin{table}[ht]",
        r"\centering",
        r"\small",
        rf"\caption{{{latex_escape(caption)}}}",
        rf"\label{{{label}}}",
        r"\resizebox{\linewidth}{!}{%",
        rf"\begin{{tabular}}{{{align}}}",
        r"\toprule",
        " & ".join(latex_escape(col) for col in cols) + r" \\",
        r"\midrule",
    ]
    for _, row in df.iterrows():
        lines.append(" & ".join(latex_escape(row[col]) for col in cols) + r" \\")
    lines.extend([r"\bottomrule", r"\end{tabular}%", r"}"])
    if notes:
        lines.append(rf"\vspace{{2mm}}\parbox{{0.95\linewidth}}{{\footnotesize {latex_escape(notes)}}}")
    lines.append(r"\end{table}")
    path.write_text("\n".join(lines) + "\n")


def make_tables(primary: pd.DataFrame, model_df: pd.DataFrame, wgs: pd.DataFrame, enrichment: pd.DataFrame) -> None:
    def manuscript_interpretation(row: pd.Series) -> str:
        pair = str(row["pair"])
        site = str(row["site"])
        centered, _, _ = parse_auc_ci(row["stratum_centered_auc_95ci"])
        if "caution" in str(row["adequacy"]):
            return "Caution; low matched support"
        if pair == "E. coli / Cipro" and centered >= 0.59:
            return "Retained residual signal"
        if pair == "E. coli / Amox-Clav" and site in {"DRIAMS-C", "DRIAMS-D"}:
            return "Near chance after matching"
        if pair == "E. coli / Amox-Clav":
            return "Weak or uncertain residual signal"
        return str(row.get("interpretation", ""))

    t1 = primary[primary["pair"].isin(["E. coli / Cipro", "E. coli / Amox-Clav"])][
        ["pair", "site", "raw_auc_95ci", "stratum_centered_auc_95ci", "matched_retention_pct", "n_matched", "permutation_p", "adequacy", "interpretation"]
    ].copy()
    t1["interpretation"] = t1.apply(manuscript_interpretation, axis=1)
    t1.columns = ["Pair", "Site", "Raw AUC", "Centered AUC", "Retention (%)", "n matched", "P perm", "Adequacy", "Interpretation"]
    t1["Retention (%)"] = t1["Retention (%)"].map(lambda x: f"{float(x):.1f}")
    t1["P perm"] = t1["P perm"].map(lambda x: "" if pd.isna(x) else f"{float(x):.3f}")
    write_latex_table(
        TABLE_DIR / "table_1_primary_audit.tex",
        "Primary background-matched audit for the E. coli ciprofloxacin and amoxicillin-clavulanic acid contrast.",
        "tab:primary-audit",
        t1,
        "Centered AUC is computed after subtracting the mean model score within each co-resistance background stratum. Caution rows are reported but not used as the main evidence.",
    )

    t2 = model_df[model_df["drug"].isin(["Cipro", "Amox-Clav"]) & model_df["site"].isin(["A-2018", "DRIAMS-C", "DRIAMS-D"])][
        ["site", "drug", "cnn_raw_auc", "cnn_centered_auc", "lgbm_raw_auc", "lgbm_centered_auc", "model_family_consensus"]
    ].copy()
    t2.columns = ["Site", "Drug", "CNN raw", "CNN centered", "LGBM raw", "LGBM centered", "Consensus"]
    for col in ["CNN raw", "CNN centered", "LGBM raw", "LGBM centered"]:
        t2[col] = t2[col].map(lambda x: "" if pd.isna(x) else f"{float(x):.3f}")
    write_latex_table(
        TABLE_DIR / "table_2_model_replication.tex",
        "Model-family replication of raw-to-background-centered attenuation.",
        "tab:model-replication",
        t2,
        "The same audit is applied to CNN and multi-task LGBM predictions.",
    )

    edges = pd.read_csv(FINAL / "table_5_top_cross_resistance_edges.csv").head(6).copy()
    edges = edges[["drug_a", "drug_b", "n_both_known", "n_rr", "rr_lift", "phi", "resistant_jaccard"]]
    edges.columns = ["Drug A", "Drug B", "n both", "n RR", "RR lift", "Phi", "R Jaccard"]
    for col in ["RR lift", "Phi", "R Jaccard"]:
        edges[col] = edges[col].map(lambda x: f"{float(x):.3f}")
    write_latex_table(
        TABLE_DIR / "table_3_cross_resistance_edges.tex",
        "Strongest co-resistance edges in the E. coli expanded panel.",
        "tab:cross-resistance",
        edges,
        "Edges are computed across isolate-level AST labels. RR lift compares observed double resistance with the expectation under independence.",
    )

    t4 = wgs.copy()
    t4["target"] = t4["target"].str.replace("_", " ", regex=False)
    t4 = t4[["target", "n", "class_1", "auc", "folds", "model"]]
    t4.columns = ["Target", "n", "positive", "AUC", "folds", "Model"]
    t4["AUC"] = t4["AUC"].map(lambda x: f"{float(x):.3f}")
    write_latex_table(
        TABLE_DIR / "table_4_public_wgs_auc.tex",
        "Public WGS-linked Basel UPEC MALDI validation.",
        "tab:public-wgs",
        t4,
        "A simple centroid-direction classifier is used to ask whether lineage and resistance labels are encoded in Bruker MALDI peak features.",
    )

    t5 = enrichment[enrichment["target"].isin(["ST131", "Ciprofloxacin_R", "Ceftriaxone_R", "ALL_TARGETS"])].copy()
    t5["target"] = t5["target"].str.replace("_", " ", regex=False)
    t5 = t5[["target", "observed_overlap_count", "top_peak_count", "null_mean_overlap_count", "fold_enrichment", "empirical_p_ge_observed"]]
    t5.columns = ["Target", "Observed", "Top peaks", "Null mean", "Fold enrichment", "P perm"]
    for col in ["Null mean", "Fold enrichment"]:
        t5[col] = t5[col].map(lambda x: f"{float(x):.2f}")
    t5["P perm"] = t5["P perm"].map(lambda x: f"{float(x):.4f}")
    write_latex_table(
        TABLE_DIR / "table_5_biomarker_enrichment.tex",
        "Mass-matched enrichment of published ST131 MALDI biomarkers among discriminative public UPEC peaks.",
        "tab:biomarker-enrichment",
        t5,
        "The null model preserves the coarse m/z-stratum distribution of selected peaks.",
    )

    guardrails = pd.DataFrame(
        [
            ["Supported", "MALDI-AMR prediction is background-sensitive, and raw AUC can overstate focal-drug signal."],
            ["Supported", "Background-centered AUC separates retained within-background ranking from background-driven attenuation."],
            ["Supported", "Public WGS-linked Bruker data show that ST131 lineage is strongly encoded in MALDI peak features."],
            ["Not claimed", "The DRIAMS CNN is definitively detecting ST131."],
            ["Not claimed", "The discriminative DRIAMS or UPEC peaks have definitive protein identities."],
            ["Not claimed", "All MALDI-AMR performance is clonal confounding."],
        ],
        columns=["Status", "Statement"],
    )
    write_latex_table(
        TABLE_DIR / "table_6_claim_guardrails.tex",
        "Claim guardrails used throughout the manuscript.",
        "tab:claim-guardrails",
        guardrails,
        "These guardrails distinguish the evaluation-framework claim from direct clone or protein-identification claims.",
    )

    comparison = pd.DataFrame(
        [(dim, prior, tripod, ours) for dim, prior, tripod, ours, _ in COMPARISON_ROWS],
        columns=["Dimension", "Prior MALDI-AMR", "TRIPOD+AI", "This framework"],
    )
    write_latex_table(
        TABLE_DIR / "table_7_framework_comparison.tex",
        "What the background-matched audit adds relative to prior work and TRIPOD+AI.",
        "tab:framework-comparison",
        comparison,
        r"Prior MALDI-AMR refers to Weis et al.\ 2022 (\textit{Nature Medicine}) and comparable publications. "
        r"TRIPOD+AI refers to the transparent reporting guideline (2024). "
        r"N/A denotes dimensions not applicable to a reporting standard.",
    )


def write_source_data(primary: pd.DataFrame, model_df: pd.DataFrame, wgs: pd.DataFrame,
                      enrichment: pd.DataFrame,
                      matrix_df: pd.DataFrame | None = None,
                      sa_summary: pd.DataFrame | None = None,
                      ecoli_bg: pd.DataFrame | None = None,
                      saureus_bg: pd.DataFrame | None = None) -> None:
    """Write figure source-data CSVs in a Nature-style layout."""
    primary.copy().to_csv(SOURCE_DIR / "source_data_fig2_primary_background_audit.csv", index=False)

    if matrix_df is not None and not matrix_df.empty:
        fig3 = figure_3_matrix_rows(matrix_df)
        fig3 = fig3[
            [
                "display_target",
                "model_column",
                "model_class",
                "model_variant",
                "scope",
                "organism",
                "drug",
                "site",
                "raw_auc",
                "centered_auc",
                "pairwise_accuracy",
                "matched_retention",
                "valid_strata",
                "adequacy_label",
                "interpretation",
                "source_path",
            ]
        ].copy()
    else:
        fig3 = model_df[
            model_df["drug"].isin(["Cipro", "Amox-Clav"])
            & model_df["site"].isin(["A-2018", "DRIAMS-C", "DRIAMS-D"])
        ].copy()
    fig3.to_csv(SOURCE_DIR / "source_data_fig3_model_family_replication.csv", index=False)

    edges = pd.read_csv(ANALYSIS / "cross_resistance_network" / "cross_resistance_edges.csv")
    edges[edges["site"].eq("ALL")].to_csv(SOURCE_DIR / "source_data_fig4_cross_resistance_edges_all_sites.csv", index=False)

    wgs.copy().to_csv(SOURCE_DIR / "source_data_fig5a_public_wgs_maldi_auc.csv", index=False)
    enrichment.copy().to_csv(SOURCE_DIR / "source_data_fig5b_proteomic_biomarker_enrichment.csv", index=False)

    manifest = [
        ["Figure 2", "source_data_fig2_primary_background_audit.csv", "Raw, matched and background-centered AUC values for the primary E. coli contrast."],
        ["Figure 3", "source_data_fig3_model_family_replication.csv", "Completed CNN/Mega, LightGBM multi-task and LightGBM single-task raw-to-centered model-class audit rows."],
        ["Figure 4", "source_data_fig4_cross_resistance_edges_all_sites.csv", "All-sites E. coli cross-resistance edges used for the phi heatmap and strongest-edge annotations."],
        ["Figure 5a", "source_data_fig5a_public_wgs_maldi_auc.csv", "Public UPEC WGS-linked MALDI peak-feature AUCs."],
        ["Figure 5b", "source_data_fig5b_proteomic_biomarker_enrichment.csv", "Published ST131 biomarker enrichment results."],
    ]
    if sa_summary is not None and ecoli_bg is not None and saureus_bg is not None:
        fig6_df = _build_three_way_data(_prep_auc_ci(primary), sa_summary, ecoli_bg, saureus_bg)
        fig6_df.to_csv(SOURCE_DIR / "source_data_fig6_three_way_decomposition.csv", index=False)
        manifest.append(["Figure 6", "source_data_fig6_three_way_decomposition.csv",
                         "Three-way AUC decomposition: raw MALDI, co-resistance-only, background-centred, for E. coli Cipro/Amox-Clav and S. aureus/Oxacillin."])

    pd.DataFrame(manifest, columns=["display_item", "file", "description"]).to_csv(SOURCE_DIR / "source_data_manifest.csv", index=False)


def main() -> None:
    ensure_dirs()
    primary    = pd.read_csv(FINAL / "table_1_primary_background_matched_audit.csv")
    model_df   = pd.read_csv(FINAL / "table_2_cnn_vs_lgbm_multi_background_audit.csv")
    wgs        = pd.read_csv(FINAL / "table_7_public_wgs_maldi_auc.csv")
    enrichment = pd.read_csv(FINAL / "table_9_published_st131_biomarker_enrichment.csv")
    matrix_path = ANALYSIS / "model_class_matrix" / "model_class_matrix.csv"
    matrix_df = pd.read_csv(matrix_path) if matrix_path.exists() else pd.DataFrame()

    sa_summary_path = ANALYSIS / "saureus_panel_oxa_background_audit" / "background_matched_audit_summary.csv"
    ecoli_bg_path   = ANALYSIS / "co_resistance_only_baseline_ecoli"   / "co_resistance_only_baseline.csv"
    saureus_bg_path = ANALYSIS / "co_resistance_only_baseline_saureus" / "co_resistance_only_baseline.csv"

    sa_summary  = pd.read_csv(sa_summary_path)  if sa_summary_path.exists()  else None
    ecoli_bg    = pd.read_csv(ecoli_bg_path)    if ecoli_bg_path.exists()    else None
    saureus_bg  = pd.read_csv(saureus_bg_path)  if saureus_bg_path.exists()  else None

    created = [
        figure_1_framework(),
        figure_2_primary_audit(primary),
        figure_3_model_replication(model_df, matrix_df),
        figure_4_cross_resistance(),
        figure_5_public_support(wgs, enrichment),
        figure_7_framework_comparison(),
    ]
    if sa_summary is not None and ecoli_bg is not None and saureus_bg is not None:
        created.append(figure_6_three_way_decomposition(primary, sa_summary, ecoli_bg, saureus_bg))
    else:
        print("  [figure_6] one or more inputs missing — check saureus_panel_oxa_background_audit, "
              "co_resistance_only_baseline_ecoli, co_resistance_only_baseline_saureus.")

    make_tables(primary, model_df, wgs, enrichment)
    write_source_data(primary, model_df, wgs, enrichment, matrix_df, sa_summary, ecoli_bg, saureus_bg)

    print("Created figures:")
    for path in created:
        print(f"  {path.relative_to(ROOT)}")
    print("Created LaTeX tables:")
    for path in sorted(TABLE_DIR.glob("*.tex")):
        print(f"  {path.relative_to(ROOT)}")
    print("Created source-data files:")
    for path in sorted(SOURCE_DIR.glob("*.csv")):
        print(f"  {path.relative_to(ROOT)}")


if __name__ == "__main__":
    main()

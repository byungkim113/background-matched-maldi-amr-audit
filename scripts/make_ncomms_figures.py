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
    "E. coli / Cipro": "E. coli / Cipro",
    "E. coli / Amox-Clav": "E. coli / Amox-Clav",
}

DRUG_ORDER = ["Cipro", "Norflox", "Amox-Clav", "CRO", "CAZ", "FEP"]


def ensure_dirs() -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    TABLE_DIR.mkdir(parents=True, exist_ok=True)


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


def marker(c: canvas.Canvas, x: float, y: float, color, shape: str, hollow: bool = False, size: float = 5) -> None:
    c.setStrokeColor(color)
    c.setFillColor(colors.white if hollow else color)
    c.setLineWidth(1.2)
    if shape == "circle":
        c.circle(x, y, size, stroke=1, fill=1)
    elif shape == "square":
        c.rect(x - size, y - size, 2 * size, 2 * size, stroke=1, fill=1)
    else:
        c.circle(x, y, size, stroke=1, fill=1)


def figure_1_framework() -> Path:
    path = FIG_DIR / "figure_1_framework.pdf"
    c = canvas.Canvas(str(path), pagesize=landscape(letter))
    w, h = landscape(letter)
    text(c, 0.55 * inch, h - 0.45 * inch, "Fig. 1 | Background-matched MALDI-AMR audit", 14, True)
    text(c, 0.55 * inch, h - 0.70 * inch, "A model-agnostic test of whether focal-drug prediction survives after controlling for co-resistance background.", 9, False, GRAY)

    y = h - 1.65 * inch
    boxes = [
        ("Model predictions", "isolate ID, site/year,\norganism, drug, label,\nprobability score"),
        ("Raw transfer", "External AUC/AUPR\nbefore background control"),
        ("Background match", "Group isolates by labels\nof the other antibiotics"),
        ("Centered audit", "Subtract each stratum's\nmean model score"),
        ("Interpretation", "retained focal signal,\npartial retention, or\nbackground-driven collapse"),
    ]
    x0 = 0.55 * inch
    box_w = 1.65 * inch
    gap = 0.33 * inch
    box_h = 1.28 * inch
    for i, (title, body) in enumerate(boxes):
        x = x0 + i * (box_w + gap)
        c.setFillColor(PALE)
        c.setStrokeColor(LIGHT_GRAY)
        c.roundRect(x, y, box_w, box_h, 8, stroke=1, fill=1)
        text(c, x + 0.12 * inch, y + box_h - 0.28 * inch, title, 9, True, DARK)
        for j, line in enumerate(body.split("\n")):
            text(c, x + 0.12 * inch, y + box_h - 0.52 * inch - j * 0.18 * inch, line, 7.5, False, GRAY)
        if i < len(boxes) - 1:
            c.setStrokeColor(MID_GRAY)
            c.setLineWidth(1.0)
            ax = x + box_w + 0.05 * inch
            ay = y + box_h / 2
            c.line(ax, ay, ax + gap - 0.10 * inch, ay)
            c.line(ax + gap - 0.15 * inch, ay + 4, ax + gap - 0.10 * inch, ay)
            c.line(ax + gap - 0.15 * inch, ay - 4, ax + gap - 0.10 * inch, ay)

    panel_label(c, 0.60 * inch, 2.25 * inch, "Audit outputs")
    outputs = [
        ("raw AUC", "apparent external performance"),
        ("matched AUC", "performance in strata with both R and S isolates"),
        ("background-centered AUC", "ranking after stratum-level score shifts are removed"),
        ("matched retention", "how much of the target comparison remains interpretable"),
        ("cross-resistance network", "the label ecology the model could exploit"),
    ]
    for i, (name, desc) in enumerate(outputs):
        yy = 1.95 * inch - i * 0.25 * inch
        c.setFillColor(BLUE if i == 2 else LIGHT_GRAY)
        c.circle(0.75 * inch, yy + 3, 3, stroke=0, fill=1)
        text(c, 0.90 * inch, yy, name, 8.5, True, DARK)
        text(c, 2.35 * inch, yy, desc, 8.5, False, GRAY)

    c.showPage()
    c.save()
    return path


def figure_2_primary_audit(primary: pd.DataFrame) -> Path:
    path = FIG_DIR / "figure_2_primary_background_audit.pdf"
    c = canvas.Canvas(str(path), pagesize=landscape(letter))
    w, h = landscape(letter)
    text(c, 0.55 * inch, h - 0.45 * inch, "Fig. 2 | Focal-drug signal after background matching", 14, True)
    text(c, 0.55 * inch, h - 0.70 * inch, "Raw external AUC is compared with stratum-centered AUC. Centering removes score shifts shared by co-resistance strata.", 9, False, GRAY)

    data = primary.copy()
    for col in ["raw_auc_95ci", "stratum_centered_auc_95ci"]:
        values = data[col].map(parse_auc_ci)
        data[col + "_mid"] = values.map(lambda t: t[0])
        data[col + "_low"] = values.map(lambda t: t[1])
        data[col + "_high"] = values.map(lambda t: t[2])

    panels = [
        ("E. coli / Cipro", BLUE, 0.55 * inch, "Ciprofloxacin retains residual within-background ranking."),
        ("E. coli / Amox-Clav", ORANGE, 5.65 * inch, "Amox-Clav weakens or collapses after background control."),
    ]
    xmin, xmax = 0.45, 0.90
    for pair, col, px, subtitle in panels:
        panel_label(c, px, h - 1.15 * inch, pair)
        text(c, px, h - 1.38 * inch, subtitle, 8, False, GRAY)
        x0, x1 = px + 1.15 * inch, px + 4.55 * inch
        y0 = 1.55 * inch
        axis_auc(c, x0, x1, y0, xmin, xmax)
        sub = data[data["pair"].eq(pair)].copy()
        sub["site_order"] = sub["site"].map({"A-2018": 0, "DRIAMS-B": 1, "DRIAMS-C": 2, "DRIAMS-D": 3})
        sub = sub.sort_values("site_order")
        for i, row in enumerate(sub.itertuples(index=False)):
            y = h - 1.88 * inch - i * 0.47 * inch
            site = row.site
            caution = "caution" in str(row.adequacy)
            right_text(c, x0 - 0.18 * inch, y - 3, site, 8, False, GRAY if caution else DARK)
            raw = getattr(row, "raw_auc_95ci_mid")
            cen = getattr(row, "stratum_centered_auc_95ci_mid")
            raw_x = xmap(raw, x0, x1, xmin, xmax)
            cen_x = xmap(cen, x0, x1, xmin, xmax)
            c.setStrokeColor(col if not caution else MID_GRAY)
            c.setLineWidth(1.2)
            c.line(raw_x, y, cen_x, y)
            marker(c, raw_x, y, col if not caution else MID_GRAY, "circle", hollow=False)
            marker(c, cen_x, y, col if not caution else MID_GRAY, "square", hollow=caution)
            delta = getattr(row, "raw_to_centered_delta")
            text(c, x1 + 0.12 * inch, y - 3, f"Delta={delta:+.3f}; ret={row.matched_retention_pct:.1f}%", 7.2, False, GRAY)
        marker(c, px + 0.10 * inch, y0 - 0.18 * inch, DARK, "circle")
        text(c, px + 0.25 * inch, y0 - 0.22 * inch, "raw", 7.5, False, GRAY)
        marker(c, px + 0.75 * inch, y0 - 0.18 * inch, DARK, "square")
        text(c, px + 0.90 * inch, y0 - 0.22 * inch, "background-centered", 7.5, False, GRAY)
    c.showPage()
    c.save()
    return path


def figure_3_model_replication(model_df: pd.DataFrame) -> Path:
    path = FIG_DIR / "figure_3_model_family_replication.pdf"
    c = canvas.Canvas(str(path), pagesize=landscape(letter))
    w, h = landscape(letter)
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


def figure_5_public_support(wgs: pd.DataFrame, enrichment: pd.DataFrame) -> Path:
    path = FIG_DIR / "figure_5_public_wgs_proteomic_support.pdf"
    c = canvas.Canvas(str(path), pagesize=landscape(letter))
    w, h = landscape(letter)
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
        rf"\begin{{tabular}}{{{align}}}",
        r"\toprule",
        " & ".join(latex_escape(col) for col in cols) + r" \\",
        r"\midrule",
    ]
    for _, row in df.iterrows():
        lines.append(" & ".join(latex_escape(row[col]) for col in cols) + r" \\")
    lines.extend([r"\bottomrule", r"\end{tabular}"])
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


def main() -> None:
    ensure_dirs()
    primary = pd.read_csv(FINAL / "table_1_primary_background_matched_audit.csv")
    model_df = pd.read_csv(FINAL / "table_2_cnn_vs_lgbm_multi_background_audit.csv")
    wgs = pd.read_csv(FINAL / "table_7_public_wgs_maldi_auc.csv")
    enrichment = pd.read_csv(FINAL / "table_9_published_st131_biomarker_enrichment.csv")

    created = [
        figure_1_framework(),
        figure_2_primary_audit(primary),
        figure_3_model_replication(model_df),
        figure_4_cross_resistance(),
        figure_5_public_support(wgs, enrichment),
    ]
    make_tables(primary, model_df, wgs, enrichment)

    print("Created figures:")
    for path in created:
        print(f"  {path.relative_to(ROOT)}")
    print("Created LaTeX tables:")
    for path in sorted(TABLE_DIR.glob("*.tex")):
        print(f"  {path.relative_to(ROOT)}")


if __name__ == "__main__":
    main()

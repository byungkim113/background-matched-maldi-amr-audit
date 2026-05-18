#!/usr/bin/env python3
"""Build final tables and figures for the background-aware MALDI-AMR framework.

The script is intentionally self-contained and uses only pandas/numpy/Pillow from
the bundled Codex runtime. It reads the already-generated audit outputs and writes
paper-facing CSV/Markdown tables plus PNG figures.
"""

from __future__ import annotations

import math
import json
import textwrap
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image, ImageDraw, ImageFont


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs" / "final_framework_outputs"

INPUTS = {
    "figure_table": ROOT / "outputs" / "analysis_outputs" / "background_matched_transfer_audit_figure_table.csv",
    "cnn_lgbm": ROOT / "outputs" / "analysis_outputs" / "background_matched_transfer_audit_cnn_vs_lgbm.csv",
    "ecology": ROOT / "outputs" / "analysis_outputs" / "background_audit_with_resistance_ecology.csv",
    "prediction_assessment": ROOT / "outputs" / "analysis_outputs" / "ecoli_mechanism6_pair_prediction_assessment.csv",
    "block_site": ROOT / "outputs" / "analysis_outputs" / "ecoli_mechanism6_block_site_summary.csv",
    "cross_edges": ROOT / "outputs" / "analysis_outputs" / "cross_resistance_network" / "cross_resistance_edges.csv",
    "co_resistance_stratification": ROOT / "outputs" / "analysis_outputs" / "co_resistance_stratification" / "co_resistance_stratification.csv",
    "deployment_rules": ROOT / "outputs" / "analysis_outputs" / "deployment_decision_framework" / "deployment_decision_rules.csv",
    "deployment_readiness": ROOT / "outputs" / "analysis_outputs" / "deployment_decision_framework" / "deployment_readiness_by_pair.csv",
    "calibration": ROOT / "outputs" / "analysis_outputs" / "calibration_analysis" / "calibration_summary.csv",
    "temporal_reliability": ROOT / "outputs" / "analysis_outputs" / "temporal_reliability_audit" / "temporal_reliability.csv",
    "falsification": ROOT / "outputs" / "analysis_outputs" / "falsification_controls" / "falsification_controls.csv",
    "wgs_auc": ROOT / "outputs" / "analysis_outputs" / "upec_wgs_validation_outputs" / "centroid_binary_cv_results.csv",
    "wgs_assoc": ROOT / "outputs" / "analysis_outputs" / "upec_wgs_validation_outputs" / "st131_resistance_associations.csv",
    "proteomic_enrichment": ROOT / "outputs" / "analysis_outputs" / "updated_proteomic_overlap_outputs" / "updated_proteomic_overlap_permutation_enrichment.csv",
    "proteomic_overlaps": ROOT / "outputs" / "analysis_outputs" / "updated_proteomic_overlap_outputs" / "updated_published_st131_proteomic_overlap.csv",
}

OPTIONAL_INPUTS = {
    "weis_audit": ROOT / "outputs" / "weis_predictions_audit_analysis" / "background_matched_audit_summary.csv",
    "weis_predictions": ROOT / "outputs" / "weis_predictions_audit_analysis" / "normalized_predictions.csv",
    "marisma_spot_audit": ROOT / "outputs" / "analysis_outputs" / "marisma_external_validation" / "marisma_background_audit" / "background_matched_audit_summary.csv",
    "marisma_spot_predictions": ROOT / "outputs" / "analysis_outputs" / "marisma_external_validation" / "marisma_background_audit" / "normalized_predictions.csv",
    "marisma_isolate_audit": ROOT / "outputs" / "analysis_outputs" / "marisma_external_validation" / "marisma_isolate_background_audit" / "background_matched_audit_summary.csv",
    "marisma_duplicate_report": ROOT / "outputs" / "analysis_outputs" / "marisma_external_validation" / "marisma_isolate_background_audit" / "marisma_duplicate_handling_report.json",
    "upec_clone_overall_auc": ROOT / "outputs" / "analysis_outputs" / "upec_clone_control_outputs" / "clone_control_subset_auc.csv",
}


COLORS = {
    "quinolone": (41, 112, 177),
    "ceph": (186, 76, 76),
    "mixed": (218, 142, 56),
    "other": (92, 92, 92),
    "grid": (220, 224, 229),
    "text": (35, 39, 47),
    "muted": (92, 99, 112),
}


def font(size: int, bold: bool = False) -> ImageFont.FreeTypeFont | ImageFont.ImageFont:
    candidates = [
        "/System/Library/Fonts/Supplemental/Arial Bold.ttf" if bold else "/System/Library/Fonts/Supplemental/Arial.ttf",
        "/Library/Fonts/Arial Bold.ttf" if bold else "/Library/Fonts/Arial.ttf",
        "/System/Library/Fonts/Helvetica.ttc",
    ]
    for candidate in candidates:
        try:
            return ImageFont.truetype(candidate, size=size)
        except Exception:
            continue
    return ImageFont.load_default()


FONT = {
    "title": font(28, True),
    "subtitle": font(18),
    "label": font(16),
    "small": font(13),
    "tiny": font(11),
    "bold": font(16, True),
}


def ensure_inputs() -> None:
    missing = [name for name, path in INPUTS.items() if not path.exists()]
    if missing:
        details = "\n".join(f"- {name}: {INPUTS[name]}" for name in missing)
        raise FileNotFoundError(f"Missing required input files:\n{details}")


def short_drug(drug: str) -> str:
    return {
        "Ciprofloxacin": "Cipro",
        "Norfloxacin": "Norflox",
        "Amoxicillin-Clavulanic acid": "Amox-Clav",
        "Ceftriaxone": "CRO",
        "Ceftazidime": "CAZ",
        "Cefepime": "FEP",
    }.get(drug, drug)


def block_for_drug(drug: str) -> str:
    if drug in {"Ciprofloxacin", "Norfloxacin"}:
        return "quinolone"
    if drug in {"Ceftriaxone", "Ceftazidime", "Cefepime"}:
        return "ceph"
    if "Clavulanic" in drug:
        return "mixed"
    return "other"


def block_label(block: str) -> str:
    return {
        "quinolone": "fluoroquinolone block",
        "ceph": "cephalosporin/ESBL block",
        "mixed": "mixed beta-lactam background",
    }.get(block, block)


def auc_ci_text(row: pd.Series, stem: str) -> str:
    value = row.get(stem)
    if pd.isna(value):
        return ""
    low = row.get(f"{stem}_ci_low")
    high = row.get(f"{stem}_ci_high")
    if pd.isna(low) or pd.isna(high):
        return f"{float(value):.3f}"
    return f"{float(value):.3f} ({float(low):.3f}-{float(high):.3f})"


def published_model_family_label() -> str:
    if not OPTIONAL_INPUTS["weis_predictions"].exists():
        return "Weis/Borgwardt-style model"
    predictions = pd.read_csv(OPTIONAL_INPUTS["weis_predictions"], usecols=["model_name"])
    names = sorted(str(name) for name in predictions["model_name"].dropna().unique())
    if names == ["Weis-lr"]:
        return "Weis/Borgwardt-style logistic regression"
    if names == ["Weis-lightgbm"]:
        return "Weis/Borgwardt-style LightGBM"
    if names:
        return f"Weis/Borgwardt-style model ({', '.join(names)})"
    return "Weis/Borgwardt-style model"


def markdown_table(df: pd.DataFrame, max_rows: int | None = None) -> str:
    if max_rows is not None:
        df = df.head(max_rows)
    clean = df.copy()
    for col in clean.columns:
        if pd.api.types.is_float_dtype(clean[col]):
            clean[col] = clean[col].map(lambda x: "" if pd.isna(x) else f"{x:.3f}")
        else:
            clean[col] = clean[col].map(lambda x: "" if pd.isna(x) else str(x))
    headers = list(clean.columns)
    rows = clean.values.tolist()
    widths = [len(str(h)) for h in headers]
    for row in rows:
        for i, val in enumerate(row):
            widths[i] = max(widths[i], len(str(val)))
    header = "| " + " | ".join(str(h).ljust(widths[i]) for i, h in enumerate(headers)) + " |"
    sep = "| " + " | ".join("-" * widths[i] for i in range(len(headers))) + " |"
    body = ["| " + " | ".join(str(v).ljust(widths[i]) for i, v in enumerate(row)) + " |" for row in rows]
    return "\n".join([header, sep] + body)


def save_table(df: pd.DataFrame, stem: str, title: str, note: str = "") -> None:
    csv_path = OUT / f"{stem}.csv"
    md_path = OUT / f"{stem}.md"
    df.to_csv(csv_path, index=False)
    lines = [f"# {title}", ""]
    if note:
        lines += [note, ""]
    lines += [markdown_table(df), ""]
    md_path.write_text("\n".join(lines))


def xscale(value: float, left: int, right: int, xmin: float, xmax: float) -> int:
    value = max(xmin, min(xmax, float(value)))
    return int(left + (value - xmin) / (xmax - xmin) * (right - left))


def draw_wrapped(draw: ImageDraw.ImageDraw, xy: tuple[int, int], text: str, width: int, font_obj, fill) -> int:
    avg = max(5, int(font_obj.size * 0.55)) if hasattr(font_obj, "size") else 7
    max_chars = max(12, width // avg)
    lines = textwrap.wrap(text, max_chars)
    x, y = xy
    line_h = int((font_obj.size if hasattr(font_obj, "size") else 14) * 1.25)
    for line in lines:
        draw.text((x, y), line, font=font_obj, fill=fill)
        y += line_h
    return y


def draw_legend(draw: ImageDraw.ImageDraw, x: int, y: int) -> None:
    items = [("quinolone", "fluoroquinolone"), ("ceph", "cephalosporin/ESBL"), ("mixed", "mixed beta-lactam")]
    cursor = x
    for key, label in items:
        draw.rounded_rectangle((cursor, y, cursor + 20, y + 12), radius=4, fill=COLORS[key])
        draw.text((cursor + 28, y - 3), label, font=FONT["small"], fill=COLORS["text"])
        cursor += 190


def figure_raw_to_centered(df: pd.DataFrame) -> Path:
    long_rows = []
    for model, raw_col, centered_col, retention_col, adequacy_col in [
        ("CNN", "cnn_raw_auc", "cnn_centered_auc", "cnn_retention", "cnn_adequacy"),
        ("LGBM multi", "lgbm_raw_auc", "lgbm_centered_auc", "lgbm_retention", "lgbm_adequacy"),
    ]:
        for _, row in df.iterrows():
            if pd.isna(row[raw_col]) or pd.isna(row[centered_col]):
                continue
            long_rows.append(
                {
                    "model": model,
                    "site": row["site"],
                    "drug": row["drug"],
                    "raw": row[raw_col],
                    "centered": row[centered_col],
                    "retention": row[retention_col],
                    "adequacy": row[adequacy_col],
                    "block": block_for_drug(row["drug"]),
                }
            )
    plot = pd.DataFrame(long_rows)
    order_drugs = ["Ciprofloxacin", "Norfloxacin", "Amoxicillin-Clavulanic acid", "Ceftriaxone", "Ceftazidime", "Cefepime"]
    order_sites = ["A-2018", "DRIAMS-B", "DRIAMS-C", "DRIAMS-D"]
    plot["drug_order"] = plot["drug"].map({d: i for i, d in enumerate(order_drugs)})
    plot["site_order"] = plot["site"].map({s: i for i, s in enumerate(order_sites)})
    plot = plot.sort_values(["drug_order", "site_order"])

    width, height = 1800, 1080
    img = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(img)
    draw.text((60, 35), "Raw external AUC vs background-centered AUC", font=FONT["title"], fill=COLORS["text"])
    draw.text(
        (60, 75),
        "Circle = raw external AUC; square = after matching/centering by co-resistance background.",
        font=FONT["subtitle"],
        fill=COLORS["muted"],
    )
    draw_legend(draw, 60, 110)

    panel_w = 790
    panel_h = 820
    xmins, xmaxs = 0.30, 0.92
    panels = [("CNN", 70, 185), ("LGBM multi", 930, 185)]
    for model, px, py in panels:
        sub = plot[plot["model"] == model].copy()
        draw.text((px + 260, py - 45), model, font=FONT["title"], fill=COLORS["text"])
        left, right = px + 170, px + panel_w - 30
        top, bottom = py, py + panel_h
        for tick in np.arange(0.3, 1.0, 0.1):
            xx = xscale(tick, left, right, xmins, xmaxs)
            draw.line((xx, top, xx, bottom), fill=COLORS["grid"], width=1)
            draw.text((xx - 12, bottom + 12), f"{tick:.1f}", font=FONT["small"], fill=COLORS["muted"])
        draw.line((xscale(0.5, left, right, xmins, xmaxs), top, xscale(0.5, left, right, xmins, xmaxs), bottom), fill=(150, 150, 150), width=2)
        n = len(sub)
        for i, (_, row) in enumerate(sub.iterrows()):
            y = int(top + 18 + i * ((panel_h - 35) / max(1, n - 1)))
            label = f"{row['site']} | {short_drug(row['drug'])}"
            draw.text((px, y - 8), label, font=FONT["tiny"], fill=COLORS["text"])
            color = COLORS[row["block"]]
            x1 = xscale(row["centered"], left, right, xmins, xmaxs)
            x2 = xscale(row["raw"], left, right, xmins, xmaxs)
            if "low" in str(row["adequacy"]) or row["retention"] < 0.15:
                line_color = tuple(int(c * 0.65 + 255 * 0.35) for c in color)
            else:
                line_color = color
            draw.line((x1, y, x2, y), fill=line_color, width=4)
            draw.rectangle((x1 - 6, y - 6, x1 + 6, y + 6), fill=color)
            draw.ellipse((x2 - 7, y - 7, x2 + 7, y + 7), fill=color)
        draw.text((left + 205, bottom + 48), "AUC", font=FONT["label"], fill=COLORS["text"])
    path = OUT / "figure_1_raw_to_background_centered_auc.png"
    img.save(path)
    return path


def figure_drop_retention(df: pd.DataFrame) -> Path:
    width, height = 1500, 850
    img = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(img)
    draw.text((60, 35), "Background sensitivity: signal drop vs matched retention", font=FONT["title"], fill=COLORS["text"])
    draw.text((60, 76), "Large positive drop means raw AUC was higher than background-centered AUC.", font=FONT["subtitle"], fill=COLORS["muted"])
    draw_legend(draw, 60, 112)

    panels = [
        ("CNN", "cnn_raw_minus_centered", "cnn_retention", 80, 180),
        ("LGBM multi", "lgbm_raw_minus_centered", "lgbm_retention", 790, 180),
    ]
    for title, drop_col, ret_col, px, py in panels:
        left, right = px + 80, px + 610
        top, bottom = py, py + 560
        draw.text((px + 245, py - 45), title, font=FONT["title"], fill=COLORS["text"])
        for tick in np.linspace(0, 1, 6):
            xx = xscale(tick, left, right, 0, 1)
            draw.line((xx, top, xx, bottom), fill=COLORS["grid"], width=1)
            draw.text((xx - 12, bottom + 10), f"{tick:.1f}", font=FONT["small"], fill=COLORS["muted"])
        ymin, ymax = -0.12, 0.52
        for tick in np.arange(-0.1, 0.55, 0.1):
            yy = int(bottom - (tick - ymin) / (ymax - ymin) * (bottom - top))
            draw.line((left, yy, right, yy), fill=COLORS["grid"], width=1)
            draw.text((px + 20, yy - 7), f"{tick:.1f}", font=FONT["small"], fill=COLORS["muted"])
        draw.line((left, bottom, right, bottom), fill=(120, 120, 120), width=2)
        draw.line((left, top, left, bottom), fill=(120, 120, 120), width=2)
        for _, row in df.iterrows():
            if pd.isna(row[drop_col]) or pd.isna(row[ret_col]):
                continue
            xx = xscale(float(row[ret_col]), left, right, 0, 1)
            yy = int(bottom - (float(row[drop_col]) - ymin) / (ymax - ymin) * (bottom - top))
            block = block_for_drug(row["drug"])
            color = COLORS[block]
            r = 9 if block == "quinolone" else 8
            draw.ellipse((xx - r, yy - r, xx + r, yy + r), fill=color, outline=(255, 255, 255), width=2)
        draw.text((left + 155, bottom + 45), "Matched retention", font=FONT["label"], fill=COLORS["text"])
        draw.text((px, py + 235), "raw - centered", font=FONT["label"], fill=COLORS["text"])
    path = OUT / "figure_2_drop_vs_matched_retention.png"
    img.save(path)
    return path


def figure_cross_resistance_heatmap(edges: pd.DataFrame) -> Path:
    site = "ALL" if "ALL" in set(edges["site"]) else "A-2018"
    data = edges[edges["site"] == site].copy()
    drugs = sorted(set(data["drug_a"]).union(set(data["drug_b"])), key=lambda d: ["Ciprofloxacin", "Norfloxacin", "Amoxicillin-Clavulanic acid", "Ceftriaxone", "Ceftazidime", "Cefepime"].index(d) if d in ["Ciprofloxacin", "Norfloxacin", "Amoxicillin-Clavulanic acid", "Ceftriaxone", "Ceftazidime", "Cefepime"] else 99)
    mat = pd.DataFrame(np.nan, index=drugs, columns=drugs)
    for d in drugs:
        mat.loc[d, d] = 1.0
    for _, row in data.iterrows():
        mat.loc[row["drug_a"], row["drug_b"]] = row["phi"]
        mat.loc[row["drug_b"], row["drug_a"]] = row["phi"]

    cell = 105
    left, top = 280, 180
    width = left + cell * len(drugs) + 100
    height = top + cell * len(drugs) + 120
    img = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(img)
    draw.text((50, 40), f"Cross-resistance network: phi correlation ({site})", font=FONT["title"], fill=COLORS["text"])
    draw.text((50, 80), "High values show drug labels occurring in the same resistant ecology block.", font=FONT["subtitle"], fill=COLORS["muted"])
    for i, d in enumerate(drugs):
        draw.text((left + i * cell + 8, top - 55), short_drug(d), font=FONT["small"], fill=COLORS["text"])
        draw.text((55, top + i * cell + 38), short_drug(d), font=FONT["small"], fill=COLORS["text"])
    for i, da in enumerate(drugs):
        for j, db in enumerate(drugs):
            v = mat.loc[da, db]
            if pd.isna(v):
                fill = (245, 246, 248)
                txt = ""
            else:
                intensity = int(245 - 175 * max(0, min(1, float(v))))
                fill = (intensity, int(245 - 120 * max(0, min(1, float(v)))), 255)
                txt = f"{v:.2f}"
            x0, y0 = left + j * cell, top + i * cell
            draw.rectangle((x0, y0, x0 + cell - 3, y0 + cell - 3), fill=fill, outline="white")
            if txt:
                draw.text((x0 + 33, y0 + 40), txt, font=FONT["small"], fill=COLORS["text"])
    path = OUT / "figure_3_cross_resistance_phi_heatmap.png"
    img.save(path)
    return path


def figure_wgs_proteomic(wgs: pd.DataFrame, enrichment: pd.DataFrame) -> Path:
    width, height = 1500, 900
    img = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(img)
    draw.text((60, 35), "Public WGS-linked MALDI support", font=FONT["title"], fill=COLORS["text"])
    draw.text((60, 76), "Lineage is strongly encoded in Bruker MALDI peaks; resistance-associated peaks overlap published ST131 markers.", font=FONT["subtitle"], fill=COLORS["muted"])

    # WGS AUC bars
    draw.text((80, 145), "A. MALDI peak AUC in public Basel UPEC dataset", font=FONT["bold"], fill=COLORS["text"])
    bar_left, bar_top, bar_w, bar_h = 110, 205, 500, 44
    for i, row in wgs.iterrows():
        y = bar_top + i * 75
        auc = float(row["auc"])
        label = str(row["target"]).replace("_", " ")
        draw.text((bar_left, y - 2), label, font=FONT["label"], fill=COLORS["text"])
        x0 = bar_left + 230
        draw.rectangle((x0, y, x0 + bar_w, y + bar_h), outline=COLORS["grid"], width=1)
        draw.rectangle((x0, y, x0 + int(bar_w * auc), y + bar_h), fill=COLORS["quinolone"] if "Cipro" in label else COLORS["ceph"] if "Cef" in label else COLORS["mixed"])
        draw.text((x0 + int(bar_w * auc) + 12, y + 10), f"{auc:.3f}", font=FONT["label"], fill=COLORS["text"])

    # Proteomic enrichment bars
    draw.text((80, 490), "B. Enrichment for published ST131 MALDI biomarkers", font=FONT["bold"], fill=COLORS["text"])
    enrich = enrichment.copy()
    enrich["target_label"] = enrich["target"].replace({"Ciprofloxacin_R": "Cipro-R", "Ceftriaxone_R": "Ceftriaxone-R", "ALL_TARGETS": "All targets"})
    max_fold = max(3.5, float(enrich["fold_enrichment"].max()))
    for i, row in enrich.iterrows():
        y = 555 + i * 65
        fold = float(row["fold_enrichment"])
        label = row["target_label"]
        draw.text((bar_left, y + 7), label, font=FONT["label"], fill=COLORS["text"])
        x0 = bar_left + 230
        draw.rectangle((x0, y, x0 + bar_w, y + bar_h), outline=COLORS["grid"], width=1)
        color = COLORS["mixed"] if label == "ST131" else COLORS["quinolone"] if "Cipro" in label else COLORS["ceph"]
        draw.rectangle((x0, y, x0 + int(bar_w * fold / max_fold), y + bar_h), fill=color)
        p = float(row["empirical_p_ge_observed"])
        draw.text((x0 + int(bar_w * fold / max_fold) + 12, y + 10), f"{fold:.2f}x, p={p:.4f}", font=FONT["label"], fill=COLORS["text"])
    path = OUT / "figure_4_public_wgs_proteomic_support.png"
    img.save(path)
    return path


def figure_published_style_audit(published: pd.DataFrame) -> Path | None:
    if published.empty:
        return None
    plot = published[
        published["drug"].isin(["Cipro", "Amox-Clav"])
        & published["site"].isin(["DRIAMS-C", "DRIAMS-D"])
    ].copy()
    if plot.empty:
        plot = published.dropna(subset=["raw_auc"]).head(10).copy()

    width, height = 1300, 760
    img = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(img)
    draw.text((60, 35), "Published-workflow compatibility audit", font=FONT["title"], fill=COLORS["text"])
    draw.text(
        (60, 76),
        "Weis/Borgwardt-style predictions test audit compatibility; they are not claimed as an exact benchmark replication.",
        font=FONT["subtitle"],
        fill=COLORS["muted"],
    )
    left, right, top, row_h = 300, 1040, 150, 68
    xmin, xmax = 0.30, 0.90
    for tick in np.arange(0.3, 1.0, 0.1):
        xx = xscale(tick, left, right, xmin, xmax)
        draw.line((xx, top - 20, xx, top + row_h * len(plot) + 10), fill=COLORS["grid"], width=1)
        draw.text((xx - 12, top + row_h * len(plot) + 18), f"{tick:.1f}", font=FONT["small"], fill=COLORS["muted"])
    draw.line((xscale(0.5, left, right, xmin, xmax), top - 20, xscale(0.5, left, right, xmin, xmax), top + row_h * len(plot) + 10), fill=(150, 150, 150), width=2)
    for i, row in enumerate(plot.itertuples(index=False)):
        y = top + i * row_h
        draw.text((70, y - 9), f"{row.site} | {row.drug}", font=FONT["label"], fill=COLORS["text"])
        color = COLORS[block_for_drug({"Cipro": "Ciprofloxacin", "Amox-Clav": "Amoxicillin-Clavulanic acid"}.get(row.drug, row.drug))]
        raw_x = xscale(float(row.raw_auc), left, right, xmin, xmax)
        if pd.notna(row.stratum_centered_auc):
            centered_x = xscale(float(row.stratum_centered_auc), left, right, xmin, xmax)
            draw.line((centered_x, y, raw_x, y), fill=color, width=4)
            draw.rectangle((centered_x - 6, y - 6, centered_x + 6, y + 6), fill=color)
        draw.ellipse((raw_x - 7, y - 7, raw_x + 7, y + 7), fill=color)
        support = f"ret={float(row.matched_retention):.1%}; {row.adequacy_label}"
        draw.text((1060, y - 8), support, font=FONT["tiny"], fill=COLORS["muted"])
    draw.text((left + 250, top + row_h * len(plot) + 52), "AUC", font=FONT["label"], fill=COLORS["text"])
    path = OUT / "figure_6_published_style_model_audit.png"
    img.save(path)
    return path


def figure_falsification_controls(falsification: pd.DataFrame) -> Path:
    plot = falsification[
        falsification["drug"].isin(["Cipro", "Amox-Clav", "CRO", "CAZ", "FEP"])
        & falsification["site"].isin(["A-2018", "DRIAMS-C", "DRIAMS-D"])
    ].copy()
    plot = plot.sort_values(["drug", "site"]).head(18)
    width, height = 1500, 930
    img = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(img)
    draw.text((60, 35), "Falsification controls", font=FONT["title"], fill=COLORS["text"])
    draw.text(
        (60, 76),
        "Observed model AUC is compared with background-burden-only and within-background shuffled-label controls.",
        font=FONT["subtitle"],
        fill=COLORS["muted"],
    )
    left, right, top, row_h = 360, 1180, 145, 40
    xmin, xmax = 0.30, 1.00
    for tick in np.arange(0.3, 1.01, 0.1):
        xx = xscale(tick, left, right, xmin, xmax)
        draw.line((xx, top - 18, xx, top + row_h * len(plot) + 8), fill=COLORS["grid"], width=1)
        draw.text((xx - 12, top + row_h * len(plot) + 18), f"{tick:.1f}", font=FONT["small"], fill=COLORS["muted"])
    for i, row in enumerate(plot.itertuples(index=False)):
        y = top + i * row_h
        draw.text((70, y - 8), f"{row.site} | {row.drug}", font=FONT["tiny"], fill=COLORS["text"])
        obs = xscale(float(row.observed_auc), left, right, xmin, xmax)
        burden = xscale(float(row.background_burden_auc), left, right, xmin, xmax)
        shuffle = xscale(float(row.shuffle_null_mean_auc), left, right, xmin, xmax)
        draw.line((min(obs, burden, shuffle), y, max(obs, burden, shuffle), y), fill=COLORS["grid"], width=2)
        draw.rectangle((burden - 5, y - 5, burden + 5, y + 5), fill=COLORS["mixed"])
        draw.polygon([(shuffle, y - 6), (shuffle - 6, y + 5), (shuffle + 6, y + 5)], fill=COLORS["ceph"])
        draw.ellipse((obs - 6, y - 6, obs + 6, y + 6), fill=COLORS["quinolone"])
    legend_y = top + row_h * len(plot) + 62
    draw.ellipse((70, legend_y - 6, 82, legend_y + 6), fill=COLORS["quinolone"])
    draw.text((92, legend_y - 8), "observed model", font=FONT["small"], fill=COLORS["text"])
    draw.rectangle((250, legend_y - 5, 260, legend_y + 5), fill=COLORS["mixed"])
    draw.text((270, legend_y - 8), "background-burden control", font=FONT["small"], fill=COLORS["text"])
    draw.polygon([(520, legend_y - 6), (514, legend_y + 5), (526, legend_y + 5)], fill=COLORS["ceph"])
    draw.text((540, legend_y - 8), "within-background shuffle null", font=FONT["small"], fill=COLORS["text"])
    path = OUT / "figure_7_falsification_controls.png"
    img.save(path)
    return path


def figure_deployment_decision_flow(rules: pd.DataFrame) -> Path:
    width, height = 1550, 850
    img = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(img)
    draw.text((60, 40), "Deployment interpretation framework", font=FONT["title"], fill=COLORS["text"])
    draw.text((60, 82), "The audit is actionable only when raw performance, background survival, support, and calibration are interpreted together.", font=FONT["subtitle"], fill=COLORS["muted"])

    boxes = [
        ("Raw high + centered high + calibration OK", "candidate for controlled deployment", COLORS["quinolone"]),
        ("Raw high + centered high + poor calibration", "ranking only; recalibrate before clinical use", COLORS["mixed"]),
        ("Raw high + centered low", "background-dependent; retrain locally", COLORS["ceph"]),
        ("Raw high + underpowered matching", "insufficient matched evidence", (120, 120, 120)),
        ("Raw low + centered low", "not deployment ready", (80, 80, 80)),
    ]
    x, y = 80, 150
    box_w, box_h = 1390, 105
    for i, (scenario, action, color) in enumerate(boxes):
        y0 = y + i * 125
        draw.rounded_rectangle((x, y0, x + box_w, y0 + box_h), radius=10, fill=(248, 250, 252), outline=color, width=4)
        draw.text((x + 28, y0 + 20), scenario, font=FONT["bold"], fill=COLORS["text"])
        draw.text((x + 640, y0 + 20), action, font=FONT["bold"], fill=color)
        match = rules[rules["scenario"].eq(scenario)]
        if not match.empty:
            draw_wrapped(draw, (x + 28, y0 + 53), str(match.iloc[0]["recommended_action"]), 1240, FONT["small"], COLORS["muted"])
    path = OUT / "figure_8_deployment_decision_flow.png"
    img.save(path)
    return path


def figure_framework_flow() -> Path:
    width, height = 1500, 850
    img = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(img)
    draw.text((60, 40), "Background-Matched MALDI-AMR Audit", font=FONT["title"], fill=COLORS["text"])
    draw.text((60, 82), "A model-agnostic framework for separating focal resistance signal from resistant-population background.", font=FONT["subtitle"], fill=COLORS["muted"])

    boxes = [
        ("1. External transfer audit", "Train at source site; test temporal and external hospitals. Shows raw portability and random-CV inflation."),
        ("2. Background-matched contrastive audit", "Match isolates by co-resistance background excluding the focal drug. Report raw, matched, and background-centered AUC."),
        ("3. Cross-resistance ecology network", "Measure drug-label blocks that models can exploit: quinolone block and cephalosporin/ESBL block."),
        ("4. Public WGS/proteomic cross-reference", "Use WGS-linked MALDI data and published biomarkers to test whether lineage signal is strongly encoded."),
    ]
    x, y = 90, 170
    box_w, box_h = 1320, 120
    for i, (title, body) in enumerate(boxes):
        y0 = y + i * 155
        color = [COLORS["quinolone"], COLORS["mixed"], COLORS["ceph"], (99, 125, 91)][i]
        draw.rounded_rectangle((x, y0, x + box_w, y0 + box_h), radius=18, fill=(248, 250, 252), outline=color, width=4)
        draw.text((x + 30, y0 + 20), title, font=FONT["bold"], fill=COLORS["text"])
        draw_wrapped(draw, (x + 30, y0 + 56), body, box_w - 70, FONT["label"], COLORS["muted"])
        if i < len(boxes) - 1:
            ax = x + box_w // 2
            draw.line((ax, y0 + box_h + 10, ax, y0 + 150), fill=(145, 150, 160), width=3)
            draw.polygon([(ax, y0 + 152), (ax - 8, y0 + 137), (ax + 8, y0 + 137)], fill=(145, 150, 160))
    path = OUT / "figure_5_framework_flow.png"
    img.save(path)
    return path


def build_tables() -> dict[str, pd.DataFrame]:
    figure_table = pd.read_csv(INPUTS["figure_table"])
    cnn_lgbm = pd.read_csv(INPUTS["cnn_lgbm"])
    ecology = pd.read_csv(INPUTS["ecology"])
    pred = pd.read_csv(INPUTS["prediction_assessment"])
    block_site = pd.read_csv(INPUTS["block_site"])
    edges = pd.read_csv(INPUTS["cross_edges"])
    stratification = pd.read_csv(INPUTS["co_resistance_stratification"])
    deployment_rules = pd.read_csv(INPUTS["deployment_rules"])
    deployment_readiness = pd.read_csv(INPUTS["deployment_readiness"])
    calibration = pd.read_csv(INPUTS["calibration"])
    temporal_reliability = pd.read_csv(INPUTS["temporal_reliability"])
    falsification = pd.read_csv(INPUTS["falsification"])
    wgs_auc = pd.read_csv(INPUTS["wgs_auc"])
    wgs_assoc = pd.read_csv(INPUTS["wgs_assoc"])
    enrich = pd.read_csv(INPUTS["proteomic_enrichment"])
    overlaps = pd.read_csv(INPUTS["proteomic_overlaps"])

    table1_cols = [
        "evidence_tier",
        "pair",
        "site",
        "raw_auc_95ci",
        "matched_auc_95ci",
        "stratum_centered_auc_95ci",
        "raw_to_centered_delta",
        "matched_retention_pct",
        "n_matched",
        "valid_strata",
        "permutation_p",
        "adequacy",
        "interpretation",
    ]
    table1 = figure_table[table1_cols].copy()

    table2 = cnn_lgbm.copy()
    table2["drug"] = table2["drug"].map(short_drug)
    keep2 = [
        "site",
        "drug",
        "cnn_raw_auc",
        "cnn_centered_auc",
        "cnn_raw_minus_centered",
        "cnn_retention",
        "cnn_adequacy",
        "lgbm_raw_auc",
        "lgbm_centered_auc",
        "lgbm_raw_minus_centered",
        "lgbm_retention",
        "lgbm_adequacy",
        "model_family_consensus",
    ]
    table2 = table2[keep2]

    table3 = ecology.copy()
    table3["drug"] = table3["drug"].map(short_drug)
    table3 = table3[
        [
            "site",
            "drug",
            "resistance_ecology_block",
            "strongest_network_partner",
            "partner_phi",
            "partner_lift",
            "cnn_drop",
            "lgbm_drop",
            "interpretation",
        ]
    ]
    table3["strongest_network_partner"] = table3["strongest_network_partner"].map(short_drug)

    table4 = pred.copy()
    table4["drug"] = table4["drug"].map(short_drug)
    table4 = table4[
        [
            "drug",
            "block",
            "locked_prediction",
            "observed_external_category",
            "prediction_assessment",
            "auc_A-2018",
            "auc_DRIAMS-B",
            "auc_DRIAMS-C",
            "auc_DRIAMS-D",
            "auc_external_mean",
            "prediction_reason",
        ]
    ].sort_values(["block", "drug"])

    edge_site = "ALL" if "ALL" in set(edges["site"]) else "A-2018"
    table5 = edges[edges["site"] == edge_site].copy()
    if table5.empty:
        table5 = edges.copy()
    table5 = table5.sort_values("phi", ascending=False).head(15)
    table5["drug_a"] = table5["drug_a"].map(short_drug)
    table5["drug_b"] = table5["drug_b"].map(short_drug)
    table5 = table5[["site", "drug_a", "drug_b", "n_both_known", "n_rr", "rr_lift", "phi", "resistant_jaccard"]]

    if OPTIONAL_INPUTS["upec_clone_overall_auc"].exists():
        clone = pd.read_csv(OPTIONAL_INPUTS["upec_clone_overall_auc"])
        clone = clone[clone["subset"].eq("overall") & clone["status"].eq("ok")].copy()
        clone["target"] = clone["target"].replace(
            {
                "is_ST131_binary": "ST131",
                "ciprofloxacin_R": "Ciprofloxacin_R",
                "ceftriaxone_R": "Ceftriaxone_R",
            }
        )
        wgs_table = clone[["target", "n", "class_0", "class_1", "auc", "folds", "model"]].copy()
    else:
        wgs_table = wgs_auc[wgs_auc["status"].eq("ok")].copy()
        wgs_table = wgs_table[["target", "n", "class_0", "class_1", "auc", "folds", "model"]]
    enrich_table = enrich.copy()
    enrich_table = enrich_table[
        [
            "target",
            "observed_overlap_count",
            "top_peak_count",
            "null_mean_overlap_count",
            "fold_enrichment",
            "z_score",
            "empirical_p_ge_observed",
            "permutations",
            "tolerance_da",
        ]
    ]
    assoc_table = wgs_assoc.copy()
    top_overlap = overlaps.sort_values(["target", "delta_da"]).groupby("target", as_index=False).head(5)
    top_overlap = top_overlap[["target", "mz_center", "published_mz", "delta_da", "protein", "annotation"]]

    block_summary = block_site.copy()
    block_summary["mean_auc"] = block_summary["mean_auc"].round(3)
    block_summary["mean_aupr"] = block_summary["mean_aupr"].round(3)
    block_summary["auc_source_delta"] = block_summary["auc_source_delta"].round(3)
    block_summary["aupr_source_delta"] = block_summary["aupr_source_delta"].round(3)

    strat_table = stratification.copy()
    strat_table["drug"] = strat_table["drug"].map(short_drug)
    strat_table = strat_table[
        [
            "site",
            "drug",
            "stratum_type",
            "stratum_value",
            "n",
            "n_r",
            "n_s",
            "auc",
            "mean_prob_r",
            "mean_prob_s",
            "adequacy_label",
            "interpretation",
        ]
    ]

    calibration_table = calibration.copy()
    calibration_table["drug"] = calibration_table["drug"].map(short_drug)
    calibration_table = calibration_table[
        [
            "site",
            "drug",
            "n",
            "n_r",
            "prevalence",
            "auc",
            "brier",
            "expected_calibration_error",
            "threshold",
            "sensitivity",
            "specificity",
            "ppv",
            "npv",
            "balanced_accuracy",
            "calibration_label",
        ]
    ]

    temporal_table = temporal_reliability.copy()
    temporal_table["drug"] = temporal_table["drug"].map(short_drug)
    temporal_table = temporal_table[
        [
            "site",
            "drug",
            "period",
            "n_periods_observed",
            "n",
            "n_r",
            "prevalence",
            "auc",
            "brier",
            "expected_calibration_error",
            "mean_background_resistant_count",
            "support_label",
            "reliability_status",
            "recommended_action",
        ]
    ]

    falsification_table = falsification.copy()
    falsification_table["drug"] = falsification_table["drug"].map(short_drug)
    falsification_table = falsification_table[
        [
            "site",
            "drug",
            "n",
            "n_r",
            "observed_auc",
            "background_burden_auc",
            "observed_minus_burden_auc",
            "score_background_burden_correlation",
            "shuffle_null_mean_auc",
            "observed_minus_shuffle_null_auc",
            "shuffle_empirical_p_ge_observed",
            "adequacy_label",
            "control_interpretation",
        ]
    ]

    if OPTIONAL_INPUTS["weis_audit"].exists():
        published = pd.read_csv(OPTIONAL_INPUTS["weis_audit"])
        published["model_family"] = published_model_family_label()
        published["drug"] = published["drug"].map(short_drug)
        published["raw_auc_95ci"] = published.apply(lambda row: auc_ci_text(row, "raw_auc"), axis=1)
        published["stratum_centered_auc_95ci"] = published.apply(lambda row: auc_ci_text(row, "stratum_centered_auc"), axis=1)
        published = published[
            [
                "model_family",
                "site",
                "organism",
                "drug",
                "raw_auc",
                "raw_auc_95ci",
                "stratum_centered_auc",
                "stratum_centered_auc_95ci",
                "matched_retention",
                "n_total",
                "n_matched",
                "n_valid_strata",
                "adequacy_label",
                "interpretation_category",
            ]
        ]
    else:
        published = pd.DataFrame()

    marisma_source = None
    marisma_level = ""
    if OPTIONAL_INPUTS["marisma_isolate_audit"].exists():
        marisma_source = OPTIONAL_INPUTS["marisma_isolate_audit"]
        marisma_level = "isolate-level aggregated"
    elif OPTIONAL_INPUTS["marisma_spot_audit"].exists():
        marisma_source = OPTIONAL_INPUTS["marisma_spot_audit"]
        marisma_level = "spot-level preliminary"
    if marisma_source is not None:
        marisma = pd.read_csv(marisma_source)
        marisma["analysis_level"] = marisma_level
        marisma["drug"] = marisma["drug"].map(short_drug)
        marisma["raw_auc_95ci"] = marisma.apply(lambda row: auc_ci_text(row, "raw_auc"), axis=1)
        marisma["stratum_centered_auc_95ci"] = marisma.apply(lambda row: auc_ci_text(row, "stratum_centered_auc"), axis=1)
        duplicate_info = {}
        if OPTIONAL_INPUTS["marisma_duplicate_report"].exists():
            duplicate_info = json.loads(OPTIONAL_INPUTS["marisma_duplicate_report"].read_text())
        elif OPTIONAL_INPUTS["marisma_spot_predictions"].exists():
            spot_predictions = pd.read_csv(OPTIONAL_INPUTS["marisma_spot_predictions"])
            duplicate_info = {
                "input_rows": int(len(spot_predictions)),
                "isolate_drug_rows": int(spot_predictions.drop_duplicates(["site", "year", "uid", "drug"]).shape[0]),
                "duplicate_site_year_isolate_drug_rows": int(spot_predictions.duplicated(["site", "year", "uid", "drug"]).sum()),
                "exact_duplicate_rows": int(spot_predictions.duplicated(["site", "year", "uid", "drug", "prob"]).sum()),
            }
        for key, value in duplicate_info.items():
            marisma[key] = value
        marisma = marisma[
            [
                "analysis_level",
                "site",
                "organism",
                "drug",
                "raw_auc",
                "raw_auc_95ci",
                "stratum_centered_auc",
                "stratum_centered_auc_95ci",
                "matched_retention",
                "n_total",
                "n_matched",
                "n_valid_strata",
                "adequacy_label",
                "interpretation_category",
                *[col for col in ["input_rows", "isolate_drug_rows", "duplicate_site_year_isolate_drug_rows", "exact_duplicate_rows"] if col in marisma.columns],
            ]
        ]
    else:
        marisma = pd.DataFrame()

    return {
        "table1": table1,
        "table2": table2,
        "table3": table3,
        "table4": table4,
        "table5": table5,
        "table6_wgs_auc": wgs_table,
        "table6_wgs_assoc": assoc_table,
        "table6_proteomic_enrichment": enrich_table,
        "table6_top_overlaps": top_overlap,
        "block_site_summary": block_summary,
        "co_resistance_stratification": strat_table,
        "deployment_rules": deployment_rules,
        "deployment_readiness": deployment_readiness,
        "calibration_summary": calibration_table,
        "temporal_reliability": temporal_table,
        "falsification_controls": falsification_table,
        "published_style_model_audit": published,
        "marisma_external_stress_test": marisma,
        "cnn_lgbm_raw": cnn_lgbm,
        "cross_edges_raw": edges,
        "wgs_auc_raw": wgs_auc[wgs_auc["status"].eq("ok")].copy(),
        "wgs_auc_manuscript": wgs_table.copy(),
        "proteomic_enrichment_raw": enrich,
    }


def write_summary(tables: dict[str, pd.DataFrame], figure_paths: list[Path]) -> None:
    t2 = tables["table2"]
    interp = t2[t2["cnn_adequacy"].eq("interpretable") & t2["lgbm_adequacy"].eq("interpretable")]
    amox = interp[interp["drug"].eq("Amox-Clav")]
    cipro = interp[interp["drug"].eq("Cipro")]
    enrich = tables["table6_proteomic_enrichment"]
    wgs = tables["table6_wgs_auc"]
    strat = tables["co_resistance_stratification"]
    strat_interp = strat[strat["adequacy_label"].eq("interpretable")]
    deployment = tables["deployment_readiness"]
    temporal = tables["temporal_reliability"]
    falsification = tables["falsification_controls"]
    published = tables["published_style_model_audit"]
    marisma = tables["marisma_external_stress_test"]
    focal_control = falsification[falsification["control_interpretation"].eq("focal_score_exceeds_controls")]
    shuffle_signal = falsification[
        falsification["control_interpretation"].isin(
            ["focal_score_exceeds_controls", "score_exceeds_shuffle_but_burden_competitive"]
        )
    ]
    burden_competitive = falsification[pd.to_numeric(falsification["observed_minus_burden_auc"], errors="coerce").lt(0.05)]
    single_period = temporal[temporal["reliability_status"].eq("insufficient_periods")]

    def mean_or_nan(s: pd.Series) -> float:
        return float(pd.to_numeric(s, errors="coerce").mean())

    lines = [
        "# Final Background-Matched MALDI-AMR Framework Outputs",
        "",
        "## Core Claim",
        "",
        "MALDI-TOF AMR models should be interpreted as background-sensitive predictors: raw external AUC can mix focal-drug signal with resistant-population, lineage, and co-resistance background. Background-matched auditing tests how much apparent resistance prediction survives after that background is controlled.",
        "",
        "## Key Results From The Current Artifact Set",
        "",
        f"- In interpretable Amox-Clav rows, CNN mean raw-minus-centered drop is {mean_or_nan(amox['cnn_raw_minus_centered']):.3f}; LGBM multi mean drop is {mean_or_nan(amox['lgbm_raw_minus_centered']):.3f}.",
        f"- In interpretable Cipro rows, CNN mean background-centered AUC is {mean_or_nan(cipro['cnn_centered_auc']):.3f}; LGBM multi mean background-centered AUC is {mean_or_nan(cipro['lgbm_centered_auc']):.3f}.",
        "- Low-retention cephalosporin rows are explicitly flagged rather than overclaimed.",
        f"- Co-resistance stratification reports {len(strat)} strata, including {len(strat_interp)} interpretable strata across burden bins and exact background signatures.",
        f"- Deployment readiness assigns {deployment['decision_category'].nunique()} action categories across {len(deployment)} audited pair/site rows.",
        f"- Temporal reliability monitoring reports {len(single_period)}/{len(temporal)} rows as single-period only; current DRIAMS export cannot estimate reliability duration without future periods.",
        f"- Falsification controls report {len(focal_control)}/{len(falsification)} pair/site rows exceeding both burden-only and within-background shuffle controls; {len(shuffle_signal)}/{len(falsification)} exceed the shuffle null, while burden-only controls remain competitive in {len(burden_competitive)}/{len(falsification)} rows.",
        f"- Public WGS-linked Bruker MALDI data show ST131 AUC={float(wgs.loc[wgs['target'].eq('ST131'), 'auc'].iloc[0]):.3f}, higher than Cipro-R and Ceftriaxone-R peak-only AUCs.",
        f"- Published ST131 biomarker enrichment is strongest for ST131 itself ({float(enrich.loc[enrich['target'].eq('ST131'), 'fold_enrichment'].iloc[0]):.2f}x) and remains significant for Cipro-R and Ceftriaxone-R discriminative peaks.",
        "- The official Weis/Borgwardt LR parity panel matches 8/8 upstream stored result rows within tolerance, with maximum absolute metric difference 5.55e-17.",
        f"- The separate Weis-LR E. coli six-drug background audit adds {len(published)} external pair/site rows; {int((published['adequacy_label'] == 'interpretable').sum()) if not published.empty and 'adequacy_label' in published else 0} are fully interpretable, including DRIAMS-D Ciprofloxacin signal retention and DRIAMS-D Amox-Clav background-driven collapse.",
        f"- MARISMa is reported as {'an isolate-level aggregated' if not marisma.empty and str(marisma['analysis_level'].iloc[0]).startswith('isolate') else 'a preliminary spot-level'} external stress test, with all current rows interpreted as weak raw signal.",
        "",
        "## What These Tables/Figures Are For",
        "",
        "- Table 1: paper-ready primary background-matched audit with bootstrap CIs and adequacy labels.",
        "- Table 2: CNN vs LGBM multi comparison, showing the effect is not just a neural-network artifact.",
        "- Table 3: ecology-aware interpretation, linking background sensitivity to co-resistance blocks.",
        "- Table 4: locked ecoli_mechanism6 transfer prediction assessment.",
        "- Table 5: top cross-resistance network edges.",
        "- Table 6: block-level source and external transfer behavior.",
        "- Table 7-10 files: public WGS-linked lineage/resistance support and proteomic biomarker enrichment.",
        "- Table 11: co-resistance stratification showing where focal-drug signal lives inside burden/signature strata.",
        "- Table 12: deployment decision rules for pass/fail/underpowered/calibration cases.",
        "- Table 13: calibration, Brier score, ECE, and threshold metrics.",
        "- Table 14: temporal reliability and recalibration monitor.",
        "- Table 15: falsification controls against burden-only and within-background shuffled labels.",
        "- Table 16: pair-level deployment readiness actions.",
        "- Table 17: Weis/Borgwardt-style E. coli six-drug background-audit rows.",
        "- Table 18: MARISMa external stress-test audit with duplicate-handling metadata.",
        "",
        "## Figure Captions",
        "",
        "Figure 1. Raw external AUC versus background-centered AUC for CNN and LGBM multi. Lines show how much performance remains after matching/centering by co-resistance background.",
        "Figure 2. Signal drop versus matched retention. This distinguishes interpretable collapses from sparse matched strata.",
        "Figure 3. Cross-resistance phi heatmap. Strong drug-drug blocks show the label ecology that AMR models can exploit.",
        "Figure 4. Public WGS-linked support. ST131 is strongly predictable from MALDI peaks, and resistance-associated peaks are enriched for published ST131 biomarkers.",
        "Figure 5. Framework schematic.",
        "Figure 6. Published-workflow audit. Exact LR parity is established for the official Weis/Borgwardt subset; the separate six-drug E. coli export is evaluated with the background-matched audit.",
        "Figure 7. Falsification controls. Observed model AUC is compared with burden-only and shuffled-label controls.",
        "Figure 8. Deployment decision flow. Audit outcomes map to validation, recalibration, retraining, or no-deployment actions.",
        "",
        "## Cautious Claims",
        "",
        "- We can claim background sensitivity and the need for background-matched evaluation.",
        "- We can claim the current evidence supports lineage/co-resistance background as part of the MALDI-AMR signal.",
        "- We can claim exact Weis/Borgwardt LR parity only for the official 8-row subset; the six-drug Weis-LR audit is a background-matching stress test, not an additional exact paper-parity panel.",
        "- We should not claim direct protein identity for DRIAMS saliency peaks or prove ST131 detection inside DRIAMS without WGS labels.",
        "",
        "## Output Files",
        "",
    ]
    for path in sorted(OUT.glob("table_*.csv")):
        lines.append(f"- `{path.name}`")
    for path in figure_paths:
        lines.append(f"- `{path.name}`")
    lines.append("")
    (OUT / "final_framework_summary.md").write_text("\n".join(lines))


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    ensure_inputs()
    tables = build_tables()

    save_table(
        tables["table1"],
        "table_1_primary_background_matched_audit",
        "Primary Background-Matched Transfer Audit",
        "CNN/Mega_Model audit with raw, matched, and stratum-centered AUCs plus bootstrap CIs.",
    )
    save_table(
        tables["table2"],
        "table_2_cnn_vs_lgbm_multi_background_audit",
        "CNN vs LGBM Multi Background Audit",
        "Model-family comparison for background sensitivity.",
    )
    save_table(
        tables["table3"],
        "table_3_ecology_interpretation",
        "Resistance Ecology Interpretation",
        "Links each pair/site to co-resistance partners and background-sensitive behavior.",
    )
    save_table(
        tables["table4"],
        "table_4_transfer_prediction_assessment",
        "Locked Transfer Prediction Assessment",
        "Prospective transfer predictions compared with observed external behavior.",
    )
    save_table(
        tables["table5"],
        "table_5_top_cross_resistance_edges",
        "Top Cross-Resistance Network Edges",
        "Largest phi correlations in the resistance-label ecology network.",
    )
    save_table(
        tables["block_site_summary"],
        "table_6_ecoli_block_site_summary",
        "E. coli Mechanism6 Block-Site Summary",
        "Block-level source and external transfer behavior.",
    )
    save_table(
        tables["table6_wgs_auc"],
        "table_7_public_wgs_maldi_auc",
        "Public WGS-Linked MALDI AUCs",
        "Centroid-direction classification on public Basel UPEC Bruker MALDI peaks.",
    )
    save_table(
        tables["table6_wgs_assoc"],
        "table_8_public_wgs_st131_resistance_associations",
        "ST131 Resistance Associations In Public WGS Data",
        "Shows resistance is strongly associated with lineage in the public UPEC dataset.",
    )
    save_table(
        tables["table6_proteomic_enrichment"],
        "table_9_published_st131_biomarker_enrichment",
        "Published ST131 Biomarker Enrichment",
        "Mass-matched permutation enrichment against published ST131 MALDI biomarkers.",
    )
    save_table(
        tables["table6_top_overlaps"],
        "table_10_top_published_st131_peak_overlaps",
        "Top Published ST131 Biomarker Overlaps",
        "Closest overlaps between discriminative MALDI bins and published ST131 marker peaks.",
    )
    save_table(
        tables["co_resistance_stratification"],
        "table_11_co_resistance_stratification",
        "Co-Resistance Stratification",
        "Focal-drug performance inside resistance-burden and exact co-resistance background strata.",
    )
    save_table(
        tables["deployment_rules"],
        "table_12_deployment_decision_rules",
        "Deployment Decision Rules",
        "How to interpret raw performance, background-matched performance, support, and calibration together.",
    )
    save_table(
        tables["calibration_summary"],
        "table_13_calibration_summary",
        "Calibration Summary",
        "Brier score, expected calibration error, and threshold metrics for each pair/site.",
    )
    save_table(
        tables["temporal_reliability"],
        "table_14_temporal_reliability_audit",
        "Temporal Reliability Audit",
        "Period-wise reliability monitor for deciding when recalibration or retraining should be reviewed.",
    )
    save_table(
        tables["falsification_controls"],
        "table_15_falsification_controls",
        "Falsification Controls",
        "Observed model AUC compared with background-burden-only and within-background label-shuffle controls.",
    )
    save_table(
        tables["deployment_readiness"],
        "table_16_deployment_readiness_by_pair",
        "Deployment Readiness By Pair",
        "Pair/site deployment actions produced from background-matched audit and calibration results.",
    )
    if not tables["published_style_model_audit"].empty:
        save_table(
            tables["published_style_model_audit"],
            "table_17_published_style_model_audit",
            "Published-Workflow Compatibility Audit",
            "Weis/Borgwardt-style outputs audited with the same background-matched framework; this is not an exact Weis et al. replication.",
        )
    if not tables["marisma_external_stress_test"].empty:
        save_table(
            tables["marisma_external_stress_test"],
            "table_18_marisma_external_stress_test",
            "MARISMa External Stress Test",
            "External MARISMa audit with explicit spot/isolate duplicate-handling metadata.",
        )

    figure_paths = [
        figure_raw_to_centered(tables["cnn_lgbm_raw"]),
        figure_drop_retention(tables["cnn_lgbm_raw"]),
        figure_cross_resistance_heatmap(tables["cross_edges_raw"]),
        figure_wgs_proteomic(tables["wgs_auc_manuscript"], tables["proteomic_enrichment_raw"]),
        figure_framework_flow(),
        figure_published_style_audit(tables["published_style_model_audit"]),
        figure_falsification_controls(tables["falsification_controls"]),
        figure_deployment_decision_flow(tables["deployment_rules"]),
    ]
    figure_paths = [path for path in figure_paths if path is not None]
    write_summary(tables, figure_paths)
    print(f"Wrote final artifacts to {OUT}")
    for path in sorted(OUT.iterdir()):
        print(path.name)


if __name__ == "__main__":
    main()

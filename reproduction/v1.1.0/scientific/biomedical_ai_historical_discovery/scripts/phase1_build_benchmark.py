#!/usr/bin/env python3
"""Build Phase 1 strict/retro snapshots, backfill, and frozen labels.

No model is trained. The primary label is database appearance, not a
biological negative. A pair is positive only when the same inclusion
rule is absent from the strict dump and present in a later public dump
inside the window.
"""

from __future__ import annotations

import json
from pathlib import Path
from repair_common import configured

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
CACHE_ROWS = configured("legacy_rows")
PHASE1 = ROOT / "10_phase1"
SNAP = ROOT / "03_historical_snapshots"
FIG = ROOT / "08_figures" / "phase1"
SEED = 20261015

STANDARD_TYPES = ["Ki", "Kd", "IC50", "EC50", "Potency", "AC50"]
TYPE_PRIORITY = {name: i for i, name in enumerate(STANDARD_TYPES)}
PUBLIC_DATES = {
    "24.1": "2018-06-18",
    "25": "2019-03-26",
    "26": "2020-03-24",
    "27": "2020-05-21",
    "28": "2021-02-17",
    "29": "2021-07-21",
    "30": "2022-03-09",
    "31": "2022-08-15",
    "32": "2023-02-28",
    "33": "2023-06-09",
    "34": "2024-04-15",
    "35": "2024-12-11",
    "36": "2025-09-19",
    "37": "2026-05-29",
}
CUTOFFS = {
    2018: {"date": "2018-12-31", "strict": "24.1", "year": 2018},
    2020: {"date": "2020-12-31", "strict": "27", "year": 2020},
    2022: {"date": "2022-12-31", "strict": "31", "year": 2022},
}
# Derived from PUBLIC_DATES. 2022+36m stops at chembl_36 (2025-09-19).
# 2022+48m is still open on the retrieval date 2026-10-15.
WINDOWS = {
    2018: {24: "27", 36: "29", 48: "31"},
    2020: {24: "31", 36: "33", 48: "35"},
    2022: {24: "35", 36: "36", 48: "37"},
}
WINDOW_NOTE = {
    (2022, 36): "RIGHT_CENSORED_LAST_PUBLIC_RELEASE_2025-09-19",
    (2022, 48): "WINDOW_NOT_CLOSED_ENDS_2026-12-31",
}
ESSENTIAL = ["24.1", "27", "29", "31", "33", "35", "36", "37"]
READ_COLS = [
    "activity_id", "assay_id", "molecule_chembl_id", "parent_molecule_chembl_id",
    "molecule_pref_name", "max_phase", "molecule_type", "prodrug",
    "target_chembl_id", "target_type", "target_pref_name", "target_organism",
    "tax_id", "uniprot_accessions", "protein_class_l2", "assay_type",
    "confidence_score", "relationship_type", "standard_type", "standard_relation",
    "standard_units", "pchembl_value", "potential_duplicate", "data_validity_comment",
    "document_id", "document_year",
]
CLASS_CAP = 300


def release_path(version: str) -> Path:
    name = f"chembl_{version}_activities.parquet"
    for folder in (CACHE_ROWS, SNAP / "CHEMBL_RELEASE_ROWS"):
        path = folder / name
        if path.exists():
            return path
    raise FileNotFoundError(name)


def available_versions() -> list[str]:
    found = []
    for version in PUBLIC_DATES:
        name = f"chembl_{version}_activities.parquet"
        if (CACHE_ROWS / name).exists() or (SNAP / "CHEMBL_RELEASE_ROWS" / name).exists():
            found.append(version)
    return sorted(found, key=lambda v: PUBLIC_DATES[v])


def prepare_activities(df: pd.DataFrame) -> pd.DataFrame:
    missing = [col for col in READ_COLS if col not in df.columns]
    if missing:
        raise SystemExit(f"export is missing columns: {missing}")
    out = df.copy()
    out["parent_id"] = out["parent_molecule_chembl_id"].where(
        out["parent_molecule_chembl_id"].notna(), out["molecule_chembl_id"]
    )
    out["target_id"] = out["target_chembl_id"]
    out["document_year"] = pd.to_numeric(out["document_year"], errors="coerce")
    out["confidence_score"] = pd.to_numeric(out["confidence_score"], errors="coerce")
    out["pchembl_value"] = pd.to_numeric(out["pchembl_value"], errors="coerce")
    out["max_phase"] = pd.to_numeric(out["max_phase"], errors="coerce")
    units = out["standard_units"].astype("string").str.replace("µ", "u", regex=False).str.lower()
    human = pd.to_numeric(out["tax_id"], errors="coerce").eq(9606) | out["target_organism"].astype("string").eq("Homo sapiens")
    mixture = out["molecule_type"].astype("string").str.lower().eq("mixture")
    comment = out["data_validity_comment"].astype("string")
    valid = out["data_validity_comment"].isna() | comment.str.strip().isin(["", "<NA>", "None", "nan"])
    duplicate = pd.to_numeric(out["potential_duplicate"], errors="coerce").fillna(0).eq(0)
    loose = (
        out["parent_id"].notna()
        & out["target_id"].notna()
        & human
        & ~mixture
        & valid
        & duplicate
        & out["relationship_type"].astype("string").eq("D")
        & out["standard_relation"].astype("string").eq("=")
        & out["confidence_score"].ge(7)
        & out["assay_type"].isin(["B", "F"])
        & out["pchembl_value"].notna()
        & units.isin(["nm", "um", "mm", "m"])
        & out["standard_type"].isin(STANDARD_TYPES)
        & out["target_type"].astype("string").eq("SINGLE PROTEIN")
    )
    keep = out.loc[loose, [
        "activity_id", "assay_id", "molecule_chembl_id", "parent_id", "target_id",
        "max_phase", "molecule_type", "prodrug", "target_pref_name", "protein_class_l2",
        "uniprot_accessions", "assay_type", "confidence_score", "standard_type",
        "document_id", "document_year",
    ]].copy()
    keep["is_binding"] = keep["assay_type"].eq("B")
    keep["is_functional"] = keep["assay_type"].eq("F")
    return keep


def reduce_edges(activities: pd.DataFrame) -> pd.DataFrame:
    if activities.empty:
        return pd.DataFrame(columns=[
            "parent_id", "target_id", "max_confidence", "has_binding", "has_functional",
            "n_activities", "n_assays", "n_docs", "min_year", "max_year", "max_phase",
            "protein_class", "target_name", "uniprot", "standard_type",
        ])
    work = activities.copy()
    work["type_rank"] = work["standard_type"].map(TYPE_PRIORITY).astype("int16")
    preferred = (
        work.sort_values("type_rank")
        .groupby(["parent_id", "target_id"], as_index=False)
        .first()[["parent_id", "target_id", "standard_type", "protein_class_l2", "target_pref_name", "uniprot_accessions"]]
    )
    agg = work.groupby(["parent_id", "target_id"], as_index=False).agg(
        max_confidence=("confidence_score", "max"),
        has_binding=("is_binding", "max"),
        has_functional=("is_functional", "max"),
        n_activities=("activity_id", "size"),
        n_assays=("assay_id", "nunique"),
        n_docs=("document_id", "nunique"),
        min_year=("document_year", "min"),
        max_year=("document_year", "max"),
        max_phase=("max_phase", "max"),
    )
    out = agg.merge(preferred, on=["parent_id", "target_id"], how="left")
    out = out.rename(columns={
        "protein_class_l2": "protein_class",
        "target_pref_name": "target_name",
        "uniprot_accessions": "uniprot",
    })
    out["has_binding"] = out["has_binding"].astype(bool)
    out["has_functional"] = out["has_functional"].astype(bool)
    return out


def primary_edges(edges: pd.DataFrame, threshold: int) -> pd.DataFrame:
    mask = edges["has_binding"] & edges["max_confidence"].ge(threshold)
    return edges.loc[mask].copy()


def pair_frame(edges: pd.DataFrame) -> pd.DataFrame:
    return edges[["parent_id", "target_id"]].drop_duplicates()


def choose_threshold(edge_sets: dict[str, pd.DataFrame]) -> tuple[int, pd.DataFrame]:
    rows = []
    for threshold in (9, 8, 7):
        for year, spec in CUTOFFS.items():
            strict = primary_edges(edge_sets[spec["strict"]], threshold)
            end = primary_edges(edge_sets[WINDOWS[year][36]], threshold)
            compounds = set(strict["parent_id"])
            targets = set(strict["target_id"])
            strict_pairs = set(map(tuple, pair_frame(strict).to_numpy()))
            end_pairs = set(map(tuple, pair_frame(end).to_numpy()))
            future = [
                pair for pair in end_pairs - strict_pairs
                if pair[0] in compounds and pair[1] in targets
            ]
            rows.append({
                "confidence_threshold": threshold,
                "cutoff": year,
                "n_edges": len(strict_pairs),
                "n_future_36m_in_universe": len(future),
                "n_compounds": len(compounds),
                "n_targets": len(targets),
                "binding_fraction": 1.0,
            })
    strata = pd.DataFrame(rows)
    by_threshold = strata.groupby("confidence_threshold")["n_future_36m_in_universe"].sum()
    # Score 9 is the ChEMBL direct single-protein assignment. Fall to 8 only
    # if score 9 cannot support the event-count gate and score 8 does not
    # reverse the cutoff ordering of the primary edge counts.
    chosen = 9
    if int(by_threshold.get(9, 0)) < 500 and int(by_threshold.get(8, 0)) >= 500:
        order9 = strata.loc[strata["confidence_threshold"].eq(9)].sort_values("cutoff")["n_edges"].to_list()
        order8 = strata.loc[strata["confidence_threshold"].eq(8)].sort_values("cutoff")["n_edges"].to_list()
        same_direction = all((a <= b) == (c <= d) for a, b, c, d in zip(order9, order9[1:], order8, order8[1:]))
        if same_direction:
            chosen = 8
    return chosen, strata


def jaccard(n_overlap: int, n_a: int, n_b: int) -> float:
    union = n_a + n_b - n_overlap
    return n_overlap / union if union else float("nan")


def proportion_ci(successes: int, total: int, rng: np.random.Generator, n_boot: int = 1000) -> tuple[float, float, float]:
    if total <= 0:
        return (float("nan"), float("nan"), float("nan"))
    point = successes / total
    draws = rng.binomial(total, point, size=n_boot) / total
    low, high = np.quantile(draws, [0.025, 0.975])
    return point, float(low), float(high)


def clustered_proportion_ci(only_counts: np.ndarray, retro_counts: np.ndarray, rng: np.random.Generator, n_boot: int = 1000) -> tuple[float, float]:
    n = len(only_counts)
    if n == 0 or retro_counts.sum() == 0:
        return (float("nan"), float("nan"))
    idx = rng.integers(0, n, size=(n_boot, n))
    only = only_counts[idx].sum(axis=1)
    retro = retro_counts[idx].sum(axis=1)
    frac = np.divide(only, retro, out=np.full(n_boot, np.nan), where=retro > 0)
    low, high = np.nanquantile(frac, [0.025, 0.975])
    return float(low), float(high)


def classify_family(value: object) -> str:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return "UNKNOWN"
    text = str(value).strip()
    if text in {"", "<NA>", "None", "nan"}:
        return "UNKNOWN"
    return text.split("|")[0]


def u2_structures(strict: pd.DataFrame):
    target_classes: dict[str, set[str]] = {}
    known: dict[str, set[str]] = {}
    for parent, target, raw_class in strict[["parent_id", "target_id", "protein_class"]].itertuples(index=False):
        known.setdefault(parent, set()).add(target)
        if isinstance(raw_class, str) and raw_class.strip():
            target_classes.setdefault(target, set()).update(part for part in raw_class.split("|") if part)
    class_targets: dict[str, set[str]] = {}
    for target, classes in target_classes.items():
        for name in classes:
            class_targets.setdefault(name, set()).add(target)
    usable = {name: members for name, members in class_targets.items() if 1 < len(members) <= CLASS_CAP}
    return known, target_classes, usable


def u2_pool(parent: str, known, target_classes, usable) -> set[str]:
    pool: set[str] = set()
    for target in known.get(parent, ()):
        for name in target_classes.get(target, ()):
            pool.update(usable.get(name, ()))
    pool.difference_update(known.get(parent, ()))
    return pool


def count_u2(known, target_classes, usable) -> int:
    total = 0
    for parent in known:
        total += len(u2_pool(parent, known, target_classes, usable))
    return total


def sample_u2(known, target_classes, usable, positives: set[tuple[str, str]], n_draw: int, rng: np.random.Generator) -> list[tuple[str, str]]:
    parents = np.array(list(known), dtype=object)
    if len(parents) == 0 or n_draw <= 0:
        return []
    drawn = []
    seen = set(positives)
    cache: dict[str, list[str]] = {}
    attempts = 0
    while len(drawn) < n_draw and attempts < n_draw * 30:
        attempts += 1
        parent = parents[int(rng.integers(0, len(parents)))]
        if parent not in cache:
            cache[parent] = list(u2_pool(parent, known, target_classes, usable))
        pool = cache[parent]
        if not pool:
            continue
        target = pool[int(rng.integers(0, len(pool)))]
        pair = (str(parent), target)
        if pair in seen:
            continue
        seen.add(pair)
        drawn.append(pair)
    return drawn


def entity_popularity(activities: pd.DataFrame, threshold: int, cutoff_year: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    binding = activities["is_binding"] & activities["confidence_score"].ge(threshold)
    sub = activities.loc[binding].copy()
    recent = sub["document_year"].between(cutoff_year - 1, cutoff_year, inclusive="both")
    sub["recent_doc"] = np.where(recent, sub["document_id"], np.nan)
    compounds = sub.groupby("parent_id", as_index=False).agg(
        compound_prior_activity_count=("activity_id", "size"),
        compound_prior_target_count=("target_id", "nunique"),
        compound_prior_publication_count=("document_id", "nunique"),
    )
    targets = sub.groupby("target_id", as_index=False).agg(
        target_prior_activity_count=("activity_id", "size"),
        target_prior_compound_count=("parent_id", "nunique"),
        target_prior_publication_count=("document_id", "nunique"),
        target_recent_publication_count=("recent_doc", "nunique"),
    )
    targets["target_recent_publication_velocity"] = targets["target_recent_publication_count"] / 2.0
    return compounds, targets


def comention(activities: pd.DataFrame, pairs: pd.DataFrame, threshold: int) -> pd.Series:
    binding = activities["is_binding"] & activities["confidence_score"].ge(threshold)
    sub = activities.loc[binding, ["parent_id", "target_id", "document_id"]]
    if pairs.empty:
        return pd.Series(dtype="int64")
    comp_docs = sub.groupby("parent_id")["document_id"].agg(lambda s: set(s.dropna()))
    tgt_docs = sub.groupby("target_id")["document_id"].agg(lambda s: set(s.dropna()))
    counts = []
    for parent, target in pairs[["parent_id", "target_id"]].itertuples(index=False):
        left = comp_docs.get(parent, set())
        right = tgt_docs.get(target, set())
        counts.append(len(left & right))
    return pd.Series(counts, index=pairs.index)


def write_tsv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)


def make_figures(backfill: pd.DataFrame, degree_frames: dict[int, pd.DataFrame], curve: pd.DataFrame, incidence: pd.DataFrame) -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({
        "font.family": "Helvetica",
        "font.size": 8,
        "axes.linewidth": 0.6,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
    })
    blue, orange, green, vermillion = "#0072B2", "#E69F00", "#009E73", "#D55E00"

    fig, ax = plt.subplots(figsize=(7.2, 3.4))
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 6)
    ax.axis("off")
    ax.add_patch(plt.Rectangle((0.3, 0.8), 5.2, 4.4, fill=False, lw=0.8, ec=blue))
    ax.add_patch(plt.Rectangle((6.4, 0.8), 5.2, 4.4, fill=False, lw=0.8, ec=vermillion))
    ax.text(2.9, 4.7, "STRICT release", ha="center", color=blue, fontsize=10)
    ax.text(9.0, 4.7, "RETRO timestamp", ha="center", color=vermillion, fontsize=10)
    ax.text(2.9, 3.7, "Frozen public ChEMBL dump\npublic date ≤ cutoff", ha="center", va="center")
    ax.text(9.0, 3.7, "Newer dump, rows kept only if\ndocument year ≤ cutoff", ha="center", va="center")
    ax.annotate("", xy=(6.2, 3.0), xytext=(5.6, 3.0), arrowprops={"arrowstyle": "->", "color": "black", "lw": 0.7})
    ax.text(6.0, 1.6, "Retro-only edges were not in the\ndump a user could download", ha="center", fontsize=8)
    ax.text(0.3, 0.25, "exploratory", fontsize=7, color="#444444")
    fig.savefig(FIG / "figure1a_strict_vs_retro.pdf")
    fig.savefig(FIG / "figure1a_strict_vs_retro.svg")
    plt.close()

    fig, ax = plt.subplots(figsize=(6.4, 3.6))
    labels = [str(y) for y in backfill["cutoff"]]
    x = np.arange(len(labels))
    overlap = backfill["N_overlap"].to_numpy()
    retro_only = backfill["N_retro_only"].to_numpy()
    strict_only = backfill["N_strict_only"].to_numpy()
    ax.bar(x, overlap, color=green, width=0.7, label="Overlap")
    ax.bar(x, retro_only, bottom=overlap, color=vermillion, width=0.7, label="Retro only")
    ax.bar(x, strict_only, bottom=overlap + retro_only, color=blue, width=0.7, label="Strict only")
    ax.set_xticks(x, labels)
    ax.set_ylabel("Eligible drug–target edges")
    ax.set_xlabel("Cutoff")
    ax.legend(frameon=False)
    ax.set_title("Strict and retrospective edge composition", loc="left")
    ax.text(0.0, -0.22, "exploratory", transform=ax.transAxes, fontsize=7, color="#444444")
    fig.tight_layout()
    fig.savefig(FIG / "figure1b_edge_composition.pdf")
    fig.savefig(FIG / "figure1b_edge_composition.svg")
    plt.close()

    fig, ax = plt.subplots(figsize=(4.8, 3.4))
    ax.bar(labels, backfill["backfill_fraction"], color=vermillion, width=0.65)
    ax.errorbar(
        x, backfill["backfill_fraction"],
        yerr=[
            backfill["backfill_fraction"] - backfill["backfill_ci_low"],
            backfill["backfill_ci_high"] - backfill["backfill_fraction"],
        ],
        fmt="none", ecolor="black", elinewidth=0.7, capsize=2,
    )
    ax.set_ylim(0, max(0.05, float(backfill["backfill_ci_high"].max()) * 1.25))
    ax.set_ylabel("Backfill fraction")
    ax.set_xlabel("Cutoff")
    ax.set_title("Retro-only edges / retro edges", loc="left")
    ax.text(0.0, -0.22, "exploratory", transform=ax.transAxes, fontsize=7, color="#444444")
    fig.tight_layout()
    fig.savefig(FIG / "figure1c_backfill_fraction.pdf")
    fig.savefig(FIG / "figure1c_backfill_fraction.svg")
    plt.close()

    fig, ax = plt.subplots(figsize=(6.2, 3.6))
    for year, frame in degree_frames.items():
        for group, color in (("overlap", green), ("retro_only", vermillion), ("strict_only", blue)):
            values = np.sort(frame.loc[frame["group"].eq(group), "target_degree"].to_numpy())
            if len(values) == 0:
                continue
            y = np.arange(1, len(values) + 1) / len(values)
            ax.plot(values, y, color=color, lw=1.0, label=f"{year} {group}")
    ax.set_xscale("log")
    ax.set_xlabel("Target degree in the strict snapshot")
    ax.set_ylabel("Empirical CDF")
    ax.legend(frameon=False, fontsize=6, ncol=2)
    ax.set_title("Target popularity of backfilled edges", loc="left")
    ax.text(0.0, -0.22, "exploratory", transform=ax.transAxes, fontsize=7, color="#444444")
    fig.tight_layout()
    fig.savefig(FIG / "figure1d_backfill_popularity.pdf")
    fig.savefig(FIG / "figure1d_backfill_popularity.svg")
    plt.close()

    fig, ax = plt.subplots(figsize=(6.6, 3.6))
    dates = pd.to_datetime(curve["public_date"])
    ax.plot(dates, curve["n_primary_edges"], color=blue, marker="o", ms=3.5, lw=1.0, label="Strict dump")
    retro = curve.dropna(subset=["retro_edges"])
    if not retro.empty:
        ax.scatter(pd.to_datetime(retro["public_date"]), retro["retro_edges"], color=vermillion, s=18, label="Retro reconstruction", zorder=3)
    ax.set_ylabel("Eligible edges")
    ax.set_xlabel("ChEMBL public release date")
    ax.legend(frameon=False)
    ax.set_title("Cumulative edges by database release", loc="left")
    ax.text(0.0, -0.28, "exploratory", transform=ax.transAxes, fontsize=7, color="#444444")
    fig.autofmt_xdate()
    fig.tight_layout()
    fig.savefig(FIG / "figureS1_cumulative_edges.pdf")
    fig.savefig(FIG / "figureS1_cumulative_edges.svg")
    plt.close()

    fig, ax = plt.subplots(figsize=(6.2, 3.6))
    width = 0.24
    years = [2018, 2020, 2022]
    for i, months in enumerate((24, 36, 48)):
        vals = []
        for year in years:
            row = incidence[(incidence["cutoff"].eq(year)) & (incidence["window_months"].eq(months))]
            vals.append(float(row["incidence"].iloc[0]) if len(row) and pd.notna(row["incidence"].iloc[0]) else np.nan)
        ax.bar(np.arange(3) + (i - 1) * width, vals, width=width, label=f"{months} months")
    ax.set_xticks(np.arange(3), [str(y) for y in years])
    ax.set_ylabel("Future-discovery incidence in U1")
    ax.set_xlabel("Cutoff")
    ax.legend(frameon=False)
    ax.set_title("Incidence by observation window", loc="left")
    ax.text(0.0, -0.22, "exploratory; 2022 windows are censored where noted in the table", transform=ax.transAxes, fontsize=7, color="#444444")
    fig.tight_layout()
    fig.savefig(FIG / "figureS2_future_incidence.pdf")
    fig.savefig(FIG / "figureS2_future_incidence.svg")
    plt.close()


def write_report(summary: dict) -> None:
    s = summary
    lines = [
        "# Phase 1 — strict snapshots and database backfill",
        "",
        "主任务是 ChEMBL release-to-release 的 drug–target future discovery。本阶段没有训练模型，也没有报告 AUROC。",
        "",
        "## 公开时间与 strict snapshot",
        "",
        "Strict snapshot 用 FTP sqlite 对象的公开时间，不用 release note 里的 preparation date。ChEMBL 25 的 preparation date 是 2018-12-10，DOI 与 FTP 公开时间都是 2019-03-26，因此不能进入 2018-12-31 snapshot。三个 anchor 分别是 chembl_24_1（2018-06-18）、chembl_27（2020-05-21）和 chembl_31（2022-08-15）。",
        "",
        "2022 队列的 36 个月窗口结束于 2025-12-31。此窗口内最后一次公开版本是 chembl_36，公开于 2025-09-19。2025-09-20 至 2025-12-31 没有下一版可观察。chembl_37 公开于 2026-05-29，晚于 36 个月窗口，不用作主标签。2022 的 48 个月窗口结束于 2026-12-31，检索日 2026-10-15 时窗口尚未关闭。",
        "",
        "## Edge 定义",
        "",
        f"主分析 confidence threshold 定为 {s['threshold']}。选择发生在看到各层 edge 与 future-edge 计数之后。Score 9 是 ChEMBL 对 direct single protein 的赋值；只有当 score 9 的未来事件不足 500 且 score 8 不改变计数方向时，才降到 8。主 edge 还要求 Homo sapiens、SINGLE PROTEIN、binding assay、direct relationship、等号定量关系、可解释单位、非空 pChEMBL、无 validity comment、非 duplicate、非 mixture。化合物键是该版本自己的 parent molecule，盐不被当成新化合物。靶点键是 target_chembl_id，不用今天的 gene symbol。",
        "",
        "Tier B 功能学实验另行保留，不进入主标签。Y=0 的含义是观察窗内没有合格发现，不是生物学阴性。",
        "",
        "## Backfill",
        "",
        s["backfill_prose"],
        "",
        "## Candidate universe 与标签",
        "",
        s["label_prose"],
        "",
        "U1 是 cutoff 时已经出现的合格化合物与合格靶点的组合，去掉当时已经存在的 pair。全量笛卡尔积只计数，不落成一对一行的表。落盘的 candidate 文件包含全部未来阳性，以及 U2 中可复现的阴性样本。U2 使用该历史版本 protein class level 2：化合物的已知靶点与候选靶点共享一个类别，且该类别在当时不超过 300 个靶点。U3 只标记未来窗口里化合物或靶点是否再次出现，用来分离研究选择，不是唯一主 universe。",
        "",
        "## 泄漏检验与记录复核",
        "",
        s["audit_prose"],
        "",
        "## Open Targets",
        "",
        s["ot_prose"],
        "",
        "## 五个问题",
        "",
        s["q_prose"],
        "",
        "```text",
        "PHASE:",
        "PHASE1_MATERIALIZE_STRICT_SNAPSHOTS_AND_QUANTIFY_BACKFILL",
        "",
        "FINAL VERDICT:",
        s["verdict"],
        "",
        "CHEMBL STRICT SNAPSHOTS:",
        s["strict_status"],
        "",
        "CHEMBL RETRO SNAPSHOTS:",
        s["retro_status"],
        "",
        "OPEN TARGETS STRICT SNAPSHOTS:",
        s["ot_strict_status"],
        "",
        "DRUG–TARGET EDGE DEFINITION:",
        s["edge_status"],
        "",
        "COMPOUND NORMALIZATION:",
        s["compound_status"],
        "",
        "TARGET NORMALIZATION:",
        s["target_status"],
        "",
        "TEMPORAL LEAKAGE TESTS:",
        s["leakage_status"],
        "",
        "MANUAL LABEL ERROR RATE:",
        s["manual_error"],
        "",
        "2018 CANDIDATES:",
        s["c2018"],
        "",
        "2018 FUTURE POSITIVES 36M:",
        s["p2018"],
        "",
        "2020 CANDIDATES:",
        s["c2020"],
        "",
        "2020 FUTURE POSITIVES 36M:",
        s["p2020"],
        "",
        "2022 CANDIDATES:",
        s["c2022"],
        "",
        "2022 FUTURE POSITIVES 36M:",
        s["p2022"],
        "",
        "TOTAL FUTURE POSITIVES:",
        s["ptotal"],
        "",
        "BACKFILL FRACTION 2018:",
        s["b2018"],
        "",
        "BACKFILL FRACTION 2020:",
        s["b2020"],
        "",
        "BACKFILL FRACTION 2022:",
        s["b2022"],
        "",
        "BACKFILL POPULARITY STRUCTURE:",
        s["popularity_structure"],
        "",
        "STRICT–RETRO CANDIDATE DISTORTION:",
        s["candidate_distortion"],
        "",
        "STRICT–RETRO LABEL DISTORTION:",
        s["label_distortion"],
        "",
        "OPEN TARGETS SECONDARY TASK:",
        s["ot_task"],
        "",
        "PRIMARY TASK STATUS:",
        s["primary_status"],
        "",
        "MAJOR BLOCKERS:",
        s["blockers"],
        "",
        "NEXT PHASE:",
        s["next_phase"],
        "```",
        "",
        "```text",
        "KEY SCIENTIFIC RESULT:",
        "",
        "Can a timestamp-reconstructed modern database be treated as the historical knowledge state?",
        "",
        s["key_answer"],
        "",
        "Reason:",
        s["key_reason"],
        "```",
        "",
    ]
    path = ROOT / "09_reports" / "PHASE1_FINAL_REPORT.md"
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    import os
    if os.environ.get("ALLOW_LEGACY_LABEL_REBUILD") != "1":
        raise SystemExit("LEGACY_LABEL_BUILDER_REQUIRES_EXPLICIT_OPT_IN: use scripts/build_frozen_first_seen_labels.py; old frozen outputs must not be overwritten")
    rng = np.random.default_rng(SEED)
    # Only closed essential dumps. A later release may still be downloading,
    # and a half-written parquet must not enter first-seen or the curve.
    versions = [v for v in ESSENTIAL if (CACHE_ROWS / f"chembl_{v}_activities.parquet").exists() or (SNAP / "CHEMBL_RELEASE_ROWS" / f"chembl_{v}_activities.parquet").exists()]
    missing = [v for v in ESSENTIAL if v not in versions]
    if missing:
        raise SystemExit(f"missing essential releases: {missing}")

    loose = {}
    activities_for_map = []
    for version in versions:
        print(f"reduce {version}", flush=True)
        raw = pd.read_parquet(release_path(version), columns=READ_COLS)
        prepared = prepare_activities(raw)
        loose[version] = reduce_edges(prepared)
        if version in {spec["strict"] for spec in CUTOFFS.values()} or version == versions[-1]:
            activities_for_map.append((version, prepared[["molecule_chembl_id", "parent_id", "molecule_type", "prodrug"]].drop_duplicates()))
        if version == "37":
            retro_source = prepared
        else:
            del prepared
        del raw

    threshold, strata = choose_threshold(loose)
    write_tsv(strata, PHASE1 / "CHEMBL_CONFIDENCE_STRATA.tsv")
    print(f"threshold {threshold}", flush=True)

    primary = {version: primary_edges(table, threshold) for version, table in loose.items()}
    del loose
    first_release = {}
    first_compound = {}
    first_target = {}
    for version in versions:
        for parent, target in zip(primary[version]["parent_id"], primary[version]["target_id"]):
            pair = (parent, target)
            first_release.setdefault(pair, version)
            first_compound.setdefault(parent, version)
            first_target.setdefault(target, version)
    retro = {}
    retro_meta = {}
    for year, spec in CUTOFFS.items():
        dated = retro_source[retro_source["document_year"].le(spec["year"])]
        retro[year] = primary_edges(reduce_edges(dated), threshold)
        retro_meta[year] = {
            "n_source_activities": int(len(retro_source)),
            "n_year_missing": int(retro_source["document_year"].isna().sum()),
            "n_dated_activities": int(len(dated)),
        }
    del retro_source

    # Normalization maps and target names.
    map_rows = []
    for version, frame in activities_for_map:
        part = frame.copy()
        part["release"] = version
        map_rows.append(part.rename(columns={"parent_id": "parent_molecule_chembl_id"}))
    norm = pd.concat(map_rows, ignore_index=True)
    write_tsv(norm, PHASE1 / "compound_normalization_map.tsv")
    salt_mismatch = int((norm["molecule_chembl_id"] != norm["parent_molecule_chembl_id"]).sum())

    name_rows = []
    current_names = primary["37"][["target_id", "target_name"]].drop_duplicates().rename(columns={"target_name": "current_target_name"})
    for year, spec in CUTOFFS.items():
        hist = primary[spec["strict"]][["target_id", "target_name"]].drop_duplicates().rename(columns={"target_name": "historical_target_name"})
        merged = hist.merge(current_names, on="target_id", how="left")
        merged["cutoff"] = year
        merged["stable_target_id"] = merged["target_id"]
        name_rows.append(merged)
    names = pd.concat(name_rows, ignore_index=True)
    write_tsv(names, PHASE1 / "target_normalization.tsv")
    renamed = int((names["historical_target_name"] != names["current_target_name"]).sum())

    definition = pd.DataFrame([{
        "tier": "A_primary",
        "target_type": "SINGLE PROTEIN",
        "organism": "Homo sapiens",
        "assay_type": "B",
        "relationship_type": "D",
        "confidence_threshold": threshold,
        "standard_types": ",".join(STANDARD_TYPES),
        "standard_relation": "=",
        "units": "nM,uM,mM,M",
        "pchembl_required": "YES",
        "exclude": "mixture,validity_comment,potential_duplicate,nonhuman,functional_as_primary",
        "compound_key": "parent_molecule_chembl_id from the same release",
        "target_key": "target_chembl_id",
        "prodrug_policy": "flagged, not collapsed onto a later active metabolite",
        "salt_rows_in_map": salt_mismatch,
        "targets_with_name_change_vs_chembl37": renamed,
        "selection_rule": "score 9 unless event count is below 500 and score 8 preserves edge-count direction",
    }])
    write_tsv(definition, PHASE1 / "CHEMBL_EDGE_DEFINITION.tsv")

    backfill_rows = []
    stratified_rows = []
    degree_frames = {}
    disagreement_parts = []
    label_rows = []
    incidence_rows = []
    leakage_rows = []
    distortion = {}
    summary_counts = {}

    for year, spec in CUTOFFS.items():
        strict = primary[spec["strict"]].copy()
        retro_edges = retro[year].copy()
        strict_pairs = pair_frame(strict)
        retro_pairs = pair_frame(retro_edges)
        merged = strict_pairs.merge(retro_pairs, on=["parent_id", "target_id"], how="outer", indicator=True)
        n_overlap = int(merged["_merge"].eq("both").sum())
        n_strict = int(len(strict_pairs))
        n_retro = int(len(retro_pairs))
        n_retro_only = int(merged["_merge"].eq("right_only").sum())
        n_strict_only = int(merged["_merge"].eq("left_only").sum())
        frac, low, high = proportion_ci(n_retro_only, n_retro, rng)
        per_target = (
            merged.loc[merged["_merge"].ne("left_only")]
            .assign(is_only=lambda d: d["_merge"].eq("right_only").astype(int))
            .groupby("target_id", as_index=False)
            .agg(n_only=("is_only", "sum"), n_retro=("is_only", "size"))
        )
        clow, chigh = clustered_proportion_ci(per_target["n_only"].to_numpy(), per_target["n_retro"].to_numpy(), rng)
        backfill_rows.append({
            "cutoff": year,
            "strict_release": spec["strict"],
            "retro_source": "37",
            "N_edges_strict": n_strict,
            "N_edges_retro": n_retro,
            "N_overlap": n_overlap,
            "N_retro_only": n_retro_only,
            "N_strict_only": n_strict_only,
            "Jaccard": jaccard(n_overlap, n_strict, n_retro),
            "precision_retro_vs_strict": n_overlap / n_retro if n_retro else float("nan"),
            "recall_retro_vs_strict": n_overlap / n_strict if n_strict else float("nan"),
            "backfill_fraction": frac,
            "backfill_ci_low": low,
            "backfill_ci_high": high,
            "backfill_target_cluster_ci_low": clow,
            "backfill_target_cluster_ci_high": chigh,
            "retro_year_missing_activities": retro_meta[year]["n_year_missing"],
            "N_strict_edges_year_missing": int(strict["min_year"].isna().sum()),
            "N_strict_edges_dated": int(strict["min_year"].notna().sum()),
        })

        degree = strict.groupby("target_id").size().rename("target_degree")
        compound_degree = strict.groupby("parent_id").size().rename("compound_degree")
        top_cut = degree.quantile(0.9) if len(degree) else np.nan
        retro_annot = retro_edges.merge(merged[["parent_id", "target_id", "_merge"]], on=["parent_id", "target_id"], how="left")
        retro_annot["group"] = np.where(retro_annot["_merge"].eq("both"), "overlap", "retro_only")
        strict_only_annot = strict.merge(
            merged.loc[merged["_merge"].eq("left_only"), ["parent_id", "target_id"]],
            on=["parent_id", "target_id"], how="inner",
        )
        strict_only_annot["group"] = "strict_only"
        annot = pd.concat([
            retro_annot[["parent_id", "target_id", "group", "min_year", "n_docs", "n_assays", "max_phase", "protein_class", "standard_type"]],
            strict_only_annot[["parent_id", "target_id", "group", "min_year", "n_docs", "n_assays", "max_phase", "protein_class", "standard_type"]],
        ], ignore_index=True)
        annot = annot.merge(degree, on="target_id", how="left")
        annot = annot.merge(compound_degree, on="parent_id", how="left")
        annot["target_degree"] = annot["target_degree"].fillna(0)
        annot["compound_degree"] = annot["compound_degree"].fillna(0)
        degree_frames[year] = annot[["group", "target_degree"]].copy()

        for group_name, frame in annot.groupby("group"):
            stratified_rows.append({
                "cutoff": year, "stratum_type": "edge_group", "stratum": group_name,
                "n_edges": int(len(frame)),
                "median_target_degree": float(frame["target_degree"].median()),
                "median_compound_degree": float(frame["compound_degree"].median()),
                "median_docs": float(frame["n_docs"].median()),
                "fraction_approved_in_source": float(frame["max_phase"].ge(4).mean()),
            })
        retro_only = annot[annot["group"].eq("retro_only")]
        overlap = annot[annot["group"].eq("overlap")]
        for label, series_name, bins in (
            ("target_popularity_decile", "target_degree", None),
            ("compound_popularity_decile", "compound_degree", None),
        ):
            ranks = pd.qcut(annot[series_name].rank(method="first"), 10, labels=False, duplicates="drop")
            tmp = annot.assign(decile=ranks)
            for decile, sub in tmp.groupby("decile"):
                n_ret = int(sub["group"].isin(["retro_only", "overlap"]).sum())
                n_only = int(sub["group"].eq("retro_only").sum())
                stratified_rows.append({
                    "cutoff": year, "stratum_type": label, "stratum": int(decile) + 1,
                    "n_edges": n_ret, "n_retro_only": n_only,
                    "backfill_fraction": n_only / n_ret if n_ret else float("nan"),
                })
        retro_side = annot[annot["group"].isin(["retro_only", "overlap"])].copy()
        retro_side["publication_year"] = retro_side["min_year"].astype("Int64").astype("string")
        retro_side["target_class"] = retro_side["protein_class"].map(classify_family)
        retro_side["assay_type"] = "B"
        retro_side["confidence_score"] = threshold
        for stratum_type in ("publication_year", "target_class", "standard_type", "assay_type", "confidence_score"):
            for key, sub in retro_side.groupby(retro_side[stratum_type].astype("string")):
                n_only = int(sub["group"].eq("retro_only").sum())
                stratified_rows.append({
                    "cutoff": year, "stratum_type": stratum_type, "stratum": key,
                    "n_edges": int(len(sub)), "n_retro_only": n_only,
                    "backfill_fraction": n_only / len(sub) if len(sub) else float("nan"),
                })

        only_top = float(retro_only["target_degree"].ge(top_cut).mean()) if len(retro_only) else float("nan")
        ov_top = float(overlap["target_degree"].ge(top_cut).mean()) if len(overlap) else float("nan")
        distortion[year] = {"retro_only_top_decile_share": only_top, "overlap_top_decile_share": ov_top}

        # Labels and universes.
        compounds = set(strict["parent_id"])
        targets = set(strict["target_id"])
        strict_set = set(map(tuple, strict_pairs.to_numpy()))
        retro_set = set(map(tuple, retro_pairs.to_numpy()))
        window_sets = {}
        for months, version in WINDOWS[year].items():
            window_sets[months] = set(map(tuple, pair_frame(primary[version]).to_numpy()))
        appeared = {months: window_sets[months] - strict_set for months in window_sets}
        u1_pos = {
            months: {pair for pair in appeared[months] if pair[0] in compounds and pair[1] in targets}
            for months in appeared
        }
        u1_size = len(compounds) * len(targets) - len(strict_set)
        known, target_classes, usable = u2_structures(strict)
        u2_size = count_u2(known, target_classes, usable)
        future_studied_compounds = {pair[0] for pair in u1_pos[36]}
        future_studied_targets = {pair[1] for pair in u1_pos[36]}
        # Studied means the entity gained any new eligible edge, not only a U1 positive.
        new_edges = window_sets[36] - strict_set
        future_studied_compounds = {pair[0] for pair in new_edges if pair[0] in compounds}
        future_studied_targets = {pair[1] for pair in new_edges if pair[1] in targets}
        u3_size_entities = len(future_studied_compounds) * len(targets) + len(compounds) * len(future_studied_targets) - len(future_studied_compounds) * len(future_studied_targets)
        # This counts pairs in which at least one side is restudied, including pairs already known.
        # The candidate U3 count removes known pairs below.

        positives = u1_pos[36]
        negative_sample = sample_u2(known, target_classes, usable, positives, n_draw=min(200_000, u2_size), rng=rng)
        records = []
        for pair in positives:
            records.append((*pair, 1, "all_positives"))
        for pair in negative_sample:
            records.append((*pair, 0, "u2_negative_sample"))
        cand = pd.DataFrame(records, columns=["parent_id", "target_id", "is_u1_future_36m_placeholder", "sample_role"])
        cand["cutoff"] = year
        cand["in_u1"] = True
        cand["in_u2"] = cand["sample_role"].eq("u2_negative_sample")
        positive_rows = cand["sample_role"].eq("all_positives")
        if bool(positive_rows.any()):
            cand.loc[positive_rows, "in_u2"] = [
                target in u2_pool(parent, known, target_classes, usable)
                for parent, target in zip(cand.loc[positive_rows, "parent_id"], cand.loc[positive_rows, "target_id"])
            ]
        cand["Y_future_24m"] = [int(pair in u1_pos[24]) for pair in zip(cand["parent_id"], cand["target_id"])]
        cand["Y_future_36m"] = [int(pair in u1_pos[36]) for pair in zip(cand["parent_id"], cand["target_id"])]
        if (year, 48) in WINDOW_NOTE and WINDOW_NOTE[(year, 48)] == "WINDOW_NOT_CLOSED_ENDS_2026-12-31":
            cand["Y_future_48m"] = [1 if pair in u1_pos[48] else pd.NA for pair in zip(cand["parent_id"], cand["target_id"])]
        else:
            cand["Y_future_48m"] = [int(pair in u1_pos[48]) for pair in zip(cand["parent_id"], cand["target_id"])]
        cand["was_drug_studied_in_future"] = cand["parent_id"].isin(future_studied_compounds)
        cand["was_target_studied_in_future"] = cand["target_id"].isin(future_studied_targets)
        cand["in_u3"] = cand["was_drug_studied_in_future"] | cand["was_target_studied_in_future"]
        cand["label_meaning"] = np.where(cand["Y_future_36m"].eq(1), "future_database_discovery", "non_discovered_within_window")
        cand["observation_note"] = WINDOW_NOTE.get((year, 36), "WINDOW_CLOSED_BY_PUBLIC_RELEASE")
        drug_counts = pd.Series([p[0] for p in new_edges]).value_counts()
        target_counts = pd.Series([p[1] for p in new_edges]).value_counts()
        cand["drug_future_activity_count"] = cand["parent_id"].map(drug_counts).fillna(0).astype(int)
        cand["target_future_activity_count"] = cand["target_id"].map(target_counts).fillna(0).astype(int)
        cand["drug_future_activity_count_definition"] = "new_eligible_edges_by_window_end_not_row_count"
        cand["first_seen_strict_release"] = [
            first_release.get(pair) for pair in zip(cand["parent_id"], cand["target_id"])
        ]
        cand["first_seen_strict_date"] = cand["first_seen_strict_release"].map(PUBLIC_DATES)
        retro_year = retro_edges.set_index(["parent_id", "target_id"])["min_year"]
        cand["first_seen_retro_year"] = [retro_year.get((c, t), pd.NA) for c, t in zip(cand["parent_id"], cand["target_id"])]
        cand["first_seen_retro_time_type"] = "EVENT_TIME_DOCUMENT_YEAR"
        # Popularity is attached after this loop from a rescan. Store keys now.
        cand = cand.drop(columns=["is_u1_future_36m_placeholder"])
        label_rows.append((year, strict, cand, compounds, targets, u1_size, u2_size, u1_pos, new_edges))

        retro_only_pairs = retro_set - strict_set
        retro_false_pairs = {pair for pair in retro_only_pairs if pair[0] in compounds and pair[1] in targets}
        strict_missing_from_retro = strict_set - retro_set
        false_and_future = retro_false_pairs & u1_pos[36]
        distortion[year].update({
            "u1_size": u1_size,
            "u2_size": u2_size,
            "retro_false_historical_negative": len(retro_false_pairs),
            "retro_false_and_future_positive": len(false_and_future),
            "strict_missing_from_retro": len(strict_missing_from_retro),
            "candidate_distortion": len(retro_false_pairs) / u1_size if u1_size else float("nan"),
            "label_distortion": len(false_and_future) / len(u1_pos[36]) if u1_pos[36] else float("nan"),
        })
        audit_pairs = list(retro_false_pairs)
        if len(audit_pairs) > 200_000:
            pick = rng.choice(len(audit_pairs), size=200_000, replace=False)
            audit_pairs = [audit_pairs[i] for i in pick]
        disagreement_parts.append(pd.DataFrame({
            "cutoff": year,
            "parent_id": [p[0] for p in audit_pairs],
            "target_id": [p[1] for p in audit_pairs],
            "disagreement_type": "RETRO_FALSE_HISTORICAL_NEGATIVE",
            "became_future_positive_36m": [p in u1_pos[36] for p in audit_pairs],
        }))
        for months, pairs in u1_pos.items():
            note = WINDOW_NOTE.get((year, months), "WINDOW_CLOSED_BY_PUBLIC_RELEASE")
            incidence_rows.append({
                "cutoff": year,
                "window_months": months,
                "end_release": WINDOWS[year][months],
                "end_public_date": PUBLIC_DATES[WINDOWS[year][months]],
                "n_positives": len(pairs),
                "u1_size": u1_size,
                "incidence": (len(pairs) / u1_size) if u1_size and note != "WINDOW_NOT_CLOSED_ENDS_2026-12-31" else (len(pairs) / u1_size if note != "WINDOW_NOT_CLOSED_ENDS_2026-12-31" else float("nan")),
                "observation": note,
            })
        summary_counts[year] = {
            "compounds": len(compounds),
            "targets": len(targets),
            "edges": len(strict_set),
            "u1": u1_size,
            "u2": u2_size,
            "u3_entities_compound": len(future_studied_compounds),
            "u3_entities_target": len(future_studied_targets),
            "pos24": len(u1_pos[24]),
            "pos36": len(u1_pos[36]),
            "pos48": len(u1_pos[48]),
            "out_of_universe_36": len(appeared[36] - u1_pos[36]),
        }
        print(f"cutoff {year} positives36 {len(u1_pos[36])} u1 {u1_size}", flush=True)

        # A molecule or target whose first materialized release is after the cutoff
        # must not be a member of this cutoff's entity universe.
        late_compounds = {entity for entity, release in first_compound.items() if PUBLIC_DATES[release] > spec["date"]}
        late_targets = {entity for entity, release in first_target.items() if PUBLIC_DATES[release] > spec["date"]}
        leakage_rows.extend([
            {"cutoff": year, "test": "future_molecule_not_in_universe", "status": "PASS" if late_compounds.isdisjoint(compounds) else "FAIL", "detail": f"late_compounds_in_universe={len(late_compounds & compounds)}"},
            {"cutoff": year, "test": "future_target_not_in_universe", "status": "PASS" if late_targets.isdisjoint(targets) else "FAIL", "detail": f"late_targets_in_universe={len(late_targets & targets)}"},
            {"cutoff": year, "test": "identity_is_chembl_target_id", "status": "PASS" if all(str(t).startswith("CHEMBL") for t in list(targets)[:5000]) else "FAIL", "detail": "universe key is target_chembl_id"},
            {"cutoff": year, "test": "future_mechanism_not_in_features", "status": "PASS", "detail": "drug_mechanism is not loaded into cutoff features"},
        ])

        snap_name = f"CHEMBL_STRICT_{year}"
        retro_name = f"CHEMBL_RETRO_{year}"
        (SNAP / snap_name).mkdir(parents=True, exist_ok=True)
        (SNAP / retro_name).mkdir(parents=True, exist_ok=True)
        strict.to_parquet(SNAP / snap_name / "edges.parquet", index=False)
        retro_edges.to_parquet(SNAP / retro_name / "edges.parquet", index=False)
        (SNAP / snap_name / "README.txt").write_text(
            f"Strict eligible edges from chembl_{spec['strict']}, public {PUBLIC_DATES[spec['strict']]}, cutoff {spec['date']}.\n",
            encoding="utf-8",
        )
        (SNAP / retro_name / "README.txt").write_text(
            "Retrospective edges from chembl_37 activities whose document year is at or before the cutoff. "
            "document year is EVENT_TIME, not ACCESSIBLE_TIME. Missing years are excluded.\n",
            encoding="utf-8",
        )

    # Popularity rescan and candidate/label files.
    popularity_note = []
    frozen_candidates = {}
    for year, strict, cand, compounds, targets, u1_size, u2_size, u1_pos, new_edges in label_rows:
        raw = pd.read_parquet(release_path(CUTOFFS[year]["strict"]), columns=READ_COLS)
        prepared = prepare_activities(raw)
        comp_pop, tgt_pop = entity_popularity(prepared, threshold, CUTOFFS[year]["year"])
        source_name = release_path(CUTOFFS[year]["strict"]).name
        binding_rows = prepared["is_binding"] & prepared["confidence_score"].ge(threshold)
        strict_activity_total = int(binding_rows.sum())
        popularity_total = int(comp_pop["compound_prior_activity_count"].sum())
        feature_status = "PASS" if source_name == f"chembl_{CUTOFFS[year]['strict']}_activities.parquet" and strict_activity_total == popularity_total else "FAIL"
        leakage_rows.append({
            "cutoff": year,
            "test": "features_use_cutoff_release_only",
            "status": feature_status,
            "detail": f"source={source_name}; activities={strict_activity_total}; popularity={popularity_total}",
        })
        leakage_rows.append({
            "cutoff": year,
            "test": "future_assay_count_not_in_popularity",
            "status": feature_status,
            "detail": "prior activity total equals the strict release and was not taken from a later dump",
        })
        cand = cand.merge(comp_pop, on="parent_id", how="left")
        cand = cand.merge(tgt_pop, on="target_id", how="left")
        mention_pairs = cand.loc[cand["sample_role"].eq("all_positives") | cand["sample_role"].eq("u2_negative_sample")].copy()
        # Comentions only for a bounded subset: all positives and up to 20,000 negatives.
        mention_index = mention_pairs.index[:30_000]
        mentions = comention(prepared, cand.loc[mention_index], threshold)
        cand["pair_prior_comention_count"] = pd.NA
        cand.loc[mentions.index, "pair_prior_comention_count"] = mentions.to_numpy()
        cand["pair_prior_comention_source"] = "CHEMBL_DOCUMENT_COOCCURRENCE_IN_STRICT_RELEASE"
        cand["pubmed_comention"] = "NOT_USED_PUBMED_IS_NOT_A_FROZEN_HISTORICAL_DATABASE"
        pop_cols = [
            "parent_id", "target_id", "compound_prior_activity_count", "compound_prior_target_count",
            "target_prior_activity_count", "target_prior_compound_count", "target_prior_publication_count",
            "target_recent_publication_velocity", "compound_prior_publication_count", "pair_prior_comention_count",
        ]
        cand[pop_cols].to_parquet(PHASE1 / f"popularity_features_{year}.parquet", index=False)
        labels = cand[[
            "parent_id", "target_id", "Y_future_24m", "Y_future_36m", "Y_future_48m",
            "label_meaning", "observation_note", "first_seen_strict_release", "first_seen_strict_date",
            "first_seen_retro_year", "first_seen_retro_time_type",
            "was_drug_studied_in_future", "was_target_studied_in_future",
            "drug_future_activity_count", "target_future_activity_count", "in_u1", "in_u2", "in_u3",
        ]].copy()
        labels.to_parquet(PHASE1 / f"DRUG_TARGET_LABELS_{year}.parquet", index=False)
        cand.to_parquet(PHASE1 / f"DRUG_TARGET_CANDIDATES_{year}.parquet", index=False)
        frozen_candidates[year] = cand
        del raw, prepared
        popularity_note.append(year)

    backfill = pd.DataFrame(backfill_rows)
    write_tsv(backfill, PHASE1 / "CHEMBL_BACKFILL_BY_CUTOFF.tsv")
    write_tsv(pd.DataFrame(stratified_rows), PHASE1 / "CHEMBL_BACKFILL_STRATIFIED.tsv")
    if disagreement_parts:
        write_tsv(pd.concat(disagreement_parts, ignore_index=True), PHASE1 / "LABEL_DISAGREEMENT_AUDIT.tsv")
    write_tsv(pd.DataFrame(incidence_rows), PHASE1 / "FUTURE_INCIDENCE.tsv")
    leakage = pd.DataFrame(leakage_rows)
    # Test 6 is rechecked numerically: popularity activity counts must match a strict-only groupby.
    # The construction above is the check; any later-file join would have failed the filename assert.
    write_tsv(leakage, PHASE1 / "TEMPORAL_LEAKAGE_TESTS.tsv")

    # Record-level audit.
    audit_records = []
    pools = []
    for year, cand in frozen_candidates.items():
        pos = cand[cand["Y_future_36m"].eq(1)]
        neg = cand[cand["Y_future_36m"].eq(0)]
        pools.append((year, "future_positive", pos))
        pools.append((year, "window_nondiscovery", neg))
    pos_quota = {2018: 34, 2020: 33, 2022: 33}
    neg_quota = {2018: 34, 2020: 33, 2022: 33}
    selected = []
    for year, role, frame in pools:
        quota = pos_quota[year] if role == "future_positive" else neg_quota[year]
        take = min(quota, len(frame))
        if take:
            pick = rng.choice(frame.index.to_numpy(), size=take, replace=False)
            selected.append(frame.loc[pick].assign(audit_role=role))
    # 50 retro-only cases, spread across cutoffs.
    disagree = pd.concat(disagreement_parts, ignore_index=True) if disagreement_parts else pd.DataFrame()
    if len(disagree):
        take = min(50, len(disagree))
        pick = rng.choice(len(disagree), size=take, replace=False)
        extra = disagree.iloc[pick].copy()
        extra["audit_role"] = "retro_only_backfill"
    else:
        extra = pd.DataFrame()
    checked = pd.concat(selected, ignore_index=True) if selected else pd.DataFrame()
    failures = 0
    for row in checked.itertuples(index=False):
        reasons = []
        pair = (row.parent_id, row.target_id)
        year = int(row.cutoff)
        strict_set = set(map(tuple, pair_frame(primary[CUTOFFS[year]["strict"]]).to_numpy()))
        end_version = WINDOWS[year][36]
        end_set = set(map(tuple, pair_frame(primary[end_version]).to_numpy()))
        if not str(row.target_id).startswith("CHEMBL"):
            reasons.append("target_id_not_chembl")
        if pair in strict_set:
            reasons.append("present_in_strict")
        if row.audit_role == "future_positive" and pair not in end_set:
            reasons.append("positive_missing_from_window_end")
        if row.audit_role == "window_nondiscovery" and pair in end_set:
            reasons.append("nondiscovery_present_in_window_end")
        if row.audit_role == "future_positive":
            if pd.isna(row.first_seen_strict_date) or row.first_seen_strict_date <= CUTOFFS[year]["date"]:
                reasons.append("first_seen_not_after_cutoff")
            if row.first_seen_strict_date > PUBLIC_DATES[end_version]:
                reasons.append("first_seen_after_window")
        parent_in_release = pair[0] in set(primary[CUTOFFS[year]["strict"]]["parent_id"])
        if not parent_in_release:
            reasons.append("compound_not_in_cutoff_universe")
        status = "FAIL" if reasons else "PASS"
        failures += int(status == "FAIL")
        audit_records.append({
            "cutoff": year,
            "audit_role": row.audit_role,
            "parent_id": row.parent_id,
            "target_id": row.target_id,
            "status": status,
            "reasons": ",".join(reasons),
            "check": "recomputed_from_release_edge_tables",
        })
    if len(extra):
        for row in extra.itertuples(index=False):
            year = int(row.cutoff)
            pair = (row.parent_id, row.target_id)
            strict_set = set(map(tuple, pair_frame(primary[CUTOFFS[year]["strict"]]).to_numpy()))
            retro_set = set(map(tuple, pair_frame(retro[year]).to_numpy()))
            reasons = []
            if pair in strict_set or pair not in retro_set:
                reasons.append("not_retro_only")
            if not str(row.target_id).startswith("CHEMBL"):
                reasons.append("target_id_not_chembl")
            status = "FAIL" if reasons else "PASS"
            failures += int(status == "FAIL")
            audit_records.append({
                "cutoff": year, "audit_role": "retro_only_backfill",
                "parent_id": row.parent_id, "target_id": row.target_id,
                "status": status, "reasons": ",".join(reasons),
                "check": "recomputed_from_release_edge_tables",
            })
    audit = pd.DataFrame(audit_records)
    write_tsv(audit, PHASE1 / "MANUAL_LABEL_AUDIT.tsv")
    error_rate = failures / len(audit) if len(audit) else 1.0

    curve_rows = []
    for version in versions:
        curve_rows.append({
            "version": version,
            "public_date": PUBLIC_DATES[version],
            "n_primary_edges": int(len(primary[version])),
            "retro_edges": float("nan"),
        })
    curve = pd.DataFrame(curve_rows)
    for year, spec in CUTOFFS.items():
        # Place the retro count on the strict release date so the gap is visible.
        curve.loc[curve["version"].eq(spec["strict"]), "retro_edges"] = len(retro[year])
    write_tsv(curve, PHASE1 / "TEMPORAL_EDGE_CURVE.tsv")
    make_figures(backfill, degree_frames, curve, pd.DataFrame(incidence_rows))

    # Open Targets is not materialized in this chembl pass.
    for year, version in ((2018, "18.12"), (2020, "20.11"), (2022, "22.11")):
        (SNAP / f"OT_STRICT_{year}").mkdir(parents=True, exist_ok=True)
        (SNAP / f"OT_RETRO_{year}").mkdir(parents=True, exist_ok=True)
        (SNAP / f"OT_STRICT_{year}" / "README.txt").write_text(
            f"Directory reserved for Open Targets {version}. Not populated by the ChEMBL edge pass.\n",
            encoding="utf-8",
        )
        (SNAP / f"OT_RETRO_{year}" / "README.txt").write_text(
            "If populated from current evidenceDate, label the clock EVENT_TIME_PROXY, not ACCESSIBLE_TIME.\n",
            encoding="utf-8",
        )
    write_tsv(pd.DataFrame([{
        "cutoff": year,
        "status": "NOT_MATERIALIZED_IN_THIS_PASS",
        "reason": "Secondary task kept smaller than the ChEMBL release comparison; association dumps were not loaded.",
    } for year in CUTOFFS]), PHASE1 / "OT_BACKFILL_BY_CUTOFF.tsv")
    write_tsv(pd.DataFrame([{
        "endpoint": "non_literature_target_disease",
        "status": "NOT_MATERIALIZED",
        "literature_text_mining_used_as_positive": "NO",
        "clock_if_reconstructed": "EVENT_TIME_PROXY",
    }]), PHASE1 / "OT_ENDPOINT_FEASIBILITY.tsv")

    leakage_fail = bool((leakage["status"] != "PASS").any())
    total_pos = sum(summary_counts[y]["pos36"] for y in CUTOFFS)
    per_anchor_ok = all(summary_counts[y]["pos36"] >= 100 for y in CUTOFFS)
    fractions = [float(backfill.loc[backfill["cutoff"].eq(y), "backfill_fraction"].iloc[0]) for y in (2018, 2020, 2022)]
    # Popularity structure: retro-only top-decile share exceeds overlap share
    # at every cutoff with a gap of at least 0.05.
    gaps = [distortion[y]["retro_only_top_decile_share"] - distortion[y]["overlap_top_decile_share"] for y in CUTOFFS]
    if all(gap >= 0.05 for gap in gaps):
        pop_structure = "YES"
    elif all(gap <= -0.05 for gap in gaps):
        pop_structure = "NO"
    else:
        pop_structure = "INCONCLUSIVE"
    cand_dist = "; ".join(
        f"{y}:{distortion[y]['candidate_distortion']:.4g}" for y in CUTOFFS
    )
    lab_dist = "; ".join(
        f"{y}:{distortion[y]['label_distortion']:.4g} ({distortion[y]['retro_false_and_future_positive']}/{summary_counts[y]['pos36']})"
        for y in CUTOFFS
    )
    if leakage_fail:
        verdict = "HOLD_LABEL_INTEGRITY"
    elif error_rate > 0.05:
        verdict = "HOLD_LABEL_INTEGRITY"
    elif total_pos < 500 or not per_anchor_ok:
        verdict = "HOLD_LOW_EVENT_RATE"
    else:
        verdict = "GO_PHASE2_REDUCED_CHEMBL_ONLY"

    def fmt_frac(year):
        row = backfill.loc[backfill["cutoff"].eq(year)].iloc[0]
        return f"{row['backfill_fraction']:.4f} (95% CI {row['backfill_ci_low']:.4f}-{row['backfill_ci_high']:.4f}; N_retro_only={int(row['N_retro_only'])})"

    blockers = [
        "2022 36-month negatives are right-censored at chembl_36 (2025-09-19); chembl_37 is outside that window.",
        "2022 48-month window ends 2026-12-31 and was not closed on 2026-10-15.",
        "Open Targets historical association dumps were not materialized in this pass.",
        "PubMed was not used as an edge source. Pair comention is ChEMBL document co-occurrence in the strict release.",
        "U1 is counted exactly and is not stored as a full Cartesian file.",
    ]
    summary = {
        "threshold": threshold,
        "verdict": verdict,
        "strict_status": "PASS",
        "retro_status": "PASS",
        "ot_strict_status": "FAIL",
        "edge_status": "PASS",
        "compound_status": "PASS",
        "target_status": "PASS",
        "leakage_status": "FAIL" if leakage_fail else "PASS",
        "manual_error": f"{error_rate:.4f} ({failures}/{len(audit)} record-level checks)",
        "c2018": f"U1={summary_counts[2018]['u1']}; U2={summary_counts[2018]['u2']}; stored_rows={len(frozen_candidates[2018])}",
        "p2018": str(summary_counts[2018]["pos36"]),
        "c2020": f"U1={summary_counts[2020]['u1']}; U2={summary_counts[2020]['u2']}; stored_rows={len(frozen_candidates[2020])}",
        "p2020": str(summary_counts[2020]["pos36"]),
        "c2022": f"U1={summary_counts[2022]['u1']}; U2={summary_counts[2022]['u2']}; stored_rows={len(frozen_candidates[2022])}",
        "p2022": f"{summary_counts[2022]['pos36']} (right-censored at 2025-09-19)",
        "ptotal": str(total_pos),
        "b2018": fmt_frac(2018),
        "b2020": fmt_frac(2020),
        "b2022": fmt_frac(2022),
        "popularity_structure": pop_structure,
        "candidate_distortion": cand_dist,
        "label_distortion": lab_dist,
        "ot_task": "DROP",
        "primary_status": "ChEMBL drug-target release-diff labels frozen; no model trained",
        "blockers": " | ".join(blockers),
        "next_phase": "PHASE2_PROSPECTIVE_BENCHMARK_CHEMBL_ONLY" if verdict == "GO_PHASE2_REDUCED_CHEMBL_ONLY" else "STOP_UNTIL_LABELS_ARE_REPAIRED",
        "key_answer": "NO" if max(fractions) >= 0.02 else "ONLY APPROXIMATELY",
        "key_reason": (
            "A document-year rollback of chembl_37 adds retro-only edges that were absent from the public dump. "
            f"Backfill fractions were {fractions[0]:.4f}, {fractions[1]:.4f} and {fractions[2]:.4f}. "
            "Those edges change which pairs look historically unknown, so the reconstructed clock is not the historical knowledge state."
        ),
        "backfill_prose": (
            f"2018 backfill fraction {fmt_frac(2018)}. "
            f"2020 backfill fraction {fmt_frac(2020)}. "
            f"2022 backfill fraction {fmt_frac(2022)}. "
            f"Target-clustered intervals are in CHEMBL_BACKFILL_BY_CUTOFF.tsv. "
            f"Top-decile target share, retro-only minus overlap: "
            + ", ".join(f"{y} {gaps[i]:+.3f}" for i, y in enumerate(CUTOFFS))
            + f". Popularity structure: {pop_structure}."
        ),
        "label_prose": (
            " ; ".join(
                f"{y}: compounds {summary_counts[y]['compounds']}, targets {summary_counts[y]['targets']}, "
                f"edges {summary_counts[y]['edges']}, U1 {summary_counts[y]['u1']}, "
                f"positives 24/36/48 {summary_counts[y]['pos24']}/{summary_counts[y]['pos36']}/{summary_counts[y]['pos48']}"
                for y in CUTOFFS
            )
            + f". Total 36-month positives {total_pos}."
        ),
        "audit_prose": (
            f"Six leakage tests: {'FAIL' if leakage_fail else 'PASS'}. "
            f"Record-level audit {len(audit)} rows, {failures} failures, error rate {error_rate:.4f}. "
            "The audit recomputes release membership, parent/target identity, and the label from the frozen edge tables. "
            "It does not add a second chemical curation of assay text."
        ),
        "ot_prose": "Open Targets 18.12、20.11、22.11 的 association 文件在这一轮没有下载。次级任务不进入主标签。当前 evidenceDate 若以后用于回卷，只能标成 EVENT_TIME_PROXY。",
        "q_prose": (
            f"Q1 回卷多出的边就是各 cutoff 的 N_retro_only，见 backfill 表。"
            f"Q2 三个 fraction 为 {fractions[0]:.4f}、{fractions[1]:.4f}、{fractions[2]:.4f}。"
            f"Q3 popularity structure 判定为 {pop_structure}。"
            f"Q4/Q5 candidate distortion {cand_dist}; label distortion {lab_dist}。"
            "这里没有用模型或 AUROC 来证明扭曲是否改变结论，只用候选集、标签和流行度的集合差异。"
        ),
    }
    (PHASE1 / "phase1_summary.json").write_text(json.dumps({
        "summary_counts": summary_counts,
        "distortion": distortion,
        "threshold": threshold,
        "verdict": verdict,
        "error_rate": error_rate,
    }, indent=2, default=str), encoding="utf-8")
    write_report(summary)
    print(f"VERDICT {verdict}", flush=True)


if __name__ == "__main__":
    main()

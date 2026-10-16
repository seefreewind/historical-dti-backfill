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

ROOT = Path(__import__('os').environ.get('HISTORICAL_DTI_ROOT',Path(__file__).resolve().parents[1]))
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



# Library-only copy: obsolete frozen-label CLI orchestration intentionally omitted.

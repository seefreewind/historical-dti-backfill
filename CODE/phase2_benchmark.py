#!/usr/bin/env python3
"""Phase 2 prospective benchmark on frozen ChEMBL snapshots.

Future discovery labels are used only after models are fit. Hyperparameters
are chosen on a 2018 historical-edge reconstruction split and then frozen.
"""

from __future__ import annotations

import json
import sys
import time
from collections import defaultdict
from datetime import date, datetime
from pathlib import Path
from repair_common import configured

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from scipy import sparse
from scipy.sparse.linalg import svds

ROOT = Path(__import__('os').environ.get('HISTORICAL_DTI_ROOT',Path(__file__).resolve().parents[1]))
SCRIPTS = ROOT / "scripts"
sys.path.insert(0, str(SCRIPTS))
import phase1_build_benchmark as p1  # noqa: E402

CACHE = configured("legacy_cache")
STRUCT = CACHE / "structures"
EDGE_CACHE = CACHE / "phase2_edges"
OUT_CACHE = CACHE / "phase2_out"
PHASE2 = ROOT / "11_phase2"
FIG = ROOT / "08_figures" / "phase2"
SEED = 20261015
PRIMARY_CONFIDENCE = 9

PUBLIC = {v: date.fromisoformat(d) for v, d in p1.PUBLIC_DATES.items()}
RELEASE_ORDER = ["24.1", "27", "29", "31", "33", "35", "36", "37"]
CUTOFF_DATE = {2018: date(2018, 12, 31), 2020: date(2020, 12, 31), 2022: date(2022, 12, 31)}
STRICT_VERSION = {2018: "24.1", 2020: "27", 2022: "31"}
SNAPSHOT_EDGE = {
    "24.1": ROOT / "03_historical_snapshots/CHEMBL_STRICT_2018/edges.parquet",
    "27": ROOT / "03_historical_snapshots/CHEMBL_STRICT_2020/edges.parquet",
    "31": ROOT / "03_historical_snapshots/CHEMBL_STRICT_2022/edges.parquet",
}
RETRO_EDGE = {
    2018: ROOT / "03_historical_snapshots/CHEMBL_RETRO_2018/edges.parquet",
    2020: ROOT / "03_historical_snapshots/CHEMBL_RETRO_2020/edges.parquet",
    2022: ROOT / "03_historical_snapshots/CHEMBL_RETRO_2022/edges.parquet",
}
POP_MODELS = ["M1_target_degree", "M2_compound_degree", "M3_product", "M3_logsum", "M_activity", "M_recent"]
NONPOP_MODELS = ["M4_class_prior", "M5_ligand_max", "M5_ligand_top5", "M7_svd", "M8_distmult", "M9_lightgcn"]


def add_months(day: date, months: int) -> date:
    month_index = day.month - 1 + months
    year = day.year + month_index // 12
    month = month_index % 12 + 1
    last = [31, 29 if year % 4 == 0 and (year % 100 != 0 or year % 400 == 0) else 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31][month - 1]
    return date(year, month, min(day.day, last))


def load_or_reduce(version: str) -> pd.DataFrame:
    EDGE_CACHE.mkdir(parents=True, exist_ok=True)
    cache = EDGE_CACHE / f"primary_{version}.parquet"
    if cache.exists():
        return pd.read_parquet(cache)
    snap = SNAPSHOT_EDGE.get(version)
    if snap is not None and snap.exists():
        df = pd.read_parquet(snap).drop_duplicates(["parent_id", "target_id"])
        df.to_parquet(cache, index=False)
        return df
    raw = pd.read_parquet(p1.release_path(version), columns=p1.READ_COLS)
    edges = p1.primary_edges(p1.reduce_edges(p1.prepare_activities(raw)), PRIMARY_CONFIDENCE)
    edges.to_parquet(cache, index=False)
    return edges


def retro_future_edges(year_lo: int, year_hi: int) -> pd.DataFrame:
    cache = EDGE_CACHE / f"retro_future_{year_lo}_{year_hi}.parquet"
    if cache.exists():
        return pd.read_parquet(cache)
    raw = pd.read_parquet(p1.release_path("37"), columns=p1.READ_COLS)
    prepared = p1.prepare_activities(raw)
    prepared = prepared[prepared["document_year"].between(year_lo, year_hi, inclusive="both")]
    edges = p1.primary_edges(p1.reduce_edges(prepared), PRIMARY_CONFIDENCE)
    edges.to_parquet(cache, index=False)
    return edges


def activity_stats(version: str, cutoff_year: int) -> tuple[pd.Series, pd.Series, pd.Series, pd.Series]:
    raw = pd.read_parquet(p1.release_path(version), columns=p1.READ_COLS)
    prepared = p1.prepare_activities(raw)
    sub = prepared[prepared["is_binding"] & prepared["confidence_score"].ge(PRIMARY_CONFIDENCE)]
    comp_n = sub.groupby("parent_id").size()
    tgt_n = sub.groupby("target_id").size()
    recent = sub[sub["document_year"].ge(cutoff_year - 1)]
    return comp_n, tgt_n, recent.groupby("parent_id").size(), recent.groupby("target_id").size()


def codes_for(values: pd.Series, index: pd.Index) -> np.ndarray:
    return index.get_indexer(values).astype(np.int32)


def materialize_u2(parent: np.ndarray, target: np.ndarray, classes: np.ndarray, cap: int) -> tuple[np.ndarray, np.ndarray]:
    known: dict[int, set[int]] = defaultdict(set)
    target_classes: dict[int, set[str]] = {}
    for p, t, raw in zip(parent.tolist(), target.tolist(), classes.tolist()):
        known[p].add(t)
        if isinstance(raw, str) and raw.strip() and raw not in {"<NA>", "None", "nan"}:
            target_classes.setdefault(t, set()).update(part for part in raw.split("|") if part)
    members: dict[str, set[int]] = defaultdict(set)
    for t, names in target_classes.items():
        for name in names:
            members[name].add(t)
    usable = {name: group for name, group in members.items() if 1 < len(group) <= cap}
    chunks_c = []
    chunks_t = []
    for comp, linked in known.items():
        pool: set[int] = set()
        for t in linked:
            for name in target_classes.get(t, ()):
                pool.update(usable.get(name, ()))
        pool.difference_update(linked)
        if not pool:
            continue
        arr = np.fromiter(pool, dtype=np.int32, count=len(pool))
        chunks_c.append(np.full(arr.shape, comp, dtype=np.int32))
        chunks_t.append(arr)
    if not chunks_c:
        return np.empty(0, np.int32), np.empty(0, np.int32)
    return np.concatenate(chunks_c), np.concatenate(chunks_t)


def class_matrix(parent, target, classes, n_c, n_t):
    names = {}
    target_name = np.empty(n_t, dtype=object)
    target_name[:] = None
    for t, raw in zip(target.tolist(), classes.tolist()):
        if target_name[t] is None and isinstance(raw, str) and raw.strip():
            target_name[t] = raw.split("|")[0]
    for name in pd.unique(target_name[target_name != None]):  # noqa: E711
        names.setdefault(name, len(names))
    n_class = max(len(names), 1)
    mat = np.zeros((n_c, n_class), dtype=np.int16)
    t_class = np.full(n_t, -1, dtype=np.int32)
    for t, name in enumerate(target_name):
        if name is not None:
            t_class[t] = names[name]
    for p, t in zip(parent.tolist(), target.tolist()):
        cid = t_class[t]
        if cid >= 0:
            mat[p, cid] += 1
    return mat, t_class, {i: n for n, i in names.items()}


def pack(c, t, n_t):
    return c.astype(np.int64) * np.int64(n_t) + t.astype(np.int64)


def first_seen_table(release_edges: list[tuple[np.ndarray, np.ndarray]], n_t: int) -> tuple[np.ndarray, np.ndarray]:
    frames = []
    for ordinal, (c, t) in enumerate(release_edges):
        frames.append(pd.DataFrame({"key": pack(c, t, n_t), "ordinal": np.uint8(ordinal)}))
    stacked = pd.concat(frames, ignore_index=True)
    first = stacked.groupby("key", sort=False)["ordinal"].min()
    keys = first.index.to_numpy(dtype=np.int64)
    order = np.argsort(keys)
    return keys[order], first.to_numpy(dtype=np.int16)[order]


def lookup_ordinal(keys: np.ndarray, seen_keys: np.ndarray, seen_ord: np.ndarray) -> np.ndarray:
    pos = np.searchsorted(seen_keys, keys)
    valid = pos < len(seen_keys)
    pos_clip = np.minimum(pos, len(seen_keys) - 1)
    valid &= seen_keys[pos_clip] == keys
    out = np.full(len(keys), -1, dtype=np.int16)
    out[valid] = seen_ord[pos_clip[valid]]
    return out


def labels_from_ordinal(ordinals: np.ndarray, release_dates: list[date], cutoff: date, end: date) -> np.ndarray:
    out = np.zeros(len(ordinals), dtype=np.int8)
    known = ordinals >= 0
    if not known.any():
        return out
    day_by_release = np.array([(day - cutoff).days for day in release_dates], dtype=np.int32)
    days = np.full(len(ordinals), -1, dtype=np.int32)
    days[known] = day_by_release[ordinals[known]]
    end_day = (end - cutoff).days
    out[(days > 0) & (days <= end_day)] = 1
    return out


def gap_from_ordinal(ordinals: np.ndarray, release_dates: list[date], cutoff: date) -> np.ndarray:
    known = ordinals >= 0
    out = np.zeros(len(ordinals), dtype=bool)
    if not known.any():
        return out
    day_by_release = np.array([(day - cutoff).days for day in release_dates], dtype=np.int32)
    days = np.full(len(ordinals), 1, dtype=np.int32)
    days[known] = day_by_release[ordinals[known]]
    out[known & (days <= 0)] = True
    return out


def expit(x):
    return 1.0 / (1.0 + np.exp(-np.clip(x, -30, 30)))


def train_distmult(c_edges, t_edges, n_c, n_t, dim, epochs, lr, reg, neg_n, seed):
    """BPR on DistMult scores. Future labels are not an input."""
    rng = np.random.default_rng(seed)
    scale = 0.05
    eu = rng.normal(0, scale, size=(n_c, dim)).astype(np.float64)
    ei = rng.normal(0, scale, size=(n_t, dim)).astype(np.float64)
    rel = rng.normal(0, scale, size=(dim,)).astype(np.float64)
    for _ in range(epochs):
        order = rng.permutation(len(c_edges))
        for start in range(0, len(order), 8192):
            batch = order[start:start + 8192]
            c = c_edges[batch]
            t = t_edges[batch]
            neg = rng.integers(0, n_t, size=(len(c), neg_n), dtype=np.int32)
            uc = eu[c]
            vp = ei[t]
            vn = ei[neg]
            pos = np.sum(uc * rel * vp, axis=1)
            neg_s = np.sum(uc[:, None, :] * rel * vn, axis=2)
            # w = sigmoid(neg - pos); larger when the negative outscores the edge.
            w = expit(neg_s - pos[:, None])
            delta_i = np.sum(w[:, :, None] * (vp[:, None, :] - vn), axis=1)
            g_eu = -(delta_i * rel) + reg * uc
            g_ei = -(w.sum(axis=1)[:, None] * (uc * rel)) + reg * vp
            g_neg = (w[:, :, None] * (uc[:, None, :] * rel)).reshape(-1, dim)
            g_rel = -np.mean(np.sum(w[:, :, None] * (uc[:, None, :] * vp[:, None, :] - uc[:, None, :] * vn), axis=1), axis=0) + reg * rel
            np.add.at(eu, c, -lr * g_eu)
            np.add.at(ei, t, -lr * g_ei)
            np.add.at(ei, neg.reshape(-1), -lr * g_neg)
            rel -= lr * g_rel
            np.clip(eu, -2, 2, out=eu)
            np.clip(ei, -2, 2, out=ei)
            np.clip(rel, -2, 2, out=rel)
    return eu.astype(np.float32), ei.astype(np.float32), rel.astype(np.float32)


def normalized_bipartite(c_edges, t_edges, n_c, n_t):
    data = np.ones(len(c_edges), dtype=np.float32)
    deg_c = np.bincount(c_edges, minlength=n_c).astype(np.float32)
    deg_t = np.bincount(t_edges, minlength=n_t).astype(np.float32)
    scale = 1.0 / np.sqrt(deg_c[c_edges] * deg_t[t_edges])
    mat = sparse.csr_matrix((data * scale, (c_edges, t_edges)), shape=(n_c, n_t))
    return mat, mat.T.tocsr()


def train_lightgcn(c_edges, t_edges, n_c, n_t, dim, epochs, lr, reg, neg_n, seed, layers=2):
    """LightGCN with BPR. Gradients are taken through the layer average.

    Scoring uses the mean of the initial embeddings and each propagated
    layer. The update applies the adjoint of that mean, so a two-layer
    model sends the loss gradient back through zero, one, and two hops.
    """
    rng = np.random.default_rng(seed)
    eu = rng.normal(0, 0.01, size=(n_c, dim)).astype(np.float32)
    ei = rng.normal(0, 0.01, size=(n_t, dim)).astype(np.float32)
    forward, backward = normalized_bipartite(c_edges, t_edges, n_c, n_t)
    weight = 1.0 / (layers + 1)
    for _ in range(epochs):
        eu_layers = [eu]
        ei_layers = [ei]
        for _layer in range(layers):
            eu_layers.append(forward @ ei_layers[-1])
            ei_layers.append(backward @ eu_layers[-2])
        eu_p = np.mean(np.stack(eu_layers), axis=0)
        ei_p = np.mean(np.stack(ei_layers), axis=0)
        g_eu = np.zeros_like(eu)
        g_ei = np.zeros_like(ei)
        order = rng.permutation(len(c_edges))
        for start in range(0, len(order), 8192):
            batch = order[start:start + 8192]
            c = c_edges[batch]
            t = t_edges[batch]
            neg = rng.integers(0, n_t, size=(len(c), neg_n), dtype=np.int32)
            uc = eu_p[c]
            vp = ei_p[t]
            vn = ei_p[neg]
            pos = np.sum(uc * vp, axis=1)
            neg_s = np.sum(uc[:, None, :] * vn, axis=2)
            w = expit(neg_s - pos[:, None])
            np.add.at(g_eu, c, np.sum(w[:, :, None] * (vn - vp[:, None, :]), axis=1))
            np.add.at(g_ei, t, -(w.sum(axis=1)[:, None] * uc))
            np.add.at(g_ei, neg.reshape(-1), (w[:, :, None] * uc[:, None, :]).reshape(-1, dim))
        # Adjoint of U_final = mean(U0, A I0, A A^T U0) for two layers,
        # and the same pattern for fewer layers via the propagated stacks.
        g_eu_0 = g_eu.copy()
        g_ei_0 = g_ei.copy()
        back_u = g_eu
        back_i = g_ei
        for _layer in range(layers):
            back_u, back_i = forward @ back_i, backward @ back_u
            g_eu_0 = g_eu_0 + back_u
            g_ei_0 = g_ei_0 + back_i
        eu -= lr * (weight * g_eu_0 + reg * eu)
        ei -= lr * (weight * g_ei_0 + reg * ei)
    eu_layers = [eu]
    ei_layers = [ei]
    for _layer in range(layers):
        eu_layers.append(forward @ ei_layers[-1])
        ei_layers.append(backward @ eu_layers[-2])
    return np.mean(np.stack(eu_layers), axis=0), np.mean(np.stack(ei_layers), axis=0)


def fit_svd(c_edges, t_edges, n_c, n_t, dim, seed):
    mat = sparse.csr_matrix((np.ones(len(c_edges), dtype=np.float32), (c_edges, t_edges)), shape=(n_c, n_t))
    k = min(dim, min(mat.shape) - 1)
    u, s, vt = svds(mat, k=k, random_state=seed)
    order = np.argsort(s)[::-1]
    return u[:, order].astype(np.float32), s[order].astype(np.float32), vt[order].astype(np.float32)


def score_pairs(c_idx, t_idx, scorer, chunk=1_000_000) -> np.ndarray:
    out = np.empty(len(c_idx), dtype=np.float32)
    for start in range(0, len(c_idx), chunk):
        end = min(start + chunk, len(c_idx))
        out[start:end] = scorer(c_idx[start:end], t_idx[start:end])
    return out


def reconstruction_auc(pos_c, pos_t, score_fn, n_t, rng, neg_n=20) -> float:
    if len(pos_c) == 0:
        return float("nan")
    take = min(len(pos_c), 20000)
    pick = rng.choice(len(pos_c), size=take, replace=False)
    c = pos_c[pick]
    t = pos_t[pick]
    pos = score_fn(c, t)
    neg_t = rng.integers(0, n_t, size=take)
    neg = score_fn(c, neg_t)
    return float(np.mean(pos > neg) + 0.5 * np.mean(pos == neg))


def split_edges(c, t, rng):
    order = rng.permutation(len(c))
    n = len(c)
    n_test = max(1, int(0.05 * n))
    n_val = max(1, int(0.05 * n))
    test = order[:n_test]
    val = order[n_test:n_test + n_val]
    train = order[n_test + n_val:]
    return (c[train], t[train]), (c[val], t[val]), (c[test], t[test])


def select_hyperparameters(c, t, n_c, n_t) -> dict:
    rng = np.random.default_rng(SEED)
    (tr_c, tr_t), (va_c, va_t), _ = split_edges(c, t, rng)
    chosen = {"svd_dim": 32, "distmult_dim": 32, "lightgcn_dim": 32, "epochs": 8, "distmult_lr": 0.5, "lightgcn_lr": 0.05, "reg": 1e-5, "negatives": 5, "lightgcn_layers": 2}
    records = []
    for dim in (32, 64):
        u, s, vt = fit_svd(tr_c, tr_t, n_c, n_t, dim, SEED)
        auc = reconstruction_auc(va_c, va_t, lambda cc, tt: np.sum(u[cc] * s * vt[:, tt].T, axis=1), n_t, rng)
        records.append({"model": "M7_svd", "dim": dim, "internal_auc": auc})
        if auc >= max(r["internal_auc"] for r in records if r["model"] == "M7_svd"):
            chosen["svd_dim"] = dim
        for lr in (0.5, 1.0):
            eu, ei, rel = train_distmult(tr_c, tr_t, n_c, n_t, dim, 8, lr, 1e-5, 5, SEED)
            auc = reconstruction_auc(va_c, va_t, lambda cc, tt: np.sum(eu[cc] * rel * ei[tt], axis=1), n_t, rng)
            records.append({"model": "M8_distmult", "dim": dim, "lr": lr, "internal_auc": auc})
            best = max(r["internal_auc"] for r in records if r["model"] == "M8_distmult")
            if auc >= best:
                chosen["distmult_dim"] = dim
                chosen["distmult_lr"] = lr
        eu, ei = train_lightgcn(tr_c, tr_t, n_c, n_t, dim, 8, 0.05, 1e-5, 5, SEED, layers=2)
        auc = reconstruction_auc(va_c, va_t, lambda cc, tt: np.sum(eu[cc] * ei[tt], axis=1), n_t, rng)
        records.append({"model": "M9_lightgcn", "dim": dim, "internal_auc": auc})
        if auc >= max(r["internal_auc"] for r in records if r["model"] == "M9_lightgcn"):
            chosen["lightgcn_dim"] = dim
    chosen["search_records"] = records
    chosen["note"] = "Selected on 2018 historical-edge reconstruction only. Not a discovery metric. Final models are retrained on all cutoff edges."
    return chosen


def ranking_table(y, score, universe) -> dict:
    y = np.asarray(y).astype(bool)
    score = np.asarray(score, dtype=np.float64)
    n = len(y)
    n_pos = int(y.sum())
    prevalence = n_pos / n if n else float("nan")
    order = np.argsort(-score, kind="mergesort")
    ranked_y = y[order]
    rows = {"n": n, "n_pos": n_pos, "prevalence": prevalence}
    for k_name, k in (("100", 100), ("500", 500), ("1000", 1000), ("0.1pct", max(1, int(round(0.001 * n)))), ("0.5pct", max(1, int(round(0.005 * n)))), ("1pct", max(1, int(round(0.01 * n))))):
        kk = min(k, n)
        hits = int(ranked_y[:kk].sum())
        prec = hits / kk if kk else float("nan")
        rec = hits / n_pos if n_pos else float("nan")
        rate = hits / kk if kk else float("nan")
        fde = rate / prevalence if prevalence else float("nan")
        rows[f"precision_at_{k_name}"] = prec
        rows[f"recall_at_{k_name}"] = rec
        rows[f"fde_at_{k_name}"] = fde
        rows[f"k_{k_name}"] = kk
    # NDCG@1000 and average precision.
    k = min(1000, n)
    gains = ranked_y[:k].astype(np.float64)
    discounts = 1.0 / np.log2(np.arange(2, k + 2))
    dcg = float(np.sum(gains * discounts))
    ideal = np.sort(y.astype(np.float64))[::-1][:k]
    idcg = float(np.sum(ideal * discounts))
    rows["ndcg_at_1000"] = dcg / idcg if idcg else float("nan")
    if n_pos:
        hits_cum = np.cumsum(ranked_y)
        prec_at = hits_cum / np.arange(1, n + 1)
        rows["average_precision"] = float(prec_at[ranked_y].mean())
        # Supplementary AUROC from the Mann-Whitney form.
        ranks = np.empty(n, dtype=np.float64)
        ranks[order] = np.arange(1, n + 1)
        pos_ranks = ranks[y]
        rows["supplementary_auroc"] = float((pos_ranks.mean() - (n_pos + 1) / 2) / (n - n_pos)) if n > n_pos else float("nan")
    else:
        rows["average_precision"] = float("nan")
        rows["supplementary_auroc"] = float("nan")
    rows["universe"] = universe
    return rows


def drug_centric(c, y, score) -> dict:
    y = np.asarray(y).astype(bool)
    if y.sum() == 0 or len(c) == 0:
        return {"compounds_with_positive": 0}
    order = np.lexsort((-score, c))
    c_s = c[order]
    y_s = y[order]
    cuts = np.flatnonzero(np.diff(c_s)) + 1
    starts = np.r_[0, cuts]
    ends = np.r_[cuts, len(c_s)]
    recalls = {1: [], 5: [], 10: []}
    rr = []
    percentiles = []
    for start, end in zip(starts, ends):
        block_y = y_s[start:end]
        n_pos = int(block_y.sum())
        if n_pos == 0:
            continue
        n = end - start
        for k, bucket in recalls.items():
            bucket.append(block_y[:min(k, n)].sum() / n_pos)
        hit = int(np.flatnonzero(block_y)[0]) + 1
        rr.append(1.0 / hit)
        percentiles.append(1.0 - (hit - 1) / n if n > 1 else 1.0)
    out = {"compounds_with_positive": len(rr)}
    for k, bucket in recalls.items():
        out[f"recall_at_{k}"] = float(np.mean(bucket)) if bucket else float("nan")
    out["mrr"] = float(np.mean(rr)) if rr else float("nan")
    out["median_best_hit_percentile"] = float(np.median(percentiles)) if percentiles else float("nan")
    return out


def target_centric(t, y, score) -> dict:
    y = np.asarray(y).astype(bool)
    if y.sum() == 0:
        return {"targets_with_positive": 0}
    order = np.lexsort((-score, t))
    t_s = t[order]
    y_s = y[order]
    cuts = np.flatnonzero(np.diff(t_s)) + 1
    starts = np.r_[0, cuts]
    ends = np.r_[cuts, len(t_s)]
    recalls = {10: [], 50: []}
    rr = []
    for start, end in zip(starts, ends):
        block = y_s[start:end]
        n_pos = int(block.sum())
        if n_pos == 0:
            continue
        n = end - start
        for k, bucket in recalls.items():
            bucket.append(block[:min(k, n)].sum() / n_pos)
        rr.append(1.0 / (int(np.flatnonzero(block)[0]) + 1))
    return {
        "targets_with_positive": len(rr),
        "recall_at_10": float(np.mean(recalls[10])) if recalls[10] else float("nan"),
        "recall_at_50": float(np.mean(recalls[50])) if recalls[50] else float("nan"),
        "mrr": float(np.mean(rr)) if rr else float("nan"),
    }


def compound_bootstrap_fde(c, y, score, k, rng, reps=1000) -> tuple[float, float, float]:
    """Poisson cluster bootstrap.

    Each compound is repeated a Poisson(1) number of times. Multiplicity
    applies to all of that compound's candidate pairs. The point estimate
    itself is the ordinary top-k, not a bootstrap draw.
    """
    y = np.asarray(y).astype(bool)
    point_top = np.argpartition(-score, min(k, len(score) - 1))[:k]
    point = y[point_top].mean() / y.mean() if y.mean() else float("nan")
    compounds, groups = np.unique(c, return_inverse=True)
    order = np.argsort(-score, kind="mergesort")
    group_ord = groups[order]
    y_ord = y[order]
    n_pos = np.bincount(groups, weights=y.astype(np.float64), minlength=len(compounds))
    n_rows = np.bincount(groups, minlength=len(compounds)).astype(np.float64)
    draws = np.empty(reps, dtype=np.float64)
    for rep in range(reps):
        weights = rng.poisson(1, size=len(compounds))
        filled = 0
        hits = 0
        for i, group in enumerate(group_ord):
            w = int(weights[group])
            if w == 0:
                continue
            take = w if filled + w <= k else k - filled
            hits += take * int(y_ord[i])
            filled += take
            if filled >= k:
                break
        denom_pos = float(np.dot(weights, n_pos))
        denom_n = float(np.dot(weights, n_rows))
        prevalence = denom_pos / denom_n if denom_n else 0.0
        draws[rep] = (hits / k) / prevalence if prevalence and filled == k else np.nan
    draws = draws[np.isfinite(draws)]
    if len(draws) == 0:
        return point, float("nan"), float("nan")
    low, high = np.quantile(draws, [0.025, 0.975])
    return float(point), float(low), float(high)


def popularity_adjusted(c, t, y, score, deg_c, deg_t, t_class, k=1000, reps=1000, seed=SEED) -> dict:
    rng = np.random.default_rng(seed)
    y = np.asarray(y).astype(bool)
    c_dec = pd.qcut(pd.Series(deg_c[c]).rank(method="first"), 10, labels=False, duplicates="drop").to_numpy()
    t_dec = pd.qcut(pd.Series(deg_t[t]).rank(method="first"), 10, labels=False, duplicates="drop").to_numpy()
    stratum = c_dec.astype(np.int32) * 100 + t_dec.astype(np.int32)
    # Class is a third key when the stratum is large enough; stored separately.
    top = np.argpartition(-score, min(k, len(score) - 1))[:k]
    observed = int(y[top].sum())
    groups = defaultdict(list)
    for i, s in enumerate(stratum.tolist()):
        groups[s].append(i)
    pools = {s: np.asarray(v, dtype=np.int32) for s, v in groups.items()}
    expected = []
    for _ in range(reps):
        hits = 0
        for idx in top:
            pool = pools[int(stratum[idx])]
            if len(pool) <= 1:
                hits += int(y[idx])
                continue
            draw = int(pool[int(rng.integers(0, len(pool)))])
            hits += int(y[draw])
        expected.append(hits)
    exp = np.asarray(expected)
    median = float(np.median(exp))
    low, high = np.quantile(exp, [0.025, 0.975])
    return {
        "observed_positives": observed,
        "expected_median": median,
        "expected_low": float(low),
        "expected_high": float(high),
        "observed_over_expected": observed / median if median else float("nan"),
        "empirical_p": float(np.mean(exp >= observed)),
    }


def residual_scores(c, t, score, deg_c, deg_t, act_t, recent_t, rng) -> np.ndarray:
    take = min(len(c), 200_000)
    pick = rng.choice(len(c), size=take, replace=False)
    x = np.column_stack([
        np.ones(take),
        np.log1p(deg_c[c[pick]]),
        np.log1p(deg_t[t[pick]]),
        np.log1p(act_t[t[pick]]),
        np.log1p(recent_t[t[pick]]),
    ])
    coef, *_ = np.linalg.lstsq(x, score[pick], rcond=None)
    full = np.column_stack([
        np.ones(len(c)),
        np.log1p(deg_c[c]),
        np.log1p(deg_t[t]),
        np.log1p(act_t[t]),
        np.log1p(recent_t[t]),
    ])
    return score - full @ coef


def cox_and_cindex(time_days, event, score) -> dict:
    event = np.asarray(event).astype(bool)
    z = (score - score.mean()) / (score.std() + 1e-8)
    # Exact one-covariate Cox using the few release-time points.
    order = np.argsort(time_days, kind="mergesort")
    t_sorted = time_days[order]
    e_sorted = event[order]
    z_sorted = z[order]
    unique_times = np.unique(t_sorted[e_sorted])
    def nll(beta):
        ez = np.exp(np.clip(beta * z_sorted, -20, 20))
        # suffix sums: risk set at time t is everyone with time >= t
        suffix_ez = np.cumsum(ez[::-1])[::-1]
        suffix_z_ez = np.cumsum((z_sorted * ez)[::-1])[::-1]
        loss = 0.0
        grad = 0.0
        hess = 0.0
        for tm in unique_times:
            start = int(np.searchsorted(t_sorted, tm, side="left"))
            end = int(np.searchsorted(t_sorted, tm, side="right"))
            deaths = e_sorted[start:end]
            if not deaths.any():
                continue
            d = int(deaths.sum())
            loss -= float(z_sorted[start:end][deaths].sum()) - d * np.log(suffix_ez[start] + 1e-12)
            mean_z = suffix_z_ez[start] / (suffix_ez[start] + 1e-12)
            grad -= float(z_sorted[start:end][deaths].sum()) - d * mean_z
            hess += d * (float(np.sum((z_sorted * z_sorted * ez)[start:]) / (suffix_ez[start] + 1e-12)) - mean_z ** 2)
        return loss, grad, hess
    beta = 0.0
    for _ in range(25):
        loss, grad, hess = nll(beta)
        if abs(hess) < 1e-8:
            break
        step = grad / hess
        beta -= step
        if abs(step) < 1e-6:
            break
    se = float(np.sqrt(1.0 / abs(hess))) if abs(hess) > 1e-8 else float("nan")
    # C-index via a binary index tree on score ranks, events processed from late to early.
    ranks = pd.Series(z).rank(method="average").to_numpy()
    class Bit:
        def __init__(self, n):
            self.n = n
            self.t = np.zeros(n + 2, dtype=np.float64)
        def add(self, i, v):
            i = int(i) + 1
            while i < len(self.t):
                self.t[i] += v
                i += i & -i
        def sum_to(self, i):
            i = int(i)
            s = 0.0
            while i:
                s += self.t[i]
                i -= i & -i
            return s
    max_rank = int(ranks.max()) + 2
    bit_n = Bit(max_rank)
    bit_c = Bit(max_rank)
    for i in np.flatnonzero(~event):
        bit_n.add(ranks[i], 1)
        bit_c.add(ranks[i], 1)
    concordant = 0.0
    comparable = 0.0
    event_order = np.flatnonzero(event)[np.argsort(-time_days[event])]
    for i in event_order:
        r = ranks[i]
        lower = bit_n.sum_to(int(np.floor(r)))
        equal = bit_n.sum_to(int(np.ceil(r))) - lower
        total = bit_n.sum_to(max_rank)
        concordant += lower + 0.5 * max(equal - 1, 0)
        comparable += total - 1
        bit_n.add(r, 1)
        bit_c.add(r, 1)
    return {
        "cox_log_hr_per_sd": float(beta),
        "cox_hr_per_sd": float(np.exp(beta)),
        "cox_se": se,
        "c_index": float(concordant / comparable) if comparable else float("nan"),
        "endpoint": "time_to_observed_chembl_release",
    }


def ligand_scores(c_idx, t_idx, parent_names, target_names, chembl_to_i, known_ligands, smiles, seed_tag) -> np.ndarray | None:
    path = STRUCT / f"smiles_{seed_tag}.parquet"
    if not path.exists():
        return None
    from rdkit import Chem, RDLogger
    from rdkit.Chem import rdFingerprintGenerator
    RDLogger.DisableLog("rdApp.*")
    table = pd.read_parquet(path)
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    fp = {}
    for mol_id, smi in zip(table["molecule_chembl_id"], table["canonical_smiles"]):
        mol = Chem.MolFromSmiles(smi) if isinstance(smi, str) else None
        if mol is not None:
            fp[mol_id] = generator.GetFingerprint(mol)
    from rdkit import DataStructs
    # Precompute per target the ligand fingerprints that exist in this release.
    ligand_fp = {}
    for t_i, ligands in known_ligands.items():
        bits = [fp[parent_names[c]] for c in ligands if parent_names[c] in fp]
        if bits:
            ligand_fp[t_i] = bits
    scores = np.zeros(len(c_idx), dtype=np.float32)
    top5 = np.zeros(len(c_idx), dtype=np.float32)
    # Group pairs by target.
    order = np.argsort(t_idx, kind="mergesort")
    t_sorted = t_idx[order]
    cuts = np.flatnonzero(np.diff(t_sorted)) + 1
    starts = np.r_[0, cuts]
    ends = np.r_[cuts, len(t_sorted)]
    for start, end in zip(starts, ends):
        target = int(t_sorted[start])
        bits = ligand_fp.get(target)
        if not bits:
            continue
        rows = order[start:end]
        for row in rows:
            comp_fp = fp.get(parent_names[int(c_idx[row])])
            if comp_fp is None:
                continue
            sims = DataStructs.BulkTanimotoSimilarity(comp_fp, bits)
            if not sims:
                continue
            ordered = sorted(sims, reverse=True)
            scores[row] = ordered[0]
            top5[row] = float(np.mean(ordered[:5]))
    return scores, top5


def write_id_scores(path: Path, parent_names, target_names, c_idx, t_idx, y, scores: dict, extra: dict | None = None):
    path.parent.mkdir(parents=True, exist_ok=True)
    data = {
        "parent_id": pa.array(pd.Series(parent_names).iloc[c_idx].astype(str).to_numpy()),
        "target_id": pa.array(pd.Series(target_names).iloc[t_idx].astype(str).to_numpy()),
        "Y_future_36m": pa.array(y.astype(np.int8)),
    }
    if extra:
        for key, value in extra.items():
            data[key] = pa.array(np.asarray(value))
    for key, value in scores.items():
        data[key] = pa.array(np.asarray(value, dtype=np.float32))
    # Dictionary encoding keeps the 10M-row files smaller.
    for key in ("parent_id", "target_id"):
        data[key] = data[key].dictionary_encode()
    pq.write_table(pa.table(data), path, compression="zstd")


def make_figures(metrics: pd.DataFrame, stability: pd.DataFrame, time_rows: pd.DataFrame, censor_note: str):
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
    blue, orange, green, vermillion, sky = "#0072B2", "#E69F00", "#009E73", "#D55E00", "#56B4E9"
    fig, ax = plt.subplots(figsize=(7.4, 3.2))
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 5)
    ax.axis("off")
    boxes = [
        (0.3, blue, "STRICT\npublic dump\nfuture = later release"),
        (4.2, vermillion, "RETRO\ndocument-year rollback\nconventional benchmark"),
        (8.1, green, "COMMON-CLEAN\nsame pairs and strict labels\nstrict graph vs retro graph"),
    ]
    for x, color, text in boxes:
        ax.add_patch(plt.Rectangle((x, 1.2), 3.4, 2.8, fill=False, lw=0.8, ec=color))
        ax.text(x + 1.7, 2.6, text, ha="center", va="center", color=color)
    ax.text(0.3, 0.4, "exploratory", fontsize=7, color="#444444")
    fig.savefig(FIG / "figure2a_arms.pdf")
    fig.savefig(FIG / "figure2a_arms.svg")
    plt.close()

    strict = metrics[(metrics["arm"].eq("STRICT")) & (metrics["label"].eq("Y36")) & (metrics["cutoff"].isin([2018, 2020]))]
    models = list(dict.fromkeys(strict["model"]))
    fig, ax = plt.subplots(figsize=(7.2, 3.8))
    x = np.arange(len(models))
    for offset, year, color in ((-0.18, 2018, blue), (0.18, 2020, orange)):
        vals = [float(strict[(strict["model"].eq(m)) & (strict["cutoff"].eq(year))]["fde_at_1000"].iloc[0]) if len(strict[(strict["model"].eq(m)) & (strict["cutoff"].eq(year))]) else np.nan for m in models]
        ax.bar(x + offset, vals, width=0.36, color=color, label=str(year))
    ax.axhline(1, color="#444444", lw=0.6)
    ax.set_xticks(x, models, rotation=55, ha="right")
    ax.set_ylabel("FDE@1000")
    ax.legend(frameon=False)
    ax.set_title("Strict historical benchmark", loc="left")
    fig.tight_layout()
    fig.savefig(FIG / "figure2b_fde1000.pdf")
    fig.savefig(FIG / "figure2b_fde1000.svg")
    plt.close()

    fig, ax = plt.subplots(figsize=(6.2, 3.6))
    focus = [m for m in ["M3_product", "M5_ligand_max", "M7_svd", "M8_distmult", "M9_lightgcn"] if m in models]
    sub = strict[strict["cutoff"].eq(2018)]
    vals = [float(sub[sub["model"].eq(m)]["fde_at_1000"].iloc[0]) for m in focus if len(sub[sub["model"].eq(m)])]
    labels = [m for m in focus if len(sub[sub["model"].eq(m)])]
    ax.barh(labels[::-1], vals[::-1], color=[orange if m in POP_MODELS else blue for m in labels][::-1])
    ax.axvline(1, color="#444444", lw=0.6)
    ax.set_xlabel("FDE@1000, 2018 strict")
    fig.tight_layout()
    fig.savefig(FIG / "figure2c_popularity_vs_models.pdf")
    fig.savefig(FIG / "figure2c_popularity_vs_models.svg")
    plt.close()

    fig, ax = plt.subplots(figsize=(6.4, 3.6))
    if "pae_ratio" in metrics.columns:
        pae = metrics[(metrics["arm"].eq("STRICT")) & (metrics["label"].eq("Y36")) & (metrics["cutoff"].eq(2018))]
        ax.barh(pae["model"], pae["pae_ratio"], color=green)
        ax.axvline(1, color="#444444", lw=0.6)
        ax.set_xlabel("Observed/expected positives in Top 1000")
    ax.set_title("Popularity-matched enrichment, 2018", loc="left")
    fig.tight_layout()
    fig.savefig(FIG / "figure2d_pae.pdf")
    fig.savefig(FIG / "figure2d_pae.svg")
    plt.close()

    fig, ax = plt.subplots(figsize=(6.2, 3.6))
    if len(stability):
        for year, color in ((2018, blue), (2020, orange)):
            sub = stability[stability["cutoff"].eq(year)]
            ax.plot(sub["k"], sub["jaccard"], marker="o", color=color, label=str(year))
        ax.set_xlabel("Top K")
        ax.set_ylabel("Jaccard, strict graph vs retro graph")
        ax.legend(frameon=False)
    ax.set_title("Common-clean rank stability", loc="left")
    fig.tight_layout()
    fig.savefig(FIG / "figure2e_topk_stability.pdf")
    fig.savefig(FIG / "figure2e_topk_stability.svg")
    plt.close()

    fig, ax = plt.subplots(figsize=(6.2, 3.6))
    if len(time_rows):
        for name, color in zip(time_rows["stratum"].unique(), [vermillion, orange, sky, blue, green]):
            sub = time_rows[time_rows["stratum"].eq(name)].sort_values("day")
            ax.plot(sub["day"], sub["cumulative"], color=color, lw=1.0, label=str(name))
        ax.set_xlabel("Days after cutoff")
        ax.set_ylabel("Cumulative observed discovery")
        ax.legend(frameon=False, fontsize=6)
    ax.set_title("Time to ChEMBL appearance, 2018", loc="left")
    ax.text(0.0, -0.22, "exploratory; time to database release, not biological onset. " + censor_note, transform=ax.transAxes, fontsize=7, color="#444444")
    fig.tight_layout()
    fig.savefig(FIG / "figure2f_time_to_discovery.pdf")
    fig.savefig(FIG / "figure2f_time_to_discovery.svg")
    plt.close()

    fig, ax = plt.subplots(figsize=(6.2, 3.6))
    harm = metrics[(metrics["arm"].eq("STRICT")) & (metrics["label"].eq("Y_HARMONIZED"))]
    if len(harm):
        models_h = list(dict.fromkeys(harm["model"]))[:8]
        x = np.arange(len(models_h))
        for offset, year, color in ((-0.25, 2018, blue), (0, 2020, orange), (0.25, 2022, vermillion)):
            vals = []
            for m in models_h:
                hit = harm[(harm["model"].eq(m)) & (harm["cutoff"].eq(year))]
                vals.append(float(hit["fde_at_1000"].iloc[0]) if len(hit) else np.nan)
            ax.bar(x + offset, vals, width=0.25, color=color, label=str(year))
        ax.axhline(1, color="#444444", lw=0.6)
        ax.set_xticks(x, models_h, rotation=55, ha="right")
        ax.set_ylabel("FDE@1000, harmonized window")
        ax.legend(frameon=False)
    ax.set_title("2022 is right-censored replication", loc="left")
    fig.tight_layout()
    fig.savefig(FIG / "figureS_2022_harmonized.pdf")
    fig.savefig(FIG / "figureS_2022_harmonized.svg")
    plt.close()


def self_test() -> None:
    rng = np.random.default_rng(0)
    n_c, n_t = 40, 15
    c = np.repeat(np.arange(n_c), 3).astype(np.int32)
    t = np.tile((np.arange(n_c) % n_t), 3).astype(np.int32)[: len(c)]
    u, s, vt = fit_svd(c, t, n_c, n_t, 8, 0)
    auc = reconstruction_auc(c, t, lambda cc, tt: np.sum(u[cc] * s * vt[:, tt].T, axis=1), n_t, rng)
    eu, ei, rel = train_distmult(c, t, n_c, n_t, 8, 8, 0.5, 1e-5, 5, 0)
    auc2 = reconstruction_auc(c, t, lambda cc, tt: np.sum(eu[cc] * rel * ei[tt], axis=1), n_t, rng)
    eu, ei = train_lightgcn(c, t, n_c, n_t, 8, 8, 0.1, 1e-4, 3, 0, layers=2)
    auc3 = reconstruction_auc(c, t, lambda cc, tt: np.sum(eu[cc] * ei[tt], axis=1), n_t, rng)
    print(f"self-test auc svd={auc:.3f} distmult={auc2:.3f} lightgcn={auc3:.3f}")
    if min(auc, auc2, auc3) < 0.7:
        raise SystemExit("model self-test failed")



# Library-only copy: obsolete frozen-label CLI orchestration intentionally omitted.

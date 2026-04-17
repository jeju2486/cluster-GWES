#!/usr/bin/env python3
"""
tree_weighted_conditional_mi.py

Compute tree-weighted pooled MI and tree-weighted conditional MI for candidate
locus pairs using:
  - a Newick tree with branch lengths,
  - a space-delimited pair file,
  - a fake FASTA with A=0 and C=1.

Main outputs:
  - per-pair TSV with:
      u, v, distance, [optional extra pair columns],
      mi_tree_pool, mi_tree_cond, delta_mi,
      mono_weight_mass, informative_weight_mass
  - sidecar TSVs:
      <out>.tip_weights.tsv
      <out>.cluster_summary.tsv

Assumptions:
  - pair file has no header
  - pair columns 1 and 2 are 0-based locus indices
  - pair column 3 is carried through as "distance"
  - fake FASTA contains only A/C characters
"""

from __future__ import annotations

import argparse
import csv
import itertools
import math
import re
import sys
from pathlib import Path
from typing import Iterator, List, Sequence, Tuple

import numpy as np

try:
    from Bio import Phylo
except ImportError as e:
    raise SystemExit(
        "Missing dependency: biopython\n"
        "Install with: pip install biopython"
    ) from e


# -----------------------------
# Utilities
# -----------------------------

def eprint(*args, **kwargs) -> None:
    print(*args, file=sys.stderr, **kwargs)


def ensure_parent_dir(path: str | Path) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)


def is_blank_or_comment(line: str) -> bool:
    s = line.strip()
    return (not s) or s.startswith("#")


# -----------------------------
# Input readers
# -----------------------------

def read_tree(tree_path: str, midpoint_root: bool = True):
    tree = Phylo.read(tree_path, "newick")
    if (not tree.rooted) and midpoint_root:
        tree.root_at_midpoint()
    return tree


def read_fake_fasta(fasta_path: str) -> Tuple[List[str], np.ndarray]:
    """
    Reads FASTA-like alignment with:
      A = 0
      C = 1

    Returns:
      tip_names: list[str]
      G: uint8 array of shape (N, L)
    """
    tip_names: List[str] = []
    seqs: List[str] = []

    current_name = None
    current_chunks: List[str] = []

    with open(fasta_path, "r", encoding="utf-8") as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_name is not None:
                    tip_names.append(current_name)
                    seqs.append("".join(current_chunks).upper())
                current_name = line[1:].strip()
                current_chunks = []
            else:
                if current_name is None:
                    raise ValueError("FASTA parse error: sequence encountered before first header.")
                current_chunks.append(line)

    if current_name is not None:
        tip_names.append(current_name)
        seqs.append("".join(current_chunks).upper())

    if not tip_names:
        raise ValueError("No sequences found in fake FASTA.")

    if len(set(tip_names)) != len(tip_names):
        raise ValueError("Duplicate FASTA headers found.")

    lengths = {len(s) for s in seqs}
    if len(lengths) != 1:
        raise ValueError(f"FASTA sequences do not all have the same length: {sorted(lengths)}")
    L = lengths.pop()
    N = len(seqs)

    G = np.empty((N, L), dtype=np.uint8)
    allowed = np.array([ord("A"), ord("C")], dtype=np.uint8)

    for i, seq in enumerate(seqs):
        b = np.frombuffer(seq.encode("ascii"), dtype=np.uint8)
        bad_mask = ~np.isin(b, allowed)
        if bad_mask.any():
            bad_pos = int(np.flatnonzero(bad_mask)[0])
            bad_char = chr(int(b[bad_pos]))
            raise ValueError(
                f"Invalid character in fake FASTA at sequence '{tip_names[i]}', "
                f"position {bad_pos}, char '{bad_char}'. Only A/C are allowed."
            )
        G[i, :] = (b == ord("C")).astype(np.uint8)

    return tip_names, G


def iter_pair_rows(pair_path: str) -> Iterator[Tuple[int, List[str]]]:
    """
    Yields (line_number, fields) for non-empty, non-comment lines.
    No header is assumed.
    """
    with open(pair_path, "r", encoding="utf-8") as fh:
        for line_no, raw in enumerate(fh, start=1):
            if is_blank_or_comment(raw):
                continue
            fields = re.split(r"\s+", raw.strip())
            if len(fields) < 3:
                raise ValueError(
                    f"Pair file line {line_no} has fewer than 3 columns: {raw.rstrip()}"
                )
            yield line_no, fields


def peek_first_pair_row(pair_path: str) -> Tuple[int, List[str]]:
    for line_no, fields in iter_pair_rows(pair_path):
        return line_no, fields
    raise ValueError("Pair file contains no data rows.")


# -----------------------------
# Alignment and tree annotation
# -----------------------------

def align_tree_and_matrix(tree, fasta_tip_names: Sequence[str], G: np.ndarray) -> Tuple[List[str], np.ndarray]:
    """
    Reorder FASTA matrix rows to match tree terminal order.
    Also annotates each terminal with _tip_index.
    """
    terminals = tree.get_terminals()
    tree_tip_names = [tip.name for tip in terminals]

    if any(name is None for name in tree_tip_names):
        raise ValueError("Tree contains tip(s) with missing names.")

    if len(set(tree_tip_names)) != len(tree_tip_names):
        raise ValueError("Duplicate tip names found in tree.")

    fasta_set = set(fasta_tip_names)
    tree_set = set(tree_tip_names)

    only_in_tree = sorted(tree_set - fasta_set)
    only_in_fasta = sorted(fasta_set - tree_set)

    if only_in_tree or only_in_fasta:
        msg = []
        if only_in_tree:
            msg.append(f"Tips present in tree but missing in FASTA: {only_in_tree[:10]}")
        if only_in_fasta:
            msg.append(f"Headers present in FASTA but missing in tree: {only_in_fasta[:10]}")
        raise ValueError("Tree/FASTA tip mismatch.\n" + "\n".join(msg))

    name_to_row = {name: i for i, name in enumerate(fasta_tip_names)}
    row_order = [name_to_row[name] for name in tree_tip_names]
    G_aligned = G[row_order, :]

    for i, tip in enumerate(terminals):
        tip._tip_index = i  # type: ignore[attr-defined]

    return tree_tip_names, G_aligned


def root_and_annotate_tree(tree, n_min: int, L_min: float) -> float:
    """
    Annotates each node with:
      _n_desc
      _subtree_len
      _subtree_len_ratio
      _tip_indices
      _eligible
      _has_eligible
    Returns:
      total tree branch length
    """
    root = tree.root

    def postorder(node) -> None:
        if node.is_terminal():
            node._n_desc = 1
            node._subtree_len = 0.0
            node._tip_indices = [node._tip_index]
            return

        n_desc = 0
        subtree_len = 0.0
        tip_indices: List[int] = []

        for child in node.clades:
            postorder(child)
            bl = float(child.branch_length or 0.0)
            n_desc += child._n_desc
            subtree_len += bl + child._subtree_len
            tip_indices.extend(child._tip_indices)

        node._n_desc = n_desc
        node._subtree_len = subtree_len
        node._tip_indices = tip_indices

    postorder(root)
    total_len = float(root._subtree_len)

    if total_len < 0:
        raise ValueError("Tree total branch length is negative, which is invalid.")

    def preorder_mark(node, parent_n: int | None, is_root: bool) -> None:
        node._subtree_len_ratio = (node._subtree_len / total_len) if total_len > 0 else 0.0
        node._eligible = (
            (not is_root)
            and (not node.is_terminal())
            and (node._n_desc >= n_min)
            and (parent_n is not None and (parent_n - node._n_desc) >= n_min)
            and (node._subtree_len_ratio >= L_min)
        )
        for child in node.clades:
            preorder_mark(child, node._n_desc, False)

    preorder_mark(root, None, True)

    def postorder_has_eligible(node) -> None:
        child_has = False
        for child in node.clades:
            postorder_has_eligible(child)
            child_has = child_has or bool(child._has_eligible)
        node._has_eligible = bool(node._eligible or child_has)

    postorder_has_eligible(root)
    return total_len


# -----------------------------
# Cluster construction
# -----------------------------

def build_clusters(tree) -> Tuple[np.ndarray, List[np.ndarray], int]:
    """
    Build final non-overlapping antichain partition.

    Rule used:
      - if a child is eligible, select it as a final cluster and do not descend
      - else if that child has eligible descendants, recurse into it
      - else take the whole child subtree as one residual final cluster
    """
    root = tree.root
    clusters: List[List[int]] = []

    def recurse(node) -> None:
        for child in node.clades:
            if bool(child._eligible):
                clusters.append(sorted(child._tip_indices))
            elif bool(child._has_eligible):
                recurse(child)
            else:
                clusters.append(sorted(child._tip_indices))

    recurse(root)

    N = len(tree.get_terminals())
    if not clusters:
        # Degenerate fallback: one cluster containing all tips
        clusters = [sorted(root._tip_indices)]

    cluster_tip_indices: List[np.ndarray] = []
    cluster_id_per_tip = np.full(N, -1, dtype=np.int32)

    for cid, tip_idx_list in enumerate(clusters):
        arr = np.array(tip_idx_list, dtype=np.int32)
        if arr.size == 0:
            continue
        cluster_tip_indices.append(arr)
        cluster_id_per_tip[arr] = cid

    if np.any(cluster_id_per_tip < 0):
        missing = np.flatnonzero(cluster_id_per_tip < 0).tolist()
        raise RuntimeError(f"Some tips were not assigned to any cluster: {missing[:10]}")

    return cluster_id_per_tip, cluster_tip_indices, len(cluster_tip_indices)


# -----------------------------
# Fair-proportion tip weights
# -----------------------------

def compute_tip_weights(tree) -> np.ndarray:
    """
    Fair-proportion novelty weight:
      q_t = sum_{e in path(root, t)} branch_length(e) / descendant_tips(e)

    Normalized to sum to N:
      a_t = N q_t / sum(q)
    """
    def preorder(node, acc: float) -> None:
        node._fp_acc = acc
        for child in node.clades:
            bl = float(child.branch_length or 0.0)
            desc = int(child._n_desc)
            add = (bl / desc) if desc > 0 else 0.0
            preorder(child, acc + add)

    preorder(tree.root, 0.0)

    terminals = tree.get_terminals()
    q = np.array([float(tip._fp_acc) for tip in terminals], dtype=float)
    N = q.size
    q_sum = float(q.sum())

    if q_sum <= 0.0:
        # Fallback for zero-length trees
        return np.ones(N, dtype=float)

    a = (N * q) / q_sum
    return a


def compute_cluster_weights(cluster_tip_indices: Sequence[np.ndarray], a_t: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    A_c = np.array([float(a_t[idx].sum()) for idx in cluster_tip_indices], dtype=float)
    N = float(a_t.size)
    W_c = A_c / N
    return A_c, W_c


# -----------------------------
# MI calculations
# -----------------------------

def weighted_counts_2x2(x: np.ndarray, y: np.ndarray, w: np.ndarray) -> np.ndarray:
    """
    Returns counts in order:
      [m00, m01, m10, m11]
    """
    x1 = (x == 1)
    y1 = (y == 1)
    x0 = ~x1
    y0 = ~y1

    m00 = float(w[x0 & y0].sum())
    m01 = float(w[x0 & y1].sum())
    m10 = float(w[x1 & y0].sum())
    m11 = float(w[x1 & y1].sum())

    return np.array([m00, m01, m10, m11], dtype=float)


def smoothed_mi_from_counts(counts: np.ndarray, alpha: float = 0.5) -> float:
    """
    counts order: [m00, m01, m10, m11]
    smoothing:
      p_ij = (m_ij + alpha) / (sum(m) + 4*alpha)
    MI in bits.
    """
    counts = np.asarray(counts, dtype=float)
    total = float(counts.sum())
    denom = total + 4.0 * alpha

    if denom <= 0:
        return 0.0

    p = (counts + alpha) / denom
    p = p.reshape(2, 2)  # rows X=0/1, cols Y=0/1

    p_i = p.sum(axis=1, keepdims=True)
    p_j = p.sum(axis=0, keepdims=True)

    mi = float(np.sum(p * np.log2(p / (p_i * p_j))))
    return mi


def compute_pair_metrics(
    u: int,
    v: int,
    G: np.ndarray,
    cluster_tip_indices: Sequence[np.ndarray],
    a_t: np.ndarray,
    A_c: np.ndarray,
    W_c: np.ndarray,
    alpha: float = 0.5,
) -> Tuple[float, float, float, float, float]:
    """
    Returns:
      mi_tree_pool, mi_tree_cond, delta_mi, mono_weight_mass, informative_weight_mass
    """
    X = G[:, u]
    Y = G[:, v]

    counts_pool = weighted_counts_2x2(X, Y, a_t)
    mi_tree_pool = smoothed_mi_from_counts(counts_pool, alpha=alpha)

    mi_tree_cond = 0.0
    mono_weight_mass = 0.0

    for idx, weight_c in zip(cluster_tip_indices, W_c):
        x_c = X[idx]
        y_c = Y[idx]
        w_c = a_t[idx]

        counts_c = weighted_counts_2x2(x_c, y_c, w_c)
        mi_c = smoothed_mi_from_counts(counts_c, alpha=alpha)
        mi_tree_cond += float(weight_c) * mi_c

        # Monostate defined on raw unweighted cluster size/state constancy
        n_c = int(idx.size)
        s_x = int(x_c.sum())
        s_y = int(y_c.sum())
        if s_x == 0 or s_x == n_c or s_y == 0 or s_y == n_c:
            mono_weight_mass += float(weight_c)

    informative_weight_mass = max(0.0, 1.0 - mono_weight_mass)
    delta_mi = mi_tree_cond - mi_tree_pool

    return (
        mi_tree_pool,
        mi_tree_cond,
        delta_mi,
        mono_weight_mass,
        informative_weight_mass,
    )


# -----------------------------
# Output helpers
# -----------------------------

def write_global_metadata(
    out_path: str,
    tip_names: Sequence[str],
    a_t: np.ndarray,
    cluster_id_per_tip: np.ndarray,
    cluster_tip_indices: Sequence[np.ndarray],
    A_c: np.ndarray,
    W_c: np.ndarray,
) -> None:
    tip_out = f"{out_path}.tip_weights.tsv"
    clu_out = f"{out_path}.cluster_summary.tsv"

    ensure_parent_dir(tip_out)
    ensure_parent_dir(clu_out)

    with open(tip_out, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["tip_index", "tip_name", "tip_weight", "cluster_id"])
        for i, (name, wt, cid) in enumerate(zip(tip_names, a_t, cluster_id_per_tip)):
            w.writerow([i, name, f"{wt:.12g}", int(cid)])

    with open(clu_out, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["cluster_id", "n_tips", "effective_mass", "cluster_weight"])
        for cid, (idx, eff, wt) in enumerate(zip(cluster_tip_indices, A_c, W_c)):
            w.writerow([cid, int(idx.size), f"{eff:.12g}", f"{wt:.12g}"])


# -----------------------------
# Main pipeline
# -----------------------------

def run_all(
    tree_path: str,
    pair_path: str,
    fasta_path: str,
    n_min: int,
    L_min: float,
    out_path: str,
    alpha: float = 0.5,
    midpoint_root: bool = True,
    keep_extra_cols: bool = True,
) -> None:
    if n_min < 1:
        raise ValueError("n_min must be >= 1")
    if not (0.0 <= L_min <= 1.0):
        raise ValueError("L_min must be in [0, 1]")
    if alpha <= 0:
        raise ValueError("alpha must be > 0")

    eprint("[1/6] Reading tree")
    tree = read_tree(tree_path, midpoint_root=midpoint_root)

    eprint("[2/6] Reading fake FASTA")
    fasta_tip_names, G = read_fake_fasta(fasta_path)
    N, L = G.shape
    eprint(f"      tips={N}, loci={L}")

    eprint("[3/6] Aligning FASTA rows to tree tip order")
    tree_tip_names, G = align_tree_and_matrix(tree, fasta_tip_names, G)

    eprint("[4/6] Root/annotate tree and build clusters")
    total_len = root_and_annotate_tree(tree, n_min=n_min, L_min=L_min)
    cluster_id_per_tip, cluster_tip_indices, K = build_clusters(tree)
    eprint(f"      total_tree_length={total_len:.12g}")
    eprint(f"      clusters={K}")

    eprint("[5/6] Computing fair-proportion tip weights and cluster weights")
    a_t = compute_tip_weights(tree)
    A_c, W_c = compute_cluster_weights(cluster_tip_indices, a_t)
    eprint(f"      sum_tip_weights={a_t.sum():.12g}")
    eprint(f"      sum_cluster_weights={W_c.sum():.12g}")

    write_global_metadata(
        out_path=out_path,
        tip_names=tree_tip_names,
        a_t=a_t,
        cluster_id_per_tip=cluster_id_per_tip,
        cluster_tip_indices=cluster_tip_indices,
        A_c=A_c,
        W_c=W_c,
    )

    eprint("[6/6] Processing pairs and writing output")
    first_line_no, first_fields = peek_first_pair_row(pair_path)
    n_extra = max(0, len(first_fields) - 3)

    out_cols = ["u", "v", "distance"]
    if keep_extra_cols:
        out_cols.extend([f"pair_col{i}" for i in range(4, 4 + n_extra)])
    out_cols.extend(
        [
            "mi_tree_pool",
            "mi_tree_cond",
            "delta_mi",
            "mono_weight_mass",
            "informative_weight_mass",
        ]
    )

    ensure_parent_dir(out_path)
    n_written = 0

    with open(out_path, "w", encoding="utf-8", newline="") as out_fh:
        writer = csv.writer(out_fh, delimiter="\t")
        writer.writerow(out_cols)

        pair_iter = iter_pair_rows(pair_path)
        all_rows = itertools.chain([(first_line_no, first_fields)], pair_iter)

        seen_first = False
        for line_no, fields in all_rows:
            # Avoid processing the first row twice because iter_pair_rows starts from top again
            if (line_no == first_line_no) and (fields == first_fields) and seen_first:
                continue
            seen_first = True

            try:
                u = int(fields[0])
                v = int(fields[1])
            except ValueError as e:
                raise ValueError(
                    f"Pair file line {line_no}: first two columns must be integers."
                ) from e

            if not (0 <= u < L) or not (0 <= v < L):
                raise IndexError(
                    f"Pair file line {line_no}: locus index out of range. "
                    f"u={u}, v={v}, valid range=[0,{L-1}]"
                )

            distance = fields[2]
            extra = fields[3:] if keep_extra_cols else []

            (
                mi_tree_pool,
                mi_tree_cond,
                delta_mi,
                mono_weight_mass,
                informative_weight_mass,
            ) = compute_pair_metrics(
                u=u,
                v=v,
                G=G,
                cluster_tip_indices=cluster_tip_indices,
                a_t=a_t,
                A_c=A_c,
                W_c=W_c,
                alpha=alpha,
            )

            row = [u, v, distance]
            if keep_extra_cols:
                row.extend(extra)
            row.extend(
                [
                    f"{mi_tree_pool:.12g}",
                    f"{mi_tree_cond:.12g}",
                    f"{delta_mi:.12g}",
                    f"{mono_weight_mass:.12g}",
                    f"{informative_weight_mass:.12g}",
                ]
            )
            writer.writerow(row)
            n_written += 1

            if n_written % 10000 == 0:
                eprint(f"      processed {n_written} pairs")

    eprint(f"Done. Wrote {n_written} pair rows to: {out_path}")
    eprint(f"Sidecars:")
    eprint(f"  {out_path}.tip_weights.tsv")
    eprint(f"  {out_path}.cluster_summary.tsv")


# -----------------------------
# CLI
# -----------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Compute tree-weighted pooled MI and tree-weighted conditional MI."
    )
    p.add_argument("--tree", required=True, help="Newick tree file with branch lengths")
    p.add_argument("--pairs", required=True, help="Whitespace-delimited pair file (no header)")
    p.add_argument("--fasta", required=True, help="Fake FASTA with A=0 and C=1")
    p.add_argument("--n-min", type=int, required=True, help="Minimum cluster size")
    p.add_argument("--L-min", type=float, required=True, help="Minimum relative subtree length in [0,1]")
    p.add_argument("--out", required=True, help="Output TSV path")
    p.add_argument("--alpha", type=float, default=0.5, help="Pseudocount alpha (default: 0.5)")
    p.add_argument(
        "--no-midpoint-root",
        action="store_true",
        help="Do not midpoint-root when tree is unrooted",
    )
    p.add_argument(
        "--drop-extra-cols",
        action="store_true",
        help="Do not carry through pair columns after distance",
    )
    return p.parse_args()


def main() -> None:
    args = parse_args()
    run_all(
        tree_path=args.tree,
        pair_path=args.pairs,
        fasta_path=args.fasta,
        n_min=args.n_min,
        L_min=args.L_min,
        out_path=args.out,
        alpha=args.alpha,
        midpoint_root=(not args.no_midpoint_root),
        keep_extra_cols=(not args.drop_extra_cols),
    )


if __name__ == "__main__":
    main()

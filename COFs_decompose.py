#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Batch CIF -> unique fragments (global dedup) with canonical SMILES (isomericSmiles=True).
SQLite-backed dedup + resume + parallel workers.
Also builds a CSV mapping each global unique fragment to the list of CIFs that contain it.

Usage:
  python COFs_decompose.py \
      --cif-dir ./cifs \
      --out-dir ./uniq_xyz \
      --record ./fragments_record.txt \
      --db ./fragments_seen.sqlite \
      --csv-out ./fragment_to_cifs.csv \
      [--repeat 3 3 3] [--scale 1.10] [--span-warn 0.75] [--workers 16] [--max-files 0] [--resume]

Dependencies: rdkit, ase, numpy
"""

import os
import sys
import csv
import argparse
from pathlib import Path
from dataclasses import dataclass
from typing import Iterable, Tuple, List, Optional, Dict, Any

import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed

# ---- Make threading polite on shared HPCs ----
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

# ASE
from ase.io import read as ase_read
from ase.neighborlist import neighbor_list
from ase.data import covalent_radii

# RDKit (quiet)
from rdkit import rdBase
rdBase.DisableLog('rdApp.*')
from rdkit import Chem
from rdkit.Chem.rdmolops import SanitizeFlags

# SQLite
import sqlite3
from contextlib import closing

# ---------------- Ring & cut logic ----------------

@dataclass
class RingData:
    bond_rings: list
    atom_rings: list
    bond_in_any_ring: set
    rings_of_sizes: dict

def compute_ring_data(mol: Chem.Mol) -> RingData:
    info = mol.GetRingInfo()
    atom_rings = [tuple(r) for r in info.AtomRings()]
    bond_rings = [tuple(r) for r in info.BondRings()]
    bond_in_any = set(b for br in bond_rings for b in br)
    sized: Dict[int, List[Tuple[int, ...]]] = {}
    for r in atom_rings:
        sized.setdefault(len(r), []).append(r)
    return RingData(bond_rings, atom_rings, bond_in_any, sized)

CENTER_TO_NEIGHBOR_PAIRS = {
    "N": {("O", "O"), ("B", "B"), ("C", "C")},
    "B": {("O", "O")},
    "Si": {("O", "O")},
    "C": {("N", "N")},
}

def _two_ring_neighbors_in_cycle(atom_idx: int, cycle: List[int]) -> Tuple[int, int]:
    n = len(cycle)
    k = cycle.index(atom_idx)
    return cycle[(k - 1) % n], cycle[(k + 1) % n]

def find_cut_bonds_by_center_first(mol, ring_sizes=(4, 5, 6), require_degree3=True):
    """
    Returns:
        selected_bonds (sorted list of bond indices), stats (dict)

    Pass 1: Center-first patterns on small rings.
    Pass 2: N–N bonds not in any 5/6-member ring. --> changed to not in small rings
    Pass 3: C–N bonds where:
            - C is in at least one 6-member ring (no requirement for 5-rings),
            - N is NOT in any 5- or 6-member ring,
            - CN bond is not in a small ring,
            - veto if N has exactly two other neighbors that are both H or both O.

    Lone-atom guard: accept cuts greedily by pass priority, skipping any cut that would
    leave either endpoint with degree 0 after the set of cuts applied so far.
    """
    ringdata = compute_ring_data(mol)
    def sym(i): return mol.GetAtomWithIdx(i).GetSymbol()

    stats = {
        "pattern_hits": 0,   # pass 1 candidates found
        "nn_hits": 0,        # pass 2 candidates found (N–N rule)
        "cn_hits": 0,        # pass 3 candidates found (C–N rule)
        "accepted": 0,       # cuts actually kept after lone-atom guard
        "skipped_to_avoid_lone": 0
    }

    # ---- Precompute ring helpers ----
    # Bonds in any ring with size in ring_sizes (for pass 1 veto)
    small_ring_bonds = set()
    for br in ringdata.bond_rings:
        if len(br) in ring_sizes:
            small_ring_bonds.update(br)

    # Atom membership in 5/6-member rings (for N veto, old behavior)
    atoms_in_5_or_6 = set()
    for sz in (5, 6):
        for ring in ringdata.rings_of_sizes.get(sz, []):
            atoms_in_5_or_6.update(ring)

    # Atom membership in 6-member rings only (for C requirement in pass 3)
    atoms_in_6 = set()
    for ring in ringdata.rings_of_sizes.get(6, []):
        atoms_in_6.update(ring)

    # Bonds in 5/6-member rings (for N–N veto)
    bonds_in_5_or_6 = set()
    for br in ringdata.bond_rings:
        if len(br) in (5, 6):
            bonds_in_5_or_6.update(br)

    # ---- Collect candidates with pass priorities (1 > 2 > 3). Earlier pass wins ties. ----
    candidates_with_priority = {}  # bidx -> priority

    # ---------------- Pass 1: center-first patterns on small rings ----------------
    for sz in ring_sizes:
        for ring in ringdata.rings_of_sizes.get(sz, []):
            # skip all-carbon cycles
            if all(sym(i) == "C" for i in ring):
                continue
            cycle = list(ring)
            for c_idx in cycle:
                c_atom = mol.GetAtomWithIdx(c_idx)
                csym = c_atom.GetSymbol()
                if csym not in CENTER_TO_NEIGHBOR_PAIRS:
                    continue
                if require_degree3 and c_atom.GetDegree() != 3:
                    continue
                left_idx, right_idx = _two_ring_neighbors_in_cycle(c_idx, cycle)
                pair = tuple(sorted((sym(left_idx), sym(right_idx))))
                if pair not in CENTER_TO_NEIGHBOR_PAIRS[csym]:
                    continue
                # third neighbor (outside the ring)
                third = next(
                    (n for n in c_atom.GetNeighbors() if n.GetIdx() not in {left_idx, right_idx}),
                    None
                )
                if third is None:
                    continue
                b = mol.GetBondBetweenAtoms(c_idx, third.GetIdx())
                if b is None:
                    continue
                bidx = b.GetIdx()
                # veto bonds that belong to any small ring (4/5/6 by default)
                if bidx in small_ring_bonds:
                    continue
                if bidx not in candidates_with_priority:
                    candidates_with_priority[bidx] = 1
                    stats["pattern_hits"] += 1

    # ---------------- Pass 2: N–N bonds NOT in any 5/6-member ring ----------------
    for b in mol.GetBonds():
        a1, a2 = b.GetBeginAtom(), b.GetEndAtom()
        if a1.GetSymbol() == "N" and a2.GetSymbol() == "N":
            bidx = b.GetIdx()
            # if bidx in bonds_in_5_or_6:
            #     continue  # only cut N–N if not in 5/6 ring
            if bidx in small_ring_bonds:
                continue  # only cut N–N if not in small ring
            if bidx not in candidates_with_priority:
                candidates_with_priority[bidx] = 2
                stats["nn_hits"] += 1

    # ---------------- Pass 3: C–N ring-asymmetry rule (C in 6-ring, N not in 5/6) ----------------
    for b in mol.GetBonds():
        a1, a2 = b.GetBeginAtom(), b.GetEndAtom()
        s1, s2 = a1.GetSymbol(), a2.GetSymbol()
        if {s1, s2} != {"C", "N"}:
            continue
        bidx = b.GetIdx()
        if bidx in small_ring_bonds:
            continue

        c_atom = a1 if s1 == "C" else a2
        n_atom = a2 if s1 == "C" else a1

        # NEW: require C in at least one 6-member ring specifically
        c_in_6 = c_atom.GetIdx() in atoms_in_6
        # keep N veto as "not in any 5 or 6 ring"
        n_in_56 = n_atom.GetIdx() in atoms_in_5_or_6

        if not (c_in_6 and (not n_in_56)):
            continue

        # veto HH or OO as the two *other* neighbors on N
        other_neighbors = [nbr for nbr in n_atom.GetNeighbors() if nbr.GetIdx() != c_atom.GetIdx()]
        if len(other_neighbors) == 2:
            sA, sB = other_neighbors[0].GetSymbol(), other_neighbors[1].GetSymbol()
            if (sA == "H" and sB == "H") or (sA == "O" and sB == "O"):
                continue

        if bidx not in candidates_with_priority:
            candidates_with_priority[bidx] = 3
            stats["cn_hits"] += 1

    # ---- Lone-atom guard: greedy accept by (priority, bond index), updating endpoint degrees ----
    deg_remaining = np.array(
        [mol.GetAtomWithIdx(i).GetDegree() for i in range(mol.GetNumAtoms())],
        dtype=int
    )
    ordered = sorted(candidates_with_priority.items(), key=lambda kv: (kv[1], kv[0]))

    final_cuts = []
    for bidx, _prio in ordered:
        b = mol.GetBondWithIdx(bidx)
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        if (deg_remaining[i] - 1) >= 1 and (deg_remaining[j] - 1) >= 1:
            final_cuts.append(bidx)
            deg_remaining[i] -= 1
            deg_remaining[j] -= 1
        else:
            stats["skipped_to_avoid_lone"] += 1

    stats["accepted"] = len(final_cuts)
    return sorted(final_cuts), stats

# ---------------- PBC graph & central-image mapping ----------------

def build_supercell(prim_atoms, repeat=(3, 3, 3)):
    return prim_atoms.repeat(repeat)

def get_central_indices(supercell, lower=1/3, upper=2/3, tol=1e-12):
    frac = supercell.get_scaled_positions()
    keep = []
    for i, (u, v, w) in enumerate(frac):
        if (lower - tol <= u < upper - tol and
            lower - tol <= v < upper - tol and
            lower - tol <= w < upper - tol):
                keep.append(i)
    return set(keep)

#was using scale=1.10 and max_r = float(radii.max()) * 2.2
def build_periodic_adjacency(supercell, scale=1.20):
    Z = supercell.get_atomic_numbers()
    cell = supercell.get_cell()
    pos  = supercell.get_positions()
    radii = np.array([covalent_radii[z] if z < len(covalent_radii) else 0.77 for z in Z])
    max_r = float(radii.max()) * 2.4
    i_idx, j_idx, Shifts = neighbor_list('ijS', supercell, max_r)
    adj: Dict[int, List[Tuple[int, Tuple[int,int,int]]]] = {}
    for i, j, S in zip(i_idx, j_idx, Shifts):
        if i >= j:
            continue
        ri, rj = radii[i], radii[j]
        cutoff = (ri + rj) * scale
        disp = (pos[j] + S @ cell) - pos[i]
        dist = np.linalg.norm(disp)
        if dist <= cutoff:
            tS = (int(S[0]), int(S[1]), int(S[2]))
            adj.setdefault(i, []).append((j, tS))
            adj.setdefault(j, []).append((i, (-tS[0], -tS[1], -tS[2])))
    return adj

def supercell_to_rdkit_central(supercell, prim_count, central_indices, adj):
    rw = Chem.RWMol()
    Z = supercell.get_atomic_numbers()
    idx_map = {}
    for sc_idx in sorted(central_indices):
        a = Chem.Atom(int(Z[sc_idx]))
        idx_map[sc_idx] = rw.AddAtom(a)

    reps = {c % prim_count: c for c in central_indices}
    for i in sorted(central_indices):
        for j, S in adj.get(i, []):
            j_cent = reps.get(j % prim_count, None)
            if j_cent is None or j_cent not in idx_map:
                continue
            ii, jj = idx_map[i], idx_map[j_cent]
            if ii >= jj:
                continue
            rw.AddBond(ii, jj, Chem.BondType.SINGLE)

    mol = rw.GetMol()
    Chem.SanitizeMol(mol, sanitizeOps=SanitizeFlags.SANITIZE_ALL)
    return mol

# ---------------- Unwrap fragment under PBC ----------------

def unwrap_fragment_by_frac(cut_mol: Chem.Mol,
                            frag_idx_tuple: Tuple[int, ...],
                            molidx_to_scidx: Dict[int,int],
                            super_pos: np.ndarray,
                            inv_cell: np.ndarray):
    EPS = 1e-12
    real_atoms = [i for i in frag_idx_tuple if cut_mol.GetAtomWithIdx(i).GetAtomicNum() > 0]
    if not real_atoms:
        return [], np.zeros((0, 3)), np.zeros((0, 3))

    raw_frac = {i: super_pos[molidx_to_scidx[i]] @ inv_cell for i in real_atoms}

    adj = {i: [] for i in real_atoms}
    real_set = set(real_atoms)
    for i in real_atoms:
        for b in cut_mol.GetAtomWithIdx(i).GetBonds():
            j = b.GetOtherAtomIdx(i)
            if j in real_set:
                adj[i].append(j)

    assigned: Dict[int, Optional[np.ndarray]] = {i: None for i in real_atoms}
    for root in real_atoms:
        if assigned[root] is not None:
            continue
        assigned[root] = raw_frac[root].copy()
        stack = [root]
        while stack:
            cur = stack.pop()
            f_cur = assigned[cur]
            for nbr in adj[cur]:
                if assigned[nbr] is None:
                    delta_raw = raw_frac[nbr] - raw_frac[cur]
                    assigned[nbr] = f_cur + (delta_raw - np.round(delta_raw))
                    stack.append(nbr)

    for k in assigned:
        a = assigned[k]
        a[np.isclose(a, 0.0, atol=EPS)] = 0.0

    allF = np.vstack([assigned[i] for i in real_atoms])
    shift = np.floor(allF.min(axis=0))
    for k in assigned:
        assigned[k] = assigned[k] - shift

    real_atoms_sorted = sorted(real_atoms)
    symbols = [cut_mol.GetAtomWithIdx(i).GetSymbol() for i in real_atoms_sorted]
    fracs   = np.vstack([assigned[i] for i in real_atoms_sorted])
    return symbols, fracs, fracs

# ---------------- Fingerprints ----------------

def frag_mol_without_dummies(cut_mol: Chem.Mol, frag_idxs: Tuple[int, ...]) -> Optional[Chem.Mol]:
    real = [i for i in frag_idxs if cut_mol.GetAtomWithIdx(i).GetAtomicNum() > 0]
    if not real:
        return None
    amap = {old: new for new, old in enumerate(sorted(real))}
    rw = Chem.RWMol()
    for old in sorted(real):
        rw.AddAtom(Chem.Atom(int(cut_mol.GetAtomWithIdx(old).GetAtomicNum())))
    present = set(real)
    for old in real:
        ai = cut_mol.GetAtomWithIdx(old)
        for b in ai.GetBonds():
            j = b.GetOtherAtomIdx(old)
            if old < j and j in present:
                rw.AddBond(amap[old], amap[j], b.GetBondType())
    sub = rw.GetMol()
    try:
        Chem.SanitizeMol(sub)
    except Exception:
        pass
    return sub

def fingerprint_smiles_isomeric(cut_mol: Chem.Mol, frag_idxs: Tuple[int, ...]) -> str:
    sub = frag_mol_without_dummies(cut_mol, frag_idxs)
    if sub is None or sub.GetNumAtoms() == 0:
        return "NA"
    try:
        smi = Chem.MolToSmiles(sub, canonical=True, isomericSmiles=True)
        return "SMI:" + smi
    except Exception:
        labels = []
        for a in sub.GetAtoms():
            s  = a.GetSymbol()
            nbs = sorted(nb.GetSymbol() for nb in a.GetNeighbors())
            labels.append(f"{s}|{','.join(nbs)}")
        labels.sort()
        return "ENVv1:" + "|".join(labels)

# ---------------- SQLite registry ----------------

def init_seen_db(db_path: Path):
    with closing(sqlite3.connect(db_path)) as con:
        con.execute("""
            CREATE TABLE IF NOT EXISTS seen (
                fp   TEXT PRIMARY KEY,
                cif  TEXT,
                xyz  TEXT,
                idx  INTEGER
            )
        """)
        con.execute("""
            CREATE TABLE IF NOT EXISTS processed (
                cif      TEXT PRIMARY KEY,
                status   TEXT,         -- 'done' or 'error'
                n_frags  INTEGER,
                n_uniqs  INTEGER,
                message  TEXT
            )
        """)
        # New: which CIFs contain which fragments (by fp)
        con.execute("""
            CREATE TABLE IF NOT EXISTS fragment_occurrence (
                fp   TEXT NOT NULL,
                cif  TEXT NOT NULL,
                PRIMARY KEY (fp, cif)
            )
        """)
        con.execute("CREATE INDEX IF NOT EXISTS idx_fragment_occurrence_fp ON fragment_occurrence(fp)")
        con.execute("PRAGMA journal_mode=WAL;")
        con.execute("PRAGMA synchronous=NORMAL;")
        con.commit()

def try_insert_fp(db_path: Path, fp: str, cif_name: str, xyz_path: str, idx: int) -> bool:
    with closing(sqlite3.connect(db_path)) as con:
        try:
            con.execute("INSERT INTO seen(fp, cif, xyz, idx) VALUES (?,?,?,?)",
                        (fp, cif_name, xyz_path, idx))
            con.commit()
            return True
        except sqlite3.IntegrityError:
            return False

def get_idx_for_fp(db_path: Path, fp: str) -> Optional[int]:
    with closing(sqlite3.connect(db_path)) as con:
        cur = con.execute("SELECT idx FROM seen WHERE fp = ?", (fp,))
        row = cur.fetchone()
        return row[0] if row else None

def upsert_fragment_occurrence(db_path: Path, fp: str, cif_name: str):
    with closing(sqlite3.connect(db_path)) as con:
        # On conflict ignore — we only care that (fp,cif) exists
        con.execute("""
            INSERT OR IGNORE INTO fragment_occurrence(fp, cif) VALUES (?,?)
        """, (fp, cif_name))
        con.commit()

def mark_processed(db_path: Path, cif_name: str, status: str, n_frags: int, n_uniqs: int, message: str):
    with closing(sqlite3.connect(db_path)) as con:
        con.execute("""
            INSERT INTO processed(cif, status, n_frags, n_uniqs, message)
            VALUES (?,?,?,?,?)
            ON CONFLICT(cif) DO UPDATE SET
                status=excluded.status,
                n_frags=excluded.n_frags,
                n_uniqs=excluded.n_uniqs,
                message=excluded.message
        """, (cif_name, status, n_frags, n_uniqs, message[:500]))
        con.commit()

def already_processed(db_path: Path, cif_name: str) -> bool:
    with closing(sqlite3.connect(db_path)) as con:
        cur = con.execute("SELECT status FROM processed WHERE cif = ?", (cif_name,))
        row = cur.fetchone()
        return bool(row and row[0] == "done")

# ---------------- Worker ----------------

def process_cif_worker(cif_path: str,
                       repeat: Tuple[int,int,int],
                       scale: float,
                       span_warn: float) -> Dict[str, Any]:
    """
    Return dict with:
      'cif': name,
      'error': str or None,
      'frags': list of {'k': int, 'fp': str, 'symbols': [..], 'cart': np.ndarray},
      'n_frags': int,
      'warns': int
    """
    name = Path(cif_path).name
    try:
        prim = ase_read(cif_path)
    except Exception as e:
        return {"cif": name, "error": f"ASE read error: {e}", "frags": [], "n_frags": 0, "warns": 0}

    try:
        supercell = build_supercell(prim, repeat)
        central_indices = get_central_indices(supercell)
        adj = build_periodic_adjacency(supercell, scale=scale)
        mol = supercell_to_rdkit_central(supercell, len(prim), central_indices, adj)

        selected_bonds, _stats = find_cut_bonds_by_center_first(mol)
        if selected_bonds:
            cut = Chem.FragmentOnBonds(mol, selected_bonds, addDummies=True)
            frag_idx_tuples = Chem.GetMolFrags(cut, asMols=False, sanitizeFrags=False)
        else:
            cut = mol
            frag_idx_tuples = [tuple(range(mol.GetNumAtoms()))]

        cell = np.array(prim.get_cell())
        inv_cell = np.linalg.inv(cell)
        super_pos = supercell.get_positions()
        central_list = sorted(central_indices)
        if mol.GetNumAtoms() != len(central_list):
            return {"cif": name, "error": "central_indices mismatch", "frags": [], "n_frags": 0, "warns": 0}
        molidx_to_scidx = {i: central_list[i] for i in range(mol.GetNumAtoms())}

        results = []
        warned = 0
        for k, frag_idxs in enumerate(frag_idx_tuples):
            fp = fingerprint_smiles_isomeric(cut, frag_idxs)
            symbols, fracs_unwrap, _ = unwrap_fragment_by_frac(
                cut, frag_idxs, molidx_to_scidx, super_pos, inv_cell
            )
            if len(symbols) == 0:
                continue
            span = (fracs_unwrap.max(axis=0) - fracs_unwrap.min(axis=0)) % 1.0
            if np.any(span > span_warn + 1e-12):
                warned += 1
            cart = fracs_unwrap @ cell
            results.append({"k": k, "fp": fp, "symbols": symbols, "cart": cart})

        return {"cif": name, "error": None, "frags": results, "n_frags": len(frag_idx_tuples), "warns": warned}
    except Exception as e:
        return {"cif": name, "error": f"processing error: {e}", "frags": [], "n_frags": 0, "warns": 0}

# ---------------- IO helpers ----------------

def write_xyz(path: Path, symbols: List[str], cart_coords: np.ndarray, title: str = ""):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        f.write(f"{len(symbols)}\n")
        f.write(title.strip() + "\n")
        for s, (x, y, z) in zip(symbols, cart_coords):
            f.write(f"{s:2s} {x:16.8f} {y:16.8f} {z:16.8f}\n")

# ---------------- Reporting ----------------

def write_fragment_to_cifs_csv(db_path: Path, csv_path: Path):
    """
    Output CSV columns:
      fragment_idx, n_cifs, cif_basenames (semicolon-separated, de-duplicated)
    """
    import csv
    with closing(sqlite3.connect(db_path)) as con, open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["fragment_idx", "n_cifs", "cif_basenames"])
        # Deduplicate (idx, cif) pairs in a subquery, THEN aggregate.
        for row in con.execute("""
            SELECT t.fragment_idx,
                   COUNT(*) AS n_cifs,
                   GROUP_CONCAT(t.cif, ';') AS cifs
            FROM (
                SELECT s.idx AS fragment_idx, fo.cif AS cif
                FROM seen s
                JOIN fragment_occurrence fo ON fo.fp = s.fp
                GROUP BY s.idx, fo.cif   -- DISTINCT (idx, cif)
            ) AS t
            GROUP BY t.fragment_idx
            ORDER BY n_cifs DESC, t.fragment_idx ASC
        """):
            fragment_idx, n_cifs, cifs = row
            writer.writerow([fragment_idx, n_cifs, cifs or ""])

# ---------------- Main ----------------

def main():
    ap = argparse.ArgumentParser(description="Loop CIFs, global dedup by isomeric SMILES (SQLite-backed), write unique XYZs + fragment->CIFs CSV.")
    ap.add_argument("--cif-dir", type=Path, required=True, help="Folder containing .cif files (non-recursive).")
    ap.add_argument("--out-dir", type=Path, default=Path("./uniq_xyz"), help="Output folder for xyz.")
    ap.add_argument("--record", type=Path, default=Path("fragments_record.txt"), help="TSV record file to write/append.")
    ap.add_argument("--db", type=Path, default=Path("fragments_seen.sqlite"), help="SQLite DB for dedup + resume.")
    ap.add_argument("--csv-out", type=Path, default=Path("fragment_to_cifs.csv"), help="CSV mapping fragment_idx -> CIF list.")
    ap.add_argument("--repeat", type=int, nargs=3, default=[3, 3, 3], help="Supercell repeat (e.g., 3 3 3).")
    ap.add_argument("--scale", type=float, default=1.20, help="Covalent-radius scaling for bonding.")
    ap.add_argument("--span-warn", type=float, default=0.75, help="Warn if fragment spans > this in any frac dim.")
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) // 2), help="Parallel workers.")
    ap.add_argument("--max-files", type=int, default=0, help="Limit number of CIFs (0=all).")
    ap.add_argument("--resume", action="store_true", help="Skip CIFs already marked as done in DB.")
    args = ap.parse_args()

    if not args.cif_dir.exists():
        print(f"[error] cif-dir not found: {args.cif_dir}", file=sys.stderr)
        sys.exit(2)

    init_seen_db(args.db)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    cif_paths = sorted(p for p in args.cif_dir.iterdir() if p.is_file() and p.suffix.lower() == ".cif")
    if args.max_files > 0:
        cif_paths = cif_paths[:args.max_files]
    if not cif_paths:
        print("[info] no CIF files found.")
        with open(args.record, "a") as rec:
            rec.write("# No CIFs found.\n")
        # Still write an empty CSV header for consistency
        write_fragment_to_cifs_csv(args.db, args.csv_out)
        return

    # Optionally skip already-done files
    if args.resume:
        cif_paths = [p for p in cif_paths if not already_processed(args.db, p.name)]
        if not cif_paths:
            print("[info] nothing to do (resume mode; all done).")
            # Even in resume no-op, refresh CSV from DB
            write_fragment_to_cifs_csv(args.db, args.csv_out)
            return

    print(f"[info] files: {len(cif_paths)} | workers: {args.workers} | out: {args.out_dir}")

    # Open record in append mode once
    rec = open(args.record, "a")
    if rec.tell() == 0:
        rec.write("# CIF\tuniq_index\tfingerprint\txyz_path\n")
        rec.flush()

    # Global unique index counter is derived from current max idx in DB
    with closing(sqlite3.connect(args.db)) as con:
        cur = con.execute("SELECT COALESCE(MAX(idx), -1) FROM seen")
        row = cur.fetchone()
        next_idx = (row[0] or -1) + 1

    submitted = 0
    completed = 0
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        futures = []
        for p in cif_paths:
            futures.append(ex.submit(
                process_cif_worker,
                str(p),
                tuple(args.repeat),
                float(args.scale),
                float(args.span_warn)
            ))
            submitted += 1

        for fut in as_completed(futures):
            res = fut.result()
            cif_name = res["cif"]
            if res["error"]:
                print(f"[fail] {cif_name}: {res['error']}", file=sys.stderr)
                mark_processed(args.db, cif_name, "error", 0, 0, res["error"])
                completed += 1
                continue

            kept_for_this_cif = 0

            # For each fragment from this CIF:
            for frag in res["frags"]:
                fp = frag["fp"]
                symbols = frag["symbols"]
                cart = frag["cart"]

                # Record the occurrence of this fragment fingerprint in this CIF (always)
                upsert_fragment_occurrence(args.db, fp, cif_name)

                out_xyz = args.out_dir / f"{Path(cif_name).stem}_uniq{next_idx:08d}.xyz"

                # Attempt to insert FP (global dedup). If succeeds -> first time seen; write XYZ.
                if try_insert_fp(args.db, fp, cif_name, str(out_xyz), next_idx):
                    title = f"{Path(cif_name).stem} uniq{next_idx:08d} | atoms={len(symbols)} | fp=isomericSMILES"
                    write_xyz(out_xyz, symbols, cart, title=title)
                    rec.write(f"{cif_name}\t{next_idx:08d}\t{fp}\t{out_xyz}\n")
                    rec.flush()
                    next_idx += 1
                    kept_for_this_cif += 1
                # else duplicate globally -> we already logged occurrence; skip writing XYZ

            mark_processed(args.db, cif_name, "done", res["n_frags"], kept_for_this_cif,
                           f"warns={res['warns']}")
            if kept_for_this_cif == 0:
                # also log to the record for quick visibility
                rec.write(f"# no-unique\t{cif_name}\n")
                rec.flush()

            completed += 1
            if completed % 200 == 0:
                print(f"[progress] {completed}/{submitted} CIFs done")

    rec.close()

    # Build the fragment -> CIFs CSV summary
    write_fragment_to_cifs_csv(args.db, args.csv_out)

    print(f"[done] Finished. Unique XYZs in: {args.out_dir}")
    print(f"[info] Record TSV: {args.record} | DB: {args.db} | Fragment→CIF CSV: {args.csv_out}")

if __name__ == "__main__":
    main()


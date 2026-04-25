"""
src/data_loader.py
==================
Dataset loaders for SG-SNIR experiments.

Supported datasets:
  1. HIV Transmission  — ICPSR #22140 DS0002  (N≈1,005, E≈6,726)
  2. p2p-Gnutella08   — SNAP                  (N=6,301,  E=20,777)
  3. Wiki-Vote        — SNAP                  (N=7,115,  E=103,689)
  4. soc-Epinions1    — SNAP                  (N=75,879, E=508,837)
  5. email-EuAll      — SNAP                  (N=265,214,E=420,045)

All datasets are natively directed — no conversion needed.
SNAP files use '#' as comment character and tab-separated node pairs.

Auto-download:
  SNAP datasets are downloaded automatically from snap.stanford.edu on first
  use and cached in  combo/data/snap/.  Subsequent runs use the local cache.
  If paper-2/data/ is available locally, files are copied from there instead
  of re-downloading (avoids redundant network traffic).
  HIV (ICPSR) must be placed manually — see HIV_TSV_PATH below.

Initial state seeding:
  seed_initial_states  — fraction-based (general use)
  seed_fixed_nodes     — fixed count of infected seeds (plan Section 3.3)
"""

import gzip
import random
import shutil
import urllib.request
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import networkx as nx
import numpy as np


# ---------------------------------------------------------------------------
# Dataset paths & download registry
# ---------------------------------------------------------------------------

# Resolve project root = combo/
_PROJECT_ROOT = Path(__file__).parent.parent

# Local cache directory — all downloaded files land here
SNAP_CACHE_DIR = _PROJECT_ROOT / "data" / "snap"
SNAP_CACHE_DIR.mkdir(parents=True, exist_ok=True)

# HIV TSV must be placed manually (ICPSR requires registration).
# Expected location: combo/data/ICPSR_22140/DS0002/22140-0002-Data.tsv
HIV_TSV_PATH = _PROJECT_ROOT / "data" / "ICPSR_22140" / "DS0002" / "22140-0002-Data.tsv"

# SNAP download registry: name → (remote_url, local_filename, gzipped?)
# All URLs use the public SNAP download endpoint.
_SNAP_REGISTRY: Dict[str, tuple] = {
    "p2p-Gnutella": (
        "https://snap.stanford.edu/data/p2p-Gnutella08.txt.gz",
        "p2p-Gnutella08.txt",
        True,
    ),
    "Wiki-Vote": (
        "https://snap.stanford.edu/data/wiki-Vote.txt.gz",
        "Wiki-Vote.txt",
        True,
    ),
    "soc-Epinions1": (
        "https://snap.stanford.edu/data/soc-Epinions1.txt.gz",
        "soc-Epinions1.txt",
        True,
    ),
    "email-EuAll": (
        "https://snap.stanford.edu/data/email-EuAll.txt.gz",
        "Email-EuAll.txt",
        True,
    ),
}


def _ensure_snap_dataset(name: str) -> Path:
    """
    Return the local path for a SNAP dataset, downloading it first if absent.

    Resolution order (stops at first hit):
      1. combo/data/snap/<filename>         — local cache hit, use immediately.
      2. ../paper-2/data/<filename>         — sibling repo copy, copy to cache.
      3. snap.stanford.edu download         — fetch, decompress .gz, cache.
    """
    url, filename, gzipped = _SNAP_REGISTRY[name]
    local_path = SNAP_CACHE_DIR / filename

    # 1. Already cached locally
    if local_path.exists():
        return local_path

    # 2. Copy from paper-2 sibling directory if available (avoids download)
    paper2_path = _PROJECT_ROOT.parent / "paper-2" / "data" / filename
    if paper2_path.exists():
        print(f"  [Data] Copying {name} from paper-2 sibling → {local_path}")
        shutil.copy2(paper2_path, local_path)
        return local_path

    # 3. Download from SNAP
    print(f"  [Data] Downloading {name} ...")
    print(f"         URL: {url}")
    tmp_path = local_path.with_suffix(".gz" if gzipped else ".tmp")
    try:
        urllib.request.urlretrieve(url, tmp_path)
        if gzipped:
            print(f"  [Data] Decompressing {tmp_path.name} ...")
            with gzip.open(tmp_path, "rb") as gz_f, open(local_path, "wb") as out_f:
                shutil.copyfileobj(gz_f, out_f)
            tmp_path.unlink()
        else:
            tmp_path.rename(local_path)
    except Exception as exc:
        if tmp_path.exists():
            tmp_path.unlink()
        raise RuntimeError(
            f"Failed to download {name} from {url}.\n"
            f"You can manually place the file at:\n  {local_path}\n"
            f"Original error: {exc}"
        )

    size_kb = local_path.stat().st_size // 1024
    print(f"  [Data] Saved {name} → {local_path} ({size_kb:,} KB)")
    return local_path


# Canonical dataset names and their roles (for experiment scripts)
PRIMARY_DATASETS   = ["HIV", "p2p-Gnutella", "Wiki-Vote"]     # Exp 1–4, 7–8
LARGE_DATASETS     = ["soc-Epinions1"]                        # Exp 1 if timing OK
SCALABILITY_ONLY   = ["email-EuAll"]                           # Exp 2, 5 only
ALL_DATASETS       = PRIMARY_DATASETS + LARGE_DATASETS + SCALABILITY_ONLY


# ---------------------------------------------------------------------------
# 1. HIV Transmission Network
# ---------------------------------------------------------------------------

def load_hiv_network(
    tsv_path: Path = HIV_TSV_PATH,
    default_weight: float = 0.05,
) -> nx.DiGraph:
    """
    Parse ICPSR #22140 DS0002 and build a directed contact network.

    Returns a DiGraph where each edge (u→v) has 'weight' = default_weight.
    Self-loops and invalid rows are dropped.
    """
    G = nx.DiGraph()

    with open(tsv_path, "r") as f:
        header = f.readline().strip().split("\t")
        try:
            id1_col = header.index("ID1")
            id2_col = header.index("ID2")
        except ValueError:
            raise ValueError("Expected columns 'ID1' and 'ID2' in TSV header.")

        for line in f:
            parts = line.strip().split("\t")
            if len(parts) <= max(id1_col, id2_col):
                continue
            try:
                u = int(parts[id1_col])
                v = int(parts[id2_col])
            except ValueError:
                continue
            if u > 0 and v > 0 and u != v:
                G.add_edge(u, v, weight=default_weight)

    print(f"[HIV] {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
    return G


# ---------------------------------------------------------------------------
# 2. SNAP datasets (p2p-Gnutella, Wiki-Vote, soc-Epinions1, email-EuAll)
# ---------------------------------------------------------------------------

def load_snap_network(
    name: str,
    default_weight: float = 0.05,
    path: Optional[Path] = None,
) -> nx.DiGraph:
    """
    Load a SNAP-format directed graph (tab-separated, '#' comment lines).

    Parameters
    ----------
    name          : One of 'p2p-Gnutella', 'Wiki-Vote', 'soc-Epinions1',
                    'email-EuAll'.
    default_weight: Edge transmission probability assigned to every edge.
    path          : Override the default path (for testing).

    Returns
    -------
    G : nx.DiGraph with 'weight' attribute on every edge.
    """
    if name not in _SNAP_REGISTRY:
        raise ValueError(
            f"Unknown SNAP dataset '{name}'. Known: {list(_SNAP_REGISTRY)}"
        )

    # Resolve path — downloads/copies automatically on first call
    fpath = path or _ensure_snap_dataset(name)

    G = nx.DiGraph()
    with open(fpath, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            try:
                u, v = int(parts[0]), int(parts[1])
            except ValueError:
                continue
            if u != v:
                G.add_edge(u, v, weight=default_weight)

    print(f"[{name}] {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
    return G


# ---------------------------------------------------------------------------
# 3. Unified loader
# ---------------------------------------------------------------------------

def load_dataset(
    name: str,
    default_weight: float = 0.05,
) -> nx.DiGraph:
    """
    Load any supported dataset by name.

    Parameters
    ----------
    name : 'HIV', 'p2p-Gnutella', 'Wiki-Vote', 'soc-Epinions1', 'email-EuAll'

    Returns
    -------
    G : nx.DiGraph
    """
    if name == "HIV":
        return load_hiv_network(default_weight=default_weight)
    elif name in _SNAP_REGISTRY:
        return load_snap_network(name, default_weight=default_weight)
    else:
        raise ValueError(
            f"Unknown dataset '{name}'. Supported: 'HIV', "
            + ", ".join(f"'{k}'" for k in _SNAP_REGISTRY)
        )


# ---------------------------------------------------------------------------
# 4. Initial state seeding
# ---------------------------------------------------------------------------

def seed_initial_states(
    G: nx.DiGraph,
    frac_I: float = 0.02,
    frac_N: float = 0.02,
    seed: int = 42,
) -> Tuple[Set, Set, Set, Set]:
    """
    Fraction-based seeding. frac_I and frac_N control the proportion of
    Infectious and Non-symptomatic nodes at t=0. Remaining nodes = S.

    Returns (initial_S, initial_N, initial_I, initial_R).
    """
    rng = random.Random(seed)
    nodes = list(G.nodes())
    rng.shuffle(nodes)

    n   = len(nodes)
    n_I = max(1, int(frac_I * n))
    n_N = max(1, int(frac_N * n))

    initial_I = set(nodes[:n_I])
    initial_N = set(nodes[n_I:n_I + n_N])
    initial_R: Set = set()
    initial_S = set(nodes[n_I + n_N:])

    print(f"  [Seed] S={len(initial_S)}, N={len(initial_N)}, I={len(initial_I)}, R=0")
    return initial_S, initial_N, initial_I, initial_R


def seed_fixed_nodes(
    G: nx.DiGraph,
    n_infected: int = 5,
    seed: int = 42,
) -> Tuple[Set, Set, Set, Set]:
    """
    Fixed-count seeding. Selects exactly `n_infected` nodes as Infectious (I)
    and places all others in Susceptible (S).  N = R = ∅.

    Seeds are chosen preferentially from nodes that have at least one
    out-neighbour pointing to another node, so that Γ(W) is non-empty
    from iteration 0. If fewer than n_infected such nodes exist, all
    nodes are eligible as a fallback.

    This matches the plan's Section 3.3: "5 randomly chosen infected nodes".
    Used in all experiments unless overridden (e.g., Experiment 8 varies count).

    Returns (initial_S, initial_N, initial_I, initial_R).
    """
    rng = random.Random(seed)
    # Prefer nodes with out-edges so Gamma(W) is non-empty at t=0
    nodes_with_out = [u for u in G.nodes() if G.out_degree(u) > 0]
    pool = nodes_with_out if len(nodes_with_out) >= n_infected else list(G.nodes())
    infected = set(rng.sample(pool, min(n_infected, len(pool))))
    susceptible = set(G.nodes()) - infected

    print(f"  [Seed] S={len(susceptible)}, N=0, I={len(infected)}, R=0")
    return susceptible, set(), infected, set()


# ---------------------------------------------------------------------------
# 5. Transmission probability helpers
# ---------------------------------------------------------------------------

def assign_fixed_weights(G: nx.DiGraph) -> nx.DiGraph:
    """
    Fixed transmission probability: p(u,v) = 1 / avg_degree(G).
    Modifies G in-place and returns it.
    """
    avg_deg = G.number_of_edges() / max(G.number_of_nodes(), 1)
    p = 1.0 / max(avg_deg, 1.0)
    nx.set_edge_attributes(G, p, "weight")
    return G


def assign_variable_weights(G: nx.DiGraph) -> nx.DiGraph:
    """
    Variable transmission probability: p(u,v) = 1 / out_degree(u).
    Nodes with out_degree = 0 get weight = 0.
    Modifies G in-place and returns it.
    """
    for u, v in G.edges():
        deg = G.out_degree(u)
        G[u][v]["weight"] = (1.0 / deg) if deg > 0 else 0.0
    return G


# ---------------------------------------------------------------------------
# 6. Synthetic networks (kept for backward compatibility)
# ---------------------------------------------------------------------------

ER_CONFIGS = {
    "ER-Small":  {"n": 1000, "p": 0.004},
    "ER-Medium": {"n": 3000, "p": 0.0013},
    "ER-Large":  {"n": 5000, "p": 0.0008},
}


def generate_er_network(
    n: int,
    p: float,
    seed: int = 42,
    default_weight: float = 0.05,
) -> nx.DiGraph:
    """Generate a directed Erdős–Rényi graph G(n,p) with uniform edge weights."""
    G = nx.erdos_renyi_graph(n, p, seed=seed, directed=True)
    nx.set_edge_attributes(G, default_weight, "weight")
    print(f"[ER] n={n}, p={p}: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
    return G

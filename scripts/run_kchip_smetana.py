#!/usr/bin/env python3
"""
Compute SMETANA detailed scores for every 2-organism permutation of the kchip
models across all non-M9 JSON media.

Usage:
    python scripts/run_kchip_smetana.py [--output-dir DIR] [--verbose] [--workers N]

Outputs one TSV per medium to the output directory, plus a combined TSV.
"""

import os
import sys
import json
import glob
import time
import argparse
import warnings
import traceback
from itertools import combinations, permutations
from math import inf
from concurrent.futures import ProcessPoolExecutor, as_completed

import cobra
import pandas as pd
from reframed import Environment
from reframed.solvers import set_default_solver

set_default_solver('gurobi')

from smetana.cobrapy import convert_cobrapy_model
from smetana.legacy import Community
from smetana.smetana import sc_score, mu_score, mp_score

# ── paths ────────────────────────────────────────────────────────────────────
MODELS_DIR = "/home/freiburger/Documents/CommScores/data/processed/kchip/models_fixed"
MEDIA_DIR  = "/home/freiburger/Documents/CommScores/data/processed/kchip/media"

COLUMNS = ["community", "medium", "receiver", "donor", "compound",
           "scs", "mus", "mps", "smetana"]


# ── model loading & caching ─────────────────────────────────────────────────

def load_all_models(models_dir, verbose=False):
    """Load all .sbml models via COBRApy and convert to reframed, cached."""
    sbml_files = sorted(glob.glob(os.path.join(models_dir, "*.sbml"))
                        + glob.glob(os.path.join(models_dir, "*.xml")))
    models = {}
    for path in sbml_files:
        name = os.path.basename(path)
        for suffix in (".genome.mdl.sbml", ".genome.mdl.xml", ".sbml", ".xml"):
            if name.endswith(suffix):
                name = name[:-len(suffix)]
                break
        if verbose:
            print(f"  Loading {name} ...", end=" ", flush=True)
        t0 = time.time()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            cm = cobra.io.read_sbml_model(path)
        cm.id = name
        rm = convert_cobrapy_model(cm)
        models[name] = rm
        if verbose:
            print(f"{time.time() - t0:.1f}s  ({len(rm.reactions)} rxns, "
                  f"{len(rm.get_exchange_reactions())} exchange)")
    return models


# ── media loading ────────────────────────────────────────────────────────────

def load_json_media(media_dir):
    """Load all non-M9 JSON media files. Returns {name: {cpd_id: bound}}."""
    media = {}
    for path in sorted(glob.glob(os.path.join(media_dir, "*.json"))):
        name = os.path.basename(path).replace(".json", "")
        if name.startswith("M9_"):
            continue
        with open(path) as f:
            media[name] = json.load(f)
    return media


def make_environment(media_dict, community):
    """Build a reframed Environment from a JSON media dict for a community.

    The JSON maps cpd IDs (e.g. 'cpd00001') to uptake bounds. In the merged
    community model, pool exchange reactions are named R_EX_M_{cpd}_e0_pool.
    """
    merged = community.merged
    pool_exchanges = set(merged.get_exchange_reactions())

    env = Environment()
    for cpd_id, bound in media_dict.items():
        rxn_id = f"R_EX_M_{cpd_id}_e0_pool"
        if rxn_id in pool_exchanges:
            env[rxn_id] = (-float(bound), inf)

    return env


# ── scoring ──────────────────────────────────────────────────────────────────

def run_pair(org1_id, org2_id, models, medium_name, media_dict):
    """Run detailed SMETANA for an unordered pair on a single medium.

    Builds one community and extracts both receiver/donor directions,
    returning rows for both (org1->org2) and (org2->org1).
    """
    comm_id = f"{org1_id}__{org2_id}"
    rm1 = models[org1_id].copy()
    rm2 = models[org2_id].copy()
    community = Community(comm_id, [rm1, rm2], copy_models=False, create_biomass=False)

    env = make_environment(media_dict, community)
    env.apply(community.merged, inplace=True, warning=False)

    # SCS — returns scores for both organisms
    try:
        scs = sc_score(community, verbose=False)
    except Exception:
        return []

    # MUS — returns uptake frequencies for both organisms
    try:
        mus = mu_score(community, verbose=False)
    except Exception:
        return []

    # MPS — returns production scores for both organisms
    try:
        mps = mp_score(community)
    except Exception:
        return []

    # Extract both directions from the single community run
    rows = []
    for receiver, donor in permutations(community.organisms, 2):
        if scs[receiver] is None or mus[receiver] is None or mps[donor] is None:
            continue

        metabolites = set(mus[receiver]) | set(mps[donor])
        for met in sorted(metabolites):
            scs_val = scs[receiver].get(donor, 0)
            mus_val = mus[receiver].get(met, 0)
            mps_val = mps[donor].get(met, 0)
            smt = scs_val * mus_val * mps_val
            rows.append((comm_id, medium_name, receiver, donor, met,
                         scs_val, mus_val, mps_val, smt))

    return rows


# ── worker for multiprocessing ───────────────────────────────────────────────

# Global state for worker processes (populated by initializer)
_worker_models = None
_worker_media = None

def _worker_init(models_dir, media_dir):
    """Initialize per-worker model and media caches."""
    global _worker_models, _worker_media
    _worker_models = load_all_models(models_dir)
    _worker_media = load_json_media(media_dir)

def _worker_run(args):
    """Run a single pair-medium job in a worker process."""
    org1, org2, medium_name = args
    try:
        return run_pair(org1, org2, _worker_models, medium_name, _worker_media[medium_name])
    except Exception as e:
        return []


# ── main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Run SMETANA on all kchip model permutations and media.")
    parser.add_argument("--output-dir", default="kchip_smetana_results",
                        help="Directory for output TSVs (default: kchip_smetana_results)")
    parser.add_argument("--verbose", "-v", action="store_true")
    parser.add_argument("--workers", "-w", type=int, default=1,
                        help="Number of parallel worker processes (default: 1)")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # Load and cache all models in the main process
    print("Loading models...")
    models = load_all_models(MODELS_DIR, verbose=args.verbose)
    model_ids = sorted(models.keys())
    print(f"  {len(model_ids)} models loaded.\n")

    # Load media
    print("Loading media...")
    all_media = load_json_media(MEDIA_DIR)
    media_names = sorted(all_media.keys())
    print(f"  {len(media_names)} non-M9 media loaded.\n")

    # Unordered pairs — each community is built once and both directions extracted
    pairs = list(combinations(model_ids, 2))
    total_jobs = len(pairs) * len(media_names)
    print(f"Running {len(pairs)} unordered pairs x {len(media_names)} media "
          f"= {total_jobs} community-medium jobs.\n")

    all_data = []
    done = 0
    failed = 0
    t_start = time.time()

    if args.workers > 1:
        # ── parallel execution ──
        jobs = [(o1, o2, med) for med in media_names for o1, o2 in pairs]

        with ProcessPoolExecutor(
            max_workers=args.workers,
            initializer=_worker_init,
            initargs=(MODELS_DIR, MEDIA_DIR)
        ) as executor:
            futures = {executor.submit(_worker_run, job): job for job in jobs}
            for future in as_completed(futures):
                done += 1
                job = futures[future]
                try:
                    rows = future.result()
                    all_data.extend(rows)
                except Exception as e:
                    failed += 1
                    print(f"  FAILED: {job[0]}+{job[1]} on {job[2]}: {e}",
                          file=sys.stderr)

                if done % 100 == 0 or done == total_jobs:
                    elapsed = time.time() - t_start
                    rate = done / elapsed if elapsed > 0 else 0
                    eta = (total_jobs - done) / rate if rate > 0 else 0
                    print(f"  Progress: {done}/{total_jobs} "
                          f"({elapsed:.0f}s elapsed, ~{eta:.0f}s remaining, "
                          f"{failed} failures)")

    else:
        # ── sequential execution ──
        for medium_name in media_names:
            media_dict = all_media[medium_name]
            medium_rows = []
            t_medium = time.time()

            for org1, org2 in pairs:
                done += 1
                if args.verbose:
                    print(f"  [{done}/{total_jobs}] {org1} + {org2} on {medium_name}",
                          flush=True)

                try:
                    rows = run_pair(org1, org2, models, medium_name, media_dict)
                    medium_rows.extend(rows)
                except Exception as e:
                    failed += 1
                    print(f"  FAILED: {org1}+{org2} on {medium_name}: {e}",
                          file=sys.stderr)

            # Write per-medium file
            if medium_rows:
                df = pd.DataFrame(medium_rows, columns=COLUMNS)
                outpath = os.path.join(args.output_dir, f"{medium_name}_detailed.tsv")
                df.to_csv(outpath, sep="\t", index=False)

            all_data.extend(medium_rows)
            elapsed = time.time() - t_start
            rate = done / elapsed if elapsed > 0 else 0
            eta = (total_jobs - done) / rate if rate > 0 else 0
            print(f"  {medium_name}: {len(medium_rows)} rows in {time.time()-t_medium:.0f}s "
                  f"({done}/{total_jobs}, ~{eta:.0f}s remaining, {failed} failures)")

    # Write per-medium files (parallel mode — group after all done)
    if args.workers > 1 and all_data:
        df_all = pd.DataFrame(all_data, columns=COLUMNS)
        for medium_name in media_names:
            df_med = df_all[df_all["medium"] == medium_name]
            if len(df_med) > 0:
                outpath = os.path.join(args.output_dir, f"{medium_name}_detailed.tsv")
                df_med.to_csv(outpath, sep="\t", index=False)

    # Write combined file
    if all_data:
        df_all = pd.DataFrame(all_data, columns=COLUMNS)
        combined_path = os.path.join(args.output_dir, "all_detailed.csv")
        df_all.to_csv(combined_path, index=False)
        print(f"\nCombined results: {combined_path} ({len(df_all)} rows)")

    elapsed = time.time() - t_start
    print(f"\nDone. {done} jobs completed in {elapsed:.0f}s, {failed} failures.")


if __name__ == "__main__":
    main()

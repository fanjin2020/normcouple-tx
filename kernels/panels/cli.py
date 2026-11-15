"""Dependency-free CLI, including progress inspection."""
import argparse
import os
from pathlib import Path


def parser():
    p = argparse.ArgumentParser(description="Paper 3: 50 pre-fixed randomized stress subsets, three paired direct effects; cached inputs only.")
    p.add_argument("--project-root", type=Path, default=Path("/home/research/git/repository/paper1_public_evidence"))
    p.add_argument("--output", type=Path, default=Path("/home/research/git/repository/paper3_random_marker_results_v1"))
    p.add_argument("--workers", type=int, default=36 if os.name != "nt" else 1)
    p.add_argument("--repetitions", type=int, default=200)
    p.add_argument("--subsets", type=int, default=50, help="Fixed at 50 for research; adjustable only with --synthetic")
    p.add_argument("--batch-size", type=int, default=5)
    p.add_argument("--subset-batch-size", type=int, default=5)
    p.add_argument("--r-timeout", type=int, default=7200)
    p.add_argument("--rscript", default="Rscript")
    for name in ("libd-dir", "seaad-dir", "gene-map", "external-cache", "input-bundle", "r-library", "old-results"):
        p.add_argument("--" + name, type=Path)
    for name in ("setup", "offline", "status", "verify", "pilot-only", "synthetic", "summarize-only"):
        p.add_argument("--" + name, action="store_true")
    return p

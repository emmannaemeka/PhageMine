#!/usr/bin/env python3
"""Outcome-independent planning scenarios for the validation panel.

This script uses only declared planning assumptions. It never reads genomes,
PhageMine outputs, benchmark labels, or validation outcomes.
"""
from __future__ import annotations

import argparse
import json
import math
from statistics import NormalDist


def mean_precision_n(sd: float, half_width: float, confidence: float = 0.95) -> int:
    z = NormalDist().inv_cdf((1 + confidence) / 2)
    return math.ceil((z * sd / half_width) ** 2)


def paired_power_n(sd: float, difference: float, power: float = 0.80, alpha: float = 0.05) -> int:
    z_alpha = NormalDist().inv_cdf(1 - alpha / 2)
    z_beta = NormalDist().inv_cdf(power)
    return math.ceil(((z_alpha + z_beta) * sd / difference) ** 2)


def scenarios() -> dict:
    precision = [
        {"sd": sd, "half_width": h, "n": mean_precision_n(sd, h)}
        for sd in (0.10, 0.15, 0.20)
        for h in (0.05, 0.075, 0.10)
    ]
    paired = [
        {"paired_sd": sd, "target_difference": d, "power": p,
         "n": paired_power_n(sd, d, p)}
        for sd in (0.10, 0.15, 0.20)
        for d in (0.10, 0.15, 0.20)
        for p in (0.80, 0.90)
    ]
    return {
        "assumptions_are_planning_only": True,
        "confidence": 0.95,
        "alpha": 0.05,
        "precision_scenarios": precision,
        "paired_scenarios": paired,
        "operational_targets": {
            "minimum_analyzable_panel": 32,
            "preferred_target_panel": 64,
            "reserve_pool": 16,
            "rationale": "32 covers the 80% paired-difference scenario at SD=0.20 and difference=0.10; 64 covers the 0.20 SD, 0.05 half-width precision scenario and exceeds the corresponding 90% paired scenario. These are not observed estimates."
        }
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=str)
    args = parser.parse_args()
    payload = scenarios()
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.output:
        with open(args.output, "w", encoding="utf-8") as handle:
            handle.write(rendered)
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()

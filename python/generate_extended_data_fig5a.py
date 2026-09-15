"""Reproduce regressions and plots for Extended Data Figure 5a."""

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from fit_gS_untreated import fit_gs
from tumor_data import load_tumor_volumes


RESULTS = Path(__file__).resolve().parents[1] / "results"


def main():
    RESULTS.mkdir(exist_ok=True)
    time, untreated, _ = load_tumor_volumes()
    estimates = []
    fig, axes = plt.subplots(2, 4, figsize=(11, 5.5), sharex=True)

    for tumour, (volume, axis) in enumerate(zip(untreated, axes.flat), start=1):
        valid = np.isfinite(time) & np.isfinite(volume) & (volume > 0)
        t = time[valid]
        v = volume[valid]
        gs, r_squared = fit_gs(t, v)
        intercept = np.mean(np.log(v)) - gs * np.mean(t)
        t_fit = np.linspace(t.min(), t.max(), 200)
        v_fit = np.exp(intercept + gs * t_fit)

        estimates.append((tumour, len(t), gs, r_squared))
        axis.scatter(t, v, s=22, color="#333333", zorder=2)
        axis.plot(t_fit, v_fit, color="#C43C39", linewidth=1.8)
        axis.set_yscale("log")
        axis.set_title(f"Tumour {tumour}: $g_S$={gs:.3f}, $R^2$={r_squared:.3f}", fontsize=9)
        axis.set_xlabel("Time (days)")
        axis.set_ylabel("Tumour volume (mm$^3$)")

    fig.tight_layout()
    fig.savefig(RESULTS / "Extended_Data_Fig_5a.pdf", bbox_inches="tight")
    fig.savefig(RESULTS / "Extended_Data_Fig_5a.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    with (RESULTS / "gS_estimates.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["tumour", "n_points", "gS_day-1", "R_squared"])
        writer.writerows(estimates)

    rates = np.array([row[2] for row in estimates])
    print(f"gS mean ± SD (n=8): {rates.mean():.6f} ± {rates.std(ddof=1):.6f} day^-1")


if __name__ == "__main__":
    main()

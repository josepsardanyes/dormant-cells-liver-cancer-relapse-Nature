"""Reproduce relapse regressions and plots for Extended Data Figure 5b."""

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import linregress

from fit_gR_relapse import fit_gr
from tumor_data import load_tumor_volumes


RESULTS = Path(__file__).resolve().parents[1] / "results"
RELAPSING_TUMOURS = (1, 2, 3, 4, 5, 6, 9, 10, 12)


def relapse_fit(time, volume):
    valid = np.isfinite(time) & np.isfinite(volume) & (volume > 0)
    t = time[valid]
    v = volume[valid]
    post = np.flatnonzero(t >= 26)
    nadir_index = post[np.argmin(v[post])]
    t_relapse = t[nadir_index:]
    v_relapse = v[nadir_index:]
    gr, r_squared = fit_gr(t, v, start_day=26)
    regression = linregress(t_relapse, np.log(v_relapse))
    prediction = np.exp(regression.intercept + gr * t_relapse)
    nrmse = 100 * np.sqrt(np.mean((prediction - v_relapse) ** 2)) / np.ptp(v_relapse)
    return t, v, t_relapse, v_relapse, gr, r_squared, nrmse, regression.intercept


def main():
    RESULTS.mkdir(exist_ok=True)
    time, _, treated = load_tumor_volumes()
    estimates = []
    fig, axes = plt.subplots(3, 3, figsize=(9, 8), sharex=True)

    for tumour, axis in zip(RELAPSING_TUMOURS, axes.flat):
        t, v, tr, vr, gr, r2, nrmse, intercept = relapse_fit(time, treated[tumour - 1])
        estimates.append((tumour, tr[0], len(tr), gr, r2, nrmse))
        t_fit = np.linspace(tr.min(), tr.max(), 200)
        axis.scatter(t, v, s=18, color="#A0A0A0", label="All observations")
        axis.scatter(tr, vr, s=24, color="#333333", label="Relapse phase", zorder=3)
        axis.plot(t_fit, np.exp(intercept + gr * t_fit), color="#2B6CB0", linewidth=1.8)
        axis.set_yscale("log")
        axis.set_title(f"Tumour {tumour}: $g_R$={gr:.3f}, $R^2$={r2:.3f}", fontsize=9)
        axis.set_xlabel("Time (days)")
        axis.set_ylabel("Tumour volume (mm$^3$)")

    fig.tight_layout()
    fig.savefig(RESULTS / "Extended_Data_Fig_5b.pdf", bbox_inches="tight")
    fig.savefig(RESULTS / "Extended_Data_Fig_5b.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    with (RESULTS / "gR_estimates.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["tumour", "start_day", "n_points", "gR_day-1", "R_squared", "NRMSE_percent"])
        writer.writerows(estimates)

    all_rates = np.array([row[3] for row in estimates])
    primary_rates = np.array([row[3] for row in estimates if row[2] > 3])
    print(f"gR mean ± SD (n=9): {all_rates.mean():.6f} ± {all_rates.std(ddof=1):.6f} day^-1")
    print(
        "gR mean ± SD for tumours with >3 points "
        f"(n=6): {primary_rates.mean():.6f} ± {primary_rates.std(ddof=1):.6f} day^-1"
    )


if __name__ == "__main__":
    main()

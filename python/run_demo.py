"""Run the two regression functions on simulated demonstration data."""

import csv
from pathlib import Path

import numpy as np

from fit_gR_relapse import fit_gr
from fit_gS_untreated import fit_gs


def load_demo_data():
    path = Path(__file__).resolve().parents[1] / "data" / "demo_regression_data.csv"
    rows = list(csv.DictReader(path.open(encoding="utf-8")))
    time = np.array([float(row["time_days"]) for row in rows])
    untreated = np.array([float(row["untreated_volume"]) for row in rows])
    relapse = np.array([float(row["relapse_volume"]) for row in rows])
    return time, untreated, relapse


if __name__ == "__main__":
    time, untreated_volume, relapse_volume = load_demo_data()
    gs, gs_r2 = fit_gs(time, untreated_volume)
    gr, gr_r2 = fit_gr(time, relapse_volume, start_day=20)
    print(f"Untreated fit: gS = {gs:.6f} day^-1, R^2 = {gs_r2:.6f}")
    print(f"Relapse fit:   gR = {gr:.6f} day^-1, R^2 = {gr_r2:.6f}")

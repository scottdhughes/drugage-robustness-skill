"""Plot helpers."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def plot_binary_heatmap(frame: pd.DataFrame, path: Path, title: str) -> None:
    if frame.shape[0] == 0 or frame.shape[1] == 0:
        frame = pd.DataFrame({"not_available": [0]}, index=["no_data"])
    data = frame.to_numpy(dtype=float)
    rows = list(frame.index)
    cols = list(frame.columns)
    height = max(3, 0.45 * len(rows) + 1.8)
    width = max(5, 0.9 * len(cols) + 1.6)
    fig, ax = plt.subplots(figsize=(width, height))
    ax.imshow(data, aspect="auto", cmap="Blues", vmin=0, vmax=1)
    ax.set_xticks(np.arange(len(cols)))
    ax.set_xticklabels(cols, rotation=45, ha="right")
    ax.set_yticks(np.arange(len(rows)))
    ax.set_yticklabels(rows)
    ax.set_title(title)
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            ax.text(j, i, f"{int(data[i, j])}", ha="center", va="center", color="black", fontsize=8)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=200)
    plt.close(fig)


def plot_null_separation(null_values: list[float], observed: float, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.hist(null_values, bins=min(20, max(8, len(null_values) // 4)), color="#9ecae1", edgecolor="#1f77b4")
    ax.axvline(observed, color="#d62728", linewidth=2, label=f"Observed = {observed:.4f}")
    ax.set_xlabel("Top-10 Mean Robustness Score")
    ax.set_ylabel("Null Rerun Count")
    ax.set_title("Observed Top-10 Robustness vs Empirical Null")
    ax.legend()
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=200)
    plt.close(fig)

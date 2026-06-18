"""Build a before/after comparison collage for raw vs normalized figures.

The collage is arranged as:
UMAP before | UMAP after
Dendrogram before | Dendrogram after
Boxplot before | Boxplot after

Default canvas is A4 landscape-sized so it can be dropped into a manuscript figure.
"""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import matplotlib.image as mpimg


def copy_if_exists(source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)


def load_image(path: Path):
    if not path.exists():
        raise FileNotFoundError(f"Missing image: {path}")
    return mpimg.imread(path)


def place_image(ax, image_path: Path, title: str) -> None:
    ax.imshow(load_image(image_path))
    ax.set_title(title, fontsize=11, weight="bold")
    ax.axis("off")


def main() -> None:
    parser = argparse.ArgumentParser(description="Create a before/after collage from existing plot PNGs.")
    parser.add_argument("--root", type=Path, required=True, help="Root visualization folder for GSE142279")
    parser.add_argument("--outdir", type=Path, required=True, help="Comparison output folder")
    parser.add_argument("--page", choices=["A4_landscape", "A4_portrait", "letter_landscape", "letter_portrait"], default="A4_landscape", help="Canvas size preset")
    args = parser.parse_args()

    presets = {
        "A4_landscape": (11.69, 8.27),
        "A4_portrait": (8.27, 11.69),
        "letter_landscape": (11.0, 8.5),
        "letter_portrait": (8.5, 11.0),
    }
    figsize = presets[args.page]

    raw_dir = args.root / "raw_pca_umap"
    norm_dir = args.root / "normalized_plots" / "plots_normalized"

    before_dir = args.outdir / "before"
    after_dir = args.outdir / "after"
    before_dir.mkdir(parents=True, exist_ok=True)
    after_dir.mkdir(parents=True, exist_ok=True)

    # Copy the source figures into the comparison folder for easy reuse.
    copy_if_exists(raw_dir / "plots_raw_counts" / "raw_counts_umap.png", before_dir / "umap_before.png")
    copy_if_exists(raw_dir / "plots_raw_counts" / "raw_counts_umap.pdf", before_dir / "umap_before.pdf")
    copy_if_exists(norm_dir / "normalized_umap.png", after_dir / "umap_after.png")
    copy_if_exists(norm_dir / "normalized_umap.pdf", after_dir / "umap_after.pdf")

    copy_if_exists(raw_dir / "plots_dendrogram" / "sample_dendrogram.png", before_dir / "dendrogram_before.png")
    copy_if_exists(raw_dir / "plots_dendrogram" / "sample_dendrogram.pdf", before_dir / "dendrogram_before.pdf")
    copy_if_exists(norm_dir / "normalized_sample_dendrogram.png", after_dir / "dendrogram_after.png")
    copy_if_exists(norm_dir / "normalized_sample_dendrogram.pdf", after_dir / "dendrogram_after.pdf")

    copy_if_exists(raw_dir / "plots_boxplot" / "boxplot_per_sample.png", before_dir / "boxplot_before.png")
    copy_if_exists(raw_dir / "plots_boxplot" / "boxplot_per_sample.pdf", before_dir / "boxplot_before.pdf")
    copy_if_exists(norm_dir / "normalized_boxplot_per_sample.png", after_dir / "boxplot_after.png")
    copy_if_exists(norm_dir / "normalized_boxplot_per_sample.pdf", after_dir / "boxplot_after.pdf")

    fig, axes = plt.subplots(3, 2, figsize=figsize, constrained_layout=True)

    place_image(axes[0, 0], before_dir / "umap_before.png", "UMAP before normalization")
    place_image(axes[0, 1], after_dir / "umap_after.png", "UMAP after normalization")

    place_image(axes[1, 0], before_dir / "dendrogram_before.png", "Dendrogram before normalization")
    place_image(axes[1, 1], after_dir / "dendrogram_after.png", "Dendrogram after normalization")

    place_image(axes[2, 0], before_dir / "boxplot_before.png", "Boxplot before normalization")
    place_image(axes[2, 1], after_dir / "boxplot_after.png", "Boxplot after normalization")

    fig.suptitle("Raw vs Normalized Sample-Level Plots", fontsize=14, weight="bold")
    fig.savefig(args.outdir / "comparison_collage.png", dpi=300, bbox_inches="tight")
    fig.savefig(args.outdir / "comparison_collage.pdf", bbox_inches="tight")
    plt.close(fig)

    print(f"Comparison collage saved to {args.outdir}")


if __name__ == "__main__":
    main()

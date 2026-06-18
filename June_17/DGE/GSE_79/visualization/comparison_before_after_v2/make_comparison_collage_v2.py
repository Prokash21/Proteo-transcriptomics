"""Build a before/after comparison collage with fixed panel sizes.

This version pads each plot image to the same canvas before assembling the collage,
so the rows and columns line up symmetrically on an A4/letter page.

Layout:
UMAP before | UMAP after
Dendrogram before | Dendrogram after
Boxplot before | Boxplot after
"""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from PIL import Image, ImageOps


PAGE_PRESETS = {
    "A4_landscape": (11.69, 8.27),
    "A4_portrait": (8.27, 11.69),
    "letter_landscape": (11.0, 8.5),
    "letter_portrait": (8.5, 11.0),
}


def copy_if_exists(source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)


def load_image(path: Path) -> Image.Image:
    if not path.exists():
        raise FileNotFoundError(f"Missing image: {path}")
    return Image.open(path).convert("RGB")


def pad_to_canvas(image: Image.Image, canvas_size: tuple[int, int]) -> Image.Image:
    canvas = Image.new("RGB", canvas_size, "white")
    fitted = ImageOps.contain(image, canvas_size, method=Image.Resampling.LANCZOS)
    offset_x = (canvas_size[0] - fitted.size[0]) // 2
    offset_y = (canvas_size[1] - fitted.size[1]) // 2
    canvas.paste(fitted, (offset_x, offset_y))
    return canvas


def save_padded_copy(source: Path, destination: Path, canvas_size: tuple[int, int]) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    padded = pad_to_canvas(load_image(source), canvas_size)
    padded.save(destination)


def place_image(ax, image_path: Path, title: str) -> None:
    ax.imshow(Image.open(image_path))
    ax.set_title(title, fontsize=11, weight="bold")
    ax.axis("off")


def main() -> None:
    parser = argparse.ArgumentParser(description="Create a before/after collage from existing plot PNGs using padded panels.")
    parser.add_argument("--root", type=Path, required=True, help="Root visualization folder for GSE142279")
    parser.add_argument("--outdir", type=Path, required=True, help="Comparison v2 output folder")
    parser.add_argument("--page", choices=list(PAGE_PRESETS.keys()), default="A4_landscape", help="Page preset for the collage")
    parser.add_argument("--dpi", type=int, default=300, help="DPI for padded panels and final figure")
    args = parser.parse_args()

    raw_dir = args.root / "raw_pca_umap"
    norm_dir = args.root / "normalized_plots" / "plots_normalized"

    before_dir = args.outdir / "before"
    after_dir = args.outdir / "after"
    panels_dir = args.outdir / "panels"
    before_dir.mkdir(parents=True, exist_ok=True)
    after_dir.mkdir(parents=True, exist_ok=True)
    panels_dir.mkdir(parents=True, exist_ok=True)

    raw_sources = {
        "umap": (raw_dir / "plots_raw_counts" / "raw_counts_umap.png", raw_dir / "plots_raw_counts" / "raw_counts_umap.pdf"),
        "dendrogram": (raw_dir / "plots_dendrogram" / "sample_dendrogram.png", raw_dir / "plots_dendrogram" / "sample_dendrogram.pdf"),
        "boxplot": (raw_dir / "plots_boxplot" / "boxplot_per_sample.png", raw_dir / "plots_boxplot" / "boxplot_per_sample.pdf"),
    }
    norm_sources = {
        "umap": (norm_dir / "normalized_umap.png", norm_dir / "normalized_umap.pdf"),
        "dendrogram": (norm_dir / "normalized_sample_dendrogram.png", norm_dir / "normalized_sample_dendrogram.pdf"),
        "boxplot": (norm_dir / "normalized_boxplot_per_sample.png", norm_dir / "normalized_boxplot_per_sample.pdf"),
    }

    for key, (png_path, pdf_path) in raw_sources.items():
        copy_if_exists(png_path, before_dir / f"{key}_before.png")
        copy_if_exists(pdf_path, before_dir / f"{key}_before.pdf")
    for key, (png_path, pdf_path) in norm_sources.items():
        copy_if_exists(png_path, after_dir / f"{key}_after.png")
        copy_if_exists(pdf_path, after_dir / f"{key}_after.pdf")

    panel_paths = {
        "umap_before": panels_dir / "umap_before.png",
        "umap_after": panels_dir / "umap_after.png",
        "dendrogram_before": panels_dir / "dendrogram_before.png",
        "dendrogram_after": panels_dir / "dendrogram_after.png",
        "boxplot_before": panels_dir / "boxplot_before.png",
        "boxplot_after": panels_dir / "boxplot_after.png",
    }

    all_pngs = [
        raw_sources["umap"][0],
        norm_sources["umap"][0],
        raw_sources["dendrogram"][0],
        norm_sources["dendrogram"][0],
        raw_sources["boxplot"][0],
        norm_sources["boxplot"][0],
    ]
    images = [load_image(path) for path in all_pngs]
    max_width = max(image.size[0] for image in images)
    max_height = max(image.size[1] for image in images)
    canvas_size = (max_width, max_height)

    save_padded_copy(raw_sources["umap"][0], panel_paths["umap_before"], canvas_size)
    save_padded_copy(norm_sources["umap"][0], panel_paths["umap_after"], canvas_size)
    save_padded_copy(raw_sources["dendrogram"][0], panel_paths["dendrogram_before"], canvas_size)
    save_padded_copy(norm_sources["dendrogram"][0], panel_paths["dendrogram_after"], canvas_size)
    save_padded_copy(raw_sources["boxplot"][0], panel_paths["boxplot_before"], canvas_size)
    save_padded_copy(norm_sources["boxplot"][0], panel_paths["boxplot_after"], canvas_size)

    fig_size = PAGE_PRESETS[args.page]
    fig, axes = plt.subplots(3, 2, figsize=fig_size, constrained_layout=True)

    place_image(axes[0, 0], panel_paths["umap_before"], "UMAP before normalization")
    place_image(axes[0, 1], panel_paths["umap_after"], "UMAP after normalization")

    place_image(axes[1, 0], panel_paths["dendrogram_before"], "Dendrogram before normalization")
    place_image(axes[1, 1], panel_paths["dendrogram_after"], "Dendrogram after normalization")

    place_image(axes[2, 0], panel_paths["boxplot_before"], "Boxplot before normalization")
    place_image(axes[2, 1], panel_paths["boxplot_after"], "Boxplot after normalization")

    fig.suptitle("Raw vs Normalized Sample-Level Plots", fontsize=14, weight="bold")
    fig.savefig(args.outdir / "comparison_collage_v2.png", dpi=args.dpi, bbox_inches="tight")
    fig.savefig(args.outdir / "comparison_collage_v2.pdf", bbox_inches="tight")
    plt.close(fig)

    print(f"Comparison collage saved to {args.outdir}")


if __name__ == "__main__":
    main()

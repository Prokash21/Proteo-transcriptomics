"""
Generate sample-level visualizations (PCA, UMAP, hierarchical tree, boxplots)
for an input matrix and a derived matrix.

Usage example:
python gse79_plot_samples.py \
  --counts path/to/counts.csv \
  --meta path/to/meta.csv \
  --outdir ./plots

The script expects a matrix CSV with genes as rows and samples as columns (or will try to detect orientation).
Metadata should have a column matching sample IDs and an optional `group` or `tissue`
column for coloring (this script prefers `tissue`).
"""
import argparse
import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
import umap
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import pdist


def load_counts(path):
    df = pd.read_csv(path, index_col=0)
    # ensure numeric matrix
    if df.select_dtypes(include=[np.number]).shape[1] < df.shape[1]:
        # maybe samples are rows, genes columns
        df = df.transpose()
    return df


def load_meta(path):
    return pd.read_csv(path, index_col=0)


def cpm_normalize(counts):
    col_sums = counts.sum(axis=0)
    cpm = counts.divide(col_sums, axis=1) * 1e6
    return cpm


def ensure_outdir(path):
    os.makedirs(path, exist_ok=True)


def plot_pca(mat, meta, outpath, title="PCA"):
    pca = PCA(n_components=2)
    X = pca.fit_transform(mat.T)
    pc = pd.DataFrame(X, index=mat.columns, columns=["PC1", "PC2"])
    var = pca.explained_variance_ratio_ * 100
    fig, ax = plt.subplots(figsize=(6,5))
    joined = pc.join(meta)
    if 'tissue' in meta.columns:
        sns.scatterplot(x='PC1', y='PC2', data=joined, hue='tissue', s=80, ax=ax)
        ax.legend(title='Tissue')
    elif 'group' in meta.columns:
        sns.scatterplot(x='PC1', y='PC2', data=joined, hue='group', s=80, ax=ax)
        ax.legend(title='Group')
    else:
        sns.scatterplot(x='PC1', y='PC2', data=pc, s=80, ax=ax)
    ax.set_xlabel(f"PC1 ({var[0]:.1f}%)")
    ax.set_ylabel(f"PC2 ({var[1]:.1f}%)")
    ax.set_title(title)
    plt.tight_layout()
    fig.savefig(outpath)
    plt.close(fig)


def plot_umap(mat, meta, outpath, title="UMAP"):
    reducer = umap.UMAP(n_neighbors=15, min_dist=0.1, random_state=0)
    embedding = reducer.fit_transform(mat.T)
    emb = pd.DataFrame(embedding, index=mat.columns, columns=["UMAP1","UMAP2"])
    fig, ax = plt.subplots(figsize=(6,5))
    joined = emb.join(meta)
    if 'tissue' in meta.columns:
        sns.scatterplot(x='UMAP1', y='UMAP2', data=joined, hue='tissue', s=80, ax=ax)
        ax.legend(title='Tissue')
    elif 'group' in meta.columns:
        sns.scatterplot(x='UMAP1', y='UMAP2', data=joined, hue='group', s=80, ax=ax)
        ax.legend(title='Group')
    else:
        sns.scatterplot(x='UMAP1', y='UMAP2', data=emb, s=80, ax=ax)
    ax.set_title(title)
    plt.tight_layout()
    fig.savefig(outpath)
    plt.close(fig)


def plot_dendrogram(mat, meta, outpath, title="Hierarchical Clustering"):
    # use correlation distance
    dist = pdist(mat.T, metric='euclidean')
    Z = linkage(dist, method='average')
    fig, ax = plt.subplots(figsize=(8,4))
    dendrogram(Z, labels=mat.columns, leaf_rotation=90, ax=ax)

    # Color sample labels with exactly two colors using tissue metadata.
    if 'tissue' in meta.columns:
        tissue_colors = {'tumor': '#d62728', 'normal': '#1f77b4'}
        for tick in ax.get_xmajorticklabels():
            sample = tick.get_text()
            tissue = str(meta.loc[sample, 'tissue']).strip().lower() if sample in meta.index else ''
            tick.set_color(tissue_colors.get(tissue, '#333333'))

        handles = [
            plt.Line2D([0], [0], color='#d62728', lw=2, label='tumor'),
            plt.Line2D([0], [0], color='#1f77b4', lw=2, label='normal')
        ]
        ax.legend(handles=handles, title='Tissue', loc='upper right')

    ax.set_title(title)
    plt.tight_layout()
    fig.savefig(outpath)
    plt.close(fig)


def plot_boxplots(mat, outpath, title="Boxplots (per sample)"):
    # mat: genes x samples
    logmat = np.log2(mat + 1)
    df = logmat
    fig, ax = plt.subplots(figsize=(max(6, 0.3*df.shape[1]),4))
    sns.boxplot(data=df, orient='v', color='lightblue')
    ax.set_xticks(range(len(df.columns)))
    ax.set_xticklabels(df.columns, rotation=90)
    ax.set_title(title)
    plt.tight_layout()
    fig.savefig(outpath)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--counts', required=True, help='Input matrix CSV (genes x samples or samples x genes)')
    parser.add_argument('--after-matrix', default=None, help='Optional matrix CSV to use for "after" plots directly')
    parser.add_argument('--meta', required=True, help='Sample metadata CSV with sample IDs as index or a column named "sample"')
    parser.add_argument('--outdir', default='plots', help='Output directory for plots')
    parser.add_argument('--input-status', choices=['raw', 'normalized'], default='raw', help='Whether input matrix is raw counts or already normalized')
    parser.add_argument('--normalization', choices=['cpm','log'], default='cpm', help='Normalization method for "after" plots')
    args = parser.parse_args()

    counts = load_counts(args.counts)
    try:
        meta = load_meta(args.meta)
    except Exception:
        # try reading without index
        meta = pd.read_csv(args.meta)
        if 'sample' in meta.columns:
            meta = meta.set_index('sample')

    # align sample names
    samples = [s for s in counts.columns if s in meta.index]
    if len(samples) == 0:
        print('No sample names matched between counts columns and metadata index.', file=sys.stderr)
    else:
        counts = counts.loc[:, samples]
        meta = meta.loc[samples]

    ensure_outdir(args.outdir)

    # Input matrix plots (titles depend on raw vs normalized input)
    before_mat = counts
    if args.input_status == 'normalized':
        input_label = 'provided normalized matrix'
    else:
        input_label = 'raw counts'

    plot_pca(before_mat, meta, os.path.join(args.outdir, 'before_pca.png'), title=f'PCA ({input_label})')
    plot_umap(before_mat, meta, os.path.join(args.outdir, 'before_umap.png'), title=f'UMAP ({input_label})')
    plot_dendrogram(before_mat, meta, os.path.join(args.outdir, 'before_dendrogram.png'), title=f'Dendrogram ({input_label})')
    plot_boxplots(before_mat, os.path.join(args.outdir, 'before_boxplot.png'), title=f'Boxplots ({input_label}, log2)')

    # After normalization
    if args.after_matrix:
        after = load_counts(args.after_matrix)
        after = after.loc[:, [s for s in samples if s in after.columns]]
    elif args.normalization == 'cpm':
        after = cpm_normalize(counts)
    else:
        after = np.log2(counts + 1)

    norm_label = f'after {args.normalization.upper()} transform'
    plot_pca(after, meta, os.path.join(args.outdir, 'after_pca.png'), title=f'PCA ({norm_label})')
    plot_umap(after, meta, os.path.join(args.outdir, 'after_umap.png'), title=f'UMAP ({norm_label})')
    plot_dendrogram(after, meta, os.path.join(args.outdir, 'after_dendrogram.png'), title=f'Dendrogram ({norm_label})')
    plot_boxplots(after, os.path.join(args.outdir, 'after_boxplot.png'), title=f'Boxplots ({norm_label})')


if __name__ == '__main__':
    main()

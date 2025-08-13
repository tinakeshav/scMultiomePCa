###############################################################################
# plot_unfiltered_multiome_qc.py
# Script to plot unfiltered multiome QC metrics.
#
# Usage: python3 plot_unfiltered_multiome_qc.py --metadatafiles ... [other args]
############

import argparse
import logging
import os
from typing import List, Optional

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

LABEL_SIZE = 8
mpl.rcParams['xtick.labelsize'] = LABEL_SIZE

def parse_args() -> argparse.Namespace:
    """
    Parse command-line arguments for the QC plotting script.
    """
    parser = argparse.ArgumentParser(
        prog='QCPlots',
        description='Plots QC from unfiltered metadata'
    )
    parser.add_argument(
        '--metadatafiles',
        nargs='+',
        required=True,
        help='List of metadata CSV files to process.'
    )
    parser.add_argument(
        '--output_csv',
        default='data/QC_percentiles_samples_22032023.csv',
        help='Output CSV file for QC percentiles.'
    )
    parser.add_argument(
        '--output_figname',
        default='figures/QC_plots_samples_22032023.png',
        help='Output filename for the QC figure.'
    )
    return parser.parse_args()

def add_log_transformed_features(df: pd.DataFrame, features: List[str]) -> pd.DataFrame:
    """
    Add log1p-transformed columns for specified features to the DataFrame.
    """
    for feature in features:
        if feature in df.columns:
            log_col = f'{feature}_log1p'
            df[log_col] = np.log1p(df[feature])
            logging.debug(f'Added log1p column: {log_col}')
        else:
            logging.warning(f'Feature {feature} not found in DataFrame columns.')
    return df

def remove_infs_and_nans(df: pd.DataFrame) -> pd.DataFrame:
    """
    Replace inf/-inf with NaN and drop rows where all values are NaN.
    """
    df_clean = df.replace([np.inf, -np.inf], np.nan)
    dropped = df_clean[df_clean.isna().any(axis=1)]
    if not dropped.empty:
        logging.info(f'Dropping rows with NaN values: {dropped.index.tolist()}')
    return df_clean.dropna(how="all")

def plot_violin_box_swarm(
    df: pd.DataFrame,
    features: List[str],
    group_var: str,
    nrows: int,
    title: str,
    figsize: tuple,
    filename: str,
    plot_dots: bool = False
) -> None:
    """
    Plot multiple violin/box/swarm plots for the given features grouped by group_var.
    """
    plt.figure(figsize=figsize)
    for idx, feature in enumerate(features, 1):
        plt.subplot(nrows, 1, idx)
        sns.violinplot(x=group_var, y=feature, data=df, scale="count", inner="quartile")
        sns.boxplot(x=group_var, y=feature, data=df, showfliers=False, showcaps=False, boxprops={'facecolor':'None'})
        if plot_dots:
            sns.swarmplot(x=group_var, y=feature, data=df, color="white", edgecolor="gray", size=1)
    plt.xticks(rotation=90)
    plt.suptitle(title)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(filename, dpi=600)
    plt.close()
    logging.info(f'Figure saved to {filename}')

def load_and_concat_metadata(files: List[str]) -> pd.DataFrame:
    """
    Load and concatenate metadata CSV files into a single DataFrame.
    """
    dataframes = []
    for f in files:
        try:
            df = pd.read_csv(f)
            dataframes.append(df)
            logging.info(f'Loaded {f} with shape {df.shape}')
        except Exception as e:
            logging.error(f'Error reading {f}: {e}')
    if not dataframes:
        raise ValueError('No metadata files could be loaded.')
    return pd.concat(dataframes, axis=0, ignore_index=True)

def main():
    args = parse_args()
    metadata = load_and_concat_metadata(args.metadatafiles)
    logging.info(f'Combined metadata shape: {metadata.shape}')
    logging.info(f'Columns: {metadata.columns.tolist()}')

    features_to_log = ["nCount_ATAC", "nFeature_ATAC", "nCount_RNA", "nFeature_RNA"]
    log_feature_names = [f"{feature}_log1p" for feature in features_to_log]
    metadata = add_log_transformed_features(metadata, features_to_log)

    other_feature_names = [
        "TSS.enrichment", "nucleosome_signal", "frip", "percent.mt"
    ] + log_feature_names

    missing_features = [f for f in other_feature_names if f not in metadata.columns]
    if missing_features:
        logging.warning(f'Missing features in metadata: {missing_features}')
    features_to_plot = [f for f in other_feature_names if f in metadata.columns]

    metadata_to_plot = remove_infs_and_nans(metadata[features_to_plot].copy())
    if 'name' in metadata.columns:
        metadata_to_plot['name'] = metadata['name']
    else:
        raise KeyError("'name' column is required in metadata for grouping.")

    plot_violin_box_swarm(
        df=metadata_to_plot,
        features=features_to_plot,
        group_var='name',
        nrows=len(features_to_plot),
        title='TNBC Primaries Core Processed',
        figsize=(10, 30),
        filename=args.output_figname
    )

    # Save quantiles for further inspection
    percentiles = metadata_to_plot.groupby('name').describe(percentiles=[0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.8, 0.9, 0.95, 0.975, 1.0])
    percentiles.to_csv(args.output_csv)
    logging.info(f'Percentiles saved to {args.output_csv}')

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        logging.error(f'An error occurred: {e}')

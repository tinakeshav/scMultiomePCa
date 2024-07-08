from plot_colour_and_composition import *
import pandas as pd
import matplotlib.pyplot as plt

def load_colours(file_path):
    colours_annot = pd.read_csv(file_path)
    return dict(zip(colours_annot['annotation'], colours_annot['colour']))


def process_cohort(mt_filename, name, date, cols_to_plot, titles, annotation_to_colour):
    print(mt_filename)
    mt = pd.read_csv(mt_filename, index_col=0)
    print(mt.head())
    for col, title in zip(cols_to_plot, titles):
        count_df = get_count_df(mt, col, annotation_to_colour)
        fname = f'/results/figures/{name}_{col}_composition_bar_plots_with_counts_{date}.pdf'
        plot_composition(count_df, fname, title)
        plt.close()


def process_intrapatient(mt_filenames, pats, cols_to_plot, titles, date, annotation_to_colour):
    for filename, pat in zip(mt_filenames, pats):
    print(filename)
    mt = pd.read_csv(filename, index_col=0)
    print(mt.head())
    for col, title in zip(cols_to_plot, titles):
        count_df = get_count_df(mt, col, annotation_to_colour)
        fname = f'/results/figures/{pat}_{col}_composition_bar_plots_with_counts_{date}.pdf'
        plot_composition(count_df, fname, title)
        plt.close()


if __name__ == '__main__':
    annotation_to_colour = load_colours('data/colours_23032023_with_samples_averaged_vals_03032023_CNVrep.csv')

    # Process full inter-patient cohort
    if run_wnn_proc:
    cohort_mt_filename = 'data/mu_concat_processed_wnn_23062022_prospect_annotations_032022-scanpy_122023_03122023.csv'
    cohort_name = 'mu_concat_processed_wnn_23062022_prospect_annotations_03122023'
    cohort_date = 'tightlayout_11032024'
    cohort_cols_to_plot = ['Patient', 'TuorLN', 'broad_clusters_122023']
    cohort_titles = ['Patient', 'Site', 'Cell type annotation']

    process_cohort(cohort_mt_filename, cohort_name, cohort_date, cohort_cols_to_plot, cohort_titles, annotation_to_colour)

    # Process intra-patient data
    intrapat_date = 'tightlayout_25022024'
    intrapat_tu_filenames = [
    'data/merged_multiome_allPat1_29092022processed_clustered_chromvar_29092022_TuONLY_01112023TK_03112023_cluster_annot_short_09122023.csv',
    'data/merged_multiome_allPat2_29092022processed_clustered_chromvar_29092022_TuONLY_01112023TK_03112023_cluster_annot_short_09122023.csv',
    'data/merged_multiome_allPat4_29092022processed_clustered_chromvar_29092022_TuONLY_01112023TK_03112023_cluster_annot_short_09122023.csv',
    'data/merged_multiome_allPat5_29092022processed_clustered_chromvar_29092022_TuONLY_01112023TK_03112023_cluster_annot_short_09122023.csv'
    ]

    intrapat_tu_names = [
    'merged_multiome_allPat1_29092022_TuONLY',
    'merged_multiome_allPat2_29092022_TuONLY',
    'merged_multiome_allPat4_29092022_TuONLY',
    'merged_multiome_allPat5_29092022_TuONLY'
    ]

    intrapat_cols_to_plot = ['name', 'TK_03112023', 'TuorLN']
    intrapat_titles = ['Sample biopsy', 'Cell type annotation', 'Site']

    process_intrapatient(intrapat_tu_filenames, 
                         intrapat_tu_names, 
                         intrapat_cols_to_plot, 
                         intrapat_titles, 
                         intrapat_date, 
                         annotation_to_colour)




    intrapat_filenames = [
    'data/merged_multiome_allPat1_29092022processed_clustered_chromvar_29092022cluster_general_confirmDotPlotsMarkers_15052023_cluster_annot_short_09122023_replaced_colname_13022024.csv',
    'data/merged_multiome_allPat2_29092022processed_clustered_chromvar_29092022cluster_general_confirmDotPlotsMarkers_15052023_cluster_annot_short_09122023_replaced_colname_13022024.csv',
    'data/merged_multiome_allPat3_29092022processed_clustered_chromvar_29092022cluster_general_confirmDotPlotsMarkers_15052023_cluster_annot_short_09122023_replaced_colname_13022024.csv',
    'data/merged_multiome_allPat4_29092022processed_clustered_chromvar_29092022cluster_general_confirmDotPlotsMarkers_15052023_cluster_annot_short_09122023_replaced_colname_13022024.csv',
    'data/merged_multiome_allPat5_29092022processed_clustered_chromvar_29092022cluster_general_confirmDotPlotsMarkers_15052023_cluster_annot_short_09122023_replaced_colname_13022024.csv'
    ]

    intrapat_names = [
    'merged_multiome_allPat1_29092022', 
    'merged_multiome_allPat2_29092022', 
    'merged_multiome_allPat3_29092022', 
    'merged_multiome_allPat4_29092022', 
    'merged_multiome_allPat5_29092022', 
    ]

    intrapat_cols_to_plot = ['name', 'cluster_general_confirmDotPlotsMarkers_15052023', 'TuorLN']
    intrapat_titles = ['Sample biopsy', 'Cell type annotation', 'Site']

    process_intrapatient(intrapat_filenames, 
                         intrapat_names, 
                         intrapat_cols_to_plot, 
                         intrapat_titles, 
                         intrapat_date, 
                         annotation_to_colour)


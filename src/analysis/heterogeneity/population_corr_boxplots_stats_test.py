import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from collections import defaultdict
from numpy import where, nan
import numpy as np
from scipy.stats import mannwhitneyu


def process_corr_file(corr_filename, annotation_to_cancer_dict, annotation_to_colour_dict):
    """
    <corr_filename> : str of filename with correlation,
                      assumes annotations are as follows e.g. 8_Endothelial
    <annotation_to_cancer_dict> : dictionary of annotations to cancer or TME
    <annotation_to_colour_dict> : dictionary of annotations to colour
    """
    mt = pd.read_csv(corr_filename, index_col = 0)
    mt_melt = mt.reset_index().melt(id_vars='index').query("index != variable")

    # get annots
    mt_melt['annot1'] = mt_melt['index'].str.split('_', 1).str[1]
    mt_melt['annot2'] = mt_melt['variable'].str.split('_', 1).str[1]

    # map annots to get cancer versus TME
    mt_melt['type1'] = mt_melt['annot1'].map(annotation_to_cancer_dict)
    mt_melt['type2'] = mt_melt['annot2'].map(annotation_to_cancer_dict)

    mt_melt['colour1'] = mt_melt['annot1'].map(annotation_to_colour_dict)
    mt_melt['colour2'] = mt_melt['annot2'].map(annotation_to_colour_dict)

    # are corrs are within ME or not
    mt_melt['TMEtoME'] = mt_melt['type1'] != mt_melt['type2']

    # dataframe that contains only the unique comparisons
    mt_melt_val_dropped = mt_melt.drop_duplicates(subset=['value'], keep='first')

    return mt_melt, mt_melt_val_dropped

def get_proportions_from_metadata(metadata, colname1 = 'cluster_annot_short', colname2 = 'name', prop_threshold = 0.9):
    """
    returns a dataframe with proportions of each value in colname2 for each value in colname1
    <metadata> : pandas dataframe
    <colname1> : str, column name to group by
    <colname2> : str, column name to group by
    <prop_threshold> : float, threshold for proportion
    """
    prop_df = metadata.groupby([colname1, colname2]).size().unstack(fill_value=0)
    prop_df = prop_df.div(prop_df.sum(axis=1), axis=0)
    return prop_df

def get_counts_from_metadata(metadata, colname1 = 'cluster_annot_short'):
    """
    returns a dataframe with counts of each value in colname 
    <metadata> : pandas dataframe
    <colname1> : str, column name to group by
    """
    counts_df = metadata.groupby([colname1]).size().unstack(fill_value=0)
    return counts_df

def process_proportions(name_prop, prop_threshold = 0.9):
    """
    processes a proportions DataFrame to identify rows where values exceed
    a threshold, classifying as 'Tu' or 'LN' based on conditions.

    <name_prop> : str, path to the CSV file containing proportions data OR a pandas dataframe of metadata
                 which can be used to obtain proportions for
    <colname2> : str, column name to group by
    <prop_threshold> : float, threshold for identifying significant proportions.

    returns dictionary mapping column names to lists of indexes meeting conditions
    or 'NA' if no indexes meet the condition.
    """
    if isinstance(name_prop, str):
        metadata = pd.read_csv(name_prop, index_col = 0, header = 0)
        prop_df = get_proportions_from_metadata(prop_df, colname2 = colname2)
    else:
        prop_df = name_prop
    prop_t = prop_df.T

    # Classify index based on presence of 'Tu'
    prop_t['TuorLN'] = prop_t.index.map(lambda idx: 'Tu' if 'Tu' in idx else 'LN')

    # Initialize dictionary to hold results for every cluster
    annots_memberships = defaultdict(list)

    for column in prop_t.columns[:-1]:  # Exclude the 'TuorLN' classification column
        # Sum values by 'TuorLN' classification
        sums_by_group = prop_t.groupby('TuorLN')[column].sum()
        tuorln_classification = 'LN' if sums_by_group.get('LN', 0) > prop_threshold else 'Tu'

        # Find indexes where the condition is met, for each column
        condition_met_indexes = prop_t.index[prop_t[column] > prop_threshold].tolist()

        # Record the index or 'NA', and classification
        annots_memberships[column].append(condition_met_indexes[0] if condition_met_indexes else 'NA')
        annots_memberships[column].append(tuorln_classification)

    return annots_memberships

def add_vals_in_category(df, colname, expected_categories, return_included = True, val = nan):
    """
    defaults to NaN values in df[colname] if any expected_categories do not exist
    <df> : pandas dataframe
    <colname> : str, should exist in pandas dataframe
    <expected_categories> : list of expected categories in df
    <return_included> : returns included categories as list
    """
    included_categories = []
    new_row = {}
    for expected_category in expected_categories:
        if expected_category not in df[colname].unique():
            new_row = {col : val if col!=colname else expected_category for col in df.columns}
            # such that we can colour NA sample as white
            new_row['sample'] = 'NA'
            if val is not nan:
                new_row['value'] = val
            new_row = pd.DataFrame(new_row, index = [10])
            df = pd.concat([df, new_row], ignore_index = True)
        else:
            included_categories.append(expected_category)
    if return_included:
        return df, included_categories
    return df


def get_boxplot_colours(included_categories, order, annotation_to_colours):
    colours_to_plot = {}
    for i in range(len(order)):
        if order[i] in included_categories:
            colours_to_plot[order[i]] = annotation_to_colours[order[i]]
        else:
            colours_to_plot[order[i]] = 'white'
    return colours_to_plot


def set_ax_boxplot_colours(ax, colours):
    """
    ax : ax object
    colours : list of colours
    """
    for i in range(len(colours)):
        try:
            ax.patches[i].set_edgecolor(colours[i])
            ax.patches[i].set_linewidth(2)
            for j in range(6*i, 6*(i+1)):
                ax.lines[j].set_color(colours[i])
                ax.lines[j].set_linewidth(1)  
        except IndexError:
            print('Skipped setting colours due to lack of actual drawn patches')


def set_ax_fig_specs(ax, ylab = 'Pearson Correlation'):
    ax.set_xlabel('')
    ax.set_xticklabels(ax.get_xticklabels(),
                       rotation=90)
    ax.set_ylabel(ylab)
    ax.get_legend().remove()
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)


def plot_norm_entropy_boxplots(entropy_file_csv, date, annotation_to_colour):
    """
    plot boxplots for normalized entropies
    """
    ent_df = pd.read_csv(entropy_file_csv)
    ent_df['proliferating'] = ent_df.G2_M_markers > 0.2
    print('getting cancer values')
    # Mann-Whitney U test between Cancer and TME groups
    cancer_values = ent_df.query("is_cancer == 'Cancer'")['normalized_entropy']
    tme_values = ent_df.query("is_cancer == 'TME'")['normalized_entropy']
    stat, p_value = mannwhitneyu(cancer_values, tme_values, alternative='two-sided')
    print(f'Mann-Whitney U test p-value between Cancer and TME is {p_value}')

    # get a slightly smaller y lim for jitter?
    y_lim = (0,1)

    plt.figure(figsize=(2.5, 6))
    ax = sns.boxplot(data = ent_df
                    , x="is_cancer"
                    , y="normalized_entropy"
                    , palette = annotation_to_colour
                    , boxprops={'alpha': 0.2}
                    , linewidth=.85
                    , zorder = 1
                    , width = 0.5
                    , order = ['Cancer', 'TME']
                    )
    sns.stripplot(data = ent_df
                    , x="is_cancer"
                    , y="normalized_entropy"
                    , hue = 'proliferating'
                    #, style = 'proliferating'
                    #, markers = ['o', '*']
                    #, palette = annotation_to_colour
                    , ax = ax
                    , zorder = 2
                    , size = 5
                    , order = ['Cancer', 'TME']
                    )
    colours = [annotation_to_colour[lab.get_text()] for lab in ax.get_xticklabels()]
    set_ax_boxplot_colours(ax, colours)
    set_ax_fig_specs(ax, ylab='Normalized entropy')
    plt.tight_layout()
    plt.ylim(y_lim)


    y_max = ent_df['normalized_entropy'].max() + 0.03  
    if p_value < 0.001:
        ax.text(0.5, y_max, '***', ha='center', va='bottom', fontsize=12)
    elif p_value < 0.01:
        ax.text(0.5, y_max, '**', ha='center', va='bottom', fontsize=12)
    elif p_value < 0.05:
        ax.text(0.5, y_max, '*', ha='center', va='bottom', fontsize=12)

    # draw line connecting tested wo groups
    x1, x2 = 0, 1  # x positions for 'Cancer' and 'TME'
    # Vertical lines
    ax.plot([x1, x1, x2, x2],
            [y_max - 0.01, y_max + 0.01, y_max + 0.01, y_max - 0.01], color='black')
    # Horizontal lines
    #ax.plot([x1, x2], [y_max, y_max], color = 'black')


    plt.savefig(f'figures/allPatsTuONLY_G2M_EntropyNormMat_boxplot_pvals_{date}.pdf', dpi=600)

def main_plot_entropy(entropy_file_csv = 'data/allPatsTuONLY_EntropyNormMat_15022025.csv', date = '15022025'):
    plot_norm_entropy_boxplots(entropy_file_csv, '15022025')

def main_plot_intra_patient_corrs(corr_filenames, pats, annotation_to_colour, annotation_to_cancer, date = '02042024'):
    # correlations obtained from workflows/Snakefile_intrapat 
    corr_filenames = [
    'data/merged_multiome_allPat1_29092022processed_clustered_chromvar_29092022_TuONLY_01112023TK_03112023_cluster_annot_short_09122023__intrapat_peaks_DESeq2vstCorrMat_NonPromoterPeaks16012024.csv',
    'data/merged_multiome_allPat2_29092022processed_clustered_chromvar_29092022_TuONLY_01112023TK_03112023_cluster_annot_short_09122023__intrapat_peaks_DESeq2vstCorrMat_NonPromoterPeaks16012024.csv',
    'data/merged_multiome_allPat3_29092022processed_clustered_chromvar_29092022cluster_general_confirmDotPlotsMarkers_15052023_cluster_annot_short_09122023_replaced_colname_13022024__intrapat_peaks_DESeq2vstCorrMat_NonPromoterPeaks16012024.csv',
    'data/merged_multiome_allPat4_29092022processed_clustered_chromvar_29092022_TuONLY_01112023TK_03112023_cluster_annot_short_09122023__intrapat_peaks_DESeq2vstCorrMat_NonPromoterPeaks16012024.csv',
    'data/merged_multiome_allPat5_29092022processed_clustered_chromvar_29092022_TuONLY_01112023TK_03112023_cluster_annot_short_09122023__intrapat_peaks_DESeq2vstCorrMat_NonPromoterPeaks16012024.csv',
    ]
    pats = [
    'merged_multiome_allPat1_29092022_TuONLY',
    'merged_multiome_allPat2_29092022_TuONLY',
    'merged_multiome_allPat3_29092022',
    'merged_multiome_allPat4_29092022_TuONLY',
    'merged_multiome_allPat5_29092022_TuONLY',
    ]
    all_within_type_dfs = []
    all_me_to_tme_dfs = []
    all_tme_to_tme_dfs = []
    figsize = (2.5, 6)

    # read them one by one
    for pat, fname in dict(zip(pats, corr_filenames)).items():

        mt_melt, mt_melt_val_dropped = process_corr_file(fname,
                                       annotation_to_cancer,
                                       annotation_to_colour)

        y_lim_min = 0 if mt_melt.value.min() > 0 else mt_melt.value.min() - 0.05
        y_lim = (y_lim_min, mt_melt_val_dropped.value.max() + 0.05)

        ###########################################################################
        # Plot to showcase ME similarity wrt TME; TME used as reference points
        ###########################################################################
        df_within_type = mt_melt_val_dropped.query("TMEtoME == False")
        all_within_type_dfs.append(df_within_type)

        # Mann-Whitney U test between Cancer and TME groups
        cancer_values = df_within_type.query("type1 == 'Cancer'")['value']
        tme_values = df_within_type.query("type1 == 'TME'")['value']
        stat, p_value = mannwhitneyu(cancer_values, tme_values, alternative='two-sided')

        # print the p-value in log file
        print(pat, fname)
        print(f'Mann-Whitney U test p-value between Cancer and TME in {pat}: {p_value}')
        print(f'Pass? {p_value < 0.05}')
        if generate_pval_plot is True:
            plt.figure(figsize=(2.5, 6))
            ax = sns.boxplot(data = df_within_type
                            , x="type1"
                            , y="value"
                            , palette = annotation_to_colour
                            , boxprops={'alpha': 0.2}
                            , linewidth=.85
                            , zorder = 1
                            , width = 0.5
                            , order = ['Cancer', 'TME']
                            )
            sns.stripplot(data = df_within_type
                            , x="type1"
                            , y="value"
                            , hue = 'annot2'
                            , palette = annotation_to_colour
                            , ax = ax
                            , zorder = 2
                            , size = 5
                            , order = ['Cancer', 'TME']
                            )
            colours = [annotation_to_colour[lab.get_text()] for lab in ax.get_xticklabels()]
            set_ax_boxplot_colours(ax, colours)
            set_ax_fig_specs(ax)
            plt.tight_layout()
            plt.ylim(y_lim)


            y_max = df_within_type['value'].max() + 0.03  # Adjust y position for stars
            if p_value < 0.001:
                ax.text(0.5, y_max, '***', ha='center', va='bottom', fontsize=12)
            elif p_value < 0.01:
                ax.text(0.5, y_max, '**', ha='center', va='bottom', fontsize=12)
            elif p_value < 0.05:
                ax.text(0.5, y_max, '*', ha='center', va='bottom', fontsize=12)

            x1, x2 = 0, 1  # x positions for 'Cancer' and 'TME'
            ax.plot([x1, x1, x2, x2],
                    [y_max - 0.01, y_max + 0.01, y_max + 0.01, y_max - 0.01], color='black')
            #ax.plot([x1, x2], [y_max, y_max], color = 'black')
            plt.savefig(f'figures/{pat}_non_prom_corrs_pvals_{date}.pdf', dpi=600)


def main_plot_intra_patient_tuln_corrs(annotation_to_colour, annotation_to_cancer, date = '02042024'):
    # correlations obtained from workflows/Snakefile_intrapat 
    corr_filenames = [
    'data/merged_multiome_allPat1_29092022processed_clustered_chromvar_29092022_TuONLY_01112023TK_03112023_cluster_annot_short_09122023__intrapat_peaks_DESeq2vstCorrMat_NonPromoterPeaks16012024.csv',
    'data/merged_multiome_allPat2_29092022processed_clustered_chromvar_29092022_TuONLY_01112023TK_03112023_cluster_annot_short_09122023__intrapat_peaks_DESeq2vstCorrMat_NonPromoterPeaks16012024.csv',
    'data/merged_multiome_allPat3_29092022processed_clustered_chromvar_29092022cluster_general_confirmDotPlotsMarkers_15052023_cluster_annot_short_09122023_replaced_colname_13022024__intrapat_peaks_DESeq2vstCorrMat_NonPromoterPeaks16012024.csv',
    'data/merged_multiome_allPat4_29092022processed_clustered_chromvar_29092022_TuONLY_01112023TK_03112023_cluster_annot_short_09122023__intrapat_peaks_DESeq2vstCorrMat_NonPromoterPeaks16012024.csv',
    'data/merged_multiome_allPat5_29092022processed_clustered_chromvar_29092022_TuONLY_01112023TK_03112023_cluster_annot_short_09122023__intrapat_peaks_DESeq2vstCorrMat_NonPromoterPeaks16012024.csv',
    ]
    pats = [
    'merged_multiome_allPat1_29092022_TuONLY',
    'merged_multiome_allPat2_29092022_TuONLY',
    'merged_multiome_allPat3_29092022',
    'merged_multiome_allPat4_29092022_TuONLY',
    'merged_multiome_allPat5_29092022_TuONLY',
    ]
    # required for obtaining counts of clusters and appropriate proportions
    metadata_filenames = [
        'data/merged_multiome_allPat1_29092022processed_clustered_chromvar_29092022MarkersCellCycle_12112024.csv',
        'data/merged_multiome_allPat2_29092022processed_clustered_chromvar_29092022MarkersCellCycle_12112024.csv',
        'data/merged_multiome_allPat4_29092022processed_clustered_chromvar_29092022MarkersCellCycle_12112024.csv',
        'data/merged_multiome_allPat5_29092022processed_clustered_chromvar_29092022MarkersCellCycle_12112024.csv'
    ]

    count_threshold = 200 # minimum number of cells in a cluster to be considered
    figsize = (2.5,5.5)

    for pat, fname, metadata_fname in zip(pats, corr_filenames, metadata_filenames):

        print(pat, fname, metadata_fname)
        # read correlation plots
        mt_melt, mt_melt_val_dropped = process_corr_file(fname,
                                                        annotation_to_cancer,
                                                        annotation_to_colour)

        annots_memberships = process_proportions(metadata_fname,
                                                 prop_threshold = 0.9)

        # set limits for correlation plots
        #y_lim_min = 0 if mt_melt.value.min() > 0 else mt_melt.value.min() - 0.05
        y_lim = (0, 1)

        # retrieve ME to ME comparisons (all annot1s or annot2s will be ME anyway)
        df_me_to_me = mt_melt_val_dropped.query("TMEtoME == False and annot1 == 'Malignant Epithelial'")

        # map each annot1 or annot2 to get 1. Site membership 2. Biopsy membership
        # if it's same site, pass Intra-site {site}, if any are NA, pass NA, if not, pass Inter-Site
        df_me_to_me['idx_site'] = [annots_memberships[v][1] for v in df_me_to_me['index'].values]
        df_me_to_me['var_site'] = [annots_memberships[v][1] for v in df_me_to_me['variable'].values]


        df_me_to_me['site'] = where(df_me_to_me['idx_site'].isna() | df_me_to_me['var_site'].isna(), np.nan,
                                    where(df_me_to_me['idx_site'] == df_me_to_me['var_site'],
                                    df_me_to_me['idx_site'], 'Inter-site')
                                    )

        # Do the same for sample: if same, pass intra-sample, if NA then NA else pass Inter-sample
        df_me_to_me['idx_sample'] = [annots_memberships[v][0] for v in df_me_to_me['index'].values]
        df_me_to_me['var_sample'] = [annots_memberships[v][0] for v in df_me_to_me['variable'].values]
        df_me_to_me['sample'] = where(df_me_to_me['idx_sample'].isna() | df_me_to_me['var_sample'].isna(), np.nan,
                                    where(df_me_to_me['idx_sample'] == df_me_to_me['var_sample'],
                                    df_me_to_me['idx_sample'], 'Inter-sample')
                                    )

    
        # order annotations based on their membership
        order = ['LN', 'Inter-site', 'Tu']
        # mann-whitney u test between site comparisons
        grouped_data = {group: df_me_to_me[df_me_to_me['site'] == group]['value'] for group in order}
        p_values, p_values_order = get_pvals(grouped_data, order)

        #import ipdb; ipdb.set_trace()
        print(pat, fname, prop, count)
        print("Mann-Whitney U test p-values:")
        for comparison, p_value in p_values.items():
            print(f"{comparison}: p-value = {p_value}")


        # filter based on counts so we get rid of low-count/quality clusters
        # note: super low count clusters will generally have very low fractions of reads mapped to them
        # and can skew the analysis for no reason. 
        # thus its best to remove them for a more appropriate comparison
        if isinstance(count, str):
            count_df = pd.read_csv(count)
        else:
            count_df = get_counts_from_metadata(metadata)

        exclude_clusters = count_df.query(f"count < {count_threshold}")["cluster_annot_short"]
        print(f"excluding clusters {exclude_clusters}")
        print(exclude_clusters)
        exclude_mask = df_me_to_me.isin(exclude_clusters.values).any(axis = 1)

        plt.figure(figsize=figsize)

        # If expected values do not exist then add val
        # Had to add this chunk as there was a single Pat5 LN cluster with only 28 cells.
        df_excluded, included_categories = add_vals_in_category(df_me_to_me[~exclude_mask],
                                                                'site',
                                                                order,
                                                                val = 0.5)
        grouped_data = {group: df_excluded[df_excluded['site'] == group]['value'] for group in order}
        p_values, p_values_order = get_pvals(grouped_data, order)

        print(f'figures/{pat}_non_prom_corrs_intersite_hueSite_filter{count_threshold}_{date}.pdf')
        print(pat, fname, prop, count)
        print("Mann-Whitney U test p-values:")
        for comparison, p_value in p_values.items():
            print(f"{comparison}: p-value = {p_value}")

        print('included categories')
        print(included_categories)
        if len(included_categories) < len(order):
            # annotation_to_colour dict will be copied and the NA val will be white
            colours_dict = get_boxplot_colours(included_categories, order, annotation_to_colour)
        else:
            colours_dict = annotation_to_colour


        #print(colours_dict)
        ax = sns.boxplot(data = df_excluded
                        , x="site"
                        , y="value"
                        , palette = colours_dict
                        , boxprops={'alpha': 0.2}
                        , linewidth=.85
                        , zorder = 1
                        , width = 0.5
                        , order = order
                        )
        sns.stripplot(data = df_excluded
                        , x="site"
                        , y="value"
                        , hue = 'site'
                        , palette = colours_dict
                        , ax = ax
                        , zorder = 2
                        , size = 4
                        , order = order
                        )
        colours = [colours_dict[lab.get_text()] for lab in ax.get_xticklabels()]
        set_ax_boxplot_colours(ax, colours)
        set_ax_fig_specs(ax)
        plt.tight_layout()
        plt.ylim(y_lim)
        add_stats_to_plots(ax, p_values_order, order, y_lim)
        plt.savefig(f'figures/{pat}_non_prom_corrs_intersite_hueSite_filter{count_threshold}_{date}.pdf', dpi = 600)
        plt.close()


        plt.figure(figsize=small_figsize)
        ax = sns.boxplot(data = df_excluded
                        , x="site"
                        , y="value"
                        , palette = colours_dict
                        , boxprops={'alpha': 0.2}
                        , linewidth=.85
                        , zorder = 1
                        , width = 0.5
                        , order = order
                        )
        sns.stripplot(data = df_excluded
                        , x="site"
                        , y="value"
                        , hue = 'site'
                        , palette = colours_dict
                        , ax = ax
                        , zorder = 2
                        , size = 4
                        , order = order
                        )
        colours = [colours_dict[lab.get_text()] for lab in ax.get_xticklabels()]
        set_ax_boxplot_colours(ax, colours)
        set_ax_fig_specs(ax)
        plt.tight_layout()
        small_ylim = (0.4, 1.0)
        plt.ylim(small_ylim)
        add_stats_to_plots(ax, p_values_order, order, small_ylim)
        plt.savefig(f'figures/{pat}_non_prom_corrs_intersite_hueSite_ylim04_filter{count_threshold}_{date}.pdf', dpi = 600)



if __name__ == '__main__':
    # Prepare a mapping from annotations to cancer types and colours
    colours_annot = pd.read_csv('data/colours_23032023_with_samples_averaged_vals_03032023_CNVrep.csv')
    annotation_to_cancer = dict(zip(colours_annot['annotation'], colours_annot['cancer']))
    annotation_to_colour = dict(zip(colours_annot['annotation'], colours_annot['colour']))
    annotation_to_colour['NA'] = annotation_to_colour[np.nan]
    date = '02042024'
    plot_norm_entropy_boxplots(entropy_file_csv = 'data/allPatsTuONLY_EntropyNormMat_15022025.csv', date = date, annotation_to_colour = annotation_to_colour)
    main_plot_intra_patient_tuln_corrs(date = date, annotation_to_colour = annotation_to_colour, annotation_to_cancer = annotation_to_cancer)
    main_plot_intra_patient_corrs(date = date, annotation_to_colour = annotation_to_colour, annotation_to_cancer = annotation_to_cancer)

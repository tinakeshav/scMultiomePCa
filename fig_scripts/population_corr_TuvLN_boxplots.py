import numpy as np
import matplotlib.pyplot as plt
from population_corr_boxplots_utils import *


date = '05042024'
plt.rcParams.update({'font.size': 8})
figsize = (2.5,5.5)

corr_filenames = [
    '/data/merged_multiome_allPat1_29092022processed_clustered_chromvar_29092022cluster_general_confirmDotPlotsMarkers_15052023_cluster_annot_short_09122023_replaced_colname_13022024__intrapat_peaks_DESeq2vstCorrMat_NonPromoterPeaks16012024.csv',
    '/data/merged_multiome_allPat2_29092022processed_clustered_chromvar_29092022cluster_general_confirmDotPlotsMarkers_15052023_cluster_annot_short_09122023_replaced_colname_13022024__intrapat_peaks_DESeq2vstCorrMat_NonPromoterPeaks16012024.csv',
    '/data/merged_multiome_allPat4_29092022processed_clustered_chromvar_29092022cluster_general_confirmDotPlotsMarkers_15052023_cluster_annot_short_09122023_replaced_colname_13022024__intrapat_peaks_DESeq2vstCorrMat_NonPromoterPeaks16012024.csv',
    '/data/merged_multiome_allPat5_29092022processed_clustered_chromvar_29092022cluster_general_confirmDotPlotsMarkers_15052023_cluster_annot_short_09122023_replaced_colname_13022024__intrapat_peaks_DESeq2vstCorrMat_NonPromoterPeaks16012024.csv'
]


proportion_filenames = [
    '/data/merged_multiome_allPat1_29092022processed_clustered_chromvar_29092022cluster_general_confirmDotPlotsMarkers_15052023_cluster_annot_short_09122023_replaced_colname_13022024_proportions_cluster_annot_short_name_13022024.csv',
    '/data/merged_multiome_allPat2_29092022processed_clustered_chromvar_29092022cluster_general_confirmDotPlotsMarkers_15052023_cluster_annot_short_09122023_replaced_colname_13022024_proportions_cluster_annot_short_name_13022024.csv',
    '/data/merged_multiome_allPat4_29092022processed_clustered_chromvar_29092022cluster_general_confirmDotPlotsMarkers_15052023_cluster_annot_short_09122023_replaced_colname_13022024_proportions_cluster_annot_short_name_13022024.csv',
    '/data/merged_multiome_allPat5_29092022processed_clustered_chromvar_29092022cluster_general_confirmDotPlotsMarkers_15052023_cluster_annot_short_09122023_replaced_colname_13022024_proportions_cluster_annot_short_name_13022024.csv'
]


pats = [
    'merged_multiome_allPat1_29092022', 
    'merged_multiome_allPat2_29092022', 
    'merged_multiome_allPat4_29092022', 
    'merged_multiome_allPat5_29092022'
]


# Filter rows / columns based on counts, as one wishes
counts_filenames = [
        '/data/merged_multiome_allPat1_29092022processed_clustered_chromvar_29092022cluster_general_confirmDotPlotsMarkers_15052023_cluster_annot_short_09122023_replaced_colname_13022024_count_cluster_annot_short_13022024.csv',
        '/data/merged_multiome_allPat2_29092022processed_clustered_chromvar_29092022cluster_general_confirmDotPlotsMarkers_15052023_cluster_annot_short_09122023_replaced_colname_13022024_count_cluster_annot_short_13022024.csv',
        '/data/merged_multiome_allPat4_29092022processed_clustered_chromvar_29092022cluster_general_confirmDotPlotsMarkers_15052023_cluster_annot_short_09122023_replaced_colname_13022024_count_cluster_annot_short_13022024.csv',
        '/data/merged_multiome_allPat5_29092022processed_clustered_chromvar_29092022cluster_general_confirmDotPlotsMarkers_15052023_cluster_annot_short_09122023_replaced_colname_13022024_count_cluster_annot_short_13022024.csv'
        ]

count_threshold = 200



# Prepare a mapping from annotations to cancer types and colours
colours_annot = pd.read_csv('/data/colours_23032023_with_samples_averaged_vals_03032023_CNVrep.csv')
annotation_to_cancer = dict(zip(colours_annot['annotation'], colours_annot['cancer']))
annotation_to_colour = dict(zip(colours_annot['annotation'], colours_annot['colour']))
annotation_to_colour['NA'] = annotation_to_colour[np.nan]

# read them one by one
for pat, fname, prop, count in zip(pats, corr_filenames, 
                                   proportion_filenames, counts_filenames):

    print(pat, fname, prop)
    # read correlation plots
    mt_melt, mt_melt_val_dropped = process_corr_file(fname, 
                                   annotation_to_cancer, 
                                   annotation_to_colour)

    annots_memberships = process_proportions(prop, 
                         prop_threshold = 0.9)

    # set limits for correlation plots
    #y_lim_min = 0 if mt_melt.value.min() > 0 else mt_melt.value.min() - 0.05
    y_lim = (0, 1)

    ###########################################################################
    # New Plots
    ##########################################################################

    # Retrieve ME to ME comparisons (all annot1s or annot2s will be ME anyway)
    df_me_to_me = mt_melt_val_dropped.query("TMEtoME == False and annot1 == 'Malignant Epithelial'")

    # Map each annot1 or annot2 to get 1. Site membership 2. Biopsy membership
    # If same site, pass Intra-site {site}, if any are NA, pass NA, if not, pass Inter-Site
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

    # We plot
    plt.figure(figsize=figsize)

    order = ['LN', 'Inter-site', 'Tu']

    ax = sns.boxplot(data = df_me_to_me
                     , x="site"
                     , y="value" 
                     , palette = annotation_to_colour
                     , boxprops={'alpha': 0.2}
                     , linewidth=.85
                     , zorder = 1
                     , width = 0.5
                     , order = order
                     )
    sns.stripplot(data = df_me_to_me
                    , x="site"
                    , y="value"
                    , hue = 'sample'
                    , palette = annotation_to_colour
                    , ax = ax
                    , zorder = 2
                    , size = 4
                    , order = order
                    )
    colours = [annotation_to_colour[lab.get_text()] for lab in ax.get_xticklabels()]
    set_ax_boxplot_colours(ax, colours)
    set_ax_fig_specs(ax)
    plt.tight_layout()
    plt.ylim(y_lim)
    plt.savefig(f'/results/figures/{pat}_non_prom_corrs_intersite_hueSample_{date}.pdf', dpi = 600)
    plt.close(); plt.close(); plt.close(); plt.close()



    plt.figure(figsize=figsize)
    ax = sns.boxplot(data = df_me_to_me
                     , x="site"
                     , y="value" 
                     , palette = annotation_to_colour
                     , boxprops={'alpha': 0.2}
                     , linewidth=.85
                     , zorder = 1
                     , width = 0.5
                     , order = order
                     )
    sns.stripplot(data = df_me_to_me
                    , x="site"
                    , y="value"
                    , hue = 'site'
                    , palette = annotation_to_colour
                    , ax = ax
                    , zorder = 2
                    , size = 4
                    , order = order
                    )
    colours = [annotation_to_colour[lab.get_text()] for lab in ax.get_xticklabels()]
    set_ax_boxplot_colours(ax, colours)
    set_ax_fig_specs(ax)
    plt.tight_layout()
    plt.ylim(y_lim)
    plt.savefig(f'/results/figures/{pat}_non_prom_corrs_intersite_hueSite_{date}.pdf', dpi = 600)
    plt.close(); plt.close(); plt.close(); plt.close()




    # Same plots but filter on counts
    count_df = pd.read_csv(count)
    exclude_clusters = count_df.query(f"count < {count_threshold}")["cluster_annot_short"]
    print(f"excluding clusters {exclude_clusters}")
    print(exclude_clusters)
    exclude_mask = df_me_to_me.isin(exclude_clusters.values).any(axis = 1)

    print(df_me_to_me[~exclude_mask])
    print(df_me_to_me[df_me_to_me['site'] == 'LN'])
    print(df_me_to_me[~exclude_mask]['site'].unique())

    plt.figure(figsize=figsize)

    order = ['LN', 'Inter-site', 'Tu']

    # adding NA values in case we've removed certain categories in the process
    # This was a problem for Patient5 as the only "LN to LN" column had 28 cells

    # If expected values do not exist then add val
    df_excluded, included_categories = add_vals_in_category(df_me_to_me[~exclude_mask], 
                                       'site', 
                                       order, 
                                       val = 0.5)
    print('included categories')
    print(included_categories)
    if len(included_categories) < len(order):
        # annotation_to_colour dict will be copied and the NA val will be white
        colours_dict = get_boxplot_colours(included_categories, order, annotation_to_colour)
    else:
        colours_dict = annotation_to_colour

    
    print(colours_dict)
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
    plt.savefig(f'/results/figures/{pat}_non_prom_corrs_intersite_hueSite_filter{count_threshold}_{date}.pdf', dpi = 600)
    plt.close(); plt.close(); plt.close(); plt.close()


    plt.figure(figsize=figsize)
    ax = sns.boxplot(data = df_excluded
                     , x="site"
                     , y="value" 
                     , palette = annotation_to_colour
                     , boxprops={'alpha': 0.2}
                     , linewidth=.85
                     , zorder = 1
                     , width = 0.5
                     , order = order
                     )
    sns.stripplot(data = df_excluded
                    , x="site"
                    , y="value"
                    , hue = 'sample'
                    , palette = annotation_to_colour 
                    , ax = ax
                    , zorder = 2
                    , size = 4
                    , order = order
                    )
    # setting ax boxplot colours needs to take care of cases 
    # where patches are not drawn
    # due to NA values
    colours = [colours_dict[lab.get_text()] for lab in ax.get_xticklabels()]
    set_ax_boxplot_colours(ax, colours)
    set_ax_fig_specs(ax)
    plt.tight_layout()
    plt.ylim(y_lim)
    plt.savefig(f'/results/figures/{pat}_non_prom_corrs_intersite_hueSample_filter{count_threshold}_{date}.pdf', dpi = 600)
    plt.close(); plt.close(); plt.close(); plt.close()

    
    print(colours_dict)
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
    plt.ylim((0, 0.4))
    plt.savefig(f'/results/figures/{pat}_non_prom_corrs_intersite_hueSite_ylim04_filter{count_threshold}_{date}.pdf', dpi = 600)
    plt.close(); plt.close(); plt.close(); plt.close()

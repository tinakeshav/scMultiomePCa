import pandas as pd
import population_corr_boxplots_utils 
import seaborn as sns
import matplotlib.pyplot as plt
from numpy import where, nan

if __name__ == '__main__':


	plt.rcParams.update({'font.size': 10})

	date = '02042024'
	date = 'fsize22_02042024'

	corr_filenames = [
	'/data/merged_multiome_allPat1_29092022processed_clustered_chromvar_29092022_TuONLY_01112023TK_03112023_cluster_annot_short_09122023__intrapat_peaks_DESeq2vstCorrMat_NonPromoterPeaks16012024.csv',
	'/data/merged_multiome_allPat2_29092022processed_clustered_chromvar_29092022_TuONLY_01112023TK_03112023_cluster_annot_short_09122023__intrapat_peaks_DESeq2vstCorrMat_NonPromoterPeaks16012024.csv',
	'/data/merged_multiome_allPat3_29092022processed_clustered_chromvar_29092022cluster_general_confirmDotPlotsMarkers_15052023_cluster_annot_short_09122023_replaced_colname_13022024__intrapat_peaks_DESeq2vstCorrMat_NonPromoterPeaks16012024.csv',
	'/data/merged_multiome_allPat4_29092022processed_clustered_chromvar_29092022_TuONLY_01112023TK_03112023_cluster_annot_short_09122023__intrapat_peaks_DESeq2vstCorrMat_NonPromoterPeaks16012024.csv',
	'/data/merged_multiome_allPat5_29092022processed_clustered_chromvar_29092022_TuONLY_01112023TK_03112023_cluster_annot_short_09122023__intrapat_peaks_DESeq2vstCorrMat_NonPromoterPeaks16012024.csv',
	'/data/merged_multiome_allPat1_29092022processed_clustered_chromvar_29092022cluster_general_confirmDotPlotsMarkers_15052023_cluster_annot_short_09122023_replaced_colname_13022024__intrapat_peaks_DESeq2vstCorrMat_NonPromoterPeaks16012024.csv',
	'/data/merged_multiome_allPat2_29092022processed_clustered_chromvar_29092022cluster_general_confirmDotPlotsMarkers_15052023_cluster_annot_short_09122023_replaced_colname_13022024__intrapat_peaks_DESeq2vstCorrMat_NonPromoterPeaks16012024.csv',
	'/data/merged_multiome_allPat4_29092022processed_clustered_chromvar_29092022cluster_general_confirmDotPlotsMarkers_15052023_cluster_annot_short_09122023_replaced_colname_13022024__intrapat_peaks_DESeq2vstCorrMat_NonPromoterPeaks16012024.csv',
	'/data/merged_multiome_allPat5_29092022processed_clustered_chromvar_29092022cluster_general_confirmDotPlotsMarkers_15052023_cluster_annot_short_09122023_replaced_colname_13022024__intrapat_peaks_DESeq2vstCorrMat_NonPromoterPeaks16012024.csv'
	]

	pats = [
	'merged_multiome_allPat1_29092022_TuONLY',
	'merged_multiome_allPat2_29092022_TuONLY', 
	'merged_multiome_allPat3_29092022', 
	'merged_multiome_allPat4_29092022_TuONLY', 
	'merged_multiome_allPat5_29092022_TuONLY',
	'merged_multiome_allPat1_29092022', 
	'merged_multiome_allPat2_29092022', 
	'merged_multiome_allPat4_29092022', 
	'merged_multiome_allPat5_29092022'
	]


    all_within_type_dfs = []
    all_me_to_tme_dfs = []
    all_tme_to_tme_dfs = []

    # read them one by one
    for pat, fname in dict(zip(pats, corr_filenames)).items():
        # Prepare a mapping from annotations to cancer types and colours
        colours_annot = pd.read_csv('data/colours_23032023_with_samples_averaged_vals_03032023_CNVrep.csv')
        annotation_to_cancer = dict(zip(colours_annot['annotation'], colours_annot['cancer']))
        annotation_to_colour = dict(zip(colours_annot['annotation'], colours_annot['colour']))

        mt_melt, mt_melt_val_dropped = process_corr_file(fname, 
                                       annotation_to_cancer, 
                                       annotation_to_colour)

        y_lim_min = 0 if mt_melt.value.min() > 0 else mt_melt.value.min() - 0.05
        y_lim = (y_lim_min, mt_melt_val_dropped.value.max() + 0.05)

        ###########################################################################
        # Plots
        ###########################################################################
        # First showcase that MEs are more similar to each other than TMEs
        # i.e. within ME comparison and within TME comparison
        ###########################################################################
        df_within_type = mt_melt_val_dropped.query("TMEtoME == False")
        all_within_type_dfs.append(df_within_type)
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
        plt.savefig(f'figures/{pat}_non_prom_corrs_{date}.pdf', dpi = 600)
        plt.close(); plt.close(); plt.close(); plt.close()

        ################################################################## 
        # Plot
        # Second, which TMEs are most corr to MEs?
        ################################################################## 
        df_me_to_tme = mt_melt.query("TMEtoME == True and annot1 == 'Malignant Epithelial'")
        all_me_to_tme_dfs.append(df_me_to_tme)
        plt.figure(figsize=(2.5, 6))
        ax = sns.boxplot(data = df_me_to_tme
                         , x="type1"
                         , y="value" 
                         , palette = annotation_to_colour
                         , boxprops={'alpha': 0.2}
                         , linewidth=.85
                         , zorder = 1
                         , width = 0.5
                         )

        sns.stripplot(data = df_me_to_tme, 
                        x="type1", 
                        y="value", 
                        hue = 'annot2', 
                        palette = annotation_to_colour,
                        ax = ax, 
                        zorder = 2,
                        size = 5)

        colours = [annotation_to_colour[lab.get_text()] for lab in ax.get_xticklabels()]
        set_ax_boxplot_colours(ax, colours)
        set_ax_fig_specs(ax)
        plt.tight_layout()
        plt.ylim(y_lim)
        plt.savefig(f'figures/{pat}_non_prom_corrs_tme_to_me_{date}.pdf', dpi = 600)
        plt.close(); plt.close()

        ################################################################## 
        # Plot
        # Third, which TMEs are most corr to TMEs?
        ################################################################## 
        df_tme_to_tme = mt_melt.query("TMEtoME == False and type1 == 'TME'")
        all_tme_to_tme_dfs.append(df_tme_to_tme)

        plt.figure(figsize=(len(df_tme_to_tme['annot1'].unique())*1.5, 6))
        ax = sns.boxplot(data = df_tme_to_tme
                         , x="annot1"
                         , y="value" 
                         , palette = annotation_to_colour
                         , boxprops={'alpha': 0.2}
                         , linewidth=.85
                         , zorder = 1
                         , width = 0.5
                         )
        sns.stripplot(data = df_tme_to_tme, 
                        x="annot1", 
                        y="value", 
                        hue = 'annot2', 
                        palette = annotation_to_colour,
                        ax = ax, 
                        zorder = 2,
                        size = 5)

        colours = [annotation_to_colour[lab.get_text()] for lab in ax.get_xticklabels()]
        set_ax_boxplot_colours(ax, colours)
        set_ax_fig_specs(ax)
        plt.tight_layout()
        plt.savefig(f'figures/{pat}_non_prom_corrs_tme_to_tme_{date}.pdf', dpi = 600)
        plt.close(); plt.close()

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from collections import defaultdict 
from numpy import where, nan

def process_corr_file(corr_filename, annotation_to_cancer_dict, annotation_to_colour_dict):
    """
    <corr_filename> : str of filename with correlation, 
                      assumes annotations are as follows e.g. 8_Endothelial
    <annotation_to_cancer_dict> : dictionary of annotations to cancer or TME
    <annotation_to_colour_dict> : dictionary of annotations to colour
    """
    mt = pd.read_csv(corr_filename, index_col = 0)
    # 'melt' it
    mt_melt = mt.reset_index().melt(id_vars='index').query("index != variable")
    
    # Get their annotations
    mt_melt['annot1'] = mt_melt['index'].str.split('_', 1).str[1]
    mt_melt['annot2'] = mt_melt['variable'].str.split('_', 1).str[1]

    # Map them to get cancer versus TME
    mt_melt['type1'] = mt_melt['annot1'].map(annotation_to_cancer_dict)
    mt_melt['type2'] = mt_melt['annot2'].map(annotation_to_cancer_dict)

    mt_melt['colour1'] = mt_melt['annot1'].map(annotation_to_colour_dict)
    mt_melt['colour2'] = mt_melt['annot2'].map(annotation_to_colour_dict)

    # Determine if corrs are within ME or not
    mt_melt['TMEtoME'] = mt_melt['type1'] != mt_melt['type2']

    # A dataframe that contains only the unique comparisons
    mt_melt_val_dropped = mt_melt.drop_duplicates(subset=['value'], keep='first')

    return mt_melt, mt_melt_val_dropped


def process_proportions(name_prop_filename, prop_threshold = 0.9):
    """
    Processes a proportions DataFrame to identify rows where values exceed
    a threshold, classifying as 'Tu' or 'LN' based on conditions.

    Parameters:
    - name_prop: Path to the CSV file containing proportions data.
    - prop_threshold: Threshold for identifying significant proportions.

    Returns:
    A dictionary mapping column names to lists of indexes meeting conditions
    or 'NA' if no indexes meet the condition.
    """
    prop_df = pd.read_csv(name_prop_filename, index_col = 0, header = 0)
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
                ax.lines[j].set_linewidth(1)  # Change linewidth
        except IndexError as e:
            print('Skipped setting colours due to lack of actual drawn patches')


def set_ax_fig_specs(ax, ylab = 'Pearson Correlation'):
    ax.set_xlabel('')
    ax.set_xticklabels(ax.get_xticklabels(), 
                       rotation=90)
    ax.set_ylabel(ylab)
    ax.get_legend().remove()
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)





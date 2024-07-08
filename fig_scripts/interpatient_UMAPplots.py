import numpy as np
import colorsys
import muon as mu
import matplotlib.pyplot as plt
import pandas as pd

def gg_colour_hue(n):
    """
    function to emulate ggplot2 colours as closely as possible in Python
    s/o to chatGPT
    returns a list of hex colours
    """
    # Generate a sequence of hues
    hues = np.linspace(15, 375, num=n+1)[:-1]  # Exclude the last element to match R's behavior
    # Convert hues to [0, 1] for colorsys
    hues_norm = hues / 360.0
    # Convert HCL (approximately using HLS) to RGB
    colours = [colorsys.hls_to_rgb(h, 0.65, 1) for h in hues_norm]  # l=0.65 approximates l=65 in R, c is approximated with s=1
    # Convert RGB from [0, 1] to hex
    colours_hex = ['#%02x%02x%02x' % tuple(int(c * 255) for c in colour) for colour in colours]
    return colours_hex

def load_metadata(file_path):
    metadata = pd.read_csv(file_path, header=0, index_col=0)
    metadata['TuorLN'] = metadata['ATAC:name'].str[5:7]
    return metadata

def load_colours(file_path):
    colours_df = pd.read_csv(file_path)
    return dict(zip(colours_df.annotation, colours_df.colour))

def generate_population_dict(metadata, n_clusters):
    return {str(i): gg_colour_hue(n_clusters+1)[i] for i in range(n_clusters+1)}

def save_umap_plot(mu_object, colour_by, palette, 
                   save_path, legend_loc='on data'):
    plt.figure(figsize=(5,5))
    mu.pl.embedding(mu_object, basis="X_wnn_umap", 
                    color=colour_by, palette=palette, 
                    legend_loc=legend_loc, legend_fontsize='small')
    plt.axis('off')
    plt.title('')
    plt.savefig(save_path, dpi=600, bbox_inches='tight')

if __name__ == '__main__':
    # Define dates
    wnn_date = '23062022'
    # Save figures with following date
    experiment_date = '25032024'

    # Load metadata and colours
    metadata_file = f'data/mu_concat_processed_wnn_{wnn_date}_prospect_annotations_032022-scanpy_122023_03122023.csv'
    metadata = load_metadata(metadata_file)

    colours_file = 'data/colours_23032023_with_samples_averaged_vals_03032023.csv'
    colours_dict = load_colours(colours_file)

    # Generate population dictionary
    n_clusters = max(metadata['leiden_wnn'])
    population_dict = generate_population_dict(metadata, n_clusters)
    print(population_dict)

    # Load WNN object
    wnn_file = f'data/mu_concat_processed_wnn_{wnn_date}.h5mu'
    pro_mu = mu.read(wnn_file)

    # Save UMAP plots
    get_pltname = lambda c: f'figures/mu_concat_processed_wnn_{wnn_date}_UMAP_{c}_{experiment_date}.png'

    # Plot with 'leiden_wnn'
    save_umap_plot(pro_mu, 'leiden_wnn', population_dict, get_pltname('leiden_wnn'))

    # Plot with other colour_by columns
    colour_by_cols = ['TuorLN', 'Patient', 'Pat_TuorLN']
    for colour_by in colour_by_cols:
    save_umap_plot(pro_mu, colour_by, colours_dict, get_pltname(colour_by), legend_loc=None)

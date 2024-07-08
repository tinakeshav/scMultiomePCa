# util function to plot a circle next to a bar plus text on top
import matplotlib.pyplot as plt
import matplotlib
from numpy import log2
import pandas as pd
#matplotlib.rcParams['font.sans-serif'] = ['Arial']
#matplotlib.rcParams['font.family'] = "sans-serif"
#plt.rc('font', family='Arial')

# retrieve metadata and get counts and colour
def get_count_df(mt_df, col, colours_dict):
    """
    returns a count to annotation df; not a plotting util
    """
    count_df = pd.DataFrame(mt_df[col].value_counts()).reset_index()
    count_df = count_df.rename(columns = {'index' : 'annotation', col : 'count'})
    count_df['colour'] = count_df['annotation'].map(colours_dict)
    return count_df

def get_text_boundaries(ax, ax_text):
    """
    returns text boundaries
    """
    bbox_display_units = ax_text.get_window_extent()
    bbox_data_units = ax.transData.inverted().transform_bbox(bbox_display_units)
    return bbox_data_units.y0, bbox_data_units.y1


def get_scaling(top_boundaries, bottom_boundaries):
    """
    returns how much bar height needs to be 'scaled'
    """
    return top_boundaries[0] - bottom_boundaries[1]


def plot_composition(df, plt_savename, title = None):
    """
    df pd dataframe with columns group,count,colour
    bars are 'normalized' to max count
    """
    # first plot the bars individually
    df_sorted = df.sort_values('count')

    # variables
    height_multiplier = 1
    fig_height = height_multiplier*(df_sorted.shape[0])
    #bar_height = 0.2/df_sorted.shape[0] # y width wrt circle
    text_y_pos_add = 0.05
    text_fontsize = 28

    
    # adding one to the shape because title will take up 1/xth of the portion
    fig, ax = plt.subplots(figsize = (7, fig_height))

    # plot circles
    # get circle x position
    max_width = max(df['count']) 
    # multiply by 0.03 to move circle 3% left of the barplot
    # and the x limit is defined based on max count
    circles = []
    circle_x_position = -max_width * 0.03

    for index, colour in enumerate(df_sorted['colour']):
        circles.append(ax.text(circle_x_position,
                       index,
                       '●', 
                       fontsize=25, 
                       color=colour, 
                       va='center', 
                       ha='center')
                       )

    fig.canvas.draw()

    texts = []
    for circle, annot, count in zip(circles, df_sorted['annotation'], df_sorted['count']):
        texts.append(ax.text(0, 
                            get_text_boundaries(ax, circle)[1],
                            #bar.get_y() + bar.get_height()*4, 
                            f"{annot} ({count:,})",
                            va = 'bottom', 
                            ha = 'left', 
                            color = 'black', 
                            fontsize = text_fontsize)
                    )

    
    if df.shape[0] > 1:
        top_circle_boundaries = get_text_boundaries(ax, circles[-1])
        bottom_circle_boundaries = get_text_boundaries(ax, circles[-2])

        print(top_circle_boundaries)
        print(bottom_circle_boundaries)

        circle_height = top_circle_boundaries[0] - bottom_circle_boundaries[1]

    else:
        circle_height = 1 


    bars = ax.barh(df_sorted['annotation'], 
                   df_sorted['count'], 
                   color = df_sorted['colour'], 
                   height = circle_height*0.1)


    fig.canvas.draw()

    # remove labels
    ax.tick_params(axis='both', 
                   which='both', 
                   left=False, 
                   bottom=False, 
                   labelleft=False, 
                   labelbottom=False)
    
    # remove spines
    [ax.spines[s].set_visible(False) for s in ax.spines]

    # set x limit so circles are visible
    #ax.set_xlim(min(-max_width * 0.05, circle_x_position), max_width * 1.1)

    if title:
        # compute exactly where the first bar is taking up
        text_boundary = get_text_boundaries(ax, texts[-1])
        text_width = text_boundary[1] - text_boundary[0]

        ax_title = ax.text(0, 
                    text_boundary[1] + text_width,
                    title,
                    fontsize=text_fontsize,
                    va = 'bottom'
                    )

    # set y lim to the most bottom position
    # and to the most top position
    ax.set_ylim(get_text_boundaries(ax, circles[0])[0], get_text_boundaries(ax, ax_title)[1]+0.05)
    plt.tight_layout()
    plt.savefig(plt_savename, dpi = 600)

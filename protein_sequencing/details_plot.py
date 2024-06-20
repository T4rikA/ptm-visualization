from collections import defaultdict
import os
import plotly.graph_objects as go
import pandas as pd

from protein_sequencing import parameters, utils,sequence_plot


def get_present_regions_cleavage(cleavage_df: pd.DataFrame):
    cleavages = cleavage_df.iloc[0:1,2:].values[0].tolist()
    cleavage_ranges = []
    for range in cleavages:
        if '-' in str(range):
            start, end = map(int, range.split('-'))
            cleavage_ranges.append((start, end))
        else:
            start = end = int(range)
            cleavage_ranges.append((start, end))

    region_ranges = []
    region_start = 1
    for region_name, region_end, region_color in parameters.REGIONS:
        region_ranges.append((region_start, region_end))
        region_start = region_end + 1

    regions_present = [False] * len(region_ranges)
    region_index = 0
    for cleavage_range in cleavage_ranges:
        while cleavage_range[0] > region_ranges[region_index][1]:
            region_index += 1
        regions_present[region_index] = True
    return regions_present

def plot_line_with_label_horizontal(fig: go.Figure, x_0: int, x_1: int, y_0: int, y_1: int, y_2: int, y_3: int, y_label: int, label: str, color):
    fig.add_trace(go.Scatter(x=[x_0, x_0, x_1, x_1],
                            y=[y_0, y_1, y_2, y_3],
                            mode='lines',
                            line=dict(color="black", width=1), showlegend=False, hoverinfo='none'))
    fig.add_annotation(x=x_1, y=y_label,
                        text=label,
                        showarrow=False,
                        textangle=-90,
                        font=dict(
                            family=parameters.FONT,
                            size=parameters.SEQUENCE_PLOT_FONT_SIZE,
                            color=color,
                            ))
    return fig

def plot_range_with_label_horizontal(fig: go.Figure, x_0_start: int, x_0_end: int, x_1: int, y_0: int, y_1: int, y_2: int, y_3: int, y_label: int, label: str):
    fig.add_trace(go.Scatter(x=[x_0_start, x_0_start, x_1, x_1, x_1, x_0_end, x_0_end],
                            y=[y_0, y_1, y_2, y_3, y_2, y_1, y_0],
                            mode='lines',
                            fill='toself',
                            line=dict(color="black", width=1), showlegend=False, hoverinfo='none'))
    fig.add_annotation(x=x_1, y=y_label,
                        text=label,
                        showarrow=False,
                        textangle=-90,
                        font=dict(
                            family=parameters.FONT,
                            size=parameters.SEQUENCE_PLOT_FONT_SIZE,
                            color=parameters.CLEAVAGE_LABEL_COLOR,
                            ))
    return fig

def plot_cleavage_labels(fig: go.Figure, cleavage_df: pd.DataFrame, pixels_per_cleavage: int, label_plot_height: int, group: str):
    cleavages = cleavage_df.iloc[0:1,2:].values[0].tolist()
    longest_label = ''
    for cleavage in cleavages:
        if utils.get_label_length(str(cleavage)) > utils.get_label_length(longest_label):
            longest_label = str(cleavage)

    group_direction = 1 if group == 'A' else -1
    cleavages_visited = 0
    last_end = parameters.REGIONS[0][1]
    last_region = 0
    for cleavage in cleavages:
        if '-' in str(cleavage):
            start, end = map(int, cleavage.split('-'))
        else:
            start = end = int(cleavage)
        if start > last_end:
            last_region += 1
            last_end = parameters.REGIONS[last_region + 1][1]
            cleavages_visited += 1
        if parameters.FIGURE_ORIENTATION == 0:
            if start == end:
                label = str(start)
                x_0_line = start * utils.PIXELS_PER_PROTEIN + utils.SEQUENCE_OFFSET
                x_1_line = cleavages_visited * pixels_per_cleavage + utils.SEQUENCE_OFFSET
                y_0_line = utils.SEQUENCE_BOUNDARIES['y1'] if group == 'A' else utils.SEQUENCE_BOUNDARIES['y0']
                y_1_line = y_0_line + 10 * group_direction
                y_3_line = y_0_line + (label_plot_height - utils.get_label_length(label)) * group_direction
                y_2_line = y_3_line - (utils.get_label_length(longest_label)-utils.get_label_length(label)+10) * group_direction
                y_label = y_3_line + utils.get_label_length(label)//2 * group_direction
                
                plot_line_with_label_horizontal(fig,
                                 x_0_line, x_1_line,
                                 y_0_line, y_1_line, y_2_line, y_3_line,
                                 y_label,
                                 label, parameters.CLEAVAGE_LABEL_COLOR)
            else:
                label = f'{start}-{end}'
                x_0_start_line = start * utils.PIXELS_PER_PROTEIN + utils.SEQUENCE_OFFSET
                x_0_end_line = end * utils.PIXELS_PER_PROTEIN + utils.SEQUENCE_OFFSET
                x_1_line = cleavages_visited * pixels_per_cleavage + utils.SEQUENCE_OFFSET
                y_0_line = utils.SEQUENCE_BOUNDARIES['y1'] if group == 'A' else utils.SEQUENCE_BOUNDARIES['y0']
                y_1_line = y_0_line + 10 * group_direction
                y_3_line = y_0_line + (label_plot_height - utils.get_label_length(label)) * group_direction
                y_2_line = y_3_line - 10 * group_direction
                y_label = y_3_line + utils.get_label_length(label)//2 * group_direction

                plot_range_with_label_horizontal(fig,
                                    x_0_start_line, x_0_end_line, x_1_line,
                                    y_0_line, y_1_line, y_2_line, y_3_line,
                                    y_label,
                                    label)
        else:
            pass

        cleavages_visited += 1

def filter_relevant_modification_sights(ptm_file: str):
    threshold = 4
    df = pd.read_csv(ptm_file)
    columns_to_keep = [col for col in df.columns if df[col].iloc[0] in parameters.MODIFICATIONS.keys()]
    df_filtered = df[columns_to_keep]
    df_values = df_filtered.iloc[2:].astype(int)
    sums = df_values.sum()
    filtered_columns = sums[sums >= threshold].index
    filtered_df = df[filtered_columns]

    result_df = pd.concat([df.iloc[:, :2], filtered_df], axis=1)
    
    return result_df



def create_details_plot(input_file: str | os.PathLike, output_path: str | os.PathLike):
    if not 'A' in parameters.INPUT_FILES.keys():
        fig = sequence_plot.create_plot(input_file, 'A')
    elif not 'B' in parameters.INPUT_FILES.keys():
        fig = sequence_plot.create_plot(input_file, 'B')
    else:
        fig = sequence_plot.create_plot(input_file)
    for group in parameters.INPUT_FILES.keys():
        match parameters.INPUT_FILES[group][0]:
            case 'Cleavage':
                cleavage_file_path = parameters.INPUT_FILES[group][1]
                cleavage_group = group
            case 'PTM':
                ptm_file_path = parameters.INPUT_FILES[group][1]
                ptm_group = group
    if parameters.FIGURE_ORIENTATION == 0:
        # TODO calculate plot width based on legend to remove magic number 25
        plot_width = parameters.FIGURE_WIDTH-utils.SEQUENCE_BOUNDARIES['x0']-25
    else:
        pass

    label_plot_height = 150
    
    if cleavage_file_path:
        cleavage_df = pd.read_csv(cleavage_file_path)
        present_regions = get_present_regions_cleavage(cleavage_df)
        cleavages = cleavage_df.iloc[0:1,2:].values[0].tolist()
        pixels_per_cleavage = plot_width // (len(cleavages) + present_regions.count(True))
        assert(pixels_per_cleavage > parameters.FONT_SIZE)
        plot_cleavage_labels(fig, cleavage_df, pixels_per_cleavage, label_plot_height, cleavage_group)

    if ptm_file_path:
        ptm_df = filter_relevant_modification_sights(ptm_file_path)
        

    
    utils.show_plot(fig, output_path)

create_details_plot(parameters.FASTA_INPUT_FILE, parameters.OUTPUT_FOLDER)
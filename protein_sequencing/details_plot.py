from collections import defaultdict
import math
import os
import plotly.graph_objects as go
import pandas as pd

from protein_sequencing import parameters, utils,sequence_plot

def get_present_regions(positions):
    ranges = []
    for range in positions:
        if '-' in str(range):
            start, end = map(int, range.split('-'))
            ranges.append((start, end))
        else:
            start = end = int(range)
            ranges.append((start, end))

    region_ranges = []
    region_start = 1
    for _, region_end, _, _ in parameters.REGIONS:
        region_ranges.append((region_start, region_end))
        region_start = region_end + 1

    regions_present = [False] * len(region_ranges)
    region_index = 0
    for range in ranges:
        while range[0] > region_ranges[region_index][1]:
            region_index += 1
        regions_present[region_index] = True
    return regions_present

def get_present_regions_cleavage(cleavage_df: pd.DataFrame):
    cleavages = cleavage_df.iloc[0:1,2:].values[0].tolist()
    return get_present_regions(cleavages)

def get_present_regions_ptm(ptm_df: pd.DataFrame):
    ptms = ptm_df.iloc[1:2,2:].values[0].tolist()
    ptms = [ptm[1:] for ptm in ptms]
    return get_present_regions(ptms)
    

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

def plot_neuropathologies_horizontal(fig, df, x_0_neuropathologies, y_0_neuropathologies, dx, dy, x_label, y_label, last_region):
    x_margin = 0
    if dx % 2 != 0:
        x_margin = 1
    fig.add_shape(type='rect',
                            x0=x_0_neuropathologies - dx//2 - x_margin,
                            y0=y_0_neuropathologies,
                            x1=x_0_neuropathologies + dx * len(df.iloc[0:1,:].columns) - dx//2,
                            y1=y_0_neuropathologies + dy * len(df.index) + 1,
                            fillcolor='grey',
                            line=dict(color='grey', width=1),
                            showlegend=False,
                            layer='below',)
    fig.add_trace(go.Heatmap(z=df,
                             x0=x_0_neuropathologies,
                             y0=y_0_neuropathologies+dy//2,
                             dx=dx, dy=dy,
                             showscale=False, hoverinfo='none',
                             xgap=1, ygap=1,
                             colorscale=[[0, '#FFFFFF'], [1, '#452A79']]))
    
    fig.add_annotation(x=x_label, y=y_label,
            text=parameters.REGIONS[last_region][3],
            showarrow=False,
            font=dict(
                family=parameters.FONT,
                size=parameters.SEQUENCE_PLOT_FONT_SIZE,
                color='black',
                ))
    return fig

def plot_cleavage_labels(fig: go.Figure, cleavage_df: pd.DataFrame, pixels_per_cleavage: int, label_plot_height: int, group: str):
    # prepare data:
    cleavages = cleavage_df.iloc[0:1,2:].values[0].tolist()
    cleavage_df.columns = cleavage_df.iloc[0]
    cleavage_df = cleavage_df.iloc[1:]

    reverse_neuropathology_mapping = {}
    for key, list in parameters.NEUROPATHOLOGIES.items():
        values = list[0]
        for value in values:
            reverse_neuropathology_mapping[value] = key

    cleavage_df.iloc[:,1] = cleavage_df.iloc[:,1].map(reverse_neuropathology_mapping)
    mean_values = cleavage_df.iloc[:,2:].astype(float).groupby(cleavage_df.iloc[:,1]).mean()

    longest_label = ''
    for cleavage in cleavages:
        if utils.get_label_length(str(cleavage)) > utils.get_label_length(longest_label):
            longest_label = str(cleavage)

    group_direction = 1 if group == 'A' else -1
    first_cleavage_in_region = 0
    cleavages_visited = 0
    last_end = parameters.REGIONS[0][1]
    last_region = 0

    y_0_line = utils.SEQUENCE_BOUNDARIES['y1'] if group == 'A' else utils.SEQUENCE_BOUNDARIES['y0']
    y_1_line = y_0_line + 10 * group_direction
    y_2_line = y_0_line + (label_plot_height - utils.get_label_length(longest_label) - 10) * group_direction

    y_0_neuropathologies = y_0_line + (label_plot_height + 10) * group_direction
    vertical_space_left = utils.get_height() - y_0_neuropathologies if group == 'A' else y_0_neuropathologies
    # offset for border around heatmap
    vertical_space_left -= 2
    # offset for label
    vertical_space_left -= utils.get_label_height() + 5
    dx = pixels_per_cleavage
    dy = vertical_space_left//len(mean_values.index)*group_direction

    for i, neuropathology in enumerate(mean_values.index):
        y_0_rect = y_0_neuropathologies + i*dy
        fig.add_shape(type='rect',
                      x0 = utils.SEQUENCE_OFFSET - dx//2 - 100,
                      x1 = utils.SEQUENCE_OFFSET - dx//2,
                      y0 = y_0_rect,
                      y1 = y_0_rect + dy,
                      fillcolor=parameters.NEUROPATHOLOGIES[neuropathology][1],
                      line=dict(width=0),
                      showlegend=False,
                      layer='below',)
        # based on https://stackoverflow.com/questions/3942878/
        red, green, blue = tuple(int(parameters.NEUROPATHOLOGIES[neuropathology][1][i:i+2], 16) for i in (1, 3, 5))

        fig.add_annotation(x=utils.SEQUENCE_OFFSET - dx//2 - 50, y=y_0_rect + dy//2,
            text=neuropathology,
            showarrow=False,
            align='center',
            font=dict(
                family=parameters.FONT,
                size=parameters.SEQUENCE_PLOT_FONT_SIZE,
                color= '#000000' if red*0.299 + green*0.587 + blue*0.114 > 130 else '#ffffff'))


    for cleavage in cleavages:
        if '-' in str(cleavage):
            start, end = map(int, cleavage.split('-'))
        else:
            start = end = int(cleavage)
        if start > last_end:
            if parameters.FIGURE_ORIENTATION == 0:
                x_0_neuropathologies = first_cleavage_in_region * pixels_per_cleavage + utils.SEQUENCE_OFFSET
                x_divider = cleavages_visited * pixels_per_cleavage + utils.SEQUENCE_OFFSET
                x_label = x_0_neuropathologies + (x_divider-x_0_neuropathologies)//2 - dx//2
                y_label = y_0_neuropathologies + len(mean_values.index)*dy + (5+utils.get_label_height()//2) * group_direction

                plot_neuropathologies_horizontal(fig, mean_values.iloc[:,first_cleavage_in_region-last_region:cleavages_visited-last_region], x_0_neuropathologies, y_0_neuropathologies, dx, dy, x_label, y_label, last_region)                
                
                fig.add_trace(go.Scatter(x=[x_divider,x_divider],
                            y=[y_0_neuropathologies, y_0_neuropathologies+len(mean_values.index)*dy],
                            mode='lines',
                            line=dict(color="black", width=3), showlegend=False, hoverinfo='none'))
            
            else:
                pass
            while start > last_end:
                last_region += 1
                last_end = parameters.REGIONS[last_region][1]
            cleavages_visited += 1
            first_cleavage_in_region = cleavages_visited
        if parameters.FIGURE_ORIENTATION == 0:
            if start == end:
                label = str(start)
                x_0_line = start * utils.PIXELS_PER_PROTEIN + utils.SEQUENCE_OFFSET
                x_1_line = cleavages_visited * pixels_per_cleavage + utils.SEQUENCE_OFFSET
                y_3_line = y_0_line + (label_plot_height - utils.get_label_length(label)) * group_direction
                y_label = y_3_line + (utils.get_label_length(label) // 2 + 5) * group_direction
                
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
                y_3_line = y_0_line + (label_plot_height - utils.get_label_length(label)) * group_direction
                y_label = y_3_line + (utils.get_label_length(label) // 2 + 5) * group_direction

                plot_range_with_label_horizontal(fig,
                                    x_0_start_line, x_0_end_line, x_1_line,
                                    y_0_line, y_1_line, y_2_line, y_3_line,
                                    y_label,
                                    label)
        else:
            pass
        cleavages_visited += 1
    # plot neuropathologies for last region
    if parameters.FIGURE_ORIENTATION == 0:
        x_0_neuropathologies = first_cleavage_in_region * pixels_per_cleavage + utils.SEQUENCE_OFFSET
        region_length = len(mean_values.iloc[0:1,first_cleavage_in_region-last_region:].columns)
        x_label = x_0_neuropathologies + (region_length * pixels_per_cleavage)//2 - dx//2
        y_label = y_0_neuropathologies + len(mean_values.index)*dy + (5+utils.get_label_height()//2) * group_direction
        plot_neuropathologies_horizontal(fig, mean_values.iloc[:,first_cleavage_in_region-last_region:], x_0_neuropathologies, y_0_neuropathologies, dx, dy, x_label, y_label, last_region)
    else:
        pass

def plot_line_with_label(fig: go.Figure, x_0: int, x_1: int, y_0: int, y_1: int, y_2: int, y_3: int, y_label: int, label: str, color):
    fig.add_trace(go.Scatter(x=[x_0, x_0, x_1, x_1],
                            y=[y_0, y_1, y_2, y_3],
                            mode='lines',
                            line=dict(color=color, width=1), showlegend=False, hoverinfo='none'))
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

def plot_ptm_labels(fig: go.Figure, ptm_df: pd.DataFrame, pixels_per_ptm: int, label_plot_height: int, group: str, second_row: bool):
    group_direction = 1 if group == 'A' else -1
    ptms = ptm_df.iloc[1:2,2:].values[0].tolist()
    label_length = utils.get_label_length(ptms[-1])
    if second_row:
        second_row_offset = pixels_per_ptm // 2

    if parameters.FIGURE_ORIENTATION == 0:
        y_0_line = utils.SEQUENCE_BOUNDARIES['y1'] if group == 'A' else utils.SEQUENCE_BOUNDARIES['y0']
        y_1_line = y_0_line + 10 * group_direction
        y_2_line = y_0_line + (label_plot_height - label_length - 10) * group_direction
        if second_row:
            y_2_line = y_0_line + (label_plot_height - 2*(label_length + 10)) * group_direction
    else:
        pass

    ptms_visited = 0
    last_end = parameters.REGIONS[0][1]
    last_region = 0
    region_boundary = False
    for i, ptm in enumerate(ptms):
        ptm_position = int(ptm[1:])
        while ptm_position > last_end:
            last_region += 1
            last_end = parameters.REGIONS[last_region][1]
            region_boundary = True
        if region_boundary:
            region_boundary = False
            ptms_visited += 1
        if parameters.FIGURE_ORIENTATION == 0:
            x_0_line = ptm_position * utils.PIXELS_PER_PROTEIN + utils.SEQUENCE_OFFSET
            x_1_line = ptms_visited * pixels_per_ptm + utils.SEQUENCE_OFFSET
            y_3_line = y_2_line + 10 * group_direction
            if second_row and i % 2 == 1:
                x_1_line = (ptms_visited-1) * pixels_per_ptm + utils.SEQUENCE_OFFSET + second_row_offset
                y_3_line = y_2_line + (label_length + 10 + 5) * group_direction
            else:
                ptms_visited += 1
            y_label = y_3_line + (utils.get_label_length(ptm)+10) // 2 * group_direction
            plot_line_with_label(fig, x_0_line, x_1_line, y_0_line, y_1_line, y_2_line, y_3_line, y_label, ptm, parameters.MODIFICATIONS[str(ptm_df.iloc[0,i+2])][1])
        else:
            pass


def filter_relevant_modification_sights(ptm_file: str, threshold: int):
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
    cleavage_file_path = None
    ptm_file_path = None
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
        ptm_df = filter_relevant_modification_sights(ptm_file_path, parameters.MODIFICATION_THRESHOLD)
        present_regions = get_present_regions_ptm(ptm_df)
        number_of_ptms = len(ptm_df.columns)
        number_of_regions = present_regions.count(True)
        second_row = False
        if (number_of_ptms + number_of_regions) * parameters.FONT_SIZE <= plot_width:
            pixels_per_ptm = plot_width // (number_of_ptms + number_of_regions)
        elif (number_of_ptms + 2*number_of_regions) * parameters.FONT_SIZE <= 2 * plot_width:
            pixels_per_ptm = plot_width // (math.ceil(number_of_ptms/2) + number_of_regions)
            second_row = True
        else:
            raise ValueError('Too many PTMs to fit in plot')
        
        plot_ptm_labels(fig, ptm_df, pixels_per_ptm, label_plot_height, ptm_group, second_row)
    
    utils.show_plot(fig, output_path)

create_details_plot(parameters.FASTA_INPUT_FILE, parameters.OUTPUT_FOLDER)
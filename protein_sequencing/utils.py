from collections import defaultdict
import os
import plotly.graph_objects as go

import numpy as np
from protein_sequencing import parameters

SEQUENCE_BOUNDARIES = (None, None, None, None)
PIXELS_PER_PROTEIN = 0
SEQUENCE_OFFSET = 0

def get_width():
    if parameters.FIGURE_ORIENTATION == 0:
        return parameters.FIGURE_WIDTH
    return parameters.FIGURE_HEIGHT

def get_height():
    if parameters.FIGURE_ORIENTATION == 0:
        return parameters.FIGURE_HEIGHT
    return parameters.FIGURE_WIDTH

def get_left_margin():
    return parameters.LEFT_MARGIN * get_width()

def get_top_margin(): 
    return parameters.TOP_MARGIN * get_height()

def get_right_margin():
    return parameters.RIGHT_MARGIN * get_width()

def get_bottom_margin():
    return parameters.BOTTOM_MARGIN * get_height()

def get_label_length(label):
    if parameters.FIGURE_ORIENTATION == 1:
        return parameters.FONT_SIZE/1.5 * len(label) + 4
    return parameters.FONT_SIZE/1.5 * len(label)

def get_label_height():
    return parameters.FONT_SIZE+parameters.FONT_SIZE/4

def separate_by_group(groups_by_position):
    group_a = defaultdict(list)
    group_b = defaultdict(list)

    for protein_position in groups_by_position.keys():
        modification_sights = groups_by_position[protein_position]
        for modification_sight in modification_sights:
            if modification_sight[2] == 'A':  # group A
                group_a[protein_position].append(modification_sight)
            else:  # group B
                group_b[protein_position].append(modification_sight)
    
    return group_a, group_b

def different_possibilities_plot(width: int, height: int, different_possibilities: list[int]):
    rectangle = np.zeros((height, width))
    for i, value in enumerate(different_possibilities):        
        rectangle[:, i] = value
    fig = go.Figure(data=go.Heatmap(z=rectangle))
    fig.show()

def get_modifications_per_position(input_file):
    with open(input_file, 'r') as f:
      rows = f.readlines()[1:3]
      modification_types = rows[0].strip().split(',')
      labels = rows[1].strip().split(',')

      modifications_by_position = defaultdict(list)
      for i, (label) in enumerate(labels):
         if label == '':
            continue
         position = int(label[1:])
         letter = label[0]
         if letter in parameters.EXCLUDED_MODIFICATIONS:
               if parameters.EXCLUDED_MODIFICATIONS[letter] is None:
                  continue
               if modification_types[i] in parameters.EXCLUDED_MODIFICATIONS[letter]:
                  continue
         if modification_types[i] not in parameters.MODIFICATIONS:
            continue
         modifications_by_position[position].append((label, modification_types[i], parameters.MODIFICATIONS[modification_types[i]][2]))
      for position, mods in modifications_by_position.items():
         modifications_by_position[position] = list(set(mods))
    return modifications_by_position

def clean_up():
    directory = 'data/tmp'

    files = os.listdir(directory)

    for file_name in files:
        file_path = os.path.join(directory, file_name)
        if os.path.isfile(file_path):
            os.remove(file_path)

def show_plot(fig, output_path):
    output_file = f"{output_path}/figure1.png"
    fig.show()
    fig.write_image(output_file)

    clean_up()

    return output_file
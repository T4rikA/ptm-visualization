from collections import defaultdict
import os
import plotly.graph_objects as go
import pandas as pd
from queue import PriorityQueue

from protein_sequencing import parameters, utils,sequence_plot


def map_positions_to_bars(positions, num_bars, num_positions):
   # calculate gaps between positions
   min_pos = 1
   max_pos = parameters.REGIONS[-1][1]
   gaps = [min_pos] + positions + [max_pos]
   for i in range(len(gaps)-1):
      gaps[i] = gaps[i+1] - gaps[i]
   gaps = gaps[:-1]
   gaps_to_position = PriorityQueue()
   for i, gap in enumerate(gaps[:-1]):
      gaps_to_position.put((gap*-1, positions[i], positions[i]))
   gaps_to_position.put((gaps[-1]*-1, -1, -1))

   number_of_gaps = num_bars - num_positions

   for i in range(number_of_gaps-1):
      largest_gap = gaps_to_position.get()
      left_gap = (largest_gap[0] // 2, -1, largest_gap[2])
      right_gap = (largest_gap[0] - left_gap[0], largest_gap[1], largest_gap[2])
      gaps_to_position.put(left_gap)
      gaps_to_position.put(right_gap)

   gaps_for_position = {}
   bars_to_positions = []
   end_empty_bars = []

   while gaps_to_position.empty() == False:
      gap, pos, left_pos = gaps_to_position.get()
      if left_pos == -1:
         end_empty_bars.append(None)
      elif pos == -1:
         gaps_for_position[left_pos] = gaps_for_position.get(left_pos, 0) +  1    

   for pos in positions:
      if pos in gaps_for_position:
         empty_bars = [None]* gaps_for_position[pos]
         gaps_for_position[pos] = 0
         bars_to_positions.extend(empty_bars)
      bars_to_positions.append(pos)
   
   bars_to_positions.extend(end_empty_bars)

   return bars_to_positions

def add_labels_to_plot(fig: go.Figure, modification_sights: dict[int, list[tuple[int, str, str]]], group: str, bar_plot_label_height: int) -> go.Figure:
   group_direction = 1 if group == 'A' else -1
   group_size = 0

   positions = []
   for protein_position in modification_sights.keys():
      for modification_sight in modification_sights[protein_position]:
         positions.append(protein_position)
         group_size += 1

   if parameters.FIGURE_ORIENTATION == 0:
      bar_plot_width = utils.get_width() - utils.get_label_length('100%') - bar_plot_label_height
      num_bars = int(bar_plot_width // parameters.BAR_WIDTH)
      assert num_bars >= group_size
      map_positions_to_bars(positions, num_bars, group_size)
   else:
      pass

   height_offset = 0
   max_offset = 0
   lines = []
   modifications_visited = 0
   for protein_position in modification_sights.keys():
      for modification_sight in modification_sights[protein_position]:
         label, modification_type, modification_group = modification_sight

         if parameters.FIGURE_ORIENTATION == 0:
            x_start_line = protein_position * utils.PIXELS_PER_PROTEIN + utils.SEQUENCE_OFFSET
            x_end_line = x_start_line
            y_start_line = utils.SEQUENCE_BOUNDARIES['y1'] if group == 'A' else utils.SEQUENCE_BOUNDARIES['y0']
            y_end_line = y_start_line + group_direction * height_offset
            plot_line_with_label(fig, x_start_line, x_end_line, y_start_line, y_end_line, parameters.MODIFICATIONS[modification_type][1])
         else:
            pass
         modifications_visited += 1


         if modifications_visited > group_size/2:
            height_offset -= 2
         else:
            height_offset += 2
            max_offset = max(max_offset, height_offset)


   

   for line in lines:
      plot_line_with_label(fig, x_start_line, x_end_line, y_start_line, y_end_line, parameters.MODIFICATIONS[modification_type][1])
         
   
   
   return fig

def plot_line_with_label(fig, x_start, x_end, y_start, y_end, color):
    fig.add_trace(go.Scatter(x=[x_start, x_end], y=[y_start, y_end], mode='lines', line=dict(color=color, width=1), showlegend=False, hoverinfo='none'))

def filter_relevant_modification_sights(modification_sights: dict[int, list[tuple[int, str, str]]], helper_file: str):
   # read csv into df
   df = pd.read_csv(helper_file)

   # only keep first two columns and columns that are in MODIFICATIONS
   columns_to_keep = list(parameters.MODIFICATIONS.keys())
   df = df[[col for col in df.columns if df[col][0] in columns_to_keep or col in ['ID', 'Neuropathology']]]

   # remove rows where Neuropathology is not in parameters.Neurpathologies
   header_rows = df.iloc[:4, :]
   df = df[df['Neuropathology'].isin(parameters.NEUROPATHOLOGIES)]
   df = pd.concat([header_rows, df], ignore_index=True)

   # only keep columns where there are occurrences for the reasearched pathalogies
   data_to_filter = df.iloc[4:, 2:]
   header_columns = df.iloc[:, :2]
   for col in data_to_filter.columns:
      data_to_filter[col] = data_to_filter[col].astype(int)
   df = df[[col for col in data_to_filter if data_to_filter[col].sum(axis=0) > 0]]
   df = pd.concat([header_columns, df], axis=1)


   # create new dict for modification sights
   relevant_modification_sights = defaultdict(list)
   for column in df.columns:
      if column in ['ID', 'Neuropathology']:
         continue
      column_modification = df[column][0]
      column_position = int(df[column][1][1:])
      relevant_modification_sights[column_position].append((df[column][1], column_modification, 'A'))
   return relevant_modification_sights, df


def create_label_plot(input_file: str | os.PathLike, output_path: str | os.PathLike):
   fig = sequence_plot.create_plot(input_file)

   label_input_file = 'data/chris/label_plot/PP-MASCOT-CellAll-All_cutoff_0-05FDR_TAU_reformat_reduced_sub_binaryCell.csv'

   modification_sights = utils.get_modifications_per_position(label_input_file)

   relevant_positions, df = filter_relevant_modification_sights(modification_sights, label_input_file)

   group_a, group_b = utils.separate_by_group(relevant_positions)
   for (group, group_dict) in [('A', group_a), ('B', group_b)]:
      fig = add_labels_to_plot(fig, group_dict, group,  100)
   

   utils.show_plot(fig, output_path)

create_label_plot(parameters.FASTA_INPUT_FILE, parameters.OUTPUT_FOLDER)
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
   end_empty_bars = []

   bars_to_positions = []
   position_to_bar_position = {}

   while gaps_to_position.empty() == False:
      gap, pos, left_pos = gaps_to_position.get()
      if left_pos == -1:
         end_empty_bars.append(None)
      elif pos == -1:
         gaps_for_position[left_pos] = gaps_for_position.get(left_pos, 0) +  1    
   
   bar_counter = 0
   for pos in positions:
      if pos in gaps_for_position:
         empty_bars = [None]* gaps_for_position[pos]
         bar_counter += gaps_for_position[pos]
         gaps_for_position[pos] = 0
         bars_to_positions.extend(empty_bars)
      position_to_bar_position[pos] = bar_counter
      bars_to_positions.append(pos)
      bar_counter += 1
   
   bars_to_positions.extend(end_empty_bars)

   return bars_to_positions, position_to_bar_position

def add_labels_to_plot(fig: go.Figure, modification_sights: dict[int, list[tuple[int, str, str]]], group: str, bar_plot_label_height: int) -> go.Figure:
   group_direction = 1 if group == 'A' else -1
   group_size = 0

   positions = []
   for protein_position in modification_sights.keys():
      for modification_sight in modification_sights[protein_position]:
         positions.append(protein_position)
         group_size += 1
   bar_positions = []
   position_to_bar_position = {}
   if parameters.FIGURE_ORIENTATION == 0:
      num_bars = (utils.get_width() - utils.SEQUENCE_BOUNDARIES["x0"]) // parameters.BAR_WIDTH
      assert num_bars >= group_size
      bar_positions, position_to_bar_position  = map_positions_to_bars(positions, num_bars, group_size)
   else:
      pass

   height_offset = 0
   max_offset = 0
   modifications_visited = 0
   # for duplicate positions
   bar_offset = 0
   for protein_position in modification_sights.keys():
      for i, modification_sight in enumerate(modification_sights[protein_position]):
         label, modification_type, modification_group = modification_sight

         if parameters.FIGURE_ORIENTATION == 0:
            # x position for protein sequence
            x_0_line = protein_position * utils.PIXELS_PER_PROTEIN + utils.SEQUENCE_OFFSET
            # x position for bar plot
            x_1_line = (position_to_bar_position[protein_position]+i+bar_offset) * parameters.BAR_WIDTH + utils.SEQUENCE_BOUNDARIES["x0"] - parameters.BAR_WIDTH//2

            y_0_line = utils.SEQUENCE_BOUNDARIES['y1'] if group == 'A' else utils.SEQUENCE_BOUNDARIES['y0']
            y_1_line = y_0_line + group_direction * height_offset
            space_above_sequence = utils.get_height()-y_0_line if group == 'A' else y_0_line
            y_3_line = y_0_line + space_above_sequence // 5 * group_direction
            y_2_line = y_3_line - 10 * group_direction
            y_label = y_3_line + (utils.get_label_length(label) // 2 + 5) * group_direction 
            plot_line_with_label(fig,
                                 x_0_line, x_1_line,
                                 y_0_line, y_1_line, y_2_line, y_3_line,
                                 y_label,
                                 parameters.MODIFICATIONS[modification_type][1],
                                 label)
         else:
            pass
         modifications_visited += 1
         bar_offset += i

      if modifications_visited > group_size/2:
         height_offset -= 2
      else:
         height_offset += 2
         max_offset = max(max_offset, height_offset)

   
   return fig

def plot_line_with_label(fig, x_0, x_1, y_0, y_1, y_2, y_3, y_label, color, label):
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
   

def filter_relevant_modification_sights(helper_file: str):
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
      column_label = df[column][1][0]
      if column_modification in parameters.EXCLUDED_MODIFICATIONS.get(column_label, []):
         continue
      column_position = int(df[column][1][1:])
      relevant_modification_sights[column_position].append((df[column][1], column_modification, 'A'))
   return relevant_modification_sights, df


def create_label_plot(input_file: str | os.PathLike, output_path: str | os.PathLike):
   fig = sequence_plot.create_plot(input_file)

   label_input_file = 'data/chris/label_plot/PP-MASCOT-CellAll-All_cutoff_0-05FDR_TAU_reformat_reduced_sub_binaryCell.csv'

   relevant_positions, df = filter_relevant_modification_sights(label_input_file)

   group_a, group_b = utils.separate_by_group(relevant_positions)
   for (group, group_dict) in [('A', group_a), ('B', group_b)]:
      fig = add_labels_to_plot(fig, group_dict, group,  100)
   

   utils.show_plot(fig, output_path)

create_label_plot(parameters.FASTA_INPUT_FILE, parameters.OUTPUT_FOLDER)
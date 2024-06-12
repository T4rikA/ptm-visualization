from collections import defaultdict
import os
import plotly.graph_objects as go
import pandas as pd
from queue import PriorityQueue

from protein_sequencing import parameters, utils,sequence_plot


def add_bar_plot(fig: go.Figure, group: str, modification_sights_all: dict[int, list[tuple[int, str, str]]], modification_sights_relevant: dict[int, list[tuple[int, str, str]]], df) -> go.Figure:
   group_direction = 1 if group == 'A' else -1
   group_size = 0
   positions = []
   for protein_position in modification_sights_all.keys():
      for modification_sight in modification_sights_all[protein_position]:
         positions.append(protein_position)
         group_size += 1
   if group_size == 0:
      return fig

   if parameters.FIGURE_ORIENTATION == 0:
      max_num_bars = (utils.get_width() - utils.SEQUENCE_BOUNDARIES["x0"]) // parameters.MIN_BAR_WIDTH
      assert max_num_bars >= group_size
      bar_width = (utils.get_width() - utils.SEQUENCE_BOUNDARIES["x0"]) // group_size
   else:
      pass

   height_offset = 0
   max_offset = 0
   modifications_visited = 0
   # TODO ask chris or henne if there is better way to do this
   bar_percentages = {neuropathology: [] for neuropathology in parameters.NEUROPATHOLOGIES.keys()}
   legends = {'A': {neuropathology: False for neuropathology in parameters.NEUROPATHOLOGIES.keys()},
              'B': {neuropathology: False for neuropathology in parameters.NEUROPATHOLOGIES.keys()}}

   for protein_position in modification_sights_all.keys():
      for modification_sight in modification_sights_all[protein_position]:
         if modification_sight not in modification_sights_relevant[protein_position]:
            for neuropathology in parameters.NEUROPATHOLOGIES.keys():
               bar_percentages[neuropathology].append(0)
            modifications_visited += 1
            continue
         label, modification_type, modification_group = modification_sight
         if parameters.FIGURE_ORIENTATION == 0:
            # x position for protein sequence
            x_0_line = protein_position * utils.PIXELS_PER_PROTEIN + utils.SEQUENCE_OFFSET
            # x position for bar plot
            x_1_line = modifications_visited * bar_width + bar_width//2 + utils.SEQUENCE_BOUNDARIES["x0"] 
            y_0_line = utils.SEQUENCE_BOUNDARIES['y1'] if group == 'A' else utils.SEQUENCE_BOUNDARIES['y0']
            y_1_line = y_0_line + group_direction * height_offset
            space_above_sequence = utils.get_height()-y_0_line if group == 'A' else y_0_line
            space_per_group = space_above_sequence // len(df["Neuropathology"].unique())
            y_3_line = y_0_line + (space_per_group  - (utils.get_label_length(label) + 5)) * group_direction
            y_2_line = y_3_line - 10 * group_direction
            y_label = y_3_line + (utils.get_label_length(label) // 2 + 5) * group_direction

            # plot line with label
            plot_line_with_label(fig,
                                 x_0_line, x_1_line,
                                 y_0_line, y_1_line, y_2_line, y_3_line,
                                 y_label,
                                 parameters.MODIFICATIONS[modification_type][1],
                                 label)
            y_bar = y_3_line + (utils.get_label_length(label) + 5) * group_direction
            # -10 for margin above and below bars
            max_bar_height = space_per_group - 10
            # plot bars
            for i, neuropathology in enumerate(parameters.NEUROPATHOLOGIES.keys()):
               modification_column = [col for col in df.columns if col not in ['ID', 'Neuropathology'] and df[col].iloc[1] == label]
               bar_df = df[['ID', 'Neuropathology'] + modification_column]
               bar_df = bar_df[bar_df['Neuropathology'] == neuropathology]
               values = bar_df.iloc[:, 2].astype(int)
               percentage = values.mean()
               height = max_bar_height * percentage
               bar_percentages[neuropathology].append(percentage)
               x0 = x_1_line - bar_width // 2
               x1 = x_1_line + bar_width // 2
               y0 = y_bar + (i * space_per_group + 5) * group_direction
               y1 = y0 + height * group_direction

               fig.add_shape(
                  type="rect",
                  x0=x0+1,
                  y0=y0,
                  x1=x1-1,
                  y1=y1,
                  line=dict(color="black", width=1),
                  fillcolor=parameters.MODIFICATIONS[modification_type][1]
               )

               # plot legends
               if not legends[modification_group][neuropathology]:
                  legends[modification_group][neuropathology] = True
                  y_group = y0
                  fig.add_annotation(x = 10,
                                     y= y_group + space_per_group//2 * group_direction,
                                     width=space_per_group,
                                     text=parameters.NEUROPATHOLOGIES[neuropathology],
                                     textangle=-90,
                                     font=dict(size=parameters.SEQUENCE_PLOT_FONT_SIZE, color="black"),
                                     showarrow=False)
                  for j in range(5):
                     y_trace = y_group + max_bar_height//4 * j * group_direction
                     fig.add_trace(go.Scatter(x=[x0, utils.get_width()],
                                              y=[y_trace, y_trace],
                                              mode='lines',
                                              line=dict(color='lightgray', width=1),
                                              showlegend=False,
                                              hoverinfo='none')
                                              )
                     if j % 2 == 0:
                        text = str(j*25) + '%'
                        fig.add_annotation(x=x0-5-utils.get_label_length(text),
                                           y=y_trace,
                                           width=space_per_group,
                                           text=text,
                                           font=dict(size=parameters.SEQUENCE_PLOT_FONT_SIZE, color="black"),
                                           showarrow=False)
                     

         else:
            pass
         modifications_visited += 1

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
   return fig

def filter_relevant_modification_sights(helper_file: str):
   # read csv into df
   df = pd.read_csv(helper_file)

   # only keep first two columns and columns that are in MODIFICATIONS
   columns_to_keep = list(parameters.MODIFICATIONS.keys())
   df = df[[col for col in df.columns if df[col][0] in columns_to_keep or col in ['ID', 'Neuropathology']]]

   # create dict for all modification sights
   all_modification_sights = defaultdict(list)
   for column in df.columns:
      if column in ['ID', 'Neuropathology']:
         continue
      column_modification = df[column][0]
      column_label = df[column][1][0]
      if column_modification in parameters.EXCLUDED_MODIFICATIONS.get(column_label, []):
         continue
      column_position = int(df[column][1][1:])
      all_modification_sights[column_position].append((df[column][1], column_modification, parameters.MODIFICATIONS[column_modification][2]))

   # remove rows where Neuropathology is not in parameters.Neurpathologies
   header_rows = df.iloc[:4, :]
   df = df[df['Neuropathology'].isin(parameters.NEUROPATHOLOGIES.keys())]
   df = pd.concat([header_rows, df], ignore_index=True)

   # filter out columns that have no modifications observed for relevant neuropathologies
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
      relevant_modification_sights[column_position].append((df[column][1], column_modification, parameters.MODIFICATIONS[column_modification][2]))
   return all_modification_sights, relevant_modification_sights, df


def create_label_plot(input_file: str | os.PathLike, output_path: str | os.PathLike):
   fig = sequence_plot.create_plot(input_file)

   label_input_file = 'data/chris/label_plot/PP-MASCOT-CellAll-All_cutoff_0-05FDR_TAU_reformat_reduced_sub_binaryCell.csv'
   all_positions, relevant_positions, df = filter_relevant_modification_sights(label_input_file)
   
   group_a_all, group_b_all = utils.separate_by_group(all_positions)
   group_a_relevant, group_b_relevant = utils.separate_by_group(relevant_positions)
   for (group_label, group_all, group_relevant) in [('A', group_a_all, group_a_relevant), ('B', group_b_all, group_b_relevant)]:
      fig = add_bar_plot(fig, group_label, group_all, group_relevant, df)
   

   utils.show_plot(fig, output_path)

create_label_plot(parameters.FASTA_INPUT_FILE, parameters.OUTPUT_FOLDER)
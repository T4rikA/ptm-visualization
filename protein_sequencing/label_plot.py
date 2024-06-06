from collections import defaultdict
import os
import plotly.graph_objects as go
import pandas as pd

from protein_sequencing import parameters, utils,sequence_plot


def add_labels_to_plot(fig: go.Figure, input_file: str | os.PathLike) -> go.Figure:
   
   return fig


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
   fig = add_labels_to_plot(fig, input_file)

   label_input_file = 'data/chris/label_plot/PP-MASCOT-CellAll-All_cutoff_0-05FDR_TAU_reformat_reduced_sub_binaryCell.csv'

   modification_sights = utils.get_modifications_per_position(label_input_file)

   relevant_positions, df = filter_relevant_modification_sights(modification_sights, label_input_file)

   group_a, group_b = utils.separate_by_group(relevant_positions)
   for (group, group_dict) in [('A', group_a), ('B', group_b)]:
      pass

   

   utils.show_plot(fig, output_path)

create_label_plot(parameters.FASTA_INPUT_FILE, parameters.OUTPUT_FOLDER)
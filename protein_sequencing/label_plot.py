import os
import plotly.graph_objects as go

from protein_sequencing import utils,sequence_plot


def add_labels_to_plot(fig: go.Figure, input_file: str | os.PathLike) -> go.Figure:
   
   return fig

    


def create_label_plot(input_file: str | os.PathLike, output_path: str | os.PathLike):
    fig = sequence_plot.create_plot(input_file)
    fig = add_labels_to_plot(fig, input_file)


    utils.show_plot(fig, output_path)
import plotly.graph_objects as go
import os
from protein_sequencing import utils
import uniprot_align
import parameters

def create_plot(input_file: str | os.PathLike) -> go.Figure:

    alignments = list(uniprot_align.get_alignment(input_file))
    max_sequence_length = 0
    for alignment in alignments:
        if len(alignment.seq) > max_sequence_length:
            max_sequence_length = len(alignment.seq)

    assert all(len(alignment.seq) == max_sequence_length for alignment in alignments)

    different_possibilities = [-1]*max_sequence_length
    for i in range(len(alignments[0].seq)):
        proteins = set()
        for alignment in alignments:
            protein = alignment.seq[i]
            proteins.add(protein)
        
        if '-' in proteins:
            if len(proteins) == 2:
                different_possibilities[i] = -1
            if len(proteins) > 2:
                different_possibilities[i] = len(proteins)-1
        else:
            different_possibilities[i] = len(proteins)

    # TODO: needed later for different starts/ endings
    count, i = 0, 0
    while i < len(different_possibilities):
        if different_possibilities[i] == 2:
            count += 1
            while i < len(different_possibilities) and different_possibilities[i+1] == 2:
                i += 1
        i += 1

    # For debugging purposes
    # different_possibilities_plot(max_sequence_length, 50, different_possibilities)

    # basis for all pixel calculations
    if parameters.FIGURE_ORIENTATION == 0:
        max_sequence_length_pixels = utils.get_width() - utils.get_left_margin() - utils.get_right_margin()
        utils.PIXELS_PER_PROTEIN = int(max_sequence_length_pixels // max_sequence_length)
        utils.SEQUENCE_OFFSET = utils.get_left_margin()
    else:
        max_sequence_length_pixels = utils.get_height() - utils.get_top_margin() - utils.get_bottom_margin()
        utils.PIXELS_PER_PROTEIN = int(max_sequence_length_pixels // max_sequence_length)
        utils.SEQUENCE_OFFSET = utils.get_top_margin()

    # calculate region boundaries in pixels
    region_boundaries = []
    region_end_pixel = utils.SEQUENCE_OFFSET
    region_start = 1
    for region_name, region_end, region_color in parameters.REGIONS:
        region_start_pixel = region_end_pixel
        region_end_pixel = region_end * utils.PIXELS_PER_PROTEIN + 1 + utils.SEQUENCE_OFFSET
        region_boundaries.append((region_name, region_start_pixel, region_end_pixel, parameters.SEQUENCE_REGION_COLORS[region_color], region_start, region_end))
        region_start = region_end + 1

    fig = create_sequence_plot(region_boundaries)
    
    return fig

def create_sequence_plot(region_boundaries: list[tuple[str, int, int, str, int, int]]) -> go.Figure:
    fig = go.Figure()
    
    width = utils.get_width()
    height = utils.get_height()

    # General Layout
    fig.update_layout(
        title="",
        width = width,
        height = height,
        xaxis=dict(range=[0, width], autorange=False),
        yaxis=dict(range=[0, height], autorange=False),
        plot_bgcolor="white",
        font_family=parameters.FONT,
    )
    fig.update_xaxes(visible=False)
    fig.update_yaxes(visible=False)

    # Legend
    x_legend = 0 if parameters.FIGURE_ORIENTATION == 0 else width/2 + parameters.SEQUENCE_PLOT_HEIGHT
    if parameters.FIGURE_ORIENTATION == 1:
        fig.add_trace(go.Scatter(x=[x_legend], y=[height], mode='text', text="<b>Legend:</b>", textposition="bottom right", showlegend=False, hoverinfo='none', textfont=dict(size=parameters.SEQUENCE_PLOT_FONT_SIZE, color="black")))
    for i, modification in enumerate(parameters.MODIFICATIONS.values()):
        y_legend = height/2 + parameters.SEQUENCE_PLOT_HEIGHT + i*utils.get_label_height()
        if parameters.FIGURE_ORIENTATION == 1:
            y_legend = height - ((i+1)*utils.get_label_height())
        fig.add_trace(go.Scatter(x=[x_legend], y=[y_legend], mode='text', text=modification[0], textposition="bottom right", showlegend=False, hoverinfo='none', textfont=dict(size=parameters.SEQUENCE_PLOT_FONT_SIZE, color=modification[1])))
    if parameters.FIGURE_ORIENTATION == 0:
        fig.add_trace(go.Scatter(x=[x_legend], y=[y_legend + utils.get_label_height()], mode='text', text="<b>Legend:</b>", textposition="bottom right", showlegend=False, hoverinfo='none', textfont=dict(size=parameters.SEQUENCE_PLOT_FONT_SIZE, color="black")))
    # Sequence
    fig = plot_sequence(fig, region_boundaries)

    return fig

def plot_sequence(fig, region_boundaries):
    for i, (region_name, region_start_pixel, region_end_pixel, region_color, region_start, region_end) in enumerate(region_boundaries):
        

        if parameters.FIGURE_ORIENTATION == 0:
            x0 = region_start_pixel
            x1 = region_end_pixel
            y0 = utils.get_height()/2 - parameters.SEQUENCE_PLOT_HEIGHT/2
            y1 = y0+parameters.SEQUENCE_PLOT_HEIGHT
        else:
            y0 = utils.get_height() - region_start_pixel
            y1 = utils.get_height() - region_end_pixel
            x0 = utils.get_width()/2 - parameters.SEQUENCE_PLOT_HEIGHT/2
            x1 = x0+parameters.SEQUENCE_PLOT_HEIGHT
            

        # Region rects
        fig.add_shape(
            type="rect",
            x0=x0,
            y0=y0,
            x1=x1,
            y1=y1,
            line=dict(color="darkgrey", width=2),
            fillcolor=region_color
        )

        # Labels
        x_label = (x0 + x1) / 2
        y_label = (y0 + y1) / 2
        fig.add_annotation(
            x=x_label,
            y=y_label,
            text=region_name,
            showarrow=False,
            font=dict(size=parameters.SEQUENCE_PLOT_FONT_SIZE, color="black"),
            textangle= 90 if parameters.FIGURE_ORIENTATION == 1 else 0
        )

        if i == 0:
            if parameters.FIGURE_ORIENTATION == 0:
                x = x0 - utils.get_label_length(str(region_start))
                y = y_label
            else:
                x = x_label
                y = y0 + utils.get_label_height()
            fig.add_annotation(
                x=x,
                y=y,
                text='1',
                showarrow=False,
                font=dict(size=parameters.SEQUENCE_PLOT_FONT_SIZE, color="gray"),
                textangle= 0
            )
    if parameters.FIGURE_ORIENTATION == 0:
        x = x1 + utils.get_label_length(str(region_end))
        y = y
    else:
        x = x
        y = y1 - utils.get_label_height()
    fig.add_annotation(
        x=x,
        y=y,
        text= region_end,
        showarrow=False,
        font=dict(size=parameters.SEQUENCE_PLOT_FONT_SIZE, color="gray"),
        textangle= 0
    )

    utils.SEQUENCE_BOUNDARIES = {'x0': x0, 'x1': x1, 'y0': y0, 'y1': y1}
    return fig

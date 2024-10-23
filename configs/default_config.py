# General Settings

# Sequence Settings
# First sequence is from (1, 44), second from (45, 73) and so on
# Region Name, Region End, Group, Region Abbreviation
# REGIONS = [
#     ('N-term', 44, 'A', 'N-term'),
#     ('N1', 73, 'B', 'N1'),
#     ('N2', 102, 'B', 'N2'),
#     ('2N4R-Tau', 150, 'A', 'Mid'),
#     ('Proline-rich region', 241, 'B', 'PRR'),
#     ('R1', 272, 'B', 'R1'),
#     ('R2', 303, 'B', 'R2'),
#     ('R3', 334, 'B', 'R3'),
#     ('R4', 371, 'B', 'R4'),
#     ('C-term', 441, 'A', 'C-term'),
# ]
REGIONS = [
    ("start", 390, "A", "start"),
    ("alpha", 432, "B", "a"),
    ("epsilon", 431, "A", "e"),
]

# Modification Settings
MODIFICATION_LEGEND_TITLE = 'PTMs'
MODIFICATIONS = {
    'Phospho': ('Phosphorylation', '#000000'),
    'Acetyl': ('Acetylation', '#93478F'),
    'Methyl': ('Methylation', '#C35728'),
    'GG': ('Ubiquitination', '#548056'),
    'Citrullination': ('Citrullination', '#FF17E3'),
}

EXCLUDED_MODIFICATIONS = {'Q': None,
                          'X': None,
                          'S': ['GG'],}

# Input Output Settings
FASTA_INPUT_FILE = 'protein_sequencing/data_preprocessing/data/uniprot_analysis/output_aligned_2/P14136_aligned.fasta'
OUTPUT_FOLDER = 'output'

# Plot Settings
FIGURE_ORIENTATION = 0  # 0 for horizontal, 1 for vertical, note figure height and width are then automatically swapped

# just change width and height to change the size of the figure not the orientation
FIGURE_WIDTH = 1500
FIGURE_HEIGHT = 1000

FONT_SIZE = 12

PTMS_TO_HIGHLIGHT = ['Phospho(S)@61', 'Citrullination(R)@242', 'GG(K)@254', 'Acetyl(K)@267', 'Methyl(K)@311']
PTM_HIGHLIGHT_LABEL_COLOR = '#cfcfcf'


# Default Parameters
FONT = 'Arial'

# Margins for sequence Plot
# TODO remove margins and auto calculate based on legend
LEFT_MARGIN = 0.065
RIGHT_MARGIN = 0.025
TOP_MARGIN = 0.065
BOTTOM_MARGIN = 0.025

# Sequence Plot
SEQUENCE_PLOT_FONT_SIZE = FONT_SIZE
SEQUENCE_PLOT_HEIGHT = 50
# works best with an even number
EXONS_GAP = 10

# Sequence Region Colors
SEQUENCE_REGION_COLORS = {
    'A': 'white',
    'B': 'lightgrey',
}
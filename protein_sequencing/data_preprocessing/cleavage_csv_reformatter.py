import pandas as pd
import numpy as np

file = "data/chris/cleavage_plot/PPc_COMPLETE_cutoff_0-05FDR_reformat_XX_C_collmean.csv"

def extract_numeric_value(value):
    return int(value[1][1:])

df = pd.read_csv(file)

new_data_2 = {}

new_data_2['ID'] = df.iloc[:, 0].tolist()
new_data_2['Neuropathology'] = df.iloc[:, 1].tolist()

# new_data = {}
# for col in df.columns[6:]:
#     second_row_value = df.loc[0, col]    
#     sights = second_row_value.split(';')
#     for sight in sights:
#         mod_pos = sight.split('@')[1]
#         mod_type = sight.split('(')[1][0]
#         mod_name = sight.split('(')[0].strip()
#         data = [mod_name, mod_type+mod_pos]
#         data.extend(df[col][1:].tolist())
#         new_data[mod_name+"."+mod_type+mod_pos] = data

new_data = {}
for col in df.columns[6:]:
    second_row_value = df.loc[0, col]    
    sights = second_row_value.split(';')
    for sight in sights:
        mod_pos = sight.split('@')[1]
        mod_type = sight.split('@')[0]
        data = [mod_pos]
        data.extend(df[col][1:].tolist())
        new_data[mod_type+"."+mod_pos] = data


#sorted_dict = dict(sorted(new_data.items(), key=lambda item: extract_numeric_value(item[1])))
new_data_2.update(new_data)
new_df = pd.DataFrame(new_data_2)
new_df.to_csv("data/chris/cleavage_plot/PPc_COMPLETE_cutoff_0-05FDR_reformat_XX_C_collmean_tarik.csv", index=False)

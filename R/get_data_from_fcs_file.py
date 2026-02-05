# import libraries
import flowkit as fk
import pandas as pd
import seaborn as sns
import numpy as np
import sklearn
import matplotlib.pyplot as plt
import re as re


############## code execution ###############

#-------------------- import fcs files -----------------------
# change path based on folder containing the fcs files whose expression data need to be extracted
data_path = "/home/rstudio/marine_data/Aggregated data for paper/Pipeline/second_batch_2025/non-algae_gated/12022025_AB_fcs_non-algae"

# import a session of fcs files
fcs_session = fk.Session(fcs_samples=data_path)

# list file IDs, ie filenames from loaded group of files
list_files = fcs_session.get_sample_ids()

#--------------------- transform data ----------------

# a loop to transform all files in the list
# biex_xform = fk.transforms.WSPBiexTransform('biex', max_value=4194304, positive=5, width=-39, negative=0) # not valid anymore, probably due to update.
biex_xform = fk.transforms.WSPBiexTransform(max_value=4194304, positive=5, width=-39, negative=0)

fcs_session.get_sample(sample_id=list_files[1])

for i in range(len(list_files)): 
  print(f"Applying transform to sample {i}")
  fcs_session.get_sample(sample_id=list_files[i]).apply_transform(biex_xform, include_scatter=True)


#------ check samples channel old vs new samples

sample_id = fcs_session.get_sample_ids()[0]  # Or use your own ID
sample = fcs_session.get_sample(sample_id=sample_id)
df_expr = sample.as_dataframe(source='raw')

df_expr.columns

type(df_expr)

df_expr_old = pd.read_csv('~/marine_data/results/test_data_all_m.csv',nrows=5)

type(df_expr_old)

df_expr_old.columns

# You have two DataFrames:
# 
# 1. df_expr_old.columns ➜ a regular (flat) index
# Index([...], dtype='object')
# 2. df_expr.columns ➜ a MultiIndex, with 2 levels:
# MultiIndex([...], names=['pnn', 'pns'])
# They both look horizontal when printed, but structurally:
# 
# df_expr_old.columns is a 1D index of strings
# 
# df_expr.columns is a 2D index: ('FSC-A', ''), etc.
# 
# why This Happens
# FlowKit sometimes uses a MultiIndex for expression matrices, especially when:
# 
# You request transformed data via get_expression_data(transformed=True)
# 
# It wants to keep both Pnn (parameter name) and Pns (parameter short name) from the FCS metadata
# 
# So this:

# MultiIndex([('FSC-A', '')], names=['pnn', 'pns'])
# means:
# 
# "FSC-A" is the pnn (parameter name, e.g. $PnN)
# 
# '' is the pns (parameter short name, e.g. $PnS — sometimes blank)

# to flatten df_expr.columns

df_expr.columns = df_expr.columns.get_level_values('pnn') 

#------ checking scatter plot of data
plt.scatter(df_expr['B1-A'], df_expr['B2-A'], s=10, alpha=0.7)
plt.xlabel('FSC-A')
plt.ylabel('SSC-A')
plt.title('Simple Scatter Plot')
plt.show()


# Note: select 10 samples with a lot of events.
#----------------- export expression data to csv -------------------
# a loop to export data from all files in the list
# use source="raw" to get untransformed data or source="xform" for transformed data or if needed source="comp" for raw but compensated data
# in the loop it also adds the Species label based on the filename, the filename format should be RCCxxxx_rest of filename.fcs, then based on indexing it takes the name preceding                                                                               the "_" character
data_list = []
for i in range(len(list_files)): 
    df = fcs_session.get_sample(sample_id=list_files[i]).as_dataframe(source="xform").droplevel(1, axis="columns")
    df["Species"] = list_files[i][:list_files[i].index("_")]
    data_list.append(df)
# check if all files are imported

# merge into one frame
con_df = pd.concat(data_list, ignore_index=True)

# write into csv
data_path_out="/home/rstudio/marine_data/Aggregated data for paper/Pipeline/second_batch_2025/data_csv_all/"
data_name = "second_batch_non_algae_m3.csv"
con_df.to_csv(data_path_out+data_name, index=False)

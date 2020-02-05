import json
import os
import pandas as pd

inp_path = os.path.expanduser("~/DATA/TF_P-values/")

for file in os.listdir(inp_path):
    table = pd.read_table(inp_path + file)
    if table.empty:
        print(inp_path + file)
        os.remove(inp_path + file)

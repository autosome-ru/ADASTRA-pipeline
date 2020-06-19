import os
import sys
import pandas as pd

dirname = sys.argv[1]
ids_set = set()
for filename in os.listdir(dirname):
    df = pd.read_table(os.path.join(dirname, filename))
    ids_set |= set(df['ID'])

print(len(ids_set))

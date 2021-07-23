import pandas as pd
import os
from zipfile import ZipFile


alignments_path = os.path.expanduser('~/AlignmentsChip')
master_list = pd.read_table('~/Configs/master-chip.txt')
print(master_list)
master_list = master_list[master_list['TF_UNIPROT_NAME'] == 'CTCF_HUMAN']
master_list.to_csv('~/Configs/CTCF_master.tsv', sep='\t', index=False)
with ZipFile(os.path.expanduser('~/CTCF_dump.zip'), 'w') as zipObj2:
    for index, row in master_list.iterrows():
        path = os.path.join(alignments_path, row['#EXP'], row['ALIGNS'] + '.table_p')
        zipObj2.write(path, os.path.join('DATA', os.path.splitext(os.path.basename(path))[0] + '.tsv'))
    zipObj2.write(os.path.expanduser('~/Configs/CTCF_master.tsv'), 'CTCF_master.tsv')

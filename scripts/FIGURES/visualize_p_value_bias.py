import sys
import os
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns

stats = pd.read_table(os.path.expanduser('~/CTCF_HUMAN_common_table.tsv'))
print(stats)

fig, ax = plt.subplots(figsize=(10, 8))
x = np.array(range(1, 18))
stats = stats[stats['mean_sBAD'] == 1]
#stats = stats[stats['ID'] == '.']
palts = stats[stats['logitp_alt'] != 1.0]
sns.barplot(x=x, y=[len(palts[(1/10**(k+1) < palts['logitp_alt']) & (palts['logitp_alt'] <= 1/10**k)].index)
                    for k in x], label='alt', ax=ax, color='C1', alpha=0.5)

prefs = stats[stats['logitp_ref'] != 1.0]
sns.barplot(x=x, y=[len(prefs[(1/10**(k+1) < prefs['logitp_ref']) & (prefs['logitp_ref'] <= 1/10**k)].index)
                    for k in x], label='ref', ax=ax, color='C0', alpha=0.5)




plt.grid(True)
plt.legend()
plt.title('ref-alt p_value on BAD=1')
plt.xlabel('x: x+1 >-log10 p_value >= x')
plt.ylabel('snp count')
plt.savefig(os.path.expanduser('~/fixed_alt/p_dist.png'))
plt.show()

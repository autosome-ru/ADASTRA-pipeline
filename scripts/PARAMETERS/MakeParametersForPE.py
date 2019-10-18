import json

parameters_path = '/home/abramov/PARAMETERS/'
parallel_parameters_path = '/home/abramov/ParallelParameters/stats/'

out_path = parallel_parameters_path + 'PE_parameters.cfg'

cell_lines = [
    'K562_myelogenous_leukemia',
    'HCT-116_colon_carcinoma',
    'MCF7_Invasive_ductal_breast_carcinoma',
    'PC3_prostate_carcinoma',
    'LoVo_colorectal_adenocarcinoma',
    'HeLa_S3_cervical_adenocarcinoma',
    ]
# excluded = [
#    '_labs_michael-snyder___biosamples_ENCBS357NWO_',
#    'GSE109481',
#    ]
excluded = []

if __name__ == "__main__":
    with open(parameters_path + 'CELL_LINES.json', 'r') as read_file:
        d = json.loads(read_file.readline())

    keys = sorted(d.keys())
    with open(out_path, 'w') as file:
        ### FIXME: this part is only for correlation
        syn_path = '../CORRELATIONanalysis/synonims.tsv'
        names = []
        with open(syn_path, 'r') as syn:
            for line in syn:
                line = line.strip('\n').split('\t')
                if line[1]:
                    GTRD_name = line[0].replace('(', '').replace(')', '').replace(' ', '_')
                    if GTRD_name not in cell_lines:
                        continue
                    COSMIC_name = line[1]
                    names.append(GTRD_name)
        ###
        for key in keys:
            ### FIXME: this part is only for correlation
            try:
                name, lab = key.split('!')
            except ValueError:
                print(key)
            if name not in names or lab in excluded:
                continue
            ###
            file.write(key + '\n')

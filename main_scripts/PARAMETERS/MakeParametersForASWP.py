import os

if __name__ == "__main__":
    Ploidy_path = '/home/abramov/Ploidy/'
    out_path = '/home/abramov/ParallelParameters/ASWP_parameters.cfg'
    
    with open(out_path, 'w') as file:
        ### this part is only for correlation
        syn_path = '../CORRELATIONanalysis/synonims.tsv'
        names = []
        with open(syn_path, 'r') as syn:
            for line in syn:
                line = line.strip('\n').split('\t')
                if line[1] and line[2]:
                    GTRD_name = line[0].replace('(', '').replace(')', '').replace(' ', '_')
                    COSMIC_name = line[1]
                    names.append(GTRD_name)
        ###
        for file_name in os.listdir(Ploidy_path):
            if not os.path.isfile(Ploidy_path + file_name):
                continue
            ###
            name = file_name.split('!')[0]
            if name not in names:
                print('Skipping, not in synonims..')
                continue
            ###
            file.write(file_name + '\n')
    
    
import os

if __name__ == "__main__":
    Correlation_path = '/home/abramov/Correlation/'
    out_path = '/home/abramov/ParallelParameters/CS_parameters.cfg'
    
    with open(out_path, 'w') as file:
        for file_name in os.listdir(Correlation_path):
            if os.path.isdir(Correlation_path + file_name):
                ### FIXME: this part is only for correlation
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
                for file_name2 in os.listdir(Correlation_path + file_name):
                    ### FIXME: this part is only for correlation
                    name = file_name2[:file_name2.rfind('_')]
                    if name not in names:
                        print('Skipping, not in synonims..')
                        continue
                    ###
            
                    file.write(file_name2 + '\n')
                break

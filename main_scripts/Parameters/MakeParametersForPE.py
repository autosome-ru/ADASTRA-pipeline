import json

if __name__ == "__main__":
    JSON_path = '/home/abramov/PLOIDYcalling/CELL_LINES.json'
    out_path = '/home/abramov/ParallelParameters/PE_parameters.cfg'
    
    with open(JSON_path, 'r') as read_file:
        d = json.loads(read_file.readline())
    keys = sorted(d.keys())
    
    with open(out_path, 'w') as file:
        for key in keys:
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
            name = key.split('!')[0]
            if name not in names:
                print('Skipping, not in synonims..')
                continue
            ###
            file.write(key + '\n')
    
    
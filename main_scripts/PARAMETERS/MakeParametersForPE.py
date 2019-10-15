import json

if __name__ == "__main__":
    JSON_path = '/home/abramov/PLOIDYcalling/CELL_LINES.json'
    out_path = '/home/abramov/ParallelParameters/PE_parameters.cfg'
    
    with open(JSON_path, 'r') as read_file:
        d = json.loads(read_file.readline())
    keys = sorted(d.keys())
    
    with open(out_path, 'w') as file:
        for key in keys:
            file.write(key + '\n')
    
    
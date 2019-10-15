import os

if __name__ == "__main__":
    Ploidy_path = '/home/abramov/Ploidy/'
    out_path = '/home/abramov/ParallelParameters/ASWP_parameters.cfg'
    
    with open(out_path, 'w') as file:
        for file_name in os.listdir(Ploidy_path):
            if not os.path.isfile(Ploidy_path + file_name):
                continue
            file.write(file_name + '\n')
    
    
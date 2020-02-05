import os

inp_path = os.path.expanduser("~/DATA/CL_DICTS/")
out_path = os.path.expanduser("~/DATA/CL_DICTS_new/")
good_path = os.path.expanduser("~/ParallelParameters/Agr_parameters.cfg")
good_one = set()
with open(good_path) as f:
    for line in f:
        good_one.add(line.strip("\n"))

for file in os.listdir(inp_path):
    print(file)
    if file.replace("_DICT.json", "") in good_one:
        os.rename(inp_path + file, out_path + file)

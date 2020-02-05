import json
import os

inp_path = os.path.expanduser("~/DATA/CL_DICTS/")

for file in os.listdir(inp_path):
    print(file)
    with open(inp_path + file) as f:
        d = json.loads(f.readline())
    if len(d.keys()) == 0:
        print(inp_path + file)
        #os.remove(inp_path + file)

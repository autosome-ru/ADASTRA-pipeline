import json

parameters_path = "/home/abramov/PARAMETERS/"
alignments_path = "/home/abramov/Alignments/"

def make_reverse_dict(dictionary):
    new_dict = {}
    for key in dictionary:
        paths = dictionary[key]
        for path in paths:
            if path.find("CTRL") == -1:
                new_dict[path] = key
    return new_dict


with open(parameters_path + "CELL_LINES.json", "r") as read_file:
    old_one = json.loads(read_file.readline())
new_one = make_reverse_dict(old_one)
with open(parameters_path + "REVERSE_CELL_LINES.json", "w") as write_file:
    json.dump(new_one, write_file, sort_keys=True)

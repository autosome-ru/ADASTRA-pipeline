import json

dict_path = "/home/abramov/PLOIDYcalling/"


def make_reverse_dict(dictionary):
    new_dict = {}
    for key in dictionary:
        values = dictionary[key]
        for value in values:
            v = value.split("/")
            if v[4] == "EXP" and v[5] != "None":
                new_dict[value] = key
    return new_dict


with open(dict_path + "CELL_LINES.json", "r") as read_file:
    old_one = json.loads(read_file.readline())
new_one = make_reverse_dict(old_one)
with open(dict_path + "REVERSE_CELL_LINES.json", "w") as write_file:
    json.dump(new_one, write_file, sort_keys=True)

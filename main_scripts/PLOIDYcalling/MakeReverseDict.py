import json


def make_reverse_dict(dictionary):
    new_dict = {}
    for key in dictionary:
        values = dictionary[key]
        for value in values:
            print(value)
            v = value.split("/")
            if v[4] == "EXP":
                new_dict[value] = key
    return new_dict


with open("/home/abramov/PLOIDYcalling/CELL_LINES.json", "r") as read_file:
    old_one = json.loads(read_file.readline())
new_one = make_reverse_dict(old_one)
with open("/home/abramov/PLOIDYcalling/REVERSE_CELL_LINES.json", "w") as write_file:
    json.dump(new_one, write_file, sort_keys=True)

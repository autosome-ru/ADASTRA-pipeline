import requests
import json

alignments_path = "/home/abramov/Alignments/"
parameters_path = "/home/abramov/PARAMETERS/"


def findLAB(enc):
    if enc.find(";") != -1:
        enc = enc.split(";")[0]
    if enc.find("wgEncode") == -1:
        r = requests.get('https://www.encodeproject.org/experiments/' + enc + '/?format=json')
        lab = json.loads(r.text)['lab']['@id']
        biosample = json.loads(r.text)['replicates'][0]['library']['biosample']['@id']
        ret = lab + "_" + biosample
        return ret.replace("/", "_")
    else:
        return 'False'


def CreatePath(line, ctrl=False):
    if ctrl:
        return alignments_path + "CTRL/" + line[10] + "/" + line[14] + ".vcf.gz"
    else:
        return alignments_path + "EXP/" + line[1] + "/" + line[0] + "/" + line[6] + ".vcf.gz"


def add_to_dict(d, key, value):
    el = d.get(key, None)
    if el:
        d[key] = el | {value}
    else:
        d[key] = {value}


def add_record(d, line, ctrl=False):
    if ctrl:
        idcs = (12, 15, 16)
    else:
        idcs = (4, 8, 9)

    path = CreatePath(line, ctrl)
    line[idcs[0]] = line[idcs[0]].replace("(", "").replace(")", "").replace(" ", "_").replace("/", "_")
    if line[idcs[2]] != "None":
        Lab = findLAB(line[idcs[2]])
        if Lab:
            key = line[idcs[0]] + '!' + Lab
            add_to_dict(d, key, path)
    elif line[idcs[1]] != "None":
        key = line[idcs[0]] + '!' + line[idcs[1]]
        add_to_dict(d, key, path)


def MakeDict(masterList):
    master = open(masterList, "r")
    d = dict()
    count = 0
    for line in master:
        if line[0] == "#":
            continue
        count += 1
        if count % 10 == 0:
            print("Made {} Experiments out of ~6120".format(count))
        if count == 20:
            break
        ln = line.strip().split("\t")
        add_record(d, ln)

        if len(ln) > 16:
            add_record(d, ln, ctrl=True)
    print("Saving Dictionary")
    for key in d:
        value = d[key]
        sorted_value = sorted(list(value))
        d[key] = sorted_value
    with open(parameters_path + "CELL_LINES.json", "w") as write_file:
        json.dump(d, write_file)
    print("Dictionary Saved")


MakeDict(parameters_path + "Master-lines.tsv")

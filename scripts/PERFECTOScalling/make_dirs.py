import sys
import os

sys.path.insert(1, "/home/abramov/ASB-Project")
from scripts.HELPERS.helpers import pack
from scripts.HELPERS.paths_for_components import PWMs_path

processed_tfs = dict()

for filename in os.listdir(PWMs_path):
    filename = filename.split(".")
    tf_name = filename[0]
    rank = int(filename[2])
    quality = filename[3]
    if tf_name in processed_tfs.keys():
        best_rank, best_quality = processed_tfs[tf_name]
        if best_rank < rank:
            continue
        elif best_rank > rank:
            processed_tfs[tf_name] = (rank, quality)
            os.replace(PWMs_path + filename, PWMs_path + tf_name + "/" + filename)
        elif best_quality < quality:
            continue
        else:
            processed_tfs[tf_name] = (rank, quality)
            os.replace(PWMs_path + filename, PWMs_path + tf_name + "/" + filename)
    else:
        if not os.path.isdir(PWMs_path + tf_name):
            os.mkdir(PWMs_path + tf_name)
        os.replace(PWMs_path + filename, PWMs_path + tf_name + "/" + filename)
        processed_tfs[tf_name] = (rank, quality)


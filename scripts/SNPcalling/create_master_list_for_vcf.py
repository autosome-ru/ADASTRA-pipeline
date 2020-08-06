import os
from scripts.HELPERS.paths_for_components import configs_path
from scripts.HELPERS.helpers import make_list_for_VCFs


if __name__ == "__main__":
    make_list_for_VCFs(os.path.join(configs_path, "MasterListForVCFs.tsv"))

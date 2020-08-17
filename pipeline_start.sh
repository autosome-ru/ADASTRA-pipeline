#!/bin/bash

get_stage_index() {
  case "$1" in
    --create_reference) return 1
      break;;
    --snp_call) return 2
      break;;
    --peak_call) return 3
      break;;
    --bad_call) return 4
      break;;
    --nb_fit) return 5
      break;;
    --p_value_count) return 6
      break;;
    --aggregate_p_values) return 7
      break;;
    *)
      echo "There is no option $1"
      return 0
	    break;;
	esac
}

njobs=$1
flag=$2
python3 construct_parameters_python.py

source HELPERS/paths_for_components.py

stage_index = $( get_stage_index $2 )

if [ "$flag" -le 1 ]; then
  bash "$scripts_path/SNPcalling/CreateReference.sh" --RefFolder "$reference_path" --RefGenome "$genome_path"
  python3 "$scripts_path/PARAMETERS/make_badmaps_dict.py"
  python3 "$scripts_path/SNPcalling/"sort_columns.py
  python3 "$scripts_path/PARAMETERS/create_initial_dirs.py"
  python3 "$scripts_path/PARAMETERS/make_aggregation_dict.py" TF
  python3 "$scripts_path/PARAMETERS/make_aggregation_dict.py" CL
fi

if [ "$flag" -le 2 ]; then
  bash "$scripts_path/"snp_calling.sh "$njobs"
fi

if [ "$flag" -le 3 ]; then
  bash "$scripts_path/"annotation.sh "$njobs"
fi

if [ "$flag" -le 4 ]; then
  bash "$scripts_path"/bad_map_est.sh "$njobs" --merge
  bash "$scripts_path"/BAD_annotation.sh "$njobs"
fi

if [ "$flag" -le 5 ]; then
  python3 "$scripts_path"/FITnoise/collect_ref_bias_statistics.py
  python3 "$scripts_path"/FITnoise/fit_negative_binom_with_weights.py
fi

if [ "$flag" -le 6 ]; then
  bash "$scripts_path"/p_value_count.sh "$njobs"
fi

if [ "$flag" -le 7 ]; then
  bash "$scripts_path"/aggregation.sh "$njobs" --forTF
  bash "$scripts_path"/aggregation.sh "$njobs" --forCL
fi


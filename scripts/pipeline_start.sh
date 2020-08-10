#!/bin/bash

get_stage_index() {
  case "$1" in
    --create_reference) return 1
      break;;
    --snp_call) return 2
      break;;
    --peak_call) return 3
      break;;
    --bad_groups) return 4
      break;;
    --bad_call) return 5
      break;;
    --nb_fit) return 6
      break;;
    --p_value_count) return 7
      break;;
    --aggregate_p_values) return 8
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

if [ "$flag" -ge 1 ]; then
  snp_calling_path="$scripts_path"/SNPcalling/
  bash "$snp_calling_path"CreateReference.sh --RefFolder "$reference_path" --RefGenome "$genome_path"
  python3 "$snp_calling_path"create_master_list_for_vcf.py
fi

if [ "$flag" -ge 2 ]; then
  bash "$scripts_path/"snp_calling.sh "$njobs"
fi

if [ "$flag" -ge 3 ]; then
  bash "$scripts_path/"annotation.sh "$njobs"
fi

if [ "$flag" -ge 4 ]; then
  python3 "$scripts_path/"PARAMETERS/MakeDictForBadMaps.py
fi

if [ "$flag" -ge 5 ]; then
  bash "$scripts_path"/bad_map_est.sh "$njobs" --merge
  bash "$scripts_path"/BAD_annotation.sh "$njobs"
fi

if [ "$flag" -ge 6 ]; then
  python3 "$scripts_path"/PARAMETERS/MakeDictForAggregation.py TF
  python3 "$scripts_path"/PARAMETERS/MakeDictForAggregation.py CL
  python3 "$scripts_path"/FITnoise/collect_ref_bias_statistics.py
  python3 "$scripts_path"/FITnoise/fit_negative_binom_with_weights.py
fi

if [ "$flag" -ge 7 ]; then
  bash "$scripts_path"/p_value_count.sh "$njobs"
fi

if [ "$flag" -ge 9 ]; then
  bash "$scripts_path"/aggregation.sh "$njobs" --forTF
  bash "$scripts_path"/aggregation.sh "$njobs" --forCL
fi


#!/bin/bash

njobs=$1
flag=$2
start_script_path="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"

previous_pwd=$PWD
cd $start_script_path

python3 "construct_parameters_python.py"

source "scripts/HELPERS/paths_for_components.py"

case "$2" in
  --create_reference) stage_index=1
    ;;
  --snp_call) stage_index=2
    ;;
  --peak_call) stage_index=3
    ;;
  --bad_call) stage_index=4
    ;;
  --nb_fit) stage_index=5
    ;;
  --p_value_count) stage_index=6
    ;;
  --aggregate_p_values) stage_index=7
    ;;
  *)
    echo "There is no option $1"
    stage_index=0
    ;;
esac

if [ "$stage_index" -le 1 ]; then
  bash "$scripts_path/create_reference.sh" -RefFolder "$reference_path" -RefGenome "$genome_path"
  adastra badmaps_dict
  adastra sort_cols
  adastra init_dirs
  adastra aggregation_dict
fi

if [ "$stage_index" -le 2 ]; then
  bash "$scripts_path/"snp_calling.sh "$njobs"
fi

if [ "$stage_index" -le 3 ]; then
  bash "$scripts_path/"annotation.sh "$njobs"
fi
echo 'Annotated successfully'
if [ "$stage_index" -le 4 ]; then
  bash "$scripts_path"/bad_map_est.sh "$njobs" --merge
  bash "$scripts_path"/BAD_annotation.sh "$njobs"
fi

if [ "$stage_index" -le 5 ]; then
  adastra collect_ref_bias
  adastra fit_neg_bin
fi

if [ "$stage_index" -le 6 ]; then
  bash "$scripts_path"/p_value_count.sh "$njobs"
fi

if [ "$stage_index" -le 7 ]; then
  bash "$scripts_path"/aggregation.sh "$njobs" --forTF
  bash "$scripts_path"/aggregation.sh "$njobs" --forCL
fi

cd $previous_pwd
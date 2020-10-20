#!/bin/bash

njobs=$1
flag=$2
start_script_path="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"

previous_pwd=$PWD
cd $start_script_path

python3 "construct_parameters_python.py"

source "scripts/HELPERS/paths_for_components.py"

case "$flag" in
  --create-reference) stage_index=1
    ;;
  --snp-call) stage_index=2
    ;;
  --peak-annotation) stage_index=3
    ;;
  --bad-call) stage_index=4
    ;;
  --nb-fit) stage_index=5
    ;;
  --pvalue-count) stage_index=6
    ;;
  --aggregate-pvalues) stage_index=7
    ;;
  *)
    echo "There is no option $2"
    stage_index=0
    ;;
esac

if [ "$stage_index" -le 1 ]; then
  if !  bash "$scripts_path/create_reference.sh" -RefFolder "$reference_path" -RefGenome "$genome_path"
  then
    echo 'Create reference failed'
    exit 1
  fi
  if ! adastra badmaps_dict
  then
    echo 'BADmaps dict failed'
    exit 1
  fi
  if ! adastra sort_cols
  then
    echo 'Sort columns failed'
    exit 1
  fi
  if ! adastra init_dirs
  then
    echo 'Create directories failed'
    exit 1
  fi
  if ! adastra aggregation_dict
  then
    echo 'Aggregation dict failed'
    exit 1
  fi
fi

if [ "$stage_index" -le 2 ]; then
  if ! bash "$scripts_path/"snp_calling.sh "$njobs"
  then
    echo 'SNPcalling failed'
    exit 1
  fi
fi

if [ "$stage_index" -le 3 ]; then
  if ! bash "$scripts_path/"annotation.sh "$njobs"
  then
    echo 'Peak annotation failed'
    exit 1
  fi
fi

if [ "$stage_index" -le 4 ]; then
  if ! bash "$scripts_path"/bad_map_est.sh "$njobs" --merge
  then
    echo 'BAD estimation failed'
    exit 1
  fi
  if ! bash "$scripts_path"/correlation_with_cosmic.sh "$njobs" --annotate
  then
    echo 'Correlation analysis failed'
    exit 1
  fi
  if ! bash "$scripts_path"/BAD_annotation.sh "$njobs"
  then
    echo 'BAD annotation failed'
    exit 1
  fi
fi

if [ "$stage_index" -le 5 ]; then
  if ! adastra collect_ref_bias
  then
    echo 'Collect statistics failed'
    exit 1
  fi
  if ! adastra fit_neg_bin
  then
    echo 'Fit negative binom failed'
    exit 1
  fi
fi

if [ "$stage_index" -le 6 ]; then
  if ! bash "$scripts_path"/p_value_count.sh "$njobs"
  then
    echo 'P-value computation failed'
    exit 1
  fi
fi

if [ "$stage_index" -le 7 ]; then
  if ! bash "$scripts_path"/aggregation.sh "$njobs" --forTF
  then
    echo 'TF aggregation failed'
    exit 1
  fi
  if ! bash "$scripts_path"/aggregation.sh "$njobs" --forCL
  then
    echo 'CL aggregation failed'
    exit 1
  fi
fi

cd $previous_pwd
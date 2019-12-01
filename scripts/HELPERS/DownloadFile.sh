#!/bin/bash

DownloadPath=$1
BamPath=$2

if ! cp :"$DownloadPath" "$BamPath"
then
  echo "Failed to download ALIGNEXP $EXP"
  exit 1
fi

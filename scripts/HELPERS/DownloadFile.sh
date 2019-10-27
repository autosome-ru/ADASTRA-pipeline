#!/bin/bash

DownloadPath=$1
BamPath=$2

if ! scp :"$DownloadPath" "$BamPath"
then
  echo "Failed to download ALIGNEXP $EXP"
  exit 1
fi

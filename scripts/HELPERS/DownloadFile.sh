#!/bin/bash

DownloadPath=$1
BamPath=$2

if ! scp :"$DownloadPath" "$BamPath"
then
  echo "Failed to download $DownloadPath in $BamPath"
  exit 1
fi

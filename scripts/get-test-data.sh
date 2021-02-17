#!/usr/bin/env bash

#
# Cell division detector look up table csv file
#
source_url="https://drive.google.com/file/d/1gWFjPxj7hIITToL3VbB2UrKBFVm_3fDl/view?usp=sharing"
destination="data/img.klb"
timeout_seconds=5;
mkdir -p data
curl -sSL -m $timeout_seconds -o $destination $source_url
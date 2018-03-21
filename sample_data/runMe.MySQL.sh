#!/bin/bash

set -ev

./__run_sample_pipeline.pl --align_assembly_config mysql.confs/alignAssembly.config --annot_compare_config mysql.confs/annotCompare.config --dbtype mysql



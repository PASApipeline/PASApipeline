#!/bin/bash

set -ev

./__run_sample_pipeline.pl --align_assembly_config sqlite.confs/alignAssembly.config --annot_compare_config sqlite.confs/annotCompare.config --dbtype sqlite



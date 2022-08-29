#!/bin/bash
for filename in gt_standard_logs/*.fasta; do
	iqtree -b 100 -nt AUTO -s $filename
done 

cd gt_standard_logs
mv *.treefile ../gt_output

cd ../gt_output
rename -e 's/.treefile/.tre/'  *.treefile
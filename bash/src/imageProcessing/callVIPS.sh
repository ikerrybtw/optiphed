#!/bin/bash

dire=$1; shift
regex=$1; shift
outdir=$1; shift
layer=$1; shift
#files=`ls $dire*$regex*svs`;
#files=`ls $dire*TCGA-[A-D]*svs`
#TCGA-[E-Z0-9]
files=`find $dire -type f -name *$regex*svs`

for x in ${files[@]}
    do
        while true; do
		runs=`ps -ef | grep 'vips' | grep -v grep | wc -l`
		[ $runs -lt 25 ] && break;
                echo $runs' VIPS processes running. Waiting to complete...'  
                sleep 20;
        done
	parent=$(dirname $x)'/';
	fileN=$(basename $x);	
	patient=${fileN%.svs}
	##Check if output file exists
        if [ -f $outdir$patient'_pyramid75.tif' ]
        then
                echo "Pyramidal tif file for $patient  exists. Skipped"
        	cmd1="echo working with existing $patient_pyramid75.tif"
	else	
		cmd1="vips im_vips2tiff "$x" "$outdir$patient"_pyramid75.tif:jpeg:75,tile:256x256,pyramid"
		echo $cmd1
	fi
	
	cmd2="vips im_copy "$outdir$patient"_pyramid75.tif:"$(($layer-1))" "$outdir$patient"_layer"$layer".tif"
	#cmd6="vips im_copy "$outdir$patient"_pyramid75.tif:5 "$outdir$patient"_layer6.tif"
	eval $cmd1 && eval $cmd2 &
	sleep 10
done

#!/bin/bash

indir=$1; shift; #TIF directory
layer=$1; shift; #'layer2';
stainType=$1; shift; #'HandE' #'IPOX'
outdir=$1; shift; ##Nuclei directory

#root='~/Projects'
root='/mnt/ix2/nandor/Projects';
cpipe='/ITHinferenceFromHandE/code/imageProcessing/'$stainType'_adaptive_'$layer'.cppipe';

files=`ls $indir/*$layer'.tif'`;
for f in $files; do
	ind=`echo $(basename $f) | sed 's/\.tif//g'`
	ind=$outdir/$ind
	mkdir $ind
	imgdirout=$ind'/cellProfilerOut';
	cp $f $ind/	
	
	##Check if output already exists
	if [ -f $imgdirout'/Nuclei.txt' ];
	   then
		continue;	
	fi

	##Mac OSX
	#sudo /Applications/CellProfiler.app/Contents/MacOS/CellProfiler -p $root$cpipe -i $ind -t tmp/ -o $imgdirout -c 1  1>log.cellProfiler 2>err.cellProfiler
	
	##tamago
	cmd="cellprofiler --jvm-heap-size=15g -p $root$cpipe -i $ind -t tmp2/ -o $imgdirout -c 1" 
	echo $cmd
	eval $cmd
done

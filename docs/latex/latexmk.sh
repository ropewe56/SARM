#!/bin/bash

basename=`basename "$0"`
dirname=`dirname "$0"`

cd $dirname

latex_input=sarm

outdir=build

latex_cmd=xelatex
bib_cmd=biber

if [ ! -d $outdir ]
then
    mkdir $outdir
fi

echo "$latex_cmd -output-directory=$outdir -synctex=1 -interaction=nonstopmode -file-line-error $latex_input"
$latex_cmd -output-directory=$outdir -synctex=1 -interaction=nonstopmode -file-line-error $latex_input

echo "$bib_cmd   $outdir/$latex_input"
$bib_cmd   $outdir/$latex_input

$latex_cmd -output-directory=$outdir -synctex=1 -interaction=nonstopmode -file-line-error $latex_input
$latex_cmd -output-directory=$outdir -synctex=1 -interaction=nonstopmode -file-line-error $latex_input

cp $outdir/sarm.pdf $dirname
cp $outdir/sarm.synctex.gz $dirname

#!/bin/bash

dirname="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

#basename=`basename "$0"`
#dirname=`dirname "$0"`

cd $dirname

latex_input=sarm

outdir=build

latex_cmd=xelatex
bib_cmd=biber

if [ ! -d $outdir ]
then
    mkdir $outdir
fi

echo $latex_cmd
t1=$(($(date +%s%N)/1000000))
$latex_cmd -output-directory=$outdir -synctex=1 -file-line-error $latex_input # > log/xelatex_stdout.log 2> log/xelatex_stderr.log
t2=$(($(date +%s%N)/1000000))
echo `expr $t2 - $t1` ms

echo "end"



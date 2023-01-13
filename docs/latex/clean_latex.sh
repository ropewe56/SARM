#!/bin/bash
dirname="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

set +e
rm ${dirname}/*.aux
rm ${dirname}/*.bcf
rm ${dirname}/*.bbl
rm ${dirname}/*.blg
rm ${dirname}/*.idx
rm ${dirname}/*.log
rm ${dirname}/*.out
rm ${dirname}/*.tdo
rm ${dirname}/*.toc
rm ${dirname}/*.run.xml
rm ${dirname}/*.synctex.gz
rm ${dirname}/build/*
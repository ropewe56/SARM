#!/bin/bash

dirname="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

cd $dirname
echo $dirname

$dirname/clean_latex.sh

echo "konsole --hold --workdir ${dirname} -e bash ./xelatex.sh"

konsole --hold --workdir ${dirname} -e bash ./xelatex.sh


#!/bin/bash

julia --math-mode=fast -O3 --check-bounds=no --inline=yes --threads 16 run_spectrum.jl

# top : show processes
# top -H  : show threads
# top -H -p <pid> . theraeds of process <pid>
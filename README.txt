SiRMS
-----
Simplex representation of molecular structure - a chemoinformatic tool for calculation of simplex descriptors

Installation
-----
Just copy all *.py and *.json files in the same directory.

Usage
-----
Basic example - returns all possible 2D simplex descriptors with vertexes labelled by element
sirms.py -in input.sdf -o output.txt
The following command returns identical results
sirms.py -in input.sdf -o output.txt -d elm -t all

Returns 2D simplex descriptors with vertexes labelled by element, which belonds for 3-11 topoligical types and don't contain hydrogen atoms. Calculation will be performed on 2 cores, if they are available.
sirms.py -in input.sdf -o output.txt -x -t extended -c 2

Help
-----
sirms.py -h

Notes
-----
1. JSON-file is a dictionary of precomputed simplexes to speed up the calculation. You can disable its usage by -n option in command line

License
-----
GPLv3
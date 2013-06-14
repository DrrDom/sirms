SiRMS
-----
Simplex representation of molecular structure - a chemoinformatic tool for calculation of simplex descriptors

Features
-----
1. Calculate 2D simplexes descriptors, which are the number of identical simplexes in a molecule. Simplex - tetraatomic fragment with fixed topology and stereoconfiguration.
2. Calculate 2D simplexes descriptors for mixtures of compounds.
3. Multi-threaded calculations of descriptors. Very useful at processing of big molecules (more than 80 atoms). Otherwise this option cannot be recommended.
4. Calculate 2D simplexes descriptors for specified fragments of molecules.

Installation
-----
Just copy all *.py and *.json files in the same directory.

Usage
-----
Basic example - returns all possible 2D simplex descriptors with vertexes labeled by element
sirms.py -in input.sdf -o output.txt
The following command returns identical results
sirms.py -in input.sdf -o output.txt -d elm -t all

Returns 2D simplex descriptors with vertexes labeled by element, which belongs for 3-11 topological types and don't contain hydrogen atoms. Calculation will be performed on 2 cores, if they are available.
sirms.py -in input.sdf -o output.txt -x -t extended -c 2

Help
-----
sirms.py -h

Notes
-----
1. You need to standardize your structures before simplex descriptors calculation.
2. JSON-file is a dictionary of pre-computed simplexes to speed up the calculation. You can disable its usage by -n option in command line

License
-----
GPLv3
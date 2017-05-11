SiRMS
-----
Simplex representation of molecular structure - a chemoinformatic tool for calculation of simplex (fragment) descriptors

Theory
-----
Simplex - tetraatomic fragment with fixed topology and stereo configuration.
Simplex descriptor - number of identical simplexes in a molecule.
From the version 1.0.0 fragments of any size are supported (not only simplexes).

Features
-----
1. Calculates descriptors for single compounds, mixtures, quasi-mixtures and reactions (all-in-one tool).
2. Supports incorporation of user-defined atomic properties (e.g. partial charges, H-bond donor/acceptor, etc) in generated descriptors.
3. Works extremely well in QSAR/QSPR modeling tasks.

Limitations
-----
This version calculates only 2D descriptors.

Installation
-----
No installation is needed. Clone or download and unzip and launch sirms.py from the command line.

Python version
-----
Only Python 3 is supported.

Usage
-----
Basic example - returns 2D simplex descriptors with vertexes labeled by element
sirms.py -in input.sdf -o output.txt
sirms.py -in input.rdf -o output.txt

Help
-----
sirms.py -h

Wiki
-----
See wiki pages for more information and examples - https://github.com/DrrDom/sirms/wiki.

Notes
-----
You need to standardize your structures before simplex descriptors calculation.

License
-----
BSD 3-clause

Reference
-----
1. Kuz’min, V. E.; Artemenko, A. G.; Polischuk, P. G.; Muratov, E. N.; Khromov, A. I.; Liahovskiy, A. V.; Andronati, S. A.; Makan, S. Y. Hierarchic System of QSAR Models (1D-4D) on the Base of Simplex Representation of Molecular Structure - Journal of Molecular Modelling, 2005, 11, 457-467. (one of the first description of single compounds representation)
2. Kuz’min, V. E.; Artemenko, A. G.; Muratov, E. N. Hierarchical QSAR technology based on the Simplex representation of molecular structure - Journal of Computer-Aided Molecular Design, 2008, Volume 22, Issue 6-7, 403-421. (single compounds representation)
3. Oprisiu, I.; Varlamova, E.; Muratov, E.; Artemenko, A.; Marcou, G.; Polishchuk, P.; Kuz'min, V.; Varnek, A. QSPR Approach to Predict Nonadditive Properties of Mixtures. Application to Bubble Point Temperatures of Binary Mixtures of Liquids. Molecular Informatics 2012, 31, 491-502. (old approach of mixtures representation)
4. Mokshyna, E.; Nedostup, V. I.; Polischuk, P. G.; Kuzmin, V. E. ‘Quasi-Mixture’ Descriptors for QSPR Analysis of Molecular Macroscopic Properties. The Critical Properties of Organic Compounds. Molecular Informatics 2014, 33, 647-654. (quasi-mixture representation)

Changes
-----
version 1.1.0 (26.07.2016)
1) new canonicalization algorithm was implemented that changed descriptors names (break compatibility with older versions)

version 1.1.1 (20.02.2017)
1) multiprocessing calculation of descriptors was implemented for single molecules and svm output format only (sparse format)
2) per atom fragmentation was added
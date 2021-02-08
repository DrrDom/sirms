# SiRMS

Simplex representation of molecular structure - a chemoinformatic tool for calculation of simplex (fragment) descriptors

#### Background

Simplex - tetraatomic fragment with fixed topology and stereo configuration.
Simplex descriptor - number of identical simplexes in a molecule.
From the version 1.0.0 fragments of any size are supported (not only simplexes).

#### Features

1. Calculate descriptors for single compounds, mixtures, quasi-mixtures and reactions (all-in-one tool).
2. Support incorporation of user-defined atomic properties (e.g. partial charges, H-bond donor/acceptor, etc) in generated descriptors.
3. Work well in QSAR/QSPR modeling tasks.

#### Limitations

This version calculates only 2D descriptors.

#### Installation

`pip install sirms`

#### Usage
Returns 2D simplex descriptors with vertexes labeled by element for single compounds  
`sirms.py -in input.sdf -o output.txt`  
Returns 2D simplex descriptors with vertexes labeled by element for reactions  
`sirms.py -in input.rdf -o output.txt`  

#### Help

`sirms.py -h`

####Wiki

See wiki pages for more information and examples - https://github.com/DrrDom/sirms/wiki.

#### Notes

Standardize structures before descriptors calculation.

#### License

BSD 3-clause

#### Reference

1. Kuz’min, V. E.; Artemenko, A. G.; Polischuk, P. G.; Muratov, E. N.; Khromov, A. I.; Liahovskiy, A. V.; Andronati, S. A.; Makan, S. Y. Hierarchic System of QSAR Models (1D-4D) on the Base of Simplex Representation of Molecular Structure - Journal of Molecular Modelling, 2005, 11, 457-467. (**one of the first descriptions of single compounds representation**)
2. Kuz’min, V. E.; Artemenko, A. G.; Muratov, E. N. Hierarchical QSAR technology based on the Simplex representation of molecular structure - Journal of Computer-Aided Molecular Design, 2008, Volume 22, Issue 6-7, 403-421. (**single compounds representation**)
3. Oprisiu, I.; Varlamova, E.; Muratov, E.; Artemenko, A.; Marcou, G.; Polishchuk, P.; Kuz'min, V.; Varnek, A. QSPR Approach to Predict Nonadditive Properties of Mixtures. Application to Bubble Point Temperatures of Binary Mixtures of Liquids. Molecular Informatics 2012, 31, 491-502. (**old approach of mixtures representation**)
4. Mokshyna, E.; Nedostup, V. I.; Polischuk, P. G.; Kuzmin, V. E. ‘Quasi-Mixture’ Descriptors for QSPR Analysis of Molecular Macroscopic Properties. The Critical Properties of Organic Compounds. Molecular Informatics 2014, 33, 647-654. (**quasi-mixture representation**)
5. Polishchuk, P.; Madzhidov, T.; Gimadiev, T.; Bodrov, A.; Nugmanov, R.; Varnek, A. Structure–reactivity modeling using mixture-based representation of chemical reactions. J Comput Aided Mol Des 2017, 31 (9), 829-839 (**reaction representation & new mixture representation**).

#### Changes

version 1.1.0 (26.07.2016)
1) New canonicalization algorithm was implemented that changed descriptors names (break compatibility with older versions)

version 1.1.1 (20.02.2017)
- multiprocessing calculation of descriptors was implemented for single molecules and svm output format only (sparse format)
- per atom fragmentation was added

version 1.2.0 (06.02.2021)
- reorganized as a Python package

version 1.2.1 (08.02.2021)
- fix README formatting
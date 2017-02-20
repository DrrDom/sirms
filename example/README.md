Files which are used for demonstration are located in `example` dir of the repository. Switch to this dir and run following scripts.  
  
Note: structures should be standardized first, because calculated descriptors are fragment ones they may vary with different representation of compounds structures.  
  
### Descriptors of single compounds  
  
1. Default call is used to calculate simplex descriptors of single compounds with atoms labelled by elements. Only tatraatomic fragments fully connected or consisting of two disconnected parts will be calculated.  
`../sirms.py -i molecules.sdf -o output01.txt`  
  
2. Calculate counts of fragments with 2 to 6 atoms labelled by element and consisting of maximum two disconnected parts. Hydrogens are suppressed. 
`../sirms.py -i molecules.sdf -o output02.txt --min_atoms 2 --max_atoms 6 -x`  
  
3. Calculated counts of fully connected fragments with 2 to 6 atoms labelled by element.  
`../sirms.py -i molecules.sdf -o output03.txt --min_atoms 2 --max_atoms 6 --max_components 1`  
  
4. Calculate counts of fragments with 2 to 6 atoms labelled by element, partial atomic charge, lipophilicity, hydrogen bonding and refractivity and consisting of up to two disconnected parts.  
  
* First let's generate atomic properties by using external program (Chemaxon cxcalc)  
`/utilities/calc_atomic_properties_chemaxon.py -i example/molecules.sdf -o example/molecules_labeled.sdf`  
Several additional fields appeared in the new sdf-file which contain values of atomic properties.  
  
* Generate counts of fragments with 2 to 6 atoms labelled by element, charge, lipophilicity, refractivity, hydrogen bonding. Fragments consist of up to two disconnected parts. Hydrogens are suppressed.  
`../sirms.py -i molecules_labeled.sdf -o output04.txt --min_atoms 2 --max_atoms 6 -a elm CHARGE HB LOGP REFRACTIVITY -x`  
  
### Quasi-mixture descriptors of single compounds  
  
1. Calculate quasi-mixture descriptors for single compounds. Fragments will contain from 2 to 6 atoms labelled by charge, lipophilicity and hydrogen bonding. Hydrogens are suppressed.  
`../sirms.py -i molecules_labeled.sdf -o output05.txt --min_atoms 2 --max_atoms 6 -a CHARGE HB LOGP -x -q`  
  
### Mixture descriptors  
  
1. Calculate mixture descriptors for mixtures specified in the external file `mixtures.txt`. Fragments will contain from 2 to 6 atoms labelled by element, charge. lipophilicity and hydrogen bonding and consisting of up to two disconnected parts. Only pair-wise interaction of components will be taken into account by default. Hydrogens are suppressed.  
`../sirms.py -i molecules_labeled.sdf -o output06.txt --min_atoms 2 --max_atoms 6 -a elm CHARGE HB LOGP -x -m mixtures.txt`  
  
2. The same as above considering components amounts as relative values and taking into account pair-wise and triple-wise interaction of components with self-association.  
`../sirms.py -i molecules_labeled.sdf -o output07.txt --min_atoms 2 --max_atoms 6 -a elm CHARGE HB LOGP -x -m mixtures.txt --max_mix_components 3 --mix_self_association --mix_type rel`  
  
3. Calculate ordered mixture descriptors when the role of components are known and the position of the component in the mixture description file determines its role. Fragments with 2 to 6 atoms with suppressed hydrogens will be generated.  
`../sirms.py -i molecules.sdf -o output08.txt --min_atoms 2 --max_atoms 6 -x -r -m mixtures.txt`  
    
### Reaction descriptors  
  
They are similar to mixture descriptors except only labelling by elements is available at this moment.  
  
1. Generate reaction descriptors concatenating reactant and product mixture descriptors. Fragments will contain 2 to 6 atoms labelled by elements.  
`../sirms.py -i e2.rdf -o output09.txt --min_atoms 2 --max_atoms 6`  

  
2. The same as above but reaction descriptors will be calculated as a difference between the mixture descriptors of products and the mixture descriptors of reactants.  
`../sirms.py -i e2.rdf -o output10.txt --min_atoms 2 --max_atoms 6 --reaction_diff`  


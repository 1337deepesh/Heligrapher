HELIGRAPHER: a tool for generating antimicrobial helices
==============

Copyright (C) 2015 Deepesh Nagarajan (1337deepesh@gmail.com)

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, version 3 of the License.

This program is distributed in the hope that it will be useful, but 
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU General Public License <http://www.gnu.org/licenses/>
for more details.

Concept
--------------

Heligrapher attempts to design a helix such that its residues intersect 
the most possible maximal-common-subgraphs from a database of antimicrobial 
helical peptide structures. 

Heligrapher reduces 3D alpha-helical antimicrobial proteins to cylindrical
2D graphs. The *graph* itself is defined from the alpha-helical structure. 
Every residue is a node, and the amino-acid is a label (G,A,V,etc).
A node can have a maximum of 4 edges:
- Residue to the LEFT (covalently bound)
- Residue to the RIGHT (covalently bound)
- Residue ABOVE (hydrogen bound, i -> i+4)
- Residue BELOW (hydrogen bound, i -> i-4)

I have created a database of known alpha-helical antimicrobial peptide structures.
I have reduced the database structures to their graphs. Heligrapher then uses simulated annealing to optimize a designed helix such that it shares the maximum 
possible subgraphs with all database graphs (maximal common subgraph optimization).

It's difficult to explain in words but I hope you get the idea.

Example Output
--------------
Use the following command (inside the Heligrapher directory):
**python Heligrapher.py -i input_data.txt -o output_folder -d Heligrapher_database -S Y -n 1000 -P Y**

**-n 1000** is just for demonstration. I recommend using **-n 10000** or higher.
An example of how the output should appear is given in the folder **example_output**

Adding to the Database
--------------

Feel free to add/remove antimicrobial helical structures to the Heligrapher 
database. It is located here:
**cd Heligrapher/Heligrapher_database/antimicrobial_helices**

Keep these points in mind:
- Remove all HETATM/ANISOU/DNA. Keep only protein coordinates ("ATOM  " only).
- Remove multiple chains (especially for NMR structures). Keep only 1 chain.
- Heligrapher only considers alpha-helices (using dssp). Therefore, you don't need
to trim your peptides of beta-hairpins or loops before use.

Troubleshooting
--------------

If Heligrapher doesn't work, it may be because **dssp** doesn't have permissions.
To fix this:
cd Heligrapher/Heligrapher_database/third_party_software
sudo chmod 755 dssp-2.0.4-linux-amd64

The **-P** argument needs **pymol** installed. Either install pymol or use **-P N**.
To install pymol:
**sudo apt-get install pymol**

I have kept dependencies to a bare minimum. I have written my own simulated 
annealing protocols, my own protein structural-manipulation algos, and my own 
graphing/plotting functions. All to spare you the trouble of searching for 
libraries (you're welcome).

Contact
--------------
If you have any queries/suggestions, contact me:
Deepesh Nagarajan: 1337deepesh@gmail.com

.. -*- mode: rst -*-

Calculates gene ontology distance for a pair of UniProt/GeneOntology Identifiers using the formula::

	2*N3/(N1+N2+2*N3)

where,::
	N1 = distance of first term from common parent
	N2 = distance of second term from common parent
	N3 = distance of common parent from the root

To run the script::
	$ python go_distance.py input.txt

Sample output looks like this (depending on the input)::

	GO:0001578 GO:0030036 0.222222222222
	O04630 O13297 0.428571428571
	O01482 Q12072 0


***Note: Pointers for the common parent query from http://biostar.stackexchange.com/questions/10961/common-parent-between-go-terms
# symmetries
symmetry music

usage: python symmetries.py GROUP_NAME POLYPHONY

#

right now we have: dihedral1, dihedral2, diherdral3, dihedral4, dihedral5, and qubit

#

POLYPHONY is an integer. beware of setting it too high!

#

also, stupidly, it will crash if POLYPHONY is > than # of transformations in ensemble

#

Ps. there's an element of spagetti code in all this to be fixed in version two, which will involve refactoring everything to handle more general cases, and using better libraries.

# 

requires pyo 

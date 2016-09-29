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

#  Brief Summary

We represent the symmetries of, say, an equilateral triangle (3 flips, 3 rotations) by transformation matrices, and use these to generate sound by deriving the frequency of the next notes from the eigenvalues of the matrix and deriving the amplitudes from the projection of the prior state of the triangle onto the corresponding eigenvector. Elements of the symmetry group are assigned to keys on a midi keyboard, and there's a sound and light show.

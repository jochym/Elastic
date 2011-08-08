Introduction
============

Elastic is based on the standard elasticity theory (see [LL]_ for the detailed
introduction) and *finite deformation* approach to the calculation of elastic 
tensor of the crystal.

Physical Principles
-------------------

I have described basic physical principles on which the code rests in my 
habilitation thesis. Here I will include slightly edited second chapter of the 
thesis introducing the method and some implementation details.

Elasticity of crystals
^^^^^^^^^^^^^^^^^^^^^^

The classical, linear theory of elasticity of crystalline materials has been formu-
lated already in the 18th and 19th century by Cauchy, Euler, Poisson, Young and
many other great mathematicians and physicists of that time. The standard text-
book formulation (e.g. classical book by Landau et al. [LL]_) can be, in prin-
ciple, directly used as a basis for numerical determination of the elastic tensor
and other mechanical properties of the crystal. Nevertheless, practical imple-
mentation of these formulas have some non-obvious aspects, worthy of explicit
presentation.
The *finite deformation* method developed and used in the presented papers [1,
2, 3, 4, 5, 6, 7] is based on the fundamental relationship between stress and strain
of the solid crystalline body with a particular symmetry. This is a simple tensor
equation, sometimes called generalised *Hook’s law* (in standard tensor notation):

.. math::
    \sigma_{\lambda\ksi} = C_{λξμν} s_{μν}

This formula simply states that the stress in the crystal :math:`\sigma_{\lambda\ksi}` is a linear function
of the strain sμν incurred by its deformation, and the elasticity tensor Cλξμν is
just a tensor proportionality coefficient. The Greek indexes run through coordin-
ates x, y, z. The elasticity tensor inherits symmetries of the crystal and has some
intrinsic symmetries of its own. Therefore, only a small number of its compon-
ents are independent. This fact leads to customary representation of this entity
in the form of the matrix with components assigned according to Voight’s nota-
tion. Thus, instead of the rank-4 three dimensional tensor we have 6 × 6 matrix
Cij where the indexes i, j = 1 . . . 6. The stress and strain tensors are represen-
ted as six-dimensional vectors. The symmetries of the elastic tensor are directly
translated into symmetries of the Cij matrix. The Voight’s notation is commonly
used in tensor calculus. For this particular case we can write it as an index as-
signment where each pair of Greek indexes is replaced with a corresponding
Latin index (i, j, k, l, m, n): xx = 1, yy = 2, zz = 3, yz = 4, xz = 5, xy = 6. While
this convention makes presentation of elastic constants much easier – since it is
just a square table of numbers – it slightly complicates algebraic procedures as
we lose the simplicity of the tensor formalism.
Every class of crystal cells implies, through its symmetry, a different num-
ber of independent parameters defining the Cij matrix. For example, the cubic
lattice has just three independent elastic constants: C11 , C12 , C44 and the matrix
looks as follows:


Less symmetric crystals have, naturally, a higher number of independent elastic
constants and lower symmetry of the Cij matrix (see [53] for full introduction to
theory of elasticity).



Implementation
--------------

Elastic is implemented as an extension module to ASE system



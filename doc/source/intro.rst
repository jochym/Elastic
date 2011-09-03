Physical Principles
===================

Elastic is based on the standard elasticity theory (see [LL]_ for the detailed
introduction) and *finite deformation* approach to the calculation of elastic 
tensor of the crystal.
I have described basic physical principles on which the code rests in my 
habilitation thesis. Here I will include slightly edited second chapter of the 
thesis introducing the method and some implementation details.

Elasticity of crystals
----------------------

The classical, linear theory of elasticity of crystalline materials has been
formulated already in the 18th and 19th century by Cauchy, Euler, Poisson,
Young and many other great mathematicians and physicists of that time. The
standard textbook formulation (e.g. classical book by Landau et al. [LL]_) can
be, in principle, directly used as a basis for numerical determination of the
elastic tensor and other mechanical properties of the crystal. Nevertheless,
practical implementation of these formulas have some non-obvious aspects,
worthy of explicit presentation. The *finite deformation* method developed and
used in the mentioned papers [TiC]_, [ZrC]_ is based on the fundamental
relationship between stress and strain of the solid crystalline body with a
particular symmetry. This is a simple tensor equation, sometimes called
generalised *Hook’s law* (in standard tensor notation):

.. math::
    \sigma_{\lambda\xi} = C_{\lambda\xi\mu\nu} s_{\mu\nu}

This formula simply states that the stress in the crystal
:math:`\sigma_{\lambda\xi}` is a linear function of the strain
:math:`s_{\mu\nu}` incurred by its deformation, and the elasticity tensor
:math:`C_{\lambda\xi\mu\nu}` is just a tensor proportionality coefficient. The
Greek indexes run through coordinates x, y, z. The elasticity tensor inherits
symmetries of the crystal and has some intrinsic symmetries of its own.
Therefore, only a small number of its components are independent. This fact
leads to customary representation of this entity in the form of the matrix with
components assigned according to Voight’s notation. Thus, instead of the rank-4
three dimensional tensor we have :math:`6 \times 6` matrix :math:`C_{ij}` where
the indexes :math:`i, j = 1 \ldots 6`. The stress and strain tensors are
represented as six-dimensional vectors. The symmetries of the elastic tensor are
directly translated into symmetries of the :math:`C_{ij}` matrix. The Voight’s
notation is commonly used in tensor calculus. For this particular case we can
write it as an index assignment where each pair of Greek indexes is replaced
with a corresponding Latin index (i, j, k, l, m, n): xx=1, yy=2, zz=3, 
yz=4, xz=5, xy=6.

While this convention makes presentation of elastic constants much easier -
since it is just a square table of numbers - it slightly complicates algebraic
procedures as we lose the simplicity of the tensor formalism. Every class of
crystal implies, through its symmetry, a different number of independent
elements in the :math:`C_{ij}` matrix. 

For example, the cubic lattice has just three independent elements in the 
elastic matrix: :math:`C_{11}, C_{12}, C_{44}`, and the matrix itself has the following shape:

.. math::
    \left[\begin{array}{cccccc}
    C_{11} & C_{12} & C_{12} & 0 & 0 & 0\\
    C_{12} & C_{11} & C_{12} & 0 & 0 & 0\\
    C_{12} & C_{12} & C_{11} & 0 & 0 & 0\\
    0 & 0 & 0 & C_{44} & 0 & 0\\
    0 & 0 & 0 & 0 & C_{44} & 0\\
    0 & 0 & 0 & 0 & 0 & C_{44}\end{array}\right]

Less symmetric crystals have, naturally, a higher number of independent elastic
constants and lower symmetry of the :math:`C_{ij}` matrix (see [LL]_ for full
introduction to theory of elasticity).

Numerical derivation of elastic matrix
--------------------------------------

Numerical derivation of the :math:`C_{ij}` matrix may be approached in many
different ways. Basically, we can employ the same methods as used
effectively in experimental work. From all experimental procedures we can select
three classes which are relevant to our discussion:

#. Based on the measured sound velocity, including various methods based on determination of lattice dynamics of the crystal.
#. Based on the strain-energy relation.
#. Based on the measured stress-strain relations for some particular, simple strains.

While the first method is frequently used in laboratory measurements, it is not
direct and is not well suited to numerical derivation. For example, you can
measure the tangent of all acoustic branches of phonon dispersion curves in
several directions to get enough data points to solve the set of equations
for most of the independent components of the :math:`C_{ij}` matrix. The
tangent of the acoustic branch is connected with the sound velocity and with
components of elastic matrix by a set of equations of the general form:

.. math::
    \varrho v_{k}^{2}=L(C_{ij})

where :math:`L(C_{ij})` is a linear combination of independent components
of elastic tensor, :math:`v_{k}` is a long-wave sound velocity in particular
direction, which is equivalent to the slope of the acoustic branch
of phonon dispersion curve in this direction, and :math:`\varrho` is crystal
density. Full set of these equations for the cubic crystal is included
in [TiC]_. Unfortunately, it is difficult and non-practical
to use this method to obtain more then few of the simplest of components,
since the numerical properties of the non-linear formulas involved
lead to the error pile-up in the results. It is particularly susceptible
to errors in long-wave sound velocities -- due to the quadratic function
in above equation. Unfortunately, these asymptotic velocities
are particularly weakly constrained by most of available computational
methods. The same formulas can also be used to obtain elastic matrix
from straight-forward sound velocity measurements. The same unfavourable
numerical properties lead to high demands on accuracy of the measurements
-- but in this case these requirements could be quite easily met in
experiment since sound velocity can be measured with very high precision. 

The second method is not practical for laboratory measurements - it is not easy
to accurately measure energy of the deformed crystal. Furthermore, the
strain-energy relation is non-linear and we need to extract a derivative of the
function -- the procedure is quite complex, needs more data points and is prone
to errors.

The third method is well suited for experimental work as well as computational
derivation of the elastic matrix. The numerical properties of the formulas --
being just a set of linear equations -- are well known and provide stable and
well-controlled error propagation. Furthermore, while the sound velocity is not
directly accessible to computational quantum mechanical methods, the stresses
induced by strains on the crystal are almost universally provided by DFT based
programs and often do not require any additional computational effort. The
comparison of these methods used for computational derivation of the elastic
matrix is included in [TiC]_, [ZrC]_. The comparison shows that the finite
deformation (stress-strain) method compares favourably to the pure
energy-derivative method. The results clearly show that the strain--stress
relationship approach described here is much better suited for computational
derivation of elastic matrix and provides lower error level than other two
methods.

.. _symmetry:

Crystal symmetry and elastic matrix derivation
----------------------------------------------

As mentioned above, the symmetry of the crystal determines the number and
position of independent components of the :math:`C_{ij}` matrix. Therefore, the
stress-strain relation is effectively modified by the symmetry of the case by a
simple fact that most, of the coefficients are not independent from one another.
We aim to derive the complete set of :math:`C_{ij}` elements from the set of
computational or experimental measurements of strain and stress tensors
:math:`s^{a}`, :math:`\sigma^{a}` where the upper Latin index a numbers a
calculation/experiment setup. In the case described here the “measurement” is a
particular computational setup with the crystal deformed in various ways in
order to provide enough data points to derive all independent components of the
:math:`C_{ij}` matrix. The set of necessary deformations can be determined by
the symmetry of the crystal and contains tetragonal and sheer deformations along
some or all axis -- as the symmetry of the case dictates. To improve the
accuracy of the results the deformations may be of different sizes (typically
0.1-1% in length or 0.1-1 degree in angle). 

Having a set of calculation data :math:`\{s^{a}, \sigma^{a}\}`, we can rewrite
generalised Hook's law to form a set of linear equations (in Voight notation for
:math:`i,j` indexes): :math:`C_{ij}s_{j}^{a}=\sigma_{i}^{a}`. This set can be
further transformed for each symmetry case to the form in which the independent
components of the :math:`C_{ij}` matrix create a vector of unknowns and the
symmetry relations and strains :math:`s_{j}^{a}` create a new equation matrix
:math:`S`. :math:`S_{ju}(s^{a})C_{u}=\sigma_{j}^{a}`. The :math:`S(s)` matrix is
a linear function of the strain vector s with all symmetry relations taken into
account. The index a runs over all data sets we have in the calculation while
index u runs over all independent components of the :math:`C_{ij}` matrix. For
the cubic crystal the above equation takes explicit form:

.. math::
    \left[\begin{array}{ccc}
    s_{1} & s_{2}+s_{3} & 0\\
    s_{2} & s_{1}+s_{3} & 0\\
    s_{3} & s_{1}+s_{2} & 0\\
    0 & 0 & 2s_{4}\\
    0 & 0 & 2s_{5}\\
    0 & 0 & 2s_{6}\end{array}\right]^{a}\left[\begin{array}{c}
    C_{11}\\
    C_{12}\\
    C_{44}\end{array}\right]=\left[\begin{array}{c}
    \sigma_{1}\\
    \sigma_{2}\\
    \sigma_{3}\\
    \sigma_{4}\\
    \sigma_{5}\\
    \sigma_{6}\end{array}\right]^{a}.

Note the a index of S and :math:`\sigma`, which creates a set of
:math:`n\times6` linear equations for 3 unknowns
:math:`\left[C_{11},C_{12},C_{44}\right]`, where n is a number of independent
calculations of stresses incurred in crystal by strains. In principle, the above
relations could be expressed in the non-symmetry specific form with either a
full set of indexes and the symmetry information encoded in the single matrix of
constant elements or even in the pure tensor formulation with the four-index
elastic tensor :math:`C` and two-index stress and strain tensors. While this
type of formulation is definitely more regular and sometimes easier to
manipulate in formal transformations, it is not very useful for numerical
calculations or writing computer code -- multi-dimensional arrays are difficult
to manipulate and are prone to many trivial notation errors. Thus, it is better
to split the general formula to crystal classes with different number of
:math:`C_{ij}` components (i.e. length of the :math:`C_{u}` vector)
and separate shape of the :math:`S` matrix. This is an approach used by Elastic. 

For example, in the orthorhombic crystal the vector of independent
:math:`C_{ij}` components has nine elements and the S matrix is a :math:`9\times6` one:

.. math::
    \left[\begin{array}{ccccccccc}
    s_{1} & 0 & 0 & s_{2} & s_{3} & 0 & 0 & 0 & 0\\
    0 & s_{2} & 0 & s_{1} & 0 & s_{3} & 0 & 0 & 0\\
    0 & 0 & s_{3} & 0 & s_{1} & s_{2} & 0 & 0 & 0\\
    0 & 0 & 0 & 0 & 0 & 0 & 2s_{4} & 0 & 0\\
    0 & 0 & 0 & 0 & 0 & 0 & 0 & 2s_{5} & 0\\
    0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 2s_{6}\end{array}\right]^{a}\left[\begin{array}{c}
    C_{11}\\
    C_{22}\\
    C_{33}\\
    C_{12}\\
    C_{13}\\
    C_{23}\\
    C_{44}\\
    C_{55}\\
    C_{66}\end{array}\right]=\left[\begin{array}{c}
    \sigma_{1}\\
    \sigma_{2}\\
    \sigma_{3}\\
    \sigma_{4}\\
    \sigma_{5}\\
    \sigma_{6}\end{array}\right]^{a}.

The elements of the matrix S have direct relation to the terms of expansion of
the elastic free energy as a function of deformation (strain tensor) F(s). For
example, the orthorhombic equation can be derived from the free energy formula
(see [LL]_ for derivation):

.. math::
    F(s)   =  \frac{1}{2}C_{11}s_{1}^{2}+
             \frac{1}{2}C_{22}s_{2}^{2}+
             \frac{1}{2}C_{33}s_{3}^{2}+ 
     C_{12}s_{1}s_{2}+C_{13}s_{1}s_{3}+C_{23}s_{2}s_{3}+  \\
     2C_{44}s_{4}^{2}+2C_{55}s_{5}^{2}+2C_{66}s_{6}^{2}

The elements of the S matrix are simply coefficients of first derivatives of the
F(s) over respective strain components. Alternatively, we can rewrite the S(s)
matrix in the compact form as a mixed derivative: 

.. math::
    S_{iu}=A\frac{\partial^{2}F}{\partial s_{i}\partial C_{u}},

where A is a multiplier taking into account the double counting of the
off-diagonal components in the free energy formula (see note at the end of the
exercises in [LL]_). The multiplier :math:`A=1` for
:math:`i \leq 4`, and :math:`1/2` otherwise. The above general formula turns out to be quite
helpful in less trivial cases of trigonal or hexagonal classes. For instance,
the hexagonal elastic free energy (see [LL]_ for rather lengthy formula) leads
to the following set of equations:

.. math::
    \left[\begin{array}{ccccc}
    s_{1} & 0 & s_{2} & s_{3} & 0\\
    s_{2} & 0 & s_{1} & s_{3} & 0\\
    0 & s_{3} & 0 & s_{1}+s_{2} & 0\\
    0 & 0 & 0 & 0 & 2s_{4}\\
    0 & 0 & 0 & 0 & 2s_{5}\\
    s_{6} & 0 & -s_{6} & 0 & 0\end{array}\right]^{a}\left[\begin{array}{c}
    C_{11}\\
    C_{33}\\
    C_{12}\\
    C_{13}\\
    C_{44}\end{array}\right]=\left[\begin{array}{c}
    \sigma_{1}\\
    \sigma_{2}\\
    \sigma_{3}\\
    \sigma_{4}\\
    \sigma_{5}\\
    \sigma_{6}\end{array}\right]^{a}.

The set of linear equations, with calculated strains and stresses inserted
into the :math:`S^{a}` matrix and :math:`\sigma^{a}` vector, could be
constructed for any crystal -- only the form of the S matrix and the length of
the :math:`C_{u}` vector will be different for each symmetry. 

The set of equations is usually over-determined. Therefore, it
cannot be solved in the strict linear-algebra sense since no exact solution
could exist. Nevertheless, this set of equations can be solved in approximate
sense -- i.e. minimising the length of the residual vector of the solution.
Fortunately, a very clever algorithm capable of dealing with just this type of
linear equations has been known for a long time. It is called Singular Value
Decomposition [SVD]_. Not only does it provide the approximate solution
minimising the residual vector of the equation but also is stable against
numerically ill-conditioned equations or equations which provide too little data
to determine all components of the solution. The SVD provides also some
indication of the quality of the obtained solution in the form of the vector of
singular values, which could be used to judge whether the solution is
well-determined. It is a well known algorithm and its implementations are
available in every self-respecting numerical linear algebra library. The
implementation used in the Elastic code is the one included in the Scientific
Python library `SciPy <http://www.scipy.org/>`_.




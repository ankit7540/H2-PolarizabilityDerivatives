# H2-PolarizabilityDerivatives : under progress
Set of data on the static and dynamic invariants of polarizability (mean polarizability and anisotropy) with programs for obtaining the derivatives of these invariants (for H<sub>2</sub>, HD and D<sub>2</sub>) for specific inter-nuclear distance defined by rovibrational wavefunction for some state.

This repository contains :
 - Internuclear distance dependent polarizability from which the two non-zero invariants are accessible for the analysis of the derivatives of and Taylor series expansion. These invariants are

Property | Definition
------------ | -------------
Mean polarizability | <img src="https://github.com/ankit7540/H2-PolarizabilityMatrixElements/blob/master/image/alpha_mp.png" data-canonical-src="https://github.com/ankit7540/H2-PolarizabilityMatrixElements/blob/master/image/alpha_mp.png" width="195" height="28" />
Polarizability anisotropy | <img src="https://github.com/ankit7540/H2-PolarizabilityMatrixElements/blob/master/image/gamma.png" data-canonical-src="https://github.com/ankit7540/H2-PolarizabilityMatrixElements/blob/master/image/gamma.png" width="155" height="28" />

The above properties are available as parameter for H<sub>2</sub> HD and D<sub>2</sub>.
 - Rovibrational wavefunctions for H<sub>2</sub>, HD and D<sub>2</sub> for v=0--2 and J=0--15. (*This repository deals with computation of the derivatives of the invariants and their Taylor series expansions. For seamlessly performing this set of computation the inter-nuclear distance is truncated to 0.5--3.0 a.u. which allows for the inclusion of wavefunctions up to v=2, J=15. )
 - A python module is included which performs the above computation for the wavelength range : 182.25 to 1320.6 nm.

**Available programs**
---
The programs for computation of matrix element (which includes cubic spline interpolation and numerical integration) are written in the Python program `pol_derivative.py`.

**Usage**
---
Clone this repository or download as a zip file. According to the program of choice, refer to the `README.md` for the repository and the `README.md`in the Python-module folder.


**Comments on numerical accuracy**
---
The definite integral calculation is usually accurate to ~1e-6 or better. However, the net numerical uncertainity in the computed matrix element is  +/- 1e-4 which includes the uncertainities introduced by the accuracy of the wavefunctions, polarizability, spline interpolation procedures and physical constants.

**Comments on the sign of the matrix element**
---
Some matrix elements computed may have -ve sign which arises due to the phase of the wavefunction. In most applications, the square of the matrix elements are needed and thus the sign maybe of no real consequence.

**Computational details
---
![integral image][img0]

where, <img src="https://github.com/ankit7540/H2-PolarizabilityMatrixElements/blob/master/image/rmin.png" data-canonical-src="https://github.com/ankit7540/H2-PolarizabilityMatrixElements/blob/master/image/rmin.png" width="45" height="15" /> = 0.2 *a.u.* and  <img src="https://github.com/ankit7540/H2-PolarizabilityMatrixElements/blob/master/image/rmax.png" data-canonical-src="https://github.com/ankit7540/H2-PolarizabilityMatrixElements/blob/master/image/rmax.png" width="45" height="15" /> = 4.48 *a.u.*

**Credits**
---
Cubic spline interpolation procedure used in FORTRAN and python codes has been adapted from Numerical Recipes in FORTRAN, William H. Press, Saul A. Teukolsky, William T. Vetterling, Brian P. Flannery, Michael Metcalf, Cambridge University Press; 2<sup>nd</sup> edition.

For evaluation of the definite integral the Adaptive Gausssian Quadrature implemented in SciPy has been used.

**References on the evaluation of the definite integral and implementation:**
- T. N. L. Patterson, Math. Comput. 22, 847 (1968)
- T. N. L. Patterson, Math. Comput. 23, 892 (1969)
- R. Piessens, E. de Doncker-Kapenga, C. Uberhuber, and D. Kahaner, Quadpack - A Sub-routine Package for Automatic Integration (Springer-Verlag Berlin Heidelberg, 1983)


Python code by Ankit Raj (NCTU, Taiwan).

---




[img0]: https://github.com/ankit7540/H2-PolarizabilityMatrixElements/blob/master/image/01-05-2018_82.png "Logo Title Text 2"
[img1]: https://github.com/ankit7540/H2-PolarizabilityMatrixElements/blob/master/image/alpha_perp.png "Logo alpha_{perp}"
[img2]: https://github.com/ankit7540/H2-PolarizabilityMatrixElements/blob/master/image/alpha_parallel.png "Logo alpha_{paralell}"
[img3]: https://github.com/ankit7540/H2-PolarizabilityMatrixElements/blob/master/image/alpha_mp.png "Logo alpha_{mp}"
[img4]: https://github.com/ankit7540/H2-PolarizabilityMatrixElements/blob/master/image/gamma.png "Logo alpha_{aniso}"


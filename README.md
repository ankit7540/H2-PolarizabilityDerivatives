# H<sub>2</sub>-PolarizabilityDerivatives
Link to the article : https://doi.org/10.1080/00268976.2019.1632950 | DOI  for code : [![DOI](https://zenodo.org/badge/166674944.svg)](https://zenodo.org/badge/latestdoi/166674944)



Set of data on the static and dynamic invariants of polarizability (mean polarizability and anisotropy) with programs for obtaining the derivatives of these invariants (for H<sub>2</sub>, HD and D<sub>2</sub>) for specific inter-nuclear distance defined by rovibrational wavefunction for some state.

This repository contains :
 - Internuclear distance dependent polarizability from which the two non-zero invariants are accessible for the analysis of the derivatives of and Taylor series expansion. These invariants are

Property | Definition
------------ | -------------
Mean polarizability | <img src="https://github.com/ankit7540/H2-PolarizabilityMatrixElements/blob/master/image/alpha_mp.png" data-canonical-src="https://github.com/ankit7540/H2-PolarizabilityMatrixElements/blob/master/image/alpha_mp.png" width="195" height="28" />
Polarizability anisotropy | <img src="https://github.com/ankit7540/H2-PolarizabilityMatrixElements/blob/master/image/gamma.png" data-canonical-src="https://github.com/ankit7540/H2-PolarizabilityMatrixElements/blob/master/image/gamma.png" width="155" height="28" />

The above properties are available as parameter for H<sub>2</sub> HD and D<sub>2</sub>.
 - Rovibrational wavefunctions for H<sub>2</sub>, HD and D<sub>2</sub> for v=0--2 and J=0--15. (**#** This repository deals with computation of the derivatives of the invariants and their Taylor series expansions. For seamlessly performing this set of computation with minimal introduction of error, the inter-nuclear distance is truncated to 0.5--3.0 a.u. for v=0-2 and 0.7--3.8 for v=3-4 a.u. )
 - A python module is included which performs the above computation for the wavelength range : 182.26 to 1320.6 nm.

**Available programs**
---
The programs for computation of matrix element (which includes cubic spline interpolation and numerical integration) are written in the Python program `pol_derivative.py`.

**Usage**
---
Clone this repository or download as a zip file. According to the program of choice, refer to the `README.md` for the repository and the `README.md`in the Python-module folder.


Computational details
---
1. Parameter (mean polarizability or anisotropy) is interpolated to the asked wavelength. This is followed by fitting it over inter-nuclear distance using a polynomial function. Fit coefficients are reported in the output.
2. The expectation value (r<sub>e</sub>) of the inter-nuclear distance is computed for specific rovibrational state of H<sub>2</sub>, HD and D<sub>2</sub> , in the electronic ground state.
3. Using the coefficients of polynomial function the derivatives of the parameter are computed at r<sub>e</sub>. Set of derivatives up to 7<sup>th</sup>  order are reported in the output.
4. Taylor series expansions are generated using the derivatives centered at r<sub>e</sub>. Subsequently, the matrix elements of the Taylor series expansions of the parameter are computed along with that of the original parameter. The set of these matrix elements are reported in the output.

- The matrix elements are computed as the following integral over the inter-nuclear distance. 
![integral image][img0]

where, <img src="https://github.com/ankit7540/H2-PolarizabilityMatrixElements/blob/master/image/rmin.png" data-canonical-src="https://github.com/ankit7540/H2-PolarizabilityMatrixElements/blob/master/image/rmin.png" width="45" height="15" /> = 0.5 *a.u.* and  <img src="https://github.com/ankit7540/H2-PolarizabilityMatrixElements/blob/master/image/rmax.png" data-canonical-src="https://github.com/ankit7540/H2-PolarizabilityMatrixElements/blob/master/image/rmax.png" width="45" height="15" /> = 3.0 *a.u.*

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
**This work has been published in the following article:**

Ankit Raj, Henryk A. Witek & Hiro-o Hamaguchi<br>
Vibration–rotation interactions in H2, HD and D2 : centrifugal distortion factors and the derivatives of polarisability invariants<br>
*Molecular Physics*<br>
DOI: https://doi.org/10.1080/00268976.2019.1632950

Citing this code : [![DOI](https://zenodo.org/badge/166674944.svg)](https://zenodo.org/badge/latestdoi/166674944)


[img0]: https://github.com/ankit7540/H2-PolarizabilityMatrixElements/blob/master/image/01-05-2018_82.png "Logo Title Text 2"
[img1]: https://github.com/ankit7540/H2-PolarizabilityMatrixElements/blob/master/image/alpha_perp.png "Logo alpha_{perp}"
[img2]: https://github.com/ankit7540/H2-PolarizabilityMatrixElements/blob/master/image/alpha_parallel.png "Logo alpha_{paralell}"
[img3]: https://github.com/ankit7540/H2-PolarizabilityMatrixElements/blob/master/image/alpha_mp.png "Logo alpha_{mp}"
[img4]: https://github.com/ankit7540/H2-PolarizabilityMatrixElements/blob/master/image/gamma.png "Logo alpha_{aniso}"


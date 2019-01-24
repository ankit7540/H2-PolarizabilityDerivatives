Python-module `pol_derivative`
----------------
`pol_derivative` is python module which performs computation for obtaining the derivatives of mean polarizability and polarizability anisotropy invariants (for H2, HD and D2) for specific inter-nuclear distance defined by rovibrational wavefunction for some state.


Requirements
----------------
Python 2.x or Python 3.x with `numpy` and `scipy` modules. For plot of fit, residual and the Taylor series expansions `matplotlib` module is required.

Usage
----------------
1. After cloning the repository and moving in the `python-module` directory, add the current folder to path allowing to import the module in the current folder (this is not needed with python3). The following commands are run under the Python interpreter environment.
    > import sys
    
    > sys.path.append("..")
     
2. Import the `pol_derivative` which should be in your current folder. Directly execute the following command when using Python3.
    > import pol_derivative
3. If all requirements are met the following output should be produced.
 ```
Give  pol_derivative.compute  command with parameters:
        pol_derivative.compute(molecule, v, J, lambda, unit of lambda, operator)
         for example:  pol_derivative.compute("H2",0,4,488,"n","mp")  
                       pol_derivative.compute("D2",1,0,"static","n","g")  
                
                molecule = for H2 enter "H2", for D2 enter "D2", for HD enter "HD" 
                v    = vibrational state, [0,2]
                J    = rotataional state, [0,15]
                lambda   = wavelength in Hartree, nm or Angstrom, for static specify "s" or "static" here
                unit of lambda =  for  Hartree           use "H" or "h"  
                                  for  nanometers        use "n" or "nm" 
                                  for  Angstrom          use "a" or "A"  
                                  if static property is asked then this parameter can be any of the three 
                Available wavelength range: 0.25 - 0.0345 Hartree;
                                            182.2534 - 1320.6769 nm; 
                                            1822.5341 - 13206.7688 Angstrom
                operator        = isotropy or mean polarizability given by "iso" or "mp" or "mean" 
                                  anisotropy or polarizability difference or gamma given by "aniso" or "g"  or "diff" 
                enable_plot     = 0 or 1 
                                  0 = do not plot  
                                  1 = enable plot, requires matplotlib module  
                                      The plots will be generated as an output.pdf file having 3 plots.
                                      First plot = fit of the parameter with polynomial
                                      Second plot = residual of the fit
                                      Third plot = Taylor series expansion of the parameter along with the 
                                                  original parameter function
...ready.
 ```
4. Use the following command to do computation of the matrix element.
    > pol_derivative.compute(molecule, v, J, lambda, unit of lambda, operator)
        
    where the parameters are described below: 
      
    - mol  =    molecule specification (for H<sub>2</sub> enter "H2", for D<sub>2</sub> enter "D2", for HD enter "HD")
    - v   =    vibrational state for the bra, vl = [0,2]
    - J   =    rotational state for the bra,  Jl = [0,15]
    - wavelength =  wavelength within the specified range ( 0.25 - 0.0345 Hartree;  182.2534 - 1320.6768  nm;  1822.5341 - 13206.7688  Angstrom ). Specify unit accordingly in the next parameter. If static polarizability is needed enter "static" or "s" here.
    - wavelength_unit = specify unit using the specifier, ( for  Hartree use "H" or "h" , for  nanometers use "n" or "nm" , for  Angstrom use "a" or "A"  )
    - operator   = property namely mean polarizability (isotropy) and anisotropy. Specify this  using the specifier. ( For  isotropy  use "iso"   or  "mp" or "mean" , for  anisotropy use "aniso" or  "g"  or "diff".
    

**Examples**
---

Outputs obtained from the command execution is shown below.

- Analysis of polarizability anisotropy (400 nm) for H<sub>2</sub> at v=0, J=0.
 
```
 pol_derivative.compute("H2", 0, 0, 400, "nm", "g", 0) ⏎
Selected wavelength in nanometer : 400.0, Hartree : 0.113908
-------------------------------------------------------------------
Analysis for the derivative of the parameter at re defined by v,J
-------------------------------------------------------------------
Molecule = H2
Parameter = ['anisotropy']
Wavelength = 400 nm
Rovibrational state : v=0, J=0
Expectation value of inter-nuclear distance, r_e = 1.448725 a.u.

(1.) Fit coefficients (scaled back using 3^{n})
[parameter vs distance was fit using poly11_sc function]
c0 = -0.0089440422
c1 = 0.0721460808
c2 = 0.1607448133
c3 = 0.5929436216
c4 = -0.0088069932
c5 = -0.1681265555
c6 = 0.3340732142
c7 = -0.2943666554
c8 = 0.1450605756
c9 = -0.0430924111
c10 = 0.0070572395
c11 = -0.0004805218

(2.) Derivatives of parameter at r_e
g0 = 2.13246483
g1 = 3.98882905
g2 = 4.6636114
g3 = 0.7455217
g4 = -9.09518003
g5 = -19.10712744
g6 = -15.181164
g7 = 37.29035527

(3.) Matrix elements of parameter using Taylor series
  expansions (at r_e) using n^{th} order derivatives
<psi_0,0| g(0) |psi_0,0> = 2.132465
<psi_0,0| g(1) |psi_0,0> = 2.132463
<psi_0,0| g(2) |psi_0,0> = 2.198224
<psi_0,0| g(3) |psi_0,0> = 2.198336
<psi_0,0| g(4) |psi_0,0> = 2.197414
<psi_0,0| g(5) |psi_0,0> = 2.197374
<psi_0,0| g(6) |psi_0,0> = 2.197366
<psi_0,0| g(7) |psi_0,0> = 2.197366
<psi_0,0| g(infty) |psi_0,0> = 2.19736
-------------------------------------------------------------------
``` 


- Analysis of mean polarizability (300 nm) for D<sub>2</sub> at v=2, J=11.

```
pol_derivative.compute("D2", 2, 11, 300, "nm", "mp", 0) ⏎
Selected wavelength in nanometer : 300.0, Hartree : 0.151878
-------------------------------------------------------------------
Analysis for the derivative of the parameter at re defined by v,J
-------------------------------------------------------------------
Molecule = D2
Parameter = ['isotropy']
Wavelength = 300 nm
Rovibrational state : v=2, J=11
Expectation value of inter-nuclear distance, r_e = 1.64561 a.u.

(1.) Fit coefficients (scaled back using 3^{n})
[parameter vs distance was fit using poly11_sc function]
c0 = 1.4007005671
c1 = 0.2808254652
c2 = 3.8774527695
c3 = -4.3600502244
c4 = 5.2389317767
c5 = -4.4707273893
c6 = 2.7755601875
c7 = -1.2448614847
c8 = 0.3907130685
c9 = -0.0816625614
c10 = 0.0101647148
c11 = -0.0005639053

(2.) Derivatives of parameter at r_e
g0 = 6.96793408
g1 = 5.71318835
g2 = 2.23045816
g3 = -2.43978379
g4 = -6.26225802
g5 = -8.79800739
g6 = 1.58604615
g7 = 43.67467196

(3.) Matrix elements of parameter using Taylor series
  expansions (at r_e) using n^{th} order derivatives
<psi_2,11| mp(0) |psi_2,11> = 6.967934
<psi_2,11| mp(1) |psi_2,11> = 6.967932
<psi_2,11| mp(2) |psi_2,11> = 7.086703
<psi_2,11| mp(3) |psi_2,11> = 7.090093
<psi_2,11| mp(4) |psi_2,11> = 7.08515
<psi_2,11| mp(5) |psi_2,11> = 7.08535
<psi_2,11| mp(6) |psi_2,11> = 7.08536
<psi_2,11| mp(7) |psi_2,11> = 7.085353
<psi_2,11| mp(infty) |psi_2,11> = 7.085359
-------------------------------------------------------------------

``` 


 
 
[f1]: http://chart.apis.google.com/chart?cht=tx&chl=\langle\psi_{v=0,J=0}|\bar{\alpha}|\psi_{v=0,J=0}\rangle
[f2]: http://chart.apis.google.com/chart?cht=tx&chl=\langle\psi_{v=2,J=1}|\gamma|\psi_{v=1,J=1}\rangle
[f3]: http://chart.apis.google.com/chart?cht=tx&chl=\langle\psi_{v=2,J=1}|\alpha_{\parallel}|\psi_{v=1,J=1}\rangle


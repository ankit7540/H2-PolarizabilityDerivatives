Python-module `pol_derivative`
----------------
`pol_derivative` is python module which performs computation for obtaining the derivatives of mean polarizability and polarizability anisotropy invariants (for H<sub>2</sub>, HD and D<sub>2</sub>) for specific inter-nuclear distance defined by rovibrational state.


Requirements
----------------
Python 2.7 or Python 3.x with `numpy` and `scipy` modules. For generating the plot of fit, residual and the Taylor series expansions `matplotlib` module is required.

Usage
----------------
The following commands are run under the Python interpreter environment.

1. After cloning the repository and moving in the `python-module` directory, add the current folder to path allowing to import the module in the current folder (this is required when using Python 2.7). 
    > import sys
    
    > sys.path.append("..")
     
2. Import the `pol_derivative` which should be in your current folder (directly execute the following command when using Python3).
    > import pol_derivative
3. If all requirements are met the following output should be produced.
 ```
Give  pol_derivative.compute  command with parameters:
        pol_derivative.compute(molecule, v, J, lambda, unit of lambda, operator, enable_plot)
         for example:  pol_derivative.compute("H2", 0, 4, 488, "n", "mp", 0)  
                       pol_derivative.compute("D2", 1, 0, "static", "n", "g", 1)  
                
                molecule = for H2 enter "H2", for D2 enter "D2", for HD enter "HD" 
                v    = vibrational state, [0,4]
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
    > pol_derivative.compute(molecule, v, J, lambda, unit of lambda, operator, enable_plot)
        
    where the parameters are described below: 
      
    - mol  =    molecule specification (for H<sub>2</sub> enter "H2", for D<sub>2</sub> enter "D2", for HD enter "HD")
    - v   =    vibrational state for the bra, vl = [0,4]
    - J   =    rotational state for the bra,  Jl = [0,15]
    - wavelength =  wavelength within the specified range ( 0.25 - 0.0345 Hartree;  182.2534 - 1320.6768  nm;  1822.5341 - 13206.7688  Angstrom ). Specify unit accordingly in the next parameter. If static polarizability is needed enter "static" or "s" here.
    - wavelength_unit = specify unit using the specifier, ( for  Hartree use "H" or "h" , for  nanometers use "n" or "nm" , for  Angstrom use "a" or "A"  )
    - operator   = property namely mean polarizability (isotropy) and anisotropy. Specify this  using the specifier. ( For  isotropy  use "iso"   or  "mp" or "mean" , for  anisotropy use "aniso" or  "g"  or "diff".
    - enable_plot = boolean, 0 or 1 control the plotting of fit, residual of fit, and the Taylor series expansions of the parameter in a pdf file.
    

**Examples**
---

Outputs obtained from the command execution is shown below. The derivatives at r<sub>e</sub> are shown as d0,..,d11. The matrix elements are listed as g(0),.., g(7) or mp(0),... , mp(7).

- Analysis of polarizability anisotropy (400 nm) for H<sub>2</sub> at v=0, J=0.
 
```
 pol_derivative.compute("H2", 0, 0, 325, "nm", "g", 0) ⏎
Selected wavelength in nanometer : 325.0, Hartree : 0.140195
-------------------------------------------------------------------
Analysis for the derivative of the parameter at re defined by v,J
-------------------------------------------------------------------
Molecule = H2
Parameter = ['anisotropy']
Wavelength = 325 nm
Rovibrational state : v=0, J=0
Expectation value of inter-nuclear distance, r_e = 1.448725 a.u.

(1.) Derivatives of parameter at r_e
d0 = 2.22004074
d1 = 4.23570551
d2 = 5.22147446
d3 = 1.6413199
d4 = -8.60629586
d5 = -21.78515885
d6 = -26.09275114
d7 = 13.60065127
d8 = 324.9030138
d9 = 1897.2137504
d10 = 673.76189162
d11 = -32564.25330537

(2.) Matrix elements of parameter using Taylor series 
  expansions (at r_e) using n^{th} order derivatives
<psi_0,0| g(0) |psi_0,0> = 2.22004
<psi_0,0| g(1) |psi_0,0> = 2.22004
<psi_0,0| g(2) |psi_0,0> = 2.29367
<psi_0,0| g(3) |psi_0,0> = 2.29391
<psi_0,0| g(4) |psi_0,0> = 2.29304
<psi_0,0| g(5) |psi_0,0> = 2.29299
<psi_0,0| g(6) |psi_0,0> = 2.29298
<psi_0,0| g(7) |psi_0,0> = 2.29298
<psi_0,0| g(infty) |psi_0,0> = 2.29298
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

(1.) Derivatives of parameter at r_e
d0 = 6.96793144
d1 = 5.71319624
d2 = 2.23059187
d3 = -2.4404971
d4 = -6.27254058
d5 = -8.73646986
d6 = 2.26621137
d7 = 39.04984885
d8 = 179.71729691
d9 = 851.18388546
d10 = 786.53322125
d11 = -30882.44496207

(2.) Matrix elements of parameter using Taylor series 
  expansions (at r_e) using n^{th} order derivatives
<psi_2,11| mp(0) |psi_2,11> = 6.96793
<psi_2,11| mp(1) |psi_2,11> = 6.96793
<psi_2,11| mp(2) |psi_2,11> = 7.08671
<psi_2,11| mp(3) |psi_2,11> = 7.0901
<psi_2,11| mp(4) |psi_2,11> = 7.08515
<psi_2,11| mp(5) |psi_2,11> = 7.08535
<psi_2,11| mp(6) |psi_2,11> = 7.08536
<psi_2,11| mp(7) |psi_2,11> = 7.08535
<psi_2,11| mp(infty) |psi_2,11> = 7.08536
-------------------------------------------------------------------
``` 

**Plots**
---
The following plots are generated when the `enable_plot` option is set to `1` in the `pol_derivative.compute` command. These plots are exported in a pdf file called `output.pdf` in the same directory as the python script. (To use the `enable_plot` option, `matplotlib` must be available.)

<a href="img0"><img src="https://github.com/ankit7540/H2-PolarizabilityDerivatives/blob/master/image/fig0.png" align="center" height="350" ></a>


<a href="img0"><img src="https://github.com/ankit7540/H2-PolarizabilityDerivatives/blob/master/image/fig1.png" align="center" height="350" ></a>

<a href="img0"><img src="https://github.com/ankit7540/H2-PolarizabilityDerivatives/blob/master/image/fig2.png" align="center" height="350" ></a>
 
[f1]: http://chart.apis.google.com/chart?cht=tx&chl=\langle\psi_{v=0,J=0}|\bar{\alpha}|\psi_{v=0,J=0}\rangle
[f2]: http://chart.apis.google.com/chart?cht=tx&chl=\langle\psi_{v=2,J=1}|\gamma|\psi_{v=1,J=1}\rangle
[f3]: http://chart.apis.google.com/chart?cht=tx&chl=\langle\psi_{v=2,J=1}|\alpha_{\parallel}|\psi_{v=1,J=1}\rangle


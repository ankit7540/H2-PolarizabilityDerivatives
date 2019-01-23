'''# Purpose : Load the matrix containing polarizability and wavefunctions, interpolate the
# polarizability invariant if needed, determine the derivtive of the invariant at r_e and finally,
# compute the respective matrix elements using Taylor series expansions of parameter at r_e.'''

# Load necessary modules
import sys
import numpy as np
from scipy import integrate
from scipy.optimize import curve_fit
#-------------------------------------------------------------------

#********************************************************************
# python function to compute the first derivative at the first point
# and the last point of an array. For this computation, the first 4
# points are used for the derivative for the first data point. Similarly
# last 4 (x,y) points are used for the derivative at the last point.

def fd_ends(x, y):
    ''' Parameter:
        x       =       xaxis of the array
        y       =       y data of the array
                (x and y arrays must have atleast 4 elements)
        Returns = first derivative at the first point and the last point
    '''
    if (len(x) < 4 or len(y) < 4):
        print("Error : x and y arrays must have 4 elements")

    subx = np.zeros(4)
    suby = np.zeros(4)

    for i in range(0, 4):
        subx[i] = x[i]
        suby[i] = y[i]

    fd1 = ((subx[1]*subx[3]+subx[2]*subx[3]+subx[1]*subx[2]-2*subx[0]*subx[1] \
          -2*subx[0]*subx[3]-2*subx[0]*subx[2]+3*subx[0]**2)/(-subx[3]+subx[0]) \
          /(-subx[1]+subx[0])/(-subx[2]+subx[0])*suby[0]-(-subx[2]+subx[0])    \
          *(-subx[3]+subx[0])/(subx[1]-subx[3])/(-subx[1]+subx[0])/(subx[1]-subx[2])\
          *suby[1]+(-subx[1]+subx[0])*(-subx[3]+subx[0])/(subx[2]-subx[3])          \
          /(subx[1]-subx[2])/(-subx[2]+subx[0])*suby[2]-(-subx[1]+subx[0])          \
          *(-subx[2]+subx[0])/(subx[2]-subx[3])/(subx[1]-subx[3])/(-subx[3]+subx[0])\
          *suby[3])

    for i in range(0, 4):
        subx[i] = x[int(i-4)]
        suby[i] = y[int(i-4)]
#        print (i, int(i-4))

    fdn = ((subx[1]-subx[3])*(subx[2]-subx[3])/(-subx[3]+subx[0])/(-subx[1] \
           +subx[0])/(-subx[2]+subx[0])*suby[0]-(-subx[3]+subx[0])*(subx[2]   \
           -subx[3])/(subx[1]-subx[3])/(-subx[1]+subx[0])/(subx[1]-subx[2])   \
           *suby[1]+(-subx[3]+subx[0])*(subx[1]-subx[3])/(subx[2]-subx[3])    \
           /(subx[1]-subx[2])/(-subx[2]+subx[0])*suby[2]-(-2*subx[0]*subx[3]  \
           -2*subx[1]*subx[3]-2*subx[2]*subx[3]+subx[0]*subx[1]+subx[0]       \
           *subx[2]+subx[1]*subx[2]+3*subx[3]**2)/(subx[2]-subx[3])/(subx[1]   \
           -subx[3])/(-subx[3]+subx[0])*suby[3])

    return(fd1, fdn)

#********************************************************************

# Spline function taken from Numerical Recipes in FORTRAN, page 109 and 110
# Numerical recipes in FORTRAN, Second Edition, Press, Teukolsky, Vetterling, Flannery
# Cambridge University Press, 1992
def spline(x, y, yp1, ypn):

    '''Parameters:
        x       =       1D vector of x-values in increasing order.
        y       =       1D vector of y-values
        yp1     =       first derivative of the interpolating function at the first segment
        ypn     =       first derivative of the interpolating function at the last segment

    '''
    nx = len(x)
    ny = len(y)

    if (nx == ny):
        n = nx
    else:
        print("Error : x and y data have different lengths in spline.")
        quit()

    u = np.zeros(n)
    y2 = np.zeros(n) # this is the output
    p = 0.0
    sig = 0.0

    if yp1 > 1e30:     # lower boundar condition 'natural'
        y2[0] = 0.0
        u[0] = 0.0
    else:              # specified first derivative
        y2[0] = -0.5
        u[0] = (3/(x[1]-x[0])) * (((y[1]-y[0])/ (x[1]-x[0])) - yp1)

    for i in range(1, n-1):     # Decomposition loop of tridiagonal algorithm. y2 and u are temporary
        sig = (x[i]-x[i-1])/(x[i+1]-x[i-1])
        p = (sig * y2[i-1]) + 2.0
        y2[i] = (sig-1.0)/p
        u[i] = ((6.0*((y[i+1]-y[i])/(x[i+1]-x[i])-(y[i]-y[i-1]) / (x[i]-x[i-1]))/(x[i+1]-x[i-1])-  \
               sig*u[i-1])/p)
        # print("first loop:",i)

    if ypn > 1e30:     # upper boundary condition 'natural'
        qn = 0.0
        un = 0.0
    else:              # specified first derivative
        qn = 0.5
        un = (3.0 /(x[n-1]-x[n-2]))*(ypn-((y[n-1]-y[n-2])/(x[n-1]-x[n-2])))

    y2[n-1] = (un-qn*u[n-2])/(qn*y2[n-2]+1.0)

    for k in range(n-2, -1, -1):           # from second last point to the second point
        y2[k] = y2[k]*y2[k+1]+u[k]        # backsubstitution loop of tridiagonal algorithm
        # print("loop 2 :",k)

    return y2


#************************************************************************
# Spline interpolation function taken from Numerical Recipes in FORTRAN, page 109 and 110
# Numerical recipes in FORTRAN, Second Edition, Press, Teukolsky, Vetterling, Flannery
# Cambridge University Press, 1992

def splint(xa, ya, y2a, x):
    ''' Parameters :
        xa      =       original x-axis 1D vector
        ya      =       original y-axis 1D vector
        y2a     =       output of the spline function
        x       =       new x axis, scalar
    '''

    nxa = len(xa)
    nya = len(ya)
    ny2a = len(y2a)

    if (nxa != nya or nxa != ny2a or nya != ny2a):
        print("Error : xa or ya or y2a have incorrect dimension(s).")
        quit()

    n = nxa

    klo = int(0)
    khi = int(n-1)
    k = int(0)
    h = 0.0
    element = 0.0

    while ((khi-klo) > 1):
        k = int((khi+klo)/2)
        element = xa[k]
#        print(element,xa[k],k,x)
        if (element > x):
            khi = k
        else:
            klo = k

    h = xa[khi] - xa[klo]

    if h == 0:
        print("Error : Bad xa input in splint")
        quit()

    a = (xa[khi]-x)/h
    b = (x-xa[klo])/h
    y = a*ya[klo]+b*ya[khi] + (((a**3)-a)*y2a[klo]+((b**3)-b)*y2a[khi])*(h**2)/6.0

    return y #returns the interpolated value

#************************************************************************

# Define the fitting function : polynomial of degree 11
def poly11(x, c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11):
    '''Polynomial function without scaling the x. Used for plotting'''
    return c0+(c1*x)+c2*(x**2)+c3*(x**3)+c4*(x**4)+c5*(x**5)+c6*(x**6)+c7*(x**\
    7)+c8*(x**8)+c9*(x**9)+c10*(x**10)+c11*(x**11)

#************************************************************************

# Define the fitting function : scaled polynomial of degree 11
def poly11_sc(x, c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11):
    '''Polynomial function with scaled x. Ensures better numerical accuracy'''
    return c0+(c1*x/3.0)+c2*((x/3.0)**2)+c3*((x/3.0)**3)+c4*((x/3.0)**4)+c5*((x/3.0)**5)+c6*\
    ((x/3.0)**6)+c7*((x/3.0)**7)+c8*((x/3.0)**8)+c9*((x/3.0)**9)+c10*((x/3.0)**10)+c11*((x/3.0)**11)

#************************************************************************

# Define the wrappend intergral function which computes the matrix element by evaluating the
#   the numerical integral

# y_parameter should be already interpolated to rwave to that the product may be computed

def compute_int(rwave, psi1, psi2, y_parameter, rMin, rMax):
    p1 = np.multiply(psi1, psi2)
    p2 = np.multiply(p1, rwave)
    p3 = np.multiply(p2, rwave)
    product = np.multiply(p3, y_parameter)

    derivative = fd_ends(rwave, product)
    secarray2 = spline(rwave, product, derivative[0], derivative[1])

    # function defining the integrand which uses the spline coef array to give interpolated values
    def integrand_ME(xpoint):
        '''integrand for the quadrature'''
        result = splint(rwave, product, secarray2, xpoint)
        return result

    res = integrate.quadrature(integrand_ME, rMin, rMax, tol=1.0e-6, vec_func=False, maxiter=1000)
    return res
    # output is a tuple, [0] = integral result, [1]=error

#************************************************************************

# Loading polarizability data and checks  ::::

# Load the polarizability data ( alpha_xx and alpha_zz)---
alpha_xx = np.loadtxt("./data/matrix_xxf.txt")
alpha_zz = np.loadtxt("./data/matrix_zzf.txt")
omega = np.loadtxt("./data/freq.txt")
distance = np.loadtxt("./data/distance.txt")
static_xx = np.loadtxt("./data/static_xx.txt")
static_zz = np.loadtxt("./data/static_zz.txt")


# check size of the arrays -------------------------------
#print("Dimensions of isotropy matrix :",alpha_xx.shape)
#print("Dimensions of anisotropy matrix :",alpha_zz.shape)
if not(alpha_xx.shape == alpha_zz.shape or len(omega) == alpha_xx.shape[0] or len(static_xx) == len(distance)
      or  len(static_zz) == len(distance)):
    print("Dimension check on polarizability data matrices or wavelength file failed.")
    quit()
else:
    print("Polarizability data dimension checked.")
    print("\n")
    omega_nm = (1e7/(omega*219474.6313702000))
    omega_A = (omega_nm*10)



    print("Give  rovibME.compute  command with parameters:")
    print("\trovibME.compute(molecule, bra_v, bra_J, ket_v, ket_J, lambda, unit of lambda, operator)")
    print('\t for example:  rovibME.compute("H2",0,2,0,4,488,"n","mp")  ')
    print('\t\t       rovibME.compute("D2",1,0,1,0,"static","n","all")  ')
    print("\t\t")

    print('\t\tmolecule = for H2 enter "H2", for D2 enter "D2", for HD enter "HD" ')
    print("\t\tv    = vibrational state, v=[0,4]")
    print("\t\tJ    = rotataional state, J=[0,15]")
    print('\t\tlambda   = wavelength in Hartree, nm or Angstrom, for static specify "s" or "static" here')
    print('\t\tunit of lambda =  for  Hartree           use "H" or "h"  ')
    print('\t\t\t          for  nanometers        use "n" or "nm" ')
    print('\t\t\t          for  Angstrom          use "a" or "A"  ')
    print('\t\t\t          if static property is asked then this parameter can be any of the three ')
    print("\t\tAvailable wavelength range: {0} - {1} Hartree;\n   \t\t\t\t\t    {2} - {3} nm; \n   \t\t\t\t\t    {4} - {5} Angstrom".
         format(round(omega[0], 4), round(omega[-1], 4), round(omega_nm[0], 4), round
         (omega_nm[-1], 4), round(omega_A[0], 4), round(omega_A[-1], 4)))
    print('\t\toperator	= alpha(perpendicular) = alpha_xx given by "xx" or "x" ')
    print('\t\t\t          alpha(parallel) = alpha_xx given by "zz" or "z" ')
    print('\t\t\t          isotropy or mean polarizability given by "iso" or "mp" or "mean" ')
    print('\t\t\t          anisotropy or polarizability difference or gamma given by "aniso" or "g"  or "diff" ')
    print('\t\t\t          for all the above use "all"   or  "All"or "ALL" ')

    print("...ready.")
#********************************************************************
# the actual function for computation of the rovibrational matrix element.
# vl, Jl, vr , Jr are numbers
# mol, wavelength unit and operator are string, hence need quotes.

def compute(mol, v, J, wavelength, wavelength_unit, operator):
    '''#  parameters:
    # mol  =    molecule (for H2 enter "H2", for D2 enter "D2", for HD enter "HD")
    # v   =    vibrational state for the bra, v = [0,4]
    # J   =    rotational state for the bra,  J = [0,15]
    # wavelength =  wavelength ( can be Hartree, nanometers or Angstrom)
    # wavelength_unit = specify unit using the specifier
                                ( for  Hartree           use "H" or "h"  )
                                ( for  nanometers        use "n" or "nm"  )
                                ( for  Angstrom          use "a" or "A"  )

    # operator   = property namely alpha_xx, alpha_zz, mean polarizability
                                   (isotropy)[\bar{alpha}], anisotropy[\gamma]
                                   Specify operator using the specifier.
                                 ( for  alpha_xx = alpha(?)         use "x"     or  "xx"  )
                                 ( for  alpha_zz          use "z"     or  "zz"  )
                                 ( for  isotropy          use "iso"   or  "mp" or "mean" )
                                 ( for  anisotropy        use "aniso" or  "g"  or "diff" )
                                 ( for  all the above     use "all"   or  "All"or "ALL"  )

    This function runs on both Python 2.7x and 3.x
    '''

    # set a dictionary for output array
    d = {'output':[0]}
    #----------------------------------------------------------------
    # interpolation function defined here and is used later.
    def interpolate2D_common(input2D, originalx, finalx):
        inputSize = input2D.shape
        tempx = np.zeros((len(originalx), 1))
        tempx[:, 0] = originalx
        col = np.zeros((len(originalx), 1))
        outputArray = np.zeros((1, inputSize[1]))
        for i in range(0, inputSize[1]):
            col[:, 0] = input2D[:, i]
            der = (0.0, 0.0)	# derivatives at first and last ends
            der = fd_ends(tempx, col)
            # print(i,der[0],der[1])
            secarray = spline(tempx, col, der[0], der[1])
            interp = splint(tempx, col, secarray, finalx)
            outputArray[:, i] = interp
        d['output'] = outputArray.T
    #----------------------------------------------------------------

    # Load the polarizability data ( alpha_xx and alpha_zz)
    alpha_xx = np.loadtxt("./data/matrix_xxf.txt")
    alpha_zz = np.loadtxt("./data/matrix_zzf.txt")
    omega = np.loadtxt("./data/freq.txt")
    dist = np.loadtxt("./data/distance.txt")
    distance = np.asarray(dist)
    omega_nm = (1e7/(omega*219474.6313702000)) # convert the original freq to nm
    static_xx = np.loadtxt("./data/static_xx.txt")
    static_zz = np.loadtxt("./data/static_zz.txt")

    # compute the isotropy(mean polarizability) and anisotropy (gamma)
    isotropy = np.absolute(2*(np.array(alpha_xx))+np.array(alpha_zz))/3
    anisotropy = np.absolute(np.array(alpha_zz)-np.array(alpha_xx))

    isotropy_static = (2*static_xx+static_zz)/3
    anisotropy_static = (static_zz-static_xx)

    # step 1: load the required wavefunctions ------------------------
    Wfn1 = "./wavefunctions/{0}v{1}J{2}_norm.txt".format(mol, v, J)
    Wfn2 = "./wavefunctions/{0}v{1}J{2}_norm.txt".format(mol, v, J)
    r_wave = "./wavefunctions/r_wave.txt"
    #print(Wfn1,Wfn2)
    if v < 0  or v > 4:
        print("Error : v value out of range. v = [0,4]. Exiting ")
        quit()

    if J < 0   or J > 15:
        print("Error : J value out of range. Jl and Jr =[0,15]. Exiting ")
        quit()

    if not (mol == "H2"  or mol == "HD" or mol == "D2"):
        print("Error : Incorrect molecule chosen. For H2 enter H2, for D2 enter D2, for HD enter HD. Use quotes. Exiting  ")
        quit()

    # Proceed to load wavefunctions.
    psi1 = np.loadtxt(Wfn1)
    psi2 = np.loadtxt(Wfn2)
    rwave = np.loadtxt(r_wave)

    #----------------------------------------------------------------
    # STATIC
    if (wavelength == "static" or wavelength == "Static" or wavelength == "s"):
        #print("\tStatic")
        n = 0
        if (operator == "mean" or operator == "mp" or operator == "iso"):
            param = isotropy_static
            name = ["isotropy (static)"]
            n = 1
        elif (operator == "diff" or operator == "g" or operator == "aniso"):
            param = anisotropy_static
            name = ["anisotropy (static)"]
            n = 1
        else:
            print("Error : Operator not correctly specified. Exiting ")
            quit()

    # DYNAMIC
    elif (isinstance(wavelength, (int, float))):
        wv = float(wavelength) # entered wavelength is a number

        if (wavelength_unit == "h" or wavelength_unit == "H"):
            omegaFinal = (1e7/(wv*219474.6313702000))
        elif (wavelength_unit == "n" or wavelength_unit == "nm"):
            omegaFinal = wv
        elif (wavelength_unit == "a" or wavelength_unit == "A"):
            omegaFinal = (wv/10)
        elif not (wavelength_unit == "h" or wavelength_unit == "H" or wavelength_unit == "n" or wavelength_unit == "nm"
                  or wavelength_unit == "a" or wavelength_unit == "A"):
            print("Message : Default unit of nm will be used.")
            omegaFinal = wv

        if omegaFinal < omega_nm[0] or omegaFinal > omega_nm[-1]:
            sys.exit("Error : Requested wavelength is out of range. Exiting ")

        print("Selected wavelength in nanometer : {0}, Hartree : {1}".format(round(omegaFinal, 6), round((1e7/(omegaFinal*219474.63137020)), 6)))

        n = 0

        if (operator == "mean" or operator == "mp" or operator == "iso"):
            param = isotropy
            name = ["isotropy"]
            n = 1
        elif (operator == "diff" or operator == "g" or operator == "aniso"):
            param = anisotropy
            name = ["anisotropy"]
            n = 1
        else:
            print("Error : Operator not correctly specified. Exiting")
            quit()

    elif not (wavelength == "static" or wavelength == "s" or isinstance(wavelength, (int, float))):
        print('Error : Incorrect specification of wavelength. Use number for dynamic property and "s" or "static" for static. Exiting')
        quit()

    #-------------------------------------------------------------------------
    for i in range(n):    # evaluation of  interpolation and integral

        # step 0: prepare parameter vector(s)
        parameter = np.zeros((len(distance), 1))
        # Static
        if (wavelength == "static" or wavelength == "Static" or wavelength == "s"):
            if not(n == 1):
                parameter = list[i]
            else:
                parameter = param
        else:

        # interpolate to the asked wavelength ----------------
            if not(n == 1):
                param = list[i]
                # print(param.shape, omegaFinal)
                interpolate2D_common(param, omega_nm, omegaFinal)
            else:
                interpolate2D_common(param, omega_nm, omegaFinal)
                # print(param.shape, omegaFinal)
            temp = d['output']
            parameter = temp[:, 0]
    #-------------------------------------------------------------------------
    # print general information

        if (wavelength == "static" or wavelength == "Static" or wavelength == "s"):
            wavelength = ""
            wavelength_unit = ""

        print("-------------------------------------------------------------------")
        print("Analysis for the derivative of the parameter at re defined by v,J")
        print("-------------------------------------------------------------------")

        print("Molecule = {0}". format(mol))
        print("Parameter = {0}". format(name))
        print("Wavelength = {0} {1}". format(wavelength, wavelength_unit))
        print("Rovibrational state : v={0}, J={1}". format(v, J))


    # Step 2 : Compute the expectation value of the inter-nuclear distance -------
        result = compute_int(rwave, psi1, psi2, rwave, 0.5, 3.0)
        re = result[0]
        print("Expectation value of inter-nuclear distance, r_e = {0} a.u.". format(round(re, 6)))
    # ----------------------------------------------------------------------------

	# Step 3 : Perform truncation of parameter and corresponding x-axis to 0.5--3.0 a.u.
        distance = distance[9:]
        distance = distance[:101]

        parameter = parameter[9:]
        parameter = parameter[:101]
        parameter = np.reshape(parameter, len(parameter))

        #print("distance array: first point={0}, last point={1} " .format(round(distance[0],9) \
        #, round(distance[-1],9)  ))
        #print("parameter array: first point={0}, last point={1} " .format(np.round(parameter[0]\
        #,9) , np.round(parameter[-1],9)  ))

    # Step 4 : Fitting the parameter over distance in the range 0.5--3.0 a.u. ----
        popt_sc, pcov_sc = curve_fit(poly11_sc, distance, parameter)

        popt = np.zeros(12)

        # fit coeffcients
        print("\n(1.) Fit coefficients (scaled back using 3^{n})\nparameter vs distance was fit using poly11_sc function")
        for j in range(len(popt_sc)):
            print("c{0} = {1}".format(j, round(popt_sc[j]/(3.0**j), 10)))
            popt[j] = popt_sc[j]/(3.0**j)    # rescale the fit coefficients

        # covariance of the fit coefficients
        #print(pcov)
    #-----------------------------------------------------------------------------


    # Step 5 : Obtain the numerical value of the derivative at re ----------------
    # These are described as g0, g1, .. g6, g7

    #  g0 is the interpolated value of parameter at re
    #  g1,.. g7 are obtained using the coeffs of polynomial fit

        gn = np.zeros(8)

        der2 = fd_ends(distance, parameter)
        secarray2 = spline(distance, parameter, der2[0], der2[1])

        gn[0] = splint(distance, parameter, secarray2, re)


        gn[1] = popt[1]+2*popt[2]*re+(3*popt[3]*re**2)+(4*popt[4]*re**3)+(5*popt[5]\
        *re**4)+(6*popt[6]*re**5)+(7*popt[7]*re**6)+(8*popt[8]*re**7)+(9*popt[\
        9]*re**8)+(10*popt[10]*re**9)+(11*popt[11]*re**10)

        gn[2] = 2*popt[2]+(2*3*popt[3]*re)+(3*4*popt[4]*re**2)+(4*5*popt[5]*re**3)+\
        (5*6*popt[6]*re**4)+(6*7*popt[7]*re**5)+(7*8*popt[8]*re**6)+(8*9*popt[\
        9]*re**7)+(9*10*popt[10]*re**8)+(10*11*popt[11]*re**9)

        gn[3] = (2*3*popt[3])+(2*3*4*popt[4]*re)+(3*4*5*popt[5]*re**2)+(4*5*6*     \
        popt[6]*re**3)+(5*6*7*popt[7]*re**4)+(6*7*8*popt[8]*re**5)+(7*8*9*    \
        popt[9]*re**6)+(8*9*10*popt[10]*re**7)+(9*10*11*popt[11]*re**8)

        gn[4] = (2*3*4*popt[4])+(2*3*4*5*popt[5]*re)+(3*4*5*6*popt[6]*re**2)+(4*5*6\
        *7*popt[7]*re**3)+(5*6*7*8*popt[8]*re**4)+(6*7*8*9*popt[9]*re**5)+(7*8\
        *9*10*popt[10]*re**6)+(8*9*10*11*popt[11]*re**7)

        gn[5] = (2*3*4*5*popt[5])+(2*3*4*5*6*popt[6]*re)+(3*4*5*6*7*popt[7]*re**2)+\
        (4*5*6*7*8*popt[8]*re**3)+(5*6*7*8*9*popt[9]*re**4)+(6*7*8*9*10*      \
        popt[10]*re**5)+(7*8*9*10*11*popt[11]*re**6)

        gn[6] = (2*3*4*5*6*popt[6])+(2*3*4*5*6*7*popt[7]*re)+(3*4*5*6*7*8*popt[8]* \
        re**2)+(4*5*6*7*8*9*popt[9]*re**3)+(5*6*7*8*9*10*popt[10]*re**4)+(6*7*\
        8*9*10*11*popt[11]*re**5)

        gn[7] = (2*3*4*5*6*7*popt[7])+(2*3*4*5*6*7*8*popt[8]*re)+(3*4*5*6*7*8*9*   \
        popt[9]*re**2)+(4*5*6*7*8*9*10*popt[10]*re**3)+(5*6*7*8*9*10*11*      \
        popt[11]*re**4)

        # Derivatives of parameter at re
        print("\n(2.) Derivatives of parameter at r_e")
        for j in range(len(gn)):
            print("g{0} = {1}".format(j, round(gn[j], 8)))

        # ----------------------------------------------------------------------------
        # ----------------------------------------------------------------------------

        # Generate Taylor series expansion of parameter at r_e using the derivatives at r_e
        #    designated at g0, g1, .. gn computed above

        expn = np.zeros(shape=(len(distance), 8))
        r_re = distance-re

        expn[:, 0] = np.full(len(distance), gn[0])
        expn[:, 1] = np.fromiter((gn[0]+gn[1]*r_re[x] for x in range(         \
        len(distance))), dtype="float")

        expn[:, 2] = np.fromiter((gn[0]+gn[1]*r_re[x]+0.5*gn[2]*(r_re[x]**2) for\
        x in range(len(distance))), dtype="float")

        expn[:, 3] = np.fromiter((gn[0]+gn[1]*r_re[x]+0.5*gn[2]*(r_re[x]**2)+ \
        (1/6)*gn[3]*(r_re[x]**3) for x in range(len(distance))), dtype="float")

        expn[:, 4] = np.fromiter((gn[0]+gn[1]*r_re[x]+0.5*gn[2]*(r_re[x]**2)+ \
        (1/6)*gn[3]*(r_re[x]**3)+(1/24)*gn[4]*(r_re[x]**4) for x in range(len \
        (distance))), dtype="float")

        expn[:, 5] = np.fromiter((gn[0]+gn[1]*r_re[x]+0.5*gn[2]*(r_re[x]**2)+ \
        (1/6)*gn[3]*(r_re[x]**3)+(1/24)*gn[4]*(r_re[x]**4)+(1/120)*gn[5]*     \
        (r_re[x]**5) for x in range(len(distance))), dtype="float")

        expn[:, 6] = np.fromiter((gn[0]+gn[1]*r_re[x]+0.5*gn[2]*(r_re[x]**2)+ \
        (1/6)*gn[3]*(r_re[x]**3)+(1/24)*gn[4]*(r_re[x]**4)+(1/120)*gn[5]*     \
        (r_re[x]**5)+(1/720)*gn[6]*(r_re[x]**6) for x in range(len(distance)))\
        , dtype="float")

        expn[:, 7] = np.fromiter((gn[0]+gn[1]*r_re[x]+0.5*gn[2]*(r_re[x]**2)+ \
        (1/6)*gn[3]*(r_re[x]**3)+(1/24)*gn[4]*(r_re[x]**4)+(1/120)*gn[5]*     \
        (r_re[x]**5)+(1/720)*gn[6]*(r_re[x]**6)+(1/5040)*gn[7]*(r_re[x]**7)   \
        for x in range(len(distance))), dtype="float")

    # ----------------------------------------------------------------------------
    # Computation of the matrix element for the Taylor series expansions of
    #   the parameter using the expn array(as columns)

        # truncate the rwave, psi1 and psi2 to r=0.5--3.0 a.u.
        rwave = rwave[:701]
        rwave = rwave[75:]

        psi1 = psi1[:701]
        psi1 = psi1[75:]

        psi2 = psi2[:701]
        psi2 = psi2[75:]

        parameter_ME = np.zeros(8)    # array to keep the matrix elements
        param_sc = np.zeros(len(rwave))   # interpolation of the expn to rwave

        # generate interpolation of the expn to rwave using cubic spline
        for j in range(len(parameter_ME)):
            param = np.zeros(len(distance))
            param = expn[:, j]

            der2 = fd_ends(distance, param)
            secarray2 = spline(distance, param, der2[0], der2[1])

            for k in range(0, len(rwave)):
                param_sc[k] = splint(distance, param, secarray2, rwave[k])

            res = compute_int(rwave, psi1, psi2, param_sc, 0.5, 3.0)
            parameter_ME[j] = res[0]

        print("\n(3.) Matrix elements of parameter using Taylor series \n  expansions (at r_e) using n^{th} order derivatives")
        for j in range(len(parameter_ME)):
            print("<psi_{0},{1}| {2}({3}) |psi_{4},{5}> = {6}".format(v, J, operator, j, v, J, round(parameter_ME[j], 6)))


        # compute the matrix element for the original parameter array (no approximation)
        der2 = fd_ends(distance, parameter)
        secarray2 = spline(distance, parameter, der2[0], der2[1])
        for j in range(0, len(rwave)):
            param_sc[j] = splint(distance, parameter, secarray2, rwave[j])

        res = compute_int(rwave, psi1, psi2, param_sc, 0.5, 3.0)
        parameter_trueME = res[0]

        print("<psi_{0},{1}| {2}(infty) |psi_{3},{4}> = {5}".format(v, J, operator, v, J, round(parameter_trueME, 6)))
        print("-------------------------------------------------------------------")

        # ----------------------------------------------------------------------------

        # OPTIONAL PLOTTING : USING  MATPLOTLIB
        plot = 1    # change value to 1 to enable plotting

        if plot == 1:

            import matplotlib.pyplot as plt
            from matplotlib.backends.backend_pdf import PdfPages

            plt.figure(0)

            txt = ("*Generated from 'polDerivative.py' on the \nGitHub Repository: H2-PolarizabilityDerivatives")
            subtxt = "v={0}, J={1}\nr$_{{e}}$={2} a.u.". format(v, J, round(re, 6))
            print(txt, subtxt)

            plt.title('Fitting of the parameter using polynomial function')
            plt.plot(distance, parameter, linewidth=5)
            plt.plot(distance, poly11(distance, *popt), linewidth=2.0)
            plt.xlabel('Inter-nuclear distance [a.u.]')
            plt.ylabel('{0}, {1} {2} {3}  [a.u.]'.format(name, wavelength, wavelength_unit, mol))
            plt.grid(True)
            plt.legend(('parameter', 'polynomial fit'), loc='upper left')

            plt.text(0.05, 0.00001, txt, fontsize=5, transform=plt.gcf().transFigure)
            #plt.text(0.05, -3.1, txt, fontsize=7, wrap=True)

            residual = parameter-poly11(distance, *popt)
            plt.figure(1)
            plt.title('Residual of the fit')
            plt.plot(distance, residual, label="residual")
            plt.xlabel('Inter-nuclear distance [a.u.]')
            plt.ylabel('residual')
            plt.grid(True)
            plt.legend(loc='upper left')
            plt.text(0.05, 0.00001, txt, fontsize=5, transform=plt.gcf().transFigure)

            plt.figure(2)
            plt.title('Taylor series expansions of parameter at $r_{e}$ using $n^{th}$-order derivatives')
            plt.plot(distance, expn[:, 0], label="0")
            plt.plot(distance, expn[:, 1], label="1")
            plt.plot(distance, expn[:, 2], label="2")
            plt.plot(distance, expn[:, 3], label="3")
            plt.plot(distance, expn[:, 4], label="4")
            plt.plot(distance, expn[:, 5], label="5")
            plt.plot(distance, expn[:, 6], label="6")
            plt.plot(distance, expn[:, 7], label="7")
            plt.plot(distance, parameter, color='black', label="original")
            plt.xlabel('Inter-nuclear distance [a.u.]')
            plt.ylabel('{0}, {1} {2} {3}  [a.u.]'.format(name, wavelength, wavelength_unit, mol))
            plt.grid(True)
            plt.legend(loc='upper left')
            plt.text(0.05, 0.00001, txt, fontsize=5, transform=plt.gcf().transFigure)
            plt.text(0.725, 0.00001, subtxt, fontsize=7, color='navy', transform=plt.gcf().transFigure)

            # exporting the plots as a PDF file
            pdf = PdfPages('output.pdf')
            nfig = plt.gcf().number
            #print(nfig)

            for fig in range(nfig+1): #loop over plots
                pdf.savefig(fig, bbox_inches="tight")
            pdf.close()

        # ----------------------------------------------------------------------------
compute("H2", 0, 0, 514.5, "nm", "g")
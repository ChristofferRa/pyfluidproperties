# -*- coding: utf-8 -*-
"""
Helper functions for iapwsif97 python implementation.

This module contains functions not directly related to IAPWSIF97 or other
fluid properties but are frequently used in the calculations.

List of functions and code structure:
    Itteration functions
        secant(func,x0,x1,a, tolerance (optional), max_itt (optional))
        bisection(func, x0,x1,a, tolerance (optional), max_itt (optional))
        hybrid_bisec_secant(func, x0, x1, a, int_tolerance, tolerance (optional), max_itt (optional))
    Polynomial functions
        gen_poly_data(x_data,beta)
        gen_poly(x,beta)
        
@author: Christoffer Rappmann, christoffer.rappmann@gmail.com
"""
import numpy as np

############################
### Itteration functions ###
############################

def secant(func,x0,x1,a, tolerance = 0.0001, max_itt = 100):
    """
    Nummerical method for solving the equation f(x) = a
    where a is the wanted solution and x is the uknown variable. 
    
    Uses modified Secant method, https://en.wikipedia.org/wiki/Secant_method
    
    Use with care, finds one solution, there may be several. 
    It does not always converge. sensitive to initial guesses.
    
    Parameters
    ----------
    func :  function        f(x), wanted function.
    x0 :    double          lower first guess for x.
    x1 :    double          upper first guess for x.
    a :     double          wanted value a for f(x) = a.

    Returns
    -------
    x :     double          Solution for x, f(x) = a.

    """
    itt = 0
    x = x0 # is initial guess is good enough and while-loop never entered
    
    while np.abs(func(x1)-a) > tolerance and itt < max_itt:
        x = (((x1-x0)/(func(x1)-func(x0)))*(a-func(x1))+x1)
        x0 = x1
        x1 = x
        
        itt = itt + 1
        if itt >= max_itt:
            print(f'\nWarning convergence target not achievied, Error = {(np.abs(func(x)-a))}')

    return x

def bisection(func, x0,x1,a, tolerance = 0.001, max_itt = 100):
    """
    Nummerical method for solving the equation f(x) = a
    where a is the wanted solution and x is the uknown variable. 
    
    Uses modified Bisection method, https://en.wikipedia.org/wiki/Bisection_method
    
    Use with care, finds one solution, there may be several. 
    The solution must be within x0 and x1.
    
    Parameters
    ----------
    func :  function        f(x), wanted function.
    x0 :    double          lower first guess for x.
    x1 :    double          upper first guess for x.
    a :     double          wanted value a for f(x) = a.

    Returns
    -------
    x :     double          Solution for x, f(x) = a.

    """
    
    itt = 0
    x = x0
    
    while np.abs(func(x)-a) > tolerance and itt < max_itt:
                        
        # Calculation mid-point to create two bisections
        x = (x1 + x0)/2
        
        # Chose the bisection within which the solution lies
        # (One negative and one positive)
        if (func(x)-a)*(func(x0)-a) > 0: # if positive both terms have the same sign. 
            # solution within bisection x - x1, replace x0 with x
            x0 = x
        else:
            # solution within bisection x0 - x, replace x1 with x
            x1 = x
            
        itt = itt + 1
        if itt >= max_itt:
            print(f'\nWarning convergence target not achievied, Error = {(np.abs(func(x)-a))}')
        
    return x

def hybrid_bisec_secant(func, x0, x1, a, int_tolerance, tolerance = 0.0001, max_itt = 100):
    """
    Nummerical method for solving the equation f(x) = a
    where a is the wanted solution and x is the uknown variable. 
    
    Uses both the Bisetion and secant method
     - Bisection method, https://en.wikipedia.org/wiki/Bisection_method
     - Secant method, https://en.wikipedia.org/wiki/Secant_method
    
    Setting an intermidiary tolerance, the bisection method provides a better initial
    guess (window) for the secant method to work whithin. A good value for the
    intermeadiate tolerance needs to be determined from case to case
    
    Use with care, finds one solution, there may be several. 
    The solution must be within x0 and x1.
    
    Parameters
    ----------
    func :  function        f(x), wanted function.
    x0 :    double          lower first guess for x.
    x1 :    double          upper first guess for x.
    a :     double          wanted value a for f(x) = a.

    Returns
    -------
    x :     double          Solution for x, f(x) = a.

    """
    
    # First itterate using Bisection method
    x = bisection(func, x0,x1,a, int_tolerance, max_itt)
    
    # Continue with the secant method
    x = secant(func, x-int_tolerance,x+int_tolerance,a, tolerance, max_itt)
    
    return x

############################
### Polynomial functions ###
############################

def poly_reg(x_data, y_data, degree):
    """
    Polynomial regression using linear algebra, Ordinary least squares
    
    https://en.wikipedia.org/wiki/Ordinary_least_squares
    https://en.wikipedia.org/wiki/Regression_analysis
    https://personal.math.ubc.ca/~pwalls/math-python/linear-algebra/linear-algebra-scipy/
    
    Estimates parameters for n-th degre polynomial for a given data-set of x and y-data

    Parameters
    ----------
    x_data : list/np-array      Independent variable, x.
    y_data : list/np-array      Dependent variable y.
    degree : int                degree of polynomial function.

    Returns
    -------
    beta :      np-array        parameters of polynomial.
    r_squared : float           R^2 of fitted polynomial

    """
    # Import dependents
    import numpy as np
    import scipy.linalg as la
    
    # Initialize
    n = len(y_data)             # number of data points
    p = degree + 1              # number of modell parameters
    X_mat = np.zeros([n,p])     # init regressor matrix
    
    beta = np.array(p)
    r_squared = 0.0
    
    if n <= p: # Check for n > p
        print('Error, number of data points, n must be larger than the number of coefficients, p. i.e. n > p')
    else:
        # create regressor matrix, X
        for j in range(0,n):
            for i in range(0,p):
                X_mat[j,i] = x_data[j]**i
    
        # Performe matrix operations to solve for beta-parameters
        beta = (la.inv(X_mat.T@X_mat)@X_mat.T)@y_data
    
        # Calculate regression standard error
        #ssr = np.sum(((y_data-X_mat@beta))**2)
        #s = (ssr/(n-p))**.5
        
        # Calculate Coefficient of determination
        # https://en.wikipedia.org/wiki/Coefficient_of_determination
        y_mean = np.mean(y_data)
        r_squared = np.sum((X_mat@beta - y_mean)**2) / np.sum((y_data - y_mean)**2)
        
    
    return beta, r_squared


def gen_poly_data(x_data,beta):
    """
    Generates list of data-points, y,  for an arbitrary polynomial y = beta(n)*x^n
    
    Parameters
    ----------
    x_data : List        independent variable, x.
    beta : list/np-array parameters of polynomial.

    Returns
    -------
    y_data : list        dependent variable y data.

    """
    
    y_data = []
    for x in x_data:
        y_temp = 0.0
        for i,b in enumerate(beta):
            y_temp = y_temp + b*x**i
        y_data.append(y_temp)
    return y_data


def gen_poly(x,beta):
    """
    Generates data-point, y,  for an arbitrary polynomial y = beta(n)*x^n
    
    Parameters
    ----------
    x : double        independent variable, x.
    beta : list/np-array parameters of polynomial.

    Returns
    -------
    y : double        dependent variable y.

    """
    
    y = 0.0
    for i,b in enumerate(beta):
        y = y + b*x**i

    return y
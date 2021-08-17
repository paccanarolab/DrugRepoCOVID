import numpy.matlib
from numpy import linalg as LA
import pandas as pd 
import numpy as np

__author__ = 'diegogaleano'
__email__  = 'diego.galeano@fgv.br'
__date__  = '19-07-2021'

    
def DecompositionAlgorithm(X, k, alphas, masks):

    # maximum number of iterations of the mult algorithm
    maxiter = 2000    
    # alternative convergence criteria
    tolx = 1e-4   
    epsilon = np.finfo(float).eps
    sqrteps = np.sqrt(epsilon);
    # variance for initialization
    variance = 0.01 
    
    # Get the dimensions
    (ndrugs,nvirus) = X.shape
       
    # initialization
    W0 = np.random.uniform(0,np.sqrt(variance),(ndrugs,k))
    H0 = np.random.uniform(0,np.sqrt(variance),(k,nvirus)) 
    
    # normalization
    H0 = np.divide(H0, np.matlib.tile(np.array([np.sqrt(np.sum(np.power(H0,2),1))]).transpose(), (1, nvirus)))
    
    # cost function
    J = list()
    for iteration in range(maxiter):
        numer = 0
        denom = 0
        
        for i in masks.keys():
            numer += alphas[i]*np.multiply(masks[i], X)
            denom += alphas[i]*np.multiply(masks[i], np.dot(W0,H0))
            
        numer = np.dot(numer, H0.transpose())
        denom = np.dot(denom, H0.transpose()) + np.spacing(numer)
        
        # update W
        W = np.maximum(0, np.multiply(W0, np.divide(numer, denom)))
        
        # Delete negative values due to machine precision.
        W.clip(min = 0)

        numer = 0
        denom = 0
        
        for i in masks.keys():
            numer += alphas[i]*np.multiply(masks[i], X)
            denom += alphas[i]*np.multiply(masks[i], np.dot(W,H0))
            
        numer = np.dot(W.transpose(), numer)
        denom = np.dot(W.transpose(), denom) + np.spacing(numer)
        
        # update H
        H = np.maximum(0, np.multiply(H0, np.divide(numer, denom)))
        
        # Delete negative values due to machine precision.
        H.clip(min = 0)

        # compute loss function
        loss = 0
        for i in masks.keys():
            loss +=  0.5*alphas[i]*LA.norm(np.multiply(masks[i], (X - np.dot(W,H))), 'fro')**2
        
        J.append(loss)

        # Get norm of difference and max change in factors
        dw = np.amax(np.abs(W-W0))/(sqrteps + np.amax(np.abs(W0)));
        dh = np.amax(np.abs(H-H0))/(sqrteps + np.amax(np.abs(H0)));
        delta = np.maximum(dw,dh)
      

        # Check for convergence
        if iteration > 1:
            if delta <= tolx:
                print('Iter', iteration, 'delta', delta)
                convergence = True
                break

        # Remember previous iteration results
        W0 = W
        H0 = np.divide(H,np.matlib.tile(np.array([np.sqrt(np.sum(np.power(H0,2),1))]).transpose(), (1, nvirus))) #normalise

   
    return [W,H,J]

	
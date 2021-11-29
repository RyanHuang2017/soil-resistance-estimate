# -*- coding: utf-8 -*-
"""
Created on Sat Nov 14 19:08:10 2020
This code is used to calculate the K value for the Modified Archimedes' Law
Theory.

@author: sichuan
"""
import numpy as np

def compInteg(func, left, right, step, mode):
    if (mode == 'midpoint'):
        node = left
        integration = 0
        while (node + step <= right):
            pieceInt = step*(func(node + step/2))
            integration = integration + pieceInt
            node = node + step

        if (node < right):
            pieceInt = (right - node)*(func((node + right)/2))
            integration = integration + pieceInt
            
    elif (mode == 'trapezoid'):
        node = left
        integration = 0
        while (node + step <= right):
            pieceInt = step*(func(node) + func(node + step))/2
            integration = integration + pieceInt
            node = node + step
            
        if (node < right):
            pieceInt = (right - node)*(func(node) + func(right))/2
            integration = integration + pieceInt
            
    elif (mode == 'simpson'):
        node = left
        integration = 0
        while (node + step <= right):
            pieceInt = step*(func(node) + 4*func(node + step/2) + func(node + step))/6
            integration = integration + pieceInt
            node = node + step
            
        if (node < right):
            pieceInt = (right - node)*(func(node) + 4*func((node + right)/2) + func(right))/6
            integration = integration + pieceInt
            
    else:
        print('Oops! this method is not found')
    return integration

def compK(phi):
    # integration parameters for eta
    left = 0.0
    right = 1.0
    step = 1.0e-3
    
    # integration parameters for psi
    a = 0.0
    b = np.pi/2
    dh = 1.0e-3
    
    beta = np.pi/4.0 - phi/2.0
    A = 2.0*(1.0 + np.sin(phi))/(1.0 - np.sin(phi)) * np.exp(np.pi * np.tan(phi))
    
    rc = lambda eta : (1.0 + 2.0*(1.0 - eta)*np.exp(np.pi/2.0 * np.tan(phi))/np.tan(beta))
    rd = lambda eta : (1.0 + 1.0*(1.0 - eta)*np.exp(np.pi/2.0 * np.tan(phi))/np.tan(beta))
    re = lambda eta : eta
    
    node = left
    integration = 0
    while (node + step <= right):
        Zint = 0.0
        fz = lambda psi : (-(1.0-(node + step/2.0))*np.cos(psi + beta)*np.exp(psi*np.tan(phi)))/(np.cos(phi)*(np.sin(beta)+(1.0 - (node + step/2.0))*np.sin(psi - beta)*np.exp(psi*np.tan(phi))))
        Zint = compInteg(fz, a, b, dh, 'midpoint')
        fk = lambda eta : eta*(((rc(eta)**(1.0 + np.tan(beta)**2))/(re(eta)*rd(eta)**(np.tan(beta)**2)))**np.sin(phi)) * np.exp(np.sin(phi)*np.tan(beta)*Zint)
        pieceInt = step*fk(node + step/2.0)
        integration = integration + pieceInt
        node = node + step
    if (node < right):
        pieceInt = (right - node)*(fk((node + right)/2))
        integration = integration + pieceInt
    
    K = A*integration
    return K
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 14 16:09:40 2013

@author: Siegfried
"""
from numpy import dot, append, zeros, arange, flipud, divide

def levinson_algorithm(r):
    a = zeros(len(r));
    k = zeros(len(r));
    
    for m in arange(0,len(r)-1):        
        alpha = -dot( flipud(r[0:m+1]), append(flipud(a[1:m+1]),1.0) );        
        mu = -dot(flipud(r[0:m+1]),append(a[1:m+1],0.0)) - r[m+1];        
        k[m] = -divide(mu,alpha);        
        a[1:m+2] = append(a[1:m+1],0.0) + k[m]*append(flipud(a[1:m+1]),1.0);
        
    a[0] = 1;
    
    return (a,k)
    

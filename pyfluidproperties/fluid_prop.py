# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 21:38:45 2022

@author: Christoffer Rappmann, christoffer.rappmann@gmail.com

Fluid thermodynamic properties superclass
"""

class fluid_prop:
    """
    Superclass containing template of basic thermodynamic properties for a fluid
    """
    
    def __init__(self, p, T):
        self.p = p # Pressure (Pa)
        self.T = T  # Temperature (K)
        
        self.v = 0
        self.h = 0
        self.u = 0
        self.s = 0
        self.cp = 0
        self.cv = 0
        self.w = 0
        self.x = 0
        self.my = 0
        self.tc = 0
        self.st = 0
    
    @staticmethod
    def v_pt(p, T):
            
        v = 0
        
        return v
    
    @staticmethod
    def h_pt(p, T):

        h = 0
        
        return h
    
    @staticmethod
    def u_pt(p, T):
        
        u = 0
        
        return u
    
    @staticmethod
    def s_pt(p, T):
        
        s = 0
        
        return s
    
    @staticmethod
    def cp_pt(p, T):
        
        cp = 0
        
        return cp
    
    @staticmethod
    def cv_pt(p, T):
        
        cv = 0
        
        return cv
    
    @staticmethod
    def w_pt(p, T):
        
        w = 0
        
        return w
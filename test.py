#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 11 07:45:58 2017

@author: jiahaozhang
"""

# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import matplotlib.pyplot as plt
import numpy as np
import math
import cmath
def dis(x,y,z):
    return math.sqrt(x**2+y**2+z**2);
def px(phi):
    re=np.zeros((ngrid,ngrid,ngrid),dtype='complex');
    for i in range(0,ngrid-1,1):
        for j in range(0,ngrid,1):
            for k in range(0,ngrid,1):
                re[i,j,k]=1/cmath.sqrt(-1)*(phi[i+1,j,k]-phi[i,j,k])/dx;
    return re;
def py(phi):
    re=np.zeros((ngrid,ngrid,ngrid),dtype='complex');
    for i in range(0,ngrid,1):
        for j in range(0,ngrid-1,1):
            for k in range(0,ngrid,1):
                re[i,j,k]=1/cmath.sqrt(-1)*(phi[i,j+1,k]-phi[i,j,k])/dy;
    return re;
def pz(phi):
    re=np.zeros((ngrid,ngrid,ngrid),dtype='complex');
    for i in range(0,ngrid,1):
        for j in range(0,ngrid,1):
            for k in range(0,ngrid-1,1):
                re[i,j,k]=1/cmath.sqrt(-1)*(phi[i,j,k+1]-phi[i,j,k])/dy;
    return re;
def xop(phi):
    re=np.zeros((ngrid,ngrid,ngrid),dtype='complex');
    xscope=np.linspace(-1*scope,scope,ngrid);
    for i in range(0,ngrid,1):
        for j in range(0,ngrid,1):
            for k in range(0,ngrid,1):
                re[i,j,k]=xscope[i]*phi[i,j,k];
    return re;
def yop(phi):
    re=np.zeros((ngrid,ngrid,ngrid),dtype='complex');
    yscope=np.linspace(-1*scope,scope,ngrid);
    for i in range(0,ngrid,1):
        for j in range(0,ngrid,1):
            for k in range(0,ngrid,1):
                re[i,j,k]=yscope[j]*phi[i,j,k];
    return re;
def zop(phi):
    re=np.zeros((ngrid,ngrid,ngrid),dtype='complex');
    zscope=np.linspace(-1*scope,scope,ngrid);
    for i in range(0,ngrid,1):
        for j in range(0,ngrid,1):
            for k in range(0,ngrid,1):
                re[i,j,k]=zscope[k]*phi[i,j,k];
    return re;
def innerprod(phi1,phi2):
    return np.vdot(np.conjugate(phi1.reshape((ngrid**3,1,1))),phi2.reshape((ngrid**3,1,1)))*dx*dy*dz;
ngrid=500;
scope=5;
xgrid=np.linspace(-1*scope,scope,ngrid);
dx=xgrid[2]-xgrid[1];
ygrid=np.linspace(-1*scope,scope,ngrid);
dy=ygrid[2]-ygrid[1];
zgrid=np.linspace(-1*scope,scope,ngrid);
dz=zgrid[2]-zgrid[1];
phiall=np.zeros((4,ngrid,ngrid,ngrid),dtype='complex');
for i in range(0,ngrid,1):
    for j in range(0,ngrid,1):
        for k in range(0,ngrid,1):
            phiall[0,i,j,k]=math.exp(-1*dis(xgrid[i],ygrid[j],zgrid[k])/2);
            phiall[1,i,j,k]=xgrid[i]*math.exp(-1*dis(xgrid[i],ygrid[j],zgrid[k])/2);
            phiall[2,i,j,k]=ygrid[j]*math.exp(-1*dis(xgrid[i],ygrid[j],zgrid[k])/2); 
            phiall[3,i,j,k]=zgrid[k]*math.exp(-1*dis(xgrid[i],ygrid[j],zgrid[k])/2); 
co1=math.sqrt(1/np.sum(np.power(np.reshape(phiall[0],(ngrid*ngrid*ngrid,1,1)),2)*dx*dy*dz));
co2=math.sqrt(1/np.sum(np.power(np.reshape(phiall[1],(ngrid*ngrid*ngrid,1,1)),2)*dx*dy*dz));
co3=math.sqrt(1/np.sum(np.power(np.reshape(phiall[2],(ngrid*ngrid*ngrid,1,1)),2)*dx*dy*dz));
co4=math.sqrt(1/np.sum(np.power(np.reshape(phiall[3],(ngrid*ngrid*ngrid,1,1)),2)*dx*dy*dz));
phiall[0]=co1*phiall[0];
phiall[1]=co2*phiall[1];
phiall[2]=co3*phiall[2];
phiall[3]=co4*phiall[3];
lz=np.zeros((4,4),dtype='complex');
lx=np.zeros((4,4),dtype='complex');
ly=np.zeros((4,4),dtype='complex');
for i in range(4):
    for j in range(4):
        print("get one")
        lz[i,j]=innerprod(phiall[i],xop(py(phiall[j]))-yop(px(phiall[j])));
        lx[i,j]=innerprod(phiall[i],yop(pz(phiall[j]))-zop(py(phiall[j])));
        ly[i,j]=innerprod(phiall[i],zop(px(phiall[j]))-xop(pz(phiall[j])));
data_file=open("re.txt","a");
data_file.write(lx);
data_file.write("\n");
data_file.write(ly);
data_file.write("\n");
data_file.write(lz);
data_file.close();

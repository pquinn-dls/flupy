import math
import numpy as np
import matplotlib.pyplot as plt

"""

Code for RGB Triangle or Circle taken from espacergb.py
http://www.f-legrand.fr/scidoc/docmml/image/niveaux/couleurs/couleurs.html
Note the french language rgb instead of rgb

I modified their colour blending in the maxwell triangle 
by introducing rescaleRGB to render it like a standard XRF image

Used to produce an triangle or circle annotation for a diagram or plot....

"""

           
def plot_rgb_triangle(Rlabel=None,Rxy=(0.35, 0.15),\
                 Glabel=None,Gxy=(0.0,0.8),\
                 Blabel=None,Bxy= (-0.35,0.15),\
                 Textcolor='black',
                 npoints = 1200,
                 fontsize=20,weight='heavy',name='arial',outline=False):
    a = 1.0/math.sqrt(3)
    if outline:
        plt.scatter([a,0,-a],[0,1,0],s=40,c=[(1,0,0),(0,1,0),(0,0,1)])
        plt.plot([a,0,-a,a],[0,1,0,0])
    img = fill_maxwell(npoints)
    plt.axis('off')
    if(Blabel!=None):
        plt.annotate(Blabel, xy=Bxy, \
                     xycoords='data',textcoords='data',\
                     horizontalalignment='center',\
                     verticalalignment='top',color=Textcolor,\
                     fontsize=fontsize,weight=weight,name=name)
    if(Glabel!=None):
        plt.annotate(Glabel, xy=Gxy, \
                     xycoords='data',textcoords='data',\
                     horizontalalignment='center',\
                     verticalalignment='top',\
                     color=Textcolor,fontsize=fontsize,weight=weight,name=name)
    if(Rlabel!=None):
        plt.annotate(Rlabel, xy=Rxy, \
                     xycoords='data',textcoords='data',\
                     horizontalalignment='center',\
                     verticalalignment='top',\
                     color=Textcolor,fontsize=fontsize,weight=weight,name=name)
    plt.imshow(img,origin='lower',extent=[-a,a,0.0,1.0])
            
           
def fill_maxwell(npoints=1200):
    """
    
    Generate a maxwell triangle with the range 
    
    
    """
    Nlines=npoints
    Ncol  =npoints
    img = np.zeros((Nlines,Ncol,4))
    dx = 2.0/(Ncol-1)
    dy = 1.0/(Nlines-1)
    for i in range(Ncol-1):
        for j in range(Nlines-1):
            x =-1.0+i*dx
            y = j*dy
            v = y
            #r = (x+1-v)
            r = (x+1.-v)/2.0
            b = 1.0-v-r
            r,v,b=_rescaleRGB(r,v,b)
            if (r>=0) and (r<=1.0) and (v>=0) and (v<=1.0) and (b>=0) and (b<=1.0):
                #t,s,l=rgb2tsl([r,v,b])
                img[j][i] = np.array([r,v,b,1.0])
            else:
                img[j][i] = np.array([1.0,1.0,1.0,0.0])

    return img
                    
            
def _rgb2tsl(rgb):
        r=rgb[0]*1.0
        v=rgb[1]*1.0
        b=rgb[2]*1.0
        somme = r+v+b
        r=r/somme
        v=v/somme
        b=b/somme
        Max = max(r,v,b)
        Min = min(r,v,b)
        C=Max-Min
        L = Max
        if L==0:
            return [0,0,0]
        S = C/L
        if C==0:
            return [0,0,L]
        if max==r:
            T = 60.0*(v-b)/C % 360
        elif max==v:
            T = 120.0+60.0*(b-r)/C
        else:
            T = 240.0+60.0*(r-v)/C
        return [T,S,L]
            
def _tsl2rgb(tsl):
        T=tsl[0]*1.0
        S=tsl[1]*1.0
        L=tsl[2]*1.0
        C=L*S
        Min = L-C
        if (T>300) and (T<=360):
            r = L
            v = Min
            b = v+C*(360.0-T)/60
        elif (T>=0) and (T<=60):
            r = L
            b = Min
            v = b+C*(T/60)
        elif (T>60) and (T<=120):
            v = L
            b = Min
            r = b+C*(120.0-T)/60
        elif (T>120) and (T<=180):
            v = L
            r = Min
            b = r+C*(T-120.0)/60
        elif (T>180) and (T<=240):
            b = L
            r = Min
            v = r+C*(240.0-T)/60
        else:
            b = L
            v = Min
            r = v+C*(T-240.0)/60
        return [r,v,b]
            
def circleTSL(npoints=1200):
    Ncol = npoints
    Nlines = npoints
    img = np.zeros((Nlines,Ncol,4))
    dx = 2.0/(Ncol-1)
    dy = 2.0/(Nlines-1)
    rad2deg = 180.0/math.pi
    for i in range(Ncol):
        for j in range(Nlines):
            x=-1.0+i*dx
            y=-1.0+j*dy
            r = math.sqrt(x*x+y*y)
            if r<1.0:
                if x==0:
                    if y>0:
                        a = 90.0
                    else:
                        a = -90.0
                else:
                    a = math.atan(y/x)*rad2deg
                if x<0:
                    a = a+180.0
                a = a % 360
                rgb = _tsl2rgb([a,r,1.0])
                img[j][i] = np.array([rgb[0],rgb[1],rgb[2],1.0])
            else:
                img[j][i] = np.array([1.0,1.0,1.0,0.0])
    return img

def plot_rgb_circle(Rlabel=None,Rxy=(0.66, 0.1),\
                 Glabel=None,Gxy=(-0.33,0.5),\
                 Blabel=None,Bxy= (-0.33, -0.50),\
                 Textcolor='black',
                 npoints = 1200,
                 fontsize=20,weight='heavy',name='arial'):
    img = circleTSL(npoints)
    plt.axis('off')
    if(Blabel!=None):
        plt.annotate(Blabel, xy=Bxy, \
                     xycoords='data',textcoords='data',\
                     horizontalalignment='center',\
                     verticalalignment='top',color=Textcolor,\
                     fontsize=fontsize,weight=weight,name=name)
    if(Glabel!=None):
        plt.annotate(Glabel, xy=Gxy, \
                     xycoords='data',textcoords='data',\
                     horizontalalignment='center',\
                     verticalalignment='top',\
                     color=Textcolor,fontsize=fontsize,weight=weight,name=name)
    if(Rlabel!=None):
        plt.annotate(Rlabel, xy=Rxy, \
                     xycoords='data',textcoords='data',\
                     horizontalalignment='center',\
                     verticalalignment='top',\
                     color=Textcolor,fontsize=fontsize,weight=weight,name=name)
    plt.imshow(img,origin='lower',extent=[-1,1,-1,1])        
            
            
def _rescaleRGB(r,g,b):
    """
    PDQ: Re-normalise the blend
    so that the max of the RGB is 1.0
    R,G,B = 0.5,0.5,0.2 becomes 1,1,0.4
    """
    cmax=max([r,g,b])
    if cmax>0.0:
        return r/cmax,g/cmax,b/cmax
    else:
        return 0,0,0



#plt.figure(figsize=(6,6))
#plot_rgb_triangle(Rlabel="Cr",Glabel="Co",Blabel="Ca",Textcolor='white')
#plt.savefig('triangle.png', transparent=True)
#plt.show()
#plot_rgb_circle(Rlabel="Cr",Glabel="Co",Blabel="Ca",Textcolor='white')
#plt.show()           

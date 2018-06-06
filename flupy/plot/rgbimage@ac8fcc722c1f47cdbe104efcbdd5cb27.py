import math
import numpy as np
import matplotlib.pyplot as plt
"""
Code for maxwell triangle taken from espaceRVB.py
http://www.f-legrand.fr/scidoc/docmml/image/niveaux/couleurs/couleurs.html
Note the french language rvb instead of rgb

I modified their colour blending in the maxwell triangle 
by introducing rescaleRGB

"""

           
def maxwell():
    a = 1.0/math.sqrt(3)
    plt.scatter([a,0,-a],[0,1,0],s=40,c=[(1,0,0),(0,1,0),(0,0,1)])
    plt.plot([a,0,-a,a],[0,1,0,0],'k-')
    plt.axis([-0.7,0.7,-0.2,1.2])
            
def point(rvb):
    somme = rvb[0]+rvb[1]+rvb[2]
    r = rvb[0]*1.0/somme
    v = rvb[1]*1.0/somme
    b = rvb[2]*1.0/somme
    plt.scatter([(r-b)/math.sqrt(3)],[v],s=50,c=(r,v,b))
            
def fill_maxwell(Rlabel=None,Glabel=None,Blabel=None,Tcolor='black'):
        Nlignes=300
        Ncol=300
        img = np.zeros((Nlignes,Ncol,4))
        dx = 2.0/(Ncol-1)
        dy = 1.0/(Nlignes-1)
        for i in range(Ncol-1):
            for j in range(Nlignes-1):
                x =-1.0+i*dx
                y = j*dy
                v = y
                #r = (x+1-v)
                r = (x+1-v)/2.0
                b = 1.0-v-r
                r,v,b=rescaleRGB(r,v,b)#print x,y,r,v,b 
                if (r>=0) and (r<=1.0) and (v>=0) and (v<=1.0) and (b>=0) and (b<=1.0):
                    #t,s,l=rvb2tsl([r,v,b])
                    img[j][i] = np.array([r,v,b,1.0])
                    #img[j][i] = numpy.array([r,v,b])
                    #img[j][i] = numpy.array([t,s,l,1.0])
                else:
                    img[j][i] = np.array([1.0,1.0,1.0,0.0])
                    #img[j][i] = numpy.array([1.0,1.0,1.0])
        a = 1.0/math.sqrt(3)

        plt.axis('off')
        if(Blabel!=None):
            plt.annotate(Blabel, xy=(-0.63, 0.0), xytext=(-0.25, 0.20), xycoords='data',textcoords='data',horizontalalignment='center',verticalalignment='top',color=Tcolor,fontsize=20)
        if(Glabel!=None):
            plt.annotate(Glabel, xy=(0.0,1.0), xytext=(0.0,0.60), xycoords='data',textcoords='data',horizontalalignment='center',verticalalignment='bottom',color=Tcolor,fontsize=20)
        if(Rlabel!=None):
            plt.annotate(Rlabel, xy=(0.63, 0.0), xytext=(0.25, 0.20), xycoords='data',textcoords='data',horizontalalignment='center',verticalalignment='top',color=Tcolor,fontsize=20)

        plt.imshow(img,origin='lower',extent=[-a,a,0.0,1.0])
            
def rvb2tsl(rvb):
        r=rvb[0]*1.0
        v=rvb[1]*1.0
        b=rvb[2]*1.0
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
            
def tsl2rvb(tsl):
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
            
def disqueTSL():
    Ncol = 300
    Nlignes = 300
    img = np.zeros((Nlignes,Ncol,4))
    dx = 2.0/(Ncol-1)
    dy = 2.0/(Nlignes-1)
    rad2deg = 180.0/math.pi
    for i in range(Ncol):
        for j in range(Nlignes):
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
                rvb = tsl2rvb([a,r,1.0])
                img[j][i] = np.array([rvb[0],rvb[1],rvb[2],1.0])
            else:
                img[j][i] = np.array([1.0,1.0,1.0,0.0])
    plt.imshow(img,origin='lower',extent=[-1,1,-1,1])        
            
            
            
def pointsTSV(dT,S):
    T = np.arange(0,360,dT)
    x = np.zeros(T.size)
    y = np.zeros(T.size)
    for k in range(T.size):
        rvb = tsl2rvb([T[k],S,1.0])
        somme = rvb[0]+rvb[1]+rvb[2]
        x[k] = (rvb[0]-rvb[2])/math.sqrt(3)/somme
        y[k] = rvb[1]/somme
    plt.plot(x,y,"k.")

def rescaleRGB(r,g,b):
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



def fill_twocolour(top=0,bottom=1,toplabel=None,bottomlabel=None):
        Nrow=300
        Ncol=300
        img = np.zeros((Nrow,Ncol,4))
        dx = 1.0/(Ncol-1)
        dy = 1.0/(Nrow-1)
        for i in range(Ncol-1):
            for j in range(Nrow-1):
                x =-1.0+i*dx
                y = j*dy
                
                r = y
                g = (1.0-r)*i*dx
                b = 0.0
                #r,g,b=rescaleRGB(r,g,b)#print x,y,r,v,b 
                if (r>=0) and (r<=1.0) and (g>=0) and (g<=1.0) and (b>=0) and (b<=1.0):
                    img[j][i] = np.array([r,g,b,1.0])
                else:
                    img[j][i] = np.array([1.0,1.0,1.0,0.0])
        a = 1.0/math.sqrt(3)

        plt.axis('off')
#         if(Blabel!=None):
#             plt.annotate(Blabel, xy=(-0.63, 0.0), xytext=(-0.63, -0.05), xycoords='data',textcoords='data',horizontalalignment='center',verticalalignment='top')
#         if(Glabel!=None):
#             plt.annotate(Glabel, xy=(0.0,1.0), xytext=(0.0,1.05), xycoords='data',textcoords='data',horizontalalignment='center',verticalalignment='bottom')
#         if(Rlabel!=None):
#             plt.annotate(Rlabel, xy=(0.63, 0.0), xytext=(0.63, -0.05), xycoords='data',textcoords='data',horizontalalignment='center',verticalalignment='top')

        plt.imshow(img,origin='lower',extent=[-a,a,0.0,1.0])



#def normalize2DArray(2darr):
#    2darr /= np.max(np.abs(2darr),axis=0)
#    image *= (255.0/image.max())
    
def rgbtrueimage(R,G,B):
    rows,cols = R.shape
    R*= (1.0/R.max())
    G*= (1.0/G.max())
    B*= (1.0/B.max())
    img = np.zeros((rows,cols,3))
    img[:,:,0]=R
    img[:,:,1]=G
    img[:,:,2]=B
    plt.imshow(img)        
    

def rgbnormimage(R,G,B):
    rows,cols = R.shape
    R*= (1.0/R.max())
    G*= (1.0/G.max())
    B*= (1.0/B.max())
    img = np.zeros((rows,cols,3))
    img[:,:,0]=R
    img[:,:,1]=G
    img[:,:,2]=B
    plt.imshow(img)        
    

    

plt.figure(figsize=(6,6))
maxwell()
#fill_twocolour(top=0,bottom=1,toplabel=None,bottomlabel=None)
fill_maxwell(Rlabel="Cr",Glabel="Co",Blabel="Ca",Tcolor='black')
plt.show()

#N=200
#R=np.random.random((N,N))
#G=np.random.random((N,N))
#B=np.random.random((N,N))
#rgbimage(R,G,B)
#plt.show()           
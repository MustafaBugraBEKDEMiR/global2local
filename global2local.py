import numpy as np  #calling neccessary libraries.NumPy and Math.
import math as m    #

x,y,z=(X,Y,Z)     #Input values

a,b=(6378137.0,6356752.3141) #Parameters of GRS80 ellipsoid

e2=(a*a-b*b)/(a*a) # Eccentrcity formula 

#By Pollard Method Algorithm#

p=np.sqrt(x*x+y*y) #auxiliary point

z0=b*z/(m.sqrt(x*x+y*y+z*z)) #the radius of curvature in the prime vertical


k=m.sqrt((x*x+y*y)+(z+e2*z0)*(z+e2*z0))

l=x/k
m=y/k
n=(z+e2*z0)/k
r=l+e2*n*n
s=l*x+m*y+(a/b)*(a/b)*n*z
t=x*x+y*y+(a/b)*(a/b)*z*z-a*a
h=(np.sqrt((s*s)-(r*t)))/r          # ellipsoidal height

if h > 1000000:                     #
    h=(s-np.sqrt((s*s)-(r*t)))/r
                                    #h should be in -1000< h < 1000000
elif h < -1000:
    h=(s+np.sqrt((s*s)-(r*t)))/r    #

z0=z-n*h
B=np.arctan((z+e2*z0)/p)            #ellipsoidal latitude
L=2*np.arctan(y/(x+p))              #ellipsoidal longitude

A=np.matrix([[-np.sin(B)*np.cos(L), -np.sin(L), np.cos(B)*np.cos(L)],
            [-np.sin(B)*np.sin(L), np.cos(L), np.cos(B)*np.sin(L)],
            [np.cos(B), 0, np.sin(B)]])


v=np.array((B,L,h))

print(A.dot(v))

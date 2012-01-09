#!/usr/bin/python
#from __future__ import division
from matplotlib.pyplot import *
from numpy import *
from matplotlib import rc
from math import sqrt

rc('text', usetex=True)
rc('font', family='serif')

n = 6;
xGridSize = 50;
x = zeros(xGridSize+1);
xn = zeros(n+1);
fx = zeros(xGridSize+1);
fxn = zeros(n+1);
hn = (1.0+1.0)/n;
h = (1.0+1.0)/xGridSize;
EN2 = 0.0;

for i in range(0, n + 1):
   xn[i] = -1.0 + i*hn;
   fxn[i]= 1.0/(25.0*(xn[i]**2)+1.0);
else:
   print "Loop 1 is finished"


for i in range(0, xGridSize + 1):
   x[i] = -1.0 + i*h;
   for j in range(1, n + 1):
      if (x[i] <= xn[j]):
         break
   fx[i] = (x[i]-xn[j-1])/(xn[j]-xn[j-1])*fxn[j] + (x[i]-xn[j])/(xn[j-1]-xn[j])*fxn[j-1];
else:
   print "Loop 2 is finished"

for i in range(0, xGridSize + 1):
   f = 1.0/(25.0*(x[i]**2)+1.0);
   EN2 = EN2 + ((fx[i]-f)/f)**2;
else:
   EN2 = sqrt(EN2) / n;
   print "EN2 = ", EN2


figure(1)
plot(x, fx, 'r', label='n = %i' %(n))
plot(xn, fxn, 'ro' )

###########

n = 8;
x = zeros(xGridSize+1);
xn = zeros(n+1);
fx = zeros(xGridSize+1);
fxn = zeros(n+1);
hn = (1.0+1.0)/n;
h = (1.0+1.0)/xGridSize;
EN2 = 0.0;

for i in range(0, n + 1):
   xn[i] = -1.0 + i*hn;
   fxn[i]= 1.0/(25.0*(xn[i]**2)+1.0);
else:
   print "Loop 1 is finished"


for i in range(0, xGridSize + 1):
   x[i] = -1.0 + i*h;
   for j in range(1, n + 1):
      if (x[i] <= xn[j]):
         break
   fx[i] = (x[i]-xn[j-1])/(xn[j]-xn[j-1])*fxn[j] + (x[i]-xn[j])/(xn[j-1]-xn[j])*fxn[j-1];
else:
   print "Loop 2 is finished"

for i in range(0, xGridSize + 1):
   f = 1.0/(25.0*(x[i]**2)+1.0);
   EN2 = EN2 + ((fx[i]-f)/f)**2;
else:
   EN2 = sqrt(EN2) / n;
   print "EN2 = ", EN2

figure(1)
plot(x, fx, 'b', label='n = %i' %(n))
plot(xn, fxn, 'bo' )


#################################

n = 10;
x = zeros(xGridSize+1);
xn = zeros(n+1);
fx = zeros(xGridSize+1);
fxn = zeros(n+1);
hn = (1.0+1.0)/n;
h = (1.0+1.0)/xGridSize;
EN2 = 0.0;

for i in range(0, n + 1):
   xn[i] = -1.0 + i*hn;
   fxn[i]= 1.0/(25.0*(xn[i]**2)+1.0);
else:
   print "Loop 1 is finished"


for i in range(0, xGridSize + 1):
   x[i] = -1.0 + i*h;
   for j in range(1, n + 1):
      if (x[i] <= xn[j]):
         break
   fx[i] = (x[i]-xn[j-1])/(xn[j]-xn[j-1])*fxn[j] + (x[i]-xn[j])/(xn[j-1]-xn[j])*fxn[j-1];
else:
   print "Loop 2 is finished"

for i in range(0, xGridSize + 1):
   f = 1.0/(25.0*(x[i]**2)+1.0);
   EN2 = EN2 + ((fx[i]-f)/f)**2;
else:
   EN2 = sqrt(EN2) / n;
   print "EN2 = ", EN2

figure(1)
plot(x, fx, 'g', label='n = %i' %(n))
plot(xn, fxn, 'go' )

##########################

n = 12;
x = zeros(xGridSize+1);
xn = zeros(n+1);
fx = zeros(xGridSize+1);
fxn = zeros(n+1);
hn = (1.0+1.0)/n;
h = (1.0+1.0)/xGridSize;
EN2 = 0.0;

for i in range(0, n + 1):
   xn[i] = -1.0 + i*hn;
   fxn[i]= 1.0/(25.0*(xn[i]**2)+1.0);
else:
   print "Loop 1 is finished"


for i in range(0, xGridSize + 1):
   x[i] = -1.0 + i*h;
   for j in range(1, n + 1):
      if (x[i] <= xn[j]):
         break
   fx[i] = (x[i]-xn[j-1])/(xn[j]-xn[j-1])*fxn[j] + (x[i]-xn[j])/(xn[j-1]-xn[j])*fxn[j-1];
else:
   print "Loop 2 is finished"

for i in range(0, xGridSize + 1):
   f = 1.0/(25.0*(x[i]**2)+1.0);
   EN2 = EN2 + ((fx[i]-f)/f)**2;
else:
   EN2 = sqrt(EN2) / n;
   print "EN2 = ", EN2

figure(1)
plot(x, fx, 'm', label='n = %i' %(n))
plot(xn, fxn, 'mo' )


#####################

n = 50;
x = zeros(n+1);
fx = zeros(n+1);
h = (1.0+1.0)/n;

for i in range(0, n + 1):
   x[i] = -1.0 + i*h;
   fx[i]= 1.0/(25.0*(x[i]**2)+1.0);
else:
   print "Loop 1 is finished"


figure(1)
plot(x, fx, 'k', label='Analytical')
xlabel(r'\textbf{x}')
ylabel(r'\textit{f(x)}',fontsize=16)


legend()
savefig('plot6.pdf')


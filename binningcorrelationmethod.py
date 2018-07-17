from sys import argv
from MMTK.Minimization import SteepestDescentMinimizer
from Scientific import N
from numpy import zeros, array, sqrt, mean, var, log, floor
from MMTK import Units

rmean=0.0
rvariance=0.0

inputfile=open(argv[1],"r")
#label=argv[2]
#dt=3.0*Units.fs
#nsteps=int(100.0*Units.ps/dt)
#ntau=int(5.0*Units.ps/dt)

nsteps=int(argv[2])
nstepsbins=int(pow(2,int(log(nsteps)/log(2))))
nbins=int(floor(log(nsteps)/log(2)))-5
corrfile=open("corr-"+str(argv[1]),"w")
#nsteps=nstepsbins

r=zeros(nsteps,float)
for counter in range(nsteps):
    r[counter]=inputfile.readline()

print "Binning Level "+str(0)+": "+str(mean(r))+" "+str(sqrt(var(r)/(nsteps-1)))+" "+str(nsteps)
corrfile.write(str(0)+" "+str(sqrt(var(r)/(nsteps-1)))+"\n")

for l in range(1,nbins):
    nsteps=nsteps/2
    rnew=zeros(nsteps,float)
    for i in range(nsteps):
        rnew[i]=0.5*(r[2*i]+r[2*i+1])
    print "Binning Level "+str(l)+": "+str(mean(rnew))+" "+str(sqrt(var(rnew)/(nsteps-1)))+" "+str(nsteps)
    corrfile.write(str(l)+" "+str(sqrt(var(rnew)/(nsteps-1)))+"\n")
    r=rnew

corrfile.close()

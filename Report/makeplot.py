from sys import argv
import numpy as np
import matplotlib.pyplot as plt

# example data
filename10 = argv[1]
filename20 = argv[2]
filename30 = argv[3]
filename50 = argv[4]

Avg_Rot = np.zeros((4),float)
Var_Rot = np.zeros((4),float)

wd_T10= open(filename10,"r")
num_lines = sum(1 for line in open(filename10))
wd_T10_rot=np.zeros((num_lines), float)
for i in range(num_lines):
    column = wd_T10.readline().split()
    wd_T10_rot[i]= float(column[4])
wd_T10.close()
Avg_Rot[0] = np.mean(wd_T10_rot)
Var_Rot[0] = np.var(wd_T10_rot)

wd_T20= open(filename20,"r")
num_lines = sum(1 for line in open(filename20))
wd_T20_rot=np.zeros((num_lines), float)
for i in range(num_lines):
    column = wd_T20.readline().split()
    wd_T20_rot[i]= float(column[4])
wd_T20.close()
Avg_Rot[1] = np.mean(wd_T20_rot)
Var_Rot[1] = np.var(wd_T20_rot)



wd_T30= open(filename30,"r")
num_lines = sum(1 for line in open(filename30))
wd_T30_rot=np.zeros((num_lines), float)
for i in range(num_lines):
    column = wd_T30.readline().split()
    wd_T30_rot[i]= float(column[4])
wd_T30.close()
Avg_Rot[2] = np.mean(wd_T30_rot)
Var_Rot[2] = np.var(wd_T30_rot)


wd_T50= open(filename50,"r")
num_lines = sum(1 for line in open(filename50))
wd_T50_rot=np.zeros((num_lines), float)
for i in range(num_lines):
    column = wd_T50.readline().split()
    wd_T50_rot[i]= float(column[4])
wd_T50.close()
Avg_Rot[3] = np.mean(wd_T50_rot)
Var_Rot[3] = np.var(wd_T50_rot)

print np.sqrt( Var_Rot)
print Avg_Rot

fig, ax = plt.subplots()
Tem=[10,20,30,50]
Eigen_cal=[4.0658,22.8049,39.7044,70.2394]
plt.xlim(0,60)
plt.plot(Tem, Eigen_cal,"ro",markersize=4)
ax.errorbar(Tem, Avg_Rot,yerr=Var_Rot)
plt.xlabel('Temperature(Kelvin)')
plt.ylabel('Energy')
plt.show()

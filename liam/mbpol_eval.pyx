cimport numpy as np
import numpy as np
import ctypes

cdef extern from "mbpol.h" namespace "x2o":
	cdef cppclass mbpol:
		mbpol()
		double operator()(int, double*, double*)
		double operator()(int, double*)

def mbpolEvalEnergy():
	print 'blah'

def mbpolEvalEnergyAndGrad(configuration, numberOfWaters):

	cdef int nw = numberOfWaters
	cdef double E 

	grd = np.zeros(len(configuration))

	def callMbpol(np.ndarray[double, ndim=1, mode="c"] coord not None, 
					np.ndarray[double, ndim=1, mode="c"] grad not None, 
					int value):
		E = mbpol()(nw, &coord[0], &grad[0])
		return E, grad

	E, grd = callMbpol(configuration,grd,nw)

	if E < -20.7662037:
		Omass = 15.9994
		Hmass = 1.00794
		H2Omass = Omass + 2*Hmass
		com1 = [configuration[0]*Omass/H2Omass + configuration[3]*Hmass/H2Omass  + configuration[6]*Hmass/H2Omass, configuration[1]*Omass/H2Omass  + configuration[4]*Hmass/H2Omass  + configuration[7]*Hmass/H2Omass,configuration[2]*Omass/H2Omass + configuration[5]*Hmass/H2Omass + configuration[8]*Hmass/H2Omass ]
		com2 = [configuration[9]*Omass/H2Omass + configuration[12]*Hmass/H2Omass + configuration[15]*Hmass/H2Omass,configuration[10]*Omass/H2Omass + configuration[13]*Hmass/H2Omass + configuration[16]*Hmass/H2Omass,configuration[11]*Omass/H2Omass + configuration[14]*Hmass/H2Omass + configuration[17]*Hmass/H2Omass ]
		dist = np.sqrt((com1[0]-com2[0])**2 + (com1[1]-com2[1])**2 + (com1[2]-com2[2])**2)

		error_file = open('errors','a')
		#error_file.write(str(E) + ' ' + str(dist) + ' ' + str(configuration) + '\n')
		error_file.write(str(E) + ' ' + str(dist)  + '\n')
		error_file.close()

	return E, grd

	# np.ndarray[double, ndim=2, mode="c"] test = configuration

	# print type(configuration)

	# E = mbpol()(nw, <double*> configuration.data)


	# E = mbpol()(nw,coord)
	# E = mbpol()(nw,&configuration[0])

	# print E

	# print type(nw)

	# print configuration

	# nw = 2
	# grd = np.zeros(18)

	# print type(configuration), type(configuration[0])

	# arr = (ctypes.c_double * len(configuration))(*configuration)

	# print arr, type(arr)

	# Ep = mbpol()(nw,arr,grd)

#cdef mbpol* test = new mbpol()

#nw = 1
#coord = np.array([10.,10.1,10.2])
#grd =  np.array([0.,0.,0.])

# cdef int nw = 2
# cdef double[18] coord
# cdef double[18] grd

# coord[0] = -1.516074336e+00
# coord[1] = -2.023167650e-01
# coord[2] = 1.454672917e+00
# coord[3] = -6.218989773e-01
# coord[4] = -6.009430735e-01 
# coord[5] = 1.572437625e+00
# coord[6] = -2.017613812e+00
# coord[7] = -4.190350349e-01
# coord[8] = 2.239642849e+00
# coord[9] = -1.763651687e+00 
# coord[10] = -3.816594649e-01 
# coord[11] = -1.300353949e+00 
# coord[12] = -1.903851736e+00
# coord[13] = -4.935677617e-01
# coord[14] = -3.457810126e-01
# coord[15] = -2.527904158e+00
# coord[16] = -7.613550077e-01
# coord[17] = -1.733803676e+00
# grd[0]=0.
# grd[1]=0.
# grd[2]=0.
# grd[3]=0.
# grd[4]=0.
# grd[5]=0.
# grd[6]=0.
# grd[7]=0.
# grd[8]=0.
# grd[9]=0.
# grd[10]=0.
# grd[11]=0.
# grd[12]=0.
# grd[13]=0.
# grd[14]=0.
# grd[15]=0.
# grd[16]=0.
# grd[17]=0.

# #cdef:
# #	double *ccoord = <double *> cord.data
# #	double *cgrd = <double *> grd.data

# cdef double blah
# blah = mbpol()(nw,coord,grd)
# print blah
# print grd[0]
# #blah = mbpol()(nw,coord,grd)
# #print blah

# gfd = np.zeros(18)


# cdef int i
# cdef double[18] tmp
# cdef double x_orig, Ep, Em
# cdef double eps = 1.0e-4

# for i in range(18):
# 	x_orig = coord[i]
# 	coord[i] = x_orig + eps
# 	Ep = mbpol()(nw,coord,tmp)
# 	coord[i] = x_orig - eps
# 	Em = mbpol()(nw,coord,tmp)
# 	gfd[i] = 0.5*(Ep-Em)/eps
# 	coord[i] = x_orig
	
# 	#print grd[i]-gfd

# print grd[0],grd[1],grd[2]
# print grd[3],grd[4],grd[5]
# print grd[6],grd[7],grd[8]
# print grd[9],grd[10],grd[11]
# print grd[12],grd[13],grd[14]
# print grd[15],grd[16],grd[17]

# print gfd

# #print 'hello'

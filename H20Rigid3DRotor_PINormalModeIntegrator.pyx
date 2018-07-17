# This module implements path integral MD integrator using normal mode coordinates
#
# Written by Konrad Hinsen
#

#cython: boundscheck=False, wraparound=False, cdivision=True

"""
Path integral MD integrator using normal-mode coordinates
"""

__docformat__ = 'restructuredtext'

from cpython.pycapsule cimport PyCapsule_GetPointer, PyCapsule_New

from libc.stdint cimport int32_t
import numpy as N
cimport numpy as N

from MMTK import Units, ParticleProperties, Features, Environment, Vector
import MMTK.PIIntegratorSupport
cimport MMTK.PIIntegratorSupport
import numbers

from MMTK.forcefield cimport energy_data
cimport MMTK.mtrand

include 'MMTK/trajectory.pxi'

cdef extern from "fftw3.h":
    ctypedef struct fftw_complex
    ctypedef void *fftw_plan
    cdef int FFTW_FORWARD, FFTW_BACKWARD, FFTW_ESTIMATE
    cdef void fftw_execute(fftw_plan p)
    cdef fftw_plan fftw_plan_dft_1d(int n, fftw_complex *data_in, fftw_complex *data_out,
                                    int sign, int flags)
    cdef void fftw_destroy_plan(fftw_plan p)

cdef extern from "stdlib.h":
    cdef double fabs(double)
    cdef double sqrt(double)
    cdef double sin(double)
    cdef double cos(double)
    cdef double exp(double)
    cdef double M_PI

cdef extern from "time.h":
    ctypedef unsigned long clock_t
    cdef clock_t clock()
    cdef enum:
        CLOCKS_PER_SEC

cdef double hbar = Units.hbar
cdef double k_B = Units.k_B

cdef bytes PLAN_CAPSULE_NAME = b'plan_capsule'

cdef void plan_capsule_destructor(object cap):
    fftw_destroy_plan(PyCapsule_GetPointer(cap, PLAN_CAPSULE_NAME))

#
# Velocity Verlet integrator in normal-mode coordinates
#
cdef class Rigid3DRotor_PINormalModeIntegrator(MMTK.PIIntegratorSupport.PIIntegrator):

    """
    Molecular dynamics integrator for path integral systems using
    normal-mode coordinates.

    The integration is started by calling the integrator object.
    All the keyword options (see documentation of __init__) can be
    specified either when creating the integrator or when calling it.

    The following data categories and variables are available for
    output:

     - category "time": time

     - category "configuration": configuration and box size (for
       periodic universes)

     - category "velocities": atomic velocities

     - category "gradients": energy gradients for each atom

     - category "energy": potential and kinetic energy, plus
       extended-system energy terms if a thermostat and/or barostat
       are used

     - category "thermodynamic": temperature

     - category "auxiliary": primitive and virial quantum energy estimators

    """

    cdef N.ndarray workspace1, workspace2
    cdef double *workspace_ptr_1
    cdef double *workspace_ptr_2
    cdef dict plans
    cdef N.ndarray denrho, denerot, denesq
    cdef double rotmove
    cdef int rotstepskip

    def __init__(self, universe, **options):
        """
        :param universe: the universe on which the integrator acts
        :type universe: MMTK.Universe
        :keyword steps: the number of integration steps (default is 100)
        :type steps: int
        :keyword delta_t: the time step (default is 1 fs)
        :type delta_t: float
        :keyword actions: a list of actions to be executed periodically
                          (default is none)
        :type actions: list
        :keyword threads: the number of threads to use in energy evaluation
                          (default set by MMTK_ENERGY_THREADS)
        :type threads: int
        :keyword background: if True, the integration is executed as a
                             separate thread (default: False)
        :type background: bool
        """
        MMTK.PIIntegratorSupport.PIIntegrator.__init__(
            self, universe, options, "Path integral normal-mode integrator")
        # Supported features: PathIntegrals
        self.features = [Features.PathIntegralsFeature]

    default_options = {'first_step': 0, 'steps': 100, 'delta_t': 1.*Units.fs,
                       'background': False, 'threads': None,
                       'frozen_subspace': None, 'actions': []}

    available_data = ['time', 'configuration', 'velocities', 'gradients',
                      'energy', 'thermodynamic', 'auxiliary']

    restart_data = ['configuration', 'velocities', 'energy']

    # The implementation of the equations of motion follows the article
    #   Ceriotti et al., J. Chem. Phys. 133, 124104 (2010)
    # with the following differences:
    # 1) All the normal mode coordinates are larger by a factor sqrt(nbeads),
    #    and the non-real ones (k != 0, k != n/2) are additionally smaller by
    #    sqrt(2).
    # 2) The spring energy is smaller by a factor of nbeads to take
    #    into account the factor nbeads in Eq. (3) of the paper cited above.
    #    The potential energy of the system is also smaller by a factor of
    #    nbeads compared to the notation in this paper.
    # 3) Velocities are used instead of momenta in the integrator.
    # 4) Eq. (18) is also used for odd n, ignoring the k = n/2 case.

    cdef cartesianToNormalMode(self, N.ndarray[double, ndim=2] x, N.ndarray[double, ndim=2] nmc,
                               Py_ssize_t bead_index, int32_t nb):
        cdef double *w1 = self.workspace_ptr_1
        cdef double *w2 = self.workspace_ptr_2
        cdef fftw_plan p
        cdef Py_ssize_t i, j
        if nb == 1:
            for i in range(3):
                nmc[i, bead_index] = x[bead_index, i]
        else:
            try:
                p = PyCapsule_GetPointer(self.plans[(FFTW_FORWARD, nb)], PLAN_CAPSULE_NAME)
            except KeyError:
                p = fftw_plan_dft_1d(nb, <fftw_complex *>w1, <fftw_complex *>w2,
                                     FFTW_FORWARD, FFTW_ESTIMATE)
                self.plans[(FFTW_FORWARD, nb)] = \
                        PyCapsule_New(p, PLAN_CAPSULE_NAME, plan_capsule_destructor)
            for i in range(3):
                for j in range(nb):
                    w1[2*j] = x[bead_index+j, i]
                    w1[2*j+1] = 0.
                fftw_execute(p)
                nmc[i, bead_index+0] = w2[0]
                for j in range(1, (nb+1)/2):
                    nmc[i, bead_index+j] = w2[2*j]
                    nmc[i, bead_index+nb-j] = w2[2*j+1]
                if nb % 2 == 0:
                    nmc[i, bead_index+nb/2] = w2[nb]

    cdef normalModeToCartesian(self, N.ndarray[double, ndim=2] x, N.ndarray[double, ndim=2] nmc,
                               Py_ssize_t bead_index, int32_t nb):
        cdef double *w1 = self.workspace_ptr_1
        cdef double *w2 = self.workspace_ptr_2
        cdef fftw_plan p
        cdef Py_ssize_t i, j
        if nb == 1:
            for i in range(3):
                x[bead_index, i] = nmc[i, bead_index]
        else:
            try:
                p = PyCapsule_GetPointer(self.plans[(FFTW_BACKWARD, nb)], PLAN_CAPSULE_NAME)
            except KeyError:
                p = fftw_plan_dft_1d(nb, <fftw_complex *>w1, <fftw_complex *>w2,
                                     FFTW_BACKWARD, FFTW_ESTIMATE)
                self.plans[(FFTW_BACKWARD, nb)] = \
                        PyCapsule_New(p, PLAN_CAPSULE_NAME, plan_capsule_destructor)
            for i in range(3):
                w1[0] = nmc[i, bead_index+0]
                w1[1] = 0.
                for j in range(1, (nb+1)/2):
                    w1[2*j] = nmc[i, bead_index+j]
                    w1[2*j+1] = nmc[i, bead_index+nb-j]
                    w1[2*nb-2*j] = w1[2*j]
                    w1[2*nb-2*j+1] = -w1[2*j+1]
                if nb % 2 == 0:
                    w1[nb] = nmc[i, bead_index+nb/2]
                    w1[nb+1] = 0.
                fftw_execute(p)
                for j in range(nb):
                    x[bead_index+j, i] = w2[2*j]/nb

    cdef void propagateOscillators(self, N.ndarray[double, ndim=2] nmc,
                                   N.ndarray[double, ndim=2] nmv,
                                   Py_ssize_t bead_index, int32_t nb, double beta, double dt):
        cdef double omega_n = nb/(beta*hbar)
        cdef double omega_k, omega_k_dt, s, c
        cdef double temp
        cdef Py_ssize_t i, k
        for i in range(3):
            nmc[i, bead_index] += dt*nmv[i, bead_index]
            for k in range(1, nb):
                omega_k = 2.*omega_n*sin(k*M_PI/nb)
                omega_k_dt = omega_k*dt
                s = sin(omega_k_dt)
                c = cos(omega_k_dt)
                temp = c*nmv[i, bead_index+k]-omega_k*s*nmc[i, bead_index+k]
                nmc[i, bead_index+k] = s*nmv[i, bead_index+k]/omega_k + c*nmc[i, bead_index+k]
                nmv[i, bead_index+k] = temp

    cdef double springEnergyNormalModes(self, N.ndarray[double, ndim=2] nmc,
                                        N.ndarray[double, ndim=1] m,
                                        N.ndarray[N.int32_t, ndim=2] bd,
                                        double beta):
        cdef Py_ssize_t i, j, k
        cdef int32_t nb
        cdef double sumsq
        cdef double omega_n, omega_k
        cdef double e = 0.
        for i in range(nmc.shape[1]):
            if bd[i, 0] == 0:
                nb = bd[i, 1]
                omega_n = nb/(beta*hbar)
                # Start at j=1 because the contribution from the centroid is zero
                for j in range(1, nb):
                    omega_k = 2.*omega_n*sin(j*M_PI/nb)
                    sumsq = 0.
                    for k in range(3):
                        sumsq += nmc[k, i+j]*nmc[k, i+j]
                    # j=nb/2 corresponds to the real-valued coordinate at
                    # the maximal frequency.
                    if nb % 2 == 0 and j == nb/2:
                        sumsq *= 0.5
                    e += m[i]*sumsq*omega_k*omega_k/nb
        return e

    cdef void applyThermostat(self, N.ndarray[double, ndim=2] v, N.ndarray[double, ndim=2] nmv,
                              N.ndarray[double, ndim=1] m, N.ndarray[N.int32_t, ndim=2] bd,
                              double dt, double beta):
        pass

    cdef void atomtocm(self, N.ndarray[double, ndim=2] x, N.ndarray[double, ndim=2] v,
                       N.ndarray[double, ndim=2] g, N.ndarray[double, ndim=1] m,
                       N.ndarray[double, ndim=2] xcm, N.ndarray[double, ndim=2] vcm,
                       N.ndarray[double, ndim=2] gcm, N.ndarray[double, ndim=1] mcm,
                       N.ndarray[N.int32_t, ndim=2] bdcm, int Nmol):

         cdef int tot_atoms,i,j,k,z,natomspmol,nbeadspmol, atom_index
         tot_atoms=0
         for i in range (Nmol):
            natomspmol=self.universe.objectList()[i].numberOfAtoms()
            # nbeadspmol is the number of beads we want our molecule COM to have. 
            # Therefore is the number of beads each atom has in the molecule.
            nbeadspmol=self.universe.objectList()[i].numberOfPoints()/natomspmol          
            
            for z in range (nbeadspmol):
                bdcm[i*nbeadspmol+z,0]=N.int32(z)
                if bdcm[i*nbeadspmol+z,0] == N.int32(0):
                    bdcm[i*nbeadspmol+z,1]=N.int32(nbeadspmol)
                mcm[i*nbeadspmol+z]=self.universe.objectList()[i].mass()/nbeadspmol
                
                
                for k in range(3):
                    xcm[i*nbeadspmol+z,k]=0.0
                    vcm[i*nbeadspmol+z,k]=0.0
                    gcm[i*nbeadspmol+z,k]=0.0
                    for j in range(natomspmol):
                        atom_index=tot_atoms+j
                        xcm[i*nbeadspmol+z,k]+=m[atom_index*nbeadspmol+z]*x[atom_index*nbeadspmol+z,k]/mcm[i*nbeadspmol+z]
                        vcm[i*nbeadspmol+z,k]+=m[atom_index*nbeadspmol+z]*v[atom_index*nbeadspmol+z,k]/mcm[i*nbeadspmol+z]
                        gcm[i*nbeadspmol+z,k]+=g[atom_index*nbeadspmol+z,k]

            tot_atoms+=natomspmol


    cdef void cmtoatom(self, N.ndarray[double, ndim=2] x, N.ndarray[double, ndim=2] v,
                       N.ndarray[double, ndim=2] g, N.ndarray[double, ndim=1] m,
                       N.ndarray[double, ndim=2] xcm, N.ndarray[double, ndim=2] vcm,
                       N.ndarray[double, ndim=2] gcm, N.ndarray[double, ndim=1] mcm,
                       int Nmol):

         #xcom is ORIGINAL center of mass!
         cdef N.ndarray[double,ndim=1] xcom
         cdef int tot_atoms,i,j,k,z,natomspmol,nbeadspmol, atom_index

         xcom=N.zeros((3,),N.float)

         tot_atoms=0

         for i in range (Nmol):
            natomspmol=self.universe.objectList()[i].numberOfAtoms()
            nbeadspmol=self.universe.objectList()[i].numberOfPoints()/natomspmol
            for z in range (nbeadspmol):
                for k in range(3):
                    xcom[k]=0.
                    for j in range(natomspmol):
                        atom_index=tot_atoms+j
                        xcom[k]+=m[atom_index*nbeadspmol+z]*x[atom_index*nbeadspmol+z,k]/mcm[i*nbeadspmol+z]
                for k in range(3):
                    for j in range(natomspmol):
                        atom_index=tot_atoms+j
                        x[atom_index*nbeadspmol+z,k]=x[atom_index*nbeadspmol+z,k]-xcom[k]+xcm[i*nbeadspmol+z,k]
                        g[atom_index*nbeadspmol+z,k]=gcm[i*nbeadspmol+z,k]*m[atom_index*nbeadspmol+z]/mcm[i*nbeadspmol+z]
                        v[atom_index*nbeadspmol+z,k]=vcm[i*nbeadspmol+z,k]

            tot_atoms+=natomspmol

##### Lori Modified the eulertocart function (General for C2V triatomic molecule (FROM WFF TO SFF
    cdef void eulertocart(self, int bindex, int molind, N.ndarray[double,ndim=2] rcm_stand, N.ndarray[double, ndim=2] x, N.ndarray[double,ndim=1] eulerangles, N.ndarray[double,ndim=2] xcm):
        cdef N.ndarray[double,ndim=2] R
        cdef N.ndarray[double,ndim=1] rcm_LN, rcm_CN, rcm_RN
        natomspmol=self.universe.objectList()[molind].numberOfAtoms()
        nbeadspmol=self.universe.objectList()[molind].numberOfPoints()/natomspmol
        R = N.zeros((3, 3), N.float)
        rcm_LN = N.zeros((3), N.float)
        rcm_CN = N.zeros((3), N.float)
        rcm_RN = N.zeros((3), N.float)
        # As eulertocart() used in MC-Rot Step, it has already looped over molnum, beadnum, and we seperate the atom index, so in this function, we don't need any loop
        atomL=0
        atomR=1
        atomC=2
        aindex = molind*natomspmol
        # For every R, xcm and rcm, the index is based on the molecule, not atom
        R = self.matpre(eulerangles)
        for m in range(3):
            for n in range(3):
                # From WFF to SFF CM coords
                rcm_LN[m]+=R[m][n]*rcm_stand[0][n]
                rcm_RN[m]+=R[m][n]*rcm_stand[1][n]
                rcm_CN[m]+=R[m][n]*rcm_stand[2][n]
        # After Rotation, the xcm doesn't change
        x[(aindex+atomL)*nbeadspmol+bindex]=xcm[molind*nbeadspmol+bindex]+rcm_LN
        x[(aindex+atomR)*nbeadspmol+bindex]=xcm[molind*nbeadspmol+bindex]+rcm_RN
        x[(aindex+atomC)*nbeadspmol+bindex]=xcm[molind*nbeadspmol+bindex]+rcm_CN

#print "rcm_LN= ", N.linalg.norm(rcm_LN), "rcm_RN= ", N.linalg.norm(rcm_RN), "rcm_CN= ", N.linalg.norm(rcm_CN)
#print "det(R)= ", N.linalg.det(R)
############################
##### Lori ADDS THE ASYMMETRIC TOP DENSITY GENERATOR (BASED ON TOBY's CODES
##### Judge the value of angles

    def within(self, double value):
        if (value > 1.0):
            value = 1.0
        elif (value < -1.0):
            value = -1.0

        return value

##############################

##### Generate Rotational Matrix (FROM MFF to SFF CM
    cdef N.ndarray[double,ndim=2] matpre(self, N.ndarray[double,ndim=1] eulang):
        cdef N.ndarray[double, ndim=2] rotmat
        rotmat = N.zeros((3,3), N.float)

        phi = eulang[0]
        theta = eulang[1]
        chi = eulang[2]

        cp = N.cos(phi)
        sp = N.sin(phi)
        ct = N.cos(theta)
        st = N.sin(theta)
        ck = N.cos(chi)
        sk = N.sin(chi)

        rotmat[0][0]=cp*ct*ck-sp*sk
        rotmat[0][1]=-cp*ct*sk-sp*ck
        rotmat[0][2]=cp*st
        rotmat[1][0]=sp*ct*ck+cp*sk
        rotmat[1][1]=-sp*ct*sk+cp*ck
        rotmat[1][2]=sp*st
        rotmat[2][0]=-st*ck
        rotmat[2][1]=st*sk
        rotmat[2][2]=ct

        return rotmat
        
##############################
##### GET R-1 from R Transpose
    def rotinv(self, N.ndarray[double,ndim=2] rotmat):
        cdef N.ndarray[double, ndim=2] roti
        roti = N.zeros((3,3), N.float)
        roti = N.transpose(rotmat)
        return roti

##############################       
##### Getting the delta Euler angles between two sets of Euler angles (in WFF)
    def deleul(self, N.ndarray[double,ndim=1] eulan1, N.ndarray[double,ndim=1] eulan2):
        cdef N.ndarray[double, ndim=2] rotmar,  rotma1,  rotma2, rotma1t
        cdef N.ndarray[double, ndim=1] eulrel   
        cdef double small
        rotma1 = N.zeros((3,3),N.float)
        rotmar = N.zeros((3,3),N.float)
        rotma2 = N.zeros((3,3),N.float)
        rotma1t = N.zeros((3,3), N.float)
        eulrel = N.zeros((3), N.float)
        small = 1.0e-10
        
        rotma2 = self.matpre(eulan2) # The second Rotation Matrix from Second Euler angles
        rotma1 = self.matpre(eulan1) # The first Rotation Matrix from First Euler angles
        
        for i in range(3):
            for j in range(3):
                rotmar[i][j] = 0.0 
                for k in range(3):
                    rotmar[i][j] = rotmar[i][j] + rotma1[k][i]*rotma2[k][j]
        
        cost = rotmar[2][2]
        cost = self.within(cost)
        thetar = N.arccos(cost)
        sint = N.sin(thetar)
        
        if (N.abs(1.0-cost) < small):
            phir = 0.0
            cchi = rotmar[0][0]
            schi = rotmar[1][0]
            cchi = self.within(cchi)
            schi = self.within(schi)
            if (schi > 0.0):
                chir = N.arccos(cchi)
            else:        
                chir = 2.0*N.pi - N.arccos(cchi)
        
        elif (N.abs(1.0+cost) < small):
            phir=0.0
            cchi = rotmar[1][1]
            schi = rotmar[0][1]
            cchi = self.within(cchi)
            schi = self.within(schi)
            if (schi > 0.0):
                chir = N.arccos(cchi)
            else:
                chir = 2.0*N.pi - N.arccos(cchi)
        else:
            cphi = rotmar[0][2]/sint
            sphi = rotmar[1][2]/sint
            cchi = -rotmar[2][0]/sint
            schi = rotmar[2][1]/sint
            cphi = self.within(cphi)
            sphi = self.within(sphi)
            cchi = self.within(cchi)
            schi = self.within(schi)
            if (sphi > 0.0):
                phir = N.arccos(cphi)
            else:
                phir = 2.0*N.pi - N.arccos(cphi)
            if (schi > 0.0):
                chir = N.arccos(cchi)
            else:
                chir = 2.0*N.pi-N.arccos(cchi)
                 
        if (phir<0.0):
            phir = 2.0*N.pi+phir
        if (chir<0.0):
            chir = 2.0*N.pi+chir
            
        phir = phir % (2.0*N.pi)
        chir = chir % (2.0*N.pi)
         
        if (cost>1.0):
            cost = 2.0 - cost
        elif (cost<-1.0):
            cost = -2.0 - cost
        
        thetar = N.arccos(cost)


        eulrel[0]=phir
        eulrel[1]=thetar
        eulrel[2]=chir
        return eulrel
###########################
#### Rot Density Generator
    def rotden(self, N.ndarray[double,ndim=1] eulan1, N.ndarray[double,ndim=1] eulan2, N.ndarray[double,ndim=1] rhomat, N.ndarray[double,ndim=1] erotmat, N.ndarray[double,ndim=1] esqmat, int istop):         
        cdef double wntok
        cdef N.ndarray[double, ndim=1] eulrel, rotinfo
        eulrel = N.zeros((3), N.float)
        rotinfo = N.zeros((3), N.float)
        eulrel = self.deleul(eulan1, eulan2)
        
        # the angles of rho are based on degree
        phi=eulrel[0]*180.0/N.pi
        theta=eulrel[1]*180.0/N.pi
        chi=eulrel[2]*180.0/N.pi
        
        jstop=0
        
        cdef int indchi, indphi, indthe, index_p0, index_p1
        cdef double rho_0, erot_0, esq_0, delchi, delphi, delthe, delch2, delch3, delph2, delph3, delth2, delth3
        # Flooring the euler angles to get the index of the angles         
        indchi=int(chi)
        indphi=int(phi)
        indthe=int(theta)

        if (indchi > 360 or indchi < 0):
            print "index of chi is out of range", indchi, chi
            indchi = 0
            jstop = 1
        elif (indphi > 360 or indphi < 0):
            print "index of phi is out of range", indphi, phi
            indphi = 0
            jstop = 1
        elif (indthe > 180 or indthe < 0):
            print "index of theta is out of range", indthe, theta
            indthe = 0
            jstop = 1

        # This index is based on chi, for each fixed theta, we have 361 phi indexes, and for each fixed phi, we have 361 chi indexes
        index_p0=(indthe*361+indphi)*361+indchi # index for point 0
        rho_0=rhomat[index_p0]
        erot_0=erotmat[index_p0]
        esq_0=esqmat[index_p0]

        if (indchi != 360):
            index_p1=(indthe*361+indphi)*361+indchi+1 # next index is just the difference by 1 degree
            # CHI part for each properties
            delchi=rhomat[index_p1]-rho_0
            delch2=erotmat[index_p1]-erot_0
            delch3=esqmat[index_p1]-esq_0

        if (indphi != 360):
            index_p1=(indthe*361+indphi+1)*361+indchi
            # PHI part for each properties
            delphi=rhomat[index_p1]-rho_0
            delph2=erotmat[index_p1]-erot_0
            delph3=esqmat[index_p1]-esq_0

        if (indthe != 180):
            index_p1=((indthe+1)*361+indphi)*361+indchi
            # THETA part for each properties
            delthe=rhomat[index_p1]-rho_0
            delth2=erotmat[index_p1]-erot_0
            delth3=esqmat[index_p1]-esq_0

        # Linear Interpolation for all three angles
        rhoden=rho_0+delchi*(chi-N.float(indchi))+delphi*(phi-N.float(indphi))+delthe*(theta-N.float(indthe))
        erotden=erot_0+delch2*(chi-N.float(indchi))+delph2*(phi-N.float(indphi))+delth2*(theta-N.float(indthe))
        esqden=esq_0+delch3*(chi-N.float(indchi))+delph3*(phi-N.float(indphi))+delth3*(theta-N.float(indthe))

        if (jstop == 1):
            print eulan1[0], eulan1[1], eulan1[2], eulan2[0], eulan2[1], eulan2[2], eulrel[0], eulrel[1], eulrel[2], phi, theta, chi
            istop =1
            print "Large matrix test error"
            raise()
        
        rotinfo[0] = rhoden
        rotinfo[1] = erotden
        rotinfo[2] = esqden
        
        return rotinfo
#############################
###### Liam wrote the codes for generating the Euler angles from SFF CM coords (ONLY FOR WATER!
###### Lori modified it for general C2V triatomic molecule and integrated with MMTK
    def EulerGen(self, double MassL, double MassR, double MassC, N.ndarray[double,ndim=1] rcm_LNew, N.ndarray[double,ndim=1] rcm_RNew, N.ndarray[double,ndim=1] rcm_CNew):
        cdef N.ndarray[double,ndim=1] aX, aY, aZ, w_new, dummy_vec, eulerangles
        cdef N.ndarray[double,ndim=2] BFF, I_New, I_New_Diag, v_new, v_inv_new, RotDeux
        cdef double phi, theta, chi

        eulerangles = N.zeros((3), N.float)
        
        phi = 0.0
        theta = 0.0
        chi = 0.0

        # Hard code the standard BFF principle axes
        # Check if it is general for all C2V triatomic molecule? 
        BFF = N.zeros((3,3), N.float)
        BFF[0][0] = 1.
        BFF[2][1] = 1.
        BFF[1][2] = 1.
        
        # Generate the New Intertia of tensor
        I_New = N.zeros( (3,3) , N.float) 

        I_New[0][0] = MassC*(rcm_CNew[1]*rcm_CNew[1] + rcm_CNew[2]*rcm_CNew[2]) + MassR*(rcm_RNew[1]*rcm_RNew[1] + rcm_RNew[2]*rcm_RNew[2]) + MassL*(rcm_LNew[1]*rcm_LNew[1] + rcm_LNew[2]*rcm_LNew[2])
        I_New[1][1] = MassC*(rcm_CNew[0]*rcm_CNew[0] + rcm_CNew[2]*rcm_CNew[2]) + MassR*(rcm_RNew[0]*rcm_RNew[0] + rcm_RNew[2]*rcm_RNew[2]) + MassL*(rcm_LNew[0]*rcm_LNew[0] + rcm_LNew[2]*rcm_LNew[2])
        I_New[2][2] = MassC*(rcm_CNew[0]*rcm_CNew[0] + rcm_CNew[1]*rcm_CNew[1]) + MassR*(rcm_RNew[0]*rcm_RNew[0] + rcm_RNew[1]*rcm_RNew[1]) + MassL*(rcm_LNew[0]*rcm_LNew[0] + rcm_LNew[1]*rcm_LNew[1])

        I_New[0][1] = (-1)*(MassC*(rcm_CNew[0]*rcm_CNew[1]) + MassR*(rcm_RNew[0]*rcm_RNew[1]) + MassL*(rcm_LNew[0]*rcm_LNew[1]))
        I_New[1][0] = I_New[0][1]
        I_New[0][2] = (-1)*(MassC*(rcm_CNew[0]*rcm_CNew[2]) + MassR*(rcm_RNew[0]*rcm_RNew[2]) + MassL*(rcm_LNew[0]*rcm_LNew[2]))
        I_New[2][0] = I_New[0][2]
        I_New[1][2] = (-1)*(MassC*(rcm_CNew[1]*rcm_CNew[2]) + MassR*(rcm_RNew[1]*rcm_RNew[2]) + MassL*(rcm_LNew[1]*rcm_LNew[2]))
        I_New[2][1] = I_New[1][2]

        w_new, v_new = N.linalg.eig(I_New)

        #  get the inverse of the eigenvectors

        v_inv_new = N.linalg.inv(v_new)

        # Generate the Diaogonalize I_New

        I_New_Diag = N.zeros( (3,3) , N.float)

        for i in range(3):
                for n in range(3):
                        for j in range(3):
                                for l in range(3):
                                        I_New_Diag[i][n] += v_inv_new[i][j]*I_New[j][l]*v_new[l][n]

        # cleanup

        for i in range(3):
                for k in range(3):
                        if ( I_New_Diag[i][k] <= 1e-10 ):
                                I_New_Diag[i][k] = 0.0

        for i in range(3):
            I_New_Diag[i][i] = 1./I_New_Diag[i][i]

        # Define the Principal Axis from I_New

        aX = N.zeros(3,N.float)
        aY = N.zeros(3,N.float)
        aZ = N.zeros(3,N.float)

        dummy_vec = N.zeros(3,N.float)

        for i in range(3):
            if (I_New_Diag[1][1] >= I_New_Diag[0][0]):
                dummy = I_New_Diag[0][0]
                I_New_Diag[0][0] = I_New_Diag[1][1]
                w_new[0] = w_new[1]
                I_New_Diag[1][1] = dummy
                w_new[1] = 1./dummy
                for j in range(3):
                    dummy_vec[j] = v_new[j][0]
                    v_new[j][0] = v_new[j][1]
                    v_new[j][1] = dummy_vec[j]
            if (I_New_Diag[2][2] > I_New_Diag[1][1]):
                dummy = I_New_Diag[1][1]
                I_New_Diag[1][1] = I_New_Diag[2][2]
                w_new[1] = w_new[2]
                I_New_Diag[2][2] = dummy
                w_new[2] = 1./dummy
                for j in range(3):
                    dummy_vec[j] = v_new[j][1]
                    v_new[j][1] = v_new[j][2]
                    v_new[j][2] = dummy_vec[j]

        for j in range(3):
            aX[j] = v_new[j][0]
            aY[j] = v_new[j][1]
            aZ[j] = v_new[j][2]

        v_inv_new = N.linalg.inv(v_new)

        # Generate the Rotation Matrix RotDeux from eigenvectors

        RotDeux = N.zeros((3,3), N.float)

        for i in range(3):
            for j in range(3):
                for k in range(3):
                    RotDeux[i][k] += BFF[i][j]*v_inv_new[j][k]

        # Check the Axis Sign

        rcm_CP = N.zeros((3), N.float)
        rcm_RP = N.zeros((3), N.float) 
        rcm_LP = N.zeros((3), N.float)
                                                                        
        for i in range(3):
                for j in range(3):
                        rcm_CP[i] += RotDeux[i][j]*rcm_CNew[j]
                        rcm_RP[i] += RotDeux[i][j]*rcm_RNew[j]
                        rcm_LP[i] += RotDeux[i][j]*rcm_LNew[j]
        # Check the sign of oxygen z-axis
        if rcm_CP[2] <= 0:
            for i in range(3):
                RotDeux[2][i] = -1.*RotDeux[2][i]
        # Check the sign of right hydrogen x-axis
        if rcm_RP[0] <= 0:
            for j in range(3):
                RotDeux[0][j] = -1.*RotDeux[0][j]
        # Check the y-axis, whether it is right handed axis, use cross product
        for j in range(3):
            RotDeux[1][j] = N.cross(RotDeux[2],RotDeux[0])[j]

        # Derive the Euler angles from Rotational Matrix
        cost=RotDeux[2][2]
        cost=self.within(cost)
        theta = N.arccos(cost)
        sint=N.sin(theta)
        if (abs(1.0-cost)<1.0e-10):
            phi=0.0
            cchi=RotDeux[0][0]
            schi=RotDeux[1][0]
            cchi=self.within(cchi)
            schi=self.within(schi)
            if (schi > 0.0):
                chi=N.arccos(cchi)
            else:
                chi=2.0*N.pi-N.arccos(cchi)
        elif (abs(1.0+cost)<1.0e-10):
            phi=0.0
            cchi=RotDeux[1][1]
            schi=RotDeux[0][1]
            cchi=self.within(cchi)
            schi=self.within(schi)
            if (schi>0):
                chi=N.arccos(cchi)
            else:
                chi=2.0*N.pi-N.arccos(cchi)
        else:
            cphi=RotDeux[2][0]/sint
            sphi=RotDeux[2][1]/sint
            cchi=-1.0*RotDeux[0][2]/sint
            schi=RotDeux[1][2]/sint
            cphi=self.within(cphi)
            sphi=self.within(sphi)
            cchi=self.within(cchi)
            schi=self.within(schi)
            if (sphi>0):
                phi=N.arccos(cphi)
            else:
                phi=2.0*N.pi-N.arccos(cphi)
            if (schi>0):
                chi=N.arccos(cchi)
            else:
                chi=2.0*N.pi-N.arccos(cchi)
         
        if (phi<0.0):
            phi = 2.0*N.pi+phi
        if (chi<0.0):
            chi = 2.0*N.pi+chi
            
        phi = phi % (2.0*N.pi)
        chi = chi % (2.0*N.pi)
         
        if (cost>1.0):
            cost = 2.0 - cost
        elif (cost<-1.0):
            cost = -2.0 - cost
        
        theta = N.arccos(cost)

        eulerangles[0] = phi
        eulerangles[1] = theta
        eulerangles[2] = chi

        return eulerangles
##############################
    def energyCalculator(self, x):
        cdef energy_data energytemp
        energytemp.gradients = NULL
        energytemp.gradient_fn = NULL
        energytemp.force_constants = NULL
        energytemp.fc_fn = NULL
        self.calculateEnergies(x, &energytemp, 0)
        return energytemp.energy

    cdef start(self):
        cdef double acceptratio, rd, sint, pot_old, pot_new, dens_old, dens_new, indexp0val, indexp1val
        cdef int t0b, t1b, t2b, t0, t1, t2, atombead, indexp0, indexp1, indexp0n, indexp1n

        cdef N.ndarray[double, ndim=2] x, v, g, dv, nmc, nmv, xcm, vcm, gcm, rcm_L, rcm_C, rcm_R, rcm_stand
        cdef N.ndarray[double, ndim=1] m, mcm
        cdef N.ndarray[double, ndim=1] rcm_standL, rcm_standC, rcm_standR
        cdef N.ndarray[N.int32_t, ndim=2] bd, bdcm
        cdef N.ndarray[double, ndim=3] ss
        cdef energy_data energy
        cdef double time, delta_t, ke, ke_nm, se, beta, temperature
        cdef double qe_prim, qe_vir, qe_cvir, qe_rot, srot
        cdef int natoms, nbeads, nsteps, step, df, cdf, nb, Nmol, Ntruemol,rotbdcount,rotbdskip, istop
        cdef Py_ssize_t i, j, k

        cdef double propthe, propphi, propchi, propcth
        cdef double rho, erot, esq
        cdef int P
        cdef N.ndarray[double, ndim=2] MCAngles
        cdef N.ndarray[double, ndim=1] MCAngprop
        cdef N.ndarray[double, ndim=1] Eulan1, Eulan2
        cdef N.ndarray[double, ndim=2] xold
        cdef N.ndarray[double, ndim=1] rhomat, erotmat, esqmat
        cdef double rotstep,ndens
        cdef int rotskipstep, nrotsteps
        rhomat = self.denrho
        ndens = 1.0*len(rhomat)
        erotmat = self.denerot
        esqmat = self.denesq
        rotstep = self.rotmove
        rotskipstep = self.rotstepskip

        # Check if velocities have been initialized
        if self.universe.velocities() is None:
            raise ValueError("no velocities")

        # Gather state variables and parameters
        configuration = self.universe.configuration()
        velocities = self.universe.velocities()
        gradients = ParticleProperties.ParticleVector(self.universe)
        masses = self.universe.masses()
        delta_t = self.getOption('delta_t')
        nsteps = self.getOption('steps')
        natoms = self.universe.numberOfAtoms()
        nbeads = self.universe.numberOfPoints()
        bd = self.evaluator_object.global_data.get('bead_data')
        pi_environment = self.universe.environmentObjectList(Environment.PathIntegrals)[0]
        beta = pi_environment.beta

        # For efficiency, the Cython code works at the array
        # level rather than at the ParticleProperty level.
        x = configuration.array
        v = velocities.array
        g = gradients.array
        m = masses.array

	# MATT-Introduce X-COM variable, number of molecules Nmol
        acceptratio = 0.0
        P=nbeads/natoms
        Nmol = len(self.universe.objectList())
        nbeads_mol = N.int32(P*Nmol)
        xcm = N.zeros((nbeads_mol, 3), N.float)
        vcm = N.zeros((nbeads_mol, 3), N.float)
        gcm = N.zeros((nbeads_mol, 3), N.float)
        mcm = N.zeros(nbeads_mol, N.float)
        dv = N.zeros((nbeads_mol, 3), N.float)
        nmc = N.zeros((3, nbeads_mol), N.float)
        nmv = N.zeros((3, nbeads_mol), N.float)
        bdcm = N.zeros((nbeads_mol,2), N.int32)
        rcm_C = N.zeros((nbeads_mol,3), N.float)
        rcm_L = N.zeros((nbeads_mol,3), N.float)
        rcm_R = N.zeros((nbeads_mol,3), N.float)

        #ROTATIONAL VARIABLES
        nrotsteps=0
        Eulan1 = N.zeros(3, N.float)
        Eulan2 = N.zeros(3, N.float)
        MCAngprop = N.zeros((3), N.float)
        MCAngles = N.zeros((3,nbeads_mol), N.float)

        # Check if there is a frozen_subspace
        subspace = self.getOption('frozen_subspace')
        if subspace is None:
            ss = N.zeros((0, nbeads_mol, 3), N.float)
            df = 3*nbeads_mol
            cdf = 3*Nmol
        else:
            ss = subspace.getBasis().array
            df = 3*nbeads-ss.shape[0]
            cdf = self.centroidDegreesOfFreedom(subspace, bdcm)

        # Initialize the plan cache.
        self.plans = {}

        # Ask for energy gradients to be calculated and stored in
        # the array g. Force constants are not requested.
        energy.gradients = <void *>g
        energy.gradient_fn = NULL
        energy.force_constants = NULL
        energy.fc_fn = NULL

        # Declare the variables accessible to trajectory actions.
        self.declareTrajectoryVariable_double(
            &time, "time", "Time: %lf\n", time_unit_name, PyTrajectory_Time)
        self.declareTrajectoryVariable_array(
            v, "velocities", "Velocities:\n", velocity_unit_name,
            PyTrajectory_Velocities)
        self.declareTrajectoryVariable_array(
            g, "gradients", "Energy gradients:\n", energy_gradient_unit_name,
            PyTrajectory_Gradients)
        self.declareTrajectoryVariable_double(
            &energy.energy,"potential_energy", "Potential energy: %lf\n",
            energy_unit_name, PyTrajectory_Energy)
        self.declareTrajectoryVariable_double(
            &ke, "kinetic_energy", "Kinetic energy: %lf\n",
            energy_unit_name, PyTrajectory_Energy)
        self.declareTrajectoryVariable_double(
            &se, "spring_energy", "Spring energy: %lf\n",
            energy_unit_name, PyTrajectory_Energy)
        self.declareTrajectoryVariable_double(
            &temperature, "temperature", "Temperature: %lf\n",
            temperature_unit_name, PyTrajectory_Thermodynamic)
        self.declareTrajectoryVariable_double(
            &qe_prim, "quantum_energy_primitive",
            "Primitive quantum energy estimator: %lf\n",
            energy_unit_name, PyTrajectory_Auxiliary)
        self.declareTrajectoryVariable_double(
            &qe_vir, "quantum_energy_virial",
            "Virial quantum energy estimator: %lf\n",
            energy_unit_name, PyTrajectory_Auxiliary)
        self.declareTrajectoryVariable_double(
            &qe_cvir, "quantum_energy_centroid_virial",
            "Centroid virial quantum energy estimator: %lf\n",
            energy_unit_name, PyTrajectory_Auxiliary)
        self.declareTrajectoryVariable_double(
            &qe_rot, "quantum_energy_rotation",
            "Rotation quantum energy estimator: %lf\n",
            energy_unit_name, PyTrajectory_Auxiliary)
        self.initializeTrajectoryActions()

        # Acquire the write lock of the universe. This is necessary to
        # make sure that the integrator's modifications to positions
        # and velocities are synchronized with other threads that
        # attempt to use or modify these same values.
        #
        # Note that the write lock will be released temporarily
        # for trajectory actions. It will also be converted to
        # a read lock temporarily for energy evaluation. This
        # is taken care of automatically by the respective methods
        # of class EnergyBasedTrajectoryGenerator.
        self.acquireWriteLock()

        # Preparation: Calculate initial half-step accelerations
        # and run the trajectory actions on the initial state.
        self.foldCoordinatesIntoBox()

        Ntruemol=0
        for i in range(Nmol):
            if (self.universe.objectList()[i].numberOfAtoms()>1):
                Ntruemol+=1
                

        for i in range (Nmol):
            natomspmol=self.universe.objectList()[i].numberOfAtoms()
            # nbeadspmol is the number of beads we want our molecule COM to have.
            # Therefore is the number of beads each atom has in the molecule.
            nbeadspmol=self.universe.objectList()[i].numberOfPoints()/natomspmol

            for z in range (nbeadspmol):
                    bdcm[i*nbeadspmol+z,0]=N.int32(z)
                    if bdcm[i*nbeadspmol+z,0] == N.int32(0):
                        bdcm[i*nbeadspmol+z,1]=N.int32(nbeadspmol)
                    mcm[i*nbeadspmol+z]=self.universe.objectList()[i].mass()/nbeadspmol
										


        # Allocate workspace for Fourier transforms
        nb_max = bdcm[:, 1].max()
        self.workspace1 = N.zeros((2*nb_max,), N.float)
        self.workspace_ptr_1 = <double *>self.workspace1.data
        self.workspace2 = N.zeros((2*nb_max,), N.float)
        self.workspace_ptr_2 = <double *>self.workspace2.data


        #Calculate Energy and Fill Gradient Vector
        self.calculateEnergies(x, &energy, 0)
        self.atomtocm(x,v,g,m,xcm,vcm,gcm,mcm,bdcm,Nmol)


        ##########################################################
        ### CALCULATE ANGLES AND FILL MCAngles (ONLY FOR WATER ###
        ##########################################################
        atomL=0    # atom Left index arbitrary L
        atomR=1     # atom Center index
        atomC=2    # atom Right index arbitrary R

        rcm_standL = N.zeros((3),N.float)
        rcm_standR = N.zeros((3), N.float)
        rcm_standC = N.zeros((3), N.float)
 
        rcm_stand = N.zeros((3,3), N.float)

        rcm_standC[0] = 0.0
        rcm_standC[1] = 0.0
        rcm_standC[2] = 0.00655617
        rcm_standR[0] = 0.07569503
        rcm_standR[1] = 0.0
        rcm_standR[2] = -0.05203206
        rcm_standL[0] = -0.07569503
        rcm_standL[1] = 0.0
        rcm_standL[2] = -0.05203206

##### Create a standard MFF for C2V molecules (xz plane, and all three angles thould be zero
#        napmol=self.universe.objectList()[0].numberOfAtoms()
#        nbpmol=self.universe.objectList()[0].numberOfPoints()/natomspmol
#        # Just choose one set of bead to calculate the standard MFF
#        xC = x[atomC*nbpmol+nbpmol-1]
#        xL = x[atomL*nbpmol+nbpmol-1]
#        xR = x[atomR*nbpmol+nbpmol-1]
#        print 'xC= ', xC, 'xL', xL, 'xR', xR
#
#        rcm_standC = N.zeros((3),N.float)
#        rcm_standR = N.zeros((3), N.float)
#        rcm_standL = N.zeros((3), N.float)
# 
#        rcmCP = N.zeros((3),N.float)
#    
#        rcm_stand = N.zeros((3,3), N.float)
#        # initial MFF but not standard one
#        
#        rcmCP = xC-xcm[0*nbpmol+nbpmol-1]
#        
#        # angle between L-C-R
#        
#        aLCR = N.arccos(self.within(N.dot(xL-xC, xR-xC)/(N.linalg.norm(xL-xC)*N.linalg.norm(xR-xC))))
#        print 'lengLC=', N.linalg.norm(xL-xC)
#        print 'lengRC=', N.linalg.norm(xR-xC)
#        print 'aLCR=', aLCR
#
#        rcm_standC[2] = N.linalg.norm(rcmCP)
#        print 'lengthC=', rcm_standC[2]
#        rcm_standR[0] = N.linalg.norm(xR-xC)*N.sin(aLCR/2.)
#        rcm_standR[2] = -1.*(N.linalg.norm(xR-xC)*N.cos(aLCR/2.)-rcm_standC[2])
#        rcm_standL[0] = -rcm_standR[0]
#        rcm_standL[2] = rcm_standR[2]
        
        rcm_stand[0] = rcm_standL
        rcm_stand[1] = rcm_standR
        rcm_stand[2] = rcm_standC
       
        # Because we seperate the atoms index, there is no atom loop
        tot_atoms=0
        for i in range(Nmol):
            natomspmol=self.universe.objectList()[i].numberOfAtoms()              #number of atoms per mol
            nbeadspmol=self.universe.objectList()[i].numberOfPoints()/natomspmol   #number of beads per mol COM
            if (natomspmol>1):
               for p in range(nbeadspmol):
                    rcm_L[i*nbeadspmol+p]=x[(tot_atoms+atomL)*nbeadspmol+p]-xcm[i*nbeadspmol+p]
                    rcm_R[i*nbeadspmol+p]=x[(tot_atoms+atomR)*nbeadspmol+p]-xcm[i*nbeadspmol+p]
                    rcm_C[i*nbeadspmol+p]=x[(tot_atoms+atomC)*nbeadspmol+p]-xcm[i*nbeadspmol+p]
                    MCAngles[:,i*nbeadspmol+p]=self.EulerGen(m[(tot_atoms+atomL)*nbeadspmol+p], m[(tot_atoms+atomR)*nbeadspmol+p], m[(tot_atoms+atomC)*nbeadspmol+p], rcm_L[i*nbeadspmol+p], rcm_R[i*nbeadspmol+p], rcm_C[i*nbeadspmol+p])
            tot_atoms+=natomspmol
        
        self.freeze(vcm, ss)

        for i in range(nbeads_mol):
            if bdcm[i, 0] == 0:
                self.fixBeadPositions(xcm, i, bdcm[i, 1])
                self.cartesianToNormalMode(xcm, nmc, i, bdcm[i, 1])

        se = self.springEnergyNormalModes(nmc, mcm, bdcm, beta)
        qe_prim = energy.energy - se + 0.5*df/beta
        #qe_vir = energy.energy - 0.5*energy.virial
        #qe_cvir = energy.energy \
        #          - 0.5*self.centroidVirial(x, g, bd) \
        #          + 0.5*cdf/beta

        ke = 0.
        for i in range(nbeads_mol):
            for j in range(3):
                dv[i, j] = -0.5*delta_t*gcm[i, j]/mcm[i]
                ke += 0.5*mcm[i]*vcm[i, j]*vcm[i, j]
        temperature = 2.*ke/(df*k_B)

        #print "Before check FFT"
        #print g
        #print gcm

        # Check FFT
        if False:
            xcm_test = N.zeros((nbeads_mol, 3), N.float)
            vcm_test = N.zeros((nbeads_mol, 3), N.float)
            for i in range(nbeads_mol):
                if bdcm[i, 0] == 0:
                    self.cartesianToNormalMode(xcm, nmc, i, bdcm[i, 1])
                    self.normalModeToCartesian(xcm_test, nmc, i, bdcm[i, 1])
                    self.cartesianToNormalMode(vcm, nmv, i, bdcm[i, 1])
                    self.normalModeToCartesian(vcm_test, nmv, i, bdcm[i, 1])
            for i in range(nbeads_mol):
                for j in range(3):
                    assert fabs(xcm[i, j]-xcm_test[i, j]) < 1.e-7
                    assert fabs(vcm[i, j]-vcm_test[i, j]) < 1.e-7

        #timep1=clock()
        #time0 = (<double> (timep1 - timep0)) / CLOCKS_PER_SEC
        #print "Initialization Time (s): ", time0
        # Main integration loop
        time = 0.

        self.trajectoryActions(0)

        #print "Before integration step"
        #print g
        #print gcm

        for step in range(nsteps):
            #timetransstart=clock()
	    
            # First application of thermostat
            self.applyThermostat(vcm, nmv, mcm, bdcm, delta_t, beta)
            # First integration half-step
            for i in range(nbeads_mol):
                for j in range(3):
                    dv[i, j] = -0.5*delta_t*gcm[i, j]/mcm[i]
                    vcm[i, j] += dv[i, j]
            # Remove frozen subspace
            self.freeze(vcm, ss)

            #print "After Apply Thermostat"
            #print g
            #print gcm

            # Conversion to normal mode coordinates
            for i in range(nbeads_mol):
                # bd[i, 0] == 0 means "first bead of an atom"
                if bdcm[i, 0] == 0:
                    self.fixBeadPositions(xcm, i, bdcm[i, 1])
                    self.cartesianToNormalMode(xcm, nmc, i, bdcm[i, 1])
                    self.cartesianToNormalMode(vcm, nmv, i, bdcm[i, 1])

            # Harmonic oscillator time propagation
            for i in range(nbeads_mol):
                # bd[i, 0] == 0 means "first bead of an atom"
                if bdcm[i, 0] == 0:
                    self.propagateOscillators(nmc, nmv, i, bdcm[i, 1], beta, delta_t)
            # Conversion back to Cartesian coordinates
            for i in range(nbeads_mol):
                # bd[i, 0] == 0 means "first bead of an atom"
                if bdcm[i, 0] == 0:
                    self.normalModeToCartesian(xcm, nmc, i, bdcm[i, 1])
                    self.normalModeToCartesian(vcm, nmv, i, bdcm[i, 1])


            # Mid-step energy calculation
            self.cmtoatom(x,v,g,m,xcm,vcm,gcm,mcm,Nmol)

            self.calculateEnergies(x, &energy, 1)
            self.atomtocm(x,v,g,m,xcm,vcm,gcm,mcm,bdcm,Nmol)

            #print "After Energy Calculation"
            #print g
            #print gcm


            # Quantum energy estimators
            se = self.springEnergyNormalModes(nmc, mcm, bdcm, beta)

            qe_prim = energy.energy - se + 0.5*df/beta
            #qe_vir = energy.energy - 0.5*energy.virial
            #qe_cvir = energy.energy \
            #          - 0.5*self.centroidVirial(x, g, bd) \
            #          + 0.5*cdf/beta

            # Second integration half-step
            for i in range(nbeads_mol):
                for j in range(3):
                    dv[i, j] = -0.5*delta_t*gcm[i, j]/mcm[i]
                    vcm[i, j] += dv[i, j]
            # Second application of thermostat

            #print "After Conversion to CM"
            #print g
            #print gcm

            self.applyThermostat(vcm, nmv, mcm, bdcm, delta_t, beta)

            # Remove frozen subspace
            self.freeze(vcm, ss)

            #CHECK
            for i in range(nbeads_mol):
                if bdcm[i, 0] == 0:
                    self.cartesianToNormalMode(vcm, nmv, i, bdcm[i, 1])


            # Calculate kinetic energy
            ke = 0.
            for i in range(nbeads_mol):
                for j in range(3):
                    ke += 0.5*mcm[i]*vcm[i, j]*vcm[i, j]
            temperature = 2.*ke/(df*k_B)
            if False:
                ke_nm = 0.
                for i in range(nbeads_mol):
                    if bdcm[i, 0] == 0:
                        for j in range(3):
                            for k in range(bdcm[i,1]):
                                if k == 0 or (bdcm[i,1] % 2 == 0 and k == bdcm[i,1]/2):
                                    ke_nm += 0.5*mcm[i]*nmv[j, i+k]*nmv[j, i+k]/bdcm[i,1]
                                else:
                                    ke_nm += mcm[i]*nmv[j, i+k]*nmv[j, i+k]/bdcm[i,1]
                assert fabs(ke-ke_nm) < 1.e-7

            self.cmtoatom(x,v,g,m,xcm,vcm,gcm,mcm,Nmol)
            pot_old=energy.energy

            #print "After Translation"
            #print g
            #print gcm

            #timerotstart=clock()
            #timetrans+=(<double> (timerotstart - timetransstart)) / CLOCKS_PER_SEC
           

            #######################################
            ### PERFORM MC RIGID BODY ROTATIONS ###
            #######################################
            if (step%rotskipstep == 0):
                nrotsteps+=1
                rotbdcount=1
                rotbdskip=1
                for stp in range(rotbdcount):
                    for t1b in range(stp%rotbdskip,P,rotbdskip):
                        atomcount=0
                        for a in range(Nmol):
                            natomspmol=self.universe.objectList()[a].numberOfAtoms()
    
                            if (natomspmol==1):
                                atomcount+=natomspmol
                                continue
    
                            #timeinitstart=clock()
                            t0b = t1b - 1
                            t2b = t1b + 1
                     
                            if (t0b < 0):  # if the middle bead is the first bead, we make the left-most bead the last bead.
                                t0b += P
                            
                            if (t2b > (P-1)):  # if the middle bead is the last bead, we make the right-most bead the first bead.   
                                t2b -= P
    
                            t0 = a*P + t0b
                            t1 = a*P + t1b
                            t2 = a*P + t2b
    	         
                            atombead = atomcount*P+t1b
                            xold = N.zeros((natomspmol,3),N.float)
                            for i in range(natomspmol):
                                for j in range(3):
                                    xold[i,j] = x[atombead+i*P,j]

                            propphi = MCAngles[0][t1] + 2.0*N.pi*rotstep*(N.random.random()-0.5)
                            propchi = MCAngles[2][t1] + 2.0*N.pi*rotstep*(N.random.random()-0.5)
                            propcth = N.cos(MCAngles[1][t1]) + rotstep*(N.random.random()-0.5)
                            
                            if (propphi<0.0):
                                propphi = 2.0*N.pi + propphi
                            if (propchi<0.0):
                                propchi = 2.0*N.pi + propchi
                            propphi = propphi % (2.0*N.pi)
                            propchi = propchi % (2.0*N.pi)

                            if (propcth > 1.0):
                                propcth = 2.0 - propcth
                            elif (propcth < (-1.0)):
                                propcth = -2.0 - propcth

                            propthe = N.arccos(propcth)

                            MCAngprop[0] = propphi
                            MCAngprop[1] = propthe
                            MCAngprop[2] = propchi

#########################################
##OLD DENSITY
                            #timepot1start=clock()
                            #timeanalysis1start=clock()
                            rho = 0.0

                            istop = 0
                            Eulan1[0]=MCAngles[0][t0]
                            Eulan1[1]=MCAngles[1][t0]
                            Eulan1[2]=MCAngles[2][t0]
                            Eulan2[0]=MCAngles[0][t1]
                            Eulan2[1]=MCAngles[1][t1]
                            Eulan2[2]=MCAngles[2][t1]
                            
                            rho = self.rotden(Eulan1, Eulan2, rhomat, erotmat, esqmat, istop)[0]
                            dens_old = rho

                            Eulan1[0]=MCAngles[0][t1]
                            Eulan1[1]=MCAngles[1][t1]
                            Eulan1[2]=MCAngles[2][t1]
                            Eulan2[0]=MCAngles[0][t2]
                            Eulan2[1]=MCAngles[1][t2]
                            Eulan2[2]=MCAngles[2][t2]
                            
                            istop=0
                            
                            rho = self.rotden(Eulan1, Eulan2, rhomat, erotmat, esqmat, istop)[0]
                            dens_old = dens_old*rho

                            if (fabs(dens_old)<(1.0e-10)):
                                dens_old=0.0
                            if (dens_old < 0.0):
                                print "Rotational Density Negative"
                                raise()
                            
###########################################
##NEW DENSITY

                            Eulan1[0] = MCAngles[0][t0]
                            Eulan1[1] = MCAngles[1][t0]
                            Eulan1[2] = MCAngles[2][t0]
                            Eulan2[0] = MCAngprop[0]
                            Eulan2[1] = MCAngprop[1]
                            Eulan2[2] = MCAngprop[2]
                            
                            istop = 0
                             
                            rho = self.rotden(Eulan1, Eulan2, rhomat, erotmat, esqmat, istop)[0]
                            dens_new = rho
                           
                            Eulan1[0] = MCAngprop[0]
                            Eulan1[1] = MCAngprop[1]
                            Eulan1[2] = MCAngprop[2]
                            Eulan2[0] = MCAngles[0][t2]
                            Eulan2[1] = MCAngles[1][t2]
                            Eulan2[2] = MCAngles[2][t2]
 
                            istop = 0
                             
                            rho = self.rotden(Eulan1, Eulan2, rhomat, erotmat, esqmat, istop)[0]
                            dens_new = dens_new*rho

                            if (fabs(dens_new)<(1.0e-10)):
                                dens_new=0.0
                            if (dens_new < 0.0):
                                print "Rotational Density Negative"
                                raise()

                            #timeconvertstart=clock()
                            self.eulertocart(t1b, a, rcm_stand, x, MCAngprop, xcm)
    
                            #timepot2start=clock()
                            pot_new=self.energyCalculator(N.asarray(x))
    
                            #timeanalysis2start=clock()
                            rd = 1.0
                            if (dens_old > (1.0e-10)):
                                rd = dens_new/dens_old
    	                    
                            rd *= exp(-(beta/P)*(pot_new-pot_old))
    	         
                            accept = False
                            if (rd > 1.0):
                                accept = True
                            elif (rd > N.random.random()):
                                accept = True
    	                
                            if accept:
                                pot_old = pot_new
                                acceptratio += 1.0
                                for co in range(3):
                                    MCAngles[co][t1]=MCAngprop[co]
                            else:
                                for i in range(natomspmol):
                                    for j in range(3):
                                        x[atombead+i*P,j]=xold[i,j]
                            
                            atomcount+=natomspmol

            #Rotational Energy should be calculated for all steps
            qe_rot=0.0
            for a in range(Nmol):
                if (self.universe.objectList()[a].numberOfAtoms() == 1):
                    continue
                srot = 0.0
                for t1b in range(P):
                    t0b = t1b - 1
                    if (t0b < 0): 
                        t0b += P
                    
                    t0 = a*P + t0b
                    t1 = a*P + t1b

                    Eulan1[0]=MCAngles[0][t0]
                    Eulan1[1]=MCAngles[1][t0]
                    Eulan1[2]=MCAngles[2][t0]
                    Eulan2[0]=MCAngles[0][t1]
                    Eulan2[1]=MCAngles[1][t1]
                    Eulan2[2]=MCAngles[2][t1]

                    istop = 0
                     
                    erot = self.rotden(Eulan1, Eulan2, rhomat, erotmat, esqmat, istop)[1]
                    esq = self.rotden(Eulan1, Eulan2, rhomat, erotmat, esqmat, istop)[2]
                    srot += erot

                srot = srot/P
                qe_rot += srot

            self.calculateEnergies(x, &energy, 0)
            self.atomtocm(x,v,g,m,xcm,vcm,gcm,mcm,bdcm,Nmol)

            # End of time step
            time += delta_t
            self.foldCoordinatesIntoBox()
            self.trajectoryActions(step+1)

        # Release the write lock.
        self.releaseWriteLock()

        # Finalize all trajectory actions (close files etc.)
        self.finalizeTrajectoryActions(nsteps)

        # Deallocate the Fourier transform workspace
        self.workspace_ptr_1 = NULL
        self.workspace_ptr_2 = NULL
        self.workspace1 = None
        self.workspace2 = None

        acceptratio/=Ntruemol*float(P*nrotsteps*rotbdcount/rotbdskip)
        print "Acceptance Ratio: ", acceptratio

#
# Velocity Verlet integrator in normal-mode coordinates
# with a Langevin thermostat
#
cdef class Rigid3DRotor_PILangevinNormalModeIntegrator(Rigid3DRotor_PINormalModeIntegrator):

    """
    Molecular dynamics integrator for path integral systems using
    normal-mode coordinates and a Langevin thermostat.

    This integrator works like PINormalModeIntegrator, but has
    an additional option "centroid_friction", which is a ParticleScalar
    (one friction constant per atom) or a plain number.

    """

    cdef N.ndarray gamma
    
    cdef void applyThermostat(self, N.ndarray[double, ndim=2] v, N.ndarray[double, ndim=2] nmv,
                              N.ndarray[double, ndim=1] m, N.ndarray[N.int32_t, ndim=2] bd,
                              double dt, double beta):
        cdef N.ndarray[double, ndim=1] g = self.gamma
        cdef int nbeads = v.shape[0]
        cdef double f, c1, c2
        cdef double omega_n, mb
        cdef Py_ssize_t i, j, k
        cdef int32_t nb
        for i in range(nbeads):
            # bd[i, 0] == 0 means "first bead of an atom"
            if bd[i, 0] == 0:
                nb = bd[i, 1]
                # Conversion to normal mode coordinates
                self.cartesianToNormalMode(v, nmv, i, nb)
                # Modify velocities
                omega_n = nb/(beta*hbar)
                mb = sqrt(nb/(beta*m[i]))
                for k in range(nb):
                    if k == 0:
                        f = g[i]
                    else:
                        f = 4.*omega_n*sin(k*M_PI/nb)
                    c1 = exp(-0.5*dt*f)
                    c2 = sqrt(1-c1*c1)
                    for j in range(3):
                        if k == 0 or (nb % 2 == 0 and k == nb/2):
                            nmv[j, i+k] = c1*nmv[j, i+k] + c2*mb*MMTK.mtrand.standard_normal()
                        else:
                            nmv[j, i+k] = c1*nmv[j, i+k] + sqrt(0.5)*c2*mb*MMTK.mtrand.standard_normal()
                # Conversion back to Cartesian coordinates
                self.normalModeToCartesian(v, nmv, i, nb)

    cdef start(self):
        friction = self.getOption('centroid_friction')
        self.denrho=self.getOption('denrho')
        self.denerot=self.getOption('denerot')
        self.denesq=self.getOption('denesq')
        self.rotmove=self.getOption('rotstep')
        self.rotstepskip=self.getOption('rotskipstep')
        if isinstance(friction, ParticleProperties.ParticleScalar):
            self.gamma = friction.array
        else:
            assert isinstance(friction, numbers.Number)
            nbeads = self.universe.numberOfPoints()
            self.gamma = N.zeros((nbeads,), N.float)+friction
        Rigid3DRotor_PINormalModeIntegrator.start(self)


#!/usr/bin/python

from __future__ import division
import numpy as np
from numpy import sin, cos, exp, sqrt, pi, abs, newaxis
import scipy
from scipy.special import j0, j1, jn
from scipy.integrate import simps

I = complex(0,1)

def energy_density(r, phi, z, alpha, beta, k, nthetas=1000):
    """
    :param r: radial coordinate
    :param z: axial coordinate
    :param phi: azimuthal angle (0 is the X axis)
    :param k: wavenumber
    :param alpha: half-angle of objective opening in radians
    :param beta: the ratio of the back-aperture radius to the beam radius
    """
    # integration variables
    dtheta = alpha / nthetas
    theta = np.linspace(0, alpha, nthetas)[newaxis,:]

    # coordinates
    v = k * np.array(r) * sin(alpha)
    u = k * np.array(z) * sin(alpha)**2

    # convenient quantities
    q = sin(theta) / sin(alpha)
    qq = exp(I * u * cos(theta) / sin(alpha)**2)
    A = sqrt(cos(theta) * exp(-2 * beta**2 * q**2))

    vq = v[:,newaxis] * q
    psi0 = simps(A * sin(theta) * (1 + cos(theta)) * j0(vq) * qq, dx=dtheta)
    psi1 = simps(A * sin(theta)**2 * j1(vq) * qq, dx=dtheta)
    psi2 = simps(A * sin(theta) * (1 - cos(theta)) * jn(2, vq) * qq, dx=dtheta)
    energy = 1 / 16 / pi * (abs(psi0)**2 + 4*abs(psi1)**2 * cos(phi)**2 + abs(psi2)**2  + 2*cos(2*phi) * np.real(psi0 * psi2.conj()))
    return energy

# length in nanometers
alpha = 1.12
beta = 0.9
k = 2*pi / 514
r = np.linspace(0, 500, 1000)
res = energy_density(r, 0.1, 0.1, alpha, beta, k)

from matplotlib import pyplot as pl
pl.plot(res)
pl.show()

def corr(initial, range, alpha, beta, k, nthetas=1000):
    """
    :type initial: array of shape ''(nr, nphi, ntheta)''
    :param initial: the initial concentration field
    :type param: tuple of shape ''(r, phi, z)''
    """
    npts = initial.shape
    coords = [np.linspace(-range[i], range[i], npts[i]) for i in range(3)]
    r,phi,z = np.mgrid[coords]
    density = energy_density(r, phi, z, alpha, beta, k, nthetas=nthetas)
    A = scipy.fftconvolve(conc, density)



import matplotlib.pyplot as plt
import numpy as np

import math
from math import cos, sin, tan, pi, atan2, sqrt

N_STATE = 4

# state = [x, y, theta, v, 1/L, k/m].Transpose
# u = [pwm, steer].Transpose


P = np.diag(np.array([1,1,0.02,0.01]))

Qt = P

def compute_jacobian(dt, state, parameter, u):

	pwm = u[0,0]
	steer = u[1,0]

	v = state[3,0]
	theta = state[2,0]


	L_inv = parameter[0,0]
	K_m = parameter[1,0]


	G = np.matrix([	[ 1, 0, -dt*v*sin(theta), dt*cos(theta)],
					[ 0, 1, dt*v*cos(theta), dt*sin(theta)],
					[ 0, 0, 1, dt*L_inv*tan(steer)],
					[ 0, 0, 0, 1]])

	H = np.matrix([ [1, 0, 0, 0],
					[0, 1, 0, 0],
					# [0, 0, 1, 0],
					[0, 0, 1, 0],
					[0, 0, 0, 1]])

	Jp = np.matrix([[v*tan(steer)*dt, 0],
					[0, pwm*dt]])

	return G, H

N = 7;
PARAMETER_PARTICLES = np.matrix(np.random.rand(2, N))
PARAMETER_PARTICLES = np.multiply(PARAMETER_PARTICLES, [[4], [100]])

pw = np.zeros(N) + 1.0/N;

E = [[0.01], [0.2]];
dt = 0.01
sigma1 = 100.0*dt*np.pi/180
sigma2 = 10.0*dt

# print(PARAMETER_PARTICLES)
topper_particle = None
topper_pw = 0

def gauss_likelihood(x, sigma):
	p = 1.0 * np.exp(-x ** 2 / (2 * sigma ** 2)) / np.sqrt(2.0 * np.pi * sigma ** 2)
	return p

c = 30
def resampling():
	"""
	best candidate re-sampling
	"""

	global PARAMETER_PARTICLES, pw, topper_particle, topper_pw

	i = np.argmax(pw)
	winner = PARAMETER_PARTICLES[:, i]

	# if(pw[i] > topper_pw):
	# 	topper_particle = winner
	# 	topper_pw = pw[i]

	PARAMETER_PARTICLES = np.repeat(winner, N, axis=1)

	# NP = N
	# NTh = 0.5 * N
	#
	# Neff = 1.0 / (np.matmul(pw, pw))  # Effective particle number
	# if Neff < NTh:
	# 	wcum = np.cumsum(pw)
	#
	# 	r = np.random.rand()/NP
	#
	# 	inds = []
	#
	# 	for ip in range(NP):
	# 		ind = 0
	# 		while r*ip/N > wcum[ind]:
	# 			ind += 1
	# 		inds.append(ind)
	#
	# 	PARAMETER_PARTICLES = PARAMETER_PARTICLES[:, inds]
	# 	pw = np.zeros(NP) + 1.0 / NP  # init weight

	
cov_evolution = np.matrix(np.ones((2,1)))
def run(cycle, u, z):
	global Qt, PARAMETER_PARTICLES, pw, c, cov_evolution


	# Particle filter for parameter
	evolution = np.random.randn(2, N);
	evolution = np.multiply(evolution, E)
	evolution = np.multiply(evolution, cov_evolution)
	PARAMETER_PARTICLES = PARAMETER_PARTICLES + evolution
	for i in range(N):
		# Predict Z
		para = PARAMETER_PARTICLES[:,i]
		# para[0,0] = 1
		pw[i] = 0.
		cov = np.matrix(np.zeros((2,1)))
		for k in range(c):
			xp, zp = cycle.step(cycle.dt, cycle.state, para, u);
			# Assign weight
			dz = zp - z
			cov = cov + dz[2:]

			pw[i] = pw[i] + gauss_likelihood(dz[0,0], sigma1) * gauss_likelihood(dz[1,0], sigma1) * gauss_likelihood(dz[2, 0], sigma1) * gauss_likelihood(dz[3, 0], sigma1)

		cov_evolution = np.fabs(cov/c/dt);
		# print(cov_evolution)


	psum = pw.sum()
	pw = pw/psum

	# print(PARAMETER_PARTICLES)
	# print(pw)

	pEst = np.matmul(PARAMETER_PARTICLES, pw).T

	resampling()

	# print(pEst)
	# print(PARAMETER_PARTICLES)
	# print("pest=> ", pEst)
	# print("p=> ", PARAMETER_PARTICLES)

	cycle.set_parameter(pEst)

	# pEst = cycle.parameter



	# EKF with new parameter estimate
	xp, zp = cycle.step(cycle.dt, cycle.state, cycle.parameter, u);

	dz = z - zp

	G, H = compute_jacobian(cycle.dt, cycle.state, cycle.parameter, u)

	Qtp = np.matmul(G, np.matmul(Qt, G.T)) + R;


	S = np.matmul(H, np.matmul(Qtp, H.T)) + Q
	S_inv = np.linalg.pinv(S)
	K = np.matmul(Qtp, np.matmul(H.T, S_inv))

	xEst = xp + np.matmul(K, dz)
	Qt = np.matmul( (np.identity(N_STATE) - np.matmul(K, H)) , Qt);

	# xEst = xp

	return xEst, Qt, pEst








import matplotlib.pyplot as plt
import numpy as np

import math
from math import cos, sin, tan, pi, atan2, sqrt

from world import World


from bicycle import Bicycle

R_SIM = np.diag([0.01, 0.01, 0.01*pi/180, 0.001, 0, 0]);
Q_SIM = np.diag([0.05, 0.05, 0.1, 0.3]);


R_MODEL = np.diag([0.01, 0.01, 0.01*pi/180, 0.001, 0, 0]);
Q_MODEL = np.diag([0.05, 0.05, 0.1, 0.3])





W = World(100, 100)

N = 1000;
control = np.random.rand(2, N) - 0.5;


state = np.matrix([[10, 20, 0, 0, 1, 1]]).T


cycle_sim = Bicycle(state = state, color="green", name="Original");
cycle_sim.process_noise(R_SIM)
cycle_sim.observation_noise(Q_SIM)
W.add_agent(cycle_sim)


state = np.matrix([[10, 20, 0, 0, 1, 0.1]]).T

cycle_model = Bicycle(state = state, color="red", name="Model");
cycle_model.process_noise(R_MODEL)
cycle_model.observation_noise(Q_MODEL)
W.add_agent(cycle_model)



for i in range(N):
	# u = control[:, i:i+1]
	u = np.matrix([[0.2],[ 0.1]])

	# Move actual system
	nX, nY = cycle_sim.move(u)
	nX, nY = cycle_model.step(cycle_model.dt, cycle_model.state, u)
	cycle_model.set(nX)

	if(i%W.frequency == 0):
		W.plot()
		# for particle in P.px:
		# 	plt.scatter((particle[0]), (particle[1]), s= 10, color='blue')





# velocity
def summary(name, i, type = "X"):
	fig1, ax = plt.subplots()


	ax.set_title(name)
	for a in W.agents:
		if a.name != "Original" and type == 'Y':
			continue
		data = a.historyX.A;
		if(type == 'Y'):
			data = a.historyY.A
		vel = data[i, 0:]
		ax.plot(vel, color=a.color, label=a.name)

	ax.legend(loc='upper right')

summary("velocity", 3)
summary("Theta", 2)
summary("Measured w", 2, type="Y")
summary("Measured acce", 3, type="Y")

plt.show()

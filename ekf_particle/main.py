import matplotlib.pyplot as plt
import numpy as np

import math
from math import cos, sin, tan, pi, atan2, sqrt, fabs

from world import World

from bicycle import Bicycle

import EKF

a = 0.1
R_SIM = np.diag([0.001, 0.001, 0.01 * pi / 180, 0.01])*a;
Q_SIM = np.diag([0.01, 0.01, 0.1 * pi / 180, 0.01])*a;

R_MODEL = np.diag([0.001, 0.001, 0.01 * pi / 180, 0.01])*a;
Q_MODEL = np.diag([0.02, 0.02, 0.2 * pi / 180, 0.015])*a
# x,y,w,a

EKF.R = R_MODEL
EKF.Q = Q_MODEL

# state = [x, y, theta, v, 1/L, k/m].Transpose
# u = [pwm, steer].Transpose
# [L_inv, K_by_M].Transpose


W = World(100, 100)
# W.load('track2.pkl')

N = 1000;
control = np.random.rand(2, N) - 0.5;

state = np.matrix([[10, 20, 0, 0]]).T
parameter = np.matrix([[1, 50]]).T

cycle_sim = Bicycle(state=state, parameter=parameter, color="green", name="Original");
cycle_sim.process_noise(R_SIM)
cycle_sim.observation_noise(Q_SIM)
W.add_agent(cycle_sim)

state = np.matrix([[10, 20, 0, 0]]).T
parameter1 = np.matrix([[0.1, 5]]).T

cycle_model = Bicycle(state=state, parameter=parameter1, color="red", name="Model");
cycle_model.process_noise(R_MODEL)
cycle_model.observation_noise(Q_MODEL)
W.add_agent(cycle_model)

target = 1
# path = W.track["path"][0:, 0:]
# L = len(W.track["path"])
vmax = 0.5

p = np.matrix(np.zeros((2, 1)));

for i in range(N):
    # u = control[:, i:i+1]
    # steer = atan2(path[target][1] - cycle_sim.state[1,0], path[target][0] - cycle_sim.state[0,0]) - cycle_sim.state[2,0]

    # steer = (steer + pi)%(2*pi) - pi;
    # if(fabs(steer) > pi/2) and sqrt((path[target][1] - cycle_sim.state[1,0])**2 + (path[target][0] - cycle_sim.state[0,0])**2) < 5:
    # 	steer = 0.;
    # 	target = (target + 1)%L;
    # steer = min(pi/2, max(-pi/2, steer))
    # u = np.matrix([[vmax - cycle_sim.state[3,0]],[ steer]])
    u = np.matrix([[(4. + 0.1 * np.random.randn() - 1.0*cycle_sim.state[3, 0] + 2.0 * sin(0.01 * i) ) * 0.05],
                   [0.0 * sin(0.004 * i) + 0.1 + 0.2 * np.random.randn()]])

    # Move actual system
    nX, nZ = cycle_sim.move(u)

    # xEst, nY = cycle_model.step(cycle_model.dt, cycle_model.state, parameter, u);
    # cycle_model.historyY = np.hstack((cycle_model.historyY, nY));

    xEst, Qt, pEst = EKF.run(cycle_model, u, nZ);

    # xEst = nX
    print(pEst)
    p = np.hstack((p, pEst))

    cycle_model.set(xEst)

    if (i % W.frequency == 0):
        W.plot()


# velocity
def summary(name, i, type="X"):
    fig1, ax = plt.subplots()

    ax.set_title(name)
    for a in W.agents:
        if a.name != "Original" and type == 'Y':
            continue
        data = a.historyX.A;
        if (type == 'Y'):
            data = a.historyY.A
        val = data[i, 0:]
        ax.plot(val, color=a.color, label=a.name)

    ax.legend(loc='upper right')

def parameter(name, d):
    fig1, ax = plt.subplots()

    ax.set_title(name)
    ax.plot(d, color="red")


# summary("X", 0)
# summary("Y", 1)
summary("V", 3)
# summary("theta", 2)
p = p.A
parameter("L_inv 1", p[0,:])
parameter("k_by_m 50", p[1,:])

# summary("Measured w", 2, type="Y")
# summary("Measured acce", 3, type="Y")

plt.show()

import scipy.integrate as integrate
import matplotlib.pyplot as plt
import numpy as np

#dT/dt = -l/T
l = 1.0/5.0 #min
t0 = 130 #F

def dT_dt(self, T):
	return -(1/5)*T

y0 = 0.0
yf = 10.0
sol = integrate.solve_ivp(dT_dt, (y0, yf), [t0,])
print([t0,])

print(sol.t)
print(sol.y)
#print(sol.t.shape, sol.y.reshape([33,]).shape)
#print(type(sol.t))

plt.plot(sol.t, np.log(sol.y.reshape([5,])))
#plt.plot(sol.y.reshape([33,]))
plt.show()
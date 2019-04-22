import scipy.integrate as integrate
import scipy.interpolate as interpolate
import matplotlib.pyplot as plt
import numpy as np

pi = np.pi
c1 = (1.0 / 0.118)
c2 = (1.0 / 0.0337)
c3 = (1.0 / 0.0170)
b1 = (41.0 / 10.0)
b2 = (-19.0 / 6.0)
b3 = -7.0
yt_initial = 0.955
yb_initial = 0.0174
ytao_initial = 0.0102
M_z = 91.2 #GeV
M_susy = 1000.0 #GeV, 1 TeV
#m_z = 91.2 GeV; M_susy = 1 TeV = 1000 GeV
t0 = np.log(M_z / M_z)
tf = np.log(M_susy / M_z)
scale = (t0, tf)
initial_conditions = [yt_initial, yb_initial, ytao_initial]



def dy_dt (t, Y): #return a vector
	#g**2; function of t
	g1 = np.sqrt((4*pi) / ((-b1/(2*pi))*t + c1))
	g2 = np.sqrt((4*pi) / ((-b2/(2*pi))*t + c2))
	g3 = np.sqrt((4*pi) / ((-b3/(2*pi))*t + c3))

	dy_dt = ((1.0 / (16.0*(pi**2)))*Y[0] * (((3.0/2.0)*(Y[0]**2 - Y[1]**2)) +
				(3.0*(Y[0]**2 + Y[1]**2)) + (Y[2]**2) - (8.0*(g3**2)) - ((9.0/4.0)*(g2**2)) - ((17.0/20.0)*(g1**2))))
	dyb_dt = ((1.0 / (16.0*(pi**2)))*Y[1] * (((3.0/2.0)*(Y[1]**2 - Y[0]**2)) +
				(3.0*(Y[0]**2 + Y[1]**2)) + (Y[2]**2) - (8.0*(g3**2)) - ((9.0/4.0)*(g2**2)) - ((1.0/4.0)*(g1**2))))
	dytao_dt = ((1.0 / (16.0*(pi**2)))*Y[2] * (((3.0/2.0)*(Y[2]**2)) +
				(3.0*(Y[0]**2 + Y[1]**2)) + (Y[2]**2) - ((9.0/4.0)*(g2**2)) - ((9.0/4.0)*(g2**2)) - ((9.0/4.0)*(g1**2)))) 

	return [dy_dt, dyb_dt, dytao_dt]
	

sol_y = integrate.solve_ivp(dy_dt, scale, initial_conditions)
#print(sol_y)
#plt.plot(sol_y.t, sol_y.y.T)
#plt.show()

dydt_sol = sol_y[0]
dybdt_sol = sol_y[1]
dytaodt_sol = sol_y[2]

def dlambda_dt (t, Y):

	def dy_dt (t):
		return 

	def dyb_dt (t):
		return

	def dytao_dt (t):
		return
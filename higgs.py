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
M_susy_stop = 1.0e16
#m_z = 91.2 GeV; M_susy = 1 TeV = 1000 GeV
#inputs for lambda
Mstop1 = M_susy #1000 GeV, 1 TeV
Mstop2 = M_susy + 100.0 #1000 GeV, 1 TeV
mu = 1000.0 #1000 GeV, 1 TeV
beta = (pi*40.0 / 180.0) #change degrees
beta = np.arctan(40.0)
At = 0 #1000 GeV, 1 TeV
t0 = np.log(M_z / M_z)
tf = np.log(M_susy / M_z)
scale_y = (t0, tf)
scale_lambda = (tf, t0)
initial_conditions_y = [yt_initial, yb_initial, ytao_initial]
v = 246 #246 GeV
m_h = 125 #GeV

#g**2; function of t
def g1 (t): 
	return np.sqrt((4*pi) / ((-b1/(2*pi))*t + c1))
def g2 (t):
	return np.sqrt((4*pi) / ((-b2/(2*pi))*t + c2))
def g3 (t):
	return np.sqrt((4*pi) / ((-b3/(2*pi))*t + c3))


def dy_dt (t, Y): #return a vector

	dy_dt = ((1.0 / (16.0*(pi**2)))*Y[0] * (((3.0/2.0)*(Y[0]**2 - Y[1]**2)) +
				(3.0*(Y[0]**2 + Y[1]**2)) + (Y[2]**2) - (8.0*(g3(t)**2)) - ((9.0/4.0)*(g2(t)**2)) - ((17.0/20.0)*(g1(t)**2))))
	dyb_dt = ((1.0 / (16.0*(pi**2)))*Y[1] * (((3.0/2.0)*(Y[1]**2 - Y[0]**2)) +
				(3.0*(Y[0]**2 + Y[1]**2)) + (Y[2]**2) - (8.0*(g3(t)**2)) - ((9.0/4.0)*(g2(t)**2)) - ((1.0/4.0)*(g1(t)**2))))
	dytao_dt = ((1.0 / (16.0*(pi**2)))*Y[2] * (((3.0/2.0)*(Y[2]**2)) +
				(3.0*(Y[0]**2 + Y[1]**2)) + (Y[2]**2) - ((9.0/4.0)*(g2(t)**2)) - ((9.0/4.0)*(g2(t)**2)) - ((9.0/4.0)*(g1(t)**2)))) 

	return [dy_dt, dyb_dt, dytao_dt]
	

sol_y = integrate.solve_ivp(dy_dt, scale_y, initial_conditions_y)
print(sol_y)
plt.plot(sol_y.t, sol_y.y.T)
plt.show()

#solutions of dy_dt
dytdt_sol = sol_y.y[0]
dybdt_sol = sol_y.y[1]
dytaodt_sol = sol_y.y[2]

#yt(Msusy) = Yt(Msusy)sin(beta)
Yt = dytdt_sol[2] / np.sin(beta)

def Yt_func(Y):
	return Y / np.sin(beta)

#x = Mstop1 / Mstop2
def F (x):
	return (2.0*np.log(x))/(x**2 - 1.0)
def G (x):
	return ((12.0*x**2)*(1.0 - (x**2) + (1.0+(x**2))*np.log(x))) / (((x**2) - 1.0)**3)

def Xt (At, mu, beta):
	return At - (mu*(1.0 / np.tan(beta)))

#correct
def delta_lambda (Mstop1, Mstop2, Yt):
	return (3.0*(Yt**4)*(2.0*((Xt(At, mu, beta)**2)/(Mstop1*Mstop2))*F(Mstop1/Mstop2) - ((1.0/6.0)*(Xt(At, mu, beta)**4)/(Mstop1**2*Mstop2**2)*G(Mstop1/Mstop2)))) / (16*(pi**2))


initial_conditions_lambda = [(1.0 / 4.0)*((g2(tf)**2) + ((3.0 / 5.0)*g1(tf)**2))*(np.cos(2.0*beta)**2) + delta_lambda(Mstop1, Mstop2, Yt),]

#x = sol_y.t; y = sol_y.y indexed; kind = linear
#interpolation returns a function; pass numbers to the function to be estimated/evaluated
yt = interpolate.interp1d(sol_y.t, dytdt_sol, kind="linear")
yb = interpolate.interp1d(sol_y.t, dybdt_sol, kind="linear")
ytao = interpolate.interp1d(sol_y.t, dytaodt_sol, kind="linear")

#correct
#Y = lambda
def dlambda_dt (t, Y): 
	return (4.0*Y*((3.0*(yt(t)**2))+(3.0*(yb(t)**2))+(ytao(t))**2)) - (9.0*Y*((1.0/5.0)*(g1(t)**2))+(g2(t)**2)) - (4.0*(3.0*(yt(t)**4)+(3.0*(yb(t)**4))+(ytao(t)**4))) + ((27.0/100.0)*(g1(t)**4)) + ((9.0/10.0)*(g2(t)**2)*(g1(t)**1)) + ((9.0/4.0)*(g2(t)**4)) + (12.0*(Y**2))
	
sol_lambda = integrate.solve_ivp(dlambda_dt, scale_lambda, initial_conditions_lambda)

print(sol_lambda)
#print(sol_lambda.t, sol_lambda.y.T)
plt.plot(sol_lambda.t, sol_lambda.y.T)
plt.show()

#closest to m_h / Mz
l = interpolate.interp1d(sol_lambda.t, sol_lambda.y, kind="linear")
#print(l(m_h / M_z))

print(np.sqrt(l(m_h / M_z)) * v)

"""
#evaluate everything I've done before on a scale of M_susy -> M_susy_stop
def lambda_on_scale(M_susy):

	Mstop1 = M_susy #GeV
	Mstop2 = M_susy + 100.0 #GeV
	t0 = np.log(M_z / M_z)
	tf = np.log(M_susy / M_z)
	scale_y = (t0, tf)
	scale_lambda = (tf, t0)

	sol_y = integrate.solve_ivp(dy_dt, scale_y, initial_conditions_y)

	dytdt_sol = sol_y.y[0]
	dybdt_sol = sol_y.y[1]
	dytaodt_sol = sol_y.y[2]

	Yt = dytdt_sol[-1] / np.sin(beta)

	yt = interpolate.interp1d(sol_y.t, dytdt_sol, kind="linear")
	yb = interpolate.interp1d(sol_y.t, dybdt_sol, kind="linear")
	ytao = interpolate.interp1d(sol_y.t, dytaodt_sol, kind="linear")

	def dlambda_dt (t, Y): 
		return (4.0*Y*((3.0*(yt(t)**2))+(3.0*(yb(t)**2))+(ytao(t))**2)) - (9.0*Y*((1.0/5.0)*(g1(t)**2))+(g2(t)**2)) - (4.0*(3.0*(yt(t)**4)+(3.0*(yb(t)**4))+(ytao(t)**4))) + ((27.0/100.0)*(g1(t)**4)) + ((9.0/10.0)*(g2(t)**2)*(g1(t)**1)) + ((9.0/4.0)*(g2(t)**4)) + (12.0*(Y**2))

	#initial_conditions_lambda = [(1.0 / 4.0)*((g2(tf)**2) + ((3.0 / 5.0)*g1(tf)**2))*(np.cos(2.0*beta)**2) + delta_lambda(Mstop1, Mstop2, Yt),]
	g1_const = 1.220397455531664369257
	initial_conditions_lambda = [(1.0 / 4.0)*((g2(tf)**2) + ((3.0 / 5.0)*g1_const**2))*(np.cos(2.0*beta)**2) + delta_lambda(Mstop1, Mstop2, Yt),]


	sol_lambda = integrate.solve_ivp(dlambda_dt, scale_lambda, initial_conditions_lambda)

	return sol_lambda


M_susy_scale = np.linspace(1000.0, 1.0e16)
lambda_solutions = []
higgs_mass = []

for i in M_susy_scale:
	print(i)
	lambda_sol = lambda_on_scale(i)
	lambda_solutions.append(lambda_sol)
	l = interpolate.interp1d(lambda_sol.t, lambda_sol.y, kind="linear")
	higgs_mass.append(np.sqrt(l(m_h / M_z)) * v)


#print(higgs_mass)
"""
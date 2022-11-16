import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simps
from sklearn.linear_model import LinearRegression
import pandas as pd

def function(x,e):
    derivative = np.zeros(2)
    derivative[0] = x[1]
    derivative[1] = -e*(x[0])
    return derivative

def rk4(IC,x_array,N,e):
	#Step Size
	h = (x_array[-1] - x_array[0])/N
	x_o = x_array[0] ; u_o = IC[0] ; v_o = IC[1]
	u_array = [u_o] ; v_array = [v_o]
	
	for i in range(N):
		psi = [u_o,v_o,x_o]
		k1 = h * function(psi,e)
		psi = [u_o + k1[0]/2,v_o + k1[1]/2,x_o + h/2]
		k2 = h * function(psi,e)
		psi = [u_o + k2[0]/2,v_o + k2[1]/2,x_o + h/2]
		k3 = h * function(psi,e)
		psi = [u_o + k3[0],v_o + k3[1],x_o + h]
		k4 = h * function(psi,e)
		
		u_o = u_o + 1/6 * (k1[0] + 2*(k2[0] + k3[0]) + k4[0])
		v_o = v_o + 1/6 * (k1[1] + 2*(k2[1] + k3[1]) + k4[1])
		
		u_array.append(u_o) ; v_array.append(v_o)
		x_o = x_o + h
		
	return u_array,v_array

def normalisation(array,x_array):
	A = 1/np.sqrt((simps((array)**2,x_array)))
	u_norm = A*array
	return u_norm

def secant_method(e,u): # phi = [phi(s_i),phi(s_i+1)]
	e2 = e[1] - ((e[1]-e[0])/(u[1]-u[0]))*u[1]
	return e2

def analytic_sol(L,x,n):
    psiodd = 0 ; psieven = 0
    if n % 2 != 0:
        k = n*(np.pi/L)
        psiodd = np.sqrt(2/L) * np.cos(k*x)
        return psiodd
    else:
        k = n * (np.pi / L)
        psieven = np.sqrt(2/L) * np.sin(k*x)
        return psieven
    
def eV(m,L,e_approx):
    h=6.63* 10**(-34)
    Eval=[]
    for i in range(1,6):
        Eval.append((i**2 * np.pi**2 * h**2)/(8* m * (L**2)))
    Eval=np.array(Eval)* 6.242 * 10**(18)
    Eigenval=e_approx*((h**2)/(8* m * (L**2))) * 6.242 * 10**(18)
    dtf1=pd.DataFrame({"Analytical Energy Value": Eval, "Eigen ENergy Val": Eigenval})
    print(dtf1)

if __name__ == "__main__":
	N = 100
	x_array = np.linspace(-1/2,1/2,N+1)
	IC = [0,1]
	e_array = np.linspace(0,250,251)
	un_array = []
	u_last = []
	u_sign = []
	u_sign = []
	energy_sign = []
	indicies = []
	indicies = []

	for i in (e_array):
		sol = rk4(IC,x_array,N,i)
		un_array.append(sol[0])
  
	for i in range(len(un_array)):
		u_last.append(un_array[i][-1])

	#u_norm = normalisation(np.array(u),x_array)

	# PLOT OF u and e
	plt.plot(e_array,u_last)
	plt.xlabel("e")
	plt.ylabel("u(e)")
	plt.title("$u_r$(e) Vs e")
	plt.grid(ls = "--")
	plt.show()
 
	#print(u_last)
 
	# SIGN CHANGE AND ENERGY
	for i in range(len(u_last)-2):
		if u_last[i]*u_last[i+1] < 0:
			u_sign.append(u_last[i])
			u_sign.append(u_last[i+1])
			indicies.append(i)
			indicies.append(i+1)
   
	energy_pair1 = indicies[0:2]
	u_pair1 = u_sign[0:2]
	energy_pair2 = indicies[2:4]
	u_pair2 = u_sign[2:4]
	energy_pair3 = indicies[4:6]
	u_pair3 = u_sign[4:6]
	energy_pair4 = indicies[6:8]
	u_pair4 = u_sign[6:8]
	energy_pair5 = indicies[8:10]
	u_pair5 = u_sign[8:10]
 
 
	energy1 = secant_method(energy_pair1,u_pair1)
	energy2 = secant_method(energy_pair2,u_pair2)
	energy3 = secant_method(energy_pair3,u_pair3)
	energy4 = secant_method(energy_pair4,u_pair4)
	energy5 = secant_method(energy_pair5,u_pair5)

	energy = np.array([energy1,energy2,energy3,energy4,energy5])
	u_final = []
	u_norm_final = []
	
	for i in energy:
		u,v = rk4(IC,x_array,N,i)
		u_final.append(np.array(u))
  
	u_final = np.array(u_final)
 
	for i in range(len(u_final)):
		u_norm = normalisation(u_final[i],x_array)
		u_norm_final.append(u_norm)

	u_norm_final = np.array(u_norm_final)
 
	fig1, ax1 = plt.subplots()
	fig2, ax2 = plt.subplots()
	fig3, ax3 = plt.subplots()
	fig4, ax4 = plt.subplots()
	fig5, ax5 = plt.subplots()
 
	ax1.plot(x_array,u_final[0],label = "Unormalised")
	ax1.plot(x_array,u_norm_final[0],label = "Nomalised")
	ax2.plot(x_array,u_final[1],label = "Unormalised")
	ax2.plot(x_array,u_norm_final[1],label = "Nomalised")
	ax3.plot(x_array,u_final[2],label = "Unormalised")
	ax3.plot(x_array,u_norm_final[2],label = "Nomalised")
	ax4.plot(x_array,u_final[3],label = "Unormalised")
	ax4.plot(x_array,u_norm_final[3],label = "Nomalised")
	ax5.plot(x_array,u_final[4],label = "Unormalised")
	ax5.plot(x_array,u_norm_final[4],label = "Nomalised")
 
	ax1.set(title = "e ="+str(energy[0]))
	ax2.set(title = "e ="+str(energy[1]))
	ax3.set(title = "e ="+str(energy[2]))
	ax4.set(title = "e ="+str(energy[3]))
	ax5.set(title = "e ="+str(energy[4]))

	ax1.legend()
	ax2.legend()
	ax3.legend()
	ax4.legend()
	ax5.legend()
	ax1.grid(ls = "--")
	ax2.grid(ls = "--")
	ax3.grid(ls = "--")
	ax4.grid(ls = "--")
	ax5.grid(ls = "--")
	plt.show()
 
	# B PART CURVE FIT B/W e & n^2
 
	nsqr = np.array(energy / np.pi ** 2)
	print("nsqr: ",nsqr)
	plt.scatter(nsqr,energy)
	model = LinearRegression()
	model.fit(nsqr.reshape((-1,1)),energy)
	ypred = model.predict(nsqr.reshape((-1,1)))
	print("slope: ",model.coef_)
	print("Actual Value Of Slope: ", np.pi**2)
	print("R_sqr: ",model.score(nsqr.reshape((-1,1)),energy))

	plt.plot(nsqr,ypred,linestyle='dashdot',color='red')
	plt.grid(ls = "--")
	plt.xlabel("$n^2$")
	plt.ylabel("$e_n$")
	plt.title("Plot of $e_n$ as a function of $n^2$")
	plt.show()

	# PART C 
 
	for i in range(1,6):
		energy_c = (i**2)*(np.pi)**2
		u,v = rk4(IC,x_array,N,energy_c)
		u_norm = normalisation(np.array(u),x_array)
		anal = analytic_sol(1,x_array,i)
		plt.scatter(x_array,(u_norm)**2,label='calculated soltion',c = "r",s = 20)
		plt.plot(x_array,(anal)**2, label='analytic solution')
		plt.legend()
		plt.title("n ="+str(i))
		plt.grid(ls = "--")
		plt.show()

	# PART D & E
 
	print("--------FOR ELECTRON WHEN WIDTH OF WELL= 5 Angstrom------")
	eV(9.11*(10**(-31)),5 * (10**(-10)),energy)
	print()

	print("--------FOR ELECTRON WHEN WIDTH OF WELL= 10 Angstrom------")
	eV(9.11*(10**(-31)),10 * (10**(-10)),energy)
	print()

	print("--------FOR ELECTRON WHEN WIDTH OF WELL= 5 Fermi------")
	eV(1.67*(10**(-27)),5 * (10**(-15)),energy)
	print()
	
	
	

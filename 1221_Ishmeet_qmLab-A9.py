import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import simps
from scipy.linalg import eigh
import scipy.special as spe

def V(x,l):
    x = x[1:-1]
    v = np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            if i == j:
                v[i,j] = ((-2)/(x[i]))+(l*(l+1))/((x[i])**2)

    return v,x

def effective_potential(x,l):
    v_eff = ((-2)/(x))+(l*(l+1))/((x)**2)
    v = -2/x
    return v_eff,v

def finite_difference(a,b,n,l):
    x_array = np.linspace(a,b,n+2)
    h = x_array[1] - x_array[0]
    
    K = np.zeros((n,n)) 

    v,x_new = V(x_array,l)

    for i in range(n):
        for j in range(n):

            if i == 0:
                if i == j:
                    K[i,j] = -2
                elif j == i+1 :
                    K[i,j] = 1
            elif i == n-2:
                if i == j:
                    K[i,j] = -2
                elif j == i-1:
                    K[i,j] = 1
            else:
                if i == j:
                    K[i,j] = -2
                elif j == i+1 or j == i-1:
                    K[i,j] = 1
    
    K_new = -K/(h**2)

    A = K_new + v

    return A,x_new

def normalisation(array,x_array):
	A = 1/np.sqrt((simps((array)**2,x_array)))
	u_norm = A*array
	return u_norm

def analytic(r,n,l):
    coeff = np.sqrt((2.0/n)**3 * spe.factorial(n-l-1) /(2.0*n*spe.factorial(n+l)))
    laguerre = spe.assoc_laguerre(2.0*r/n,n-l-1,2*l+1)
    
    return coeff * np.exp(-r/n) * (2.0*r/n)**l * laguerre


if __name__ == "__main__":
    a = 10**(-14)
    b = 250
    n = 1000
    V_eff_list = [] ; V_nor_list = []
    
    fig1,ax1 = plt.subplots()

    A,x_new = finite_difference(a,b,n,0)
       
    e,v = eigh(A)
    print("\n---------------------------------------------------------------")
    print("ENERGY EIGENVALUES FOR l = 0 :-")
    print("---------------------------------------------------------------\n")
    for i in range(10):
        print("\nEnergy eigenvalue for n = "+str(i+1)+" and l = 0 is :",e[i],"\n")
    print("------------------------------------------------------------------\n")
    
    # FIRST TEN ENERGY EIGENVALUES FOR l = 1,2 :-
    
    A1,x_new = finite_difference(a,b,n,1)
    A2,x_new = finite_difference(a,b,n,2)
       
    e1,v1 = eigh(A1)
    e2,v2 = eigh(A2)
    
    print("\n---------------------------------------------------------------")
    print("ENERGY EIGENVALUES FOR l = 1 :-")
    print("---------------------------------------------------------------\n")
    for i in range(10):
        print("\nEnergy eigenvalue for n = "+str(i+1)+" and l = 1 is :",e1[i],"\n")
    print("------------------------------------------------------------------\n")
    
    print("\n---------------------------------------------------------------")
    print("ENERGY EIGENVALUES FOR l = 2 :-")
    print("---------------------------------------------------------------\n")
    for i in range(10):
        print("\nEnergy eigenvalue for n = "+str(i+1)+" and l = 2 is :",e2[i],"\n")
    print("------------------------------------------------------------------")
    
    # PLOT OF V_EFF AND V
    
    for i in range(1,4):
        v_eff,v_nor = effective_potential(x_new,i)
        V_eff_list.append(v_eff) ; V_nor_list.append(v_nor)

        
    for i in range(len(V_eff_list)):
        ax1.plot(x_new,V_eff_list[i],label = "Veff for l ="+str(i))
        ax1.plot(x_new,V_nor_list[i],label = "V for l ="+str(i))
    ax1.set(xlabel = "r",ylabel = "$V_{eff}$ or V")
    ax1.grid(ls = "--")
    ax1.legend()
    plt.show()
    
    
    # PLOT OF FIRST FOUR WAVE FUNCTION FOR l = 0
    
    n = 500
    A,x_new = finite_difference(10**(-14),50,n,0)
    
    e,v = eigh(A)
    
    #  ANALYTIC WAVE FUNCTION

    anal_0 = [] ; anal_1 = [] ; anal_2 = []
    for i in range(1,5):
        R0 = analytic(x_new,i,0)
        R1 = analytic(x_new,i,1)
        R2 = analytic(x_new,i,2)
        anal_0.append(R0)
        anal_1.append(R1)
        anal_2.append(R2)

    
    fig1,ax1 = plt.subplots()
    fig2,ax2 = plt.subplots()
    fig3,ax3 = plt.subplots()
    fig4,ax4 = plt.subplots()
    
    ax1.plot(x_new,normalisation(v[:,0],x_new),lw = 2,c = "b",label = "Calculated")
    ax1.scatter(x_new,x_new*anal_0[0],s = 10,c = "r",label = "Analytic")
    ax2.plot(x_new,-normalisation(v[:,1],x_new),lw = 2,c = "b",label = "Calculated")
    ax2.scatter(x_new,x_new*anal_0[1],s = 10,c = "r",label = "Analytic")
    ax3.plot(x_new,normalisation(v[:,2],x_new),lw = 2,c = "b",label = "Calculated")
    ax3.scatter(x_new,x_new*anal_0[2],s = 10,c = "r",label = "Analytic")
    ax4.plot(x_new,-normalisation(v[:,3],x_new),lw = 2,c = "b",label = "Calculated")
    ax4.scatter(x_new,x_new*anal_0[3],s = 10,c = "r",label = "Analytic")
    ax1.set(xlabel = "r",ylabel = "$rR_{nl}$",title = "Wave function for n = 1 and l = 0")
    ax2.set(xlabel = "r",ylabel = "$rR_{nl}$",title = "Wave function for n = 2 and l = 0")
    ax3.set(xlabel = "r",ylabel = "$rR_{nl}$",title = "Wave function for n = 3 and l = 0")
    ax4.set(xlabel = "r",ylabel = "$rR_{nl}$",title = "Wave function for n = 4 and l = 0")
    ax1.grid(ls = "--")
    ax2.grid(ls = "--")
    ax3.grid(ls = "--")
    ax4.grid(ls = "--")
    ax1.legend()
    ax2.legend()
    ax3.legend()
    ax4.legend()
    plt.show()
    
    # PLOT OF RADIAL PROBABILITY DENSITY 
    
    A1,x_new = finite_difference(10**(-14),50,n,1)
    A2,x_new = finite_difference(10**(-14),50,n,2)
    
    e1,v1 = eigh(A1)
    e2,v2 = eigh(A2)
    
    fig3,ax3 = plt.subplots()
    fig4,ax4 = plt.subplots()
    fig5,ax5 = plt.subplots()
    
    ax3.plot(x_new,(normalisation(v[:,0],x_new))**2,lw = 2,c = "b",label = "Calculated for n = 1 and l = 0")
    ax3.scatter(x_new,(x_new*anal_0[0])**2,s = 15,c = "r",label = "Analytic for n = 1 and l = 0")
    ax4.plot(x_new,(normalisation(v[:,1],x_new))**2,lw = 2,c = "b",label = "Calculated for n = 2 and l = 0")
    ax4.scatter(x_new,(x_new*anal_0[1])**2,s = 15,c = "r",label = "Analytic for n = 2 and l = 0")
    ax4.plot(x_new,(normalisation(v1[:,0],x_new))**2,lw = 2,c = "maroon",label = "Calculated for n = 2 and l = 1")
    ax4.scatter(x_new,(x_new*anal_1[1])**2,s = 15,c = "green",label = "Analytic for n = 2 and l = 1")
    ax5.plot(x_new,(normalisation(v[:,2],x_new))**2,lw = 2,c = "b",label = "Calculated for n = 3 and l = 0")
    ax5.scatter(x_new,(x_new*anal_0[2])**2,s = 15,c = "r",label = "Analytic for n = 3 and l = 0")
    ax5.plot(x_new,(normalisation(v1[:,1],x_new))**2,lw = 2,c = "maroon",label = "Calculated for n = 3 and l = 1")
    ax5.scatter(x_new,(x_new*anal_1[2])**2,s = 15,c = "orange",label = "Analytic for n = 3 and l = 1")
    ax5.plot(x_new,(normalisation(v2[:,0],x_new))**2,lw = 2,c = "black",label = "Calculated for n = 3 and l = 2")
    ax5.scatter(x_new,(x_new*anal_2[2])**2,s = 15,c = "gold",label = "Analytic for n = 3 and l = 2")
    ax3.set(xlabel = "r",ylabel = "$(rR_{nl})^2$",title = "Probability density for n = 1")
    ax4.set(xlabel = "r",ylabel = "$(rR_{nl})^2$",title = "Probability density for n = 2")
    ax5.set(xlabel = "r",ylabel = "$(rR_{nl})^2$",title = "Probability density for n = 3")
    ax3.grid(ls = "--")
    ax4.grid(ls = "--")
    ax5.grid(ls = "--")
    ax3.legend()
    ax4.legend()
    ax5.legend()
    plt.show()
    
    
    
    
    
    
import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigh
from scipy.integrate import simps

def finite_difference(a,b,n):
    x_array = np.linspace(a,b,n)
    h = x_array[1] - x_array[0]
    K = np.zeros((n,n))
    v = np.zeros((n,n)) 

    for i in range(n):
        for j in range(n):

            if i == j:
                v[i,j] = 0 

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

    return A,x_array

def normalisation(array,x_array):
	A = 1/np.sqrt((simps((array)**2,x_array)))
	u_norm = A*array
	return u_norm

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

if __name__ == "__main__":
    a = -1/2 ; b = 1/2 ; n = 100
    states = 10

    A,x_array = finite_difference(a,b,n)
    
    e,v = eigh(A)
    

    for i in range(states):
        plt.plot(x_array,normalisation(v[:,i],x_array))
        plt.xlabel("x")
        plt.ylabel("$\Psi$(x)")
        plt.title("Wave Function for n ="+str(i+1))
        plt.grid(ls = "--")
        plt.savefig(f"A8_{i+1}.png")
        plt.show()
    
    fig1,ax1 = plt.subplots(2,2)
    fig2,ax2 = plt.subplots(2,2) 
        
    anal = [analytic_sol(1,x_array,i) for i in range(1,5)]
        
    ax1[0,0].plot(x_array,anal[0], label='analytic solution')
    ax1[0,0].scatter(x_array,normalisation(v[:,0],x_array),label = "Finite difference method",s =10,c = "r")
    ax1[0,1].plot(x_array,anal[1], label='analytic solution')
    ax1[0,1].scatter(x_array,normalisation(v[:,1],x_array),label = "Finite difference method",s =10,c = "r")
    ax1[1,0].plot(x_array,anal[2], label='analytic solution')
    ax1[1,0].scatter(x_array,normalisation(v[:,2],x_array),label = "Finite difference method",s =10,c = "r")
    ax1[1,1].plot(x_array,anal[3], label='analytic solution')
    ax1[1,1].scatter(x_array,normalisation(v[:,3],x_array),label = "Finite difference method",s =10,c = "r")
    
    ax2[0,0].plot(x_array,(anal[0])**2, label='analytic solution')
    ax2[0,0].scatter(x_array,(normalisation(v[:,0],x_array))**2,label = "Finite difference method",s = 10,c = "r")
    ax2[0,1].plot(x_array,(anal[1])**2, label='analytic solution')
    ax2[0,1].scatter(x_array,(normalisation(v[:,1],x_array))**2,label = "Finite difference method",s = 10,c = "r")
    ax2[1,0].plot(x_array,(anal[2])**2, label='analytic solution')
    ax2[1,0].scatter(x_array,(normalisation(v[:,2],x_array))**2,label = "Finite difference method",s = 10,c = "r")
    ax2[1,1].plot(x_array,(anal[3])**2, label='analytic solution')
    ax2[1,1].scatter(x_array,(normalisation(v[:,3],x_array))**2,label = "Finite difference method",s = 10,c = "r")
        
    for i in range(2):
        for j in range(2):
            ax1[i,j].set(xlabel = "x",ylabel = "$\Psi$(x)")
            ax2[i,j].set(xlabel = "x",ylabel = "|$\Psi$(x)|$^{2}$")
            ax1[i,j].grid(ls = "--")
            ax2[i,j].grid(ls = "--")
            ax1[i,j].legend()
            ax2[i,j].legend()
    fig1.suptitle("Wave functions for n = 1,2,3,4")
    fig2.suptitle("Probability density for n = 1,2,3,4")
    plt.show()
        
    

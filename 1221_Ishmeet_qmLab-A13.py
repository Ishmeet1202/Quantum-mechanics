import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import simps
from scipy.linalg import eigh
import pandas as pd

def V(x,alpha):
    x = x[1:-1]
    v = np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            if i == j:
                v[i,j] = x[i]**2 + (2/3)*(alpha*(x[i])**3)

    return v,x

def finite_difference(a,b,n,alpha):
    x_array = np.linspace(a,b,n+2)
    h = x_array[1] - x_array[0]
    
    K = np.zeros((n,n)) 

    v,x_new = V(x_array,alpha)

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

def graph(type,n,alpha):
    v_list = [] ; e_glist = []
    fig1,ax1 = plt.subplots()
    for i in alpha:
        A,x_new = finite_difference(-10,10,200,i)
        e,v = eigh(A)
        v_list.append(v[:,n])
        e_glist.append(e[n])
    if type == 'pd':
        for i in range(len(alpha)):
            ax1.plot(x_new,(normalisation(v_list[i],x_new))**2,label = "$\\alpha$ = "+str(alpha[i]))
            ax1.set(xlabel = "x",ylabel = "|$\psi$(x)|$^{2}$",title = "Probability density for n = "+str(n))
            ax1.grid(ls = "--")
            ax1.legend()
        plt.show()
    else:
        for i in range(len(alpha)):
            ax1.plot(x_new,normalisation(v_list[i],x_new),label = "$\\alpha$ = "+str(alpha[i]))
            ax1.set(xlabel = "x",ylabel = "$\psi$(x)",title = "Wave function for n = "+str(n))
            ax1.legend()
            ax1.grid(ls = "--")
        plt.show()
    
if __name__ == "__main__":
    a = -10
    b = 10
    n = 200
    alpha = [0,10**(-0),10**(-1),10**(-2),10**(-3),10**(-4)]
        
    x = np.linspace(a,b,n)
    for i in alpha:
        v = x**2 + (2/3)*(i*(x)**3)
        plt.plot(x,v,label='alpha= '+str(i))
    plt.grid()
    plt.legend()
    plt.xlabel('x')
    plt.ylabel('V(x)')
    plt.title('Potential VS x')
    
    for j in range(11):
        En_calc = []
        En_analy = []
        alpha_list = []
        for i in alpha:
            A,x_new = finite_difference(a,b,n,i)
            e,v = eigh(A)
            analytic= (2*j + 1) - (1/8)*(i**2)*(15*(2*j + 1)**2 + 7)
            En_calc.append(e[j])
            alpha_list.append(i)
            En_analy.append(analytic)
        print("\n-----------------------------------------------------------")
        print('Energy eigen values For n =',j)
        print("-----------------------------------------------------------\n")
        print(pd.DataFrame({'Alpha':alpha,'Calculated':En_calc,'Analytic': En_analy}))
        print("-----------------------------------------------------------\n")

    states = [i for i in range(11)]
    En=[]
    for j in alpha:
        E=[]
        for i in range(11):
            A,x_new = finite_difference(a,b,n,j)
            e,v = eigh(A)
            E.append(e[i])
        En.append(E)
    
    fig, axes = plt.subplots(2,3, figsize=(15,10))
    
    axes[0][0].plot(states,En[0],label = "$\\alpha$ = 0")
    axes[0][1].plot(states,En[1],label = "$\\alpha$ = 1")
    axes[0][2].plot(states,En[2],label = "$\\alpha$ = $10^{-1}$")
    axes[1][0].plot(states,En[3],label = "$\\alpha$ = $10^{-2}$")
    axes[1][1].plot(states,En[4],label = "$\\alpha$ = $10^{-3}$")
    axes[1][2].plot(states,En[5],label = "$\\alpha$ = $10^{-4}$")
    
    for i in range(2):
        for j in range(3):
            axes[i][j].set(xlabel = "n(state)",ylabel = "Energy")
            axes[i][j].grid(ls = "--")
            axes[i][j].legend()
    fig.suptitle("n(States) Vs Energy of state at different $\\alpha$")
    plt.show()
        
    alpha = [0,0.01,0.1]
    
    for i in range(5):
        graph('wf',i,alpha)
        graph('pd',i,alpha)
    
    
   
        
            
    
    
        
        
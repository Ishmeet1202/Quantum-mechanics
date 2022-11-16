import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simps
from scipy.special import hermite
import pandas as pd

def V(x):
    return (x**2)/2

def numerov_method(x_array,ic,n,delta_e):
    u0 = ic[0] ; u1 = ic[1]
    h = (x_array[1] - x_array[0])
    V_array = [V(i) for i in x_array]
    V_array = np.array(V_array)
    
    if n%2 != 0:
        e = n + 0.5 + delta_e
        alpha = 2*(-V_array + e)
        ci = 1 + (((h ** 2) / 12) * alpha)
        u = [u0,u1]
        for i in range(1,len(x_array)-1):
            u2 = (((12 - 10*ci[i]) * u[i]) - (ci[i-1] * u[i-1])) / (ci[i+1])
            u.append(u2)
        
        # NORMALIZATION
        A = 1/np.sqrt((simps((np.array(u))**2,x_array)))
        u_norm = A*np.array(u)
    
    else:
        e = n + 0.5 + delta_e
        alpha = 2*(-V_array + e)
        ci = 1 + (((h ** 2) / 12) * alpha)
        u = [u0,(u0*(12 - 10*ci[0]))/2*ci[1]]
        for i in range(1,len(x_array)-1):
            u2 = (((12 - 10*ci[i]) * u[i]) - (ci[i-1] * u[i-1])) / (ci[i+1])
            u.append(u2)
            
        # NORMALIZATION
        A = 1/np.sqrt((simps((np.array(u))**2,x_array)))
        u_norm = A*np.array(u)
            
    return u_norm

def parity(array,n):
    u1 = array
    x_array = np.linspace(-x_max,x_max,(N*2)-1)
    if n%2 == 0:
        u2 = u1
        #print(u2)
        u2 = np.delete(u2,0)
        u2 = u2[::-1]
    else:
        u2 = (-1)*u1
        u2 = np.delete(u2,0)
        u2 = u2[::-1]
            
    u = np.concatenate((u2,u1),axis = None)
    
    return x_array,u

def analytic_solution(x,n):
    Hn = hermite(n)
    psi = (np.exp(-((x**2)/2))) * Hn(x)
    
    # NORMALISATION
    A = 1/np.sqrt((simps((psi)**2,x)))
    psi_norm = A*psi
    
    return psi_norm

def energy(n,delta_e,omega):
    # CALCULATED ENERGY
    ec = [] ; ea = []
    hbar = 1.055 * 10**(-34)
    for i in range(len(n)):
        e_c = (n[i] + 0.5 + delta_e)*hbar*omega
        E_c = e_c / 1.6 * 10**(-19)
        ec.append(E_c)
    for i in range(len(n)):
        e_a = (n[i] + 0.5)*hbar*omega
        E_a = e_a / 1.6  * 10**(-19)
        ea.append(E_a)
        
    return ec,ea

if __name__ == "__main__":
    x_min = 0
    x_max = 4
    ic_o = [0,0.001]
    ic_e = [1,0]
    n = [0,1,2,3]
    delta_e = [10**(-2),10**(-4),10**(-6),10**(-8)]
    N = 2000
    x_array = np.linspace(x_min,x_max,N)
    u1_array = [] ; u2_array = []
    
    # PART A
    
    for i in range(len(delta_e)):
        ic = ic_e
        ui = numerov_method(x_array,ic,0,delta_e[i])
        x,uf = parity(ui,0)
        u1_array.append(uf)
        
            
    for i in range(len(n)):
        if n[i]%2 == 0:
            ic = ic_e
            ui = numerov_method(x_array,ic,n[i],10**(-6))
            x,uf = parity(ui,n[i])
            if n[i] == 1 or n[i] == 0:
                uf = uf
            else:
                uf = -uf
            u2_array.append(uf)
        
        else:
            ic = ic_o
            ui = numerov_method(x_array,ic,n[i],10**(-6))
            x,uf = parity(ui,n[i])
            if n[i] == 1 or n[i] == 0:
                uf = uf
            else:
                uf = -uf
            u2_array.append(uf)
        
    u1_array = np.array(u1_array)
    u2_array = np.array(u2_array)
    
    psi_1 = analytic_solution(x,0)
    psi_2 = analytic_solution(x,1)
    psi_3 = analytic_solution(x,2)
    psi_4 = analytic_solution(x,3)
    
    fig4,ax4 = plt.subplots(2,2)
    ax4[0,0].plot(x,u1_array[0],label = "Numerov method")
    ax4[0,0].plot(x,psi_1,label = "Analytic solution")
    ax4[0,1].plot(x,u1_array[1],label = "Numerov method")
    ax4[0,1].plot(x,psi_1,label = "Analytic solution")
    ax4[1,0].plot(x,u1_array[2],label = "Numerov method")
    ax4[1,0].plot(x,psi_1,label = "Analytic solution")
    ax4[1,1].plot(x,u1_array[3],label = "Numerov method")
    ax4[1,1].plot(x,psi_1,label = "Analytic solution")
    ax4[0,0].set(xlabel = "x",ylabel = "|$\Psi$(x)|",title = "Wave function for n = 0 at $\Delta$e = 10$^{-2}$")
    ax4[0,1].set(xlabel = "x",ylabel = "|$\Psi$(x)|",title = "Wave function for n = 0 at $\Delta$e = 10$^{-4}$")
    ax4[1,0].set(xlabel = "x",ylabel = "|$\Psi$(x)|",title = "Wave function for n = 0 at $\Delta$e = 10$^{-6}$")
    ax4[1,1].set(xlabel = "x",ylabel = "|$\Psi$(x)|",title = "Wave function for n = 0 at $\Delta$e = 10$^{-8}$")
    for i in range(2):
        for j in range(2):
            ax4[i,j].grid(ls = "--")
            ax4[i,j].legend()
    plt.show()
       
    fig1,ax1 = plt.subplots()
    fig2,ax2 = plt.subplots()
    fig3,ax3 = plt.subplots()
    
    ax1.plot(x,u2_array[1],label = "Numerov method")
    ax1.plot(x,psi_2,label = "Analytic solution")
    ax2.plot(x,u2_array[2],label = "Numerov method")
    ax2.plot(x,psi_3,label = "Analytic solution")
    ax3.plot(x,u2_array[3],label = "Numerov method")
    ax3.plot(x,psi_4,label = "Analytic solution")
    ax1.set(xlabel = "x",ylabel = "$\Psi$(x)",title = "Wave function for n = 1 at $\Delta$ = 10$^{-6}$")
    ax2.set(xlabel = "x",ylabel = "$\Psi$(x)",title = "Wave function for n = 2 at $\Delta$ = 10$^{-6}$")
    ax3.set(xlabel = "x",ylabel = "$\Psi$(x)",title = "Wave function for n = 3 at $\Delta$ = 10$^{-6}$")
    ax1.grid(ls = "--")
    ax2.grid(ls = "--")
    ax3.grid(ls = "--")
    ax1.legend()
    ax2.legend()
    ax3.legend()
    plt.show()
    
    # B PART
    
    u1 = (u2_array[0])**2 ; u2 = (u2_array[1])**2 ; u3 = (u2_array[2])**2 ; u4 = (u2_array[3])**2
    
    fig1,ax1 = plt.subplots()
    fig2,ax2 = plt.subplots()
    fig3,ax3 = plt.subplots()
    fig4,ax4 = plt.subplots()
    
    ax1.plot(x,u1,label = "Numerov method")
    ax1.plot(x,(psi_1)**2,label = "Analytic solution")
    ax2.plot(x,u2,label = "Numerov method")
    ax2.plot(x,(psi_2)**2,label = "Analytic solution")
    ax3.plot(x,u3,label = "Numerov method")
    ax3.plot(x,(psi_3)**2,label = "Analytic solution")
    ax4.plot(x,u4,label = "Numerov method")
    ax4.plot(x,(psi_4)**2,label = "Analytic solution")
    ax1.set(xlabel = "x",ylabel = "|$\Psi$(x)|$^{2}$",title = "Probability density curve for n = 0")
    ax2.set(xlabel = "x",ylabel = "|$\Psi$(x)|$^{2}$",title = "Probability density curve for n = 1")
    ax3.set(xlabel = "x",ylabel = "|$\Psi$(x)|$^{2}$",title = "Probability density curve for n = 2")
    ax4.set(xlabel = "x",ylabel = "|$\Psi$(x)|$^{2}$",title = "Probability density curve for n = 3")
    ax1.grid(ls = "--")
    ax2.grid(ls = "--")
    ax3.grid(ls = "--")
    ax4.grid(ls = "--")
    ax1.legend()
    ax2.legend()
    ax3.legend()
    ax4.legend()
    plt.show()
    
    # PART C
    
    ec,ea = energy(n,10**(-6),5.5 * 10**(14))
    data = {'Energy state (n)':n,'Calculated energy (in eV)':ec,'Analytic energy (in eV)':ea}
    print("\n",pd.DataFrame(data),"\n")
    
    
    
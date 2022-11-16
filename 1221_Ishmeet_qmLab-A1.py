import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def energy(n):
    en = -13.6/n**2
    return en

def frequency_classical(n):
    f_cln = (m*(e**4))/(32*(np.pi**3)*(epsilon**2)*(hbar**3)*(n**3))
    return f_cln

def frequency_quantum(n):
    f_qn = (m*(e**4)*(2*n -1))/(64*(np.pi**3)*(epsilon**2)*(hbar**3)*((n**2)*(n-1)**2))
    return f_qn

def energy_graph():
    E = [] ;  n = []
    for j in range(1,11):
        E.append(energy(j))
        n.append(j)
    for i in range(len(E)):
        plt.axhline(y=E[i], xmin=0.1, xmax=0.8, linewidth=2, label='n= '+ str(n[i]),c = "r")
        plt.ylabel('Energy (in eV)')
        plt.legend()
    plt.title("Energy level plot")
    plt.grid(ls = "--")
    plt.show()

def tolerance():
    p = 0.5
    n_list = []
    error_list = []
    freq_cl = []
    freq_qn = []
    while True:
        n = 10**(p)
        error = np.fabs((frequency_classical(n)- frequency_quantum(n))/frequency_quantum(n))
        if error < 10**(-5):
            error_list.append(error)
            n_list.append(n)
            freq_cl.append(frequency_classical(n))
            freq_qn.append(frequency_quantum(n))
            break
        else:
            error_list.append(error)
            n_list.append(n)
            freq_cl.append(frequency_classical(n))
            freq_qn.append(frequency_quantum(n))
            p = p + 0.5
            
    return np.array(n_list),np.array(error_list),np.array(freq_cl),np.array(freq_qn)

if __name__ == "__main__":
    
    epsilon = 8.85 * 10**(-12)
    m = 9.1 * 10**(-31)
    e = -1.6 * 10**(-19)
    hbar = 1.054 * 10**(-34)
    
    sol = tolerance()
    
    n_list = np.log(sol[0])
    error_list =  sol[1]
    freq_cl = sol[2]
    freq_qn = sol[3]
    
    data = {'n':sol[0],'Frequency classical':sol[2],'Frequency quantum':sol[3],'Relative Error':sol[1]}
    print(pd.DataFrame(data))
    
    energy_graph()
    
    plt.plot(n_list,error_list)
    plt.xlabel("ln(n)")
    plt.ylabel("Relative error")
    plt.title("RELATIVE ERROR Vs ln(n)")
    plt.grid(ls = "--")
    plt.show()
    
    
    
    

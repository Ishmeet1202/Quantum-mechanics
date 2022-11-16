import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simps

def V(x):
    return (x**2)/2

def numerov_method(x_array,ic,nodes,e):
    u0 = ic[0] ; u1 = ic[1]
    h = (x_array[1] - x_array[0])
    V_array = [V(i) for i in x_array]
    V_array = np.array(V_array)
    
    if nodes%2 != 0:
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

def energy_eigenvalue(u_in,nodes,e_max,e_min):
    while np.fabs(e_max - e_min) >= 10**(-10):
        count = 0
        for i in range(len(u_in)-1):
            if u_in[i]*u_in[i+1] < 0:
                count = count + 1
                print(count)
        if count <= int((nodes)/2):
            e_min = (e_max + e_min)/2
        else:
            e_max = (e_max + e_min)/2

        e_new = (e_min + e_max)/2

        if nodes%2  == 0:
            u_in = numerov_method(x_array,ic_e,nodes,e_new)
        else:
            u_in = numerov_method(x_array,ic_o,nodes,e_new)

    return u_in,e_new

def classical_tp(n):
    return np.sqrt(2*n +1)
        
if __name__ == "__main__":
    x_max = 4 ; x_min = 0
    N = 2000
    x_array = np.linspace(x_min,x_max,N)
    e_max = 10 ; e_min = 0
    ic_o = [0,-(x_array[1]-x_array[0])]
    ic_e = [1,0]
    nodes = 0

    if nodes%2  == 0:
        u_in = numerov_method(x_array,ic_e,nodes,e_max)
    else:
        u_in = numerov_method(x_array,ic_o,nodes,e_max)

    u_in, e_new = energy_eigenvalue(u_in,nodes,e_max,e_min)

    print(e_new)

    x,u_f = parity(u_in,nodes)

    plt.plot(x,u_f)
    #plt.plot(x_array,u_in)
    plt.show()


import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import simps
from scipy.linalg import eigh

def V(x,l,alpha):
    x = x[1:-1]
    v = np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            if i == j:
                v[i,j] = ((-2)/(x[i]))*(np.exp(-x[i]/alpha))+(l*(l+1))/((x[i])**2)

    return v,x

def finite_difference(a,b,n,l,alpha):
    x_array = np.linspace(a,b,n+2)
    h = x_array[1] - x_array[0]
    
    K = np.zeros((n,n)) 

    v,x_new = V(x_array,l,alpha)

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

if __name__ == "__main__":
    a = 10**(-14)
    b = 10
    n = 200
    v_list = [] ; e_glist = []
    alpha = [2,5,10,20,100]
    
    fig1,ax1 = plt.subplots()
    fig2,ax2 = plt.subplots()
    fig3,ax3 = plt.subplots()
    
    print("\n----------------------------------------------------------------------------")
    print("For n = 1:")
    print("----------------------------------------------------------------------------\n")
    
    for i in alpha:
        A,x_new = finite_difference(a,b,n,0,i)
        
        e,v = eigh(A)
        
        print("Bound state energy eigen value exists for alpha =",i,": ",e[0],"\n")
        
        v_list.append(v[:,0])
        e_glist.append(e[0])
    
    for i in range(len(alpha)):    
        ax1.plot(x_new,normalisation(v_list[i],x_new),label = "$\\alpha$ ="+str(alpha[i]))
        ax2.plot(x_new,(normalisation(v_list[i],x_new))**2,label = "$\\alpha$ ="+str(alpha[i]))
    ax1.set(xlabel = "x",ylabel = "$\psi$(x)",title = "Wave function for n = 1 and l = 0")
    ax2.set(xlabel = "x",ylabel = "|$\psi$(x)|$^{2}$",title = "Probability density for n = 1 and l = 0")
    ax3.plot(alpha,e_glist,c = "b")
    ax3.scatter(alpha,e_glist,c = "r")
    ax3.set(xlabel = "$\\alpha$",ylabel = "Ground state energy at different $\\alpha$",title = "Ground state energy Vs $\\alpha$")
    ax1.legend()
    ax2.legend()
    ax1.grid(ls = "--")
    ax2.grid(ls = "--")
    ax3.grid(ls = "--")
    plt.show()
    
    
    
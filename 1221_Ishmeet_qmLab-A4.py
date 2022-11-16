import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate

# initial conitions:
u_values = [0, 0.001]


def V(x):
    V = 0
    return V

def numerov(a,b,e,i_cs, n, V):

    x = np.linspace(a, b, n)
    u0, u1 = i_cs[0], i_cs[1]

    Vx = []
    for i in x:
        Vx.append(V(i))
    Vx = np.array(Vx)
    alpha = -(Vx - e)
    #alpha=np.ones(len(x))*50
    h=x[1]-x[0]
    print("h is:", h)
    ci = 1 + (((h ** 2) / 12) * alpha)
    #print(ci)
    p = 0
    u = [u0, u1]
    for i in range(1,len(x)-1):
        p+=1
        u2 = (((12 - 10*ci[i]) * u[i]) - (ci[i - 1] * u[i-1])) / (ci[i + 1])
        u.append(u2)
        if p==len(x)-2:
            break
    return x,np.array(u)

#print(numerov(200, u_values, 100))
plt.plot(np.linspace(-1/2, 1/2, 500),numerov(-1/2,1/2,40, u_values, 500,V)[1])
plt.xlabel('x')
plt.ylabel('u')
plt.grid()
plt.show()


#EQ. 3

u_values = [1, 1.001]

def ics(n):
    x=np.linspace(0,5,n)
    h=x[1]-x[0]
    u_0=1
    u_1= 1+ ((h**2)/2) + 3*((h**4)/24)

    return [u_0,u_1]

def V2(x):
    V = x**2
    return V

plt.plot(np.linspace(0, 5, 4),numerov(0,5,-1,ics(4),4,V2)[1])
plt.title('u vs x for N=4')
plt.xlabel('x')
plt.ylabel('u')
plt.grid()
plt.show()

for i in range(1,7):
    plt.plot(np.linspace(0, 5, 2**i), numerov(0,5,-1, ics(2**i), 2**i, V)[1], label='for grid points= '+str(2**i))
plt.legend()
plt.title('u VS x Fofr Different Values Of k')
plt.xlabel('x')
plt.ylabel('u')
plt.grid()
plt.show()




x=numerov(0,1,-1,ics(100),20,V2)[0]
u=[1,0]

def f(u,x):
    return (u[1], (1+x**2)*u[0])


y0=[1,0]
us=integrate.odeint(f,y0,x)

ys=us[:,0]

plt.plot(x,numerov(0,1,-1,ics(100),20,V2)[1],label='calculated value by numerov')
plt.scatter(x, ys, c='g', label='value by scipy inbuilt func.')
plt.grid()
plt.xlabel('x')
plt.ylabel('u')
plt.title('Comparing Values Calculated By Numerov & Scipy Func.')
plt.legend()
plt.show()



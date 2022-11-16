import numpy as np
import scipy.integrate as integrate
import matplotlib.pyplot as plt

def a_k(b,x,k):
    A=1/np.sqrt(2*b)
    a=(1/np.sqrt(2*np.pi))*integrate.simps(A*np.exp(-1j*k*x),x)
    return a

def psi(a,b,x,k,t):
    #a=a_k(b,x,k)
    p=[]
    for i in x:
        psi=(1/np.sqrt(2*np.pi))*integrate.simps(a*np.exp(1j*(k*i - (k**2)*t)),k)
        p.append(psi)
    return p
    
if __name__ == "__main__":

	k=np.linspace(-10,10,200)
	x=np.linspace(-1,1,200)    #Put LImit according to b
	A=[]

	#A- ii
	for i in k:
		a=a_k(1,x,i)
		A.append(a)
		plt.scatter(i,a,c='b')
	plt.ylabel('Momentum')
	plt.xlabel('x')
	plt.title("For t=0")
	plt.grid()
	plt.show()

	#A-i

	q=psi(np.array(A),1,x,k,0)
	plt.plot(x,np.array(q)**2)
	plt.ylabel('Ψ square')
	plt.xlabel('x')
	plt.title("Probability Density For t=0")
	plt.grid()
	plt.show()

	#B
	t=np.arange(0,2,0.1)
	for i in t:
		q=psi(np.array(A),1,x,k,i)
		plt.plot(x,np.array(q)**2)
	plt.xlabel('x')
	plt.ylabel('Ψ square')
	plt.title("probability Density For Different Time Intervals")
	plt.grid()
	plt.show()

	#C
	x=np.linspace(-1/2,1/2,200)
	for i in t:
		q=psi(np.array(A),1,x,k,i)
		plt.plot(x,np.array(q)**2)
	plt.legend()
	plt.xlabel('x')
	plt.ylabel('Ψ square')
	plt.title("probability Density For |x|<b/2")
	plt.grid()
	plt.show()




#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate
from scipy.special import eval_hermite
import pandas as pd
from scipy.stats import linregress
from scipy.optimize import fsolve


def potential_dimless(x_arr,x0=None):  
    if x0==None:
        return 0.5*np.square(x_arr)
    else:
        return 0.5*np.square(x_arr-(x0*np.ones(len(x_arr))))
    
def normalization_factor(y_arr,x_arr):  #function to normalise wavefunction
    new_arr=np.square(y_arr)
    k=scipy.integrate.simps(new_arr,x_arr)
    c=np.sqrt(1/k)
    new_arr2=c*y_arr
    return new_arr2,c

def analytic_HO(x_arr,n):               #analytic solution of Harmonic Oscillator 
    ans_arr=np.zeros((len(x_arr)))
    for i in range(len(x_arr)):
        ans_arr[i]=eval_hermite(n,x_arr[i])*np.exp(-0.5*(x_arr[i]**2))
    ans_arr=normalization_factor(ans_arr,x_arr)[0]
    return ans_arr

def numerov_method_HO(u0,x_arr,alpha,e,du0=1): #u0 is value at first point x_arr is input array alpha is differential eqn function
    #e is energy of the particle du0 is rise of the function 1=positive derivative -1 will mean negative derivative
    h=abs(x_arr[1]-x_arr[0])
    u_arr=np.zeros((len(x_arr)))
    u_arr[0]=u0
    c_arr=np.ones((len(x_arr)))+np.multiply(alpha(x_arr,e),((h**2)/12))
    if u0==0:
        if du0==1:
            u_arr[1]=h
        else:
            u_arr[1]=-h
    else:
        u_arr[1]=((6-5*c_arr[1])*u_arr[0])/c_arr[1]
    for i in range(1,len(x_arr)-1,1):
        u_arr[i+1]=(1/c_arr[i+1])*(((12-10*c_arr[i])*u_arr[i])-(c_arr[i-1]*u_arr[i-1]))
    return u_arr

def alpha(x_arr,e):
    return 2*(e*np.ones(len(x_arr))-potential_dimless(x_arr))

def parity_operator(u_arr,key):
    u_arr2=np.zeros((len(u_arr)-1))
    if key==0:     #even
        u_arr2=u_arr[1:]
    else:          #odd
        u_arr2=-1*u_arr[1:]
    return np.append(u_arr2[::-1],u_arr)

def nodes_finder(u_arr,key=None):
    count=0
    if key==None:        #find number of node
        for i in range(len(u_arr)-1):
            if (u_arr[i]*u_arr[i+1])<0:
                count=count+1
        return count
    else:                #find the first point at which function changed from positive to negative
        for i in range(len(u_arr)-1):    
            if (u_arr[i]*u_arr[i+1])<0:
                return i      #used for finding classical turning points
            
def matching_function(u_arr1,u_arr2):  #matches forward and backward array
    u_l=u_arr1[-1]
    u_r=u_arr2[0]
    c=(u_r/u_l)
    u_n=np.multiply(u_arr1,c)
    return np.append(u_n,u_arr2[1:])

def derivative(p_arr,h,key=None,index=None):#calculates the derivative of given array
    p2_arr=np.zeros(len(p_arr)) 
    for i in range(len(p_arr)):
        if i==0:                #forward difference for first point
            p2_arr[i]=(p_arr[1]-p_arr[0])/h
        elif i==(len(p_arr)-1): #backward difference for last point
            p2_arr[i]=(p_arr[-1]-p_arr[-2])/h
        else:                   
            if key==None:
                p2_arr[i]=(p_arr[i+1]-p_arr[i])/(h) #forward difference
            else:
                p2_arr[i]=(p_arr[i]-p_arr[i-1])/(h) #backward difference
    if index==None:   #returns full arrat
            return p2_arr
    else:             #returns value of derivative at a particular index
        return p2_arr[index]
    
def phi(e,u0,du0,cl):
    h=x_arr[1]-x_arr[0]                                    #step size
    u_arr_1=numerov_method_HO(u0,x_arr_1,alpha,e,du0)      #integrating forward
    u_arr_2=numerov_method_HO(0,np.flip(x_arr_2),alpha,e)  #integrating backward
    u_arr=matching_function(u_arr_1,np.flip(u_arr_2))      #combined answer after rescaling
    u_arr=np.sqrt(2)*normalization_factor(u_arr,x_arr)[0]  #normalising wavefunction
    duL=derivative(u_arr,h,index=cl)                       #calculates fwd difference at classical turning point
    duR=derivative(u_arr,h,key=1,index=cl)                 #calculates backward difference at classical t.p
    return (duL-duR)                                       #difference in derivative

def phi2(e,u0,du0,cl):
    h=x_arr[1]-x_arr[0]          #step size
    u_arr_1=numerov_method_HO(u0,x_arr_1,alpha,e,du0)     #integrating forward
    u_arr_2=numerov_method_HO(0,np.flip(x_arr_2),alpha,e) #integrating backward
    u_arr=matching_function(u_arr_1,np.flip(u_arr_2))     #combined answer after rescaling
    u_arr=np.sqrt(2)*normalization_factor(u_arr,x_arr)[0]  #normalising wavefunction
    c_arr=np.ones((len(x_arr)))+np.multiply(alpha(x_arr,e),((h**2)/12)) 
    ans=(1/h)*(u_arr[cl+1]+u_arr[cl-1]-((12*c_arr[cl]-10)*u_arr[cl])) #calculates difference in derivative(indirect method)
    return (ans)

x0=0            
x_max=5         
N=40*(x_max-x0) 

x_arr=np.linspace(x0,x_max,N,float) #this array will be used for calculation
x_full=parity_operator(x_arr,1)     #this array will be used for plotting it is from -x_max to x_max

n_nodes_list=range(0,6,1)           #nodes/eigenstates at which wavefuntion will be calculated
 
eigen_val_list=np.zeros((len(n_nodes_list)))        #initialising array to store eigenvalue
u_matrix=np.zeros((len(n_nodes_list),len(x_full)))  #initialising matrix to store wavefunction(numerical)
u_num_mat=np.zeros((len(n_nodes_list),len(x_full))) #initialising matrix to store wavefunction(analytic)

tol=10**-19  #tolerance for derivative matching
tol2=10**-3  #tolerance for guess value
k=0          #counter for node list
wave_cl=np.zeros((len(n_nodes_list)))   
der_cl=np.zeros((len(n_nodes_list))) 
for n_nodes in n_nodes_list: 
    E_min=min(potential_dimless(x_arr))  
    E_max=max(potential_dimless(x_arr))
    count=1                 #counter for number of iteration
    max_iter=100              
    e=0                     #initialising energy value 
    while count<max_iter:   
        e=(E_max+E_min)/2   
        n_half=int(n_nodes/2)
        du0=0
        if n_nodes%4==0:
            u0=1
        elif n_nodes%4==1:
            u0=0
            du0=1
        elif n_nodes%4==2:
            u0=-1
        else:
            u0=0
            du0=-1
        cl=nodes_finder(alpha(x_arr,e),key=1)
        x_arr_1=x_arr[:(cl)+1]
        x_arr_2=x_arr[(cl):]
        u_arr_1=numerov_method_HO(u0,x_arr_1,alpha,e,du0)
        u_arr_2=numerov_method_HO(0,np.flip(x_arr_2),alpha,e)
        u_arr=matching_function(u_arr_1,np.flip(u_arr_2))
        n_f=nodes_finder(u_arr)
        
        if n_f<n_half or n_f==n_half:
            E_min=e
        else:
            E_max=e
        if abs(E_max-E_min)<tol2:
            e_cor=fsolve(phi2,e,(u0,du0,cl))    
            e=e_cor
            cl=nodes_finder(alpha(x_arr,e),key=1)
            x_arr_1=x_arr[:(cl)+1]
            x_arr_2=x_arr[(cl):]
            u_arr_1=numerov_method_HO(u0,x_arr_1,alpha,e,du0)
            u_arr_2=numerov_method_HO(0,np.flip(x_arr_2),alpha,e)
            u_arr=matching_function(u_arr_1,np.flip(u_arr_2))
            
            u1_arr=np.sqrt(2)*normalization_factor(u_arr,x_arr)[0]
            wave_cl[k]=u1_arr[cl]
            der_cl[k]=derivative(u1_arr,x_arr[1]-x_arr[0],index=cl)
            
            eigen_val_list[k]=e_cor
            u_matrix[k,:]=normalization_factor(parity_operator(u_arr,(k%2)),x_full)[0]
            u_num_mat[k,:]=analytic_HO(x_full,k)
            break
        count=count+1
    k=k+1
    
eigen_val_ana=np.zeros((len(eigen_val_list)))
for i in range(len(eigen_val_list)):
    eigen_val_ana[i]=i+0.5
    
table_mat=np.column_stack((n_nodes_list,eigen_val_list,eigen_val_ana))
df=pd.DataFrame(table_mat,columns=["Energy State","Energy value(numerical)","Energy value(analytic)"])
print(df)

#print(wave_cl,der_cl)
table_mat2=np.column_stack((wave_cl,der_cl))
df2=pd.DataFrame(table_mat2,columns=["Wavefn at cl","Derivative at cl"])
print(df2)

for i in range(k):
    plt.title("Wavefunctions")
    plt.scatter(x_full,u_matrix[i,:],marker=".",label="Energy state "+str(i)+"(numerical)")
    plt.plot(x_full,u_num_mat[i,:],label="Energy state "+str(i)+"(analytic)")
    plt.xlabel("x")
    plt.ylabel("u(x)")
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.grid()
plt.show()

for i in range(k):
    plt.title("Energy State= "+str(i))
    plt.scatter(x_full,u_matrix[i,:],marker=".",label="wavefunction")
    plt.scatter(x_full,np.square(u_matrix[i,:]),marker=".",label="probability density")
    i1=int((x_max-np.sqrt((2*i)+1))*(N/x_max))
    i2=int((x_max+np.sqrt((2*i)+1))*(N/x_max))
    y_c=normalization_factor(0.5*np.square(x_full[i1:i2]),x_full[i1:i2])[0]
    plt.plot(x_full[i1:i2],y_c,label="classical")
    plt.vlines(x =[-1*np.sqrt(2*i+1), np.sqrt(2*i+1)], ymin = min(y_c),ymax=max(y_c),
           colors = 'purple',linestyle=':',
           label = 'classical allowed region')
    plt.legend(loc='best')
    plt.grid()
    plt.show()

for i in range(k):
    plt.scatter(x_full,np.square(u_matrix[i,:]),marker=".",label="Energy state "+str(i)+"(numerical)")
    plt.plot(x_full,np.square(u_num_mat[i,:]),label="Energy state "+str(i)+"(analytic)")
    plt.xlabel("x")
    plt.ylabel("u(x)")
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.title("Probability Density Energy State")
plt.grid()
plt.show()


# In[14]:


x0=0
x_max=2
N=40*(x_max-x0)

x_arr=np.linspace(x0,x_max,N,float) #this array will be used for calculation
x_full=parity_operator(x_arr,1)     #this array will be used for plotting it is from -x_max to x_max

n_nodes_list=range(0,2,1)
eigen_val_list=np.zeros((len(n_nodes_list)))
u_matrix=np.zeros((len(n_nodes_list),len(x_full)))
u_num_mat=np.zeros((len(n_nodes_list),len(x_full)))

tol=10**-19
tol2=10**-3
k=0
for n_nodes in n_nodes_list: 
    E_min=min(potential_dimless(x_arr))
    E_max=max(potential_dimless(x_arr))
    count=1
    max_iter=100
    e=0
    while count<max_iter:
        e=(E_max+E_min)/2
        n_half=int(n_nodes/2)
        du0=0
        if n_nodes%4==0:
            u0=1
        elif n_nodes%4==1:
            u0=0
            du0=1
        elif n_nodes%4==2:
            u0=-1
        else:
            u0=0
            du0=-1
        cl=nodes_finder(alpha(x_arr,e),key=1)
        x_arr_1=x_arr[:(cl)+1]
        x_arr_2=x_arr[(cl):]
        u_arr_1=numerov_method_HO(u0,x_arr_1,alpha,e,du0)
        u_arr_2=numerov_method_HO(0,np.flip(x_arr_2),alpha,e)
        u_arr=matching_function(u_arr_1,np.flip(u_arr_2))
        n_f=nodes_finder(u_arr)
        
        if n_f<n_half or n_f==n_half:
            E_min=e
        else:
            E_max=e
        if abs(E_max-E_min)<tol2:
            e_cor=fsolve(phi2,e,(u0,du0,cl))    
            e=e_cor
            cl=nodes_finder(alpha(x_arr,e),key=1)
            x_arr_1=x_arr[:(cl)+1]
            x_arr_2=x_arr[(cl):]
            u_arr_1=numerov_method_HO(u0,x_arr_1,alpha,e,du0)
            u_arr_2=numerov_method_HO(0,np.flip(x_arr_2),alpha,e)
            u_arr=matching_function(u_arr_1,np.flip(u_arr_2))
            
            eigen_val_list[k]=e_cor
            u_matrix[k,:]=normalization_factor(parity_operator(u_arr,(k%2)),x_full)[0]
            u_num_mat[k,:]=analytic_HO(x_full,k)
            break
        count=count+1
    k=k+1
    
eigen_val_ana=np.zeros((len(eigen_val_list)))
for i in range(len(eigen_val_list)):
    eigen_val_ana[i]=i+0.5
    
table_mat=np.column_stack((n_nodes_list,eigen_val_list,eigen_val_ana))
df=pd.DataFrame(table_mat,columns=["Energy State","Energy value(numerical)","Energy value(analytic)"])
print(df)

for i in range(k):
    plt.title("Wavefunctions")
    plt.scatter(x_full,u_matrix[i,:],marker=".",label="Energy state "+str(i)+"(numerical)")
    plt.plot(x_full,u_num_mat[i,:],label="Energy state "+str(i)+"(analytic)")
    plt.xlabel("x")
    plt.ylabel("u(x)")
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
#plt.legend(loc='lower right')
plt.grid()
plt.show()

for i in range(k):
    title3="Energy State= "+str(i)+" for x_max=2"
    plt.title(title3)
    plt.scatter(x_full,u_matrix[i,:],marker=".",label="wavefunction")
    plt.scatter(x_full,np.square(u_matrix[i,:]),marker=".",label="probability density")
    i1=int((x_max-np.sqrt((2*i)+1))*(N/x_max))
    i2=int((x_max+np.sqrt((2*i)+1))*(N/x_max))
    y_c=normalization_factor(0.5*np.square(x_full[i1:i2]),x_full[i1:i2])[0]
    plt.plot(x_full[i1:i2],y_c,label="classical")
    plt.vlines(x =[-1*np.sqrt(2*i+1), np.sqrt(2*i+1)], ymin = min(y_c),ymax=max(y_c),
           colors = 'purple',linestyle=':',
           label = 'classical allowed region')
    plt.legend(loc='best')
    plt.grid()
    plt.savefig("C:\\Users\\vyash\\Desktop\\YASH\\quantum mechanics\\practical_qm\\A7\\"+title3)
    plt.show()

for i in range(k):
    plt.scatter(x_full,np.square(u_matrix[i,:]),marker=".",label="Energy state "+str(i)+"(numerical)")
    plt.plot(x_full,np.square(u_num_mat[i,:]),label="Energy state "+str(i)+"(analytic)")
    plt.xlabel("x")
    plt.ylabel("u(x)")
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.title("Probability Density Energy State")
plt.grid()
plt.show()


# In[12]:


x0=0
x_max=10
N=40*(x_max-x0)

x_arr=np.linspace(x0,x_max,N,float) #this array will be used for calculation
x_full=parity_operator(x_arr,1)     #this array will be used for plotting it is from -x_max to x_max

n_nodes_list=range(0,6,1)
eigen_val_list=np.zeros((len(n_nodes_list)))
u_matrix=np.zeros((len(n_nodes_list),len(x_full)))
u_num_mat=np.zeros((len(n_nodes_list),len(x_full)))

tol=10**-19
tol2=10**-3
k=0
for n_nodes in n_nodes_list: 
    E_min=min(potential_dimless(x_arr))
    E_max=max(potential_dimless(x_arr))
    count=1
    max_iter=100
    e=0
    while count<max_iter:
        e=(E_max+E_min)/2
        n_half=int(n_nodes/2)
        du0=0
        if n_nodes%4==0:
            u0=1
        elif n_nodes%4==1:
            u0=0
            du0=1
        elif n_nodes%4==2:
            u0=-1
        else:
            u0=0
            du0=-1
        cl=nodes_finder(alpha(x_arr,e),key=1)
        x_arr_1=x_arr[:(cl)+1]
        x_arr_2=x_arr[(cl):]
        u_arr_1=numerov_method_HO(u0,x_arr_1,alpha,e,du0)
        u_arr_2=numerov_method_HO(0,np.flip(x_arr_2),alpha,e)
        u_arr=matching_function(u_arr_1,np.flip(u_arr_2))
        n_f=nodes_finder(u_arr)
        
        if n_f<n_half or n_f==n_half:
            E_min=e
        else:
            E_max=e
        if abs(E_max-E_min)<tol2:
            e_cor=fsolve(phi2,e,(u0,du0,cl))    
            e=e_cor
            cl=nodes_finder(alpha(x_arr,e),key=1)
            x_arr_1=x_arr[:(cl)+1]
            x_arr_2=x_arr[(cl):]
            u_arr_1=numerov_method_HO(u0,x_arr_1,alpha,e,du0)
            u_arr_2=numerov_method_HO(0,np.flip(x_arr_2),alpha,e)
            u_arr=matching_function(u_arr_1,np.flip(u_arr_2))
            
            eigen_val_list[k]=e_cor
            u_matrix[k,:]=normalization_factor(parity_operator(u_arr,(k%2)),x_full)[0]
            u_num_mat[k,:]=analytic_HO(x_full,k)
            break
        count=count+1
    k=k+1
    
eigen_val_ana=np.zeros((len(eigen_val_list)))
for i in range(len(eigen_val_list)):
    eigen_val_ana[i]=i+0.5
    
table_mat=np.column_stack((n_nodes_list,eigen_val_list,eigen_val_ana))
df=pd.DataFrame(table_mat,columns=["Energy State","Energy value(numerical)","Energy value(analytic)"])
print(df)


for i in range(k):
    plt.title("Wavefunctions")
    plt.scatter(x_full,u_matrix[i,:],marker=".",label="Energy state "+str(i)+"(numerical)")
    plt.plot(x_full,u_num_mat[i,:],label="Energy state "+str(i)+"(analytic)")
    plt.xlabel("x")
    plt.ylabel("u(x)")
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
#plt.legend(loc='lower right')
plt.grid()
plt.show()

for i in range(k):
    title4="Energy State= "+str(i)+" for x_max=10"
    plt.title(title4)
    plt.scatter(x_full,u_matrix[i,:],marker=".",label="wavefunction")
    plt.scatter(x_full,np.square(u_matrix[i,:]),marker=".",label="probability density")
    i1=int((x_max-np.sqrt((2*i)+1))*(N/x_max))
    i2=int((x_max+np.sqrt((2*i)+1))*(N/x_max))
    y_c=normalization_factor(0.5*np.square(x_full[i1:i2]),x_full[i1:i2])[0]
    plt.plot(x_full[i1:i2],y_c,label="classical")
    plt.vlines(x =[-1*np.sqrt(2*i+1), np.sqrt(2*i+1)], ymin = min(y_c),ymax=max(y_c),
           colors = 'purple',linestyle=':',
           label = 'classical allowed region')
    plt.legend(loc='best')
    plt.grid()
    plt.savefig("C:\\Users\\vyash\\Desktop\\YASH\\quantum mechanics\\practical_qm\\A7\\"+title4)
    plt.show()

for i in range(k):
    plt.scatter(x_full,np.square(u_matrix[i,:]),marker=".",label="Energy state "+str(i)+"(numerical)")
    plt.plot(x_full,np.square(u_num_mat[i,:]),label="Energy state "+str(i)+"(analytic)")
    plt.xlabel("x")
    plt.ylabel("u(x)")
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.title("Probability Density Energy State")
plt.grid()
plt.show()


# In[ ]:





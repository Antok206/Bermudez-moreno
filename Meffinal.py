# -*- coding: utf-8 -*-
"""
Created on Fri Jul  7 11:08:33 2023

@author: antob
"""

from pylab import *
import time
import numpy as np
import scipy
import scipy.linalg
import math


a = 0. # extremo inferior del intervalo
b = 2*pi # extremo superior del intervalo
alpha = 0 # condicion de contorno en a
beta = 0 # condicion de contorno en b
N = 10# numero de particiones menos 1
nodos = np.linspace(a, b, N+1)


h = (b-a)/(N-1) #paso de malla

# Definimos la matriz A como en los apuntes y le aplicamos el metodo LU.
def matrizA(N):
    A = np.zeros([N, N]) # matriz del sistema
    for i in range(N-1):
        A[i, i] = 4*h/3 + 2/h # diagonal
        A[i, i+1] = 2*h/3 - 1/h # superdiagonal
        A[i+1, i] = 2*h/3 - 1/h # subdiagonal
        A[N-1, N-1] = 2. # ultimo elemento de la diagonal
    return A

print(matrizA(N))



#def f(x):
 #   return x

#B = np.zeros(N)
#x = np.linspace(a,b, N)

#for i in range(N):
 #   B[i] = f(x[i])

#print(B)
#Vamos a aplicar el algoritmo de Bérmudez Moreno para comprobar la convergencia en el caso


#Definimos las funciones base de nuestro subespacio V_h
def funcion_a_trozos(x,puntos, j):
    if x <= puntos[j-1]:
        return 0
    elif puntos[j-1] < x < puntos[j]:
        return (x - puntos[j-1]) / h
    elif puntos[j] < x < puntos[j + 1]:
        return (puntos[j+1] - x) / h
    elif x == puntos[j]:
        return 1
    else:
        return 0


#Creamos nuestra función u que transforma el vector U a una función multiplicando 
#cada componente por las funciones de la base
def u(z,lista,nodos):
   u = 0
   for i in range(N):
       u += lista[i]*funcion_a_trozos(z, nodos, i)
   return u

print(u(3,nodos,nodos))
#definimos las derivadas de las funciones base
def derivada_a_trozos(x,puntos,j):
    if x <= puntos[j-1]:
        return 0
    elif puntos[j-1] < x < puntos[j]:
        return 1 / h
    elif puntos[j] < x < puntos[j + 1]:
        return -1 / h
    elif x == puntos[j]:
        return 1
    else:
        return 0
print(derivada_a_trozos(1, nodos, 1))

#definimos la derivada de la función u        
def du(z,lista,nodos):
   u = 0
   for i in range(len(lista)):
       u += lista[i]*derivada_a_trozos(z, nodos, i)
   return u
  
print(du(3,nodos,nodos))       
            
#Construyamos nuestra función p
def funcion(r, y, lb, p):
    # Define la función para la cual quieres encontrar la raíz
    return r**(p-1) + (r-1)/(lb*abs(y)**(p-2))

def derivada(r,y,lb,p):
    # Define la derivada de la función
    return (p-1)*r**(p-2) + 1/(lb*(abs(y)**(p-2)))

def metodo_newton(f, f_derivada, x0, tolerancia, max_iter, p0, lb, y):
    # Implementa el método de Newton
    x = x0
    iteracion = 0

    while abs(f(x,y,lb,p0)) > tolerancia and iteracion < max_iter:
        x = x - f(x,y,lb,p0) / f_derivada(x,y,lb,p0)
        iteracion += 1

    if abs(f(x,y,lb,p0)) <= tolerancia:
        return x
    else:
        return None

#Vamos a fijar parámetros
def p_0(z):
    return 0
tolerancia = 0.0001
max_iter = 100
p0 = 2
x0 = 1

raiz = metodo_newton(funcion, derivada, x0, tolerancia, max_iter,p0,1,2)
if raiz is not None:
    print("La raíz encontrada es:", raiz)
else:
    print("El método de Newton no convergió.") 


def p(z,i,u_inicial,du_inicial,lista,nodos,lb,p_inicial):
    if i == 1:
        return p_inicial(z)
    else:
        return (1- metodo_newton(funcion, derivada, x0, tolerancia, max_iter,p0,lb,u_inicial(z,lista,nodos)))*(du_inicial(z,lista,nodos) + lb*p(z,i-1,u_inicial,du_inicial,lista,nodos,lb,p_inicial))/lb
    
print(p(3,2,u,du,nodos,nodos,1,p_0))

def BM(a,b,nodos,f,p0,u0,lb,eps,Nmax):
     N = len(nodos)
     A = matrizA(N)
     h = (b-a)/(N-1)
     error = eps + 1
     #nodos = np.linspace(a,b, N)
     #i = 1
     B = np.zeros(N)
     for i in range(N):
         B[i] = h*f(nodos[i])
         print(B)
     contador = 1
     
     uprimera = u0((a+b)/2)
     U = np.zeros(N)
     while(error>=eps and i< Nmax ):
             
        P = np.zeros(N)
        P[0] = 0
        P[N-1] = 0
        for i in range(1,N-1):
             P[i] = (p(nodos[i-1],contador,u,du,U,nodos,lb,p0) - p(nodos[i+1],contador,u,du,U,nodos,lb,p0))/2
        B += P
        U = solve(A,B) #Hemos aplicado el método de elementos finitos para desarrollar la matriz A y el vector B
        error = u((a+b)/2,U,nodos) - uprimera
        uprimera = u((a+b)/2,U,nodos)
        contador +=1
     if contador == Nmax:
         return Nmax
     else:
        y = u(nodos,U,nodos)
        plot(nodos,y)
        # Personalizar el gráfico
        plt.xlabel('x')
        plt.ylabel('solución aproximada')
        
            
        

def f_inicial(z):
    return sin(z)


BM(a,b,nodos,f_inicial,p_0,p_0,1,0.0001,200)
         
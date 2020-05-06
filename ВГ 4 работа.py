import math as mt
import scipy as sp
from scipy import integrate
import numpy as np
from sympy import *
C = np.array ([[1,0,0],
               [0,0,0],
               [-0.000484165217061,-0.000000000338846075704,0.00000243934736621]])
S = np.array ([[0,0,0],
               [0,0,0],
               [0,0.00000000146306108906,-0.00000140030429947]])
n =[0,1,2]
k =[0,1,2]
GM=398600441500000
G=6.6743*10**(-11)
H=327376.34*10**(-8)
a=6378136.46
print ("Часть 2")
print ("1) Масса планеты M")
M=GM/G*C[0,0]
print ("M =",M)
print ("2) Координаты центра масс")
x0 = float(mt.sqrt(1*mt.factorial(n[1]+k[1])/(2*(2*n[1]+1)*mt.factorial(n[1]-k[1])))*C[1,1]*a)
y0= float(mt.sqrt(1*mt.factorial(n[1]+k[1])/(2*(2*n[1]+1)*mt.factorial(n[1]-k[1])))*S[1,1]*a)
z0 = float(1/mt.sqrt(2*n[1]+1)*C[1,0]*a)
print ("x0 =" , x0 ,"\n" "y0 =", y0 ,"\n" "z0 =",z0)
print ("3) Разности осевых моментов инерции")
A = float(mt.sqrt(5)*M*a**2*((1-1/H)*C[2,0]-C[2,2]/mt.sqrt(3)))
B = float(mt.sqrt(5)*M*a**2*((1-1/H)*C[2,0]+C[2,2]/mt.sqrt(3)))
Cc = float(-C[2,0]*mt.sqrt(5)*M*a**2/H)
print ("C-A",Cc-A,"\nC-B",Cc-B,"\nB-A",B-A)
print ("4) Максимальный момент инерции Ma^2")
print ("Ma^2 = ", M*a**2)
print ("5) Отношение разностей осевых моментов к максимальному моменту инерции")
print ("(C-A)/Ma^2 =",(Cc-A)/M*a**2,"\n(C-B)/Ma^2 =",(Cc-B)/M*a**2,"\n(B-A)/Ma^2 =",(B-A)/M*a**2)
print ("6) Осевые и центробежные моменты инерции (элементы тензора инерции)")
e= float(5/mt.sqrt(15)*C[2,1]*M*a**2)
E=-e
d= float(5/mt.sqrt(15)*S[2,1]*M*a**2)
D=-d
f= float(5/mt.sqrt(15)*S[2,2]*M*a**2)
F=-f
print ("Осевые моменты инерции\nA =", A,"\nB =", B,"\nC =",Cc)
print ("Центробежные моменты инерции\n-E =", E,"\n-D =", D,"\n-F =",F)
I  = np.array ([[A,F,E],
                [F,B,D],
                [E,D,Cc]])
print("Тензор инерции\n", I)                
print ("7) Безразмерный момент инерции I* для всех осевых моментов инерции")
print ("I*(A) =",A/(M*a**2),"\nI*(B) =",B/(M*a**2),"\nI*(C) =",Cc/(M*a**2))
print ("Часть 3")
I  = np.array ([[A,F,E],
                [F,B,D],
                [E,D,Cc]])
print("Тензор инерции\n", I)
MI, Ncos = np.linalg.eig(I)
print("Главные моменты инерции \n", MI)
print("Собственные векторы тензора инерции(направляющие косинусы)\n", Ncos)
print("Направление главных осей инерции")
la1 = mt.degrees(mt.atan(Ncos[1, 0] / Ncos[0, 0]))
fi1 = mt.degrees(mt.atan(Ncos[2, 0] / mt.sqrt(Ncos[0, 0] ** 2 + Ncos[1, 0] ** 2)))
la2 = mt.degrees(mt.atan(Ncos[1, 1] / Ncos[0, 1]))
fi2 = mt.degrees(mt.atan(Ncos[2, 1] / mt.sqrt(Ncos[0, 1] ** 2 + Ncos[1, 1] ** 2)))
la3 = mt.degrees(mt.atan(Ncos[1, 2] / Ncos[0, 2])) + 180
fi3 = mt.degrees(mt.atan(Ncos[2, 2] / mt.sqrt(Ncos[0, 2] ** 2 + Ncos[1, 2] ** 2)))
print("Широта оси X %1.10f," %fi1,"Долгота оси X %1.10f" %la1,"\nШирота оси Y %1.10f," %fi2,"Долгота оси Y %1.10f" %la2,"\nШирота оси Z %1.10f," %fi3,"Долгота оси Z %1.10f" %la3)
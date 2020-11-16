% Peña Alarcón Raquel
% Práctica 2 
% Algebra lineal numérica e Interpolación


%%Ejercicio 1
A= [0.1 0.2 0.3; 0.4 0.5 0.6 ; 0.7 0.8 0.9]
cond(A)

det(A) %%Determinante de A

b=[0.1 ; 0.3; 0.5]

x=A\b %%Resolviendo el sistema por Eliminación Gaussiana

x=inv(A)*b  %%Resolviendo el sistema por otro método

cond(A,1) %%Número de condición norma 1
cond(A,2) %%Norma 2
cond(A,inf) %%Norma infinita

%%Ejercicio 2
function [x,U] = gaussianasimple (A,b)

 
 


% Pe�a Alarc�n Raquel
% Pr�ctica 2 
% Algebra lineal num�rica e Interpolaci�n


%%Ejercicio 1
A= [0.1 0.2 0.3; 0.4 0.5 0.6 ; 0.7 0.8 0.9]
cond(A)

det(A) %%Determinante de A

b=[0.1 ; 0.3; 0.5]

x=A\b %%Resolviendo el sistema por Eliminaci�n Gaussiana

x=inv(A)*b  %%Resolviendo el sistema por otro m�todo

cond(A,1) %%N�mero de condici�n norma 1
cond(A,2) %%Norma 2
cond(A,inf) %%Norma infinita

%%Ejercicio 2
function [x,U] = gaussianasimple (A,b)

 
 


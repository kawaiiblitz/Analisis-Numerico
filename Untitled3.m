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
%%a)

A=[21 67 88 73; 76 63 7 20; 0 85 56 54; 19.3 43 30.2 29.4]
b=[141; 109; 218; 93.7]

function [x,U] = gaussianasimple (A,b)
[n,m]=size(A); %%Registra en los valores n y m el tama�o de la matriz
U=single(A);
b=single(b);
U(1:n,m+1)=single(b);  %%La matriz aumentada
for j=1:m %%Eliminaci�n gaussiana
    for i=j+1: m
        U(i,j:m+1)=single(U(i,j:m+1)+single(U(j,j:m+1)*single(single(-U(i,j))/U(j,j))));
    end
end

b=single(U(1:n, m+1)); %%Asigna a b la �ltima columna de la matriz aumentada
U(:, m+1)=[];          %%Elimina la �ltima columna
U=single(U);
x=single(zeros(n,1));
for i=n: -1:1           %%Sustituci�n hacia atr�s, inicia en n, termina en 1, decrementos de -1
    x(i,1)=single(single(b(i)-U(i,i:n)*x(i:n,1))/U(i,i));
end
end

gasussianasimple(A,b)


%%b)
A=[21 67 88 73; 76 63 7 20; 0 85 56 54; 19.3 43 30.2 29.4]
b=[141; 109; 218; 93.7]

x=gaussianasimple(A,b)   %%Ocupando precisi�n simple

r=double(double(b)-double(double(A)*double(x)));
r=single(r)


%%c)
z=gaussianasimple(A,r)

x=x+z %%Soluci�n mejorada


%%d)



%%Ejercicio 3
%%a)
for k=1:10;
    e=(10^(-2*k));
    A=[e 1; 1 1 ];
    b=[1+e; 2];
    
    B=[A b];    %%Eliminaci�n Gaussiana, matriz aumentada
    B(2,:)=(-B(2,1)/B(1,1))*B(1,:)+B(2,:);  %% El pivote inicia en (1,1) en la matriz A, hacemos el primer 0, s�lo se hace un pivoteo.
    B(1,:)=(-B(1,2)/B(2,2))*B(2,:)+B(1,:);  %%El cero est� en la entrada (1,2)
    B(1,:)=(1/B(1,1))*B(1,:);
    B(2,:)=(1/B(2,2))*B(2,:);
    x=B(: , 3)
end

%%b)
for k=1:10;
    e=(10^(-2*k));
    A=[e 1; 1 1 ];
    b=[1+e; 2];
    
    B=[A b];    %%Eliminaci�n Gaussiana, matriz aumentada
    B(2,:)=(-B(2,1)/B(1,1))*B(1,:)+B(2,:);  %% El pivote inicia en (1,1) en la matriz A, hacemos el primer 0, s�lo se hace un pivoteo.
    B(1,:)=(-B(1,2)/B(2,2))*B(2,:)+B(1,:);  %%El cero est� en la entrada (1,2)
    B(1,:)=(1/B(1,1))*B(1,:);
    B(2,:)=(1/B(2,2))*B(2,:);
    x=B(:,3)


%%Para mejorar la soluci�n utilizamos residuales
    r=b-A*x;
    z=[A x]; %%La matriz aumentada
    z(2,:)=(-z(2,1)/z(1,1))*z(1,:)+z(2,:);
    z(1,:)=(-z(1,2)/z(2,2))*z(2,:)+z(1,:);
    z(1,:)=(1/z(1,1))*z(1,:);
    z(2,:)=(1/z(2,2))*z(2,:);
    x=z(:,3);
    X=z+x
end


%%Ejercicio 4
for i=21:25
    format long
    
    e=sqrt(2^(-i));
    e=single(e)
    
    A=[1 1+e; 1-e 1];
    b=[1+(e*(1+e));1];
    A=single(A)
    b=single(b)
    
    NumeroCondicion=cond(A);
    NumeroCondicion=single(NumeroCondicion)
    
    x=gaussianasimple(A,b)
    
end

    
    




    
















 
 


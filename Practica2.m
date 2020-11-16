% Peña Alarcón Raquel
% Práctica 2 
% Algebra lineal numérica e Interpolación


%%Ejercicio 1
A= [0.1 0.2 0.3; 0.4 0.5 0.6 ; 0.7 0.8 0.9]
cond(A)

det(A) %%Determinante de A
b=[0.1 ; 0.3; 0.5]

x=A\b %%Resolviendo el sistema método de Backlash
x=inv(A)*b  %%Resolviendo el sistema por otro método

%%Resolviendo el sistema por Eliminación Gaussiana
function [x,U] = gaussianasimple (A,b)
[n,m]=size(A); %%Registra en los valores n y m el tamaño de la matriz
U=single(A);
b=single(b);
U(1:n,m+1)=single(b);  %%La matriz aumentada
for j=1:m %%Eliminación gaussiana
    for i=j+1: m
        U(i,j:m+1)=single(U(i,j:m+1)+single(U(j,j:m+1)*single(single(-U(i,j))/U(j,j))));
    end
end

b=single(U(1:n, m+1)); %%Asigna a b la última columna de la matriz aumentada
U(:, m+1)=[];          %%Elimina la última columna
U=single(U);
x=single(zeros(n,1));   %%Genera un vector de ceros
for i=n: -1:1           %%Sustitución hacia atrás, inicia en n, termina en 1, decrementos de -1
    x(i,1)=single(single(b(i)-U(i,i:n)*x(i:n,1))/U(i,i));
end
end

x=gaussianasimple(A,b)


cond(A,1) %%Número de condición norma 1
cond(A,2) %%Norma 2
cond(A,inf) %%Norma infinita


%%Ejercicio 2
%%a)

A=[21 67 88 73; 76 63 7 20; 0 85 56 54; 19.3 43 30.2 29.4]
b=[141; 109; 218; 93.7]

function [x,U] = gaussianasimple (A,b)
[n,m]=size(A); %%Registra en los valores n y m el tamaño de la matriz
U=single(A);
b=single(b);
U(1:n,m+1)=single(b);  %%La matriz aumentada
for j=1:m %%Eliminación gaussiana
    for i=j+1: m
        U(i,j:m+1)=single(U(i,j:m+1)+single(U(j,j:m+1)*single(single(-U(i,j))/U(j,j))));
    end
end

b=single(U(1:n, m+1)); %%Asigna a b la última columna de la matriz aumentada
U(:, m+1)=[];          %%Elimina la última columna
U=single(U);
x=single(zeros(n,1));   %%Genera un vector de ceros
for i=n: -1:1           %%Sustitución hacia atrás, inicia en n, termina en 1, decrementos de -1
    x(i,1)=single(single(b(i)-U(i,i:n)*x(i:n,1))/U(i,i));
end
end

x=gaussianasimple(A,b)


%%b)
A=[21 67 88 73; 76 63 7 20; 0 85 56 54; 19.3 43 30.2 29.4];
b=[141; 109; 218; 93.7];

x=gaussianasimple(A,b)   %%Ocupando precisión simple

r=double(double(b)-double(double(A)*double(x))); %%Cálculo del residual

r=single(r)


%%c)
z=gaussianasimple(A,r)

x=x+z %%Solución mejorada


%%d)
for i=1:1:10;
    r=double(double(b)-double(double(A)*double(x))); %%Cálculo del residual
    r=single(r)
    z=gaussianasimple(A,r)
    x=x+z
end




%%Ejercicio 3
%%a)
for k=1:10;
    e=(10^(-2*k));
    A=[e 1; 1 1 ];
    b=[1+e; 2];
    
    B=[A b];    %%Eliminación Gaussiana, matriz aumentada
    B(2,:)=(-B(2,1)/B(1,1))*B(1,:)+B(2,:);  %% El pivote inicia en (1,1) en la matriz A, hacemos el primer 0, sólo se hace un pivoteo.
    B(1,:)=(-B(1,2)/B(2,2))*B(2,:)+B(1,:);  %%El cero está en la entrada (1,2)
    B(1,:)=(1/B(1,1))*B(1,:);
    B(2,:)=(1/B(2,2))*B(2,:);
    x=B(: , 3)
end

%%b) 

for k=1:10;
    e=(10^(-2*k));
    A=[e 1; 1 1 ];
    b=[1+e; 2];
    
    B=[A b];    %%Eliminación Gaussiana, matriz aumentada
    B(2,:)=(-B(2,1)/B(1,1))*B(1,:)+B(2,:);  %% El pivote inicia en (1,1) en la matriz A, hacemos el primer 0, sólo se hace un pivoteo.
    B(1,:)=(-B(1,2)/B(2,2))*B(2,:)+B(1,:);  %%El cero está en la entrada (1,2)
    B(1,:)=(1/B(1,1))*B(1,:);
    B(2,:)=(1/B(2,2))*B(2,:);
    x=B(: , 3)



%%Para mejorar la solución utilizamos residuales
    r=b-A*x;
    z=[A r]; %%La matriz aumentada
    z(2,:)=(-z(2,1)/z(1,1))*z(1,:)+z(2,:);
    z(1,:)=(-z(1,2)/z(2,2))*z(2,:)+z(1,:);
    z(1,:)=1/(z(1,1))*z(1,:);
    z(2,:)=1/(z(2,2))*z(2,:);
    d=z(:,3);
    X=d+x
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
   
   B=[A b];    %%Eliminación Gaussiana, matriz aumentada
    B(2,:)=(-B(2,1)/B(1,1))*B(1,:)+B(2,:);  %% El pivote inicia en (1,1) en la matriz A, hacemos el primer 0, sólo se hace un pivoteo.
    B(1,:)=(-B(1,2)/B(2,2))*B(2,:)+B(1,:);  %%El cero está en la entrada (1,2)
    B(1,:)=(1/B(1,1))*B(1,:);
    B(2,:)=(1/B(2,2))*B(2,:);
    x=B(: , 3)
end


   
%%Ejercicio 5
disp('La dimensión de la matriz a calcular');
disp('                                            ');
n=single(input('Ingresar el valor de la dimensión n='));
disp('                                            ');
disp('Ingresar los valores de la matriz a calcular entrada por entrada');
disp('La matriz debe ser tri-diagonal');
disp('                                            ');

for i=1:n     %%Ciclo que llenará la matriz A
    for j=1:n
        disp('La entrada:');
        fprintf('(%d,%d)',i,j);
        disp(' ')
        A(i,j)=single(input(' '));
        j=j+1;
    end
    i=i+1;
end

disp('La matriz A es:');
disp(A);

disp('Ingresar los valores del vector de resultados b');
for i=1:n %%Ciclo que llenará el vector b
    disp('La entrada:');
    fprintf('(%d,%d)',i,1);
    disp(' ');
    b(i,1)=single(input(' '));
    i=i+1;
end

disp('El vector b es:');
disp(b);
A=single(A);
B=single([A b]); %%Matriz aumentada

%%Pivoteo parcial
for i=1:n-1
    if abs(B(1+1,i))>abs(B(i,i))    %%Permutación si ocurre que la primer entrada del renglón n+1 es mayor que la primer
        c=B(i,:);                  %%entrada del renglón n.
        B(i,:)=B(i+1,:);
        B(1+1,:)=c;
    end
    B
    m=(-B(i+1,i)/B(i,i))*B(i,:);
    B(i+1,:)=B(i+1,:)+m;            %%Nos da el cero debajo del pivote
end

disp('        ');

%%Obtener la solución
for j=1:n-1
    for k=j:n-1
        if (B(n-j+1,n-j+1)~=0 && B(n-k,n-j+1)~=0)
            m=(-B(n-k,n-j+1)/B(n-j+1,n-j+1))*B(n-j+1,:);
            B(n-k,:)=B(n-k,:)+m;
        end
    end
end

disp('Obtuvimos la matriz B:');
disp(B);
disp('La solución al sistema tri-diagonal es');

for i=1:n
    x=B(i,i);   %%Asigna a c la diagonal de la matriz B
    y=B(i,n+1); %%Asigna a y el vector de resultados b
    z=y/x;
    disp(z);
    i=i+1;
end


    
        





%%Ejercicio 6
%%a
t=[0 1 4 9 16 25 36 49 64]
y=[0 1 2 3 4 5 6 7 8 ]

x=0:.1:64
polagrange= polyfit(t,y,8)          %%Polinomio de Lagrange vía Matlab
lagrange= polyval(polagrange,x)    %%Polinomio evaluado en [0,64]
y1=sqrt(t)

plot(t,y,'*',x,lagrange,'b-',t,y1,'m-')
axis([0 65 -5 50])
legend('Valores exactos de y', 'Pol.de Lagrange grado 8', 'y=sqrt(t)')


%%Haciendo más pequeña la escala de t:
t=[0 1 2.25 4 6.25 9 12.25 16 20.25 25 30.25 36 42.25 49 56.25 64]
y=[0 1 1.5 2 2.5 3 3.25 4 4.5 5 5.5 6 6.5 7 7.5 8 ]

x=0:.1:64
polagrange= polyfit(t,y,8)          %%Polinomio de Lagrange vía Matlab
lagrange= polyval(polagrange,x)    %%Polinomio evaluado en [0,64]
y1=sqrt(t)

plot(t,y,'*',x,lagrange,'b-',t,y1,'m-')
axis([0 65 -5 10])
legend('Valores exactos de y', 'Pol.de Lagrange grado 8', 'y=sqrt(t)')

%%c)
t=[0 1 4 9 16 25 36 49 64]
y=[0 1 2 3 4 5 6 7 8 ]
x=0:.1:64

y1=sqrt(x)
y2=spline(t,y,x)
polagrange= polyfit(t,y,8)          
y3= polyval(polagrange,x)
norm(y3-y1)  %%Norma entre lagrange y la función
norm(y2-y1)  %%Norma entre spline cúbico y la función



%%d)
t=[0 1 4 9 16 25 36 49 64]
y=[0 1 2 3 4 5 6 7 8 ]

x=0:.1:1
polagrange= polyfit(t,y,8)          %%Polinomio de Lagrange vía Matlab
lagrange= polyval(polagrange,x)    %%Polinomio evaluado en [0,64]
y1=sqrt(t)

plot(t,y,'*',x,lagrange,'b-',t,y1,'m-')
axis([0 1 -0.2 1.2])
legend('Valores exactos de y', 'Pol.de Lagrange grado 8', 'y=sqrt(t)')
norm(lagrange,y1)


%%Ejercicio 7
%%a)
t=[ 1 2 3 4 5]
y=[ 1 1 2 6 24]


x=0:.1:5
polagrange= polyfit(t,y,4)          %%Polinomio de Lagrange vía Matlab
lagrange= polyval(polagrange,x)    %%Polinomio evaluado en [1,5]
y1= factorial(t-1)

plot(t,y,'*',x,lagrange,'b-',t,y1,'m-')
axis([-1 6 -1 26])
legend('Valores exactos de y', 'Pol.de Lagrange grado 4', 'y=(t-1)!')

%%c)
t=[ 1 2 3 4 5]
y=[ 1 1 2 6 24]
x=1:.1:5

y1= gamma(x) %%Aquí utilizamos la función gamma, pues hubo problemas con la función definida anteriormente para números negativos
y2=spline(t,y,x)
polagrange= polyfit(t,y,4)          
y3= polyval(polagrange,x)
norm(y3-y1)  %%Norma entre lagrange y la función
norm(y2-y1)  %%Norma entre spline cúbico y la función
    

%%d)
t=[ 1 2 3 4 5]
y=[ 1 1 2 6 24]


x=1:.1:2.5
polagrange= polyfit(t,y,4)          %%Polinomio de Lagrange vía Matlab
lagrange= polyval(polagrange,x)    %%Polinomio evaluado en [1,2]
y1= factorial(t-1)

plot(t,y,'*',x,lagrange,'b-',t,y1,'m-')
axis([.5 2.5 .5 2.5])
legend('Valores exactos de y', 'Pol.de Lagrange grado 4', 'y=(t-1)!')


%%Ejercicio 8 
%%a)
x=[1900: 10: 1980]
y=[76212168 92228496 106021537 123202624 132164569 151325798 179323175 203302031 226542199]

for i=1:9
    y1(i)=x(i);
    i=i+1;
end
Vander1=vander(y1) %%Matriz de Vandermonde

for i=1:9
    y2(i)=(x(i)-1900);
    i=i+1;
end
Vander2=vander(y2)
    
for i=1:9
    y3(i)=(x(i)-1940);
    i=i+1;
end
Vander3=vander(y3)


for i=1:9
    y4(i)=((x(i)-1940)/40);
    i=i+1;
end
Vander4=vander(y4)

%%Número de condición para cada matriz de Vandermonde
cond(Vander1)
cond(Vander2)
cond(Vander3)
cond(Vander4)

%%b)
u=Vander4\y' %%Pues es la del número de condición mas pequeño
u1=u'
xx=sym('xx');
yy=sym('yy');
yy=u1(1)*xx.^8+u1(2)*xx.^7+u1(3)*xx.^6+u1(4)*xx.^5+u1(5)*xx.^4+u1(6)*xx.^3+u1(7)*xx.^2+u1(8)*xx.^1+u1(9)*xx.^0
horner=horner(yy)

z=1900:1:1980;
z=(z-1940)/40;
p1=subs(horner,z);

plot(y4,y,' -',z,p1,'b*');
title('Interpolación con polinomio de grado 8');
xlabel('Años');
ylabel('Valores del polinomio');
legend('Valores reales', 'Polinomio de interplación de grado 8');


%%d)
residuo=(1990-1940)/40
res1=subs(horner, residuo)

format long
spline(y4,y,residuo)

%%e)
x=[1900      1910      1920 1930 1940 1950 1960 1970 1980]
y=[76212168 92228496 106021537 123202624 132164569 151325798 179323175 203302031 226542199]


lag=polyfit(x,y,8);
 
figure(1)
hold on
plot(x,lag,'m')
hold off







    
    



    
















 
 


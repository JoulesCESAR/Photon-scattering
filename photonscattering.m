%definimos los parametros
us=1;
ua=0;
d=2;
g=0.6;

%definimos el numero de fotones
Nfon=1000;

%contador de los q salieron
k=1;

%inicializamos
x=0;
y=0;
z=0;
xf=0;
yf=0;
zf=0;

%inicializamos
med=0;
lrfsuma=0;
%numero total de dispersiones
disptotal=0;

%iteracion de fotones
for j=1:Nfon

    %inicializamos las posiciones finales;
    x(1)=0;
    y(1)=0;
    z(1)=0;

    ux(1)=0;
    uy(1)=0;
    uz(1)=1;

    %dstotal suma los ds para cada foton
    dstotal=0;
    disp=1;

    %correccion del ds
    dsc=0;

    %CALCULAMOS PARA LA 1º DISPERSION

    %calculamos el ds
    epsilon=rand ();
    ds=-log(epsilon)/us;

    %sumamos los ds
    dstotal=dstotal+ds;

    % calculemos el teta
    epsilon=rand ();
    costeta=(1/(2*g))*(1+g*g-((1-g*g)/(1-g+2*g*epsilon))^2);
    teta=acos(costeta);

    %calculemos el phi
    epsilon=rand ();
    phi=2*pi*epsilon;

    %los cálculos para la 1º dispersion
    ux(2)=sin(teta)*cos(phi);
    uy(2)=sin(teta)*sin(phi);
    uz(2)=cos(phi);
    x(2)=x(1)+ds*ux(2);
    y(2)=y(1)+ds*uy(2);
    z(2)=z(1)+ds*uz(2);

    %desde la 2º posicion, 1º dispersion
    i=2;

    %iteracion de dispersiones
    while (0<=z(i) & z(i)<=d)

        %calculemos el ds
        epsilon=rand ();
        ds=-log(epsilon)/us;

        %sumamos los ds
        dstotal=dstotal+ds;

        % calculemos el teta
        epsilon=rand ();
        costeta=(1/(2*g))*(1+g*g-((1-g*g)/(1-g+2*g*epsilon))^2);
        teta=acos(costeta);

        %calculemos el phi
        epsilon=rand ();
        phi=2*pi*epsilon;

        %los cálculos para la nueva posiciones
        %UNITARIOS
        ux(i+1)=sin(teta)*(ux(i)*uz(i)*cos(phi)-uy(i)*sin(phi))*(1-uz(i)*uz(i))^(-1/2)+ux(i)*cos(teta);
        uy(i+1)=sin(teta)*(uz(i)*uy(i)*cos(phi)+ux(i)*sin(phi))*(1-uz(i)*uz(i))^(-1/2)+uy(i)*cos(teta);
        uz(i+1)=-sin(teta)*cos(phi)*(1-uz(i)*uz(i))^(1/2)+uz(i)*cos(teta);
        %vector
        x(i+1)=x(i)+ds*ux(i+1);
        y(i+1)=y(i)+ds*uy(i+1);
        z(i+1)=z(i)+ds*uz(i+1);

        %avanzamos al siguiente paso
        i=i+1;
        disp=disp+1;
    end

    %hacemos la correcion de la posicion final
    if (z(i)>d)
        %inicializamos
        zme=0;
        zma=0;
        Rma=0;
        Rme=0;
        cx=0;
        cy=0;

        %definimos los elementos del triangulo
        zme=z(i)-d;
        zma=z(i)-z(i-1);

        %el Rma y Rme son los catetos del triangulo principal
        Rma=(x(i)*x(i)+y(i)*y(i))^(1/2);
        Rme=(Rma*zme)/zma;

        %del triangulo secundario sacamos cx, cy, dsc
        cy=(y(i)*Rme)/Rma;
        cx=(x(i)*Rme)/Rma;
        dsc=(zme*ds)/zma;

        %posiciones finales corregidas
        x(i)=x(i)-cx;
        y(i)=y(i)-cy;

        %grabamos en un vector las posiciones finales
        k=k+1;
        xf(k)=x(i);
        yf(k)=y(i);
        zf(k)=z(i);
    end

    %calculo el clm de cada foton y hago la suma acumulada de med, va sin
    %correccion
    med=med+dstotal/disp;
    
    %numero de dispersiones
    disptotal=disptotal+disp;
    
    %aplico la correcion al dstotal
    dstotal=dstotal-dsc;

    %calculo de lrf de cada foton y hago la suma acumulada de lrf, va con
    %correccion
    lrfsuma=lrfsuma+dstotal;
end

%numero medio de dispersiones
ndismedio=disptotal/Nfon

%camino libre medio
clm=med/Nfon

%longitud recorrida media
lrm=lrfsuma/Nfon

%graficamos los puntos de llegada
plot(xf,yf,'.')

%ahora realizamos el histograma
%para esto usamos a los vectores xf, yf que tienen las posiciones finales

%PROBAMOS LO SIGUIENTE
%vamos a dividir el plano en 100*100 celdas
celdas=100;
%el plano tiene ancho y largo 10*d

%tenemos los extremos:
extizq=-2*d;
extsup=2*d;
enx=extizq;
eny=extsup;
dx=d/celdas;
dy=d/celdas;

%izq a der
%for r=1:celdas
%estoy en la coordenada x: enx
%   enx=enx+dx;

%sup a inf
%  for s=1:celdas
%estoy en la coordenada y: eny
%     eny=eny+dy;


matriz=[xf;yf];
options = statset('Display','final');
%obj = gmdistribution.fit(matriz,2,'Options',options);

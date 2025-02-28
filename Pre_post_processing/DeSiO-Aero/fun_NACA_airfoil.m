function [airfoil] = fun_NACA_airfoil(m,p,c,t,x)

% function of uncambered line of symmetrical airfoil NACA-00xx:
yt = 5*t*( 0.2969*sqrt(x/c) - 0.1260*x/c - 0.3516*(x/c)^2 + 0.2843*(x/c)^3 - 0.1015*(x/c)^4 );

% function of cmber line of cambered airfoil NACA-mpxx:
if x>=0 && x<=p*c
    yc = m/p^2 * ( 2*p*(x/c) - (x/c)^2 );
elseif x>=p*c && x<=1*c
    yc = m/(1-p)^2 * ( (1-2*p) + 2*p*(x/c) - (x/c)^2 );
end

% derivertive of camber line, to get slope rel. to uncambered line
if  x>=0 && x<=p*c
    tan_theta = 2*m/p^2 * ( p - x/c );
elseif x>=p*c && x<=1*c
    tan_theta = 2*m/(1-p)^2 * ( p - x/c );
end
theta = atan(tan_theta);

if p == 0; yc    = 0; end
if p == 0; theta = 0; end

% function of upper and lower airfoil profil lines
xu =  x - yt*sin(theta); xl =  x + yt*sin(theta);
yu = yc + yt*cos(theta); yl = yc - yt*cos(theta);

airfoil.yt = yt*c; 
airfoil.yc = yc*c; 
airfoil.xu = xu; 
airfoil.xl = xl; 
airfoil.yu = yu*c;
airfoil.yl = yl*c;

return
function g = SpherePoints(NPore,sigma,doPlot, filename)
%
% This routine distributes NPore pores of radius a in a Fibonacci Spiral,
% fixes Nfix points on the boundary of the pore and one point in the center of the pore
% and then dynamically distributes Nmove additional points. We put a point
% in the center of the pore to make sure the surface is strictly convex.
%
% We use fmincon to distribute the points by minimizing the pairwise energy
% with a set of fixed points [Z] and a set of moving points [X]
% 
%  E = sum(1/(h+dij^2) -1/h) 
%
% where h is a small mollifying parameter and 
% dij is the distance between the i^th and j^th point.
%
% Finally, we triangulate by calling createPolyGeom
%
%   Input parameters:
%       NPore = <number of pores>
%       a= <pore radius>
%   Other parameters:
%
%   Flags
%       doPlot= < Plot flag> use 1 for plotting

    doPlot = 0; % Default is no plots

% Use fmincon to minimize a pairwise interaction function on a sphere with a
% set of fixed points [Z] and a set of moving points [X]

% **** PORE GEOMETRY ****
% We assume the radius of the sphere is 1.
% We measure the pore radius, a,  as the length of a chord connecting 
% the center of the pore (on the surface of the sphere) to a point 
% on the boundary of the pore. This yields the following:
%
%    a= pore radius (measured as a chord)
%    atheta = angle subtended by pore = 2*asin(a/2)
%     *** note that a =2*sin(atheta/2) ***
%    znorth = height of a pore boundary centered on the north pole
%           = cos(atheta) = 1-a^2/2
%    rnorth = radius of a pore boundary centered on the north pole
%           = sin(atheta) = a*sqrt(1-a^2/4)
%    pore area = 2*pi*(1-cos(atheta) = pi*a^2
%    sigma = area fraction per pore = pi*a^2/(4*pi)=a^2/4


if (nargin == 0) % If no input use these default parameters
    NPore = 5;   % Number of pores (must be odd)
    sigma = 0.02; % Specifying the area fraction  
    doPlot = 1;  % Plotting is on
    Nfix = 10; % Number of points on the edge of a pore
    Nmove = 2500;  % Points for triangulation on the surface of the sphere
    filename = "gfile_" + Nfix + "fix_" + Nmove + "move";
end

MPore = floor(NPore/2);   % This line and next assures the number of pore is odd
NPore = 2*MPore+1; 

a = sqrt(4*sigma/NPore); % Pore radius

%
% Compute griding automatically
%

delta=2*a; % typical gridsize

Nfix = 10; % Number of points on the edge of a pore               
Nmove = round(8*pi/sqrt(3)/delta^2);  % Points for triangulation on the surface of the sphere
              % = (number of traingles)/2 = ((4*pi)/(sqrt(3)*delta^2/4))/2
              
%
% Fibonacci Lattice Set up 
%
thPore =  pi/2 - asin(2*(-MPore:1:MPore)/NPore); % Declination
phiPore =  4*pi*(-MPore:1:MPore)/(1+sqrt(5));   % Polar angle

polar=linspace(0,2*pi,Nfix+1); % points distributed on edge of a pore
polar=polar(1:Nfix); % Remove redundant last point

Xc(1,:) = (a.*sqrt(1-a.^2/4)).* cos(polar); % x-coordinates of a pore at the north pole
Xc(2,:) = (a.*sqrt(1-a.^2/4)).* sin(polar); % y-coordinates of a pore at the north pole
Xc(3,:) = (1-a.^2/2)*ones(size(polar));     % z-coordinates of a pore at the north pole
% Add Xc center here!!!
newPoint = [0,0,1]';
Xc = [Xc newPoint];

X = [];  %Pre-allocation for the fixed points on the pore boundary

CoordsPore = zeros(3,NPore);

for j = 1:length(thPore)
        B = radrotz(phiPore(j))*radroty(thPore(j))*Xc; %rotate pores into place
        X = [X B]; %Append pore boundary points to X list
        CoordsPore(:,j)= B(:,end);
end

% Distribute moving points initially as a Fibonacci spiral also

M = floor(Nmove/2); % Assure the number of moving points is odd
N = 2*M+1;

th =  pi/2 - asin(2*(-M:1:M)/N);  % declination angle
phi =  4*pi*(-M:1:M)/(1+sqrt(5)); % polar angle

x(3,:) = sin(th).*cos(phi);
x(2,:) = sin(th).*sin(phi);
x(1,:) = cos(th);

% This is our best guess at optimal parameters to run fmincon
options = optimoptions('fmincon','Algorithm','interior-point','MaxFunctionEvaluations',1e5,...
    'TolFun',1e-4,'TolX',1e-4,'SpecifyObjectiveGradient',true); % run interior-point algorithm

%xs is a vector of moving points of the form [x locations, y locations,z locations]
tic
xs = fmincon(@(x) pairwise(x,X),x(:),[],[],[],[],[],[],@(x) mycon(x),options); %minimize
toc
xs = reshape(xs,3,numel(xs)/3); %reshape as 3-vectors

%append the fixed points on edge of the pore
x = [xs(1,:), X(1,:)]; 
y = [xs(2,:), X(2,:)];
z = [xs(3,:), X(3,:)];

disp('SpherePoints: Minimization Complete')

% Use the PolyGeom routine to create the polyhedral data structure

g = createPolyGeom(x(:),y(:),z(:));
g.PoreCenters = CoordsPore;

disp('SpherePoints: createPolyGeom Complete')

% cart2sph puts the grid points into spherical coordinates

[phiIn,ele,~] = cart2sph(g.P(:,1),g.P(:,2),g.P(:,3));
thetaIn = pi/2 -ele;  % Fixes an unfortunate definition of the elevation angle
indABS = [];

%
%  Test to see if the incenter of a triangle is inside a pore
%  if it is, mark it as absorbing.
% 
g.BC = 0*g.BC;
g.a = a;

for j = 1:length(thPore)  
    d2 = 4*( sin(thetaIn)*sin(thPore(j)) .* sin(0.5*(phiIn-phiPore(j))).^2  + sin(0.5*(thetaIn-thPore(j))).^2 );
    indCurrentPore = find(d2<a^2);
    g.BC(indCurrentPore) = j;
    indABS = [indABS; indCurrentPore];
end

save(filename, 'g');
if (doPlot)
    %  Only plot if doPlot=1
    close all;
    color1 = [0.4 0.6 0.7];
    
    indREFLECT = setdiff(1:size(g.tri,1),indABS);
    
    figure('color','w')
    hold on;
    trisurf(g.tri(indREFLECT,:),g.X(:,1),g.X(:,2),g.X(:,3),'FaceColor', 'w', 'faceAlpha', 0.8);
    trisurf(g.tri(indABS,:),g.X(:,1),g.X(:,2),g.X(:,3),'FaceColor',color1,'faceAlpha', 0.8);
    %plot3(X(1,:),X(2,:),X(3,:),'or');
    xlabel('x','interpreter','latex');
    ylabel('y','interpreter','latex');
    axis off;
    hold off;
end

% Pairwise Energy Functional
function [f,g] = pairwise(x,z)

% Pairwise energy functional
n=numel(x)/3;

% Relative move/fixed weights
a=1/n^2; %move-move energy
b=1/n^2; %move-fixed energy

% Moving points as 3 vector
x = reshape(x,3,n);

% Move-move distances^2
d2xx = (x(1,:) - x(1,:)').^2 + (x(2,:) - x(2,:)').^2 + (x(3,:) - x(3,:)').^2;
% Move-fixed distances^2
d2xz = (z(1,:) - x(1,:)').^2 + (z(2,:) - x(2,:)').^2 + (z(3,:) - x(3,:)').^2;

% inverse square, but remove diagonal
I=eye(n,n);
e2xx = 1./(d2xx+I)-I;

e2xz = 1./d2xz;

% factor of 1/2 removes double counting 
f = a*1/2*sum(sum(e2xx))+b*sum(sum(e2xz));

% Gradient
g = zeros(size(x));

 g(1,:) = a*2*sum((x(1,:) - x(1,:)') .* e2xx.^2,2)+b*2*sum((z(1,:) - x(1,:)') .* e2xz.^2,2);
 g(2,:) = a*2*sum((x(2,:) - x(2,:)') .* e2xx.^2,2)+b*2*sum((z(2,:) - x(2,:)') .* e2xz.^2,2);
 g(3,:) = a*2*sum((x(3,:) - x(3,:)') .* e2xx.^2,2)+b*2*sum((z(3,:) - x(3,:)') .* e2xz.^2,2);


function [c,ceq] = mycon(x)
% Constrained to be on the surface of a unit sphere
X = reshape(x,3,numel(x)/3);
c = [];
ceq = sum(X.^2)-1;
ceq = ceq(:);

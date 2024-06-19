function g = EllipsoidPointsFixed(NPoly, NRing, Rpore, Zpore, Requator, delta, doPlot,filename)

addpath('functions');

%%% This routine distributes (2*NRing +1)*NPoly+2 points in a deterministic fashion to create a surface.
% 
% We triangulate by calling createPolyGeom
%
%   Input parameters:
%       NPoly = <number of points in each hoop>
%       NRing = There are 2*NRing +1 slices in parametrization
%       Rpore = Pore radius
%       Zpore = Pore height
%       Requator = Equator radius
%       delta =Pore center/equator desingularization parameter
%
%   Flags
%       doPlot= < Plot flag> use 1 for plotting
%

if (nargin == 0) % If no input use these default parameters
    NPoly = 12;
    NRing = 5;
    Rpore = 1;
    Zpore = 1;
    Requator = 5;
    delta = 0.01;
    %
    doPlot = 1; 
    filename = "EllipseTestDefault";
end

    ZporeReg=Zpore*(1+delta);
    RequatorReg=Requator*(1+delta);
    b=Zpore/sqrt(1-(Rpore/RequatorReg)^2);

% Angular distribution of points
    polar=linspace(0,2*pi,NPoly+1); % points distributed on edge of a pore
    polar=polar(1:NPoly);           % Remove redundant last point
    shift=(polar(2)-polar(1))/2;    % Need to shift every other ring

% Determine ring heights
   NZ=[-NRing:NRing];
   Zslice=Zpore*NZ/NRing;
   Rslice=RequatorReg*sqrt(1-(Zslice/b).^2);

% Now build surface points   
   [Rgrid,Thetagrid]=ndgrid(Rslice,polar); 
   [Zgrid,~]=ndgrid(Zslice,polar); 
   [NZgrid,~]=ndgrid(NZ,polar);
   Shiftgrid=shift*NZgrid;
   Npts=numel(Rgrid);

% Back to Cartesians
    Xgrid=Rgrid.*cos(Thetagrid+Shiftgrid);
    Ygrid=Rgrid.*sin(Thetagrid+Shiftgrid);
     
% Unpack points
    x = [0 reshape(Xgrid,1,Npts) 0]; 
    y = [0 reshape(Ygrid,1,Npts) 0];
    z = [-ZporeReg reshape(Zgrid,1,Npts) ZporeReg];
 
%%% Use the PolyGeom routine to create the polyhedral data structure
   g = createPolyGeom(x(:),y(:),z(:));


%%% Mark absorbing facets 
%  Test to see if the incenter of a triangle is inside a pore if it is, mark it as absorbing.
    indABS = find( abs(g.P(:,3)) > Zpore );
    g.BC = 0*g.BC;
    g.BC(indABS) = 1;
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
    axis equal;
    view([98 31])

end


return

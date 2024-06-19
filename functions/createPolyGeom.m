function g = createPolyGeom(x,y,z)
% This routine creates a data structure for the polyhedral surface. The vector of points (x,y,z) are the vertices of the
% polyhedron
doPlot = 0;

if (nargin == 0)
    %if no input arguments run with a test case of a Fibonacci Lattice
    N = 40;
    M = floor(N/2);
    N = 2*M+1;
    % Fibonacci Lattice Points
    th =  pi/2 - asin(2*(-M:1:M)/N);
    phi =  4*pi*(-M:1:M)/(1+sqrt(5));
    % Spherical to cartesian
    x = sin(th).*cos(phi);  y = sin(th).*sin(phi);  z = cos(th);
    close all
    doPlot = 1;
end


% Create Data Structure and Determine Triangulation

g.X = [x(:) y(:) z(:)];             % Points for Triangulation
g.tri = convhulln(g.X);             % Convex Hull - creates triangulation
g.Tr = triangulation(g.tri,g.X);    % Triangles 
g.P = incenter(g.Tr);               % Incenters
g.n = -faceNormal(g.Tr);            % Outward surface Normals
g.e1 = zeros(size(g.P));            % (e1,e2) is local orthogonal basis 
g.e2 = zeros(size(g.P));
g.NumTri = size(g.tri,1);           % Number of triangles (facets)
g.M = cell(size(g.NumTri,1));       % Matrix M for on surface test
g.BC = ones(g.NumTri,1);            % Boundary conditions
[g.CMC, circum_radii] = circumcenter(g.Tr); % Circumcenters and Circumradii
g.maxRad = max(circum_radii);       % Maximum circumradius
g.tP1 = zeros(g.NumTri, 3);         % Triangle vector side 1
g.tP2 = zeros(g.NumTri, 3);         % Triangle vector side 2
% Coding for g.BC (determined in SpherePoints)
% 0 - reflecting.
% 1 - absorbing.

for j = 1:g.NumTri
    
    % Get an orthogonal basis (e1,e2) in the plane each triangle.
    
    temp = g.X(g.tri(j,1),:)-g.P(j,:);
    g.e1(j,:) = temp/norm(temp);
    g.e2(j,:) = cross(g.e1(j,:),g.n(j,:));
    
    % Create matrix to test whether given point is in triangle
    % The point P is in T if for u = M*P, u(1)>0, u(2) >0 and u(1)+u(2)<1.
    
    q2 = g.X(g.tri(j,2),:) - g.X(g.tri(j,1),:);
    q3 = g.X(g.tri(j,3),:) - g.X(g.tri(j,1),:);
    
    g.M{j} = [dot(q3,q3) -dot(q2,q3); -dot(q2,q3) dot(q2,q2)] * [q2;q3]/( dot(q2,q2) *dot(q3,q3) - dot(q2,q3)^2 );
    g.tP1(j,:) = g.M{j}(1,:);
    g.tP2(j,:) = g.M{j}(2,:);
end

if (doPlot)
% Plot the surface
    s = 2;

    trisurf(g.tri(:,:),g.X(:,1),g.X(:,2),g.X(:,3),'FaceColor', 'w', 'faceAlpha', 0.8);
    axis equal;
    hold on;
% Plot the local coordinate axes
    
    quiver3(g.P(s,1),g.P(s,2),g.P(s,3),g.n(s,1),g.n(s,2),g.n(s,3),0.5, 'color','r','linewidth',2);
    quiver3(g.P(s,1),g.P(s,2),g.P(s,3),g.e1(s,1),g.e1(s,2),g.e1(s,3),0.5, 'color','b','linewidth',2);
    quiver3(g.P(s,1),g.P(s,2),g.P(s,3),g.e2(s,1),g.e2(s,2),g.e2(s,3),0.5, 'color','m','linewidth',2);
%
    xlabel('$x$','Interpreter','latex','fontsize',16);
    ylabel('$y$','Interpreter','latex','fontsize',16);
    set(gcf,'color','w');
    set(gca,'fontsize',16);
    axis off;
    
% %
% %    Plot incircles
% %      
% 
%   theta=linspace(0,2*pi,20);
% 
%   for iTri= 1:g.NumTri
% 
%         q1 = g.P(iTri,:)-g.X(g.tri(iTri,1),:);
%         q2 = g.P(iTri,:)-g.X(g.tri(iTri,2),:);
%         q3 = g.P(iTri,:)-g.X(g.tri(iTri,3),:);
%         
%         q4 = g.X(g.tri(iTri,3),:)-g.X(g.tri(iTri,1),:);
%         q5 = g.X(g.tri(iTri,1),:)-g.X(g.tri(iTri,2),:);
%         q6 = g.X(g.tri(iTri,2),:)-g.X(g.tri(iTri,3),:);
%         
%         rad1 = sqrt(sum(cross(q4,q1).^2,2))./sqrt(sum(q4.^2,2));
%         rad2 = sqrt(sum(cross(q5,q2).^2,2))./sqrt(sum(q5.^2,2));
%         rad3 = sqrt(sum(cross(q6,q3).^2,2))./sqrt(sum(q6.^2,2));
%         
%         r= min([rad1 rad2 rad3],[],2); % note rad1 rad2 and rad3 are equal to O(eps)
%         
%         C =  repmat(r.*cos(theta),[3,1])' .*repmat(g.e1(iTri,:),[numel(theta),1]) + repmat(r.*sin(theta),[3,1])' .*repmat(g.e2(iTri,:),[numel(theta),1]);
%         plot3(g.P(iTri,1)+C(:,1),g.P(iTri,2)+C(:,2),g.P(iTri,3)+C(:,3),'-r','linewidth',3)
%   end
%  %

hold off;

end

 %
 %   Compute largest triangle size
 %
   smax=0; 
   for iTri= 1:g.NumTri
         s1 = sum((g.X(g.tri(iTri,3),:)-g.X(g.tri(iTri,1),:)).^2);
         s2 = sum((g.X(g.tri(iTri,1),:)-g.X(g.tri(iTri,2),:)).^2);
         s3 = sum((g.X(g.tri(iTri,2),:)-g.X(g.tri(iTri,3),:)).^2);
         smax=max([smax,s1,s2,s3]);
   end
   dmax=sqrt(smax);
   thetamax=2*asin(dmax/2);
   g.thetamax = thetamax;

end


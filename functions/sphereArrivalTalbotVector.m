function [output]= sphereArrivalTalbotVector(Router)

%function [thetacout,tout]= sphereArrivalTalbotVector(Router)

% Vectorized code for pre-computing sphere reinjection.
% Assumes a inner sphere radius of 1.
% Called when a particle reaches a sphere of a radius R
% This code determines if a particle escapes and
%   (i) returns  ts=0 if it escapes otherwise
%   (ii) returns thetac= the azimuthal displacement, and ts=the elapsed time.

% NOTES: For the output returned:
% (1) Output is for unit inner sphere - must be scaled appropriately.
% (2) The polar angle must subsequently chosen randomly from 0 to 2*pi.

%%% Plot flag (for debugging) - doPlot=1 means plot
doPlot = 0;
%%% Default values - again primarily for debugging purpose
if (nargin == 0)
    Router = 3;
    doPlot = 1;
end

%%%  Code parameters
nTalbot = 36;  % Talbot quadrature points
Upts=400;  % # points in (0,1/R) which are mapped to reinjection times
Jmax = 40; % Number of Legendre modes
thetapts = 400; % Number of theta points on cdf interpolation grid
thetacpts = 200; % Number of theta sample

%%% Implementation of Talbot integration:
% Talbot discretization points for midpoint rule
thTalbotvec = [-nTalbot+0.5:nTalbot-0.5]'*pi/nTalbot;
% Modified Talbot parameters
sg = 1.2244; mu = 1.0034; nu = 0.5290; bt = 0.6407;
%%% Details of coding the inverse Laplace transform
% s = @(theta,t) (n/t)*(-sg + mu*theta.*cot(bt*theta) +1i*nu*theta);
% ds = @(theta,t) (n/t)*(mu*cot(bt*theta) - mu*bt*theta.*csc(bt*theta).^2 + 1i*nu);
%%% Inv LT of G(s) is g(t) = real(sum(1/(2*1i*n)*exp(sTal*time).*G(sTal).*dsTal)/(time)^2);

%%% Grid for precomputed times of reinjection
%   Note that in main code must implement probability of escape is 1-1/R.
dUpts=1/Upts;
Uvec =dUpts*([1:Upts]-1/2); % Eliminate edge cases via midpoint
tvec=0.25 *(Router-1)^2 ./ (erfcinv(Uvec).^2);

%
% Calculate the chi's for the Legendre expansion of the impact point
%
Jvec=[1:Jmax]';
[J,thTalbot,T] = ndgrid(Jvec,thTalbotvec,tvec);
% Precompute Bessels on Talbot grid points
sTal=nTalbot./T.*(-sg + mu*thTalbot.*cot(bt*thTalbot) +1i*nu*thTalbot);
dsTal=nTalbot./T.*(mu*cot(bt*thTalbot)- mu*bt*thTalbot.*csc(bt*thTalbot).^2 ...
    + 1i*nu);
% Very important! The factor of 1/sqrt(Router) is due to spherical Bessel def'n
G = 1/sqrt(Router)*besselk(J+0.5,sqrt(sTal)*Router)./besselk(J+0.5,sqrt(sTal));
%%% Talbot sum
Talbotsummand=exp(sTal.*T).*G.*dsTal/(2*1i*nTalbot);
chi_j=squeeze(real(sum(Talbotsummand,2)));
%%% Now computed cdf on a regular theta array
dtheta=pi/thetapts;
theta=dtheta*([1:thetapts]-1/2)';  % Use midpoints of intervals
[Theta,Tvec]=ndgrid(theta,tvec);
%%%
rho_ts = ((Router-1)./(2*Router*sqrt(pi)* Tvec.^(3/2)).* exp(-0.25 * (Router-1)^2./Tvec));
LP = legendrePVeryFast(Jmax+1,cos(theta)'); % Legendre array
LPdiff = LP(1:end-2,:) - LP(3:end,:);          % Legendre difference for sum
% Sum the chi's * Legendre factor (kk) as a vector product
cdf = 0.5*(1 - cos(Theta))+ 0.5*(LPdiff'*chi_j)./(rho_ts);
% The variable cdfinterp is 1-cdf which allows us to resolve the tail
% of the extreme arrivals as it asymptotes to zero, not one.
cdfinterp=0.5*(1 + cos(Theta))- 0.5*(LPdiff'*chi_j)./(rho_ts);


if (doPlot)
    clf
    figure('color','w'); ax = gca;
    plot(theta,cdf(:,1:50:400),'linewidth',2)
    set(gca,'fontsize',16);
    xlim([0 pi]);
    xticks([0 pi/2 pi]);
    xticklabels({'$0$', '$\pi/2$', '$\pi$'});
    yticks([0 0.2 0.4 0.6 0.8 1.0]);
    %xlabel('$\theta$','interpreter','latex'
    ylabel('$P[ \theta_{*} < \theta ]$','Interpreter','latex');
    xlabel('$\theta$','Interpreter','latex');
    ax.TickLabelInterpreter = 'latex';
    %exportgraphics(figure(12),'MulipleCDF.png','resolution',300);
end

%%% Finally, interpolate thetac values
dV=1/thetacpts;
Vvec =dV*([1:thetacpts]-1/2);

[thetacout,tout]=ndgrid(Vvec,tvec);

for k=1:Upts
    thetacout(:,k)=interp1(cdfinterp(:,k),theta,1-Vvec);
end


output=[thetacout(:),tout(:)];

end

% if (doPlot)
%     %
%     %  PLOT ROUTINE - First define some colors
%     %
%     lightGrey = 0.7*[1 1 1];
%     vlightGrey = 0.95*[1 1 1];
%     %
%     %   Figure 11 is point of impact on sphere
%     %
%     figure(1)
%     [x1,y1,z1] = sphere(40);
%     set(gcf,'color','w');
%     axis off;
%
%     hFig = surface(x1,y1,z1,'FaceColor',vlightGrey,'EdgeColor',lightGrey);
%     hFig.FaceAlpha = 0.5;
%     hFig.LineWidth = 0.2;
%
%     z = 2*pi*rand();
%     hold on
%     plot3(sin(thetac)*cos(z),sin(thetac)*sin(z),cos(thetac),'r.','markersize',32)
%     hold off
%
%     xlabel('x');
%     ylabel('y');
%     axis([-1 1 -1 1]);
%     view([2,2,2]);
%
%     figure(2)
%     %
%     % CDF for theta distribution
%     %
%     theta = linspace(0,pi,2^10);
%     cdf=F_cdf(theta,ts,chis,Router);
%     %
%     % If ts << 1 the distribution should be the same as arrival on the tangent
%     % plane to the north poll.
%     %
%     lfac=1;  %conversion factor between theta and distance from the polar axis
%
%     planeapprox= 1-exp(-(lfac*theta).^2/(4*ts));
%     hold on
%     plot(theta,F_cdf(theta,ts,chis,Router),'DisplayName','Full CDF');
%     %plot(theta,planeapprox,'DisplayName','Plane Approx');
%     legend
%     hold off
%     % mincdf=min(cdf) %%For debug should be zero.
%     % maxcdf=max(cdf) %%For debug should be one.
% end

function p = legendrePVeryFast(m,x)

%
%  Test for up to N=15 versus in-built Legendre and
%  Answer agres to within 5 x10^(-15)
%

p = zeros(m+1,length(x));

p(1,:) = 1;
p(2,:) = x;

for j = 3:m+1
    l = j-1;
    p(j,:) = (1/l) * ( (2*l-1) * x.* p(j-1,:) - (l-1)* p(j-2,:));
end

end
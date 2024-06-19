function PlanarArrivalGH

clear all; close all;

addpath('functions');

% This code uses KMC to compute the capture time histogram
% for a configuration of one or more pores. It can compare
% the results to asymptotics for:
%  (i)  Capacitance - compared to cpature probability
%  (ii) Small pore limit - computed either asymptotically or
%      as an inverse Laplace transform
% There are two types of initial conditions
%   (i)  Random on a hemisphere (Metric = 'Cap')
%   (ii) Point release (Metric =''Asymp')

% Two examples choices of pore configurations from paper.

%Ex = 'OnePore'; % One pore with capacticane calculation.
Ex = 'SixPore'; % Six pores with FPT calculation.

%%% Plot, Output, Parallelization Flags %%%
HistogramPlotFlag = true;
IterationPlotFlag = true;
%   rng(1);  % For deterministic random numbers

%%% LaTeX Axis Labels and titles
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');

%%% Initialization
np = 1e6;   % Number of particles
D = 1;   % Diffusion constant

%%% Initial condition Ex determines number of pores
switch Ex
    case 'OnePore'
        % Calculate capacitance
        Pore.x = [ 0 ];
        Pore.y = [ 0 ];
        Pore.r = [ 1 ];

        Metric = 'Cap';
        %% For capacitance, place all points on hemisphere.
        r0 = 5; th_start = 2*pi*rand(np,1);
        P.z = (-1 + 2*rand(np,1));
        P.x = r0*sqrt(1-P.z.^2) .* cos(th_start);
        P.y = r0*sqrt(1-P.z.^2) .* sin(th_start);
        P.z = r0*abs(P.z);

    case 'SixPore'
        Nr = 5; th = linspace(pi/2,3*pi/2,Nr);

        Pore.x = [ cos(th)  15];
        Pore.y = [ sin(th)  0.0];
        Pore.r = 10*[ 0.001*ones(1,Nr) 0.1];

        %% All particles start at the same point.
        P.x0 = 0;   P.y0 = 0.0;   P.z0 = 0;

        Metric = 'Asymp';

        P.x = P.x0*ones(np,1);
        P.y = P.y0*ones(np,1);
        P.z = P.z0*ones(np,1);
        P.D = D;

end

P.D = D;
P.t = zeros(np,1);       % Initialize times of capture.
Pore.n = numel(Pore.x);  % Number of Pores in Configuration.

% Facilitate vector calculation of particle/pore distances.
xP = repmat(Pore.x,[np,1]);
yP = repmat(Pore.y,[np,1]);
aP = repmat(Pore.r,[np,1]);

% Iteration clock and max number of iterations.
nIterMax = 10000;
nIter = 1;
IterHistory=zeros(nIterMax,1);

% Tolerance level for numerical error on the hemisphere escape time.
eps_htt = 1e-6;

% Increase trap radius by machine precision to avoid stuck particles.
eps_reg = 10*eps('double');  % eps=machine precision

% Generating the hemisphere transit time table
[taus, cdf_tau] = build_htt_table(eps_htt);

% Parameter definitions
p_large = 1-2*nthroot(eps_htt, 3);

% The index of particles still free.
indFree = 1:np;    % Index of free particles.
nFree = np;
IterHistory(1)=nFree;

% All particles start away from surface. Let S be a hemisphere centered at
% the origin with radius R0 that encloses the entire planar arrangment of pores.
% This is the smallest sphere that can be projected onto.
% Therefore we test for possibility of escape all points that are more that
% a radial distance R= Rfac*R0 from the origin.

Rfac = 4; %1000;

% R0 - radius of sphere enclosing all targets.
R0 = max(sqrt((Pore.x.^2 + Pore.y.^2))+ Pore.r)+eps_reg;
R = Rfac*R0;

% Pretabulate (hemi)sphere reinsertion data
outDraw = sphereArrivalTalbotVector(Rfac);
drawTh = outDraw(:,1); drawT = outDraw(:,2); numT = numel(drawTh);

repeat = true;

while (repeat)

    % At this stage all particles are in the bulk.
    % Project each to the plane z = 0.
    % Calculate the time and locations of impact on the plane.

    p = rand(nFree, 1);
    t = (1/(4*D)) * P.z(indFree).^2 .* (erfcinv(p).^(-2));
    r_sample = sqrt( 4*D*t .* log(1 ./ (1-rand(nFree, 1))) );
    th = 2*pi*rand(nFree, 1);

    P.x(indFree) = P.x(indFree) + r_sample.*cos(th);
    P.y(indFree) = P.y(indFree) + r_sample.*sin(th);
    P.z(indFree) = 0;
    P.t(indFree) = P.t(indFree) + t;

    % At this stage, all free particles are on the plane.
    % Divide into those that could escape and those that might be absorbed.

    d2 = P.x(indFree).^2 + P.y(indFree).^2;
    indCouldEsc = find(d2 > R^2);
    indCouldAbs = find(d2 <= R^2);

    indEsc = [];

    if ~isempty(indCouldEsc)
        UEsc = rand(size(indCouldEsc));

        indEsc = indFree(indCouldEsc(UEsc > 1./Rfac));  % Escaping
        indRem = indFree(indCouldEsc(UEsc <= 1./Rfac)); % Reinjected

        % Update any escaping particles with times to zero (removed from sim).

        if ~isempty(indEsc), P.t(indEsc) = 0; end

        % Those that remain need a new time and position.

        if ~isempty(indRem)

            % Update Positions and time.
            % The new radii is (1/Rfac) the previous.

            Rd = (1/Rfac)*sqrt(d2(indCouldEsc(UEsc <= 1./Rfac)));

            % Update times.

            indDraw = randi(numT,length(indRem),1);
            P.t(indRem) = P.t(indRem) + (Rd.^2/D) .* drawT(indDraw);
            thEsc = drawTh(indDraw);

            % Get spherical coordinates

            [Thd,Phid,~] = cart2sph(P.x(indRem),P.y(indRem),P.z(indRem));
            Phid = pi/2 - Phid;
            Z = 2*pi*rand(size(indRem))';

            ZaxisVec=repmat([0;0;1],1,length(indRem));
            V2= rotzVec(Thd,...
                rotyVec(Phid,...
                rotzVec(Z,...
                rotyVec(thEsc,...
                ZaxisVec))))*spdiags(Rd,0,length(Rd),length(Rd));
            P.x(indRem) = V2(1,:)';
            P.y(indRem) = V2(2,:)';
            P.z(indRem) = abs(V2(3,:))';

        end

    end

    % Create an indicator array for particles captured this step
    indCap = [];

    if ~isempty(indCouldAbs)
        % Get distances to the pores
        indCouldAbs = indFree(indCouldAbs);

        r = sqrt( (xP(indCouldAbs,:) - repmat(P.x(indCouldAbs),[1,Pore.n])).^2 + ...
            (yP(indCouldAbs,:) - repmat(P.y(indCouldAbs),[1,Pore.n])).^2 );

        d = r - aP(indCouldAbs,:);

        indRem = find(sum(d < eps_reg , 2) == 0 );    % The index of particles which remain free.

        hr = min(d,[],2);  % Calculate closest distance to a pore.

        hr = hr(indRem);
        indRem = indCouldAbs(indRem);   % Indices of remaining particles.
        indCap = setdiff(indCouldAbs,indRem);   % Indices of captured particles.

        if (~isempty(indRem))

            % Advance the non absorbed planar points to hemispheres.
            p = rand(numel(indRem), 1);
            tau = -1*ones(numel(indRem), 1);
            tau(p > p_large) = log(2 ./ (1-p(p>p_large)) );
            %tau(p < p_small) = -pi^2/ (2*lambertw(-1, -(pi/8)*p(p < p_small).^2));
            tau(tau < 0) = interp1(cdf_tau, taus, p(tau < 0));

            t = (1/pi^2) * hr.^2 .* tau/D;

            %iterations_to_bind(indFree) = iterations_to_bind(indFree) + 1;
            P.t(indRem) = P.t(indRem) + t;

            % Choose a random point on the hemisphere
            th = 2*pi*rand(numel(indRem) , 1);
            phi = acos(2*rand(numel(indRem) , 1) - 1);

            P.x(indRem) = P.x(indRem) + hr.*sin(phi).*cos(th);
            P.y(indRem) = P.y(indRem) + hr.*sin(phi).*sin(th);
            P.z(indRem) = abs(hr.*cos(phi));

        end

    end

    indFree = setdiff(indFree,[indCap(:); indEsc(:)] );
    nFree = numel(indFree);
    IterHistory(nIter)=nFree;

    nIter = nIter + 1;

    repeat = ((nFree>0) || (nIter>nIterMax));

    if (nIter > nIterMax)
        disp('Max iterations exceeded.');
    end

end

%%% Plots and metrics below %%%

% Histogram of capture pdf

if (HistogramPlotFlag)
    figure('color','w');
    hold on
    num_bins = 50;

    t = logspace(-2,8,num_bins);
    times = P.t(P.t>0); frac = length(times)/np;
    histogram(times, t, 'Normalization', 'probability','DisplayName','Capture Times');

    if strcmp(Metric,'Asymp')

        plot_fac = 5; t = logspace(-2,8,plot_fac*num_bins);
        totalcdf = sum(FulldistPlanarExact(t,P,Pore),1);
        tplot = 0.5*(t(2:end) + t(1:end-1));

        plot(tplot,plot_fac*diff(totalcdf)/frac,'-r','linewidth',3,DisplayName='Asymptotic density');
        
    end
    set(gca,'XScale','log'); set(gca,'Fontsize',16);
    xlabel('time')
    lgd = legend;
    lgd.FontSize = 16;
    lgd.TextColor = 'black';
    lgd.Interpreter = 'latex';
    hold off
end

if strcmp(Metric,'Cap')

    %%% Estimate of the capacitance from the initial radius
    %%% and the capture percentage compted via KMC
    CaptureProbKMC = length(times)/np;
    Cap = r0*CaptureProbKMC;

    %%% Estimate standard deviation.
    %%% For np particles with a probability of capture p

    stderr = sqrt((1-CaptureProbKMC)/(CaptureProbKMC*np));

    % Output Capacitance(KMC), Capacitance (Asymptotic)
    fprintf(1,['Capacitance (KMC) = %4.3e,' ...
        'Standard error = %4.3e \n'],Cap, stderr);

end

if(IterationPlotFlag)
    figure('color','w')
    hold on
    plot(IterHistory,'+')
    set(gca, 'YScale', 'log'); set(gca,'fontsize',16)
    xlabel('\# Iterations', 'FontSize', 16)
    ylabel('\# particles free', 'FontSize', 16)
    hold off
end




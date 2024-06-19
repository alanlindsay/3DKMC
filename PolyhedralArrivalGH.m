function PolyhedraArrivalGH

% rng(1); % For deterministic runs

addpath('functions');

close all

%% Initialization
np = 1e5; % Number of particles.
D = 1;  % Diffusion constant.
phi = pi/2; r = 2.5; % Start point (spherical).
x0 = r*cos(phi); y0 = 0; z0 = r*sin(phi); % Start point (Cartesian).

%%% Flags for metrics
plotSurface = 1;
plotHistogram = 1;
plotSplit = 1;
printIters = 1;

%% Choose a domain switch
     domain = 'sphere';
    %domain = 'ellipsoid';

%%% Compute triangulated surface based on domain switch
switch domain
    case 'sphere'

        %% This code creates the geometry from scratch
        %  NPore = 51; % Number of receptors.
        %  sigma = 0.1; % Surface absorbing fraction.
        %  g = SpherePoints(NPore,sigma,0,'SphereGeometry.mat');

         % A precomputed geometry.
          load gfile_10fix_2500move.mat
       
    case 'ellipsoid'

        %   At present pores are hardwired to two at N/S pole
        %   NPoly = <number of points in each hoop>
        %   NRing = There are 2*NRing +1 slices in parametrization
        %   Rpore = Pore radius
        %   Zpore = Pore height
        %   Requator = Equator radius
        %   delta =Pore center/equator desingularization parameter

        NPoly = 12;
        NRing = 5;
        Rpore = 1; a = Rpore;
        Zpore = 1;

        Requator = 4;

        delta = 0.01;

        doPlot = 0;  % Plotting is on

        %%% Uncomment line below to generate geometry
        g = EllipsoidPointsFixed(NPoly, NRing, Rpore, Zpore, Requator, delta, doPlot,'EllipseTest.mat');
end

% Radius of largest ball enclosing all surfaces.
L = 1.1*max( sqrt(g.X(:,1).^2 + g.X(:,2).^2 + g.X(:,3).^2 ));

% Build out intial condition
x = x0*ones(np,1);
y = y0*ones(np,1);
z = z0*ones(np,1);

% Iteration clock and max number of iterations.
nIterMax = 1000000;
nIter = 0;

% Tolerance level for numerical error on the hemisphere escape time.
epsHemi = 1e-4;

% Lower bound for radial displacement on surface
% This deals with points caught on an edge
epsRad = 1e-6;

% Generating the hemisphere transit time table
[taus, cdf_tau] = build_htt_table(epsHemi);

% Asymptotic bound  for hemisphere transit time
p_large = 1-2*nthroot(epsHemi, 3);

%%% The number of interations and times till capture.
%   its_to_bind = zeros(np,1);
time_to_bind = zeros(np,1);
face_bind = zeros(np,1); % pre-allocate face capture indexstorage

%%% The index of particles still free.
indFree = 1:np;    % Index of free particles.
nFree = np;        % Number of free particles

% All particles start away from surface.
% Let S be a sphere of radius R that encloses the surface
% For a cube  R= \sqrt{3} L  which encloses the entire cube.
% This is the smallest sphere that can be projected onto.
% Therefore we test for possibility of escape
% all points that are more that Rfac times this radius away.

Rfac = 3;

outDraw = sphereArrivalTalbotVector(Rfac);
drawTh = outDraw(:,1); drawT = outDraw(:,2); numT = numel(drawTh);

R = Rfac*L;

repeat = true;

% Prep matrix for in-facet calculation. M_facet is 2x3xNumTri.
M_facet = reshape(cell2mat(g.M),2,3,g.NumTri);
M_facet1 =  squeeze(M_facet(1,:,:));
M_facet2 =  squeeze(M_facet(2,:,:));

%%% Main KMC loop

while (repeat)

    %%% If any points are outside the large sphere,
    %%% Decide if they escape to infinity or reinject

    d2 = x(indFree).^2 + y(indFree).^2 + z(indFree).^2;

    indCouldEsc = find(d2 > R^2);

    if ~isempty(indCouldEsc)

        UEsc = rand(size(indCouldEsc));

        indEsc = indFree(indCouldEsc(UEsc > 1./Rfac));
        indRem = indFree(indCouldEsc(UEsc <= 1./Rfac));

        % Those escaping have zero arrival time.

        if ~isempty(indEsc)
            time_to_bind(indEsc) = 0;
        end

        %%% Update positions of reinjected points
        if ~isempty(indRem)
            indDraw = randi(numT,length(indRem),1);

            % The new radii is reduced by a factor Rfac.

            Rd = (1/Rfac)*sqrt(d2(indCouldEsc(UEsc <= 1./Rfac)));

            % Update times.

            time_to_bind(indRem) = time_to_bind(indRem) + (Rd.^2/D).*drawT(indDraw);

            % Get spherical coordinates

            thEsc = drawTh(indDraw);
            [Thd,Phid,~] = cart2sph(x(indRem),y(indRem),z(indRem));
            Phid = pi/2 - Phid;
            Z = 2*pi*rand(size(indRem))';

            ZaxisVec=repmat([0;0;1],1,length(indRem));
            temp = rotyVec(thEsc,ZaxisVec);
            temp = rotzVec(Z,temp);
            temp = rotyVec(Phid,temp);
            temp = rotzVec(Thd,temp);
            V2=temp*spdiags(Rd,0,length(Rd),length(Rd));

            x(indRem) = V2(1,:)';
            y(indRem) = V2(2,:)';
            z(indRem) = V2(3,:)';


        end
        indFree = setdiff(indFree,indEsc);
        %indFree(indCouldEsc(tEsc == 0)) = [];
        nFree = numel(indFree);

    end

    %%% At this stage each free particle is in the bulk.
    %%% Calculate the distance of each particle to the closest face.
    %%% Note this is the slowest part of the calculation

    if (~isempty(indFree))

        % Calculate the distances to each incenter.

        V = [x(indFree), y(indFree), z(indFree)];

        % Calculate signed distance to each plane.

        dV2 = (g.n * V' - repmat(sum(g.n .* g.P,2),[1,nFree]))';

        % Choose the plane with the smallest positive distance.

        [minD,indPlane] = max(dV2,[],2);

        % Project to the plane.

        p = rand(nFree, 1);

        t = (1/(4*D)) * minD.^2 .* (erfcinv(p).^(-2));
        r_sample = sqrt( 4*D*t .* log(1 ./ (1-rand(nFree, 1))) );
        % Avoid edge trap problems by taking minimum radial jump.
        r_sample = max(r_sample,epsRad);

        th = 2*pi*rand(nFree, 1);

        time_to_bind(indFree) = time_to_bind(indFree) + t;
        V = V - g.n(indPlane,:) .* repmat(minD,[1,3]) ...
            + repmat(r_sample.*cos(th),[1,3]) .*g.e1(indPlane,:) ...
            + repmat(r_sample.*sin(th),[1,3]).*g.e2(indPlane,:);

        % Vectorized test to see if point lands in a face.

        tp1 = dot(M_facet1(:,indPlane),(V-g.X(g.tri(indPlane,1),:))')';
        tp2 = dot(M_facet2(:,indPlane),(V-g.X(g.tri(indPlane,1),:))')';
        indFace = find( (tp1 > 0).* (tp2 >0) .*(tp1 + tp2 < 1) );

        indReflect = indFace(g.BC(indPlane(indFace)) == 0);
        indAbsorb = indFace(g.BC(indPlane(indFace)) ~= 0);

        if ~isempty(indAbsorb)
            % Note the face of arrival for each absorbed particle.
            iA = indFree(indAbsorb);
            face_bind(iA) = indPlane(indAbsorb);
        end

        % fraction of particles on reflecting faces
        fracR = numel(indReflect)/numel(indFree);

        while ~isempty(indReflect)

            % Face is reflecting. Project particle onto hemisphere.

            it = indPlane(indReflect);
            iR = indFree(indReflect);

            p = rand(size(indReflect));
            tau = -ones(size(p));
            tau(p > p_large) = log(2 ./ (1-p(p>p_large)) );
            tau(tau < 0) = interp1(cdf_tau, taus, p(tau < 0));

            % The largest possible radius for hemisphere.

            q1 = V(indReflect,:)-g.X(g.tri(it,1),:);
            q2 = V(indReflect,:)-g.X(g.tri(it,2),:);
            q3 = V(indReflect,:)-g.X(g.tri(it,3),:);

            q4 = g.X(g.tri(it,3),:)-g.X(g.tri(it,1),:);
            q5 = g.X(g.tri(it,1),:)-g.X(g.tri(it,2),:);
            q6 = g.X(g.tri(it,2),:)-g.X(g.tri(it,3),:);

            rad1 = sqrt(sum(cross(q4,q1).^2,2))./sqrt(sum(q4.^2,2));
            rad2 = sqrt(sum(cross(q5,q2).^2,2))./sqrt(sum(q5.^2,2));
            rad3 = sqrt(sum(cross(q6,q3).^2,2))./sqrt(sum(q6.^2,2));

            rad = min([rad1 rad2 rad3],[],2);

            % fprintf(1,'Hemisphere radius: %4.3e \n', rad);
            % DIFFUSION CONSTANT matters here.
            time_to_bind(iR) = time_to_bind(iR) + (1/pi^2) * (rad.^2 /D) .* tau;

            % Choose a random point on the hemisphere.
            thH = 2*pi*rand(size(indReflect));
            zH = rand(size(indReflect));

            new_x = rad .* sqrt(1-zH.^2).*cos(thH);
            new_y = rad .* sqrt(1-zH.^2).*sin(thH);
            new_z = rad .* zH;

            V(indReflect,:) = V(indReflect,:) + g.n(it,:).*repmat(new_z,[1,3]) + ...
                g.e1(it,:).*repmat(new_x,[1,3]) +  g.e2(it,:).*repmat(new_y,[1,3]);

            % Here we repeat the projection between bulk and facet until
            % another facet or bulk is reached.

            p = rand(numel(iR), 1);

            t = (1/(4*D)) * new_z.^2 .* (erfcinv(p).^(-2));
            r_sample = sqrt( 4*D*t .* log(1 ./ (1-rand(numel(iR), 1))) );
            % Avoid edge trap problems by taking minimum radial jump.
            r_sample = max(r_sample,epsRad);

            th = 2*pi*rand(numel(iR), 1);

            time_to_bind(iR) = time_to_bind(iR) + t;
            V(indReflect,:)= V(indReflect,:) - g.n(it,:) .* repmat(new_z,[1,3]) + repmat(r_sample.*cos(th),[1,3]) .*g.e1(it,:) + repmat(r_sample.*sin(th),[1,3]).*g.e2(it,:);

            % Check again for facet or bulk.

            tp1 = dot(M_facet1(:,it),(V(indReflect,:)-g.X(g.tri(it,1),:))')';
            tp2 = dot(M_facet2(:,it),(V(indReflect,:)-g.X(g.tri(it,1),:))')';

            indReflect( ~((tp1 > 0) .* (tp2 >0) .* (tp1 + tp2 < 1)) ) = [];

        end
        %%% Progress print loop
        if (printIters && mod(nIter,100) == 0)
            fprintf(1,'Remaining Particles : %g, Fraction Reflecting = %3.2f \n',nFree,100*fracR);
        end

        % Now remove particles that have reached absorbing surfaces.

        x(indFree) = V(:,1);
        y(indFree) = V(:,2);
        z(indFree) = V(:,3);

        indFree(indAbsorb) = [];
        nFree = numel(indFree);

    end

    repeat = ((nFree>0) && (nIter <= nIterMax) );

    if (printIters && nIter > nIterMax)
        disp('Max iterations exceeded.');
    end

    nIter = nIter + 1;

end

%%% End of main KMC loop

%%% Surface Plotter
if (plotSurface)
    figure('color','w')
    hold on;
    trisurf(g.tri(g.BC==0,:),g.X(:,1),g.X(:,2),g.X(:,3),'FaceColor', 'white', 'faceAlpha', 0.6);
    trisurf(g.tri(g.BC~=0,:),g.X(:,1),g.X(:,2),g.X(:,3),'FaceColor', 'blue', 'faceAlpha', 0.6);
    axis equal; axis off;
    hold off;
end

%%% Capacitance assuming point release

t = time_to_bind(time_to_bind>0); % Capture times
p = length(t)/np;
C = sqrt(x0^2 + y0^2 + z0^2)*p;
err = sqrt((1-p)/(p*np));
fprintf(1,'Surface capacitance =  %8.6f, KMC error = %8.6f, Capture Prob = %5.4f. \n',C,err,p)

%%% Histogram of capture distribution.

if (plotHistogram)

    figure('color','w');
    hold on;
    bins = logspace(-2,7,20);
    histogram(t, bins, 'Normalization', 'probability','DisplayName','Capture Times');
    xscale log;  ax = gca; xlim([10^(-2) 10^6]);
    xlabel('time'); ylabel('Capture distribution');
    ax.YAxis.TickLabelFormat = '%,.2f';
    ax.TickLabelInterpreter = 'latex';
    ax.FontSize = 16;
    hold off;

end

%%% Plot the dynamic fluxes between the top and bottom of the cell.
if (plotSplit)
    
    t_capture = sort(time_to_bind(time_to_bind>0)); % Capture times
    face_capture = face_bind(time_to_bind>0); % Capture faces
    split_capture = sign(g.P(face_capture,3));% Capture flag +/- 1
    top_capture =   t_capture(split_capture==1);
    bot_capture = t_capture(split_capture==-1);
    %
    Ncap=numel(t_capture);
    Ntop=numel(top_capture);
    Nbot=numel(bot_capture);

    % Empirical CDF of top and bottom captures.
    top_cdf= (1:Ntop)/Ncap;
    bot_cdf= (1:Nbot)/Ncap;
    %
    capture_ratio=Nbot/Ntop;
    %
    figure('color','w')
    hold on
    h1 = plot(top_capture,top_cdf,'linewidth',3,'DisplayName','Top capture');
    h2 = plot(bot_capture,bot_cdf,'linewidth',3,'DisplayName','Bottom capture');
    xscale log;
    ax = gca;
    xlabel('time','FontSize',16,'Interpreter','latex');
    ylabel('Probability','FontSize',16,'Interpreter','latex');
    set(gca, 'YTick', [0:0.1:1])
    ylim([0,1]); xlim([10^(-2) 10^6]);
    ax.YAxis.TickLabelFormat = '%,.1f';
    ax.TickLabelInterpreter = 'latex';
    ax.FontSize = 16; lgd = legend;
    hold off
end
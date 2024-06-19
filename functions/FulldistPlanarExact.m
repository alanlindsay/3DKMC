function out = FulldistPlanarExact(t,P,Pore)

%close all;

if (nargin == 0)
    load data/PlanarArrivalN1e6.mat;
    
    times = P.t(P.t>0);
    D = 1;
    p = numel(times)/numel(P.t);
    
    bins = 200;
    bin_fac = 5;
    t = linspace(log10(1+min(times)), log10(1+max(times)), bin_fac*bins);
end

D = P.D;
c = (2*Pore.r')/pi;
cdf0 = zeros(length(c),length(t));
cdf1 = zeros(length(c),length(t));

R = sqrt( (Pore.x-P.x0).^2 + (Pore.y-P.y0).^2 + P.z0^2 );

d = sqrt( (Pore.x-Pore.x').^2 + (Pore.y-Pore.y').^2 );

for j = 1:Pore.n
    cdf0(j,:) = (c(j)/R(j))*erfc(R(j)./sqrt(4*D*t));
    cdf1(j,:) = (c(j)^2/R(j)) *exp(-R(j)^2 ./(4*D*t))./sqrt(pi*D*t);
    for k = 1:Pore.n          
        if (k ~= j) 
            cdf1(j,:) = cdf1(j,:) - (c(j)*c(k)/(d(j,k)*R(k)))*erfc( (R(k)+d(j,k))./sqrt(4*D*t) );
        end
    end
end

%plot(t,sum(flux,1))
out = cdf0 + 1*cdf1;
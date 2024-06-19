function [taus, cdf_tau] = build_htt_table(epsilon)
    % BUILD_HTT_TABLE Build the lookup table for the hemisphere transit time
    % function
    
    %%% Choose min_tau as machine eps 
    % (was Choose our minimum so that the first cdf value is 10^(-200))
    min_tau = 0.005; 
    % Choose our maximum based on our asymptotic approximations
    max_tau = log(2/epsilon);
    
    % Choose table size based on linear interpolation error
    step_size = sqrt(epsilon);
    taus = min_tau:step_size:max_tau;
    table_size = length(taus);
    
    cdf_tau = zeros(table_size, 1);
    
    % Choose cutoff point for theta function inversion
    cutoff = 1;
    cutoff_i = round((1/step_size)*(cutoff-min_tau));
    
    %
    % Smallest term is roughly exp(-num_terms^2)
    %
    num_terms = 25;
    
    for i = 1:cutoff_i
        t = taus(i);
        iexpf = @(x) exp(-(x+.5)^2*(pi^2/t));
        iexps = arrayfun(iexpf, num_terms:-1:0);
        invf = 2*sqrt(pi)*(1/sqrt(t))*sum(iexps);
        cdf_tau(i) = invf;
    end
    for i = cutoff_i:table_size
        t = taus(i);
        expf = @(n) (-1)^n * exp(-n^2*t);
        exps = arrayfun(expf, num_terms:-1:1);
        f = 1+2*sum(exps);
        cdf_tau(i) = f;
    end
end
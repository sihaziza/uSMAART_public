
%%

function [w,Wall] = estimateFilterReg(x1, x2, dn, dn_overlap, varargin)
%     if(isempty(oneway)) oneway = false; end
        
    if(mod(dn,2)) dn = dn+1; end
    
    options = defaultOptions();
    if(~isempty(varargin))
        options = getOptions(options, varargin);
    end
    
    if(isfinite(options.max_delay))
        options.max_phase = ...
            min([linspace(0, options.max_delay*pi, dn/2), pi, ...
                 linspace(options.max_delay*pi, 0, dn/2) ], ...
                options.max_phase)';
        options.max_phase(options.max_phase>pi) = pi;
    end
    
    nt = length(x1);
    L = dn + 1;
    T = [-(L/2+1):1:(L/2-2)];
    W = hann(L);
    nsteps =  1:dn_overlap:(nt-dn) ;
    Wall = nan([L, length(nsteps)]);
    
    fs = ifftshift(linspace(0,1,L)-1/2); 
    
    U1fall = nan([L, length(nsteps)]);
    U2fall = nan([L, length(nsteps)]);
    for i_n = 1:length(nsteps)
        n0 = nsteps(i_n);
        uz1 = x1(n0:(n0+dn));
        uz2 = x2(n0:(n0+dn));
                
        U1fall(:, i_n) = fft(W.*uz1,L); 
        U2fall(:, i_n) = fft(W.*uz2,L);
    end

    u1f = sum(U1fall.*conj(U2fall),2);
    u2f = sum(U2fall.*conj(U2fall),2);

    s = u1f./(u2f+options.eps*mean(u2f));

    fs = linspace(0,1,length(s));
    if(~isempty(options.fref))
        [~,ind_fref] = min(abs(fs-options.fref));
        options.max_amp = options.max_amp_rel*abs(s(ind_fref));
    end
    
    if(~isempty(options.max_amp))
        a_large = (abs(s) > abs(options.max_amp));
        in_flim = (fs'<=options.flim_max & fs'<=0.5) | (flip(fs')<=options.flim_max & flip(fs')<=0.5);
        s(a_large&in_flim) = options.max_amp*(s(a_large&in_flim)./abs(s(a_large&in_flim)));
    end

    % projecting s(f) on max_phase(f) direction if phase is too big
    phase_too_large_p = angle(s) >  options.max_phase;
    phase_too_large_n = angle(s) < -options.max_phase;
    s(phase_too_large_p) = abs(s(phase_too_large_p)).*...
        max(cos(angle(s(phase_too_large_p))-options.max_phase(phase_too_large_p)), 0).*...
        exp(1.i*options.max_phase(phase_too_large_p));
    s(phase_too_large_n) = abs(s(phase_too_large_n)).*...
        max(cos(-angle(s(phase_too_large_n))-options.max_phase(phase_too_large_n)), 0).*...
        exp(-1.i*options.max_phase(phase_too_large_n));
%     s(phase_too_large) = 0;

    w = real(ifft(s));
    w = [w(ceil(L/2):L); w(1:floor(L/2)-1)]; w = [w(2:end); w(1)]; 
%     w = fftshift(w);
end

function options = defaultOptions()
    options.eps = 1e-8; % can be used to suppres everything below poisson noise level

    options.max_amp = []; % maximimal spectral amplitude of a filter (= the spectral amplitude of a \delta-filter v(t): v(0)=max_amp, v(t~=t0)=0)
    options.flim_max = 1; % maximal frequency to apply amplitude correction, \in [0,0.5]

    options.fref = []; % reference frequency to determine max_amp, overwrites max_amp, \in [0,0.5]
    options.max_amp_rel = 1.2; % scaling of the options.fref amplitude determine max_amp
    
    options.max_phase = pi; % max phase delay (at every frequency), [-\pi,pi]
    options.max_delay = Inf; % max time delay (relative: max_delay(s)*fps(Hz)), overwrites options.max_phase
end

function [w,Wall] = estimateFilter(x1, x2, dn, dn_overlap, eps, reg_func, max_phase)
        
    if(mod(dn,2)) dn = dn+1; end

    nt = length(x1);
    Wall = [];

    for n0 = 1:dn_overlap:(nt-dn)
        x1in = x1(n0:(n0+dn));
        x2in = x2(n0:(n0+dn));
        
        [t,w1, s] = deconvolution(x1in, x2in, eps, reg_func, max_phase);
        Wall = [Wall, w1];
    end
    
    w = mean(Wall, 2);
end
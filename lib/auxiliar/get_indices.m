% 
ind.fbb(1)   = 1;
ind.fbb(2)   = ind.fbb(1)-1   + n_fb_max;
%
ind.t0       = ind.fbb(2) + 1;
%
ind.ToF(1) = ind.t0 + 1;
ind.ToF(2) = ind.ToF(1) + n_fb_max;
%


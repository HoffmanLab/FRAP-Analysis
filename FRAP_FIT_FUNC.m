function [f]  = FRAP_FIT_FUNC(x,a,b,c)
f = b-(b-a)*exp(-c*x);
end
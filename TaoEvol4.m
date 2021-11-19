function [pl,ql,yl,xl,t]=TaoEvol4(F,p0,q0,y0,x0,t0,T,omega)
% 'omega' is the parameter that contrals the error.
l=4;
gammal=1/(2-2^(1/(l+1)));
[pl,ql,yl,xl,t]=TaoEvol2(F,p0,q0,y0,x0,t0,T*gammal,omega);
[pl,ql,yl,xl,t]=TaoEvol2(F,pl,ql,yl,xl,t,(1-2*gammal)*T,omega);
[pl,ql,yl,xl,t]=TaoEvol2(F,pl,ql,yl,xl,t,gammal*T,omega);
end
function [p,q,tspan]=HamiltonTaoSymp(H,p0,q0,t0,tspan,omega,emi)
F=@(pp,qq,tt) HamiltonEvol(H,pp,qq,tt,10^(-9));

[p,q,tspan]=TaoSymp(F,p0,q0,t0,omega,tspan,emi);
end
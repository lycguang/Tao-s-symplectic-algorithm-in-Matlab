function dd=PartialDiff(H,x,n,emi)
% emi=10^(-5);
H1=H(x);
x(n,1)=x(n,1)+emi;
H2=H(x);
dd=(H2-H1)/emi;
end
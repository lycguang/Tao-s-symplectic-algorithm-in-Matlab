function F=HamiltonEvol(H,p,q,t,emiT)
% This function get the "F(x)" in evolution equation "dx/dt=F(x)" from the
% Hamiltonian 'H' with the initial point 'p,q,t'. 'emiT' is the
% time-interval for partial derivative. To maintain accuracy, 'emiT' should
% be set samller than the time-period we use.
NFreedom=length(p);
x=cat(1,p,q,t);
F=zeros(2*NFreedom,1);
for i=1:NFreedom
    F(i,1)=-PartialDiff(H,x,NFreedom+i,emiT);
    F(i+NFreedom,1)=PartialDiff(H,x,i,emiT);
end
end
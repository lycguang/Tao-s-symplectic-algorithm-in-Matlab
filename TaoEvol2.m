function [pl,ql,yl,xl,t]=TaoEvol2(F,p0,q0,y0,x0,t0,T,omega)
% 'omega' is the parameter that contrals the error.
NFreedom=length(q0);

Delta1=F(y0,q0,t0);
p=p0+T/2*Delta1(1:NFreedom,1);
x=x0+T/2*Delta1((1+NFreedom):(2*NFreedom),1);
Delta2=F(p,x,t0);
q=q0+T/2*Delta2((1+NFreedom):(2*NFreedom),1);
y=y0+T/2*Delta2(1:NFreedom,1);

R(1)=sin(2*omega*T);
R(2)=cos(2*omega*T);
q1=(q+x+R(2)*(q-x)+R(1)*(p-y))/2;
p1=(p+y-R(1)*(q-x)+R(2)*(p-y))/2;
x1=(q+x-R(2)*(q-x)-R(1)*(p-y))/2;
y1=(p+y+R(1)*(q-x)-R(2)*(p-y))/2;
q=q1;
p=p1;
x=x1;
y=y1;
clear q1 p1 x1 y1

Delta2=F(p,x,t0+T/2);
ql=q+T/2*Delta2((1+NFreedom):(2*NFreedom),1);
yl=y+T/2*Delta2(1:NFreedom,1);
Delta1=F(yl,ql,t0+T/2);
pl=p+T/2*Delta1(1:NFreedom,1);
xl=x+T/2*Delta1((1+NFreedom):(2*NFreedom),1);

t=t0+T;
end
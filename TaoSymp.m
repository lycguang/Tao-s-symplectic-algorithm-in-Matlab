function [p,q,tspan]=TaoSymp(F,p0,q0,t0,omega,tspan,emi)
T=tspan(1,2)-tspan(1,1);
T=chooseT(F,p0,q0,t0,T,omega,emi),

p=zeros(length(p0),length(tspan));
q=zeros(length(q0),length(tspan));
% Deltaqx=zeros(length(tspan),length(q0));
% Deltapy=zeros(length(tspan),length(q0));

for i=1:length(tspan)
    if i==1
        jceil=ceil((tspan(1,i)-t0)/T);
        tj=t0;
        pj=p0;qj=q0;
        yj=p0;xj=q0;
    else
        jceil=ceil((tspan(1,i)-tspan(1,i-1))/T);
        tj=tspan(1,i-1);
    end
    for j=1:jceil
        if j==jceil
            if abs(tspan(1,i)-tj)<T/10^8
            	tj=tspan(1,i);
            else
                [pj,qj,yj,xj,tj]=TaoEvol2(F,pj,qj,yj,xj,tj,tspan(1,i)-tj,omega);
            end
        else
            [pj,qj,yj,xj,tj]=TaoEvol2(F,pj,qj,yj,xj,tj,T,omega);
        end
    end
    p(:,i)=(pj+yj)/2;
    q(:,i)=(qj+xj)/2;
%     Deltaqx(i,:)=qj-xj;
%     Deltapy(i,:)=pj-yj;
%     i,
end
end

function T=chooseT(F,p0,q0,t0,T0,omega,emi)
T1=T0;T2=T1/2;
[p1,q1,y1,x1,~]=TaoEvol4(F,p0,q0,p0,q0,t0,T1,omega);
[p2,q2,y2,x2,~]=TaoEvol4(F,p0,q0,p0,q0,t0,T2,omega);
[p2,q2,y2,x2,~]=TaoEvol4(F,p2,q2,y2,x2,t0+T2,T2,omega);
emi1=max(abs(cat(2,p1-p2,q1-q2,y1-y2,x1-x2)));
% i=0;
while emi1>emi
    T1=T2;
    T2=T1/2;
    p1=p2;
    q1=q2;
    [p2,q2,y2,x2,~]=TaoEvol4(F,p0,q0,p0,q0,t0,T2,omega);
    [p2,q2,y2,x2,~]=TaoEvol4(F,p2,q2,y2,x2,t0+T2,T2,omega);
    emi1=max(abs(cat(2,p1-p2,q1-q2,y1-y2,x1-x2)));
%     i=i+1,
end
T=T1;
end

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

function [pl,ql,yl,xl,t]=TaoEvol4(F,p0,q0,y0,x0,t0,T,omega)
% 'omega' is the parameter that contrals the error.
l=4;
gammal=1/(2-2^(1/(l+1)));
[pl,ql,yl,xl,t]=TaoEvol2(F,p0,q0,y0,x0,t0,T*gammal,omega);
[pl,ql,yl,xl,t]=TaoEvol2(F,pl,ql,yl,xl,t,(1-2*gammal)*T,omega);
[pl,ql,yl,xl,t]=TaoEvol2(F,pl,ql,yl,xl,t,gammal*T,omega);
end
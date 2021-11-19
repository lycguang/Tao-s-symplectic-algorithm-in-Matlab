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
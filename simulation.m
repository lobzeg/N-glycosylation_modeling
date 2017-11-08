function [hm,complex,hybrid] = simulation(b,Vm,Km,concentration,NumSim,NumAdj,NE)
    
ChangeCompartmentRate = 1 / 300;
div = 0.95;
complex = 0;
hm = 0;

for i = 1:4    
    for j = 1:NE
        s(j,i) = 100;
    end;    
    ss(i) = 100 * (NE - 1);
end;

for i = 1:NumAdj
    [s,ss] = adjustment(b,s,ss,Vm,Km,ChangeCompartmentRate,concentration,div,NE);
end;

for i = 1:4
    for j = 1:NE
        s(j,i) = s(j,i) / ss(i);
    end;    
end;

for i = 1:NumSim
    a = golgi_sim_full(b,s,Vm,Km,ChangeCompartmentRate,concentration,NE);
    if a(2) + a(3) + a(4)>1 && a(7) == 0
        hm = hm + 1;
    end;
    if a(2) + a(3) + a(4) == 0
        complex = complex + 1;
    end;
end;

hm = hm / NumSim;
complex = complex / NumSim;
hybrid = 1 - hm - complex;
end
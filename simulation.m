function [hm,complex,hybrid] = simulation(b,v_m,k_m,concentration,num_sim,num_adj,n_e)
    
change_compartment_rate = 1 / 300;
div = 0.95;
complex = 0;
hm = 0;

% set up initial concentrations (values prportional to conentration - VPTC) of glycans, available to each of enzymes

for i = 1:4    
    for j = 1:n_e
        s(j,i) = 100;
    end;    
    ss(i) = 100 * (n_e - 1);
end;

% adjust VPTC

for i = 1:num_adj
    [s,ss] = adjustment(b,s,ss,v_m,k_m,change_compartment_rate,concentration,div,n_e);
end;

% normalize VPTC

for i = 1:4
    for j = 1:n_e
        s(j,i) = s(j,i) / ss(i);
    end;    
end;

% run the model and record results on a features of the interest

for i = 1:num_sim
    a = golgi_sim_full(b,s,v_m,k_m,change_compartment_rate,concentration,n_e);
    if a(2) + a(3) + a(4)>1 && a(7) == 0
        hm = hm + 1;
    end;
    if a(2) + a(3) + a(4) == 0
        complex = complex + 1;
    end;
end;

% normalize results

hm = hm / num_sim;
complex = complex / num_sim;
hybrid = 1 - hm - complex;
end

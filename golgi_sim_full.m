function [a] = golgi_sim_full(a,s,Vm,Km,ChangeCompartmentRate,concentration,NE)

compartment = a(11);
sia = 0;

while compartment < 5
    
    [e, MM, Msum] = available_reactions(a,s,Vm,Km,ChangeCompartmentRate,concentration,NE,sia);

    [a,sia] = perform_reaction(MM,Msum,a,sia,NE);
    compartment=a(11);
end;
end
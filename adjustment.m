function [sx,ssx] = adjustment(a,sx,ssx,Vm,Km,ChangeCompartmentRate,concentration,div,NE)

compartment = a(11);
sia = 0;

for i = 1:4
    for j = 1:NE
        s(j,i) = sx(j,i) / ssx(i);
    end;
end;

while compartment < 5
    
    % Determine possible reactions and likelihood of each of them
    [e, MM, Msum] = available_reactions(a,s,Vm,Km,ChangeCompartmentRate,concentration,NE,sia);

    % Time until the next reaction
    if MM > 0
        time = -log(rand()) / MM;
    end;
    
    % Count time until the next reaction for all possible reactions
    if e(1) == 1
        sx(1,compartment) = sx(1,compartment) + time;
        sx(2,compartment) = sx(2,compartment) + time;
        sx(3,compartment) = sx(3,compartment) + time;
    end;
    
    if e(2) == 1
        sx(1,compartment) = sx(1,compartment) + time;
        sx(2,compartment) = sx(2,compartment) + time;
        sx(3,compartment) = sx(3,compartment) + time;
    end;
    
    if e(3) == 1
        sx(1,compartment) = sx(1,compartment) + time;
        sx(2,compartment) = sx(2,compartment) + time;
        sx(3,compartment) = sx(3,compartment) + time;
    end;
    
    for j = 4:9
        if e(j) == 1
            sx(j,compartment) = sx(j,compartment) + time;
        end;
    end;
    
    if e(14) == 1
        sx(14,compartment) = sx(14,compartment) + time;
    end;            
    
    if e(10) == 1
        sx(10,compartment) = sx(10,compartment) + time;
        sx(11,compartment) = sx(11,compartment) + time;
        sx(12,compartment) = sx(12,compartment) + time;
        sx(13,compartment) = sx(13,compartment) + time;
    end;                        
    if e(11) == 1
        sx(10,compartment) = sx(10,compartment) + time;
        sx(11,compartment) = sx(11,compartment) + time;
        sx(12,compartment) = sx(12,compartment) + time;
        sx(13,compartment) = sx(13,compartment) + time;
    end;                        
    if e(12) == 1
        sx(10,compartment) = sx(10,compartment) + time;
        sx(11,compartment) = sx(11,compartment) + time;
        sx(12,compartment) = sx(12,compartment) + time;
        sx(13,compartment) = sx(13,compartment) + time;
    end;                        
    if e(13) == 1
        sx(10,compartment) = sx(10,compartment) + time;
        sx(11,compartment) = sx(11,compartment) + time;
        sx(12,compartment) = sx(12,compartment) + time;
        sx(13,compartment) = sx(13,compartment) + time;
    end;
    
    if e(15) == 1
        sx(15,compartment) = sx(10,compartment) + time;
        sx(16,compartment) = sx(11,compartment) + time;
        sx(17,compartment) = sx(12,compartment) + time;
        sx(18,compartment) = sx(13,compartment) + time;
    end;
    if e(16) == 1
        sx(15,compartment) = sx(10,compartment) + time;
        sx(16,compartment) = sx(11,compartment) + time;
        sx(17,compartment) = sx(12,compartment) + time;
        sx(18,compartment) = sx(13,compartment) + time;
    end;
    if e(17) == 1
        sx(15,compartment) = sx(10,compartment) + time;
        sx(16,compartment) = sx(11,compartment) + time;
        sx(17,compartment) = sx(12,compartment) + time;
        sx(18,compartment) = sx(13,compartment) + time;
    end;
    if e(18) == 1
        sx(15,compartment) = sx(10,compartment) + time;
        sx(16,compartment) = sx(11,compartment) + time;
        sx(17,compartment) = sx(12,compartment) + time;
        sx(18,compartment) = sx(13,compartment) + time;
    end;
        
    if e(19) == 1
        sx(19,compartment) = sx(10,compartment) + time;
        sx(20,compartment) = sx(11,compartment) + time;
        sx(21,compartment) = sx(12,compartment) + time;
        sx(22,compartment) = sx(13,compartment) + time;
    end;
    if e(20) == 1
        sx(19,compartment) = sx(10,compartment) + time;
        sx(20,compartment) = sx(11,compartment) + time;
        sx(21,compartment) = sx(12,compartment) + time;
        sx(22,compartment) = sx(13,compartment) + time;
    end;
    if e(21) == 1
        sx(19,compartment) = sx(10,compartment) + time;
        sx(20,compartment) = sx(11,compartment) + time;
        sx(21,compartment) = sx(12,compartment) + time;
        sx(22,compartment) = sx(13,compartment) + time;
    end;
    if e(22) == 1
        sx(19,compartment) = sx(10,compartment) + time;
        sx(20,compartment) = sx(11,compartment) + time;
        sx(21,compartment) = sx(12,compartment) + time;
        sx(22,compartment) = sx(13,compartment) + time;
    end;

    ssx(compartment) = ssx(compartment) + time;

    [a,sia] = perform_reaction(MM,Msum,a,sia,NE);
    compartment = a(11);
end;

for i = 1:4                    
    for j = 1:NE
        sx(j,i) = sx(j,i) * div;                        
    end;
    ssx(i) = ssx(i) * div;
end;
end
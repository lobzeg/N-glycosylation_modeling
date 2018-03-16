function [a] = golgi_sim_full(a,s,v_m,k_m,change_compartment_rate,concentration,n_e)

compartment = a(11);
sia = 0;

% Perform reactions until reaction "leave 4th ompartment" happens
while compartment < 5
    
    % Determine possible reactions and likelihood of each of them
    [e, m_m, m_sum] = available_reactions(a,s,v_m,k_m,change_compartment_rate,concentration,n_e,sia);

    % Perform a random reaction from available ones
    [a,sia] = perform_reaction(m_m,m_sum,a,sia,n_e);
    compartment=a(11);
end;
end

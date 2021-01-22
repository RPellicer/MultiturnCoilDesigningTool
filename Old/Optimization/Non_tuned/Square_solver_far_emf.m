function s = Square_solver_far_emf(x)
    Nwires = x(1);
    percentwire = x(2);
    
    s = Square_coil_struc(Nwires,percentwire);
end
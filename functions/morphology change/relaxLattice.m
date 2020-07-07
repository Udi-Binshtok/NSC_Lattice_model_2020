function g = relaxLattice(g,n,resolution)
% relaxLattice will change the morphology of the lattice by minimizing the mechanical energy of the lattice
%   g is the geometric structure of the lattice
%   n is the number of times the function runs
    for i = 1:n
        ve = extractverts(g);
        dE = denergy(ve,g);
        noE = norm(dE);
        if noE > 1E-10
% %             dE = 0.001*dE/noE;
% %             dE = 0.01*dE/noE;
% %             dE = 0.02*dE/noE;
            dE = resolution*dE/noE;
            ve = ve - dE';
            g = insertverts(ve,g);
        end
    end

end
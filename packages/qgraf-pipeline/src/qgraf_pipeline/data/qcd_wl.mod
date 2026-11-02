% QCD Wilson Loop model - quadrilateral (4 distinct sides)
%
% 4 worldline fermions z1-z4 for 4 sides of the Wilson loop
% 4 auxiliary scalars phi1-phi4 for corner vertices (external only)
%
% Propagators
[A, A, +]
[c, cbar, -]
[q, qbar, -]
[z1, z1bar, -]
[z2, z2bar, -]
[z3, z3bar, -]
[z4, z4bar, -]
[phi1, phi1, +]
[phi2, phi2, +]
[phi3, phi3, +]
[phi4, phi4, +]
%
% QCD vertices
[A, A, A]
[A, A, A, A]
[cbar, A, c]
[qbar, A, q]
%
% Wilson line - gluon coupling (one per side)
[z1bar, A, z1]
[z2bar, A, z2]
[z3bar, A, z3]
[z4bar, A, z4]
%
% Corner vertices: phi_i connects side i to side i+1
% Loop direction: z1 -> z2 -> z3 -> z4 -> z1
[z1bar, phi1, z2]
[z2bar, phi2, z3]
[z3bar, phi3, z4]
[z4bar, phi4, z1]

% This subscript is just for examining the feasibility of M matrix.
% 2022.12.28
% Since the two SRAs are symmetric about the central line, we just need to
% examine one side.
clc;
n_eig=length(M.L);
eg_val=struct('L',zeros(n_eig,1));
eg_val.L=simplify(eig(M.L));
det_L=simplify(det(M.L));
trace_L=simplify(trace(M.L));
function H=H_1Ds_inter(tneg,SO)
%hopping form ~ +tneg*a_id*a_j; basis={pu,pd,hu,hd}=(ph)x(spin)x(x); SO_x~(\sum\sig) \cdot p_y
%used for constructing ribbon's inter-slab H1
%or the (upper) sub-diagonal blocks of the tridiagonal block structure of the Hamiltonian of quasi-1D system 
N=max([length(tneg) size(SO,2)]);
if length(tneg)==1, tneg=tneg*ones(1,N); end
if size(SO,2)==1,   SO=SO*ones(1,N); end
tneg=tneg(:)'; 
Htu=spdiags(tneg',0,N,N);
Hsoxy=0.5*spdiags(SO(2,:)'+i*SO(1,:)',0,N,N);
Hsozu=0.5*spdiags(i*SO(3,:)',0,N,N);
Hph=sparse(2*N,2*N);

Hp=  [Htu+Hsozu Hsoxy; -Hsoxy' Htu-Hsozu]; %always anti-Hermitian
H=[Hp,  Hph;
   -Hph', -Hp];
end
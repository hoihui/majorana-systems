function H=H_1Ds_LR(tneg,SO,B,DS,DT)
% LR=long-range
%hopping form ~ +tneg*a_id*a_j; basis={pu,pd,hu,hd}=(ph)x(spin)x(sites); SO_y~(\sum\sig) \cdot p_y
%shift_on_site_energy: due to discretization (mu->mu+2*tneg); i.e. measures mu from band bottom
%INPUTS:
% tneg: (n+1)x1 negative chemical potential / hopping
% SOx: nx3 coefficients of PAULI MATRICES coupled to the momentum along wire
%      =\alpha/(a)=sqrt(Eso*t) where Eso=2m\alpha^2 (arXiv:1206.1276)
%     Hence a scaling a->ba, k->k/b is mitigated by \alpha->b\alpha, m->m/sqrt(b)
%     i.e. band shape is the same (up to a scaling) so long as m\alpha^2 is kept the same
% B: 3xN magnetic field acting on PAULI MATRICES
% DeltaS: (n+1)x1 Singlet pairing
% DeltaT: nx1 Triplet pairing
%at least 1 parameter must reflect the size of the system
%OUTPUTS:
% sparse H whose u's and v's at x for eigen mode n is given by:
%u_u = EVEC(x,n) for x in [1,...,N]
%u_d = EVEC(N+x,n)
%v_u = EVEC(2*N+x,n)
%v_d = EVEC(3*N+x,n)
N=size(B,2);
Htu=spdiags(ones(N,1)*tneg(1),0,N,N);
for ii=2:size(tneg,1)
    Htu=Htu+spdiags(ones(N,1)*[tneg(ii) tneg(ii)],[-ii+1 ii-1],N,N);
end
% disp(full(Htu))
Hsoxy=spdiags([],[],N,N);
Hsozu=spdiags([],[],N,N);
for ii=1:size(SO,1)
    Hsoxy=Hsoxy+spdiags(ones(N,1)*[-SO(ii,2)-1i*SO(ii,1),SO(ii,2)+1i*SO(ii,1)],[-ii ii],N,N);
    Hsozu=Hsozu+spdiags(ones(N,1)*[-1i*SO(ii,3),1i*SO(ii,3)],[-ii ii],N,N);
end
Hphud=spdiags(ones(N,1)*DS(1),0,N,N);
for ii=2:size(DS,1)
    Hphud=Hphud+spdiags(ones(N,1)*[DS(ii),DS(ii)],[-ii+1 ii-1],N,N);
end
Hphxy=spdiags([],[],N,N);
Hphzu=spdiags([],[],N,N);
for ii=1:size(DT,1)
    Hphxy=Hphxy+spdiags(ones(N,1)*[-DT(ii,2)-1i*DT(ii,1),DT(ii,2)+1i*DT(ii,1)],[-ii ii],N,N);
    Hphzu=Hphzu+spdiags(ones(N,1)*[-1i*DT(ii,3),1i*DT(ii,3)],[-ii ii],N,N);
end

Hzu=spdiags(B(3,:).',0,N,N);    Hzud=spdiags(B(1,:).'-1i*B(2,:).',0,N,N);
Hz=[Hzu Hzud ; Hzud' -Hzu];
Hp=  [Htu+Hsozu Hsoxy; Hsoxy' Htu-Hsozu] + Hz; %particle
Hh= -[Htu+Hsozu Hsoxy; Hsoxy' Htu-Hsozu] + Hz; %hole
Hph= [Hphud+Hphzu Hphxy; Hphxy' Hphud-Hphzu];
H=[Hp,  Hph;
   Hph', Hh];

end
% H=sparse(H);
% Particle-Hole checked by:
%  sigy=[0 -i;i 0]; PH=kron(sparse(kron(sigy,sigy)),speye(length(H)/4)); PH*conj(H)+H*PH
% Time-Reversal, if exists, is checked by:
%  sigy=[0 -i;i 0]; TR=kron(sparse(kron(speye(2),i*sigy)),speye(length(H)/4)); TR*conj(H)-H*TR
% Transform to Majorana basis (c_Au c_Ad c_Bu c_Bd)x(sites) by:
%  MB=kron(sparse([1 0 0 -1;0 1 1 0;i 0 0 i;0 i -i 0]/sqrt(2)),speye(length(H)/4)); HA=MB*PH*MB'/(i);
% Transform to (sites)x(ph)x(spin) , useful for scattering matrix comp:
%  T=spalloc(4,4*N,4);T(1,N)=1;T(2,2*N)=1;T(3,3*N)=1;T(4,4*N)=1;U=T;for ii=1:N-1,U=[circshift(T,[0 -ii]);U];end; H=U*H*U';
%  Further, to (sites)x(Majorana):
%  U=kron(speye(length(H)/4),sparse([1 0 0 -1;0 1 1 0;i 0 0 i;0 i -i 0]/sqrt(2));end; HA=U*H*U'/(i/4);

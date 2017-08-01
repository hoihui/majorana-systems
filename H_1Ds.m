function H=H_1Ds(tneg,mu,SO,B,DeltaS,bc,shift_on_site_energy)
%hopping form ~ +tneg*a_id*a_j; basis={pu,pd,hu,hd}=(ph)x(spin)x(sites); SO_y~(\sum\sig) \cdot p_y
%shift_on_site_energy: due to discretization (mu->mu+2*tneg); i.e. measures mu from band bottom
%INPUTS:
% tneg: 1x1 or 1xN negative hopping matrix element=hbar^2/(2ma^2)
%       2nd row:NNN hopping
% mu: 1x1 or 1xN chemical potential
% SOx: 3x1 or 3xN coefficients of PAULI MATRICES coupled to the momentum along wire
%      =\alpha/(a)=sqrt(Eso*t) where Eso=2m\alpha^2 (arXiv:1206.1276)
%     Hence a scaling a->ba, k->k/b is mitigated by \alpha->b\alpha, m->m/sqrt(b)
%     i.e. band shape is the same (up to a scaling) so long as m\alpha^2 is kept the same
% B: 3x1 or 3xN magnetic field acting on PAULI MATRICES
% DeltaS: ?x1 or ?xN Singlet pairing
%         1st row: on-site s-wave; 2nd row: d-wave..
% shift_on_site_energy: whether to measure from band bottom
% lambda: extra 3x1 or 3xN spin-dependent hopping (.5*sigma*p^2) to gap out modes at k=pi for TI
%         always shifted such that 1-cos[k] 
%at least 1 parameter must reflect the size of the system
%OUTPUTS:
% sparse H whose u's and v's at x for eigen mode n is given by:
%u_u = EVEC(x,n) for x in [1,...,N]
%u_d = EVEC(N+x,n)
%v_u = EVEC(2*N+x,n)
%v_d = EVEC(3*N+x,n)
N=max([length(tneg) length(mu) size(SO,2) size(B,2) length(DeltaS)]);
if size(tneg,2)==1, tneg=tneg*ones(1,N); end
if length(mu)==1,   mu=mu*ones(1,N); end
if size(SO,2)==1,   SO=SO*ones(1,N); end
if size(B,2)==1, B=B*ones(1,N); end
if size(DeltaS,2)==1, DeltaS=DeltaS*ones(1,N); end
if nargin<6, bc=1; end
if nargin<7, shift_on_site_energy=0; end
% if nargin<8, lambda=zeros(3,N); end    %still also use nargin>=8 criterion below
if shift_on_site_energy==1,mu=mu + (tneg+circshift(tneg,[1 1]));end
if length(bc)==1, bc=bc*ones(1,4);
elseif length(bc)==2, bc=kron(bc(:).',ones(1,2)); end
% Htu=diag(tneg(1:N-1),1)+diag(tneg(1:N-1),-1)-diag(mu);
% Hsoxy=0.5*(i*diag(SO(1,1:N-1),1)-i*diag(SO(1,1:N-1),-1)+... %hopping d to u (opposite signs of going backward (1st term) / forward (2nd term))
%              diag(SO(2,1:N-1),1)-  diag(SO(2,1:N-1),-1));   % the u to d part is automatically handled by its h.c.
% Hsozu=0.5*(i*diag(SO(3,1:N-1),1)-i*diag(SO(3,1:N-1),-1));   %hopping u to u (opposite signs of going backward (1st term) / forward (2nd term))
% Hzu=diag(Vzee(3,:));    Hzud=diag(Vzee(1,:)-i*Vzee(2,:));
% Hph=diag([DeltaS DeltaS]);
tneg=tneg.'; DeltaS=DeltaS.';
Htu=spdiags([tneg(:,1) -mu' circshift(tneg(:,1),1)],[-1 0 1],N,N);
for ii=2:size(tneg,2)
    Htu=Htu+spdiags([tneg(:,ii) circshift(tneg(:,ii),ii)],[-ii ii],N,N);
end
Hsoxy=0.5*spdiags([-SO(2,:).'-1i*SO(1,:).',circshift(SO(2,:).'+1i*SO(1,:).',1)],[-1 1],N,N); %factor of i needed to make sure Hermiticity, factor of 0.5 due to finite difference with 2 steps
Hsozu=0.5*spdiags([-1i*SO(3,:).',circshift(1i*SO(3,:).',1)],[-1 1],N,N);
Hzu=spdiags(B(3,:).',0,N,N);    Hzud=spdiags(B(1,:).'-1i*B(2,:).',0,N,N);
Hph=spdiags([DeltaS(:,1);DeltaS(:,1)],0,2*N,2*N);
for ii=2:size(DeltaS,2)
    Hphud=0.5*spdiags([DeltaS(:,ii) circshift(DeltaS(:,ii),ii-1)],[-ii+1 ii-1],N,N);
    Hph=Hph+kron(sparse([1 0;0 1]),Hphud);
end
% if nargin>=8
%     if size(lambda,2)==1,   lambda=lambda*ones(1,N); end
%     Hzud=Hzud+0.5*spdiags([lambda(1,:).'-1i*lambda(2,:).',circshift(lambda(1,:).'-1i*lambda(2,:).',[1 1])],[-1 1],N,N); %spin-dependent p^2 to gap out k=pi
%     Hzud=Hzud-0.5*spdiags(lambda(1,:).'-1i*lambda(2,:).'+circshift(lambda(1,:).'-1i*lambda(2,:).',[1 1]),0,N,N);  % on-site shift
%     Hzu=Hzu+0.5*spdiags([lambda(3,:).',circshift(lambda(3,:).',[1 1])],[-1 1],N,N);
%     Hzu=Hzu-0.5*spdiags(lambda(3,:).'+circshift(lambda(3,:).',[1 1]),0,N,N);
% end
% if any(bc)==1
%     Htu(1,N)=bc(1)*tneg(N); Htu(N,1)=bc(1)*tneg(N);
%     Hsoxy(1,N)=bc(1)*0.5*(-SO(2,N)-1i*SO(1,N)); Hsoxy(N,1)=bc(1)*0.5*(SO(2,N)+1i*SO(1,N));
%     Hsozu(1,N)=-bc(1)*1i*0.5*SO(3,N); Hsozu(N,1)=bc(1)*1i*0.5*SO(3,N);
%     Hzud(1,N)=bc(1)*(0.5*(lambda(1,N)-1i*lambda(2,N)));
%     Hzud(N,1)=bc(1)*(0.5*(lambda(1,N)-1i*lambda(2,N)));
%      Hzu(1,N)=bc(1)*(0.5*lambda(3,N));
%      Hzu(N,1)=bc(1)*(0.5*lambda(3,N));
% end
Hz=[Hzu Hzud ; Hzud' -Hzu];
Hp=  [Htu+Hsozu Hsoxy; Hsoxy' Htu-Hsozu] + Hz; %particle
Hh= -[Htu+Hsozu Hsoxy; Hsoxy' Htu-Hsozu] + Hz; %hole
H=[Hp,  Hph;
   Hph', Hh];
if any(bc)==1
    %TODO:DeltaS(:,2), tneg(2,:)
    H(    1,  N)= bc(1)*tneg(N); H(  N,    1)= bc(1)'*tneg(N);   %puN hops to pu1
    H(2*N+1,3*N)=-bc(3)*tneg(N); H(3*N,2*N+1)=-bc(3)'*tneg(N);   %huN hops to hu1
    H(  N+1,2*N)= bc(2)*tneg(N); H(2*N,  N+1)= bc(2)'*tneg(N);   %pdN hops to pd1
    H(3*N+1,4*N)=-bc(4)*tneg(N); H(4*N,3*N+1)=-bc(4)'*tneg(N);   %hdN hops to hd1
    H(    1,  N)=H(    1,  N)-bc(1)*1i*0.5*SO(3,N); H(  N,    1)=H(  N,    1)+bc(1)'*1i*0.5*SO(3,N);  %puN SOz to pu1
    H(2*N+1,3*N)=H(2*N+1,3*N)+bc(3)*1i*0.5*SO(3,N); H(3*N,2*N+1)=H(3*N,2*N+1)-bc(3)'*1i*0.5*SO(3,N);  %huN SOz to hu1
    H(  N+1,2*N)=H(  N+1,2*N)+bc(2)*1i*0.5*SO(3,N); H(2*N,  N+1)=H(2*N,  N+1)-bc(2)'*1i*0.5*SO(3,N);  %pdN SOz to pd1
    H(3*N+1,4*N)=H(3*N+1,4*N)-bc(4)*1i*0.5*SO(3,N); H(4*N,3*N+1)=H(4*N,3*N+1)+bc(4)'*1i*0.5*SO(3,N);  %hdN SOz to hd1
    H(    1,2*N)=H(    1,2*N)+bc(1)*0.5*(-SO(2,N)-1i*SO(1,N)); H(2*N,    1)=H(2*N,    1)+bc(1)'*0.5*(-SO(2,N)+1i*SO(1,N));  %pdN SOxy to pu1
    H(2*N+1,4*N)=H(2*N+1,4*N)-bc(3)*0.5*(-SO(2,N)-1i*SO(1,N)); H(4*N,2*N+1)=H(4*N,2*N+1)-bc(3)'*0.5*(-SO(2,N)+1i*SO(1,N));  %hdN SOxy to hu1
    H(  N+1,  N)=H(  N+1,  N)+bc(2)*0.5*( SO(2,N)-1i*SO(1,N)); H(  N,  N+1)=H(  N,  N+1)+bc(2)'*0.5*( SO(2,N)+1i*SO(1,N));  %puN SOxy to pd1
    H(3*N+1,3*N)=H(3*N+1,3*N)-bc(4)*0.5*( SO(2,N)-1i*SO(1,N)); H(3*N,3*N+1)=H(3*N,3*N+1)-bc(4)'*0.5*( SO(2,N)+1i*SO(1,N));  %huN SOxy to hd1
%     H(    1,  N)=H(    1,  N)+bc(1)*(0.5*lambda(3,N));    H(  N,    1)=H(  N,    1)+bc(1)'*(0.5*lambda(3,N));  %puN Lambdaz to pu1
%     H(2*N+1,3*N)=H(2*N+1,3*N)+bc(3)*(0.5*lambda(3,N));    H(3*N,2*N+1)=H(3*N,2*N+1)+bc(3)'*(0.5*lambda(3,N));  %huN Lambdaz to hu1
%     H(  N+1,2*N)=H(  N+1,2*N)-bc(2)*(0.5*lambda(3,N));    H(2*N,  N+1)=H(2*N,  N+1)-bc(2)'*(0.5*lambda(3,N));  %pdN Lambdaz to pd1
%     H(3*N+1,4*N)=H(3*N+1,4*N)-bc(4)*(0.5*lambda(3,N));    H(4*N,3*N+1)=H(4*N,3*N+1)-bc(4)'*(0.5*lambda(3,N));  %hdN Lambdaz to hd1
%     H(    1,2*N)=H(    1,2*N)+bc(1)*(0.5*(lambda(1,N)-1i*lambda(2,N))); H(2*N,    1)=H(2*N,    1)+bc(1)'*(0.5*(lambda(1,N)+1i*lambda(2,N)));  %pdN Lambdaxy to pu1
%     H(2*N+1,4*N)=H(2*N+1,4*N)+bc(3)*(0.5*(lambda(1,N)-1i*lambda(2,N))); H(4*N,2*N+1)=H(4*N,2*N+1)+bc(3)'*(0.5*(lambda(1,N)+1i*lambda(2,N)));  %hdN Lambdaxy to hu1
%     H(  N+1,  N)=H(  N+1,  N)+bc(2)*(0.5*(lambda(1,N)+1i*lambda(2,N))); H(  N,  N+1)=H(  N,  N+1)+bc(2)'*(0.5*(lambda(1,N)-1i*lambda(2,N)));  %puN Lambdaxy to pd1
%     H(3*N+1,3*N)=H(3*N+1,3*N)+bc(4)*(0.5*(lambda(1,N)+1i*lambda(2,N))); H(3*N,3*N+1)=H(3*N,3*N+1)+bc(4)'*(0.5*(lambda(1,N)-1i*lambda(2,N)));  %huN Lambdaxy to hd1
end

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
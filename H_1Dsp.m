function H=H_1Dsp(tneg,mu,SO,Vzee,DeltaS,DeltaP,pbc,shift_on_site_energy)
%hopping form ~ +tneg*a_id*a_j; basis={pu,pd,hu,hd}
%SO_y~(\sum(rows)\sig) \cdot p_y; extract terms ~ py in p-wave pairing
%shift_on_site_energy: due to discretization (mu->mu+tneg)
N=length(tneg);
if nargin<7, pbc=1; end
if nargin<8, shift_on_site_energy=0; end
if shift_on_site_energy==1
    mu=mu + (tneg+circshift(tneg,[1 1]));
end
% Htu=diag(tneg(1:N-1),1)+diag(tneg(1:N-1),-1)-diag(mu);
% Hsoxy=0.5*(i*diag(SO(1,1:N-1),1)-i*diag(SO(1,1:N-1),-1)+... %hopping d to u (opposite signs of going backward (1st term) / forward (2nd term))
%              diag(SO(2,1:N-1),1)-  diag(SO(2,1:N-1),-1));   % the u to d part is automatically handled by its h.c.
% Hsozu=0.5*(i*diag(SO(3,1:N-1),1)-i*diag(SO(3,1:N-1),-1));   %hopping u to u (opposite signs of going backward (1st term) / forward (2nd term))
% Hzu=diag(Vzee(3,:));    Hzud=diag(Vzee(1,:)-i*Vzee(2,:));
% Hph=diag([DeltaS DeltaS]);
tneg=tneg(:); mu=mu(:); DeltaS=DeltaS(:);  %sparse construction
Htu=spdiags([tneg -mu circshift(tneg,1)],[-1 0 1],N,N);
Hsoxy=0.5*spdiags([-SO(2,:)'-i*SO(1,:)',circshift(SO(2,:)'+i*SO(1,:)',1)],[-1 1],N,N);
Hsozu=0.5*spdiags([-i*SO(3,:)',circshift(i*SO(3,:)',1)],[-1 1],N,N);
Hzu=spdiags(Vzee(3,:)',0,N,N);    Hzud=spdiags(Vzee(1,:)'-i*Vzee(2,:)',0,N,N);
Hphud=0.5*spdiags([DeltaP(3,:)',2*DeltaS,-circshift(DeltaP(3,:)',1)],[-1 0 1],N,N); %upper-left block of pairing block
Hphdu=0.5*spdiags([-DeltaP(3,:)',2*DeltaS,circshift(DeltaP(3,:)',1)],[-1 0 1],N,N); %lower-right block of pairing block
Hphdd=0.5*spdiags([DeltaP(2,:)'-i*DeltaP(1,:)',-circshift(DeltaP(2,:)'-i*DeltaP(1,:)',1)],[-1 1],N,N); %upper-right block of pairing block %dangerous
Hphuu=0.5*spdiags([-DeltaP(2,:)'-i*DeltaP(1,:)',-circshift(-DeltaP(2,:)'-i*DeltaP(1,:)',1)],[-1 1],N,N); %lower-left block of pairing block %dangerous

if pbc==1
    Htu(1,N)=tneg(N); Htu(N,1)=tneg(N);
    Hsoxy(1,N)=-0.5*(i*SO(1,N)+SO(2,N)); Hsoxy(N,1)=0.5*(i*SO(1,N)+SO(2,N));
    Hsozu(1,N)=-i*0.5*SO(3,N); Hsozu(N,1)=i*0.5*SO(3,N);
    Hphud(1,N)=0.5*DeltaP(3,N); Hphud(N,1)=-0.5*DeltaP(3,N);
    Hphdu(1,N)=0.5*DeltaP(3,N); Hphdu(N,1)=-0.5*DeltaP(3,N);
    Hphdd(1,N)= 0.5*(DeltaP(2,N)-i*DeltaP(1,N));  Hphdu(N,1)=-0.5*(DeltaP(2,N)-i*DeltaP(1,N));
    Hphuu(1,N)= 0.5*(-DeltaP(2,N)-i*DeltaP(1,N)); Hphdu(N,1)=-0.5*(-DeltaP(2,N)-i*DeltaP(1,N));
end
Hph=[Hphud Hphuu; Hphdd Hphdu];
Hz=[Hzu Hzud ;    Hzud' -Hzu];
Hp=  [Htu+Hsozu Hsoxy; Hsoxy' Htu-Hsozu] + Hz; %particle
Hh= -[Htu+Hsozu Hsoxy; Hsoxy' Htu-Hsozu] + Hz; %hole
H=[Hp,  Hph;
   Hph', Hh];
% H=sparse(H);
end
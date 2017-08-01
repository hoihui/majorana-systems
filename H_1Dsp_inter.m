function H=H_1Dsp_inter(tneg,SO,DeltaP)
%hopping form ~ +tneg*a_id*a_j; basis={pu,pd,hu,hd}
%used for constructing ribbon's inter-slab H1
N=length(tneg);
tneg=tneg(:); 
Htu=spdiags(tneg,0,N,N);
Hsoxy=0.5*spdiags(SO(2,:)'+i*SO(1,:)',0,N,N);
Hsozu=0.5*spdiags(i*SO(3,:)',0,N,N);
Hphud=0.5*spdiags(DeltaP(3,:)',0,N,N); %upper-left block of pairing block
Hphdu=0.5*spdiags(DeltaP(3,:)',0,N,N); %lower-right block of pairing block
Hphdd=0.5*spdiags(-DeltaP(2,:)'+i*DeltaP(1,:)',0,N,N); %upper-right block of pairing block %dangerous
Hphuu=0.5*spdiags( DeltaP(2,:)'+i*DeltaP(1,:)',0,N,N); %lower-left block of pairing block %dangerous

Hph=[Hphud Hphuu; Hphdd Hphdu];

Hp=  [Htu+Hsozu Hsoxy; -Hsoxy' Htu-Hsozu]; %always anti-Hermitian
H=[Hp,  Hph;
   -Hph', -Hp];
end
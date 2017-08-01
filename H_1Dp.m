function [H0,H1]=H_1Dp(tnegx,tnegy,mu,Dpx,Dpy,bc)
N=max([length(tnegx) length(tnegy) length(mu) length(Dpx) length(Dpy)]);
if length(tnegx(:))==1, tnegx=tnegx*ones(1,N); end
if length(tnegy(:))==1, tnegy=tnegy*ones(N,1); end
if length(mu)==1, mu=mu*ones(1,N); end
if length(Dpx(:))==1, Dpx=Dpx*ones(1,N); end
if length(Dpy(:))==1, Dpy=Dpy*ones(N,1); end
if nargin<6, bc=1; end
if length(bc)==1, bc=bc*ones(1,2);end
tnegx=tnegx.';tnegy=tnegy(:); mu=mu(:); Dpx=Dpx.'; Dpy=Dpy(:);  %sparse construction
H0p=spdiags([tnegx(:,1) -mu circshift(tnegx(:,1),1)],-1:1,N,N);
for ii=2:size(tnegx,2)
    H0p=H0p+spdiags([tnegx(:,ii) circshift(tnegx(:,ii),ii)],[-ii ii],N,N);
end
H0ph=spdiags([-Dpx(:,1),circshift(Dpx(:,1),1)],[-1 1],N,N); %hopping forward/backward always have -ve signs
for ii=2:size(Dpx,2)
    H0ph=H0ph+spdiags([-Dpx(:,ii) circshift(Dpx(:,ii),ii)],[-ii ii],N,N);
end
H0=[H0p,H0ph;H0ph',-H0p];
if any(bc)==1 %TODO for NNN and so on
    H0(  1,  N)= bc(1)*tnegx(N);           H0(  N,  1)= bc(1)'*tnegx(N);
    H0(N+1,2*N)=-bc(2)*tnegx(N);           H0(2*N,N+1)=-bc(2)'*tnegx(N); 
    H0(  1,2*N)=-sqrt(bc(1)*bc(2))*Dpx(N); H0(2*N,  1)=-sqrt(bc(1)*bc(2))*Dpx(N)'; 
    H0(  N,N+1)= sqrt(bc(1)*bc(2))*Dpx(N); H0(N+1,  N)= sqrt(bc(1)*bc(2))*Dpx(N)';
end
if nargout==2
    H1p=spdiags(tnegy,0,N,N);
    H1ph=spdiags(Dpy,0,N,N);
    H1=[H1p,H1ph;-H1ph',-H1p];
    %Considering H between x and (x+1) site: Cdag_k(x+1)*H1*C_kx, where C_k={c_k,cdag_-k} 
    %For each k, attach exp(-ik) to C_kx.
    %The upper right term is iDp*cdag_k(x+1)*cdag_-k(x) * exp(-ik)    (A)
    %whose H.c. is -iDp*c_-k(x)*c_k(x+1)*exp(ik) = iDp*c_k(x+1)*c_-k(x)*exp(ik)
    %To assemble terms of Cdag_k(x+1), need to consider -k components together
    % resulting in a term  iDp*c_-k(x+1)*c_k(x) * exp(-ik) whose coefficient is
    % identical to A, hence the form of diagonal parts of H1 above
end

end
% Contructing full H for 2D:
% Hod=kron(spdiags(ones(Ny,1),1,Ny,Ny),H1); H=kron(speye(L),H0)+(Hod+Hod');
% Particle-Hole checked by:
%  PH=kron([0 1;1 0],speye(length(H0)/2)); PH*conj(H0)+H*PH
% Transform to Majorana basis (c_Au c_Ad c_Bu c_Bd)x(sites) by:
%  MB=kron(sparse([1 1;i -i]),speye(length(H)/2)); HA=MB*H*MB'/(i/4);
% Transform to (sites)x(ph)x(spin) , useful for scattering matrix comp: TODO
%  T=spalloc(4,4*N,4);T(1,N)=1;T(2,2*N)=1;T(3,3*N)=1;T(4,4*N)=1;U=T;for ii=1:N-1,U=[circshift(T,[0 -ii]);U];end; H=U*H*U';
%  Further, to (sites)x(Majorana):
%  U=kron(speye(length(H)/4),sparse([1 0 0 -1;0 1 1 0;i 0 0 i;0 i -i 0]/sqrt(2));end; HA=U*H*U'/(i/4);
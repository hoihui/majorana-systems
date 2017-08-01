function [H,f]=H_1Ds_Fourier_TI(v,mu,B,DeltaS,nband,lambda)
% basis=(ph)x(spin)x(k); use the same Fourier basis for p & h st. uniform Delta -> diagonal 
%INPUTS:
% mu: 1xN or 1x1 chemical potential 
% SOx: 3x1 coefficients of spin matrices coupled to the momentum along wire
% B: 3xN or 3x1 magnetic field
% DeltaS: 1xN or 1x1 pairing
% lambda: spin-dependent hopping to gap out band edge for TI system (to compare with old codes)
%at least 1 parameter must reflect the size of the system
%OUTPUTS:
% sparse H whose u's and v's at x for eigen mode n is given by:
%f=f(x,:);
%u_u = f*EVEC(1:nbands,n)
%u_d = f*EVEC(nbands+1:2*nbands,n)
%v_d = f*EVEC(2*nbands+1:3*nbands,n)
%-v_u = f*EVEC(3*nbands+1:4*nbands,n)
N=max([length(mu) size(v,2) size(B,2) length(DeltaS)]);
vm=mean(v,2);vfluct=v-vm*ones(1,size(v,2));
if nargin<5||nband==0, nband=N; end
if length(mu)==1,   mu=mu*ones(1,N); end
if size(B,2)==1, B=B*ones(1,N); end
if length(DeltaS)==1, DeltaS=DeltaS*ones(1,N); end
nband=min(nband,N); nband2=ceil((nband-1)/2);
mu=mu(:).'; DeltaS=DeltaS(:).';
p1=sparse([0 1;1 0]);p2=sparse([0 -1i;1i 0]);p3=sparse([1 0;0 -1]);
nlist=round((-N+1.5)/2):round((N-0.5)/2);
if nband<N
    if nargin>=6
        [~,ind]=sort(abs(0.5*sin(2*nlist*pi/N)*sqrt(sum(sum(v.^2)))-mean(mu)));
    else
        [~,ind]=sort(abs(0.5*(2*nlist*pi/N)*sqrt(sum(sum(v.^2)))-mean(mu)));
    end
    nlist=nlist(ind);    nlist=nlist(1:nband);%E(nbands) %current cutoff
end
nMat=repmat(nlist.',[1 N]); xMat=repmat(1:N,[nband 1]);
f=sqrt(1/N)*exp(1i*2*pi*nMat.*xMat./N);    conjf=conj(f); %f(bands,x)
if nargin>=6,factor=sin(2*nlist(:)*pi/N); else factor=2*nlist(:)*pi/N;end
Hv=spdiags([(vm(1)+i*vm(2))*factor  vm(3)*factor (vm(1)-i*vm(2))*factor;
            (vm(1)+i*vm(2))*factor -vm(3)*factor (vm(1)-i*vm(2))*factor],[-nband 0 nband],2*nband,2*nband);
if any(vfluct(1,:))
    pMat=ones(nband,1)*vfluct(1,:); 
    Hdiag=zeros(nband,2*nband-1);
    Hdiag(:,nband)=-0.5*i*(sum(circshift(pMat,[0 1]).*conjf.*circshift(f,[0 1]),2)- ...
                           sum(circshift(pMat,[0 0]).*conjf.*circshift(f,[0 -1]),2));
    for shift=1:nband2
        tem=-0.5*i*(sum(circshift(pMat,[0 1]).*conjf.*circshift(f,[shift 1]),2)- ...
            sum(circshift(pMat,[0 0]).*conjf.*circshift(f,[shift -1]),2));
        Hdiag(:,nband-shift)=circshift(tem,-shift);
        Hdiag(:,nband+shift)=conj(tem);
        Hdiag(:,shift)=conj(tem);
        Hdiag(:,end+1-shift)=circshift((tem),-shift);
    end
    Hv=Hv+kron([0 1;1 0],spdiags(Hdiag,-nband+1:nband-1,nband,nband));
end
if any(vfluct(2,:))
    pMat=ones(nband,1)*vfluct(2,:); 
    Hdiag=zeros(nband,2*nband-1);
    Hdiag(:,nband)=-0.5i*(sum(circshift(pMat,[0 1]).*conjf.*circshift(f,[0 1]),2)- ...
                           sum(circshift(pMat,[0 0]).*conjf.*circshift(f,[0 -1]),2));
    for shift=1:nband2
        tem=-0.5i*(sum(circshift(pMat,[0 1]).*conjf.*circshift(f,[shift 1]),2)- ...
            sum(circshift(pMat,[0 0]).*conjf.*circshift(f,[shift -1]),2));
        Hdiag(:,nband-shift)=circshift(tem,-shift);
        Hdiag(:,nband+shift)=conj(tem);
        Hdiag(:,shift)=conj(tem);
        Hdiag(:,end+1-shift)=circshift((tem),-shift);
    end
    Hv=Hv+kron([0 -1i;1i 0],spdiags(Hdiag,-nband+1:nband-1,nband,nband));
end
if any(vfluct(3,:))
    pMat=ones(nband,1)*vfluct(3,:); 
    Hdiag=zeros(nband,2*nband-1);
    Hdiag(:,nband)=-0.5*i*(sum(circshift(pMat,[0 1]).*conjf.*circshift(f,[0 1]),2)- ...
                           sum(circshift(pMat,[0 0]).*conjf.*circshift(f,[0 -1]),2));
    for shift=1:nband2
        tem=-0.5*i*(sum(circshift(pMat,[0 1]).*conjf.*circshift(f,[shift 1]),2)- ...
            sum(circshift(pMat,[0 0]).*conjf.*circshift(f,[shift -1]),2));
        Hdiag(:,nband-shift)=circshift(tem,-shift);
        Hdiag(:,nband+shift)=conj(tem);
        Hdiag(:,shift)=conj(tem);
        Hdiag(:,end+1-shift)=circshift((tem),-shift);
    end
    Hv=Hv+kron([1 0;0 -1],spdiags(Hdiag,-nband+1:nband-1,nband,nband));
end
Hv=(Hv+Hv')/2;
if length(unique(mu))==1                  %uniform mu
    Hmu=kron(speye(2),spdiags(mu.',0,nband,nband));
else
    pMat=ones(nband,1)*mu;
    Hdiag=zeros(nband,2*nband-1);
    Hdiag(:,nband)=sum(pMat.*conjf.*f,2);
    for shift=1:nband2
        elements=sum(pMat.*conjf.*circshift(f,[shift 0]),2);
        Hdiag(:,nband-shift)=circshift(elements,-shift);
        Hdiag(:,nband+shift)=conj(elements);
        Hdiag(:,shift)=conj(elements);
        Hdiag(:,end+1-shift)=circshift(elements,-shift);
    end
    Hmu=spdiags(Hdiag,-nband+1:nband-1,nband,nband);
    Hmu=kron(speye(2),(Hmu+Hmu')/2);
end
if length(unique(DeltaS))==1                 %uniform DeltaS
    HD=kron(speye(2),spdiags(DeltaS.',0,nband,nband));
else
    pMat=ones(nband,1)*DeltaS; 
    Hdiag=zeros(nband,2*nband-1);
    Hdiag(:,nband)=sum(pMat.*conjf.*f,2);
    for shift=1:nband2
        elements=sum( pMat.*conjf.*circshift(f,[shift 0]),2);
        conjelements=sum( pMat.*f.*circshift(conjf,[shift 0]),2);
        Hdiag(:,nband-shift)=circshift(elements,-shift);
        Hdiag(:,nband+shift)=conjelements;
        Hdiag(:,shift)=conjelements;
        Hdiag(:,end+1-shift)=circshift(elements,-shift);
    end
%         Hdiag(:,nband)=sum(pMat.*conjf.*conjf,2);
%         for shift=1:nband2
%             elements=sum( pMat.*conjf.*circshift(conjf,[shift 0]),2);
%             Hdiag(:,nband-shift)=circshift((elements),-shift);
%             Hdiag(:,nband+shift)=circshift((elements),0);
%             Hdiag(:,shift)=(elements);
%             Hdiag(:,end+1-shift)=circshift(elements,-shift);
%         end
    HD=kron(speye(2),spdiags(Hdiag,-nband+1:nband-1,nband,nband));
end
HZ=spdiags([],[],2*nband,2*nband);
if length(unique(B(1,:)))==1                  %uniform Bx
    HZ=HZ+kron(p1,spdiags(B(1,:).',0,nband,nband));
else
    pMat=ones(nband,1)*B(1,:);
    Hdiag=zeros(nband,2*nband-1);
    Hdiag(:,nband)=sum(pMat.*conjf.*f,2);
    for shift=1:nband2
        elements=sum(pMat.*conjf.*circshift(f,[shift 0]),2);
        Hdiag(:,nband-shift)=circshift(elements,-shift);
        Hdiag(:,nband+shift)=conj(elements);
        Hdiag(:,shift)=conj(elements);
        Hdiag(:,end+1-shift)=circshift(elements,-shift);
    end
    HZ=HZ+kron(p1,spdiags(Hdiag,-nband+1:nband-1,nband,nband));
end
if length(unique(B(2,:)))==1                  %uniform By
    HZ=HZ+kron(p2,spdiags(B(2,:).',0,nband,nband));
else
    pMat=ones(nband,1)*B(2,:);
    Hdiag=zeros(nband,2*nband-1);
    Hdiag(:,nband)=sum(pMat.*conjf.*f,2);
    for shift=1:nband2
        elements=sum(pMat.*conjf.*circshift(f,[shift 0]),2);
        Hdiag(:,nband-shift)=circshift(elements,-shift);
        Hdiag(:,nband+shift)=conj(elements);
        Hdiag(:,shift)=conj(elements);
        Hdiag(:,end+1-shift)=circshift(elements,-shift);
    end
    HZ=HZ+kron(p2,spdiags(Hdiag,-nband+1:nband-1,nband,nband));
end
if length(unique(B(3,:)))==1                  %uniform Bz
    HZ=HZ+kron(p3,spdiags(B(3,:).',0,nband,nband));
else
    pMat=ones(nband,1)*B(3,:);
    Hdiag=zeros(nband,2*nband-1);
    Hdiag(:,nband)=sum(pMat.*conjf.*f,2);
    for shift=1:nband2
        elements=sum(pMat.*conjf.*circshift(f,[shift 0]),2);
        Hdiag(:,nband-shift)=circshift(elements,-shift);
        Hdiag(:,nband+shift)=conj(elements);
        Hdiag(:,shift)=conj(elements);
        Hdiag(:,end+1-shift)=circshift(elements,-shift);
    end
    HZ=HZ+kron(p3,spdiags(Hdiag,-nband+1:nband-1,nband,nband));
end
HZ=(HZ+HZ')/2;
if nargin>=6
    HZ=HZ+kron(p1,spdiags(lambda(1)*cos(2*nlist*pi/N).',0,nband,nband));
    HZ=HZ+kron(p1,spdiags(-lambda(1)*ones(N,1),0,nband,nband));
    HZ=HZ+kron(p2,spdiags(lambda(2)*cos(2*nlist*pi/N).',0,nband,nband));
    HZ=HZ+kron(p2,spdiags(-lambda(2)*ones(N,1),0,nband,nband));
    HZ=HZ+kron(p3,spdiags(lambda(3)*cos(2*nlist*pi/N).',0,nband,nband));
    HZ=HZ+kron(p3,spdiags(-lambda(3)*ones(N,1),0,nband,nband));
end

H=kron(p3,Hv-Hmu)+kron([0 1;0 0],HD)+kron([0 0;1 0],HD')+kron(speye(2),HZ);
% sigy=kron([0 -i;i 0],speye(nband));
% P=Hv-Hmu+HZ;
% Hole=-sigy*conj(Hv-Hmu+HZ)*sigy;
% H=kron(sparse([1 0;0 0]),P)+kron(sparse([0 0;0 1]),Hole)+kron(sparse([0 1;0 0]),HD)+kron(sparse([0 0;1 0]),HD');
f=f.';
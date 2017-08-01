function [H,f]=H_1Ds_Fourier(tneg,mu,SOx,B,DeltaS,bc,shift_mu,nband,lambda)
% basis=(ph)x(spin)x(x)
%INPUTS:
% tneg: 1x1 negative hopping matrix element
% mu: 1xN or 1x1 chemical potential
% SOx: 3x1 coefficients of spin matrices coupled to the momentum along wire
% B: 3xN or 3x1 magnetic field
% DeltaS: 1xN or 1x1 pairing
% shift_on_site_energy: whether to measure from band bottom
% lambda: spni-dependent hopping to gap out band edge for TI system
%at least 1 parameter must reflect the size of the system
%OUTPUTS:
% sparse H whose u's and v's at x for eigen mode n is given by:
%f=f(x,:); 
%u_u = f*EVEC(1:nbands,n)
%u_d = f*EVEC(nbands+1:2*nbands,n)
%v_d = conj(f)*EVEC(2*nbands+1:3*nbands,n)
%-v_u = conj(f)*EVEC(3*nbands+1:4*nbands,n)
N=max([length(tneg) size(SOx,2) length(mu) size(B,2) length(DeltaS)]);
vm=mean(SOx,2);vfluct=SOx-vm*ones(1,size(SOx,2));
if nargin<6, bc=0; end
if nargin<7, shift_mu=0; end
if nargin<8||nband==0, nband=N; end
if nargin<9, lambda=0; end
if shift_mu==1,mu = mu + tneg + circshift(tneg,[1 1]);end
nband=min(nband,N);nband2=ceil((nband-1)/2);
p1=sparse([0 1;1 0]);p2=sparse([0 -1i;1i 0]);p3=sparse([1 0;0 -1]);
if bc==0
%     nlist=[]; %if cutoff is supplied
%     for n=1:N, if abs(2*tneg*cos(n*pi/(N+1))-mean(mu))<=co
%             nlist=[nlist n];
%     end;end
%     nbands=length(nlist);
    nlist=1:N;
    [~,ind]=sort(abs(2*mean(tneg)*cos(nlist*pi/(N+1))-mean(mu))); %treat SO as perturbations
    nlist=nlist(ind);    nlist=nlist(1:nband);
    if nargout>1 || length(unique(tneg))~=1 || length(unique(mu))~=1 || length(unique(DeltaS))~=1 || length(unique(SOx(:,:,1)))~=1|| length(unique(SOx(:,:,2)))~=1|| length(unique(SOx(:,:,3)))~=1|| length(unique(B(:,:,1)))~=1|| length(unique(B(:,:,2)))~=1|| length(unique(B(:,:,3)))~=1
        nMat=repmat(nlist.',[1 N]); xMat=repmat(1:N,[nband 1]);
        f=sqrt(2/(N+1))*sin(pi*nMat.*xMat./(N+1)); %f(bands,x)
    end
    if length(unique(tneg))==1
        HKE=kron(speye(2),spdiags(2*tneg(1)*cos(nlist*pi/(N+1)).',0,nband,nband));
    else
        KEMat=ones(nband,1)*tneg(:).';
        HKEdiag=zeros(nband,2*nband-1);
        HKEdiag(:,nband)=sum(circshift(KEMat,[0 1]).*f.*[zeros(nband,1),f(:,1:end-1)],2)+ ...
                          sum(circshift(KEMat,[0 0]).*f.*[f(:,2:end) zeros(nband,1)],2);
        for shift=1:nband2
            tem=sum(circshift(KEMat,[0 1]).*f.*circshift([zeros(nband,1),f(:,1:end-1)],shift),2)+ ...
                sum(circshift(KEMat,[0 0]).*f.*circshift([f(:,2:end) zeros(nband,1)],shift),2);
            HKEdiag(:,nband-shift)=circshift(tem,-shift);
            HKEdiag(:,nband+shift)=conj(tem);
            HKEdiag(:,shift)=conj(tem);
            HKEdiag(:,end+1-shift)=circshift((tem),-shift);
        end
        HKE=spdiags(HKEdiag,-nband+1:nband-1,nband,nband);
        HKE=kron(speye(2),(HKE+HKE')/2);
    end
    if (length(unique(SOx(1,:)))==1)&& (length(unique(SOx(2,:)))==1)&& (length(unique(SOx(3,:)))==1)
        Hso=zeros(2*nband,2*nband);
        for i1=1:nband
            n1=nlist(i1);
            for i2=i1+1:nband
                n2=nlist(i2);
                factor=2*(1-(-1)^mod(n1+n2,2))/(N+1)*sin(n1*pi/(N+1))*sin(n2*pi/(N+1))/(cos(n1*pi/(N+1))-cos(n2*pi/(N+1)));
                Hso(i1,nband+i2)=-i*0.5*(SOx(1)-i*SOx(2))*factor; %pu-pd
                Hso(nband+i1,i2)=-i*0.5*(SOx(1)+i*SOx(2))*factor;%pd-pu
                Hso(i1,i2)= -i*0.5*SOx(3)*factor; %pu-pu
                Hso(nband+i1,nband+i2)=-i*0.5*(-SOx(3))*factor;%pd-pd
            end
        end
        Hso=Hso+Hso';
    else
        Hso=spdiags([],[],2*nband,2*nband);
        if any(SOx(1,:))
            pMat=ones(nband,1)*SOx(1,:);
            Hdiag=zeros(nband,2*nband-1);
            Hdiag(:,nband)=-0.5*i*(sum(circshift(pMat,[0 1]).*f.*[zeros(nband,1),f(:,1:end-1)]-...
                                           circshift(pMat,[0 0]).*f.*[f(:,2:end) zeros(nband,1)],2));
            for shift=1:nband2
                tem=-0.5*i*(sum(circshift(pMat,[0 1]).*f.*circshift([zeros(nband,1),f(:,1:end-1)],shift)-...
                                circshift(pMat,[0 0]).*f.*circshift([f(:,2:end) zeros(nband,1)],shift),2));
                Hdiag(:,nband-shift)=circshift(tem,-shift);
                Hdiag(:,nband+shift)=conj(tem);
                Hdiag(:,shift)=conj(tem);
                Hdiag(:,end+1-shift)=circshift((tem),-shift);
            end
            Hso=Hso+kron([0 1;1 0],spdiags(Hdiag,-nband+1:nband-1,nband,nband));
        end
        if any(SOx(2,:))
            pMat=ones(nband,1)*SOx(2,:);
            Hdiag=zeros(nband,2*nband-1);
            Hdiag(:,nband)=-0.5*i*(sum(circshift(pMat,[0 1]).*f.*[zeros(nband,1),f(:,1:end-1)]-...
                                           circshift(pMat,[0 0]).*f.*[f(:,2:end) zeros(nband,1)],2));
            for shift=1:nband2
                tem=-0.5*i*(sum(circshift(pMat,[0 1]).*f.*circshift([zeros(nband,1),f(:,1:end-1)],shift)-...
                                circshift(pMat,[0 0]).*f.*circshift([f(:,2:end) zeros(nband,1)],shift),2));
                Hdiag(:,nband-shift)=circshift(tem,-shift);
                Hdiag(:,nband+shift)=conj(tem);
                Hdiag(:,shift)=conj(tem);
                Hdiag(:,end+1-shift)=circshift((tem),-shift);
            end
            Hso=Hso+kron([0 -i;i 0],spdiags(Hdiag,-nband+1:nband-1,nband,nband));
        end
        if any(SOx(3,:))
            pMat=ones(nband,1)*SOx(3,:);
            Hdiag=zeros(nband,2*nband-1);
            Hdiag(:,nband)=-0.5*i*(sum(circshift(pMat,[0 1]).*f.*[zeros(nband,1),f(:,1:end-1)]-...
                                           circshift(pMat,[0 0]).*f.*[f(:,2:end) zeros(nband,1)],2));
            for shift=1:nband2
                tem=-0.5*i*(sum(circshift(pMat,[0 1]).*f.*circshift([zeros(nband,1),f(:,1:end-1)],shift)-...
                                circshift(pMat,[0 0]).*f.*circshift([f(:,2:end) zeros(nband,1)],shift),2));
                Hdiag(:,nband-shift)=circshift(tem,-shift);
                Hdiag(:,nband+shift)=conj(tem);
                Hdiag(:,shift)=conj(tem);
                Hdiag(:,end+1-shift)=circshift((tem),-shift);
            end
            Hso=Hso+kron([1 0;0 -1],spdiags(Hdiag,-nband+1:nband-1,nband,nband));
        end
        Hso=(Hso+Hso')/2;
    end
    if length(unique(mu))==1                  %uniform mu
        Hmu=kron(speye(2),spdiags(mu(1)*ones(N,1),0,nband,nband));
    else
        pMat=ones(nband,1)*mu(:).';
        Hdiag=zeros(nband,2*nband-1);
        Hdiag(:,nband)=sum(pMat.*f.*f,2);
        for shift=1:nband2
            elements=sum(pMat.*f.*circshift(f,[shift 0]),2);
            Hdiag(:,nband+shift)=elements;
            Hdiag(:,shift)=elements;
            Hdiag(:,nband-shift)=circshift(elements,-shift);
            Hdiag(:,end+1-shift)=circshift(elements,-shift);
        end
        Hmu=spdiags(Hdiag,-nband+1:nband-1,nband,nband);
        Hmu=kron(speye(2),(Hmu+Hmu')/2);
    end
    if length(unique(DeltaS))==1              %uniform DeltaS
        HD=kron(speye(2),spdiags(DeltaS(1)*ones(N,1),0,nband,nband));
    else
        pMat=ones(nband,1)*DeltaS(:).';
        Hdiag=zeros(nband,2*nband-1);
        Hdiag(:,nband)=sum(pMat.*f.*f,2);
        for shift=1:nband2
            elements=sum( pMat.*f.*circshift(f,[shift 0]),2);
            Hdiag(:,nband+shift)=elements;
            Hdiag(:,shift)=elements;
            Hdiag(:,nband-shift)=circshift(elements,-shift);
            Hdiag(:,end+1-shift)=circshift(elements,-shift);
        end
        HD=kron(speye(2),spdiags(Hdiag,-nband+1:nband-1,nband,nband));
    end
    HZ=spdiags([],[],2*nband,2*nband);
    if length(unique(B(1,:)))==1              %uniform Bx
        HZ=HZ+kron(p1,spdiags(B(1)*ones(N,1),0,nband,nband));
    else
        pMat=ones(nband,1)*B(1,:);
        Hdiag=zeros(nband,2*nband-1);
        Hdiag(:,nband)=sum(pMat.*f.*f,2);
        for shift=1:nband2
            elements=sum(pMat.*f.*circshift(f,[shift 0]),2);
            Hdiag(:,nband+shift)=elements;
            Hdiag(:,shift)=elements;
            Hdiag(:,nband-shift)=circshift(elements,-shift);
            Hdiag(:,end+1-shift)=circshift(elements,-shift);
        end
        HZ=HZ+kron(p1,spdiags(Hdiag,-nband+1:nband-1,nband,nband));
    end
    if length(unique(B(2,:)))==1              %uniform By
        HZ=HZ+kron(p2,spdiags(B(2)*ones(N,1),0,nband,nband));
    else
        pMat=ones(nband,1)*B(2,:);
        Hdiag=zeros(nband,2*nband-1);
        Hdiag(:,nband)=sum(pMat.*f.*f,2);
        for shift=1:nband2
            elements=sum(pMat.*f.*circshift(f,[shift 0]),2);
            Hdiag(:,nband+shift)=elements;
            Hdiag(:,shift)=elements;
            Hdiag(:,nband-shift)=circshift(elements,-shift);
            Hdiag(:,end+1-shift)=circshift(elements,-shift);
        end
        HZ=HZ+kron(p2,spdiags(Hdiag,-nband+1:nband-1,nband,nband));
    end
    if length(unique(B(3,:)))==1              %uniform Bz
        HZ=HZ+kron(p3,spdiags(B(3)*ones(N,1),0,nband,nband));
    else
        pMat=ones(nband,1)*B(3,:);
        Hdiag=zeros(nband,2*nband-1);
        Hdiag(:,nband)=sum(pMat.*f.*f,2);
        for shift=1:nband2
            elements=sum(pMat.*f.*circshift(f,[shift 0]),2);
            Hdiag(:,nband+shift)=elements;
            Hdiag(:,shift)=elements;
            Hdiag(:,nband-shift)=circshift(elements,-shift);
            Hdiag(:,end+1-shift)=circshift(elements,-shift);
        end
        HZ=HZ+kron(p3,spdiags(Hdiag,-nband+1:nband-1,nband,nband));
    end
    if nargin>=9
        HZ=HZ+kron(p1,spdiags(lambda(1)*cos(nlist*pi/(N+1)).',0,nband,nband));
        HZ=HZ+kron(p1,spdiags(-lambda(1)*ones(N,1),0,nband,nband));
        HZ=HZ+kron(p2,spdiags(lambda(2)*cos(nlist*pi/(N+1)).',0,nband,nband));
        HZ=HZ+kron(p2,spdiags(-lambda(2)*ones(N,1),0,nband,nband));
        HZ=HZ+kron(p3,spdiags(lambda(3)*cos(nlist*pi/(N+1)).',0,nband,nband));
        HZ=HZ+kron(p3,spdiags(-lambda(3)*ones(N,1),0,nband,nband));
    end
    HZ=(HZ+HZ')/2;
end
if bc==1 || bc==-1
    nlist=(round((-N+1.5)/2):round((N-0.5)/2))-atan2(imag(bc),real(bc))/2/pi;
    if any(tneg) || (nargin>=9 && length(lambda)==3)
        [~,ind]=sort(abs(sin(2*nlist*pi/N)*sqrt(sum(mean(SOx,2).^2))+2*mean(tneg)*cos(2*nlist*pi/N)-mean(mu))); %treats normal-hopping + SO together
    else
        [~,ind]=sort(abs((2*nlist*pi/N)*sqrt(sum(mean(SOx,2).^2))-mean(mu)));
    end
    nlist=nlist(ind);    nlist=nlist(1:nband);%E(nbands) %current cutoff
    if nargout>1 || length(unique(tneg))~=1 || length(unique(mu))~=1 || length(unique(DeltaS))~=1 || length(unique(SOx(:,:,1)))~=1|| length(unique(SOx(:,:,2)))~=1|| length(unique(SOx(:,:,3)))~=1|| length(unique(B(:,:,1)))~=1|| length(unique(B(:,:,2)))~=1|| length(unique(B(:,:,3)))~=1
        nMat=repmat(nlist.',[1 N]);    xMat=repmat(1:N,[nband 1]);
        f=sqrt(1/N)*exp(1i*2*pi*nMat.*xMat./N);    conjf=conj(f); %f(bands,x)
    end
    if length(unique(tneg))==1.
        HKE=kron(speye(2),spdiags(2*tneg(1)*cos(2*nlist(:)*pi/N),0,nband,nband));
    else
        KEMat=ones(nband,1)*tneg; 
        HKEdiag=zeros(nband,2*nband-1);
        HKEdiag(:,nband)=sum(circshift(KEMat,[0 1]).*conjf.*[bc*f(:,end) f(:,1:end-1)],2)+ ...
                          sum(circshift(KEMat,[0 0]).*conjf.*[f(:,2:end) bc*f(:,1)],2);
        for shift=1:nband2
            tem=sum(circshift(KEMat,[0 1]).*conjf.*circshift([bc*f(:,end) f(:,1:end-1)],[shift 0]),2)+...
                sum(circshift(KEMat,[0 0]).*conjf.*circshift([f(:,2:end) bc*f(:,1)],[shift 0]),2);
            HKEdiag(:,nband-shift)=circshift(tem,-shift);
            HKEdiag(:,nband+shift)=conj(tem);
            HKEdiag(:,shift)=conj(tem);
            HKEdiag(:,end+1-shift)=circshift((tem),-shift);
        end
        HKE=spdiags(HKEdiag,-nband+1:nband-1,nband,nband);
        HKE=kron(speye(2),(HKE+HKE')/2);
    end
    if any(tneg) || (nargin>=9 && length(lambda)==3),factor=sin(2*nlist(:)*pi/N); else factor=2*nlist(:)*pi/N;end
    Hso=spdiags([-(vm(1)+i*vm(2))*factor -vm(3)*factor -(vm(1)-i*vm(2))*factor;
                 -(vm(1)+i*vm(2))*factor  vm(3)*factor -(vm(1)-i*vm(2))*factor],[-nband 0 nband],2*nband,2*nband);    
    if any(vfluct(1,:))
        pMat=ones(nband,1)*vfluct(1,:); 
        Hdiag=zeros(nband,2*nband-1);
        Hdiag(:,nband)=-0.5*i*(sum(circshift(pMat,[0 1]).*conjf.*[bc*f(:,end) f(:,1:end-1)],2)- ...
                               sum(circshift(pMat,[0 0]).*conjf.*[f(:,2:end) bc*f(:,1)],2));
        for shift=1:nband2
            tem=-0.5*i*(sum(circshift(pMat,[0 1]).*conjf.*circshift([bc*f(:,end) f(:,1:end-1)],[shift 0]),2)- ...
                sum(circshift(pMat,[0 0]).*conjf.*circshift([f(:,2:end) bc*f(:,1)],[shift 0]),2));
            Hdiag(:,nband-shift)=circshift(tem,-shift);
            Hdiag(:,nband+shift)=conj(tem);
            Hdiag(:,shift)=conj(tem);
            Hdiag(:,end+1-shift)=circshift((tem),-shift);
        end
        Hso=Hso+kron([0 1;1 0],spdiags(Hdiag,-nband+1:nband-1,nband,nband));
    end
    if any(vfluct(2,:))
        pMat=ones(nband,1)*vfluct(2,:); 
        Hdiag=zeros(nband,2*nband-1);
        Hdiag(:,nband)=-0.5*i*(sum(circshift(pMat,[0 1]).*conjf.*[bc*f(:,end) f(:,1:end-1)],2)- ...
                               sum(circshift(pMat,[0 0]).*conjf.*[f(:,2:end) bc*f(:,1)],2));
        for shift=1:nband2
            tem=-0.5*i*(sum(circshift(pMat,[0 1]).*conjf.*circshift([bc*f(:,end) f(:,1:end-1)],[shift 0]),2)- ...
                sum(circshift(pMat,[0 0]).*conjf.*circshift([f(:,2:end) bc*f(:,1)],[shift 0]),2));
            Hdiag(:,nband-shift)=circshift(tem,-shift);
            Hdiag(:,nband+shift)=conj(tem);
            Hdiag(:,shift)=conj(tem);
            Hdiag(:,end+1-shift)=circshift((tem),-shift);
        end
        Hso=Hso+kron([0 -i;i 0],spdiags(Hdiag,-nband+1:nband-1,nband,nband));
    end
    if any(vfluct(3,:))
        pMat=ones(nband,1)*vfluct(3,:); 
        Hdiag=zeros(nband,2*nband-1);
        Hdiag(:,nband)=-0.5*i*(sum(circshift(pMat,[0 1]).*conjf.*[bc*f(:,end) f(:,1:end-1)],2)- ...
                               sum(circshift(pMat,[0 0]).*conjf.*[f(:,2:end) bc*f(:,1)],2));
        for shift=1:nband2
            tem=-0.5*i*(sum(circshift(pMat,[0 1]).*conjf.*circshift([bc*f(:,end) f(:,1:end-1)],[shift 0]),2)- ...
                sum(circshift(pMat,[0 0]).*conjf.*circshift([f(:,2:end) bc*f(:,1)],[shift 0]),2));
            Hdiag(:,nband-shift)=circshift(tem,-shift);
            Hdiag(:,nband+shift)=conj(tem);
            Hdiag(:,shift)=conj(tem);
            Hdiag(:,end+1-shift)=circshift((tem),-shift);
        end
        Hso=Hso+kron([1 0;0 -1],spdiags(Hdiag,-nband+1:nband-1,nband,nband));
    end
    Hso=(Hso+Hso')/2;
    if length(unique(mu))==1                  %uniform mu
        Hmu=kron(speye(2),spdiags(mu(1)*ones(N,1),0,nband,nband));
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
    if length(unique(DeltaS))==0.1              %DeltaS - Always need to compute overlaps because of the wrong choice of basis
        HD=kron(speye(2),spdiags(DeltaS(1)*ones(N,1),0,nband,nband));
    else
        pMat=ones(nband,1)*(DeltaS);
        Hdiag=zeros(nband,2*nband-1);
        Hdiag(:,nband)=sum(pMat.*conjf.*conjf,2);
        for shift=1:nband2
            elements=sum( pMat.*conjf.*circshift(conjf,[shift 0]),2);
            Hdiag(:,nband-shift)=circshift((elements),-shift);
            Hdiag(:,nband+shift)=circshift((elements),0);
            Hdiag(:,shift)=(elements);
            Hdiag(:,end+1-shift)=circshift(elements,-shift);
        end
%         Hdiag(:,nband)=sum(pMat.*conjf.*f,2);
%         for shift=1:nband2
%             elements=sum( pMat.*conjf.*circshift(f,[shift 0]),2);
%             conjelements=sum( pMat.*f.*circshift(conjf,[shift 0]),2);
%             Hdiag(:,nband-shift)=circshift(elements,-shift);
%             Hdiag(:,nband+shift)=conjelements;
%             Hdiag(:,shift)=conjelements;
%             Hdiag(:,end+1-shift)=circshift(elements,-shift);
%         end
        HD=kron(speye(2),spdiags(Hdiag,-nband+1:nband-1,nband,nband));
    end
    HZ=spdiags([],[],2*nband,2*nband);
    if length(unique(B(1,:)))==1                  %uniform Bx
        HZ=HZ+kron(p1,spdiags(B(1)*ones(N,1),0,nband,nband));
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
        HZ=HZ+kron(p2,spdiags(B(2)*ones(N,1),0,nband,nband));
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
        HZ=HZ+kron(p3,spdiags(B(3)*ones(N,1),0,nband,nband));
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
    if nargin>=9 && length(lambda)==3
        HZ=HZ+kron(p1,spdiags(lambda(1)*cos(2*nlist*pi/N).',0,nband,nband));
        HZ=HZ+kron(p1,spdiags(-lambda(1)*ones(N,1),0,nband,nband));
        HZ=HZ+kron(p2,spdiags(lambda(2)*cos(2*nlist*pi/N).',0,nband,nband));
        HZ=HZ+kron(p2,spdiags(-lambda(2)*ones(N,1),0,nband,nband));
        HZ=HZ+kron(p3,spdiags(lambda(3)*cos(2*nlist*pi/N).',0,nband,nband));
        HZ=HZ+kron(p3,spdiags(-lambda(3)*ones(N,1),0,nband,nband));
    end
end

sigy=kron([0 -i;i 0],speye(nband));
P=HKE+Hso-Hmu+HZ;H=kron(sparse([1 0;0 0]),P)+kron(sparse([0 0;0 1]),-sigy*conj(P)*sigy)+kron(sparse([0 1;0 0]),HD)+kron(sparse([0 0;1 0]),HD');
% H=kron(p3,HKE+Hso-Hmu)+kron([0 1;0 0],HD)+kron([0 0;1 0],HD')+kron(speye(2),HZ);
% if bc==1 || bc==-1 %Transform to PH sym, if nband = N
%     [r,c]=find(nlist'*ones(1,nband)+ones(nband,1)*nlist==0);r=r';
%     if length(r)<nband, miss=find([c;0]~=(1:nband)',1);r=[r(1:miss-1) miss r(miss:end)];end
%     o=[1:2*nband r+2*nband r+3*nband];
%     H=H(o,o);
% end
if nargout>1
    f=f.';
end
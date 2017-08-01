function [H,f]=H_1Ds_Fourier_TI_multi(vc,mu0c,muxc,muyc,muzc,Bc,DeltaSc,nband)
% basis=(ph)x(spin)x(k)
%INPUTS:
% vc: 3nxN coefficients of spin matrices coupled to the momentum along wire; where n = # channels
% mu0c: nxnxN symmetric scattering M that couples to Identity Matrix in spin space
% muxc: nxnxN antisymmetric imag M that couples to sigma_x
% muyc: nxnxN antisymmetric imag M that couples to sigma_y
% muzc: nxnxN antisymmetric imag M that couples to sigma_z
% Bc: 3nxN or 3nx1 magnetic field
% DeltaSc: nxN or nx1 pairing
%at least 1 parameter must reflect the size of the system
%OUTPUTS:
% sparse H whose u's and v's at x for eigen mode n is given by:
%f=f(:,x).';
%u_u = f*EVEC(1:nbands,n)
%u_d = f*EVEC(nbands+1:2*nbands,n)
%v_d = f*EVEC(2*nbands+1:3*nbands,n)
%-v_u = f*EVEC(3*nbands+1:4*nbands,n)
% vc,mu0c(:),Bc(:),DeltaSc(:)
N=max([size(vc,2) size(mu0c,3) size(Bc,2) size(DeltaSc,2)]);
Nch=max([size(vc,1)/3 size(mu0c,1) size(Bc,1)/3 size(DeltaSc,1)]);
if nband==0, nband=N; end
if size(Bc,2)==1, Bc=Bc*ones(1,N); end
if size(DeltaSc,2)==1, DeltaSc=DeltaSc*ones(1,N); end
nband=min(nband,N); nband2=ceil((nband-1)/2);
p1=sparse([0 1;1 0]);p2=sparse([0 -1i;1i 0]);p3=sparse([1 0;0 -1]);
Hchannels=cell(Nch,Nch); fchannels=cell(Nch,1);

for chi=1:Nch
    nlist=round((-N+1.5)/2):round((N-0.5)/2);
    v=vc(3*chi-2:3*chi,:);
    vm=mean(v,2);
    vfluct=v-vm*ones(1,size(v,2));
    mu0=reshape(mu0c(chi,chi,:),1,N);
    B=Bc(3*chi-2:3*chi,:);
    DeltaS=DeltaSc(chi,:);
    if nband<N
        [~,ind]=sort(abs(0.5*(2*nlist*pi/N)*sqrt(sum(sum(v.^2)))-mean(mu0)));
        nlist=nlist(ind);    nlist=nlist(1:nband);%E(nbands) %current cutoff
    end
    nMat=repmat(nlist.',[1 N]); xMat=repmat(1:N,[nband 1]);
    f=sqrt(1/N)*exp(1i*2*pi*nMat.*xMat./N);    conjf=conj(f); %f(bands,x)
    factor=2*nlist(:)*pi/N;
    Hv=spdiags([-(vm(1)+i*vm(2))*factor -vm(3)*factor -(vm(1)-i*vm(2))*factor;
                -(vm(1)+i*vm(2))*factor  vm(3)*factor -(vm(1)-i*vm(2))*factor],[-nband 0 nband],2*nband,2*nband);
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
        Hv=Hv+kron([0 -i;i 0],spdiags(Hdiag,-nband+1:nband-1,nband,nband));
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
    if length(unique(mu0))==1                  %uniform mu
        Hmu=kron(speye(2),spdiags(mu0.',0,nband,nband));
    else
        pMat=ones(nband,1)*mu0;
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
    H=kron(p3,Hv-Hmu)+kron([0 1;0 0],HD)+kron([0 0;1 0],HD')+kron(speye(2),HZ);
    Hchannels{chi,chi}=H;
    fchannels{chi,1}=f;
end
for chi=1:Nch
    for chj=chi+1:Nch
        fj=fchannels{chj,1}; conjfi=conj(fchannels{chi,1});
        %fi=fchannels{chi,1}; conjfj=conj(fj);
        mu0=reshape(mu0c(chi,chj,:),1,N);
        mux=reshape(muxc(chi,chj,:),1,N);
        muy=reshape(muyc(chi,chj,:),1,N);
        muz=reshape(muzc(chi,chj,:),1,N);
        B=Bc(3*chi-2:3*chi,:);
        if nnz(mu0)==0
            Hmu=sparse(2*nband,2*nband);
        else
            pMat=ones(nband,1)*mu0;
            Hdiag=zeros(nband,2*nband-1);
            Hdiag(:,nband)=sum(pMat.*conjfi.*fj,2);
            for shift=1:nband2
                elements=sum(pMat.*conjfi.*circshift(fj,[shift 0]),2);
                conjelements=sum(pMat.*fj.*circshift(conjfi,[shift 0]),2);
                Hdiag(:,nband-shift)=circshift(elements,-shift);
                Hdiag(:,nband+shift)=conjelements;
                Hdiag(:,shift)=conjelements;
                Hdiag(:,end+1-shift)=circshift(elements,-shift);
            end
            Hmu=spdiags(Hdiag,-nband+1:nband-1,nband,nband);
            Hmu=kron(speye(2),Hmu);
        end
        if nnz(mux)==0
            Hmux=sparse(2*nband,2*nband);
        else
            pMat=ones(nband,1)*mux;
            Hdiag=zeros(nband,2*nband-1);
            Hdiag(:,nband)=sum(pMat.*conjfi.*fj,2);
            for shift=1:nband2
                elements=sum(pMat.*conjfi.*circshift(fj,[shift 0]),2);
                conjelements=sum(pMat.*fj.*circshift(conjfi,[shift 0]),2);
                Hdiag(:,nband-shift)=circshift(elements,-shift);
                Hdiag(:,nband+shift)=conjelements;
                Hdiag(:,shift)=conjelements;
                Hdiag(:,end+1-shift)=circshift(elements,-shift);
            end
            Hmux=spdiags(Hdiag,-nband+1:nband-1,nband,nband);
            Hmux=kron(p1,Hmux);
        end
        if nnz(muy)==0
            Hmuy=sparse(2*nband,2*nband);
        else
            pMat=ones(nband,1)*muy;
            Hdiag=zeros(nband,2*nband-1);
            Hdiag(:,nband)=sum(pMat.*conjfi.*fj,2);
            for shift=1:nband2
                elements=sum(pMat.*conjfi.*circshift(fj,[shift 0]),2);
                conjelements=sum(pMat.*fj.*circshift(conjfi,[shift 0]),2);
                Hdiag(:,nband-shift)=circshift(elements,-shift);
                Hdiag(:,nband+shift)=conjelements;
                Hdiag(:,shift)=conjelements;
                Hdiag(:,end+1-shift)=circshift(elements,-shift);
            end
            Hmuy=spdiags(Hdiag,-nband+1:nband-1,nband,nband);
            Hmuy=kron(p2,Hmuy);
        end
        if nnz(muz)==0
            Hmuz=sparse(2*nband,2*nband);
        else
            pMat=ones(nband,1)*muz;
            Hdiag=zeros(nband,2*nband-1);
            Hdiag(:,nband)=sum(pMat.*conjfi.*fj,2);
            for shift=1:nband2
                elements=sum(pMat.*conjfi.*circshift(fj,[shift 0]),2);
                conjelements=sum( pMat.*fj.*circshift(conjfi,[shift 0]),2);
                Hdiag(:,nband-shift)=circshift(elements,-shift);
                Hdiag(:,nband+shift)=conjelements;
                Hdiag(:,shift)=conjelements;
                Hdiag(:,end+1-shift)=circshift(elements,-shift);
            end
            Hmuz=spdiags(Hdiag,-nband+1:nband-1,nband,nband);
            Hmuz=kron(p3,Hmuz);
        end

        H=-kron([1 0;0 0],Hmu+Hmux+Hmuy+Hmuz)-kron([0 0;0 1],-Hmu-Hmux-Hmuy-Hmuz);
        Hchannels{chi,chj}=H;
        Hchannels{chj,chi}=H';
    end
end
H=cell2mat(Hchannels);
f=cell2mat(fchannels).';
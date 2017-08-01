set(0,'DefaultFigureWindowStyle','docked')
%% 1D TDOS vs E (PRL.108.067001)
 Nv=[200 200];  bc=0;  nbands=100; N=sum(Nv);
spread=0.01; % the Lorentzian spread for DOS
  t=10;  alp=[0;5*sqrt(0.1*t);0];
  D=[0.5 0.5];
 mu=[.5   .5];  shift_mu=1;
  B=[1 1;0 0;0 0];
 sd=0.0;
  muv=[];Bv=[];Dv=[];
for ii=1:length(Nv)
      z=zeros(1,Nv(ii)); 
    muv=[muv mu(ii)*ones(1,Nv(ii))];
     Bv=[Bv B(:,ii)*ones(1,Nv(ii))];
     Dv=[Dv D(:,ii)*ones(1,Nv(ii))];
end
Ev=linspace(-2,2,1001); runs=1;
DE=zeros(length(Ev),runs); 
for run=1:runs
    noise=sd*randn(1,N);
%     H=H_1Ds(-t,mu,alp,B,del,bc,shift_mu);
    H=H_1Ds_Fourier(-t,muv,alp,Bv,Dv,bc,shift_mu,nbands);
    D=eig(full(H));
%     for ii=1:length(D),DE(:,run)=DE(:,run)+(spread/pi./(spread^2+(Ev(:)-D(ii)).^2));end
    DE(:,run)=sum(spread/pi./(spread^2+(Ev(:)*ones(1,length(D))-ones(length(Ev),1)*D(:).').^2),2);
end
figure(1), plot(Ev,mean(DE,2));
%% 1D Local DOS / Majorana Polarization: Contour vs (E,pos)
 Nv=[100]; bc=0;  n=100; %number of bands in Fourier
spread=0.01;
  t=1;  alp=[0;1;0]*2;
  D=[1];
 mu=[0];  shift_mu=1;
  B=[0;0;2];
 sd=0;
 N=sum(Nv); muv=[];Bv=[];Dv=[];
for ii=1:length(Nv)
      z=zeros(1,Nv(ii)); 
    muv=[muv mu(ii)*ones(1,Nv(ii))];
     Bv=[Bv B(:,ii)*ones(1,Nv(ii))];
     Dv=[Dv D(:,ii)*ones(1,Nv(ii))];
end

Ev=linspace(-1.5,1.5,1001); runs=1;
Pv=1:N;
DOS=zeros(length(Pv),length(Ev),runs); incV=0;
MP=zeros(length(Pv),length(Ev),runs);
for run=1:runs
    H=H_1Ds(-t,muv,alp,Bv,Dv,bc,shift_mu); [V,D]=eig(full(H)); D=diag(D);
    DErun=zeros(length(Pv),length(Ev));
    MPrun=zeros(length(Pv),length(Ev));
    for ii=1:length(D)
        DErun=DErun+(abs(V(Pv,ii)).^2+abs(V(N+Pv,ii)).^2+incV*abs(V(2*N+Pv,ii)).^2+incV*abs(V(3*N+Pv,ii)).^2)*(spread/pi./(spread^2+(Ev-D(ii)).^2));
        MPrun=MPrun+(2*real(conj(V(N+Pv,ii)).*V(2*N+Pv,ii)-conj(V(Pv,ii)).*V(3*N+Pv,ii)))*(spread/pi./(spread^2+(Ev-D(ii)).^2)); % x
    end
    DOS(:,:,run)=DErun;
    MP(:,:,run)=MPrun;
end
DOS=mean(DOS,3); DOS=min(DOS,2);
figure(1),contourf(Ev,Pv,DOS),shading flat,colorbar
%% 1D Local DOS / Majorana Polarization vs E
 Nv=[100]; bc=0;  n=100; %number of bands in Fourier
pos=[1]; spread=0.01;
  t=1;  alp=[0;1;0]*2;
  D=[1];
 mu=[0];  shift_mu=1;
  B=[0;0;2];
 sd=0;
 N=sum(Nv); muv=[];Bv=[];Dv=[];
for ii=1:length(Nv)
      z=zeros(1,Nv(ii)); 
    muv=[muv mu(ii)*ones(1,Nv(ii))];
     Bv=[Bv B(:,ii)*ones(1,Nv(ii))];
     Dv=[Dv D(:,ii)*ones(1,Nv(ii))];
end

Ev=linspace(-1.5,1.5,1001); runs=1;
DOS=zeros(length(pos),length(Ev),runs); incV=0;
MP=zeros(length(pos),length(Ev),runs);
parfor run=1:runs
    noise=sd*randn(1,N); 
     H=H_1Ds(-t,muv+noise,alp,Bv,Dv,bc,shift_mu); [V,D]=eig(full(H)); D=diag(D);
%     [H,f]=H_1Ds_Fourier(-t,muv+noise,alp,Bv,Dv,bc,shift_mu,n); [V,D]=eig(full(H)); D=diag(D);
    DErun=zeros(length(pos),length(Ev));
    MPrun=zeros(length(pos),length(Ev));
    for ii=1:length(D)
        DErun=DErun+(abs(V(pos,ii)).^2+abs(V(N+pos,ii)).^2+incV*abs(V(2*N+pos,ii)).^2+incV*abs(V(3*N+pos,ii)).^2)*(spread/pi./(spread^2+(Ev-D(ii)).^2));
        MPrun=MPrun+(2*real(conj(V(N+pos,ii)).*V(2*N+pos,ii)-conj(V(pos,ii)).*V(3*N+pos,ii)))*(spread/pi./(spread^2+(Ev-D(ii)).^2)); % x
%         DErun=DErun+(abs(f(pos,:)*V(1:n,ii)).^2+abs(f(pos,:)*V(n+1:2*n,ii)).^2+incV*abs(f(pos,:)*V(2*n+1:3*n,ii)).^2+incV*abs(f(pos,:)*V(3*n+1:4*n,ii)).^2)*(spread/pi./(spread^2+(Ev-D(ii)).^2));
%         MPrun=MPrun+(2*real(conj(f(pos,:)*V(n+1:2*n,ii)).*(f(pos,:)*V(2*n+1:3*n,ii))-conj(f(pos,:)*V(1:n,ii)).*(f(pos,:)*V(3*n+1:4*n,ii))))*(spread/pi./(spread^2+(Ev-D(ii)).^2)); % x
    end
    DOS(:,:,run)=DErun;
    MP(:,:,run)=MPrun;
end
DOS=mean(DOS,3);
% MP=mean(MP,3);
% for ii=1:length(pos)
%     figure(ii),plot(Ev,[DOS(ii,:);MP(ii,:)])
% end
figure(1),plot(Ev,DOS(1,:)),axis tight,axis([xlim 0 1])
%% 1D pos-avged LDoS / MP: Contour on (E,parameter)
Ev=linspace(-0.15,0.15,301); parv=linspace(0,0.4,100);%arXiv:1203.2643
Ev=linspace(-2,2,401); parv=linspace(0,3,200);
spread=0.01;
DOS=zeros(length(Ev),length(parv));incV=0;
parfor pari=1:length(parv)
    p=parv(pari);
    Nv=[80 20];n=100;pos=81:85;t=1;alp=[0;.2;0];D=[0.3 0];mu=[p p];shift_mu=1;B=[0.4 0.4;0 0;0 0];%arXiv:1203.2643
    Nv=300;n=100;pos=150;t=5;alp=[0;.2;0];D=1;mu=1.25;shift_mu=1;B=[p;0;0];%arXiv:1203.2643
    
    N=sum(Nv);muv=[];Bv=[];Dv=[];
    for ii=1:length(Nv)
        muv=[muv mu(ii)*ones(1,Nv(ii))];
         Bv=[Bv B(:,ii)*ones(1,Nv(ii))];
         Dv=[Dv D(:,ii)*ones(1,Nv(ii))];
    end
    [H,f]=H_1Ds_Fourier(-t,muv,alp,Bv,Dv,bc,shift_mu,n); [V,D]=eig(full(H)); D=diag(D);
    col=zeros(length(pos),length(Ev));
    for ii=1:length(D)
        col=col+(abs(f(pos,:)*V(1:n,ii)).^2+abs(f(pos,:)*V(n+1:2*n,ii)).^2+incV*abs(f(pos,:)*V(2*n+1:3*n,ii)).^2+incV*abs(f(pos,:)*V(3*n+1:4*n,ii)).^2)*(spread/pi./(spread^2+(Ev-D(ii)).^2));
%         col=col+(2*real(conj(f(pos,:)*V(n+1:2*n,ii)).*(f(pos,:)*V(2*n+1:3*n,ii))-conj(f(pos,:)*V(1:n,ii)).*(f(pos,:)*V(3*n+1:4*n,ii))))*(spread/pi./(spread^2+(Ev-D(ii)).^2)); %MPx
    end
    DOS(:,pari)=mean(col,1).';
end
contourf(parv,Ev,DOS),shading flat,colormap hot
%% 1D Wavefn / Majorana Polarization vs pos
 Nv=[200]; bc=0; N=sum(Nv); n=N/2;
 Ei=[200 201];
  t=10;  alp=[0;5;0];
  D=[4];
 mu=[3];  shift_mu=1;
  B=[5.5];
 sd=0.1;
 Nv=[80 20]*1; N=sum(Nv);n=N;Ei=2*N;t=1;alp=[0;.2;0];D=[0.3 0];mu=[0 0];shift_mu=1;B=[0 0;0 0;0.4 0.4];sd=0;%arXiv:1203.2643
  muv=[];Bv=[];Dv=[];
for ii=1:length(Nv)
    muv=[muv mu(ii)*ones(1,Nv(ii))];
     Bv=[Bv B(:,ii)*ones(1,Nv(ii))];
     Dv=[Dv D(:,ii)*ones(1,Nv(ii))];
end

runs=1;Nv=1:N;
DOS=zeros(length(Ei),N,runs); incV=1;
MP=zeros(length(Ei),N,runs);
for run=1:runs
    noise=sd*randn(1,N); 
%      H=H_1Ds(-t,muv+noise,alp,Bv,Dv,bc,shift_mu); [V,D]=eig(full(H)); D=diag(D);
    [H,f]=H_1Ds_Fourier(-t,muv+noise,alp,Bv,Dv,bc,shift_mu,n); [V,D]=eig(full(H)); D=diag(D);D(Ei)
    DErun=zeros(length(Ei),N);
    MPrun=zeros(length(Ei),N);
%    DErun=DErun+(abs(V(Nv,Ei)).^2+abs(V(N+Nv,Ei)).^2+incV*abs(V(2*N+Nv,Ei)).^2+incV*abs(V(3*N+Nv,Ei)).^2)*(spread/pi./(spread^2+(Ev-D(ii)).^2));
%    MPrun=MPrun+(2*real(conj(V(N+Nv,Ei))*V(2*N+Nv,Ei)-conj(V(Nv,Ei))*V(3*N+Nv,Ei)))*(spread/pi./(spread^2+(Ev-D(ii)).^2)); % x
    DErun=DErun+(abs(f(Nv,:)*V(1:n,Ei)).^2+abs(f(Nv,:)*V(n+1:2*n,Ei)).^2+incV*abs(f(Nv,:)*V(2*n+1:3*n,Ei)).^2+incV*abs(f(Nv,:)*V(3*n+1:4*n,Ei)).^2).';
    MPrun=MPrun+(2*real(conj(f(Nv,:)*V(n+1:2*n,Ei)).*(f(Nv,:)*V(2*n+1:3*n,Ei))-conj(f(Nv,:)*V(1:n,Ei)).*(f(Nv,:)*V(3*n+1:4*n,Ei)))).'; % x
    DOS(:,:,run)=DErun;
    MP(:,:,run)=MPrun;
end
DOS=mean(DOS,3);
MP=mean(MP,3);
for ii=1:length(Ei)
    figure(ii),plot(Nv,[DOS(ii,:);MP(ii,:)])
end
%% 1D Conductance by Meir-Wingreen vs E
  N=[0 200];  bc=0; n=sum(N)/2; %number of bands in Fourier
pos=[1 50]; GAMMA=2; Dynes=0.000;
  t=10;  alp=[0;2*sqrt(0.1*t);0];
Del=[0.5 0.5];
 mu=[.5   .5];  shift_mu=1;
  B=[1 1;0 0;0 0];
 sd=0.0;
  muv=[];Bv=[];Dv=[];
for ii=1:length(N)
    muv=[muv mu(ii)*ones(1,N(ii))];
     Bv=[Bv B(:,ii)*ones(1,N(ii))];
     Dv=[Dv Del(:,ii)*ones(1,N(ii))];
end
N=sum(N);
Ev=linspace(-2*Del(1),2*Del(1),1001); runs=1;
DE=zeros(length(pos),length(Ev),runs); 
for run=1:runs
    noise=sd*randn(1,N); 
    H=H_1Ds(-t,muv+noise,alp,Bv,Dv,bc,shift_mu); [V,D]=eig(full(H)); D=diag(D);
%     [H,f]=H_1Ds_Fourier(-t,muv+noise,alp,Bv,Dv,bc,shift_mu,n); [V,D]=eig(full(H)); D=diag(D);
    DErun=zeros(length(pos),length(Ev));
    for ii=1:length(D)
        u2=repmat(abs(V(pos,ii)).^2     + abs(V(N+pos,ii)).^2,[1 length(Ev)]);
        v2=repmat(abs(V(2*N+pos,ii)).^2 + abs(V(3*N+pos,ii)).^2,[1 length(Ev)]);
%         u2=repmat(abs(f(pos,:)*V(1:n,ii)).^2      +abs(f(pos,:)*V(n+1:2*n,ii)).^2,[1 length(Ev)]);
%         v2=repmat(abs(f(pos,:)*V(2*n+1:3*n,ii)).^2+abs(f(pos,:)*V(3*n+1:4*n,ii)).^2,[1 length(Ev)]);
%         DErun=DErun+2*GAMMA*(GAMMA*(u2+v2)+Dynes).*(u2./((repmat(Ev,[length(pos),1])-D(ii)).^2+(GAMMA*(u2+v2)+Dynes).^2));
        DErun=DErun-2*GAMMA*imag(u2./(repmat(Ev,[length(pos),1])-D(ii)+1i*(GAMMA*(u2+v2)+Dynes)));
    end
    DE(:,:,run)=DErun;
end
DE=sum(DE,3)/size(DE,3);
for ii=1:length(pos)
    figure(ii),plot(Ev,DE(ii,:))
end
%% 1D Conductance by Fisher-Lee vs E
  N=[200];
  t=10; t_couple=1; t_lead=5;zplus=1e-2i; alp=[0;2;0];
  D=[1];
 mu=[1];  shift_mu=1;
  B=[2;0;0];
  muv=[];Bv=[];Dv=[];
for ii=1:length(N)
    muv=[muv mu(ii)*ones(1,N(ii))];
     Bv=[Bv B(:,ii)*ones(1,N(ii))];
     Dv=[Dv D(:,ii)*ones(1,N(ii))];
end
N=sum(N);z=zeros(N);
H=H_1Ds(-t,muv,alp,Bv,Dv, 0,shift_mu);
Ev=linspace(-2,2,201);
Tv=zeros(size(Ev));
parfor Ei=1:length(Ev)
    E=Ev(Ei)
    rad=4*t_lead^2-E^2;
    if (rad>0.0)
        selfenergy=(E-1i*sqrt(rad)) *t_couple^2/(2*t_lead^2);
    else
        selfenergy=(E-sign(E)*sqrt(-rad)) *t_couple^2/(2*t_lead^2);
    end
    Sig=zeros(size(H));
    Sig(1,1)=1;Sig(N,N)=1;Sig(N+1,N+1)=1;Sig(2*N,2*N)=1;
    Sig(2*N+1,2*N+1)=1;Sig(3*N,3*N)=1;Sig(3*N+1,3*N+1)=1;Sig(4*N,4*N)=1;
    G=inv((E+zplus)*eye(size(H))-H-Sig*selfenergy);
    vlead=-2*imag(selfenergy);
    t=i*sqrt(vlead)*G([N,2*N],[1,N+1])*sqrt(vlead);
    Tv(Ei)=trace(t'*t);
end
figure(1),h=plot(Ev,Tv,'b-');
%% Josephson by threading flux (PRL.105,077001)
  N=[1 50 1];
  t=2;alp=[0; 1; 0];
  D=[0 1 0];
 mu=[1 0 1];  shift_mu=1;
  B=[1.25 1.25 1.25];
  muv=[];Bv=[];Dv=[];
for ii=1:length(N)
      z=zeros(1,N(ii)); 
    muv=[muv mu(ii)*ones(1,N(ii))];
     Bv=[Bv [B(ii)*ones(1,N(ii)); z; z]];
     Dv=[Dv D(ii)*ones(1,N(ii))];
end
N=sum(N);

kv=linspace(0,2*pi,51);neig=4;Ev=zeros(neig,length(kv)); 
parfor ki=1:length(kv)
    k=kv(ki);
    H=H_1Ds(-t,muv,alp,Bv,Dv,[exp(1i*k/2) exp(-1i*k/2)],1);
    Ev(:,ki)=sort(eigs(sparse(H),neig,'sm'));
end
figure(1),plot(kv,Ev,'black'); axis tight
%% Josephson by Supercondutors' phases on OBC or PBC (PRL.105,077001)
  N=[25 0 25];  bc=1; n=25; %number of bands in Fourier
  t=2;
  alp=[0 0 0; 1 1 -1; 0 0 0];
  D=[1 0 1];
 mu=[0 1 0];  shift_mu=1;
  B=[1.25 1.25 1.25];
  muv=[];Bv=[];Dv=[];alpv=[];
for ii=1:length(N)
      z=zeros(1,N(ii)); 
    muv=[muv mu(ii)*ones(1,N(ii))];
     Bv=[Bv [B(ii)*ones(1,N(ii)); z; z]];
     Dv=[Dv D(ii)*ones(1,N(ii))];
   alpv=[alpv alp(:,ii)*ones(1,N(ii))];
end
N=sum(N);

Ev=[];kv=linspace(0,2*pi,51);
for k=kv
    Dphase=[ones(1,N/2)  exp(1i*k)*ones(1,N/2)];
%     H=H_1Ds_Fourier(-t,muv,alp,Bv,Dv.*Dphase,bc,shift_mu,n); %could add _Fourier
    H=H_1Ds(-t,muv,alpv,Bv,Dv.*Dphase,bc,shift_mu);
    Ev=[Ev,sort(eigs(sparse(H),4,'sm'))];
end
figure(2),plot(kv,Ev,'black'); axis tight
%% 1D Pfaffian (arXiv:1102.3440)
  N=[200]; n=74;
  t=10;  alp=[0;2;0];
  D=[3];
 mu=[4];  shiftmu=1;
for B=[4.9,5.1]
  muv=[];Bv=[];Dv=[];
for ii=1:length(N)
      z=zeros(1,N(ii)); 
    muv=[muv mu(ii)*ones(1,N(ii))];
     Bv=[Bv [B(ii)*ones(1,N(ii)); z; z]];
     Dv=[Dv D(ii)*ones(1,N(ii))];
end
H=H_1Ds(-t,muv,alp,Bv,Dv, 0,shiftmu); disp('smallest eigenvalue in open wire:');disp(eigs(H,1,'sm'))
MB=kron(sparse([1 0 0 -1;0 1 1 0;i 0 0 i;0 i -i 0]/sqrt(2)),speye(n));
H=H_1Ds_Fourier(-t,muv,alp,Bv,Dv, 1,shiftmu,n);  A0=MB*H*(MB')/(i/4); 
H=H_1Ds_Fourier(-t,muv,alp,Bv,Dv,-1,shiftmu,n);  Api=MB*H*(MB')/(i/4); 
disp('Q='); disp(sign(pfaffian_hessenberg(full(A0))*pfaffian_hessenberg(full(Api)))) %opposite sign => Topological
end
% tic,pfaffian_LTL(A),toc %wrong answer
% tic,pfaffian_householder(A),toc %13.738042 seconds
% tic,pfaffian_hessenberg(A),toc %3.387200
%% 1D Reflection Matrix Det from Transfer Matrix (arXiv:1108.0419, PRB.78.045118App)
  N=[201];
  t=1;  alp=[0;2;0];
  D=[1];
 mu=[1];  shift_mu=1;
  B=[1.41;0;0];
  muv=[];Bv=[];Dv=[];
for ii=1:length(N)
    muv=[muv mu(ii)*ones(1,N(ii))];
     Bv=[Bv B(:,ii)*ones(1,N(ii))];
     Dv=[Dv D(:,ii)*ones(1,N(ii))];
end
H=H_1Ds(-t,muv,alp,Bv,Dv, 0,shift_mu);
UH=spalloc(4*N,4*N,4*N);for ii=0:N-1,UH(4*ii+1,ii+1)=1;UH(4*ii+2,N+ii+1)=1;UH(4*ii+3,2*N+ii+1)=1;UH(4*ii+4,3*N+ii+1)=1;end; H2=UH*H*UH'; %(site)x(ph)x(spin)
U=[zeros(4) eye(4); eye(4) zeros(4)]; step=10;
for ii=1:step:N-step
    M=eye(8);
    for n=ii:ii+step-1
        tn=[-t alp(2)/2 0 0;
            -alp(2)/2 -t 0 0;
            0 0 t -alp(2)/2;
            0 0 alp(2)/2 t];
        hn=[2*t-muv(n) Bv(1,n) Dv(n) 0;
            Bv(1,n) 2*t-muv(n) 0 Dv(n);
            Dv(n) 0 muv(n)-2*t Bv(1,n);
            0 Dv(n) Bv(1,n) muv(n)-2*t];
%         tn=H(4*(n-1)+1:4*(n-1)+4,4*n+1:4*n+4);
%         hn=H(4*(n-1)+1:4*(n-1)+4,4*(n-1)+1:4*(n-1)+4);
        M11=spalloc(4,4,0);M12=tn;M21=-inv(tn');M22=-inv(tn')*hn;
        Mn=[M11 M12;M21 M22];%     Sigy=kron([0 -i;i 0],eye(4)); Sigy*M'*Sigy-inv(M)  %current conservation
        M=Mn*M;
    end
    R=kron([1 1;i -i]/sqrt(2),eye(4)); Y=inv(R)*M*R; %Rotate to have sigz-type conservation PRB.78.045118
    Y11=Y(1:4,1:4);Y12=Y(1:4,5:8);Y21=Y(5:8,1:4);Y22=Y(5:8,5:8);
    Un=[-inv(Y22)*Y21 inv(Y22); Y11-Y12*inv(Y22)*Y21 Y12*inv(Y22)];
    a1=Un(1:4,1:4);b1=Un(1:4,5:8);c1=Un(5:8,1:4);d1=Un(5:8,5:8);
    a2=U(1:4,1:4); b2=U(1:4,5:8); c2=U(5:8,1:4); d2=U(5:8,5:8);
    U(1:4,1:4)=a1+b1*inv(eye(4)-a2*d1)*a2*c1;
    U(1:4,5:8)=b1*inv(eye(4)-a2*d1)*b2;
    U(5:8,1:4)=c2*inv(eye(4)-d1*a2)*c1;
    U(5:8,5:8)=d2+c2*inv(eye(4)-d1*a2)*d1*b2;%    U'*U %always unitary
end
R=U(1:4,1:4); det(R)
%% 1D Fidelity across a QPT point (arXiv:1302.5492)
  N=[100];
  t=10;  alp=[0;2;0];
 Ds=[3];
 mu=[4];  shiftmu=1;
 Bv=linspace(4.5,5.5,101);  dB=0.01;
 Fv=[];
for B=Bv
    H1=H_1Ds(-t*ones(1,N),mu,alp,[B+dB/2;0;0],Ds,0,shiftmu);
    H2=H_1Ds(-t*ones(1,N),mu,alp,[B-dB/2;0;0],Ds,0,shiftmu);
    [V1,D1]=eig(full(H1)); V1=V1(:,diag(D1)'<-0.001);
    [V2,D2]=eig(full(H2)); V2=V2(:,1:size(V1,2));
    f=V1'*V2;
    Fv=[Fv abs(det(f))];
end
plot(Bv,Fv)
%% 1D Berry Curvature ? (arXiv:1207.1104)
  N=10; dphi=0.1; dth=0.00001;dv=0.00001;
  t=0.5;  alp=[0;1;0]; beta=1;
  D=[0.7];
 mu=[0];  shiftmu=0;
  B=1;
  muv=[];Bv=[];Dv=[];
 vv=linspace(-pi,pi,51);
 Fv=[];
for ii=1:length(N)
      z=zeros(1,N(ii)); 
    muv=[muv mu(ii)*ones(1,N(ii))];
     Bv=[Bv [B(ii)*ones(1,N(ii)); z; z]];
     Dv=[Dv D(ii)*ones(1,N(ii))];
end
taux=kron(sparse(kron([0 -i;i 0],[0 -i;i 0])),speye(N));
syty=-kron(sparse(kron([0 1;1 0],speye(2))),speye(N));
for v=vv
    Hphi =H_1Ds(-t,muv,alp,Bv,Dv,exp(+i*dphi),shiftmu)+sin(v)*taux+beta*(2-2*cos(v))*syty;
    Hmphi=H_1Ds(-t,muv,alp,Bv,Dv,exp(-i*dphi),shiftmu)+sin(v)*taux+beta*(2-2*cos(v))*syty;
    Hth =H_1Ds(-t,muv,alp,Bv,Dv,1,shiftmu)+sin(v+dv)*taux+beta*(2-2*cos(v+dv))*syty;
    Hmth=H_1Ds(-t,muv,alp,Bv,Dv,1,shiftmu)+sin(v-dv)*taux+beta*(2-2*cos(v-dv))*syty;
    [V,~]=eig(full(Hphi));  psiphi1=V(:,1:end/2);
    [V,~]=eig(full(Hmphi)); psiphi2=V(:,1:end/2);
    [V,~]=eig(full(Hth));   psith1=V(:,1:end/2);
    [V,~]=eig(full(Hmth));  psith2=V(:,1:end/2);
    overlap=det(psith1'*psiphi1)-det(psith2'*psiphi1)-det(psith1'*psiphi2)+det(psith2'*psiphi2);
    Fv=[Fv imag(overlap)];
end
plot(vv,Fv)
%% 1D Pfaffian vs. Reflection Matrix (arXiv:1102.3440)
N=2001; N_for_Pf=20;
t=3;alp=[0;1;0];mur=1.2*randn(1,N);del=1;Vzee=[sqrt(2)+.45;0;0];
mu0v=linspace(0,25,101);Qv=[];Qpfv=[];
for mu0=mu0v
    mu=mu0+mur;
    H=H_1Ds(-t,mu,alp,Vzee,del, 0,shiftmu);
    UH=spalloc(4*N,4*N,4*N);for ii=0:N-1,UH(4*ii+1,ii+1)=1;UH(4*ii+2,N+ii+1)=1;UH(4*ii+3,2*N+ii+1)=1;UH(4*ii+4,3*N+ii+1)=1;end; H=UH*H*UH'; %(site)x(ph)x(spin)
    U=[zeros(4) eye(4); eye(4) zeros(4)]; step=10;
    for ii=1:step:N-step
        M=eye(8);
        for n=ii:ii+step-1
            tn=H(4*(n-1)+1:4*(n-1)+4,4*n+1:4*n+4);
            hn=H(4*(n-1)+1:4*(n-1)+4,4*(n-1)+1:4*(n-1)+4);
            M11=spalloc(4,4,0);M12=tn;M21=-inv(tn');M22=-inv(tn')*hn;
            Mn=[M11 M12;M21 M22];%     Sigy=kron([0 -i;i 0],eye(4)); Sigy*M'*Sigy-inv(M)  %current conservation
            M=Mn*M;
        end
        R=kron([1 1;i -i]/sqrt(2),eye(4)); Y=inv(R)*M*R; %Rotate to have sigz-type conservation PRB.78.045118
        Y11=Y(1:4,1:4);Y12=Y(1:4,5:8);Y21=Y(5:8,1:4);Y22=Y(5:8,5:8);
        Un=[-inv(Y22)*Y21 inv(Y22); Y11-Y12*inv(Y22)*Y21 Y12*inv(Y22)];
        a1=Un(1:4,1:4);b1=Un(1:4,5:8);c1=Un(5:8,1:4);d1=Un(5:8,5:8);
        a2=U(1:4,1:4); b2=U(1:4,5:8); c2=U(5:8,1:4); d2=U(5:8,5:8);
        U(1:4,1:4)=a1+b1*inv(eye(4)-a2*d1)*a2*c1;
        U(1:4,5:8)=b1*inv(eye(4)-a2*d1)*b2;
        U(5:8,1:4)=c2*inv(eye(4)-d1*a2)*c1;
        U(5:8,5:8)=d2+c2*inv(eye(4)-d1*a2)*d1*b2;%    U'*U %always unitary
        if abs(abs(det(U(1:4,1:4)))-1)<1e-10
            break
        end
    end
    Qv=[Qv real(det(U(1:4,1:4)))];
    MB=kron(sparse([1 0 0 -1;0 1 1 0;i 0 0 i;0 i -i 0]/sqrt(2)),speye(N_for_Pf));
    H=H_1Ds(-t,mu0+mur(1:N_for_Pf),alp,Vzee,del, 1,shiftmu);  A0=MB*H*(MB')/(i/4); 
    H=H_1Ds(-t,mu0+mur(1:N_for_Pf),alp,Vzee,del,-1,shiftmu);  Api=MB*H*(MB')/(i/4); 
    Qpfv=[Qpfv sign(pfaffian_householder(full(A0))*pfaffian_householder(full(Api)))];
end
plot(mu0v,[Qv;Qpfv])

%% 1D OBC Excitation E vs Parameter
parv=linspace(-8,8,96);neig=10;
Ev=zeros(neig,length(parv));opts.tol=5E-3;
parfor pari=1:length(parv)
    par=parv(pari);
%     Nv=ones(1,1000);N=sum(Nv);n=200;t=1;Ds=1*exp(i*(1:N)*par);u=sqrt(4*1/2*t);alp=[0;u;0];mu=0*Nv;shift_mu=1;B=.9*Nv; %arXiv:1110.6193 Phase Gradient
    %Nv=2000;t=12;Ds=[0;1.];alp=[0;0;4];mu=par;shift_mu=0;B=0; %arXiv:1211.0338 d-wave
    Nv=600;t=6;Ds=[0;2];alp=[0;4;0];mu=par;shift_mu=0;B=[0;0;0]; %arXiv:1110.4575 d-wave
    N=sum(Nv);muv=[];Bv=[];Dv=[];
    for ii=1:length(Nv)
          z=zeros(1,Nv(ii)); 
        muv=[muv mu(ii)*ones(1,Nv(ii))];
         Bv=[Bv  B(:,ii)*ones(1,Nv(ii))];
         Dv=[Dv Ds(:,ii)*ones(1,Nv(ii))];
    end
    H=H_1Ds(-t,muv,alp,Bv,Dv,0,shift_mu);
%     H=H_1Ds_Fourier(-t,muv,alp,Bv,Dv,0,shift_mu,n);
%    D=sort(eigs(H,2*neig,0,opts));D2=D(end/2+1:end);
    D=eig(full(H));D2=D(end/2+1:end/2+neig);
    Ev(:,pari)=D2;
end
plot(parv,Ev,'black'),axis([min(parv) max(parv) 0 0.7]),axis tight
%% 1D Floquet Quasienergies (PRL.106.220402)
N=40; %70
Ds=2; mu1=-1; mu2=-3;  Bx=1; ka=pi/8;
t0=-cos(ka);
SO_z=sin(ka)*2;
  t=t0*ones(1,N);
H1=H_1Ds(t,mu1,[SO_z;0;0],[0;0;Bx],Ds,0,0);
H2=H_1Ds(t,mu2,[SO_z;0;0],[0;0;Bx],Ds,0,0);

Tv=0.3:0.02:1.6;  Ev=zeros(length(H1),length(Tv));
for Ti=1:length(Tv)
    T=Tv(Ti);
    U=expm(-i*H2*T/2)*expm(-i*H1*T/2);
    Ev(:,Ti)=sort(imag(log(eig(U))));
end
figure(1),subplot(4,1,[2 3 4]),plot(Tv,Ev,'black','LineWidth',0.01),axis([0.2 1.6 -3.2 3.2]);
%%  . Topological Charges (TODO)
Tv=0.2:0.001:1.6;
CloseReal=zeros(4,length(Tv)); %any roots close to real axis
for Ti=1:length(Tv)
    T=Tv(Ti);
    H10=[-2*cos(ka)-mu1 Bx Ds 0;
         Bx -2*cos(ka)-mu1 0 Ds;
         Ds 0 +2*cos(ka)+mu1 Bx;
         0 Ds Bx +2*cos(ka)+mu1];
    H20=[-2*cos(ka)-mu2 Bx Ds 0;
         Bx -2*cos(ka)-mu2 0 Ds;
         Ds 0 +2*cos(ka)+mu2 Bx;
         0 Ds Bx +2*cos(ka)+mu2];
    H1pi=[+2*cos(ka)-mu1 Bx Ds 0;
          Bx +2*cos(ka)-mu1 0 Ds;
          Ds 0 -2*cos(ka)+mu1 Bx;
          0 Ds Bx -2*cos(ka)+mu1];
    H2pi=[+2*cos(ka)-mu2 Bx Ds 0;
          Bx +2*cos(ka)-mu2 0 Ds;
          Ds 0 -2*cos(ka)+mu2 Bx;
          0 Ds Bx -2*cos(ka)+mu2];
    U0=expm(-i*H20*T/2)*expm(-i*H10*T/2);
    Upi=expm(-i*H2pi*T/2)*expm(-i*H1pi*T/2);
%    M0=logm(U0);    Mpi=logm(Upi);
%    Q0pi(Ti)=pfaff(real(i*M0))*pfaff(real(i*Mpi)); %not real?
    if any(abs(imag(eig(U0)))<0.005)
        CloseReal(1,Ti)=1;
    end
    if any(abs(imag(eig(Upi)))<0.005)
        CloseReal(2,Ti)=2;
    end
    if any(abs(imag(sqrt(eig(U0))))<0.005)
        CloseReal(3,Ti)=3;
    end
    if any(abs(imag(sqrt(eig(Upi))))<0.005)
        CloseReal(4,Ti)=4;
    end
end
figure(1),subplot(4,1,1),plot(Tv,CloseReal),axis([0.2 1.6 0 4.2]);
%% Floquet on SO
%% 1D BTK N-S junction (Blonder-Tinkham-Klapwijk: coherent tunneling - reflection only)
NN=100;NS=200;  N=NN+NS;
a0=0.1;
del_g=0;del_sc=0.5;  mu_g=5;mu_sc=5;    Ev=0:0.02:1.1*del_sc;
  t_g=5;t_sc=t_g;   alp_g=0.0/a0; alp_sc=alp_g;
del=[del_g*ones(1,NN) del_sc*ones(1,NS)];
 mu=[ mu_g*ones(1,NN)  mu_sc*ones(1,NS)]; mu(NN)=1*mu(NN);
Vzee=[   0*ones(1,NN)      0*ones(1,NS);         0*ones(1,NN)      0*ones(1,NS);         0*ones(1,NN)      0*ones(1,NS)];
  t=[  t_g*ones(1,NN)   t_sc*ones(1,NS)];
alp=[zeros(1,N);
     alp_g*ones(1,NN) alp_sc*ones(1,NS);
     zeros(1,N)];
H=H_1Ds(-t,mu,alp,Vzee,del,0,1);
b=zeros(4*N+4,2);con=[];
for E=Ev
    ke=acos(( E+mu_g-2*t_g)/(-2*t_g))/a0;
    kh=acos((-E+mu_g-2*t_g)/(-2*t_g))/a0; %??
    M=[H-E*eye(4*N) zeros(4*N,4)];
    M(3*N+1,:)=[]; M(2*N+1,:)=[]; M(N+1,:)=[]; M(1,:)=[];
    M=[M; 1             zeros(1,4*N-1) -1  0  0  0]; b(4*(N-1)+1,:)=[1 0];%1p,up 
    M=[M; zeros(1,N)   1 zeros(1,3*N-1) 0 -1  0  0]; b(4*(N-1)+2,:)=[0 1];%1p,down
    M=[M; zeros(1,2*N) 1 zeros(1,2*N-1) 0  0 -1  0]; %1hu
    M=[M; zeros(1,3*N) 1 zeros(1,N-1)   0  0  0 -1]; %1hd

    M=[M; 0             1 zeros(1,4*N-2) -exp(-i*ke*a0)  0  0  0]; b(4*N+1,:)=[exp(i*ke*a0) 0];%2pu
    M=[M; zeros(1,N+1)   1 zeros(1,3*N-2) 0 -exp(-i*ke*a0)  0  0]; b(4*N+2,:)=[0 exp(i*ke*a0)];%2pd
    M=[M; zeros(1,2*N+1) 1 zeros(1,2*N-2) 0  0 -exp(+i*ke*a0)  0];%2hu
    M=[M; zeros(1,3*N+1) 1 zeros(1,N-2)   0  0  0 -exp(+i*ke*a0)];%2hd
    s=M\b;
    con=[con,abs(s(4*N+3,1))^2+abs(s(4*N+4,1))^2+abs(s(4*N+3,2))^2+abs(s(4*N+4,2))^2]
end
plot(Ev,con)
%% 1D BTK N-S-N junction (coherent tunneling)
NN=50;NS=50;  N=NN+NS;NN2=NN/2;
a0=0.1;
del_g=0;del_sc=0.5;  mu_g=5;mu_sc=5;    Ev=0:0.1:5*del_sc;
  t_g=5;t_sc=t_g;   alp_g=0.0/a0; alp_sc=alp_g;
del=[del_g*ones(1,NN2) del_sc*ones(1,NS) del_g*ones(1,NN2)];
 mu=[ mu_g*ones(1,NN2)  mu_sc*ones(1,NS)  mu_g*ones(1,NN2)]; mu(NN2)=1*mu(NN2);
Vzee=[ 0*ones(1,NN2) 0*ones(1,NS) 0*ones(1,NN2);
       0*ones(1,NN2) 0*ones(1,NS) 0*ones(1,NN2);
       0*ones(1,NN2) 0*ones(1,NS) 0*ones(1,NN2)];
  t=[  t_g*ones(1,NN2)   t_sc*ones(1,NS) t_g*ones(1,NN2)];
alp=[zeros(1,N);
     alp_g*ones(1,NN2) alp_sc*ones(1,NS) alp_g*ones(1,NN2);
     zeros(1,N)];
H=H_1Ds(-t,mu,alp,Vzee,del,0,1);
b=zeros(4*N+8,2);con=[];
for E=Ev
    ke=acos(( E+mu_g-2*t_g)/(-2*t_g))/a0;
    kh=acos((-E+mu_g-2*t_g)/(-2*t_g))/a0; %??
    kesc=acos(( sqrt(E^2+del_sc^2)+mu_g-2*t_g)/(-2*t_g))/a0;
    khsc=acos(( sqrt(E^2+del_sc^2)+mu_g-2*t_g)/(-2*t_g))/a0;
    M=[H-E*eye(4*N) zeros(4*N,8)];
    M(4*N,:)=[];M(3*N+1,:)=[];M(3*N,:)=[]; M(2*N+1,:)=[];M(2*N,:)=[]; M(N+1,:)=[];M(N,:)=[]; M(1,:)=[];
    M=[M; 1             zeros(1,4*N-1) -1  0  0  0  0 0 0 0]; b(4*(N-1)+1,:)=[1 0];%1p,up 
    M=[M; zeros(1,N)   1 zeros(1,3*N-1) 0 -1  0  0  0 0 0 0]; b(4*(N-1)+2,:)=[0 1];%1p,down
    M=[M; zeros(1,2*N) 1 zeros(1,2*N-1) 0  0 -1  0  0 0 0 0]; %1hu
    M=[M; zeros(1,3*N) 1 zeros(1,N-1)   0  0  0 -1  0 0 0 0]; %1hd

    M=[M; 0             1 zeros(1,4*N-2) -exp(-i*ke*a0)  0  0  0  0 0 0 0]; b(4*N+1,:)=[exp(i*ke*a0) 0];%2pu
    M=[M; zeros(1,N+1)   1 zeros(1,3*N-2) 0 -exp(-i*ke*a0)  0  0  0 0 0 0]; b(4*N+2,:)=[0 exp(i*ke*a0)];%2pd
    M=[M; zeros(1,2*N+1) 1 zeros(1,2*N-2) 0  0 -exp(+i*ke*a0)  0  0 0 0 0];%2hu
    M=[M; zeros(1,3*N+1) 1 zeros(1,N-2)   0  0  0 -exp(+i*ke*a0)  0 0 0 0];%2hd
    
    M=[M; zeros(1,N-1)   1 zeros(1,3*N)  0 0 0 0  -exp(+i*ke*a0) 0  0  0]; %Npu
    M=[M; zeros(1,2*N-1) 1 zeros(1,2*N)  0 0 0 0  0 -exp(+i*ke*a0)  0  0]; %Npd
    M=[M; zeros(1,3*N-1) 1 zeros(1,N)    0 0 0 0  0  0 -exp(-i*ke*a0)  0]; %Nhu
    M=[M; zeros(1,4*N-1) 1               0 0 0 0  0 0  0  -exp(-i*ke*a0)]; %Nhd
    
    M=[M; zeros(1,N-2)   1 zeros(1,3*N+1) 0 0 0 0  -1 0  0  0];
    M=[M; zeros(1,2*N-2) 1 zeros(1,2*N+1) 0 0 0 0  0 -1  0  0];
    M=[M; zeros(1,3*N-2) 1 zeros(1,N+1)   0 0 0 0  0  0 -1  0];
    M=[M; zeros(1,4*N-2) 1 0              0 0 0 0  0  0  0 -1];
    s=M\b;
    con=[con,abs(s(4*N+3,1))^2+abs(s(4*N+4,1))^2+abs(s(4*N+3,2))^2+abs(s(4*N+4,2))^2]
    %con=[con,abs(s(4*N+5,1))^2+abs(s(4*N+6,1))^2+abs(s(4*N+5,2))^2+abs(s(4*N+6,2))^2]
end
plot(Ev,con)

%% 2D Total DOS
Ny=10;Nx=200;                spread=0.001;
t=10; Ds=1; 
SOx=[2;0;0];SOy=[0;2;0];
B=zeros(Ny,Nx,3);
B(:,:,3)=6*ones(Ny,Nx);
mu=-4.37*t;
Ev=linspace(-0.4,0.4,1001);  runs=1;
DE=zeros(length(Ev),runs); 
for run=1:runs
%     noise=sd*randn(1,N);
    H=H_2Ds_Fourier(-t,mu,SOx,SOy,B,Ds,[0 0],100);
    D=eig(full(H));
    DE(:,run)=sum(spread/pi./(spread^2+(Ev(:)*ones(1,length(D))-ones(length(Ev),1)*D(:).').^2),2);
end
figure(1), plot(Ev,mean(DE,2));
%% 2D Local DOS by Brute force
Ny=10;Nx=200;N=Nx*Ny;
ix=1;iy=[3 4];%position at which DOS desired; multiple iy indicates summation over it (in the same layer)
t=10;Ds=1;SOx=[0;4;0];SOy=[4;0;0];B=[0;0;6];mu=-32.65*ones(Ny,Nx);spread=1E-3;Ev=linspace(-0.4,0.4,201);
pause(3)
mu=rand(Ny,1)*ones(1,Nx);
H=H_2Ds(-t,mu,SOx,SOy,B,Ds,[0 0]);
DOSv=zeros(size(Ev));
[V,D]=eig(full(H)); D=diag(D);
u_u = reshape(V(1:N,:),Ny,Nx,4*N);
u_d = reshape(V(N+1:2*N,:),Ny,Nx,4*N);
v_d = reshape(V(2*N+1:3*N,:),Ny,Nx,4*N);
v_u =-reshape(V(3*N+1:4*N,:),Ny,Nx,4*N);
for ii=1:length(D)
    DOSv=DOSv+sum((abs(u_u(iy,ix,ii)).^2+abs(u_d(iy,ix,ii)).^2+abs(v_u(iy,ix,ii)).^2+abs(v_d(iy,ix,ii)).^2)*(spread/pi./(spread^2+(Ev-D(ii)).^2)));
end
figure(1),plot(Ev,DOSv/2),axis([xlim 0 1])
%% . by recursive Green (arXiv:1011.6371)
Ny=11;Nx=500;
t=10; mu=-25.6*ones(Ny,Nx);Ds=1; SOx=[0;2;0];SOy=[-2;0;0];Ds=1;B=[0;0;1.5];
[H0,H1]=H_2Ds(-t,mu,SOx,SOy,B,Ds,[0 0]); %assumes translation symmetry along x; but recursive Green does not need that
ix=1;iy=1:Ny;
zplus=.01i;
Ev=linspace(-2*Ds,2*Ds,101);
parfor Ei=1:length(Ev)
    E=Ev(Ei);
    Gs=zeros(4*Ny,4*Ny,Nx);%G within n (G_nn)
    GR=Gs; %G within n (G_nn) when all slices to the right are considered, where n ranges from Nx to 1
    GRNn=Gs; %G between slice Nx and n when all slices to the right are considered, where n ranges from Nx to 1
    GL=Gs;
    GL1n=Gs;
    Gs(:,:,1)=inv((E+zplus)*eye(4*Ny)-H0);  % Leftmost slice's G
    for n=2:Nx-1
        Gs(:,:,n)=inv((E+zplus)*eye(4*Ny)-H0);
    end
    Gs(:,:,Nx)=inv((E+zplus)*eye(4*Ny)-H0);  % Rightmost slice's G
    GL(:,:,1)=Gs(:,:,1);   % GL1n(:,:,1)=Gs(:,:,1); %GL1n and GRNn can be omitted if for LDoS
    GR(:,:,Nx)=Gs(:,:,Nx); % GRNn(:,:,Nx)=Gs(:,:,Nx);
    for n=1:ix-1 %only need to propagate FORWARD to the slice of interest (otherwise 1:Nx-1)
        GL(:,:,n+1)  =(eye(4*Ny)-Gs(:,:,n+1)*H1*GL(:,:,n)*H1')\Gs(:,:,n+1); %By Dyson's. See Metalidis's "Electronic Transport in Mesoscopic Systems"
%         GL1n(:,:,n+1)=GL1n(:,:,n)*H1*GL(:,:,n+1);
    end
    for n=Nx-1:-1:ix %only need to propagate BACKWARD to the slice of interest (otherwise Nx-1:-1:1)
        GR(:,:,n)  =(eye(4*Ny)-Gs(:,:,n)*H1'*GR(:,:,n+1)*H1)\Gs(:,:,n);
%         GRNn(:,:,n)=GRNn(:,:,n+1)*H1*GR(:,:,n);
    end
	Gnn=(eye(4*Ny)-GL(:,:,ix)*H1'*GR(:,:,ix+1)*H1)\GL(:,:,ix);
    DOSv(Ei)=-1/2/pi*imag(trace(Gnn([iy Ny+iy 2*Ny+iy 3*Ny+iy],[iy Ny+iy 2*Ny+iy 3*Ny+iy])));
end
figure(2),clf,plot(Ev,DOSv),axis([xlim 0 1])

%% 2D Ribbon Pfaffian: Contour on (B,mu)  (arXiv:1011.6371,arXiv:1111.4656)
W=5;t=10; SOx=[0;5;0];SOy=[-5;0;0];Ds=1; %arXiv:1011.6371
Bzvar=linspace(0,9,300);Bvar=[zeros(2,length(Bzvar));Bzvar]; muvar=linspace(-40,-20,300);
% W=10;t=1; SOx=[0.1*t;0;0];SOy=[0;0.1*t;0];Ds=0.1*t; %arXiv:1111.4656
% Bzvar=linspace(0,t,100);Bvar=[zeros(2,length(Bzvar));Bzvar]; muvar=linspace(0,3*t,200);
% W=10;t=1; SOx=[0.26*t;0;0]*2;SOy=[0;0.26*t;0]*2;Ds=.1*t; %arXiv:1302.2736
% Bxvar=linspace(0.2*t,0.8*t,100);Bvar=[Bxvar;zeros(size(Bxvar));0.08*t*ones(size(Bxvar))]; muvar=linspace(0.2*t,1.3*t,200)-4*t;
H1=H_1Ds_inter(-t*ones(1,W),SOy);
Q=zeros(length(muvar),length(Bvar));
parfor Bi=1:size(Bvar,2)
    B=Bvar(:,Bi);
    Qcol=zeros(length(muvar),1);
    for mui=1:length(muvar)
        mu=muvar(mui);     
        H0=H_1Ds(-t*ones(1,W),mu,SOx,B,Ds,0,0);
        MB=kron(sparse([1 0 0 -1;0 1 1 0;i 0 0 i;0 i -i 0]/sqrt(2)),speye(W));
        H=H0+H1+H1';  A0=MB*H*(MB')/(i/4); 
        H=H0-H1-H1'; Api=MB*H*(MB')/(i/4);
        Qcol(mui)=sign(pfaffian_hessenberg(full(A0))*pfaffian_hessenberg(full(Api)));
    end
    Q(:,Bi)=Qcol;
end
figure(1),contourf(sum(Bvar),muvar,Q<0),colormap bone
%% 2D Ribbon Spectrum vs k
N=30;n=0;t=ones(1,N);mu=-4;Ds=1;SOx=[1;0;0];SOy=[0;-1;0];B=[0;0;1.6];shiftmu=0;yax=[-2.9 2.9];%arXiv:0901.4693 Fig1c
% N=500;n=50;t=ones(1,N)*.5;mu=-4*mean(t);u=sqrt(6*mean(t));SOx=[u;0;0];SOy=[0;-u;0];shiftmu=0;%arXiv:1205.3185 Fig3?
% lfl=2.5;Bza=2;Dsa=1;Ds=2*Dsa*cos(lfl*u/mean(t)*(1:N)/2).^2; B=[zeros(2,N);2*Bza*sin(lfl*u/mean(t)*(1:N)/2).^2];

H0=H_1Ds_Fourier(-t,mu,SOx,B,Ds,0,shiftmu,n);
H1=H_1Ds_inter(-t(:,1:length(H0)/4),SOy);
kv=linspace(-pi,pi,51);Ev=zeros(length(H0),length(kv));
for ki=1:length(kv)
    k=kv(ki);
    H=exp(i*k)*H1+H0+exp(-i*k)*H1';
    Ev(:,ki)=eig(full(H));
end
figure(1),plot(kv,Ev,'black'); axis([-pi pi yax])
%% Transverse band bottom (arXiv:1011.6371)
W=10;
t=10; SO=2; Ds=1;  Bvar=linspace(0,8,21);
  z=zeros(1,W);
 tv=t*ones(1,W);
Dsv=Ds*ones(1,W);
SOv=[SO*ones(1,W); z; z];
H1=H_1Ds_inter(-tv,SOv);
 
bandbottom=[]; %i.e. eigenvalues at k_x=k_y=0 (considering pbc).
for Bz=Bvar
    Vzv=[z;z;Bz*ones(1,W)];
    SOv=[z; SO*ones(1,W); z];
    H0=H_1Ds(-tv,0,SOv,Vzv,Ds*0,0,0);
    H=H0+H1+H1';
    subband=sort(eig(full(H)));
    bandbottom=[bandbottom subband];
end
plot(Bvar,bandbottom,'black'),axis([0 8 -4*t -2*t]);hold on
%% . minigap in real space using pbc in x
W=10;
t=10*ones(1,W); Ds=1; 
SOx=[2;0;0];SOy=[0;2;0];
Bvar=linspace(0,8,101);
H1=H_1Ds_inter(-t,SOy);
bandbottom=[];
for Bz=Bvar
    B=[0;0;Bz];
    H0=H_1Ds(-t,0,SOx,B,0,0,0);
    H=H0+H1+H1';
    subband=sort(eig(full(H)));
    bandbottom=[bandbottom subband];
end

muvar=linspace(-40,-20,101);
Gap=zeros(length(muvar),length(Bvar));
parfor Bi=1:length(Bvar)
    Bz=Bvar(Bi);
    B=[0;0;Bz];
    MUv=zeros(length(muvar),1);
    for mui=1:length(muvar)
        mu=muvar(mui);     
        if mod(nnz(bandbottom(:,Bi)<-20 & bandbottom(:,Bi)>mu),2)==1
            H0=H_1Ds(-t,mu,SOx,B,Ds,0,0);
            MB=kron(sparse([1 0 0 -1;0 1 1 0;i 0 0 i;0 i -i 0]/sqrt(2)),speye(W));
            H=H0+H1+H1';  A0=MB*H*(MB')/(i/4); 
            H=H0-H1-H1'; Api=MB*H*(MB')/(i/4); 
            Q=sign(pfaffian_hessenberg(full(A0))*pfaffian_hessenberg(full(Api)));
            if Q==-1
                minE=10;
                for k=linspace(0,pi,101)
                    H=H0+H1*exp(i*k)+exp(-i*k)*H1';
                    D=min(abs(eig(full(H))));
                    minE=min([minE,D]);
                end
                MUv(mui)=minE;
            end
        end
    end
    Gap(:,Bi)=MUv;
end
contourf(Bvar,muvar,Gap,'LineColor','none'),colormap([1 1 1;jet]);hold on

%% 2D Reflection Matrix Det from recursive Green: contour on (B,mu)
L=200;zplus=1E-6i;
W=10;t=10;t_lead=t; t_couple=t; SOx=[0;2;0];SOy=[2;0;0];Ds=1; %arXiv:1011.6371
Bzvar=linspace(0,8,100);Bvar=[zeros(2,length(Bzvar));Bzvar]; muvar=linspace(-40,-20,100);
H1=H_1Ds_inter(-t*ones(1,W),SOy);
Q=zeros(length(muvar),length(Bvar));
parfor Bi=1:size(Bvar,2)
    B=Bvar(:,Bi);
    Qcol=zeros(length(muvar),1);
    for mui=1:length(muvar)
        mu=muvar(mui);
        H0=H_1Ds(-t*ones(1,W),mu,SOx,B,Ds,0,0);
        H0L=H_1Ds(-t*ones(1,W),mu*0,SOx,B*0,Ds*0,0,0);
        [VT,D]=eig(full(H0L));D=diag(D);
        Gs=cell(L,1);%G within n (G_nn)
        GR=Gs; %G within n (G_nn) when all slices to the right are considered, where n ranges from Nx to 1
        GRNn=GR; %G between slice Nx and n when all slices to the right are considered, where n ranges from Nx to 1
        GL=Gs;
        GL1n=Gs;
        ck=(D-0-0*zplus)./(2*t_lead);ka=acos(ck);
        sigL=diag(-t_couple*exp(i.*ka));
        sigR=diag(-t_couple*exp(i.*ka));
        Gs{1}=inv((0+zplus)*eye(size(H0))-VT'*H0*VT-sigL);  % Leftmost slice's G
        for n=2:L-1
            Gs{n}=inv((0+zplus)*eye(size(H0))-VT'*H0*VT);
        end
        Gs{L}=inv((0+zplus)*eye(size(H0))-VT'*H0*VT-sigR);  % Rightmost slice's G
        GL{1}=Gs{1};  GL1n{1}=Gs{1};
        GR{L}=Gs{L};  GRNn{L}=Gs{L};
        for n=1:L-1
            H1t=VT'*H1*VT;  %H1_{n+1<-n} in the basis VT
            GL{n+1}  =(eye(size(H0))-Gs{n+1}*H1t*GL{n}*H1t')\Gs{n+1}; %By Dyson's. See Metalidis's "Electronic Transport in Mesoscopic Systems"
            GL1n{n+1}=GL1n{n}*H1t*GL{n+1};
        end
        for n=L-1:-1:1
            H1t=VT'*H1*VT;  %H1_{n+1<-n} in the basis VT
            GR{n}  =(eye(size(H0))-Gs{n}*H1t'*GR{n+1}*H1t)\Gs{n}; %By Dyson's. See Metalidis's "Electronic Transport in Mesoscopic Systems"
            GRNn{n}=GRNn{n+1}*H1t*GR{n};
        end
        vL=-2*imag(sigL);    vR=-2*imag(sigR);
        H1t=VT'*H1*VT;
        G11=(eye(size(H0))-GL{1}*H1t'*GR{1+1}*H1t)\GL{1};
        r=i*sqrt(vL)*G11*sqrt(vL)-eye(size(H0));
        Qcol(mui)=det(r);
    end
    Q(:,Bi)=Qcol;
end
figure(3),contourf(sum(Bvar),muvar,Q<0),colormap bone
%% 2D Reflection Matrix
W=1;L=100; zplus=1E-5i;
t=25; t_lead=t; t_couple=.3*t; SOx=[0;5;0];SOy=[5;0;0];Ds=1;
B=[0;0;2]; mu=sqrt(3)-.1-2*t;
H1 =H_1Ds_inter(-t*ones(1,W),SOx);
H0 =H_1Ds(-t*ones(1,W),mu,SOy,B,Ds,0,0);
H0L=H_1Ds(-t*ones(1,W),.1-2*t,SOy*1,B*1,Ds*0,0,0);
[VT,D]=eig(full(H0L));D=diag(D);
Gs=cell(L,1);%G within n (G_nn)
GR=Gs; %G within n (G_nn) when all slices to the right are considered, where n ranges from Nx to 1
GRNn=GR; %G between slice Nx and n when all slices to the right are considered, where n ranges from Nx to 1
GL=Gs;
GL1n=Gs;
E=0;
ck=(D-E-0*zplus)./(2*t_lead);ka=acos(ck);
sigL=diag(-t_couple*exp(i.*ka));
sigR=diag(-t_couple*exp(i.*ka));
Gs{1}=inv((E+zplus)*eye(size(H0))-VT'*H0*VT-sigL);  % Leftmost slice's G
for n=2:L-1
    Gs{n}=inv((E+zplus)*eye(size(H0))-VT'*H0*VT);
end
Gs{L}=inv((E+zplus)*eye(size(H0))-VT'*H0*VT-sigR);  % Rightmost slice's G
GL{1}=Gs{1};  GL1n{1}=Gs{1};
GR{L}=Gs{L};  GRNn{L}=Gs{L};
for n=1:L-1
    H1t=VT'*H1*VT;  %H1_{n+1<-n} in the basis VT
    GL{n+1}  =(eye(size(H0))-Gs{n+1}*H1t*GL{n}*H1t')\Gs{n+1}; %By Dyson's. See Metalidis's "Electronic Transport in Mesoscopic Systems"
    GL1n{n+1}=GL1n{n}*H1t*GL{n+1};
end
for n=L-1:-1:1
    H1t=VT'*H1*VT;  %H1_{n+1<-n} in the basis VT
    GR{n}  =(eye(size(H0))-Gs{n}*H1t'*GR{n+1}*H1t)\Gs{n}; %By Dyson's. See Metalidis's "Electronic Transport in Mesoscopic Systems"
    GRNn{n}=GRNn{n+1}*H1t*GR{n};
end
vL=-2*imag(sigL);    vR=-2*imag(sigR);
H1t=VT'*H1*VT;
G11=(eye(size(H0))-GL{1}*H1t'*GR{1+1}*H1t)\GL{1};
t=i*sqrt(vL)*GRNn{1}*sqrt(vR); tee=t(1:2*W,1:2*W); teh=t(1:2*W,2*W+1:4*W);
r=i*sqrt(vL)*G11*sqrt(vL)-eye(size(H0)); ree=r(1:2*W,1:2*W); reh=r(1:2*W,2*W+1:4*W);
rhh=r(2*W+1:4*W,2*W+1:4*W); rhe=r(2*W+1:4*W,1:2*W);
[trace(ree'*ree),trace(rhe'*rhe),trace(tee'*tee),trace(teh'*teh)]
2-trace(ree'*ree)+trace(rhe'*rhe)
%% 2D Conductance from Reflection Matrix by Recursive Green??
Ly=1; Lx=150; G0=1E-3i;
t=10; t_lead=1*t; t_couple=1*t; SOx=[0;5;0];SOy=[5;0;0];Ds=1;B=[0;0;1];
Ev=linspace(-2*Ds,2*Ds,101); mu=1; shift_mu=1;
H0lead=H_1Ds(-t_lead*ones(1,Ly),mu,SOy*0,B*0,Ds*0,0,shift_mu);
H0s=H_1Ds(-t*ones(1,Ly),mu,SOy,B,Ds,0,shift_mu);
H1s=H_1Ds_inter(-t*ones(1,Ly),SOx);
NT=size(H0s,1)-0;
[VT,D]=eig(full(H0lead));VT=VT(:,1:NT);D=diag(D);D=D(1:NT);
Gs=cell(Lx,1);%G within n (G_nn)
GR=Gs; %G within n (G_nn) when all slices to the right are considered, where n ranges from Nx to 1
GRNn=GR; %G between slice Nx and n when all slices to the right are considered, where n ranges from Nx to 1
GL=Gs;
GL1n=Gs;
Tv=zeros(size(Ev));Gv=zeros(size(Ev));
for Ei=1:length(Ev)
    E=Ev(Ei);
	ck=(D-E-G0*1i)./(2*t_lead);ka=acos(ck);
	sigL=diag(-t_couple*exp(i.*ka));sigR=diag(-t_couple*exp(i.*ka));
    Gs{1}=inv((E+G0*1i)*eye(NT)-VT'*H0s*VT-sigL);  % Leftmost slice's G
    for n=2:Lx-1
        Gs{n}=inv((E+G0*1i)*eye(NT)-VT'*H0s*VT);
    end
    Gs{Lx}=inv((E+G0*1i)*eye(NT)-VT'*H0s*VT-sigR);  % Rightmost slice's G
    GL{1}=Gs{1};    GL1n{1}=Gs{1};
    GR{Lx}=Gs{Lx};  GRNn{Lx}=Gs{Lx};
    H1t=VT'*H1s*VT;  %H1_{n+1<-n} in the basis VT
    for n=1:Lx-1
        GL{n+1}  =(eye(NT)-Gs{n+1}*H1t*GL{n}*H1t')\Gs{n+1}; %By Dyson's. See Metalidis's "Electronic Transport in Mesoscopic Systems"
        GL1n{n+1}=GL1n{n}*H1t*GL{n+1};
    end
    for n=Lx-1:-1:1
        GR{n}  =(eye(NT)-Gs{n}*H1t'*GR{n+1}*H1t)\Gs{n}; %By Dyson's. See Metalidis's "Electronic Transport in Mesoscopic Systems"
        GRNn{n}=GRNn{n+1}*H1t*GR{n};
    end
	GamL=1i*(sigL-sigL');GamR=1i*(sigR-sigR');
    GamL=GamL(1:end/2,1:end/2);  GamR=GamR(1:end/2,1:end/2);
	G11=(eye(NT)-GL{1}*H1t'*GR{1+1}*H1t)\GL{1};
    g=G11(1:end/2,1:end/2); f=G11(1:end/2,end/2+1:end);
    Gv(Ei)=real(trace(-2*GamL*imag(g)-GamL*g*GamL*g'+GamL*f*GamL*f'));
end
figure(2),plot(Ev,Gv)

%% 2D Conductance from Reflection Matrix by Recursive Green: contour on (mu,E)
W=5;L=50; zplus=1E-5i;
t=25; t_lead=1*t; t_couple=.3*t; SOx=[0;5;0];SOy=[5;0;0];Ds=1;B=[0;0;2];%arXiv:1212.5879
Evar=linspace(0,0.25,60); muvar=linspace(-96,-86,120);
C=zeros(length(Evar),length(muvar));

parfor mui=1:length(muvar)
    mu=muvar(mui);
    Ccol=zeros(length(Evar),1);
    H1 =H_1Ds_inter(-t*ones(1,W),SOx);
    H0 =H_1Ds(-t*ones(1,W),mu,SOy,B,Ds,0,0);
    H0L=H_1Ds(-t_lead*ones(1,W),mu*1,SOy*1,B*1,Ds*0,0,0);
    [VT,D]=eig(full(H0L));D=diag(D);
    Gs=cell(L,1);%G within n (G_nn)
    GR=Gs; %G within n (G_nn) when all slices to the right are considered, where n ranges from Nx to 1
    GRNn=GR; %G between slice Nx and n when all slices to the right are considered, where n ranges from Nx to 1
    GL=Gs;
    GL1n=Gs;
    for Ei=1:length(Evar)
        E=Evar(Ei);
        ck=(D-E-0*zplus)./(2*t_lead);ka=acos(ck);
        sigL=diag(-t_couple*exp(i.*ka));
        sigR=diag(-t_couple*exp(i.*ka));
        Gs{1}=inv((E+zplus)*eye(size(H0))-VT'*H0*VT-sigL);  % Leftmost slice's G
        for n=2:L-1
            Gs{n}=inv((E+zplus)*eye(size(H0))-VT'*H0*VT);
        end
        Gs{L}=inv((E+zplus)*eye(size(H0))-VT'*H0*VT-sigR);  % Rightmost slice's G
        GL{1}=Gs{1};  GL1n{1}=Gs{1};
        GR{L}=Gs{L};  GRNn{L}=Gs{L};
        for n=1:L-1
            H1t=VT'*H1*VT;  %H1_{n+1<-n} in the basis VT
            GL{n+1}  =(eye(size(H0))-Gs{n+1}*H1t*GL{n}*H1t')\Gs{n+1}; %By Dyson's. See Metalidis's "Electronic Transport in Mesoscopic Systems"
            GL1n{n+1}=GL1n{n}*H1t*GL{n+1};
        end
        for n=L-1:-1:1
            H1t=VT'*H1*VT;  %H1_{n+1<-n} in the basis VT
            GR{n}  =(eye(size(H0))-Gs{n}*H1t'*GR{n+1}*H1t)\Gs{n}; %By Dyson's. See Metalidis's "Electronic Transport in Mesoscopic Systems"
            GRNn{n}=GRNn{n+1}*H1t*GR{n};
        end
        vL=-2*imag(sigL);    vR=-2*imag(sigR);
        H1t=VT'*H1*VT;
        G11=(eye(size(H0))-GL{1}*H1t'*GR{1+1}*H1t)\GL{1};
        T=i*sqrt(vL)*GRNn{1}*sqrt(vR); tee=T(1:2*W,1:2*W); teh=T(1:2*W,2*W+1:4*W);
        R=i*sqrt(vL)*G11*sqrt(vL)-eye(size(H0)); ree=R(1:2*W,1:2*W); reh=R(1:2*W,2*W+1:4*W);
        rhh=R(2*W+1:4*W,2*W+1:4*W); rhe=R(2*W+1:4*W,1:2*W);
        Ccol(Ei)=2*W-trace(ree'*ree-rhe'*rhe+tee'*tee-teh'*teh);
    end
    C(:,mui)=Ccol;
end
figure(3),contourf(muvar,Evar,C),colormap jet,shading flat
%% 2D OBC Eigen vs. Parameter
% neig=12;n=800;pvar=linspace(0,27,30);%arXiv:1202.5057
neig=6;n=500;parv=linspace(-96,-86,120);%arXiv:1212.5879
Ev=zeros(neig,length(parv));
tic,parfor pari=1:length(parv)
%     Ly=30; Lx=600;t=160;Ds=1;mu=-4*t;u=40;SOx=[0;1;0]*u;SOy=[0;0;0]*u; B=zeros(Ly,Lx,3);B(:,:,3)=pvar(pari)*ones(Ly,Lx);%arXiv:1202.5057
    Ly=5;Lx=10*Ly;Ds=1;mu=parv(pari);t=25;SOx=[0;5;0];SOy=[5;0;0];B=zeros(Ly,Lx,3);B(:,:,1)=2;
    H=H_2Ds_Fourier(-t,mu,SOx,SOy,B,Ds,[0 0],n);
%     D=eig(full(H)); D=D(end/2-neig/2+1:end/2+neig/2);
    D=sort(eigs(H,neig,0));
    Ev(:,pari)=D;
end,toc
figure(1),plot(parv,Ev,'black'), axis tight, axis([xlim 0 0.5])
%% 2D wavefunction of (near) zero energy bound states (Fig2b/3c)
Ly=30; Lx=600;t=160;Ds=1;mu=-4*t;u=40;SOx=[0;1;0]*u; SOy=[0;0;0]*u; %arXiv:1202.5057
n=400; neig=3;
B=zeros(Ly,Lx,3);    B(:,:,3)=ones(Ly,Lx)*(18);
[H,f]=H_2Ds_Fourier(-t,mu,SOx,SOy,B,Ds,[0 0],n);
[V,D]=eigs(H,2*neig,0);diag(D)
u_u = f*V(1:n,:);
u_d = f*V(n+1:2*n,:);
v_u = conj(f)*V(2*n+1:3*n,:);
v_d = conj(f)*V(3*n+1:4*n,:);
wf=reshape(sqrt(abs(u_u).^2+abs(u_d).^2+abs(v_u).^2+abs(v_d).^2),Ly,Lx,2*neig);%shouldn't have sqrt
figure(2)
for ii=1:neig
	subplot(neig,1,ii),contourf(wf(:,:,neig-ii+1)+wf(:,:,neig+ii)),shading flat,colormap(flipud(jet))
end

%% Single Impurity Shiba Energy check
N=100;
t=-5;mu=-2;bc=1;
JSv=linspace(0,30,301);
BSEv=JSv;
parfor JSi=1:length(JSv)
    Bv=zeros(3,N);
    Bv(1,1)=JSv(JSi);
    H=H_1Ds(t,mu,[0;0;0],Bv,1,bc,0);
    D=sort(abs(eig(H)));
    BSEv(JSi)=D(1);
end
alpha=pi*1/abs(pi*2*t*sqrt(1-mu^2/4/t^2))*JSv;
figure(1),plot(JSv,[BSEv;abs(1-alpha.^2)./(1+alpha.^2)])

%% 2D Magnetic Helix (arXiv:1303.6363 Fig.2b)
Nx=100;Ny=11;N=Nx*Ny;n=N;
tneg=-2.34;
mu=2.12;
SOx=[0;0;0]; SOy=[0;0;0];
del=ones(Ny,Nx); bc=[0 0]; %y,x
B=zeros(Ny,Nx,3);absB=4;
theta=2*(1:Nx)*pi/3;
B((Ny+1)/2,:,1)=absB*cos(theta);
B((Ny+1)/2,:,2)=absB*sin(theta);
neig=1;
% tic,[Hk,f]=H_2Ds_Fourier(tneg,mu,SOx,SOy,B,del,bc,n);toc
H=H_2Ds(tneg,mu,SOx,SOy,B,del,bc);
tic,[V,D]=eigs(full(H),2*neig,0);toc;
diag(D)
u_u = V(1:n,:);
u_d = V(n+1:2*n,:);
v_d = V(2*n+1:3*n,:);
v_u =-V(3*n+1:4*n,:);
wf=reshape(abs(u_u).^2+abs(u_d).^2+abs(v_u).^2+abs(v_d).^2,Ny,Nx,2*neig);
figure(2)
for ii=1:neig
	subplot(neig,1,ii),contourf(wf(:,:,neig-ii+1)+wf(:,:,neig+ii)),shading flat,colormap(flipud(jet))
end
%MB=kron(sparse([1 0 0 -1;0 1 1 0;i 0 0 i;0 i -i 0]/sqrt(2)),speye(length(Hk)/4));
%HkA=MB*Hk*MB'/(i/4);
%pfaffian_hessenberg(full(HkA))
%% Fig.2a: energies (non-Self-consistent)
Nx=96;Ny=19;N=Nx*Ny; n=800;
tneg=-2.34;
mu=2.12;
SOx=[0;0;0]; SOy=[0;0;0];
del=ones(Ny,Nx); bc=[0 0]; %y,x
theta=2*(1:Nx)*pi/3;
neig=2*n;
Bv=linspace(1,9,51); EMat=zeros(neig,length(Bv));
parfor Bi=1:length(Bv)
    absB=Bv(Bi);
    B=zeros(Ny,Nx,3);
    B((Ny+1)/2,:,1)=absB*cos(theta);
    B((Ny+1)/2,:,2)=absB*sin(theta);
    H=H_2Ds_Fourier(tneg,mu,SOx,SOy,B,del,bc,n);
    D=eig(full(H));
    EMat(:,Bi)=sort(D(D>0));
end
figure(n/100),plot(Bv,EMat,'b'),axis([xlim 0 1.2])

%% Ferromagnetic Impurities on Rashba surface (E vs B)
Nx=100;Ny=21;N=Nx*Ny; n=N;
tneg=-3;
mu=4*tneg+7;
SOx=[0;1;0]; SOy=[-1;0;0];so=2.5;
del=ones(Ny,Nx); bc=[0 0]; %y,x
Bt=zeros(Ny,Nx,3);
Bt((Ny+1)/2,:,1)=1;
Bv=linspace(0,10,220); EMat=zeros(4*N,length(Bv));
H0=full(H_2Ds(tneg,mu,so*SOx,so*SOy,0*Bt,del,bc));H0=(H0+H0')/2;
HB=full(H_2Ds(tneg,mu,so*SOx,so*SOy,Bt,del,bc))-H0;HB=(HB+HB')/2;
tic,parfor Bi=1:length(Bv)
    %H=H_2Ds(tneg,mu,so*SOx,so*SOy,Bv(Bi)*Bt,del,bc);
    D=eig(H0+Bv(Bi)*HB);
    EMat(:,Bi)=sort(D);
end,toc
figure(1),plot(Bv,EMat,'b'),axis([xlim -1.1 1.1])
%% wavefunction
Nx=100;Ny=11;N=Nx*Ny; n=N;
tneg=-3;
mu=4*tneg+7;
SOx=[0;1;0]; SOy=[-1;0;0];so=2.5;
del=ones(Ny,Nx); bc=[0 0]; %y,x
Bt=zeros(Ny,Nx,3);
Bt((Ny+1)/2,:,1)=1; B=8;
H=full(H_2Ds(tneg,mu,so*SOx,so*SOy,B*Bt,del,bc));
tic,[V,D]=eigs(full(H),2,0);toc;
diag(D)
u_u = V(1:n,:);
u_d = V(n+1:2*n,:);
v_d = V(2*n+1:3*n,:);
v_u =-V(3*n+1:4*n,:);
wf=reshape(abs(u_u).^2+abs(u_d).^2+abs(v_u).^2+abs(v_d).^2,Ny,Nx,2);
figure(2),contourf(wf(:,:,1)+wf(:,:,2)),shading flat,colormap(flipud(jet))
%% pfaffian vs (B,so)
Nx=2;Ny=39;N=Nx*Ny;
tneg=-3;
mu=4*tneg+7;
SOx=[0;1;0]; SOy=[1;0;0];
del=ones(Ny,Nx); bc=[0 0]; %y,x
Bt=zeros(Ny,Nx,3);
Bt((Ny+1)/2,:,1)=1;
sov=linspace(0,5,60);Bv=linspace(0,10,60);
Q=zeros(length(sov),length(Bv));
parfor Bi=1:length(Bv)
    B=Bv(Bi);
    Qcol=zeros(length(sov),1);
    for soi=1:length(sov)
        so=sov(soi);     
        [H0,H1]=H_2Ds(tneg,mu,so*SOx,so*SOy,B*Bt,del,bc);
        MB=kron(sparse([1 0 0 -1;0 1 1 0;i 0 0 i;0 i -i 0]/sqrt(2)),speye(Ny));
        H=H0+H1+H1';  A0=MB*H*(MB')/(i/4); 
        H=H0-H1-H1'; Api=MB*H*(MB')/(i/4);
        Qcol(soi)=sign(pfaffian_hessenberg(full(A0))*pfaffian_hessenberg(full(Api)));
    end
    Q(:,Bi)=Qcol;
end
figure(3),contourf(Bv,sov,Q<0),colormap bone
%% minigap in topological region vs (B,so)
gap=zeros(length(sov),length(Bv));
kxv=linspace(0,pi,60);
parfor Bi=1:length(Bv)
    B=Bv(Bi);
    gapcol=zeros(length(sov),1);
    for soi=1:length(sov)
        if Q(soi,Bi)<0
            so=sov(soi);     
            [H0,H1]=H_2Ds(tneg,mu,so*SOx,so*SOy,B*Bt,del,bc);
            cgap=1;
            for kx=kxv
                H=H0+H1*exp(1i*kx)+H1'*exp(-1i*kx);
                D=sort(abs(eigs(H,1,0)));
                cgap=min(cgap,D(1));
            end
            gapcol(soi)=cgap;
        end
    end
    gap(:,Bi)=gapcol;
end
figure(3),hold off,contourf(Bv,sov,gap)
%% E vs k
so=2.5;B=8;
kv=linspace(0,pi,250);Ev=zeros(4*Ny,length(kv));
[H0,H1]=H_2Ds(tneg,mu,so*SOx,so*SOy,B*Bt,del,bc);
H0=full(H0);H1=full(H1);
parfor ki=1:length(kv)
    k=kv(ki);
    H=H0/2+H1*exp(1i*k);
    Ev(:,ki)=eig(H+H');
end
figure(4),plot(kv,Ev,'black'),axis tight

%% 1D Long-Range
N=200;nut2=0.1;tp=1;D0=.0001;kFa=4;Gamma=1.2;xi=1E10;t0=1;mu0=.2;
Z=1+nut2/D0;
n0vec=(0:N-1)';n1vec=(1:N-1)';
n0vec=(0:1)';n1vec=(1)';
tnvec=-cos(n0vec*kFa).*exp(-n0vec/xi)./(n0vec*kFa);
tnvec(1)=-mu0;tnvec(2)=tnvec(2)-t0;
DSvec=sin(n0vec*kFa).*exp(-n0vec/xi)./(n0vec*kFa);
DSvec(1)=1;
SOxvec=2*tp/(1+tp^2)*exp(-n1vec/xi).*(cos(n1vec*kFa)+n1vec*kFa.*sin(n1vec*kFa))./(n1vec.^2*kFa);
SOvec=[0*SOxvec,SOxvec,0*SOxvec];
DTxvec=2*tp/(1+tp^2)*exp(-n1vec/xi).*(sin(n1vec*kFa)-n1vec*kFa.*cos(n1vec*kFa))./(n1vec.^2*kFa);
DTvec=[0*DTxvec,DTxvec,0*DTxvec];
H=H_1Ds_LR(tnvec,SOvec,[0;0;Gamma]*ones(1,N),DSvec,DTvec);
[V,D]=eigs(H,1,'sm')
%% OBC Excitation E vs Magnetic term
Gammav=linspace(1.5,2.5,51);neig=20;
Ev=zeros(neig,length(Gammav));opts.tol=5E-3;
parfor Gammai=1:length(Gammav)
    H=H_1Ds_LR(tnvec,SOvec,[0;0;Gammav(Gammai)]*ones(1,N),DSvec,DTvec);
    D=eig(full(H));D2=D(end/2+1:end/2+neig);
%     D2=eigs(H,neig,'sm');
    Ev(:,Gammai)=D2;
end
plot(Gammav,Ev,'black'),axis tight
%% OBC Excitation E vs chemical potential
muv=linspace(-7,4,51);neig=20;
Ev=zeros(neig,length(muv));opts.tol=5E-3;
parfor mui=1:length(muv)
    H=H_1Ds_LR([-muv(mui);tnvec(2:end)],SOvec,[0;0;Gamma]*ones(1,N),DSvec,DTvec);
    D=eig(full(H));D2=D(end/2+1:end/2+neig);
    Ev(:,mui)=D2;
end
plot(muv,Ev,'black'),axis tight
%% MF wavefunction
H=H_1Ds_LR([-.1;tvec(2:end)],SOvec,[0;0;1.4]*ones(1,N),DSvec,DTvec);
[V,D]=eigs(H,1,'sm');eigs(H,5,'sm')
figure(1),plot(sum(reshape(abs(V),N,4),2))
figure(2),plot(log(sum(reshape(abs(V),N,4),2)))
%% LDOS
N=200;nut2=0.1;tp=1;D0=.0001;kFa=4;Gamma=1.5;xi=1E10;t0=1;mu0=.2;
spread=0.02;
Z=1+nut2/D0;
n0vec=(0:N-1)';n1vec=(1:N-1)';
% n0vec=(0:1)';n1vec=(1)';
tnvec=-cos(n0vec*kFa).*exp(-n0vec/xi)./(n0vec*kFa);
tnvec(1)=-mu0;tnvec(2)=tnvec(2)-t0;
DSvec=sin(n0vec*kFa).*exp(-n0vec/xi)./(n0vec*kFa);
DSvec(1)=1;
SOxvec=2*tp/(1+tp^2)*exp(-n1vec/xi).*(cos(n1vec*kFa)+n1vec*kFa.*sin(n1vec*kFa))./(n1vec.^2*kFa);
SOvec=[0*SOxvec,SOxvec,0*SOxvec];
DTxvec=2*tp/(1+tp^2)*exp(-n1vec/xi).*(sin(n1vec*kFa)-n1vec*kFa.*cos(n1vec*kFa))./(n1vec.^2*kFa);
DTvec=[0*DTxvec,DTxvec,0*DTxvec];
H=H_1Ds_LR(tnvec,SOvec,[0;0;Gamma]*ones(1,N),DSvec,DTvec);
Ev=linspace(-1,1,1001); runs=1;
Pv=1:N;
DOS=zeros(length(Pv),length(Ev)); incV=0;
    [V,D]=eig(full(H)); D=diag(D);
    for ii=1:length(D)
        DOS=DOS+(abs(V(Pv,ii)).^2+abs(V(N+Pv,ii)).^2)*(spread/pi./(spread^2+(Ev-D(ii)).^2));
    end
% DOS=min(DOS,2);
figure(1),contourf(Ev,Pv,DOS),shading flat,colorbar
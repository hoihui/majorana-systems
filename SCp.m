set(0,'DefaultFigureWindowStyle','docked')
%kFa=acos(mu/2t), vF=2t*a*sqrt(1-mu^2/4/t^2), gap=2Dp*sin(kFa)=2Dp*sqrt(1-mu^2/4/t^2)
%% 1D wire DOS (arXiv:1202.1293)
%Eq.(2): Int~D(c_x*c_{x+1}+hc)
N=100;
t=1;mu=0;Dp=0.5; topological=abs(mu)<2*t
H=H_1Dp(-t*ones(1,N),0,mu,Dp,0,0);
%[V,D]=eigs(H,10,0);
[V,D]=eig(full(H));
D=diag(D)
[n,xout]=hist(D,size(V,2));
figure(1), plot(xout,n); 
% . wavefunction (edge state)
nth=0;
D(size(V,2)/2+nth)
figure(2), plot(1:N,sum(reshape(abs(V(:,size(V,2)/2+nth)).^2,N,2),2))
%% 1D Local DOS
 N=400; bc=1;
pos=1:N; spread=0.01;
t=3;mu=0;Dp=.5; topological=abs(mu)<2*t
sd=0;
Ev=linspace(-2,2,201); runs=1;
DOS=zeros(length(pos),length(Ev),runs);
muv=mu*ones(1,N);muv(1)=mu+10;
parfor run=1:runs
    noise=sd*randn(1,N); 
    H=H_1Dp(-t*ones(1,N),0,noise+muv,Dp,0,bc);
    [V,D]=eig(full(H)); D=diag(D);
    DErun=zeros(length(pos),length(Ev));
    for ii=1:length(D)
        DErun=DErun+(abs(V(pos,ii)).^2+0*abs(V(N+pos,ii)).^2)*(spread/pi./(spread^2+(Ev-D(ii)).^2));
    end
    DOS(:,:,run)=DErun;
end
DOS=mean(DOS,3);
figure(1),plot(Ev,DOS(1:2:10,:)),axis tight,axis([xlim 0 1])
figure(2),contourf(pos(1:20),Ev,min(.6,DOS(1:20,:)')),shading flat
%% . Majorana property: self-conjugacy
PH=kron([0 1;1 0],speye(N));
nth=0; %only nth=0 is majorana
V1=V(:,size(V,2)/2-nth); V2=V(:,size(V,2)/2+1+nth); 
[V1 V2 PH*conj(V1) abs(PH*conj(V1))-abs(V2)]
%% Pfaffian for Finite wire
N=100;
t=1;mu=2.2;Dp=0.5; disp(strcat('topological=',num2str(abs(mu)<2*t)));
H1=H_1Dp(-t*ones(1,N),0,mu,Dp,0,1);
H2=H_1Dp(-t*ones(1,N),0,mu,Dp,0,-1);
MB=kron(sparse([1 1;i -i]/sqrt(2)),speye(length(H1)/2));
HA1=MB*H1*MB'/i;
HA2=MB*H2*MB'/i;
disp(strcat('Q=',num2str(sign(pfaffian_hessenberg(full(HA1))*pfaffian_hessenberg(full(HA2))))));

%% 2D Ribbon (PRB.84.054532 Fig.1)
%Eq.(1) Int~c_y*c_{y+1}+ i * c_x*c_{x+1} +hc
N=30; %100
t=1;  Dp=0.5; mu=-2; %1a:-2 ; 1b:2
[H0,H1]=H_1Dp(-t*ones(1,N),-t,mu,Dp,i*Dp,0);
Ev=[];kv=linspace(-pi,pi,101);
for k=kv
    H=exp(i*k)*H1+H0+exp(-i*k)*H1';
    Ev=[Ev,sort(real(eig(full(H))))];
end
figure(1),plot(kv,Ev,'black'); axis([-pi,pi,0,max(max(Ev))])
%% Dp-mu Phase Diagram by Pfaffian for 2D Ribbon (Fig.2-3)
N=8;
t=1;
Dpvar=linspace(0,1.2,100);
muvar=linspace(-4,0,100);
Q=zeros(length(muvar),length(Dpvar));
parfor Dpi=1:length(Dpvar)
    Dp=Dpvar(Dpi);
    Qcol=zeros(length(muvar),1);
    for mui=1:length(muvar)
        mu=muvar(mui);     
        [H0,H1]=H_1Dp(-t*ones(1,N),-t,mu,Dp,i*Dp,0);
        MB=kron(sparse([1 1;i -i]/sqrt(2)),speye(length(H0)/2));
        H=H0+H1+H1';  A0=MB*H*(MB')/(i); 
        H=H0-H1-H1'; Api=MB*H*(MB')/(i); 
        Qcol(mui)=sign(pfaffian_hessenberg(full(A0))*pfaffian_hessenberg(full(Api)));
    end
    Q(:,Dpi)=Qcol;
end
figure(2),contourf(muvar,Dpvar,Q'<0),colormap bone
%% Finite 2D first 10 excitation eigenvalues (Fig.4-5)
t=1; mu=-1; Dp=0.1; %Fig.4:vary mu ; Fig.5:vary Dp
L=200;n=10;
[H0,H1]=H_1Dp(-t*ones(1,n),-t,mu,Dp,i*Dp,0);
Hod=kron(spdiags(ones(L,1),1,L,L),H1);
H=kron(speye(L),H0)+(Hod+Hod');
D=eigs(H,20,0);D=sort(D(D>0))/Dp;
scatter(1:10,D)

%% 2D finite (PRL 105,227003)
%Eq.(1) Int~c_y*c_{y+1}+ i * c_x*c_{x+1} +hc 
%treat each slab of length W as a subsystem; L copies
t=1; mu=-2; Dp=0.1;
L=200;W=100;
[H0,H1]=H_1Dp(-t*ones(1,W),-t,mu,Dp,i*Dp,0);
Hod=kron(spdiags(ones(L,1),1,L,L),H1);
H=kron(speye(L),H0)+(Hod+Hod');
[V,D]=eigs(H,4,0);D
%[V,D]=eig(full(H));min(abs(diag(D)))
v=reshape(abs(V(:,1)).^2,2*W,L);
wf1=v(1:W,:)+v(W+1:2*W,:);
v=reshape(abs(V(:,2)).^2,2*W,L);
wf2=v(1:W,:)+v(W+1:2*W,:);
figure(1),contourf(sqrt(wf1+wf2)') %shouldn't have sqrt
 shading flat,colormap(flipud(summer))
%% . E vs mu
t=1; Dp=0.1; 
L=100;W=10;neig=7;
muv=linspace(-4,0,201);
Ev=zeros(neig,length(muv));
parfor mui=1:length(muv)
    mu=muv(mui);
    [H0,H1]=H_1Dp(-t*ones(1,L),-t,mu,Dp,i*Dp,0);
    Hod=kron(spdiags(ones(W,1),1,W,W),H1);
    H=kron(speye(W),H0)+(Hod+Hod');
    D=eigs(H,neig*2,0);D=sort(D);
    Ev(:,mui)=D(end/2+1:end);
end
figure(2),plot(muv,Ev/Dp,'black'),axis([min(muv) max(muv) 0 1]);
% transverse subbands
% send L to infinity. has a ribbon structure. band bottom at k=0
mu=0;
[H0,H1]=H_1Dp(-t*ones(1,W),-t,mu,Dp,i*Dp,0);
H=H0+H1+H1';
subband=eig(full(H));subband=subband(subband<0);
for te=subband'
    figure(2),line([te ; te],[0 ; 1],'LineStyle','--')
end
%% . disorder TODO
% 
Vz=6;    Vzv=[z;z;Vz*ones(1,N)];    H0=H_1Ds(-tv,z,SOv,Vzv,Dsv,0,0);
kv=-pi:0.02:pi;E=[];
for k=kv
    H=H0+H1*exp(i*k)+H1'*exp(-i*k);
    E=[E eig(full(H))];
end
plot(kv,E)

%% Floquet Driven phase with NNN pairing (arXiv:1211.2498 Fig.2)
N=200; mu=-10; t1=1; Dp1=1; Dp2=2.5; T=0.2;
t2v=linspace(-5,5,101);
Ev=zeros(2*N,length(t2v));
parfor t2i=1:length(t2v)
    t2=t2v(t2i);
    H1=H_1Dp(-[t1;t2]*ones(1,N),0,mu,[i*Dp1;Dp2]*ones(1,N),0,0);
    H2=H_1Dp(-[t1;t2]*ones(1,N),0,mu,[Dp1;i*Dp2]*ones(1,N),0,0);
    U=expm(-i*H2*T/2)*expm(-i*H1*T/2);
    Ev(:,t2i)=sort(imag(log(eig(U))));
end
figure(1), plot(t2v,Ev,'black','LineWidth',0.01)
parfor t2i=1:length(t2v)
    t2=t2v(t2i);
    H0=H_1Dp(-[t1;t2]*ones(1,N),0,mu,[Dp1;i*Dp2]*ones(1,N),0,0);
    Ev(:,t2i)=sort(eig(full(H0)));
end
figure(2), plot(t2v,Ev,'black','LineWidth',0.01)
%% . Fig.3
set(0,'DefaultFigureWindowStyle','docked')
N=200; mu=-10; t1=1; t2=1; Dp1=1; Dp2=2.5; 
H1=H_1Dp(-[t1;t2]*ones(1,N),0,mu,[i*Dp1;Dp2]*ones(1,N),0,0);
H2=H_1Dp(-[t1;t2]*ones(1,N),0,mu,[Dp1;i*Dp2]*ones(1,N),0,0);
cmap=[linspace(0,1,50)',linspace(0,1,50)',linspace(1,1,50)';linspace(1,1,50)',linspace(1,0,50)',linspace(1,0,50)'];
T=0.2;U=expm(-i*H2*T/2)*expm(-i*H1*T/2);Heff=logm(U)/(-i*T);[r,c,s]=find(real(Heff));
figure(1),scatter(r,c,1,s); set(gca,'YDir','reverse');colormap(cmap);
T=2.0;U=expm(-i*H2*T/2)*expm(-i*H1*T/2);Heff=logm(U)/(-i*T);[r,c,s]=find(real(Heff));
figure(2),scatter(r,c,1,s); set(gca,'YDir','reverse');colormap(cmap);
[r,c,s]=find(real(H1));
figure(3),scatter(r,c,1,s); set(gca,'YDir','reverse');colormap(cmap);
[r,c,s]=find(imag(H1));
figure(4),scatter(r,c,1,s); set(gca,'YDir','reverse');colormap(cmap);

%% Finite 2D: Bulk & edge excitation vs. W (arXiv:1110.4062)
kfx=100; %Should choose L as large as possible; mu as close to -4 as possible
t=1; mu=-3.6;  kfa=sqrt(mu/t+4); Dp=t*kfa/kfx;
L=500; neig=7;
Wv=1:round(8*pi/kfa)
Ev=zeros(neig,length(Wv));
parfor Wi=1:length(Wv)
    W=Wv(Wi);
    [H0,H1]=H_1Dp(-t*ones(1,W),-t,mu,Dp,i*Dp,0);
    Hod=kron(spdiags(ones(L,1),1,L,L),H1);
    H=kron(speye(L),H0)+(Hod+Hod');
    D=sort(eigs(H,neig*2,0));
    Ev(:,Wi)=D(end/2+1:end);
end
plot(Wv*kfa/pi,Ev/Dp,'black'), axis([0,8,0,1])

%% Shiba State
N=100;
t=5;mu=-2;bc=1; Dp=.5;
muiv=exp(linspace(0,5,101));muiv=[-fliplr(muiv) muiv];
BSEv=muiv;
parfor muii=1:length(muiv)
    muv=mu*ones(1,N);
    muv(1)=muv(1)+muiv(muii);
    H=H_1Dp(-t,0,muv,Dp,0,bc);
    D=sort(eig(H));
    BSEv(muii)=D(101);
end
figure(1),plot(muiv,BSEv)
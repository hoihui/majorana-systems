function [H0,H1]=H_2Ds(tneg,mu,SOx,SOy,B,del,bc)
% basis=(ph)x(spin)x(x)
%INPUTS:
% tneg: 1x1 negative hopping matrix element
% SOx, SOy: 3x1 coefficients of spin matrices coupled to the momenta
% mu,Bx,By,Bz,DeltaS: NyxNx or 1x1 chemical potential, magnetic field, pairing
%                     where at least one must be NyxNx
% bc: first along y, then along x
% co: cutoff energy
%OUTPUTS:
% sparse H whose u_up at (x,y) for eigen mode n is given by:
%(ph)x(spin)x(contiguous x)
Ny=max([size(tneg,1) size(mu,1) size(B,1) size(del,1)]);
Nx=max([size(tneg,2) size(mu,2) size(B,2) size(del,2)]);
if nargin<7, bc=[0 0]; end
if numel(B)==3, oB=B;B=zeros(Ny,Nx,3);B(:,:,1)=oB(1);B(:,:,2)=oB(2);B(:,:,3)=oB(3);end
if numel(mu)==1, mu=mu*ones(Ny,Nx);end
if numel(del)==1, del=del*ones(Ny,Nx);end
if nargout==2 %must have translational symmetry along x
    H0=H_1Ds(tneg,mu(:,1).',SOy,[B(:,1,1)';B(:,1,2)';B(:,1,3)'],del(:,1).',bc(1),0);
    H1=H_1Ds_inter(tneg*ones(1,Ny),SOx);
else
    H0=[];
    for nx=1:Nx
        H0=blkdiag(H0,H_1Ds(tneg,mu(:,nx).',SOy,[B(:,nx,1)';B(:,nx,2)';B(:,nx,3)'],del(:,nx).',bc(1),0));
    end
    H1=H_1Ds_inter(tneg*ones(1,Ny),SOx);
    H=kron(spdiags(ones(Nx,1),1,Nx,Nx),H1);
    if bc(2)==1,H(end-Ny*4+1:end,1:Ny*4)=H1;end
    Hx=H+H'+H0;
    o=[];for ii=0:Nx-1,o=[o (1:Ny)+4*Ny*ii];end;o=[o o+Ny o+2*Ny o+3*Ny];
    H0=Hx(o,o);
end
%                                 [dx,dy,dz,dy0, err] = sddir(L,Lden,Lsd,p,...
%                                     d,v,vfrm,At,DAt,dense, R,K,y,y0,b, pars)
% SDDIR  Direction decomposition for Ye-Todd-Mizuno self-dual embedding.
%   Here, p is the direction p = dx+dz. If p=[], then assume p=-v,
%   the "affine scaling" direction.
%
% **********  INTERNAL FUNCTION OF SEDUMI **********
%
% See also sedumi

function [dx,dy,dz,dy0, err,v_2,v_4] = sddir(L,Lden,Lsd,pv,...
    d,v,vfrm,At,DAt,dense, R,K,y,y0,b, pars,pMode)
%
% This file is part of SeDuMi 1.1 by Imre Polik and Oleksandr Romanko
% Copyright (C) 2005 McMaster University, Hamilton, CANADA  (since 1.1)
%
% Copyright (C) 2001 Jos F. Sturm (up to 1.05R5)
%   Dept. Econometrics & O.R., Tilburg University, the Netherlands.
%   Supported by the Netherlands Organization for Scientific Research (NWO).
%
% Affiliation SeDuMi 1.03 and 1.04Beta (2000):
%   Dept. Quantitative Economics, Maastricht University, the Netherlands.
%
% Affiliations up to SeDuMi 1.02 (AUG1998):
%   CRL, McMaster University, Canada.
%   Supported by the Netherlands Organization for Scientific Research (NWO).
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
% 02110-1301, USA
%

% ------------------------------------------------
% dy0 = v'*pv,    ADp = A*D(d)*pv
% ------------------------------------------------
switch pMode,
    case 1,
        % Spectral values w.r.t. vfrm
        dy0 = (vfrm.lab'*pv) / R.b0;
        pv = frameit(pv,vfrm.q,vfrm.s,K);       % Let p = VFRM * p.
        v_2=[];
        v_4=[];
    case 2,
       if strcmp(pars.version,'sedumi')
%          Affine scaling
         pv=-v;
         dy0=-y0;
        v_2=[];
        v_4=[];
       else 
         v = qreshape(v,1,K);
         v_2=zeros(length(v),1);
         v_3=zeros(length(v),1);
         v_4=zeros(length(v),1);
         x_i=0;
         li=0;
         v_2(1:K.l)=v(1:K.l).^2;
         v_3(1:K.l)=v(1:K.l).^3;
         v_4(1:K.l)=v(1:K.l).^4;
         x_i=x_i+K.l;
         li=li+K.l;
         lambda=eigK(v,K);
         for i=1:length(K.q)
              lambda_1=lambda(li+1);
              lambda_2=lambda(li+2);
            if norm(v(x_i+2:x_i+K.q(i)))~=0
               frame=[1 1;-v(x_i+2:x_i+K.q(i))/norm(v(x_i+2:x_i+K.q(i))) v(x_i+2:x_i+K.q(i))/norm(v(x_i+2:x_i+K.q(i)))]/sqrt(2);
            else 
               frame=[1 1;zeros(K.q(i)-1,1) zeros(K.q(i)-1,1)]/sqrt(2); 
            end 
            v_2(x_i+1:x_i+K.q(i))=lambda_1^2*frame(:,1)+lambda_2^2*frame(:,2);
            v_3(x_i+1:x_i+K.q(i))=lambda_1^3*frame(:,1)+lambda_2^3*frame(:,2);
            v_4(x_i+1:x_i+K.q(i))=lambda_1^4*frame(:,1)+lambda_2^4*frame(:,2);
            x_i=x_i+K.q(i);
            li=li+2;
         end 
        norm_v2=norm(v_2);
        pv =-v_3/norm_v2;
        dy0 = (v'*pv) / R.b0;
        pv=qreshape(pv,0,K);
        v = qreshape(v,0,K);
        clear v_3
       end
    case 3,
        dy0 = (v'*pv) / R.b0;
        v_2=[];
        v_4=[];
    case 4
         v = qreshape(v,1,K);
         v_2=zeros(length(v),1);
         v_inv=zeros(length(v),1);
         x_i=0;
         li=0;
         v_inv(1:K.l)=1./v(1:K.l);
         v_2(1:K.l)=v(1:K.l).^2;
         x_i=x_i+K.l;
         li=li+K.l;
         lambda=eigK(v,K);
         for i=1:length(K.q)
              lambda_1=lambda(li+1);
              lambda_2=lambda(li+2);
            if norm(v(x_i+2:x_i+K.q(i)))~=0
               frame=[1 1;-v(x_i+2:x_i+K.q(i))/norm(v(x_i+2:x_i+K.q(i))) v(x_i+2:x_i+K.q(i))/norm(v(x_i+2:x_i+K.q(i)))]/sqrt(2);
            else 
               frame=[1 1;zeros(K.q(i)-1,1) zeros(K.q(i)-1,1)]/sqrt(2); 
            end 
            v_inv(x_i+1:x_i+K.q(i))=(1/lambda_1)*frame(:,1)+(1/lambda_2)*frame(:,2);
            v_2(x_i+1:x_i+K.q(i))=lambda_1^2*frame(:,1)+lambda_2^2*frame(:,2);
            x_i=x_i+K.q(i);
            li=li+2;
         end 
        mu=norm(v)^2/(K.l + 2 * length(K.q) + K.rLen + K.hLen);
        pv =mu*v_inv-v;
        dy0 = (v'*pv) / R.b0;
        pv=qreshape(pv,0,K);
        v = qreshape(v,0,K);
        v_4=[];
        clear v_inv
        
end
% ------------------------------------------------------------
% Solve  AP(d)A' * ypv = dy0 * Rb + AD(dy0*DRc - pv)
% and let xpv = (dy0*DRc - pv) - DA'*ypv.
% This yields AD*xpv = err.b-dy0*Rb.
% ------------------------------------------------------------
[dy,dx,err.kcg,err.b] = wrapPcg(L,Lden,At,dense,d, DAt,K,...
    dy0*R.b, dy0*Lsd.DRc - pv, pars.cg,min(1,y0) * R.maxRb);
% ------------------------------------------------------------
% Solve denom * rdx0 = y0*(DRc'*xpv+Rb'*ypv)-(err.b'*y)
% Here, rdx0 = deltax0/x0.
% ------------------------------------------------------------
rdx0 = (y0*(Lsd.DRc'*dx+R.b'*dy) - err.b'*y) / Lsd.denom;
% ------------------------------------------------------------
% Set dy -= rdx0 * ysd
%     dx = rdx0 * xsd - dx
% error dRb := rdx0 * Lsd.b - err.b
% ------------------------------------------------------------
dy = dy - rdx0 * Lsd.y;
dx = rdx0 * Lsd.x - dx;             % Now, AD*dx = deltax0*b + dy0*Rb...
err.b = rdx0 * Lsd.b - err.b;       % ... + err.b
err.maxb = norm(err.b,inf);
dx(1) = rdx0 * v(1);
dz = pv - dx;

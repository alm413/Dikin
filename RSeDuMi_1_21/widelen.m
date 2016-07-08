%              [t,wr,w] = widelen(xc,zc,y0, dx,dz,dy0,d2y0, maxt,pars,K,alpha)
%
% WIDELEN  Computes approximate wide-region neighborhood step length.
% Does extensive line search only if it pays, that is the resulting
% rate will be at most twice the best possible rate, and the step-length
% at least half of the best possible step length.
%
% **********  INTERNAL FUNCTION OF SEDUMI **********
%
% See also sedumi

function [t,wr,w,ind,Dual_Imp] = widelen(xc,zc,y0, dx,dz,dy0,d2y0,maxt,pars,K,v_2,v_4)
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

thetaSQR = pars.theta^2;
ix = K.mainblks;
% --------------------------------------------------
% For descent directions (dy0/y0 < 0), let fullt be
% s.t. y0 + fullt * dy0 + fullt^2 * min(d2y0,0) = 0     (NB: y0 is duality gap)
% This is used to compute the desired accuracy in t (when rate -> 0).
% --------------------------------------------------

if strcmp(pars.version,'rsedumi')
    dx = qreshape(dx,1,K);
    dz = qreshape(dz,1,K);
    xc = qreshape(xc,1,K);
    dxdz=zeros(length(dx),1);
    dxdz(1:K.l)=dx(1:K.l).*dz(1:K.l);
    x_i=K.l;
    for i=1:length(K.q)
       dxdz(x_i+1)=(dx(x_i+1:x_i+K.q(i))'*dz(x_i+1:x_i+K.q(i)))/sqrt(2);
       dxdz(x_i+2:x_i+K.q(i))=(dx(x_i+1)*dz(x_i+2:x_i+K.q(i))+dz(x_i+1)*dx(x_i+2:x_i+K.q(i)))/sqrt(2);
       x_i = x_i+K.q(i);
    end  
    tau=pars.tau;
    switch  pars.neigh
        case 1
           accept=false;
           steplength=maxt;
           rho=[0.99 .9 .8 .7 .6 .5 .4 .3 .2 .1];
           i=1;
           while  ~accept && i<=length(rho)
               lambda=eigK(v_2-(steplength*v_4/norm(v_2))+steplength^2*dxdz,K);
               kappa=max(lambda)/min(lambda);
               if abs(kappa)>=(tau+1e-12) || min(lambda)<=1e-12
                   steplength=rho(i)*maxt;
               else    
                   accept=true; 
               end   
               i=i+1;
           end  
           if accept
               Dual_Imp=abs((steplength*xc'*(dx+dz))/(norm(xc)^2));
               if Dual_Imp<=0.001
                    ind=3;
               elseif  Dual_Imp<=0.1
                    ind=0;
               else 
                    ind=1;
               end 
                    
           else  
               Dual_Imp=0;
               ind=0; 
           end  
         
         case 2
           accept=false;
           alpha_prime = 1;
           steplength = 1;
           id = eyeK(K);
           best_prox=inf;
           mu=norm(xc)^2/(K.l + 2 * length(K.q) + K.rLen + K.hLen);
           rho=[.9 .8 .7 .6 .5 .4 .3 .2 .1];
           i=1;
           while ~accept 
               lambda=eigK((1-alpha_prime)*v_2+alpha_prime*mu*id+alpha_prime^2*dxdz,K);
               prox=max(lambda)/min(lambda);
               if min(lambda)>=1e-12
                  if prox+1e-12 <= best_prox && i <=length(rho)
                      steplength=alpha_prime;
                      best_prox=prox;
                      alpha_prime=rho(i);
                  else 
                      accept=true; 
                  end 
               else 
                  alpha_prime=rho(i);
               end 
               i=i+1;
           end
           %prox_imp=best_prox-1;
%            prox_imp=(init_prox-best_prox)/init_prox;
%            if prox_imp>=0.5+1e-12
%                ind=1;
%            else
%                ind=1;
%                pars.tau=pars.tau * 2;
%            end
           Dual_Imp=0;
           ind=1;
    end    
    dx = qreshape(dx,0,K);
    dz = qreshape(dz,0,K);
    xc = qreshape(xc,0,K);
    tM=steplength;
    xM = xc + tM*dx;
    zM = zc + tM*dz;
  % LORENTZ
    wM.tdetx = tdet(xM,K);
    wM.tdetz = tdet(zM,K);
    detxz = wM.tdetx .* wM.tdetz / 4;
    if isempty(K.q)
        lab2q = zeros(0,1);
    else
        halfxz = (xM(ix(1):ix(2)-1).*zM(ix(1):ix(2)-1)...
            + ddot(xM(ix(2):ix(3)-1),zM,K.qblkstart)) / 2;
        tmp = halfxz.^2 - detxz;
        if tmp > 0
            lab2q = halfxz + sqrt(tmp);
        else
            lab2q = halfxz;
        end
    end
    % PSD
    wM.ux = psdfactor(xM,K);
    wM.s = psdscale(wM.ux,zM,K);
    % ALL:
    wM.lab = [xM(1:K.l).*zM(1:K.l); detxz ./ lab2q; lab2q; psdeig(wM.s,K)];
    [deltaM,hM,alphaM] = iswnbr(wM.lab, thetaSQR);
        w = wM;
        t = tM;
        wr.h=hM;
        wr.alpha=alphaM;
        wr.delta=deltaM;
else
    
   if dy0 < -1E-5 * y0
    if d2y0 < 0   % Concave decreasing: gap(fullt) < 0 where y0+fullt*dy0=0.
        fullt = 2 * y0 /(-dy0 + sqrt(dy0^2 - 4*y0*d2y0));  % < y0/abs(dy0).
    else
        fullt = y0 / (-dy0);
    end 
    if fullt <= 0
        error('Assertion failed: fullt <= 0')
    end 
   else 
    fullt = 2*maxt;
   end 
tR = min(maxt, fullt);
if tR < 0
    error('Assertion failed: tR >= 0')
end
% ------------------------------------------------------------
% rate = (y0+t*dy0)/y0 = 1 - (t/fullt) = (fullt - t)/fullt.
% It can never be better than (fullt - tR)/fullt.
% The next criterion ensures that we're never worse than twice best rate'
% ------------------------------------------------------------
  t = 0.0;
  ntry = 0;  % do loop at least once
  while (t < 0.5 * tR) | ( (fullt-tR) + (1e-7 * fullt) < (tR - t) ) | ntry==0
    ntry = 1;
    if tR == maxt                       % Bisection
        tM = 0.1 * t + 0.9 * tR;
    else
        tM = 0.5 * (t + tR);
    end
    xM = xc + tM*dx;
    zM = zc + tM*dz;
    % --------------------------------------------------
    % Let w = D(xM)*zM, compute lambda(w).
    % --------------------------------------------------
    % LORENTZ
    wM.tdetx = tdet(xM,K);
    wM.tdetz = tdet(zM,K);
    detxz = wM.tdetx .* wM.tdetz / 4;
    if isempty(K.q)
        lab2q = zeros(0,1);
    else
        halfxz = (xM(ix(1):ix(2)-1).*zM(ix(1):ix(2)-1)...
            + ddot(xM(ix(2):ix(3)-1),zM,K.qblkstart)) / 2;
        tmp = halfxz.^2 - detxz;
        if tmp > 0
            lab2q = halfxz + sqrt(tmp);
        else
            lab2q = halfxz;
        end
    end
    % PSD
    wM.ux = psdfactor(xM,K);
    wM.s = psdscale(wM.ux,zM,K);
    % ALL:
    wM.lab = [xM(1:K.l).*zM(1:K.l); detxz ./ lab2q; lab2q; psdeig(wM.s,K)];
    [deltaM,hM,alphaM] = iswnbr(wM.lab, thetaSQR);
    if (deltaM <= pars.beta) | ((tM < fullt / 10) & (deltaM < 1))
        w = wM;
        t = tM;
        wr.h=hM;
        wr.alpha=alphaM;
        wr.delta=deltaM;
    else
        tR = tM;
    end
  end
  ind=[];
  Dual_Imp=[];
end
% ------------------------------------------------------------
% In the nasty case that we could not make any progress, (t=0),
% take last middle value.
% ------------------------------------------------------------
if strcmp(pars.version,'sedumi')
 if t == 0
    w = wM;
    t = tM;
    wr.h=hM;
    wr.alpha=alphaM;
    wr.delta=deltaM;
 end
end
wr.desc = 1;       % always descent direction.
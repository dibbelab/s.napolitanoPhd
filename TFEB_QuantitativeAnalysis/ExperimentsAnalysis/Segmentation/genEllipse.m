%
% genEllipse(mu,cv,c,bFill)
% return the ellipse with center mu, covariance cv and isocontour c
% if bFill is 0, only returns the ellipse boundary points
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright (c) 2016, Drexel University
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
% * Redistributions of source code must retain the above copyright notice, this
%   list of conditions and the following disclaimer.
% 
% * Redistributions in binary form must reproduce the above copyright notice,
%   this list of conditions and the following disclaimer in the documentation
%   and/or other materials provided with the distribution.
% 
% * Neither the name of PixelRep nor the names of its
%   contributors may be used to endorse or promote products derived from
%   this software without specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [pts_xy,bw] = genEllipse(mu_xy,cv,mask,c,bFill)   % mask add by me (SN)

if nargin<4 %3
    c=1;
end

if nargin<5 %4
    bFill=1;
end
% bw=false(1000);
bw=false(size(mask));
    
theta=linspace(0,2*pi,1000);
[evec, eval]=eig(cv);
[eval, order] = sort(diag(eval),'descend');
evec=evec(:,order);

a1 = sqrt(eval(1))*c*evec(:,1);
a2 = sqrt(eval(end))*c*evec(:,2);

x = (mu_xy(1) + a1(1)*sin(theta) + a2(1)*cos(theta));
y = (mu_xy(2) + a1(2)*sin(theta) + a2(2)*cos(theta));
x = max(x, 1);
y = max(y, 1);
x = min(x, size(mask,2));   %1000);
y = min(y, size(mask,1));   %1000);

pts_xy=[x' y'];
idx=sub2ind(size(bw),round(y),round(x));
bw(idx)=1;

if bFill
    bw = imfill(bw,'holes');
    [r c]=find(bw);
    pts_xy=[c r];
end
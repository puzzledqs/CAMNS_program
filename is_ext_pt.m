%=====================================================================
% Programmers: 
% Wing-Kin Ma, E-mail: wkma@ieee.org
% Tsung-Han Chan, E-mail: chantsunghan@gmail.com
% Date: Dec. 26, 2007
% -------------------------------------------------------
% Reference: 
% T.-H. Chan, W.-K. Ma, C.-Y. Chi, and Y. Wang, ``A convex analysis 
% framework for blind separation of non-negative sources," submitted to
% IEEE Trans. Signal Process., July 2007. Accepted subject to minor
% revision, Dec 2007.
%======================================================================
% function flag= is_ext_pt(C,d,alpha,tol)
%  Verify if alpha is an extreme point of the polyhedral set
%     { a | Ca+ d => 0 }
%  It essentially implements Lemma 6.
%---Inputs----
% (C,d) is the set of polyhedron parameters
% alpha is the point to be tested
% tol specifies the numerical tolerance (say, 1e-3)
%---Output----
% flag= 1 if alpha is an extreme point; and flag= 0 otherwise.
%=======================================================================

function flag= is_ext_pt(C,d,alpha,tol)

[L D]=size(C); 
T=C(find(abs(C*alpha+d)<tol),:); % the set of row vectors of C satisfying c_i*alpha + d = 0 
%-------------check if the set T has N-1 linear independent vectors----------------
sin_value=abs(svd(T)); 
temp=sin_value./sum(sin_value);
ran=sum((temp>tol));
flag=(ran==D); % if the set T has N-1 linear independent vectors, flag=1; otherwise, flag=0.

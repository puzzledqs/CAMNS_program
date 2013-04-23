%=====================================================================
% Programmers: 
% Wing-Kin Ma, E-mail: wkma@ieee.org
% Tsung-Han Chan, E-mail: chantsunghan@gmail.com
% Date: Sept. 03, 2009
% -------------------------------------------------------
% Reference: 
% T.-H. Chan, W.-K. Ma, C.-Y. Chi, and Y. Wang, ``A convex analysis 
% framework for blind separation of non-negative sources," 
% IEEE Trans. Signal Process., vol. 56, no. 10, pp. 5120-5134, Oct. 2008.
%
%======================================================================
% A practical implementation of the CAMNS-LP method
%
% function hS = CAMNS_LP(X,N)
%======================================================================
% Output: 
% hS is the L-by-N extracted soruce matrix, where L is the data length.
%---------------------------------------------------
% Inputs:
% X is the L-by-M observation matrix, where M is the number of
% observations.
% N is the number of sources. 
%========================================================================

function hS=CAMNS_LP(X,N)

%----------- Define default parameters------------------
TOL_LP= 1e-3; % tolerance for (small) numerical errors in LP
TOL_EXT= 1e-6; % tolerance for extreme-point validation
TOL_ZEROS= 1e-6; % tolerance for eliminating zero observation points

%-----------(Modification 1) Remove All Zero vectors in Data Set X--------
% It is not meaningful to analyze the zero mixtures. 
[L,M]=size(X); % the given data set
index=find(sum(abs(X'))>=TOL_ZEROS);
Xn=X(index,:); 
LL=length(index);

%-----------Affine Set Fitting [Proposition 1]--------------
d= Xn*ones(M,1)/M;
% For computational efficiency, we use SVD of R instead of EVD of R*R'.
[C,Sigma,V]= svds(Xn-d*ones(1,M),N-1,'L'); 

%--------LP Extreme-Point Finding Algorithm [Table 1]---------------
%------------Step 1
el= 0; % the number of the extreme point estimation
Q1= zeros(LL,1); 
hS= [];
lp_cnt= 0; % the counts of the number of LPs performed 
while el<N %----------Step 7 (begin)
    %--------------Step 2
    w= randn(LL,1); w2= Q1'*w; % Randomly generate a vector w
    r= w- Q1*w2; 
    %--------------Step 3
    b= -C'*r;
    pars.fid= 0; % Make SeDuMi quiet. Set pars.fid= 1 if you want SeDuMi to be verbose
    A= -C'; c= d; K.l= LL;
    tic; [x1, alpha1]= sedumi(A,b,c,K,pars); ttime= toc; lp_cnt= lp_cnt+1; % solve the min LP
    fprintf('%dth LP: running time= %2.5f\n',lp_cnt,ttime);
    %---------------------------------------
    tic; [x2, alpha2]= sedumi(A,-b,c,K,pars); ttime= toc; lp_cnt= lp_cnt+1; % solve the max LP
    fprintf('%dth LP: running time= %2.5f\n',lp_cnt,ttime);
    %-------------Step 4
    if el==0
        % (Modification 2) To play safe, employ the extreme point
        % validation (Lemma 6) to check the obtained optimal solutions.
        % We will reject the solution if it's not an extreme point (not commonplace by our experience).
        if is_ext_pt(C,d,alpha1,TOL_EXT),
            hS= [ hS C*alpha1+d ]; fprintf('Find a new extreme pt (minimizing LP).\n');
        end;
        if is_ext_pt(C,d,alpha2,TOL_EXT),
            hS= [ hS C*alpha2+d ]; fprintf('Find a new extreme pt (maximizing LP).\n');
        end;
    else
        p_star= abs(r'*(C*alpha1+d));
        q_star= abs(r'*(C*alpha2+d));
        % (Modification 3) There may have some small numerical error with 
        % the LP solution. Hence p_star or q_star may not be exactly zero,
        % even though they are supposed to be zero theoretically. A threshold is
        % in place to decide the acceptance/rejection of the obtained solutions.
        % Also, extreme point validation is employed just to play safe.
        if p_star/(norm(r)*norm(C*alpha1+d))>= TOL_LP;  
            if is_ext_pt(C,d,alpha1,TOL_EXT),
                hS= [ hS C*alpha1+d ]; fprintf('Find a new extreme pt (minimizing LP).\n');
            end;
        end
        if q_star/(norm(r)*norm(C*alpha2+d))>= TOL_LP;  
            if is_ext_pt(C,d,alpha1,TOL_EXT),
                hS= [ hS C*alpha2+d ]; fprintf('Find a new extreme pt (maximizing LP).\n');
            end;
        end;   
    end;
    %------------Step 5
    el= size(hS,2); % Update el to be the number of columns of hS
    %------------Step 6
    if el > 0, [Q1,R]= qr(hS,0); end; % Apply thin QR decomposition
end %------------Step 7 (end)
%--------------when el>N (Modification 4)----------------
% When the number of obtained extreme points happens to be greater than 
% the number of true sources (due to numerical errors or violation of
% assumptions), we truncate by selecting the optimal solution at the 
% last run that yields a higher optimal value
if el>N 
    hS=(p_star> q_star)*hS(:,1:1:N)+(p_star< q_star)*hS(:,[1:N-1 N+1]);
end
Y = zeros(L,N);
Y(index,:) = hS;
hS = Y;



function [zhat L ztilde Utilde data] = F_sensor_CRSNC_approxnt(A, k, MAXITER, l)
% function [zhat L ztilde Utilde data] = sens_sel_randomized_half_approxnt(A, k, MAXITER, l)
% Solves the problem
%	maximize log det (sum_{i=1}^m z_i a_i a_i^T) + kappa sum_{i=1}^m(log(z_i)+ log(1-z_i))
%	subject to sum(z) = k
%			   0 <= z_i <= 1, i=1,..., m
% variable z in R^m
% problem parameters kappa (>0), a_1, ..., a_m in R^n
%
% see paper Sensor Selection via Convex Optimization
% www.stanford.edu/~boyd/papers/sensor_selection.html
%
% Nov 2007 Siddharth Joshi & Stephen Boyd

% Newton's method parameters
%MAXITER  = 1000;
NT_TOL = 1e-3;
GAP = 1.005;
% Backtracking line search parameters
alpha = 0.1;
beta = 0.5;

[m n] = size(A);
z = ones(m,1)*(k/m); % initialize
g = zeros(m,1);
ones_m = ones(m,1);
kappa = log(GAP)*n/m; 
% guarantees GM of lengths of semi-axes of ellipsoid corresponding to 
% ztilde <= 1.01 from optimal

% fprintf('\nIter.  Step_size  Newton_decr.  Objective  log_det log_det_dis time \n');
AzA=A'*diag(z)*A;
fz = -log(det(A'*diag(z)*A)) - kappa*sum(log(z) + log(1-z));

% fprintf('   0\t  -- \t     --   %10.3f  %10.3f\n', -fz, log(det(A'*diag(z)*A)));
iflg=0;
lorg=l;
lh=round(l/2);
iflgconverge=max(floor(m/l),1);
tic
idata=0;
for iter=1:MAXITER
%      l=lorg;
%      if(mod(iter,10)==0)
%          l=m;
%      end
        [~,ind]=sort(z);
        S=zeros(l,m);
        Stmp = speye(m);
        S(1:lh  ,:) = Stmp(ind(end-lh+1:end),:);
        Stmp = Stmp(ind(randperm(m-lh)),:);
        S(lh+1:l,:) = Stmp(1:(l-lh),:);
        ones_l=ones(l,1);
        B = S*A;
        zz= S*z;
        W = inv(AzA);
        V = B*W*B';
        g = -diag(V)- kappa * (1./zz - 1./(1-zz));   
        ZR   = (1./(zz.^2) + 1./((1-zz).^2));
        %size(ZR)
        %size(S)
        RR =  diag(ZR);
        %RR =  S .* repmat(ZR',m,1) * S';    
        H = (V.*V) + kappa*RR ;

    R = chol(H);
    Hinvg = (R\(R'\g));
    Hinv1 = (R\(R'\ones_l));
    dS = -Hinvg + ((ones_l'*Hinvg) / (ones_l'*Hinv1))*Hinv1;
    dz= S'* dS;
    deczi = find(dz < 0);
    inczi = find(dz > 0);
    s = min([1; 0.99*[-z(deczi)./dz(deczi) ; (1-z(inczi))./dz(inczi)]]);

    while (1)
        zp = z + s*dz;
        AzAp=AzA+s*B'*diag(dS)*B;
        fzp = -log(det(AzAp)) - kappa*sum(log(zp) + log(1-zp));
        %fz
        %fz+alpha*s*g'*dS
        %fzp
        if (fzp <= fz + alpha*s*g'*dS)
            break;
        end
        s = beta*s;
    end
    z = zp; fz = fzp; AzA=AzAp;
    zsort=sort(z); thres=zsort(m-k); zhat=(z>thres);
    fzhat= -log(det(A(zhat,:)'*A(zhat,:)));
    if mod(iter,10)==0 
    idata=idata+1;
    time10=toc;
    % fprintf('%4d %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f \n', iter, s, -g'*dS/2, -fz, log(det(AzA)),-fzhat,time10);
    data(idata,:)=[iter -fz log(det(AzA)) -fzhat time10];
    end 
    if(-g'*S*dz/2 <= NT_TOL)
        iflg=iflg+1;
        %break;
    else
        iflg=0;
    end
    if(iflg == iflgconverge)
        break
    end
end
    
zsort=sort(z); thres=zsort(m-k); zhat=(z>thres);
L = log(det(A'*diag(zhat)*A));
ztilde = z; 
Utilde = log(det(A'*diag(z)*A)) + 2*m*kappa;

time10=toc;
% fprintf('%4d %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f \n', iter, s, -g'*dS/2, -fz, log(det(AzA)),-fzhat,time10);
idata=idata+1;
data(idata,:)=[iter -fz log(det(AzA)) -fzhat time10];

% toc
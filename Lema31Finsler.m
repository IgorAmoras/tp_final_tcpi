function out = Lema31Finsler(A,Bu,Bw,C,Du,Dw, varargin)
lambda = 1;
gamma = sdpvar;

if nargin > 6
    options = struct(varargin{:});
    if(isfield(options, 'gamma'))
        gamma = options.gamma;
    end
    if(isfield(options, 'delta'))
        lambda = options.delta;
    end
end

nx = size(A,1);
nu = size(Bu,2);
nul = size(Bu, 1);
nw = size(Bw,2);
nwl = size(Bw,1);
ndc = size(Dw, 2);
ndl = size(Dw,1);

W = sdpvar(nx,nx,'symmetric');
Z = sdpvar(nu,nx,'full');
M = sdpvar(nx,nx, 'full');
X = sdpvar(nx,nx, 'full');

LMIs = (W >= 0);

Y = sdpvar(nu, nul, 'full');
T = blkvar;
T(1,1) = A*M + M'*A' + Bu*Y + Y'*Bu';
T(2,1) = X + lambda*A*M + lambda*Bu*Y-M';
T(3,1) = Bw';
T(4,1) = C*M + Du*Y;

T(2,2) = -lambda*(M + M');
T(3,2) = lambda*Bw';

T(3,3) = -gamma*eye(ndc);
T(4,3) = Dw;

T(4,4) = -gamma*eye(ndl);

T = sdpvar(T);

LMIs = LMIs + (T <= 0);


obj = gamma;

options=sdpsettings('verbose',1,'warning',1,'solver','sedumi');

out.sol=solvesdp(LMIs,obj,options);
warning('off','YALMIP:strict');
res=min(checkset(LMIs));

if out.sol.problem == 0
    out.feas = 1;
else
    out.feas = 0;
end

out.W = double(W);
out.K = (double(Y) * double(M)^-1);
out.gamma = sqrt(double(gamma));
out.res = res;

end
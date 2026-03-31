function verifyMFPTAnisotropic

close all;

global o;

newcolors = [220, 200, 0; 161, 218, 180; 65, 182, 196; 44, 127, 184;
37,  52, 148; 8,  29,  88; 0,   0,   0] / 255;
colororder(newcolors)

cols = get(gca,'colororder');


o.rho = 0.0; o.R0 = 6;

o.x0 = o.R0; o.y0 = 0; o.z0 = 0; o.np = 1e4;
o.d = 3; % spatial dimension
o.z = 1; % radial (z=1) or tangential (z= -1).
o.pm = 1; % outward (+1) or inward motion (-1)

o.int_tol = 1e-12;

x = linspace(o.rho,o.R0,250);

%alpha_range = [ 0.99 0.5 0.25 -0.25 -0.5,-0.75];
alpha_range = [0.99 0.5];

num_cases = 6;
do_SSA = 0;

Df = (0.06)^2;
ndash = 6;
lineStyles = repmat({'-'}, 1, ndash);
lineStyles{end} = ':';

str = {
    '$\hat \gamma (r) = -\hat r, $ $\kappa_1(r)= 0.02 r$',...
    '$\hat \gamma (r) = -\hat r, $ $\kappa_2(r)= 0.1 e^{-0.5 r }$',...
    '$\hat \gamma (r) = +\hat r, $ $\kappa_2(r)=0.1 e^{-0.5 r }$',...
    '$\hat \gamma (r) = +\hat r, $ $\kappa_1(r)=0.02 r$',...
    '$\hat \gamma (r) = +\hat r, $ $\kappa\to\infty$ $\mbox{ballistic}$',...
    '$D = \sigma^2/3 \mu$ $\mbox{diffusive}$'};




%for i = [4 5]
for i = [1 2 3 4 5 6]

    if i==6
        do_SSA = 1;
        o.sigma = 0.06; o.mu = 1;
        o.pm = 1; 
        o.kap = @(r) 0*r;   
        if o.d == 2, o.k = @(x,y) 0*x; end
        if o.d == 3, o.k = @(x,y,z) 1e-3 + 0*x; end
    elseif i==2
        o.pm = -1;
        do_SSA =1;
        o.sigma =0.06; o.mu = 1;        
        A = 0.1; lambda = 2;
        o.kap = @(r) A*exp(-r/lambda);   
        if o.d == 2, o.k = @(x,y) A*exp(-sqrt(x.^2+y.^2)/lambda); end
        if o.d == 3, o.k = @(x,y,z) A*exp(-sqrt(x.^2+y.^2+z.^2)/lambda); end
    elseif i==3
        o.pm = 1;
        do_SSA =1;
        o.sigma = 0.06; o.mu = 1;
         A = 0.1; lambda = 2;
         o.kap = @(r) A*exp(-r/lambda);   
         if o.d == 2, o.k = @(x,y) A*exp(-sqrt(x.^2+y.^2)/lambda); end
         if o.d == 3, o.k = @(x,y,z) A*exp(-sqrt(x.^2+y.^2+z.^2)/lambda); end
    elseif i==1
        o.pm =-1; 
        do_SSA =1;
        o.sigma = 0.06; o.mu = 1; lambda = 50;
        o.kap = @(r) r/lambda;     
        if o.d == 2, o.k = @(x,y) sqrt(x.^2+y.^2)/lambda; end
        if o.d == 3, o.k = @(x,y,z) sqrt(x.^2+y.^2+z.^2)/lambda; end
    elseif i==4
        o.pm = +1;
        do_SSA =1;
        o.sigma = 0.06; o.mu = 1; lambda = 50;
        %o.sigma = 0.006; o.mu = 0.01; lambda = 50;
        o.kap = @(r) r/lambda;     
        if o.d == 2, o.k = @(x,y) sqrt(x.^2+y.^2)/lambda; end
        if o.d == 3, o.k = @(x,y,z) sqrt(x.^2 + y.^2 + z.^2)/lambda; end
    else
        o.pm= 1;
        do_SSA =1;
        o.sigma = 0.06; o.mu = 1;
        o.kap = @(r) 2e2 + 0*r;     
        if o.d == 2, o.k = @(x,y) 1e2 + 0*x; end
        if o.d == 3, o.k = @(x,y,z) 1e2 + 0*x; end
    end

    o.D = o.sigma^2/(2*o.mu);
    if (o.d == 2)
        o.alpha = @(x) besseli(2,o.kap(x))./besseli(0,o.kap(x));
        o.beta = @(x) besseli(1,o.kap(x))./besseli(0,o.kap(x));
    end
    if (o.d == 3)  
        o.D = (2/3)*o.D;
        o.alpha = @(r) alpha3D(o.kap(r));
        o.beta = @(r) beta3D(o.kap(r));
    end
   o.F = @(r) (r/o.rho).^( (o.d-1)*(1-o.z*o.alpha(o.rho))/(1+o.z*o.alpha(o.rho)) ) .* exp( integral(@(s) ((o.d-1)./s).*( (1-o.z*o.alpha(s))./(1+o.z*o.alpha(s)) - (1-o.z*o.alpha(o.rho))./(1+o.z*o.alpha(o.rho)) ) + ... 
       + o.pm*((o.z+1)/2)*o.d*o.mu*o.beta(s)/(o.sigma*(1+o.z*o.alpha(s)) ), o.rho,r,'RelTol',0,'AbsTol',o.int_tol) );

   % o.F = @(r) exp( integral(@(s) ((o.d-1)./s).*( (1-o.z*o.alpha(s))./(1+o.z*o.alpha(s))) + o.pm*((o.z+1)/2)*o.d*o.mu*o.beta(s)/(o.sigma*(1+o.z*o.alpha(s)) ), o.rho,r) );

    % Reflecting inside, absorbing outside.
    o.BC = [0 1 1 0];

    small = o.sigma/(o.mu*o.R0);

    %h = besseli(2,k(0))/besseli(0,k(0));
    %MFPT_Ex = @(r) R0^2 * (1- (r/R0).^2)./(D*(2));

    guess = @(x) [100;100];
    solinit = bvpinit(x,guess);
    opts = bvpset('RelTol',1e-8,'AbsTol',1e-8);
    sol = bvp4c(@bvpfcn, @bcfcn, solinit,opts);
    figure(2);
    hold on
    sol.y(1,end) = 1e-4;
    qq{i} = plot(sol.x, sol.y(1,:),'color',cols(i,:), ...
        'LineStyle',  lineStyles{i}, 'linewidth',3,'Displayname', str{i});
    xlim([o.rho o.R0])

    % r_range = linspace(o.rho,o.R0,12);
    % T_closed = zeros(size(r_range));
    % for k = 1:length(r_range)
    %     T_closed(k) = integral(@(s) intgH(s), r_range(k),o.R0)/o.D;
    % end
    % 
    % plot(r_range,T_closed,'kx','linewidth',3);
    % 
    % 
    % drawnow
    %return

    %ax = gca;
    %ax.YAxis.TickLabelFormat = '%,.1f';

    MFPT_BVP = sol.y(1,1);

    fprintf(1,'Case = %g, Small = %g\n',i,small);

    %plot(rc, Ta(rc),'ko','displayname','exact');

    title_str = '$\sigma = 0.06, \mu = 1$, $d=3$';


    set(gcf,'color','w'); set(gca,'FontName', 'Times New Roman', 'fontsize',28);
    title(title_str, 'Interpreter', 'latex');
    xlabel('$r$','interpreter','latex', 'fontsize', 28)
    ylabel('$T(r)$','interpreter','latex','rotation', 90, 'fontsize', 28)
    ytickformat('%1.1f')
    yticks = [10.^-3, ... % 1, 2, 5
           10.^-1, ... % 10, 20, 50
           10.^1, ... % 100, 200, 500
           10.^3, ... % 1000, 2000, 5000
           10.^5];
    set(gca, 'YTick', yticks);
    yscale('log'); ylim([1e-3 1e5])

    %return
    %yyaxis right
    %plot(sol.x,alpha(sol.x),'-r')
    %set(gca,'Interpreter',Latex)
    %exportgraphics(figure(2),'AnnulusBothAbs.png','resolution',300)
    %exportgraphics(figure(2),'AnnulusInnerAbs.png','resolution',300)
    %plot(x,MFPT_Ex(x),'kx', DisplayName = "Constant");

    %load MFPTData.mat
    %load MFPTDataInnerExit2.mat
    %scatter(MFPTData.r,MFPTData.T,'k', 'square', 'filled','displayname','Particle');


    %set(gca,'YScale','log')
    if (do_SSA)
 
    if o.z == 1 % radial
        if (o.d==2)
            o.nu1 = @(x,y) o.pm*x./sqrt( x.^2 + y.^2 );
            o.nu2 = @(x,y) o.pm*y./sqrt( x.^2 + y.^2 );
        end

        if (o.d ==3)
            o.nu1 = @(x,y,z) o.pm*x./sqrt( x.^2 + y.^2 + z.^2 );
            o.nu2 = @(x,y,z) o.pm*y./sqrt( x.^2 + y.^2 + z.^2 );
            o.nu3 = @(x,y,z) o.pm*z./sqrt( x.^2 + y.^2 + z.^2 );
        end
    else % tangential
        o.nu1 = @(x,y)  o.pm*y./sqrt( x.^2 + y.^2 );
        o.nu2 = @(x,y) -o.pm*x./sqrt( x.^2 + y.^2 );
    end

    % x = -1 + zeros(10000,1);
    % y = 0 + zeros(10000,1);
    % z = zeros(10000,1);
    % 
    % % single-modal Fisher - concentration alignment
    % 
    % q = @(x,y,z,n1,n2,n3) (o.k(x,y,z).*exp(o.k(x,y,z).*(o.nu1(x,y,z).*n1 + o.nu2(x,y,z).*n2 + o.nu3(x,y,z).*n3))./(4*pi*sinh(o.k(x,y,z)))).*(sqrt(1-z.^2));
    % qM = @(x,y,z) (o.k(x,y,z).*exp(o.k(x,y,z))./(4*pi*sinh(o.k(x,y,z)))).*(sqrt(1-z.^2));
    % [th,ph] = sampleFisher(x,y,z,q,qM);
    % 
    % mu = [o.nu1(x,y,z) o.nu2(x,y,z) o.nu3(x,y,z) ];
    % X = rvmf3_pairwise(o.k(x,y,z), mu);
    % 
    % figure;
    % 
    % th_bins = linspace(0,2*pi,40);
    % z_bins = linspace(-1,1,20);
    % 
    % histogram2(th,cos(ph),th_bins,z_bins);
    % 
    % return


    fac = 0.1:0.1:0.9;
    run_MFPT = zeros(length(fac),1);
    for j = 1:length(fac)
            o.x0 = fac(j)*o.R0;
            if (o.d == 2), run_MFPT(j) = run_SSA(o); end
            if (o.d == 3), run_MFPT(j) = run_SSA_3D(o); end
    end
    hold on
    hh{i} = plot(o.R0*fac,run_MFPT,'s','markersize',12,'LineWidth',2, ...
        'markeredgecolor',cols(i,:),'markerfacecolor','none', ...
        'DisplayName','Particle Simulations');
    drawnow
    hold off;
    end

end
if (do_SSA)
    lgd = legend([qq{1} qq{2} qq{3} qq{4} qq{5} qq{6} hh{1}]);
    lgd.Box = 'off';
    lgd.FontSize=20;
else
    lgd = legend([qq{1} qq{2} qq{3} qq{4} qq{5} qq{6}]);
    lgd.Box = 'off';
    lgd.FontSize=20;
end
lgd.Interpreter = 'latex';
lgd.Location = 'southwest';

if (o.d==2)
    exportgraphics(figure(2),'MFPT_radial_d=2.png','resolution',300);
end

if(o.d==3)
    exportgraphics(figure(2),'MFPT_radial_d=3.png','resolution',300);
end

return

function out = intgH(s)
global o
out = zeros(size(s));

for j = 1:length(s)
    %out(j) =  intH2(s(j))/o.F(s(j));
    out(j) =  integral(@(z)  intH1(z) ,0,s(j))/o.F(s(j));
    %
end


function out = intH1(s)
global o
out = zeros(size(s));

for j = 1:length(s)
    out(j) = o.F(s(j))/(1+ o.z*o.alpha(s(j)));
end

function out = intH2(s)
global o
out = zeros(size(s));
for j = 1:length(s)
    out(j) = integral(@(z)  intH1(z) ,0,s(j));
end


function dydx = bvpfcn(x,y)

global o
dydx = zeros(2,1);
if x>0
    dydx = [y(2)
        -1/(o.D*(1+o.z*o.alpha(x))) - (y(2)/x)*(o.d-1)*(1-o.z*o.alpha(x))/(1+o.z*o.alpha(x)) - o.pm*y(2)*((o.z+1)/2)*o.d*o.mu*o.beta(x)/(o.sigma*(1+o.z*o.alpha(x))) ];
else
    dydx = [y(2)
        -1/( o.D*( 1+o.z*o.alpha(x) + (o.d-1)*(1-o.z*o.alpha(x)) ) ) ];
end

function res = bcfcn(ya,yb)
global o
res = [ (o.BC(1)*ya(1) + o.BC(2)*ya(2))
        (o.BC(3)*yb(1) + o.BC(4)*yb(2))];

function out = run_SSA(o)

% single-modal VonMises - concentration alignment
q = @(x,y,ang) exp(o.k(x,y).*(o.nu1(x,y).*cos(ang) + o.nu2(x,y).*sin(ang)) )./(2*pi*besseli(0,o.k(x,y)));
qM = @(x,y) exp(o.k(x,y))./(2*pi*besseli(0,o.k(x,y)));

% Trap;

func = @(x,y) (x).^2 + (y).^2 - o.R0^2;

nFree = o.np;
indFree = 1:nFree;
t = zeros(nFree,1);
x = o.x0*ones(nFree,1);
y = o.y0*ones(nFree,1);
th = sampleVonMises(x,y,q,qM);

while ( nFree>0 )
    tjump = (-1/o.mu) * log(rand(nFree,1));
    t(indFree) = t(indFree) + tjump;

    x(indFree) = x(indFree) + o.sigma*tjump.*cos(th(indFree));
    y(indFree) = y(indFree) + o.sigma*tjump.*sin(th(indFree));
    % Check for boundary crossing.

    indEsc = find(func(x(indFree),y(indFree))> 0);

    if (~isempty(indEsc))
        indC = indFree(indEsc);
        xc = x(indC); yc = y(indC); thc = th(indC);
        sc = zeros(size(xc)); tj = tjump(indEsc);

        for i = 1:length(indC)
            sc(i) = bisect(@(z) func(xc(i)+z*o.sigma*cos(thc(i)),yc(i)+z*o.sigma*sin(thc(i))),-tj(i),0);
        end
        t(indC) = t(indC)+sc;
        x(indC) = x(indC) + o.sigma*sc.*cos(thc); y(indC) = y(indC) + o.sigma*sc.*sin(thc);
    end
    indFree(indEsc) = [];
    nFree = numel(indFree);

    th(indFree) = sampleVonMises(x(indFree),y(indFree),q,qM);

end

out = mean(t);

function out = run_SSA_3D(o)

% single-modal Fisher - concentration alignment
q = @(x,y,z,n1,n2,n3) (o.k(x,y,z).*exp(o.k(x,y,z).*(o.nu1(x,y,z).*n1 + o.nu2(x,y,z).*n2 + o.nu3(x,y,z).*n3))./(4*pi*sinh(o.k(x,y,z))));
qM = @(x,y,z) (o.k(x,y,z).*exp(o.k(x,y,z))./(4*pi*sinh(o.k(x,y,z))));

% Trap;

func = @(x,y,z) x.^2 + y.^2 + z.^2 - o.R0^2;

nFree = o.np;
indFree = 1:nFree;
t = zeros(nFree,1);
x = o.x0*ones(nFree,1);
y = o.y0*ones(nFree,1);
z = o.z0*ones(nFree,1);

X = sampleFisher(x,y,z,q,qM,o);

%X = rvmf3_pairwise(o.k(x,y,z), [o.nu1(x,y,z) o.nu2(x,y,z) o.nu3(x,y,z) ]);

while ( nFree>0 )
    tjump = (-1/o.mu) * log(rand(nFree,1));
    t(indFree) = t(indFree) + tjump;

    x(indFree) = x(indFree) + o.sigma*tjump.*X(indFree,1);
    y(indFree) = y(indFree) + o.sigma*tjump.*X(indFree,2);
    z(indFree) = z(indFree) + o.sigma*tjump.*X(indFree,3);
    % Check for boundary crossing.

    indEsc = find(func(x(indFree),y(indFree),z(indFree))> 0);

    if (~isempty(indEsc))
        indC = indFree(indEsc);
        xc = x(indC); yc = y(indC); zc = z(indC); 
        X_c = X(indC,:);
        
        sc = zeros(size(xc)); tj = tjump(indEsc);

        for i = 1:length(indC)
            sc(i) = bisect(@(z) func(xc(i) + z*o.sigma*X_c(i,1),...
                                        yc(i) + z*o.sigma*X_c(i,2),...
                                        zc(i) + z*o.sigma*X_c(i,3)),-tj(i),0);
        end

        t(indC) = t(indC) + sc;
        
        x(indC) = x(indC) + o.sigma*sc.*X_c(:,1); 
        y(indC) = y(indC) + o.sigma*sc.*X_c(:,2);
        z(indC) = z(indC) + o.sigma*sc.*X_c(:,3);

    end
    indFree(indEsc) = [];
    nFree = numel(indFree);

    X(indFree,:) = sampleFisher(x(indFree),y(indFree),z(indFree),q,qM,o);

end

out = mean(t);

function out = sampleVonMises(x,y,q,qM)
n = length(x); indS = 1:n;
M = qM(x,y);

out = zeros(size(x));

while (~isempty(indS))
    r = 2*pi*rand(length(indS),1); u = rand(length(indS),1);
    indAccept = find(u < q(x(indS),y(indS),r)./M(indS) );
    out(indS(indAccept)) = r(indAccept);
    indS(indAccept) = [];
end

function X = sampleFisher(x,y,z,q,qM,o)
n = length(x); indS = 1:n;
M = qM(x,y,z);

X = zeros(n,3);

while (~isempty(indS))

    Z = randn(numel(indS),3);
    Z = Z ./ vecnorm(Z,2,2);

    u = rand(length(indS),1);
    indAccept = find(u < q(x(indS),y(indS),z(indS),Z(:,1),Z(:,2),Z(:,3))./M(indS));
    X(indS(indAccept),:) = Z(indAccept,:);
    indS(indAccept) = [];

end


function out = alpha3D(s)
    
    out = zeros(size(s));
    z = s(s~=0);
    out(s~=0) = (z.^2 - 3*z.*coth(z) + 3)./(z.^2);

function out = beta3D(s)
    
    out = zeros(size(s));
    z = s(s~=0);
    out(s~=0) = (z.*coth(z) - 1)./(z);

    function X = rvmf3_pairwise(kappa, mu)
% Pairwise vMF sampler on S^2
% kappa : Nx1 (or 1xN) vector of concentrations (>=0) or scalar
% mu    : (N x 3) or (3 x 1) mean directions (will be normalized)
% Output:
%   X   : N x 3, one sample for each (kappa(i), mu(i,:))

kappa = kappa(:);
N = numel(kappa);

if nargin < 2 || isempty(mu)
    mu = [0;0;1];
end

if isvector(mu) && numel(mu)==3
    mu = mu(:)'/norm(mu);
    mu = repmat(mu, N, 1);
else
    assert(all(size(mu)==[N,3]), 'mu must be (3x1) or (N x 3)');
    mu = mu ./ vecnorm(mu,2,2);
end

% Helper vector 'a' row-wise
a = zeros(N,3);
useZ = abs(mu(:,3)) < 0.999;
a(useZ,3) = 1;                 % [0 0 1]
a(~useZ,1) = 1;                % [1 0 0]

b1 = a - sum(a.*mu,2).*mu;
b1 = b1 ./ vecnorm(b1,2,2);
b2 = cross(mu, b1, 2);

u   = rand(N,1);
phi = 2*pi*rand(N,1);

X = zeros(N,3);

mask0 = (kappa==0);
if any(mask0)
    Z = randn(sum(mask0),3);
    Z = Z ./ vecnorm(Z,2,2);
    X(mask0,:) = Z;
end

maskp = ~mask0;
if any(maskp)
    kap = kappa(maskp);
    up  = u(maskp);
    ph  = phi(maskp);

    eps = exp(-2*kap);
    w   = 1 + (1./kap).*log(up + (1-up).*eps);

    s  = sqrt(max(0, 1 - w.^2));
    cp = cos(ph); sp = sin(ph);

    b1p = b1(maskp,:); b2p = b2(maskp,:); mup = mu(maskp,:);

    X(maskp,:) = (s.*cp).*b1p + (s.*sp).*b2p + w.*mup;
end

X = X ./ vecnorm(X,2,2);

function fzero = bisect(f,a,b)
left = a;   % Left end point.
right = b;  % Right end Point.
mid = (a+b)/2;   % Mid point.

j = 0;

f_left = f(left);
f_right = f(right);
f_mid = f(mid);

j = j+3;

% Check there is a sign change on the given interval, quit otherwise.
if ( f_left*f_right >= 0 )
   error('No sign change on (a,b)'); 
end
% Continue until the left and right endpoints are 1e-5 apart.
while ( abs(right-left) > 1e-5 )    
    if ( f_left*f_mid < 0 )        
        % The root is on (left, mid)        
          right = mid;
       % f_right = f_mid;        
     else 
      % The root is on [mid, right)        
        left = mid;
       % f_left = f_mid;          
    end    
    % Update the values of mid and f_mid.    
    mid = 0.5*(left + right);
    f_mid = f(mid);  
    j = j+1;
    
end
% Output the approximation of the root.

fzero = 0.5*(left+right);

clear all
close all


%% Parameters
tic

% catch/slip adhesion
ctot = 300; %/um^2          cadherin density
fs   = .0168; % nN          slip force
fc   = .003471; % nN        catch force
phis = -1.778; %            slip parameter
phic = 2.942; %             catch parameter
kon  = 1.2; %               cadherin binding rate
kb   = 0.5;  % nN/um        cadherin bond stiffness
Fmax = 22e-3; % nN          force for rupture

% contractility
alpha0= 3; % /kPa           feedback parameter in absence of bonds
beta  = 27.7; % /kPa        chemical stiffness
alphac= 17; % /kPa          feedback parameter with bonds          
kcb   = 75; %               cadherin dependence
Kc    = 1; % kPa            cell stiffness
kr    = .4e-4; % kPa/s      rate of SF remodeling
rho0  = 1; % kPa            quiescent myosin density

% polymerization
kp   = 2.5e-3; % kPa/s      rate
rhox = 1.75;  % kPa         contractility dependence
pmax = -4.725; % kPa        maximum polymerization stress

% lattice parameters
m   = 21; %                 number junctions
ang = 30; %                 set cell shape
Lx  = 10; % um              cell radius
Lg0 = 10; % um              initial gap length
Lj0 = 1; % um               junction length after formation
Km = .1; % kPa              membrane element stiffness
Ko = .01;   %               junction stiffness when open kPa
xj  = zeros(m,1); %         initially no junctions
A0 = 11.5/m; % um^2         csa

%% simulation setup
iter = 5;         
tmin = 0;          % start time
tmax = 1.6e4;      % end overall solution at time tmax
xtol = 1e-5;       % tolerance
tout = [];
xout = [];
xjvec = [];
Xo0   = zeros(3*m,1);
xtol2v= zeros(3*m,1);

% initial conditions for m junctions
for j = 1:m
    Xo0(3*j-2) = 0;     % initial cb
    Xo0(3*j-1) = rho0;  % initial rho
    Xo0(3*j)   = -rho0;    % initial polymerization
    
    xtol2v(3*j-2) = 5e-4;  % error tolerances
    xtol2v(3*j-1) = 1e-3;
    xtol2v(3*j)   = 1e-3;
end

% create lattice
[n,elem,BC0s, BC1s] = make_lattice(m, ang, Lx, Lg0);
n0 = n;

% plot undeformed configuration
h = figure(1);
plot_lattice(n0,elem,xj,m,1, 0,ctot,0)

% Calculate reference element configs
elem_L0 = elem_ref(n, elem);
elem_L00 = elem_L0;

% Non-uniform lattice element stiffness
[s1, ~] = size(elem);
kvec(1:m) = Kc.*A0./elem_L0(1);  % fiber
kvec((m+1):(2*m-1)) = Km.*ones(m-1,1).*A0./elem_L0(m+1);  % cortex
kvec((2*m):(3*m-1)) = Ko.*A0./elem_L0(end);    % gap/junction

% lattice properties
props = zeros(m,2);
props(1:m,1) = alphac;
props(1:m,2) = pmax;

% Apply negligible deformation to aid convergance
d = 1e-6*Lx;
for i = 1:m
    x0 = n(end-m+i,1);
    y0 = n(end-m+i,2);
    new_val = [x0,y0+d];
    n(end-m+i,:) = new_val;
end

%% ode solver 
inc = 0;
inc_max = iter*m;
while inc < inc_max

tspan=[tmin tmax];  % max simulated time for call
options1 = odeset('RelTol',xtol,'AbsTol',xtol2v, 'Events', @(tt1,xx1) events1(tt1,xx1,kvec,n, elem, elem_L0,...
    BC0s,BC1s,A0,ctot,Lj0,kb,Ko,Fmax, xj, m,tmin)); % may need to include abstol
% call main function 
[tt1, xx1] = ode15s(@(tt1,xx1) xjunODE(tt1,xx1,kb,Kc,alpha0,beta,...
    rho0,Fmax,Lj0,kr,kp,rhox,Lg0,kcb,Ko,kvec,n, elem, elem_L0, BC0s,BC1s, A0,xj,m,props,...
    fs, fc, phis, phic, ctot, kon),...
    tspan,Xo0,options1);

% Accumulate output.  This could be passed out as output arguments.
tout = [tout; tt1];
xout = [xout; xx1];
xjvec = [xjvec;ones(length(tt1),length(xj)).*xj'];

% set conditions for next ODE call
tmin   = tout(end);
X = ['ode break, inc=',num2str(inc),'/',num2str(inc_max),', time=',num2str(tmin)];
disp(X);

% overall end at some time, tmax
if (tmin>=tmax)
    break
end
for j = 1:m
    Xo0(3*j-2) = xx1(end,3*j-2);   %  Pb 
    Xo0(3*j-1) = xx1(end,3*j-1);   %  rho
    Xo0(3*j)   = xx1(end,3*j);     %  polymerization
end

% change condition open/closed
xjw = 10.*ones(length(xj),1);
while ~isequal(xjw,xj)
    xjw = xj;
    [elem_L0, n, xj] = events_update(tout(end),xout(end,:),kvec,n, elem, elem_L0, BC0s,BC1s,A0,...
        ctot,Lj0,kb,Ko,Fmax, xj,elem_L00,n0, m);
end


inc = inc+1;

end
toc


%% Post simulation 
tx = 0;
h = figure(2);
for i = 1:length(tout)
    
    for j = 1:m
        if xjvec(i,j)==0
            elem_L0(end-m+j) = elem_L00(end-m+j); % no junction
            kvec((end-m+j)) = Ko*A0/elem_L0(end-m+j);
            n(end-2*m+j,:) = n0(end-2*m+j,:);

        elseif xjvec(i,j) == 1
            elem_L0(end-m+j) = (Lj0); 
            xcb = xout(i,3*j-2)*ctot;   
            xKj   = Lj0*kb*xcb; 
            kvec(end-m+j) = xKj*A0/elem_L0(end-m+j);
            v0 = n(end-m+j,:); v1 = n(end-2*m+j,:); v = v1-v0;
            u = v./norm(v);
            n(end-2*m+j,:) = v0 + Lj0*u;
        end       
        xf(j) = (xout(i,3*j-1)+xout(i,3*j)).*A0;
    end
    
    
    dt = tout(i)-tx;
    val = 533;
    if length(tout)<val
        val = length(tout);
    end
    if (dt >=50)

       [n2, eps] = quasi_stat_soln(kvec,n, elem, elem_L0, BC0s,BC1s,A0, xf, tout(i));
        clf(h)

        plot_lattice(n2,elem,xjvec(i,:),m,2,xout(i,:),ctot,1); hold all 
        pbaspect([1 .866 1]);

    ln0 = n(end-m,:);
    for j=1:m
        epsc(j) = eps(j);
        sig(j)  = Kc.*epsc(j) + xout(i,3*j-1) + xout(i,3*j); 
        xdj(j) = elem_L0(end-m+j)*(eps(end-m+j));
        xFj(j) = kb.*xdj(j);
        ln(j,:) = n(j+1,:);
        ld(j) = sqrt((ln(j,1)-ln0(1))^2 + (ln(j,2)-ln0(2))^2);
    end

        tx = tout(i);
        pause(0.1);
        
    end
    tout(i)/tout(end);
end   


for i=m:-3:17

hfig1=figure(3); set(hfig1,'color','w'); hold all
plot(tout./3600-1.75,xout(:,3*i-2).*ctot, '-','linewidth',4); 
box on; haxsY=gca;
set(haxsY,'tickdir','out','linewidth',3,'fontsize',22,...
    'FontName','Arial','FontWeight','bold')
xlabel('time (hrs)');
ylabel('bound cadherin density');
pbaspect([1.5 1 1])
xlim([0 1.5])
ylim([-1 35])

hfig1=figure(4); set(hfig1,'color','w'); hold all
plot(tout./3600-1.75,xout(:,3*i-1), '-','linewidth',4); 
box on; haxsY=gca;
set(haxsY,'tickdir','out','linewidth',3,'fontsize',22,...
    'FontName','Arial','FontWeight','bold')
xlabel('time (hrs)');
ylabel('contractility');
pbaspect([1.5 1 1])
xlim([0 1.5])
ylim([1.1 1.7])

hfig1=figure(5); set(hfig1,'color','w'); hold all
plot(tout./3600-1.75,-xout(:,3*i), '-','linewidth',4); 
box on; haxsY=gca;
set(haxsY,'tickdir','out','linewidth',3,'fontsize',22,...
    'FontName','Arial','FontWeight','bold')
xlabel('time (hrs)');
ylabel('polymerization');
pbaspect([1.5 1 1])
xlim([0 1.5])
ylim([2.5 2.7])
set(gca, 'YTick', [2.5 2.6 2.7])

end
% hfig1=figure(6); set(hfig1,'color','w'); hold all; % rupture force
% plot([0 10],[22 22], 'k:','linewidth',5);
% plot(ld,xFj.*1e3, 'o-','linewidth',5,'color', [0.6350, 0.0780, 0.1840]);
% box on; haxsY=gca;
% set(haxsY,'tickdir','out','linewidth',3,'fontsize',25,...
%     'FontName','Arial','FontWeight','bold')
% pbaspect([1.5 1 1])
% xlim([0 10])
% ylim([-20 30])

%% functions

% main modelling function
function dxdt = xjunODE(time,x,kb,Kc,alpha0,beta,rho0,Fmax,Lj0,kr,kp,rhox,Lg0,kcb,Ko,...
    kvec,n, elem, elem_L0, BC0s,BC1s, A0, xj,m, props,fs, fc, phis, phic, ctot, xkon)
%     time
    dxdt   = zeros(3*m,1); % initialize arrays
    Pb     = zeros(m,1);
    rho    = zeros(m,1);
    sigp   = zeros(m,1);
    alphac = zeros(m,1);
    pmax   = zeros(m,1);
    xcb    = zeros(m,1);
    f      = zeros(m,1);
    xf     = zeros(m,1);
    xKj    = zeros(m,1);
    epsc   = zeros(m,1);
    epsj   = zeros(m,1);
    dj     = zeros(m,1);
    koff   = zeros(m,1);
    kon    = zeros(m,1);
    xalpha = zeros(m,1);
    sigc   = zeros(m,1);
    
    % start  
    for i=1:m
        Pb(i)   = x(3*i-2);  % no. bonds
        rho(i)  = x(3*i-1);  % contractility (myosin density)
        sigp(i) = x(3*i);    % polymerization (extension)

        alphac(i)  = props(i,1);
        pmax(i)= props(i,2);
        xcb(i) = Pb(i)*ctot; % actual number of junction bonds  
        kon(i) = xkon;

        % check if junction is open or closed
        if xj(i) == 0
            kvec(end-m+i) = Ko*A0/elem_L0(end-m+i);
        elseif xj(i) == 1
            xKj(i)   = Lj0*kb*xcb(i); % Effective junction stiffness = bonds x individual bond stiffness
            kvec(end-m+i) = xKj(i)*A0/elem_L0(end-m+i);
        end

       % calculate new nodal coords
         xf(i) = (rho(i)+sigp(i))*A0;
    end

    [~, eps] = quasi_stat_soln(kvec,n, elem, elem_L0, BC0s,BC1s,A0, xf,time);

    for i = 1:m
        epsc(i) = eps(i);
        epsj(i) = eps(end-m+i);
        dj(i)   = epsj(i)*elem_L0(end-m+i);
        if dj(i) < 0
            dj(i) = 0.0;
        end

      % adhesion rate coefficient
        f(i) = kb*dj(i);

        if (xj(i)==0)
          kon(i) = 0;
        end
    
        koff(i) = exp(f(i)/fs - phis) + exp(phic - f(i)/fc);

        xalpha(i) = alpha0 + alphac(i)*(xcb(i)/(kcb+xcb(i)));  
        sigc(i)  = Kc*epsc(i) + rho(i);  

        dxdt(3*i-2) =  (1-Pb(i))*kon(i) - Pb(i)*koff(i);                      % d(cb)/dt
        dxdt(3*i-1) = -kr*(epsc(i) + beta*(rho(i)-rho0) - xalpha(i)*sigc(i)); % d(rho)/dt
        dxdt(3*i)   =  -kp*( 1 - sigp(i)/pmax(i) - rho(i)/(rhox+rho(i)) );    % d(sigp)/dt
    end

end

function [n,elem,BC0s,BC1s] = make_lattice(m, ang, Lx, Lg0)
    % define lattice

    n = [];
    elem = [];
    BC0s = [];
    BC1s = [];
          
    n(end+1,:) = [0 0]; 
    BC0s(end+1,:) = [1,0];
    BC0s(end+1,:) = [1,90];   
    n(end+1,:) = [-Lx 0]; 

    Ly = Lx*tan(deg2rad(ang));
    mx = (m-1)/2;
    dphi = ang/(m-1);
    phi = 0;
    for i = 1:m-1
        phi = phi + dphi;
        x = Lx*cos(deg2rad(phi));
        y = Lx*sin(deg2rad(phi));
        n(i+2,:) = [-x y]; 
    end
    
%     vr(1) = cos(deg2rad(ang)); vr(2) = sin(deg2rad(ang));
%     for i=1:mx
%         y = i*Ly/mx;
%         n(mx+i+2,:) = n((2+(m-1)/2),:) + y.*vr; 
%     end
    
    for i = 1:m
        elem(i,:) = [1,i+1];  %fibers
        if i>1  
            elem(m+i-1,:) = [i,i+1]; % cortex
        end
        
        x = n(i+1,1);
        y = n(i+1,2);
        angf = rad2deg(atan(y/x));
%         if i==1 || i == m
        BC0s(end+1,:) = [i+1,-angf+90];
%         end
        
    end
    
    v0 = n(1,:);
    for i=1:m
            v1 = n(i+1,:);
            v = v1-v0;
            vn = norm(v);
            u = v./vn;      
            xm = (v0(2)-v1(2))/(v0(1)-v1(1)); c = v1(2) - xm*v1(1);
            n(m+i+1,:) = [-(Lx+Lg0) (xm*(-(Lx+Lg0))+c)];
            elem(2*m-1+i,:) = [i+1 m+i+1];
            BC0s(end+1,:) = [m+i+1,0];
            BC0s(end+1,:) = [m+i+1,90];
    end
    
    BC1s(:,1) = [2 m+1];
    BC1s(:,2) = [m+1 2*m-1];
    
    xr = min(n(:,1));
    yr = n(end,2);
    n(:,1) = n(:,1)-xr;
    n(:,2) = n(:,2)-yr;

end

function out = plot_lattice(n,elem,xj,m,fig, xout,ctot, z)
    % plot points, elements, and BCs
    
    [dim,~] = size(n);
    [dim2,~] = size(elem);   
    
    hfig1=figure(fig); set(hfig1,'color','w'); hold all
    haxsY=gca;
%     set(haxsY,'tickdir','out','linewidth',3,'fontsize',15,...
%         'FontName','Arial','FontWeight','bold')
    set(gca,'color','none')
    set(gca,'XColor', 'none','YColor','none')
    pbaspect([1 1 1])
    box on
    
    rang = 120;
    nrot1 = vrot2D(n', rang);
    nrot2 = n;
    nrot2(:,1) = -nrot2(:,1); 
    nrot3 = vrot2D(nrot2', rang);   
    nrot4 = vrot2D(nrot2', 2*rang);
    nrot5 = vrot2D(n', 2*rang);
    
    xlim([nrot2(1,1) n(1,1)])
    ylim([n(1,2) max(nrot1(2,:))]);
    
    xg = [58/255 150/255 58/255];
    xr = [102/255, 0, 0];
    for j=1:dim2-m
        m1 = elem(j,1);
        m2 = elem(j,2);
        x = [n(m1,1) n(m2,1)];
        y = [n(m1,2) n(m2,2)];
        plot(x,y,'color',xr,'linewidth',1.25);
        
        x = [nrot2(m1,1) nrot2(m2,1)];
        y = [nrot2(m1,2) nrot2(m2,2)];
        plot(x,y,'color',xr,'linewidth',1.25);
        
        x = [nrot1(1,m1) nrot1(1,m2)];
        y = [nrot1(2,m1) nrot1(2,m2)];
        plot(x,y,'color',xr,'linewidth',1.25);
        
        x = [nrot3(1,m1) nrot3(1,m2)];
        y = [nrot3(2,m1) nrot3(2,m2)];
        plot(x,y,'color',xr,'linewidth',1.25);
        
        x = [nrot4(1,m1) nrot4(1,m2)];
        y = [nrot4(2,m1) nrot4(2,m2)];
        plot(x,y,'color',xr,'linewidth',1.25);
        
        x = [nrot5(1,m1) nrot5(1,m2)];
        y = [nrot5(2,m1) nrot5(2,m2)];
        plot(x,y,'color',xr,'linewidth',1.25);
    end
    if (z==1)
        for i = 1:m-1
            if xj(i)==1
                xcb = xout(3*i-2)*ctot;
                al  = (xcb-5)/25; 

                m1 = elem(dim2-m+i,1);
                m2 = elem(dim2-m+i,2);
                x = [n(m1,1) n(m2,1)];
                y = [n(m1,2) n(m2,2)];
                p=plot(x,y,'color',xg,'linewidth',1.5); p.Color(4) = al;

                x = [nrot2(m1,1) nrot2(m2,1)];
                y = [nrot2(m1,2) nrot2(m2,2)];
                p=plot(x,y,'color',xg,'linewidth',1.5); p.Color(4) = al;

                x = [nrot1(1,m1) nrot1(1,m2)];
                y = [nrot1(2,m1) nrot1(2,m2)];
                p=plot(x,y,'color',xg,'linewidth',1.5); p.Color(4) = al;

                x = [nrot3(1,m1) nrot3(1,m2)];
                y = [nrot3(2,m1) nrot3(2,m2)];
                p=plot(x,y,'color',xg,'linewidth',1.5); p.Color(4) = al;

                x = [nrot4(1,m1) nrot4(1,m2)];
                y = [nrot4(2,m1) nrot4(2,m2)];
                p=plot(x,y,'color',xg,'linewidth',1.5); p.Color(4) = al;

                x = [nrot5(1,m1) nrot5(1,m2)];
                y = [nrot5(2,m1) nrot5(2,m2)];
                p=plot(x,y,'color',xg,'linewidth',1.5); p.Color(4) = al;
            end
        end
        i=21;
        if xj(i)==1       
                xcb = xout(3*i-2)*ctot;
                al  = (xcb-5)/25; 

                m1 = elem(dim2-m+i,1);
                m2 = elem(dim2-m+i,2);
                x = [n(m1,1) n(m2,1)];
                y = [n(m1,2) n(m2,2)];
                p=plot(x,y,'color',xg,'linewidth',1.5); p.Color(4) = al; 

                x = [nrot4(1,m1) nrot4(1,m2)];
                y = [nrot4(2,m1) nrot4(2,m2)];
                p=plot(x,y,'color',xg,'linewidth',1.5); p.Color(4) = al;

                x = [nrot5(1,m1) nrot5(1,m2)];
                y = [nrot5(2,m1) nrot5(2,m2)];
                p=plot(x,y,'color',xg,'linewidth',1.5); p.Color(4) = al;
        end
    end
    
    for i=1:dim-m
        msiz = 2;
        plot(n(i,1),n(i,2),'ko','markersize', msiz,'markerfacecolor','k');
        plot(nrot1(1,i),nrot1(2,i),'ko','markersize', msiz,'markerfacecolor','k');
        plot(nrot2(i,1),nrot2(i,2),'ko','markersize', msiz,'markerfacecolor','k');
        plot(nrot3(1,i),nrot3(2,i),'ko','markersize', msiz,'markerfacecolor','k');
        plot(nrot4(1,i),nrot4(2,i),'ko','markersize', msiz,'markerfacecolor','k');
        plot(nrot5(1,i),nrot5(2,i),'ko','markersize', msiz,'markerfacecolor','k');
    end
    

    out = 1;

end

function [n_new,eps] = quasi_stat_soln(kvec, n, elem, elem_L0, BC0s,BC1s,A0, xf,time)
    
% Calculate equilibrium solution given a system of points and reference state
    xs_thrs = 1e-6;
    F_thrs = (1e-4)*A0;
    iter_max = 100;
    
    [dim, ~] = size(n);
    const_rw = [];
    
    n_s = reshape(n',1,2*dim)';

    iter = 0;
    while iter < iter_max
        iter = iter+1;
        Jn_s = xJ(kvec,n_s,elem,elem_L0,A0, BC1s);
        [Fn_s,eps] = xFn(kvec,n_s,elem,elem_L0,A0, xf,BC1s);
        
        F_err = sum(Fn_s.^2);
        if F_err < F_thrs
            break;
        end
        
        [Jn_s, Fn_s] = apply_BC(Jn_s, Fn_s, BC0s);

        mdelta_n = Jn_s\Fn_s;
        
        n_s_i1 = n_s - mdelta_n; 

        n_dif = sum((-mdelta_n).^2);
        if n_dif < xs_thrs
            break;
        end
        
        n_s = n_s_i1;
    end

    if iter == iter_max
        X = ['Diverged at time=', num2str(time),'s'];
        error(X);
    end
    
    n_new = reshape(n_s,2,dim)';
        
end

function elem_L = elem_ref(n, elem)
    % Calculate element lengths
    
    % initialize arrays
    [dim, ~] = size(elem);
    elem_L   = zeros(dim,1);
      
    for ii = 1:dim
        i = elem(ii,1);
        j = elem(ii,2);
        
        elem_L(ii) = norm(n(i,:) - n(j,:));
    end
    
end

function [out,eps] = xFn(kvec,n,elem,elem_L0,A0,xf, BC1s)
    % Calculate the total force vector acting on each node
    % n - array containing all nodal coords
    % elem - all element numbers and connecting nodes
    % elem_L0 - element reference lengths
    
    dim  = length(n);
    Fn = zeros(dim/2,2);
    [dim2,~] = size(elem);
    m = length(xf);
    eps = zeros(dim2,1);

    for ii = 1:dim2
        i = elem(ii,1); j = elem(ii,2);
        if (ii <= m)
            f = xf(ii);
        else
            f = 0.0;
        end
        k = kvec(ii);
        [Fij,eps(ii)] = F_ij(k, n((2*i-1):(2*i)), n((2*j-1):(2*j)), elem_L0(ii),f);

        [~, sz] = size(BC1s);
        a = 1; b = 1;
        for jj = 1:sz
            if BC1s(jj,2) == ii
                if BC1s(jj,1) == i
                    a = 2;
                elseif BC1s(jj,1) == j
                    b = 2;
                end
            end
        end
        Fn(i,:) = Fn(i,:) + a*Fij';
        Fn(j,:) = Fn(j,:) - b*Fij';
        
    end

    out = reshape(Fn',dim,1);
    
end

function [out, eps] = F_ij(k, n_i, n_j, L0, f)
    % Calculate the force in the element connecting node i and node j
    % k - element stiffness
    % n - vector containing x and y coords
    % L0 - initial length of element connecting nodes
   
    z = 1e-8; % zero threshold
    
    n_ij = n_i - n_j;
    L_ij = norm(n_i - n_j); % current length between nodes
    eps  = (L_ij-L0)/L0;
    
    if abs(L_ij) < z
        out = 0.0./L_ij.*n_ij; 
    else
        out = (-k * (L_ij - L0)/L_ij).*n_ij + (-f/L_ij).*n_ij; % with force
    end
    
end

function Jn = xJ(kvec, n, elem, elem_L0,A0, BC1s)
    % Calculate Jacobian (stiffness matrix)
    
    dim  = length(n);
    Jn = zeros(dim,dim);
    [dim2,~] = size(elem);
    
    % loop over every element
    for ii = 1:dim2
        i = elem(ii,1); j = elem(ii,2);
        k = kvec(ii);
        dF = dF_ij_dn_i(k, n((2*i-1):(2*i)), n((2*j-1):(2*j)), elem_L0(ii));
        
        [~, sz] = size(BC1s);
        a = 1; b = 1;
        for jj = 1:sz
            if BC1s(jj,2) == ii
                if BC1s(jj,1) == i
                    a = 2;
                elseif BC1s(jj,1) == j
                    b = 2;
                end
            end
        end
        Jn((2*i-1):(2*i),(2*j-1):(2*j)) = -a*dF;
        Jn((2*i-1):(2*i),(2*i-1):(2*i)) = Jn((2*i-1):(2*i),(2*i-1):(2*i)) + a*dF;
        
        Jn((2*j-1):(2*j),(2*i-1):(2*i)) = -b*dF;
        Jn((2*j-1):(2*j),(2*j-1):(2*j)) = Jn((2*j-1):(2*j),(2*j-1):(2*j)) + b*dF;
        
    end
    

end

function out = dF_ij_dn_i(k, n_i, n_j, L0)
    % Calculate the derivative of the force with respect to change in node coords
    % 2x2 matrix
    
    z = 1e-8; % zero threshold
    
    n_ij = n_i - n_j;
    L_ij = norm(n_i - n_j); % current length between nodes
    n_ij_o = n_ij*n_ij';  
    
    if abs(L_ij) < z
        out = -k * eye(2); 
    else
        out = -k*eye(2) -k*L0.*(n_ij_o)./(L_ij^3) + k*L0.*(eye(2))./L_ij ;  
    end
    

end

function [Jout, Fout] = apply_BC(Jin, Fin, BC0s)

    Jout = Jin;
    Fout = Fin;
    
    [sz, ~] = size(BC0s);
    for ii = 1:sz
        ix = 2*BC0s(ii,1)-1;
        iy = 2*BC0s(ii,1);
        
        J01 = Jout(ix:iy,:); 
        [~,sz2] = size(J01);
        for jj = 1:sz2/2
            Jm0 = J01(:,(2*jj-1):2*jj);     
            Jm0r = mrot2D(Jm0,BC0s(ii,2));
            Jm0r(1,:) = 0.0;
            Jm0 = mrot2D(Jm0r,-BC0s(ii,2));
            J01(:,(2*jj-1):2*jj) = Jm0;
        end
        Jout(ix:iy,:) = J01;
        
        J02 = Jout(:,ix:iy); 
        for jj = 1:sz2/2
            Jm0 = J02((2*jj-1):2*jj,:);         
            Jm0r = mrot2D(Jm0,BC0s(ii,2));
            Jm0r(:,1) = 0.0;
            Jm0 = mrot2D(Jm0r,-BC0s(ii,2));
            J02((2*jj-1):2*jj,:) = Jm0;           
        end
        Jout(:,ix:iy) = J02;
        
        J0d = Jout(ix:iy,ix:iy);        
        J0dr = mrot2D(J0d,BC0s(ii,2));
        J0dr(1,1) = 1.0;
        J0d = mrot2D(J0dr,-BC0s(ii,2));
        Jout(ix:iy,ix:iy) = J0d;
        
        
        % forces 
        F0 = [Fout(ix);Fout(iy)];
        F0r = vrot2D(F0,BC0s(ii,2));
        F0r(1,1) = 0.0;
        F0 = vrot2D(F0r,-BC0s(ii,2));
        Fout(ix:iy) = F0; 
        
    end

end

function out = vrot2D(vec,ang)

    rads = deg2rad(ang);
    rotm = [cos(rads) -sin(rads); sin(rads) cos(rads)];
    out = rotm*vec;

end

function out = mrot2D(xmat,ang)

    rads = deg2rad(ang);
    rotm = [cos(rads) -sin(rads); sin(rads) cos(rads)];
    out = rotm*(xmat)*rotm';
end

%Events
function [elem_L0_new, n_new, xj_new] = events_update(t,x,kvec,n, elem, elem_L0, BC0s,BC1s,A0,...
    ctot,Lj0,kb,Ko,Fmax, xj,elem_L00,n0, m)
% Locate the time when cells come into contact
% and stop integration.
    
    % initialize
    elem_L0_new = elem_L0;
    n_new       = n;
    xj_new      = xj;
    
    Pb     = zeros(m,1);
    rho    = zeros(m,1);
    sigp   = zeros(m,1);
    xcb    = zeros(m,1);
    xf     = zeros(m,1);
    xKj    = zeros(m,1);
    dj     = zeros(m,1);
    Lx     = zeros(m,1);
    Fj     = zeros(m,1);
    
    for i=1:m
        Pb(i)   = x(3*i-2);  % no. bonds
        rho(i)  = x(3*i-1);  % contractility (myosin density)
        sigp(i) = x(3*i);  % polymerization (extension)

        xcb(i) = Pb(i)*ctot; % actual number of junction bonds
        xKj(i)   = Lj0*kb*xcb(i); % Effective junction stiffness = bonds x individual bond stiffness
    
        if xj(i) == 0
            kvec(end-m+i) = Ko*A0/elem_L0(end-m+i);
            xf(i) = (rho(i)+sigp(i))*A0; 
        elseif xj(i) == 1
            kvec(end-m+i) = xKj(i)*A0/elem_L0(end-m+i);
            xf(i)= (rho(i)+sigp(i))*A0; 
        end
    end
        
    [~, eps] = quasi_stat_soln(kvec,n, elem, elem_L0, BC0s,BC1s,A0, xf); 
    
    for i=1:m
        if xj(i) == 0
            Lx(i) = elem_L0(end-m+i)*(1+eps(end-m+i));
            % if protrusion extends enough, form junction

            if Lx(i) <= Lj0
                xj_new(i) = 1;
                elem_L0_new(end-m+i) = (Lj0); % form junction
                v0 = n(end-m+i,:); v1 = n(end-2*m+i,:); v = v1-v0;
                u = v./norm(v);
                n_new((end-2*m+i),:) = v0 + Lj0*u;
                X = ['junction ', num2str(i), ' form'];
                disp(X);
            else
                xj_new(i) = 0;
            end
                    
        elseif xj(i) == 1
    
            dj(i) = elem_L0(end-m+i)*(eps(end-m+i));
            Fj(i) = kb*dj(i);
        
            % if protrusion force exceeds crit, break junction
            if Fj(i) >= Fmax
                xj_new(i) = 0;
                elem_L0_new(end-m+i) = elem_L00(end-m+i); % no junction
                n_new((end-2*m+i),:) = n0(end-2*m+i,:);
                X = ['junction ', num2str(i), ' break'];
                disp(X);
            else
                xj_new(i) = 1;
            end
        
        end
    end

end 

function [value,isterminal,direction] = events1(t,x,kvec,n, elem, elem_L0, BC0s,BC1s,A0,...
    ctot,Lj0,kb,Ko,Fmax, xj, m,tmin)
% Locate the time when cells come into contact
% and stop integration.

    xk     = ones(2*length(xj),1);    
    Pb     = zeros(m,1);
    rho    = zeros(m,1);
    sigp   = zeros(m,1);
    xcb    = zeros(m,1);
    xf     = zeros(m,1);
    xKj    = zeros(m,1);
    dj     = zeros(m,1);
    Lx     = zeros(m,1);
    Fj     = zeros(m,1);
    
        for i=1:m
            Pb(i)   = x(3*i-2);  % no. bonds
            rho(i)  = x(3*i-1);  % contractility (myosin density)
            sigp(i) = x(3*i);  % polymerization (extension)

            xcb(i) = Pb(i)*ctot; % actual number of junction bonds
            xKj(i) = Lj0*kb*xcb(i); % Effective junction stiffness = bonds x individual bond stiffness

            if xj(i) == 0
                kvec(end-m+i) = Ko*A0/elem_L0(end-m+i);
                xf(i) = (rho(i)+sigp(i))*A0; 
            elseif xj(i) == 1
                kvec(end-m+i) = xKj(i)*A0/elem_L0(end-m+i);
                xf(i)= (rho(i)+sigp(i))*A0; 
            end
        end

        [~, eps] = quasi_stat_soln(kvec,n, elem, elem_L0, BC0s,BC1s,A0, xf, tmin); 

        for i=1:m
            if xj(i) == 0
                Lx(i) = elem_L0(end-m+i)*(1+eps(end-m+i));
                % if protrusion extends enough, form junction
                xk(2*i-1) = Lx(i) - Lj0;
                xk(2*i)   = Fmax - 0;

            elseif xj(i) == 1

                dj(i) = elem_L0(end-m+i)*(eps(end-m+i));
                Fj(i) = kb*dj(i);
                % if protrusion force exceeds crit, break junction
%                 if Fj(i) >= Fmax
%                     xk(i) = 0;
%                 else
%                     xk(i) = 1;
%                 end
                xk(2*i-1) = 1;
                xk(2*i)   = Fmax - Fj(i);

            end
        end

    value      = xk;   % detect gap exceeded
    isterminal = ones(2*length(xj),1);   % stop the integration
    direction  = -1.*ones(2*length(xj),1);   % direction dependent

end 

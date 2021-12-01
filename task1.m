clear

%% setup
load_constants
load dh_fun.mat
load TestTrack

options = optimoptions('fmincon', 'SpecifyConstraintGradient', true, ...
                        'SpecifyObjectiveGradient', true);

interp_const = 5;

[Xmin, Xmax, goals] = get_track_bounds_and_goals(TestTrack, interp_const);


num_itrs = 25;
ub = [repmat([1000; 1000; 1000; 1000; 1000; 1000], num_itrs, 1);repmat([0.5; 5000], (num_itrs-1), 1)];
lb = [repmat([-1000; -1000; -1000; -1000; -1000; -1000], num_itrs, 1);repmat([-0.5; -5000], (num_itrs-1), 1)];


%% main loop
x0 = rand(1,6*num_itrs+2*(num_itrs-1));
x0(1,1:6) = [287,5,-176,0,2,0];
for i = 1:length(goals.center)
    curr_goal_x = goals.center(1, i);
    curr_goal_y = goals.center(2, i);
    curr_goal_h = goals.heading(i);
    
    cf = @(z) costfun(z, curr_goal_x, curr_goal_y, curr_goal_h);
    nc = @(z) nonlcon(z, Xmin(i), Xmax(i), dhf);
    z = fmincon(cf, x0, [], [], [], [], lb' ,ub' ,nc, options);
    Y0 = reshape(z(1:6*num_itrs), 6, num_itrs)';
    U = reshape(z(6:num_itrs:end),2,num_itrs-1);
    x0 = Y0(end,:);
end





%% define model equations

    



%% function defs

function [Xmin, Xmax, goals] = get_track_bounds_and_goals(TestTrack, interp_const)
%     This function takes the given bounds of the racing track and convert them
%     to 
    
    lb_world = TestTrack.bl;
    ub_world = TestTrack.br;
    track_heading = TestTrack.theta;
    cline = TestTrack.cline;
    
    lb_world_new(1, :) = interp1(1:1:length(lb_world(1, :)), lb_world(1, :), 1:1/interp_const:length(lb_world(1, :)));
    lb_world_new(2, :) = interp1(1:1:length(lb_world(2, :)), lb_world(2, :), 1:1/interp_const:length(lb_world(2, :)));
    ub_world_new(1, :) = interp1(1:1:length(ub_world(1, :)), ub_world(1, :), 1:1/interp_const:length(ub_world(1, :)));
    ub_world_new(2, :) = interp1(1:1:length(ub_world(2, :)), ub_world(2, :), 1:1/interp_const:length(ub_world(2, :)));
    track_heading_new = interp1(1:1:length(track_heading), track_heading, 1:1/interp_const:length(track_heading));
    cline_new_x = interp1(1:1:length(cline), cline(1,:), 1:1/interp_const:length(cline));
    cline_new_y = interp1(1:1:length(cline), cline(2,:), 1:1/interp_const:length(cline));
    cline_new = [cline_new_x;cline_new_y];
    
    
    
    lb_car = rotmat(-track_heading_new)*reshape(lb_world_new, [], 1);
    ub_car = rotmat(-track_heading_new)*reshape(ub_world_new, [], 1);
    Xmin = reshape(lb_world_new, 2, []);
    Xmin = Xmin(1,:);
    
    Xmax = reshape(ub_world_new, 2, []);
    Xmax = Xmax(1,:);
    
    goals.center = cline_new;
    goals.heading = track_heading_new;
end


function [g,h,dg,dh]=nonlcon(z, ub_track, lb_track, dh_fun)
    num_itrs = 25;


    states = z(1:6*num_itrs);
    x = states(1:6:6*num_itrs);
    u = states(2:6:6*num_itrs);
    y = states(3:6:6*num_itrs);
    v = states(4:6:6*num_itrs);
    phi = states(5:6:6*num_itrs);
    r = states(6:6:6*num_itrs);
   
    inputs = z(6*num_itrs+1:end);
    delta = inputs(1:2:end);
    Fx = inputs(2:2:end);
    
    dx = zeros((num_itrs-1)*6,1);
    h = zeros(6*num_itrs,1);
    dg = zeros(2*num_itrs, 6*(num_itrs)+2*(num_itrs-1));
    g = zeros(2*num_itrs,1);
    dh = zeros(6*num_itrs, 6*(num_itrs)+2*(num_itrs-1));

    %% inequality constraint g
    idx = 1;
    for i = 1:2:2*num_itrs
        g(i) = x(idx)*cos(-phi(idx))-y(idx)*sin(-phi(idx))-ub_track;
        g(i+1) = -(x(idx)*cos(-phi(idx))-y(idx)*sin(-phi(idx)))-lb_track;
        idx = idx+1;
    end
 
    %% jacobian of inequality constraint g, dg
    idx = 1;
    for i = 1:2:2*num_itrs
        col_idx = (idx-1)*6+1:(idx-1)*6+6;
        
        dg(i, col_idx) = [ cos(phi(idx)), 0,  sin(phi(idx)), 0, y(idx)*cos(phi(idx)) - x(idx)*sin(phi(idx)), 0];
        dg(i+1, col_idx) = [-cos(phi(idx)), 0, -sin(phi(idx)), 0, x(idx)*sin(phi(idx)) - y(idx)*cos(phi(idx)), 0];
        idx = idx+1;
    end
    
    
    %% obtain dx
    for i=1:num_itrs-1
        row_idx = (i-1)*6+1:(i-1)*6+6;
        dx(row_idx, 1) = bike([x(i) u(i) y(i) v(i) phi(i) r(i)], [delta(i), Fx(i)]);
    end
    
    %% obtain h
    h(1:6,1) = [x(1); u(1); y(1); v(1); phi(1); r(1)];
    for i=2:num_itrs
        h((i-1)*6+1:(i-1)*6+6, 1) = [x(i) - x(i-1) - 0.01*dx((i-2)*6+1);
                                     u(i) - u(i-1) - 0.01*dx((i-2)*6+2);
                                     y(i) - y(i-1) - 0.01*dx((i-2)*6+3);
                                     v(i) - v(i-1) - 0.01*dx((i-2)*6+4);
                                     phi(i) - phi(i-1) - 0.01*dx((i-2)*6+5);
                                     r(i) - r(i-1) - 0.01*dx((i-2)*6+6)];
    end
    
   %% obtain dh
    
    dh = dh_fun(z);

     dg = dg';
     dh = dh';
end

function [J, dJ] = costfun(z, goal_x, goal_y, goal_h)

    num_itrs = 25;
    
    states = z(1:6*num_itrs);
    inputs = z(1:6*num_itrs+1:end);
    
    
    
    x = states(1:6:6*num_itrs);
    u = states(2:6:6*num_itrs);
    y = states(3:6:6*num_itrs);
    v = states(4:6:6*num_itrs);
    phi = states(5:6:6*num_itrs);
    r = states(6:6:6*num_itrs);
    
    delta = inputs(1:2:end);
    Fx = inputs(2:2:end);
    
    
    
    J = 0;
    for i=1:length(x)
        J = J + (x(i)-goal_x)^2 + (y(i)-goal_y)^2 + (phi(i)-goal_h)^2;
    end
    
    
    idx = 1;
    for i=1:6:num_itrs*6
        dJ(i) = 2*(x(idx)-goal_x);
        dJ(i+1) = 0;
        dJ(i+2) = 2*(y(idx)-goal_y);
        dJ(i+3) = 0;
        dJ(i+4) = 2*(phi(idx)-goal_h);
        dJ(i+5) = 0;
        idx = idx+1;
    end

    dJ = [dJ zeros(1, 98)];

    
end


function dzdt=bike(x,U)
%constants
Nw=2;
f=0.01;
Iz=2667;
a=1.35;
b=1.45;
By=0.27;
Cy=1.2;
Dy=0.7;
Ey=-1.6;
Shy=0;
Svy=0;
m=1400;
g=9.806;


%generate input functions
delta_f=U(1);
F_x=U(2);

%slip angle functions in degrees
a_f=rad2deg(delta_f-atan2(x(4)+a*x(6),x(2)));
a_r=rad2deg(-atan2((x(4)-b*x(6)),x(2)));

%Nonlinear Tire Dynamics
phi_yf=(1-Ey)*(a_f+Shy)+(Ey/By)*atan(By*(a_f+Shy));
phi_yr=(1-Ey)*(a_r+Shy)+(Ey/By)*atan(By*(a_r+Shy));

F_zf=b/(a+b)*m*g;
F_yf=F_zf*Dy*sin(Cy*atan(By*phi_yf))+Svy;

F_zr=a/(a+b)*m*g;
F_yr=F_zr*Dy*sin(Cy*atan(By*phi_yr))+Svy;

F_total=sqrt((Nw*F_x)^2+(F_yr^2));
F_max=0.7*m*g;

if F_total>F_max
    
    F_x=F_max/F_total*F_x;
  
    F_yr=F_max/F_total*F_yr;
end

%vehicle dynamics
dzdt= [x(2)*cos(x(5))-x(4)*sin(x(5));...
          (-f*m*g+Nw*F_x-F_yf*sin(delta_f))/m+x(4)*x(6);...
          x(2)*sin(x(5))+x(4)*cos(x(5));...
          (F_yf*cos(delta_f)+F_yr)/m-x(2)*x(6);...
          x(6);...
          (F_yf*a*cos(delta_f)-F_yr*b)/Iz];
end

function R = rotmat(h, sparse_flag)
% Given an angle in radians h, produce the 2-by-2 rotation matrix R that
% can rotate a point in \R^2 by the angle h about the origin
%
% If h is a vector (1-by-n or n-by-1), return a 2n-by-2n matrix in which
% each 2-by-2 block on the diagonal corresponds to each entry of h. This
% matrix is returned as a double if h has up to 100 entries, and a sparse
% matrix otherwise. Sparse output can also be forced by the argument
% sparse_flag.
    if nargin < 2
        sparse_flag = false ;
    end

    n = length(h) ;
    if n == 1
        R = [cos(h) -sin(h) ; sin(h) cos(h)] ;
    else
        h = h(:)' ;
        ch = cos(h) ; sh = sin(h) ;
        Rtall = [ch ; sh ; -sh ; ch] ;
        Rlong = reshape(Rtall,2,[]) ;
        if n <= 100 && ~sparse_flag
            % for smallish matrices, don't return sparse
            Rbig = repmat(Rlong, length(h),1) ;
            Icell = repmat({ones(2)},1,size(Rbig,2)/2,1) ;
            II = blkdiag(Icell{:}) ;
            R = Rbig.*II ;
        else
            io = 1:2:(2*n) ; % odd indices
            ie = 2:2:(2*n) ; % even indices
            
            % make row indices
            rio = repmat(io,2,1) ;
            rie = repmat(ie,2,1) ;
            rstack = [rio(:) rie(:)]' ;
            r = rstack(:) ;
            
            % make column indices
            c = repmat(1:2*n,2,1) ;
            c = c(:) ;
            
            % create sparse matrix
            R = sparse(r,c,Rlong(:)) ;
        end
    end
end

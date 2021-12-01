clear

load_constants
load dh_fun.mat
load TestTrack

interp_const = 5;
[Xmin, Xmax] = get_bounds(TestTrack, interp_const);



%% define model equations





%% function defs

function [Xmin, Xmax] = get_bounds(TestTrack, interp_const)
%     This function takes the given bounds of the racing track and convert them
%     to 
    
    lb_world = TestTrack.bl;
    ub_world = TestTrack.br;
    track_heading = TestTrack.theta;
    
    lb_world_new(1, :) = interp1(1:1:length(lb_world(1, :)), lb_world(1, :), 1:1/interp_const:length(lb_world(1, :)));
    lb_world_new(2, :) = interp1(1:1:length(lb_world(2, :)), lb_world(2, :), 1:1/interp_const:length(lb_world(2, :)));
    ub_world_new(1, :) = interp1(1:1:length(ub_world(1, :)), ub_world(1, :), 1:1/interp_const:length(ub_world(1, :)));
    ub_world_new(2, :) = interp1(1:1:length(ub_world(2, :)), ub_world(2, :), 1:1/interp_const:length(ub_world(2, :)));
    track_heading_new = interp1(1:1:length(track_heading), track_heading, 1:1/interp_const:length(track_heading));
    
    
    
    lb_car = rotmat(-track_heading_new)*reshape(lb_world_new, [], 1);
    ub_car = rotmat(-track_heading_new)*reshape(ub_world_new, [], 1);
    Xmin = reshape(lb_world_new, 2, []);
    Xmin = Xmin(1,:);
    
    Xmax = reshape(ub_world_new, 2, []);
    Xmax = Xmax(1,:);

end


function [g,h,dg,dh]=nonlcon(z, ub, lb)
    states = z(1:726);
    x = states(1:6:726);
    u = states(2:6:726);
    y = states(3:6:726);
    v = states(4:6:726);
    phi = states(5:6:726);
    r = states(6:6:726);
   
    inputs = z(727:end);
    delta = inputs(1:2:end);
    Fx = inputs(2:2:end);
    
    dx = zeros(120*6,1);
    h = zeros(726,1);
    dh = zeros(363, 966);
    dg = zeros(242, 966);
    g = zeros(242,1);
    b = 1.5 ; 
    L = 3 ;

    %% inequality constraint g
    idx = 1;
    for i = 1:2:242
        g(i) = x(idx)*cos(-phi(idx))-y(idx)*sin(-phi(idx))-ub;
        g(i+1) = -(x(idx)*cos(-phi(idx))-y(idx)*sin(-phi(idx)))-lb;
        idx = idx+1;
    end
 
    %% jacobian of inequality constraint g, dg
    idx = 1;
    for i = 1:2:242
        col_idx = (idx-1)*6+1:(idx-1)*6+6;
        
        dg(i, col_idx) = [ cos(h1), 0,  sin(h1), 0, y1*cos(h1) - x1*sin(h1), 0];
        dg(i+2, col_idx) = [-cos(h1), 0, -sin(h1), 0, x1*sin(h1) - y1*cos(h1), 0];
        idx = idx+1;
    end
    
    %% obtain dx
    for i=1:120
        row_idx = (i-1)*3+1:(i-1)*3+6;
        dx(row_idx, 1) = dzdt([x(i) u(i) y(i) v(i) phi(i) r(i)], [delta(i), Fx(i)]);
    end
    
    %% obtain h
    h(1:6,1) = [x(1); u(1); y(1); v(1); phi(1); r(1)];
    for i=2:121
        h((i-1)*6+1:(i-1)*6+6, 1) = [x((i-1)*6+1) - x((i-2)*6+1) - 0.01*dx((i-2)*6+1);
                                     x((i-1)*6+2) - x((i-2)*6+2) - 0.01*dx((i-2)*6+2);
                                     x((i-1)*6+3) - x((i-2)*6+3) - 0.01*dx((i-2)*6+3);
                                     x((i-1)*6+4) - x((i-2)*6+4) - 0.01*dx((i-2)*6+4);
                                     x((i-1)*6+5) - x((i-2)*6+5) - 0.01*dx((i-2)*6+5);
                                     x((i-1)*6+6) - x((i-2)*6+6) - 0.01*dx((i-2)*6+6)];
    end
    
   %% obtain dh
    
    dh = dh_fun(z);
    
    % size of g must be 121 x 1 (no.of time steps);
    % size of dg must be 603 x 121_fun = Transpose(no. of time steps x no. of elements in 'z');
    % size of h must be 363  1 ((no. of time steps * no. of states) x 1)
    % size of dh must be 603 x 363 = Transpose((no. of time steps * no. of states) x no. of elements in 'z') ;
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

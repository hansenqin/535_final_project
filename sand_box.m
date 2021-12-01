clear

%% setup
load_constants
load dh_fun.mat
load TestTrack
global num_itrs;
num_itrs = 25;

options = optimoptions('fmincon', 'SpecifyConstraintGradient', true, ...
                        'SpecifyObjectiveGradient', true);

interp_const = 5;

[Xmin, Xmax, goals] = get_track_bounds_and_goals(TestTrack, interp_const);


b = 1.5 ; 
L = 3 ;
dt=0.01;
nsteps = 121;

%remember the format for z is as follows:
%z=[x0 y0 th0 x1 y1 th1 ... xn yn thn u0 d0 ... u(n-1) d(n-1)]';
    
%3.1
ub = [repmat([8;3; 2*pi], 121,1);repmat([4; 0.5],120,1)];
lb = [repmat(-[1;3;-2*pi],121,1);repmat(-[0; 0.5],120,1)];

%3.4
%%%%%%%%%%%%%%%% no need to change these lines  %%%%%%%%%%%%%%%%%%%%%%
options = optimoptions('fmincon','SpecifyConstraintGradient',true,...
                       'SpecifyObjectiveGradient',true) ;
x0=zeros(1,5*nsteps-2);
cf=@costfun;
nc=@nonlcon;
z=fmincon(cf,x0,[],[],[],[],lb',ub',nc,options);

Y0=reshape(z(1:3*nsteps),3,nsteps)';
U=reshape(z(3*nsteps+1:end),2,nsteps-1);
u=@(t) [interp1(0:dt:119*dt,U(1,:),t,'previous','extrap');...
        interp1(0:dt:119*dt,U(2,:),t,'previous','extrap')];
[T1,Y1]=ode45(@(t,x) odefun(x,u(t)),[0:dt:120*dt],[0 0 0]);
[T2,Y2]=ode45(@(t,x) odefun(x,u(t)),[0:dt:120*dt],[0 0 -0.01]);
plot(Y0(:,1),Y0(:,2),Y1(:,1),Y1(:,2),Y2(:,1),Y2(:,2))
theta = 0:0.01:2*pi;
hold on
plot((0.7*cos(theta)+3.5),(0.7*sin(theta)-0.5))
hold on
plot(0,0,'x');
legend('fmincon trajectory','ode45 trajectory using x0 = [0;0;0]',...
    'ode45 trajectory using x0 = [0;0;-0.01]','Buffered Obstacle','Start');
ylim([-2,2]);
xlim([-1,8]);
xlabel('x');
ylabel('y');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%3.2
function [g,h,dg,dh]=nonlcon(z, lb_track, ub_track, x0)
    global num_itrs
    states = z(1:363);
    x = states(1:3:end);
    y = states(2:3:end);
    h = states(3:3:end);

    inputs = z(364:end);
    u = inputs(1:2:end);
    delta = inputs(2:2:end);
    
    dx = zeros(121*3,1);
    h = zeros(363,1);
    dh = zeros(603, 363);
    dg = zeros(603, 121);
    g = zeros(121,1);
    b = 1.5 ; 
    L = 3 ;

    idx = 1;
    for i = 1:2:2*num_itrs
        g(i) = x(idx)*cos(-h(idx))-y(idx)*sin(-h(idx))-ub_track;
        g(i+1) = -(x(idx)*cos(-h(idx))-y(idx)*sin(-h(idx)))+lb_track;
        idx = idx+1;
    end
    
    idx = 1;
    for i = 1:2:2*num_itrs
        col_idx = (idx-1)*6+1:(idx-1)*6+6;
        dg(i, col_idx) = [ cos(h(idx)), sin(h(idx)), y(idx)*cos(h(idx)) - x(idx)*sin(h(idx))];
        dg(i+1, col_idx) = [-cos(h(idx)), -sin(h(idx)), x(idx)*sin(h(idx)) - y(idx)*cos(h(idx))];
        idx = idx+1;
    end
    
    for i=1:120
        dx((i-1)*3+1:(i-1)*3+3, 1) = odefun([x(i), y(i), h(i)], [u(i), delta(i)]);
    end
    
    
    h(1:3,1) = [z(1); z(2); z(3)];
    for i=2:121
        h((i-1)*3+1:(i-1)*3+3, 1) = [x((i-1)*3+1) - x((i-2)*3+1) - 0.01*dx((i-2)*3+1);
                                     y((i-1)*3+2) - y((i-2)*3+2) - 0.01*dx((i-2)*3+2);
                                     h((i-1)*3+3) - h((i-2)*3+3) - 0.01*dx((i-2)*3+3)];
    end

    dh(1:3, 1:3) = eye(3);
    u_idx = 1;
    x_idx = 1;
    for i = 4:3:363
        dh(x_idx:x_idx+5, i) = [-1;
                                   0;
                                   0.01*(u(u_idx)*sin(x(x_idx+2))+b/L*u(u_idx)*tan(u(u_idx+1))*cos(x(x_idx+2)));
                                   1;
                                   0;
                                   0];
                               
        dh((u_idx)+363:(u_idx+1)+363, i) = [-0.01*(cos(x(x_idx+2))-b/L*tan(u(u_idx+1))*sin(x(x_idx+2)));
                                            -0.01*(-b/L*u(u_idx)*sec(u(u_idx+1))^2*sin(x(x_idx+2)))];
                                         
                                         
        dh(x_idx:x_idx+5, i+1) = [0;
                                   -1;
                                   -0.01*(u(u_idx)*cos(x(x_idx+2))-b/L*u(u_idx)*tan(u(u_idx+1))*sin(x(x_idx+2)));
                                    0;
                                    1
                                    0];
                               
        dh((u_idx)+363:(u_idx+1)+363, i+1) = [-0.01*(sin(x(x_idx+2))+b/L*tan(u(u_idx+1))*cos(x(x_idx+2)));
                                              -0.01*(b/L*u(u_idx)*sec(u(u_idx+1))^2*cos(x(x_idx+2)))];
                                         
        
        dh(x_idx:x_idx+5, i+2) = [0;
                                    0;
                                   -1;
                                    0;
                                    0;
                                    1];
                               
        dh((u_idx)+363:(u_idx+1)+363, i+2) = [-0.01/L*tan(u(u_idx+1));
                                             -0.01/L*u(u_idx)*sec(u(u_idx+1))^2];
                                         
                 
                               
       x_idx = x_idx+3;
       u_idx = u_idx+2;
    end
    
    % size of g must be 121 x 1 (no.of time steps);
    % size of dg must be 603 x 121 = Transpose(no. of time steps x no. of elements in 'z');
    % size of h must be 363  1 ((no. of time steps * no. of states) x 1)
    % size of dh must be 603 x 363 = Transpose((no. of time steps * no. of states) x no. of elements in 'z') ;
end

%3.3
function [J, dJ] = costfun(z)
  x = z(1:363);
    u = z(364:end);
    dJ = zeros(603,1);
    J1 = 0;
    J2 = 0;

   for i = 1:121
       J1 = J1 + (x((i-1)*3+1)-7)^2 + x((i-1)*3+2)^2+x((i-1)*3+3)^2;
   end
   
   for i = 1:120
       J2 = J2 + u((i-1)*2+1)^2+u((i-1)*2+2)^2;
   end
   
   J = J1+J2;
   
 
   for i=1:3:363
       dJ(i) = 2*x(i)-14;
       dJ(i+1) = 2*x(i+1);
       dJ(i+2) = 2*x(i+2);
   end
   
   for i=364:2:603
       dJ(i) = 2*u(i-363);
       dJ(i+1) = 2*u(i-363+1);
   end
   
   dJ = dJ';
   

end

function [dx] = odefun(x,u)
    b = 1.5 ; 
    L = 3 ;
    dx = [u(1)*cos(x(3))-(b/L)*u(1)*tan(u(2))*sin(x(3)) ;...
          u(1)*sin(x(3))+(b/L)*u(1)*tan(u(2))*cos(x(3)) ;...
          u(1)*tan(u(2))/L] ;
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
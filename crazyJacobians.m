clear
load_constants

num_itrs = 25;


X = sym('x', [1 num_itrs],'real');
U = sym('u', [1 num_itrs],'real');
Y = sym('y', [1 num_itrs],'real');
V = sym('v', [1 num_itrs],'real');
H = sym('h', [1 num_itrs],'real');
R = sym('r', [1 num_itrs],'real');

Delta = sym('delta', [1 num_itrs-1]);
Fx = sym('Fx', [1 num_itrs-1]);

Z1 = reshape([X;U;Y;V;H;R], [], 1);
Z2 = reshape([Delta;Fx], [], 1);

Z = [Z1;Z2];

%% calculating dg explicitly 
idx = 1;
% g = [];
for i = 1:2:2*num_itrs
    g_(i) = X(idx)*cos(-H(idx))-Y(idx)*sin(-H(idx));
    g_(i+1) = -(X(idx)*cos(-H(idx))-Y(idx)*sin(-H(idx)));
    idx = idx+1;
end
    
dg = jacobian(g_, Z);


%% calculating dh explicitly 
h = [X(1) U(1) Y(1) V(1) H(1) R(1)]';

idx = 2;
for i = 7:6:6*num_itrs

    %generate input functions
    delta_f=Delta(idx-1);
    F_x=Fx(idx-1);

    %slip angle functions in degrees
    a_f=rad2deg(delta_f-atan2(V(idx-1)+a*R(idx-1), U(idx-1)));
    a_r=rad2deg(-atan2((V(idx-1)-b*R(idx-1)), U(idx-1)));

    %Nonlinear Tire Dynamics
    phi_yf=(1-Ey)*(a_f+Shy)+(Ey/By)*atan(By*(a_f+Shy));
    phi_yr=(1-Ey)*(a_r+Shy)+(Ey/By)*atan(By*(a_r+Shy));

    F_zf=b/(a+b)*m*g;
    F_yf=F_zf*Dy*sin(Cy*atan(By*phi_yf))+Svy;

    F_zr=a/(a+b)*m*g;
    F_yr=F_zr*Dy*sin(Cy*atan(By*phi_yr))+Svy;

    F_total=sqrt((Nw*F_x)^2+(F_yr^2));
    F_max=0.7*m*g;
    
    X_dot = U(idx-1)*cos(H(idx-1))-V(idx-1)*sin(H(idx-1));
    U_dot = (-f*m*g+Nw*F_x-F_yf*sin(delta_f))/m+V(idx-1)*R(idx-1);
    Y_dot = U(idx-1)*sin(H(idx-1))+V(idx-1)*cos(H(idx-1));
    V_dot = (F_yf*cos(delta_f)+F_yr)/m-U(idx-1)*R(idx-1);
    H_dot = R(idx-1);
    R_dot = (F_yf*a*cos(delta_f)-F_yr*b)/Iz;
    
    h(i) = X(idx) - X(idx-1) - X_dot*0.02;
    h(i+1) = U(idx) - U(idx-1) - U_dot*0.02;
    h(i+2) = Y(idx) - Y(idx-1) - Y_dot*0.02;
    h(i+3) = V(idx) - V(idx-1) - V_dot*0.02;
    h(i+4) = H(idx) - H(idx-1) - H_dot*0.02;
    h(i+5) = R(idx) - R(idx-1) - R_dot*0.02;
    idx = idx+1;
end

dh = vpa(jacobian(h, Z),3);

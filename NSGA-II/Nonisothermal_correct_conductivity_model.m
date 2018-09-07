%% Non-isothermal CRTM constant velocity
% code for the non-isothermal filling of a circular mould using a
% constant velocity press mechanism
% It is assumed throughout the code that the part is thin enough to ignore
% through thickness flows. Temperature of the top and bottom moulds are
% assumed to be the same, and therefore symmetric (about the part
% mid-surface) temperature and cure profile develops.
%% Inputs---all in SI units
% clc
% tic
function [total_time_saved,dry_comp_force,inj_force,wet_comp_force]=Nonisothermal_correct_conductivity_model(h_dot_dry,Pinj,injection_h,h_dot,TEMPr,TEMPu)
% % % is_it_higher_level =0;
% Mould cavity inputs 
diameter = 1 ; 
radius = diameter/2;
inner_diameter = 20e-3;
inner_radius = inner_diameter/2;
initial_h = 1.0e-2;
injection_h =injection_h*(10^-2);
final_h =0.75e-2;

% process inputs (pre-curing phase)
h_dot_initial = h_dot_dry*(10^-3)/60; % to SI units
h_dot = h_dot*(10^-3)/60;
Pinj = Pinj*10^6;
TEMPl=TEMPu;% injection temperature of the resin
TEMPf=TEMPu; % initial temperature of the fibre mat same as patens

% fibre properties
reinforcement_layers =21; % number of layers ----> 21 for Vf = 0.5 and 15 for Vf = 0.35
areal_weight = 450; % g/m2 per fibre mat(layer)
densityf = 2580; % density of glass fibres
unstressed_Vf = 0.2;
unstressed_h = (areal_weight/1000)*reinforcement_layers/(densityf*unstressed_Vf);
Cpf=670; % heat capacity of fibreglass, obtained from lin et al.
Kf=0.67; % thermal conductivity taken from lin et al.
% constants for pereability computation
A = 3.07e-8;
B = -12.97;

% Resin properties (resin used is epoxy)
densityr=1087; % density of resin
Cpr=1260; % heat capacity of resin
Kr=0.168; % thermal conductivity of resin
deltaH=4.641*10^8;
Hr=deltaH/densityr; % heat of reaction

% constants for cure kinetics
A1=0.5963;
A2=57526.44;
E1=21514.78;
E2=49435.55;
m1=0.5874;
m2=3.2;
gasconst=8.314462; % universal gas constant 'R'

% constants for chemorheological model
Au=5.787e-11;
Eu=56492;
alphag=0.63; % gel point of the epoxy resin
a=0.4202;

% SUPG perturabation values
% gama = 2/3.873;
gama = 0;
perturbation = gama/2;

% Details for cure simulation
% % heat_rate_1 = 30; % K/min
% % heat_rate_1 = heat_rate_1/60; % K/sec
% % dwell_temperature = 273 + 153;
% % dwell_time = 0; % seconds
% % heat_rate_2 = 30; % K/min
% % heat_rate_2 = heat_rate_2/60; % K/sec
% % max_temperature_allowed = 273 + 227; % flash point of hardner
%% Axisymmetric grid generation
nz = 5; % number of nodes above/below mid-surface mesh for through thickness 
% heat conduction problem
nr =50; % number of radial elements 
mesh=ones(1,nr+1); % nr elements corresponds to nr+1 nodes
for i=2:nr+1
    mesh(i)=mesh(i-1)+1;
end
ele_length=(radius - inner_radius)/nr;
% define nodal radial coordinates
R=zeros(1,nr+1);
R(1)=inner_radius;
for i=2:nr+1
    R(i)=R(i-1)+ele_length;
end
Area=zeros(1,nr);
for i=1:nr
    Area(i)=pi*(R(i+1)^2-R(i)^2);
end
%% Declare variables that are constantly tracked
clamping_force = [];
injection_gate_pressure = [];
injection_time = 0;
wet_compaction_time = 0;
time = [];
total_time=0;
%% Dry compaction phase
steps_for_dc_phase = 10;
dry_compaction_time = (initial_h - injection_h)/h_dot_initial;
fillfactor = zeros(1,nr);
P = zeros(1,nr+1);
Vf = unstressed_Vf*unstressed_h/initial_h;
front = 1;
isstatic = 0; % is the mould static
clamping_force = [clamping_force force_computer(P , Vf , fillfactor, Area, front ,...
    R , isstatic)];
injection_gate_pressure = [injection_gate_pressure 0];
time = [time 0];
h = initial_h;
delta_t = (initial_h-injection_h)/h_dot_initial;
dt=delta_t/steps_for_dc_phase;
for i = 1:steps_for_dc_phase
    dh = h_dot_initial * dt;
    h=h-dh;
    Vf = unstressed_Vf*unstressed_h/h;
    total_time = total_time+dt;
    clamping_force = [clamping_force force_computer(P , Vf , fillfactor, Area, front ,...
        R , isstatic)];
    injection_gate_pressure = [injection_gate_pressure 0];
    time = [time total_time];
end
dry_comp_force = max(clamping_force);
% first relaxation phase of 0 second
isstatic = 1;
clamping_force = [clamping_force force_computer(P , Vf , fillfactor, Area, front ,...
        R , isstatic)];
injection_gate_pressure = [injection_gate_pressure 0];
total_time = total_time+0;
time = [time total_time];
%% RTM phase
h = injection_h;
total_resin_volume = pi*(radius^2 - inner_radius^2)*final_h*(1 - ...
    unstressed_h*unstressed_Vf/final_h);
fillfactor(1) = 0.52;
Volume = zeros(1,nr);
Volume(1)=0.52*Area(1)*h;
dVolume=zeros(1,nr);
Qin=zeros(1,nr);Qout=zeros(1,nr);

TEMP=TEMPf *ones(2*nz+1,nr+1); % temperature distribution matrix; initially assumed to be at fibre temperature
TEMP(2:2*nz,1)=TEMPr; % temperature at the inlet nodes is constant at the injection temperature
% of the resin
nodal_spacing=h/(2*nz+1); % spacing between the nodes in the through thickness direction

alpha=zeros(2*nz+1,nr+1); % declaration of matrix representing degree of cure at each node

viscosity_dist=zeros(2*nz+1,nr+1); % to take care of viscosity distribution across the thickness of the mould
for i=1:2*nz+1
    for j=1:2
        viscosity_dist(i,j)=Au*exp(Eu/(gasconst*TEMP(i,j)));
    end
end

viscosity=zeros(1,nr+1); % intialize gap averaged nodal viscosity
viscosity(1)=sum(viscosity_dist(:,1))/(2*nz+1);
viscosity(2)=sum(viscosity_dist(:,2))/(2*nz+1);

time_stored=total_time;
alpha_stored = zeros(2*nz+1,1);
velocity=zeros(1,nr+1); % initialization of nodal superficial velocities

phi = 1 - Vf; % porosity
rhoCp=phi*densityr*Cpr+(1-phi)*densityf*Cpf; % product of the effective density 
% and thermal heat capacity of the lumped fibre resin system

keff= Kr*Kf/(phi*Kf+(1-phi)*Kr); % effective stagnant thermal conductivity of the lumped
% resin fibre system

perm = A*exp(B*Vf);
injected_resin_volume = 0;
spot_temp=[];
front = 1;
inj_force = 0;
while (injected_resin_volume < total_resin_volume)&&(max(max(alpha)) < alphag)
    while (front<(nr+1))&&(fillfactor(front)>=0.5)
        front=front+1;
    end
    K=zeros(front,front);
    F=zeros(front,1);
    F(1)=Pinj;
    for i=1:front-1
        n1=mesh(i);
        n2=mesh(i+1);
        u1=viscosity(i);u2=viscosity(i+1); % retrieve nodal viscosity values
        u_avg = (u1 + u2)/2;
        r1 = R(i); r2 = R(i+1);
        % apply gaussian quadrature
        kef = (2*pi*h*perm/ele_length)*(r1+r2)/(2*u_avg);
        Ke=[kef -kef;-kef kef]; % element stiffness matrix
        if(i==1)
            Kestored=Ke;
        end
        K(n1,n1)=K(n1,n1)+Ke(1,1);
        K(n1,n2)=K(n1,n2)+Ke(1,2);
        K(n2,n1)=K(n2,n1)+Ke(2,1);
        K(n2,n2)=K(n2,n2)+Ke(2,2);
    end
    K(1,:)=0;K(front,:)=0;
    K(1,1)=1;
    K(front,front)=1;
    pressure=K\F;
    P=[pressure' zeros(1,nr+1-front)]; 
    
    discharge=Kestored(1,:)*[P(1) P(2)]';
    dt=Area(front-1)*h*phi/discharge; 
    if (injected_resin_volume + Area(1)*h*phi> total_resin_volume)
        dt = (total_resin_volume - injected_resin_volume)/discharge;
    end
    injection_time = injection_time + dt;
    total_time = total_time+dt;
    injected_resin_volume = injected_resin_volume + dt*discharge;
    
    
    Vout=discharge*dt;
    for i=1:nr
        Vin=Vout/phi; % Vin is the effective volume computed by factoring in the presence of fibres
        if (fillfactor(i)<1)
            emptyvolume=Area(i)*h*(1-fillfactor(i));
            if (Vin>=emptyvolume)
                fillfactor(i)=1;
                Vout=phi*(Vin-emptyvolume);
                Volume(i)=Area(i)*h;
                dVolume(i)=emptyvolume;
            end
            if (Vin<emptyvolume)
                fillfactor(i)=fillfactor(i)+Vin/(Area(i)*h);
                Vout=0;
                Volume(i)=Volume(i)+Vin;
                dVolume(i)=Vin;
            end
        else
            Vout=phi*Vin;
            Volume(i)=Area(i)*h;
            dVolume(i)=0;
        end
        Qin(i)=Vin*phi; Qout(i)=Vout;
    end
    
    % determination of the nodal superficial velocities
    for i=1:front
        velocity(i)=discharge/(2*pi*R(i)*h);
    end
    
    if (front>2)
            if (alpha(:,front)==zeros(2*nz+1,1))
                travel_time=total_time-time_stored;
                for i=1:nz+1
                    alpha(i,front)=runge_kutta(alpha_stored(i,1),TEMPf,travel_time);
                end
                for i = nz+2 : 2*nz + 1
                    alpha(i,front) = alpha(2*nz +2 -i,front);
                end
                time_stored=total_time;
                alpha_stored=alpha(:,front);
           end
    end 
    
     % apply crank nicolson method for through thickness heat conduction
    
        for j=1:front-1
            coeff_matrix=zeros(nz+2,nz+2);
            const_vector=zeros(nz+2,1);
            coeff_matrix(1,1)=1;const_vector(1)=TEMPu;
            for i=2:nz + 1
                Fo=(keff/rhoCp)*(dt/nodal_spacing^2)*0.5; % half of 
                % fourier's number
                coeff_matrix(i,i-1)=-Fo;
                coeff_matrix(i,i)=1+2*Fo;
                coeff_matrix(i,i+1)=-Fo;
                K1=A1*exp(-E1/(gasconst*(TEMP(i,j))));
                K2=A2*exp(-E2/(gasconst*(TEMP(i,j))));
                R_alpha=(K1+K2*alpha(i,j)^m1)*(1-alpha(i,j))^m2;
                omega=(R_alpha*phi*densityr/rhoCp)*dt*Hr;
                const_vector(i)=TEMP(i-1,j)*Fo+TEMP(i+1,j)*Fo+TEMP(i,j)*(1-2*Fo)+omega;
            end
            coeff_matrix(nz+2,nz+2) = 1;
            coeff_matrix(nz+2,nz) = -1;
            temp = coeff_matrix\const_vector; % updated temperature at given 
            % column due to through thickness heat conduction
            for i = nz+3 : 2*nz + 1
                temp(i) = temp(2*nz + 2 - i);
            end
            
            % place updated temperature vector into the existing TEMP
            % matrix
            
            TEMP(:,j)=temp;
        end    
        
        % apply fourth order Runge Kutta for curing at nodes

        for j=1:front
            for i=1:nz+1
                alpha(i,j)=runge_kutta(alpha(i,j),TEMP(i,j),dt);
            end
            for i = nz+2 : 2*nz+1
                alpha(i,j) = alpha(2*nz+2-i,j);
            end
        end
    
        % apply streamline upwind petrov galerkin method for solving the
        % in-plane convection-diffusion equation
        
        for i=2:nz+1
            coeff_matrix=zeros(front,front);
            const_vector=zeros(front,1);
            coeff_matrix(1,1)=1;const_vector(1)=TEMPr;
            coeff_matrix(front,front)=1;const_vector(front)=TEMPf;
            for j=2:front-1
                              
                % solve for the integrals containing the spacial
                % derivatives
                
                % first solve for coefficient of (Tj-Tj-1)

                % integrate using gaussian quadrature
                
                cf1 = 0;
                r1 = R(j-1); r2 = R(j);
                for t = -0.774597 : 0.774597 : 0.774597 % 3 sample points are chosen as the polynomial is of degree 4.
                    if (t == 0)
                        gauss_weight = 8/9;
                    else
                        gauss_weight = 5/9;
                    end
                    w = (((t*ele_length + r1+r2)/2)-r1)/ele_length;
                    v = (velocity(j-1)*(r2-((t*ele_length + r1+r2)/2))/ele_length)...
                        +(velocity(j)*(((t*ele_length + r1+r2)/2)-r1)/ele_length);
                    cf1 = cf1 + gauss_weight *(w * Cpr * densityr * v ...
                        * ((t*ele_length + r1+r2)/2)^2 + ((t*ele_length + r1+r2)/2)^2*...
                        (1/ele_length)*(densityr*Cpr*v*gama*ele_length/2+keff))/2;
                end
                    
                % next solve for coefficient of (Tj+1-Tj)

                cf2 = 0;
                r1 = R(j) ; r2 = R(j+1);
                for t = -0.774597 : 0.774597 : 0.774597 % 3 sample points are chosen as the polynomial is of degree 4.
                    if (t == 0)
                        gauss_weight = 8/9;
                    else
                        gauss_weight = 5/9;
                    end
                    w = (r2 - ((t*ele_length + r1+r2)/2))/ele_length;
                    v = (velocity(j)*(r2-((t*ele_length + r1+r2)/2))/ele_length)...
                        +(velocity(j+1)*(((t*ele_length + r1+r2)/2)-r1)/ele_length);
                    cf2 = cf2 + gauss_weight *(w * Cpr * densityr * v ...
                        * ((t*ele_length + r1+r2)/2)^2 + ((t*ele_length + r1+r2)/2)^2*...
                        (-1/ele_length)*(densityr*Cpr*v*gama*ele_length/2+keff))/2;
                end
                    
                    
                % solve for integral containing the time derivative
                
                % first solve for coefficient of dTj-1/dt and part of the
                % coefficient of dTj/dt
                
                cf_previous = 0;
                cf_j = 0;
                r1 = R(j-1) ; r2 = R(j);
                for t = -0.774597 : 0.774597 : 0.774597 % 3 sample points are chosen as the polynomial is of degree 4.
                    if (t == 0)
                        gauss_weight = 8/9;
                    else
                        gauss_weight = 5/9;
                    end
                    w = (((t*ele_length + r1+r2)/2)-r1)/ele_length;
                    N2 = w;
                    N1 = (r2 - ((t*ele_length + r1+r2)/2))/ele_length;
                    cf_previous = cf_previous + gauss_weight*(w+perturbation)*rhoCp*...
                        ((t*ele_length + r1+r2)/2)^2*N1*ele_length/2;
                    cf_j = cf_j + gauss_weight*(w+perturbation)*rhoCp*((t*ele_length + r1+r2)/2)^2 ...
                        *N2*ele_length/2;
                end
                
                
                % next solve for the cofficient of dTj+1/dt and remaining
                % part of the coefficient of dTj/dt
                
                cf_succesive = 0;
                r1 = R(j); r2 = R(j+1);
                for t = -0.774597 : 0.774597 : 0.774597 % 3 sample points are chosen as the polynomial is of degree 4.
                    if (t == 0)
                        gauss_weight = 8/9;
                    else
                        gauss_weight = 5/9;
                    end
                    w = (r2 - ((t*ele_length + r1+r2)/2))/ele_length;
                    N2 = (((t*ele_length + r1+r2)/2)-r1)/ele_length;
                    N1 = w;
                    cf_succesive = cf_succesive + gauss_weight*(w-perturbation)*rhoCp*...
                        ((t*ele_length + r1+r2)/2)^2*N2*ele_length/2;
                    cf_j = cf_j + gauss_weight*(w-perturbation)*rhoCp*((t*ele_length + r1+r2)/2)^2 ...
                        *N1*ele_length/2;
                end
                 
                % generate row of the coefficient matrix and the constant
                % vector corresponding to the node currently under
                % consideration. first order accurate backward difference
                % for time discretization is done to reduce computer storage.
                
                coeff_matrix(j,j-1)=cf_previous/dt-cf1;
                coeff_matrix(j,j)=cf_j/dt+cf1-cf2;
                coeff_matrix(j,j+1)=cf_succesive/dt+cf2;
                const_vector(j)=(cf_previous*TEMP(i,j-1)/dt)+(cf_j*TEMP(i,j)/dt)+(cf_succesive*TEMP(i,j+1)/dt);
            end
            temp=coeff_matrix\const_vector; %final updated temperature at the
            % particular time step for the current row of nodes
            
            % place updated temperature vector into the existing TEMP
            % matrix
            
            TEMP(i,1:front)=temp';
        end
        for i=nz+2:2*nz+1
            TEMP(i,:)=TEMP(2*nz+2-i,:);
        end
        
        % apply streamline upwind petrov galerkin method for solution of
        % the species transport equation
        
        for i=1:nz+1
            coeff_matrix=zeros(front,front);
            const_vector=zeros(front,1);
            coeff_matrix(1,1)=1;const_vector(1)=0;
            for j=2:front
                % solve for integrals containing the spatial derivative of
                % the degree of cure
                
                % first solve for the coefficient of alphaj-alphaj-1

                r1 = R(j-1); r2 = R(j);
                cf1 = 0;
                for t = -0.57735 : 2*0.57735 : 0.57735
                    gauss_weight = 1;
                    w=(((t*ele_length + r1+r2)/2)-r1)/ele_length;
                    v = (velocity(j-1)*(r2-((t*ele_length + r1+r2)/2))/ele_length)...
                        +(velocity(j)*(((t*ele_length + r1+r2)/2)-r1)/ele_length);
                    cf1 = cf1+ gauss_weight*(w+perturbation)*((t*ele_length + r1+r2)/2)*v/2;
                end
                
                if (j~=front)
                    % next solve for coefficient of alphaj+1-alphaj
                    
                    r1 = R(j); r2 = R(j+1);
                    cf2 = 0;
                    for t = -0.57735 : 2*0.57735 : 0.57735
                        gauss_weight = 1;
                        w = (r2 - ((t*ele_length + r1+r2)/2))/ele_length;
                        v = (velocity(j)*(r2-((t*ele_length + r1+r2)/2))/ele_length)...
                            +(velocity(j+1)*(((t*ele_length + r1+r2)/2)-r1)/ele_length);
                        cf2 = cf2+ gauss_weight*(w-perturbation)*((t*ele_length + r1+r2)/2)*v/2;
                    end                    
                else
                    cf2=0;
                end
                
                % solve integrals containing time derivative of the degree
                % of cure
                
                % fisrt solve for coefficient of d alphaj-1/dt and part of
                % the coefficient of d alphaj/dt 
             
                r1 = R(j-1); r2 = R(j);
                cf_previous = 0;
                cf_j = 0;
                for t = -0.57735 : 2*0.57735 : 0.57735
                    gauss_weight = 1;
                    w = (((t*ele_length + r1+r2)/2)-r1)/ele_length;
                    N2 = w;
                    N1 = (r2 - ((t*ele_length + r1+r2)/2))/ele_length;
                    cf_previous = cf_previous + gauss_weight*(w+perturbation)*phi ...
                        *((t*ele_length + r1+r2)/2)*N1*ele_length/2;
                    cf_j = cf_j + gauss_weight*(w+perturbation)*phi*((t*ele_length + r1+r2)/2)...
                        *N2*ele_length/2;
                end
                
                if (j~=front)
                    % next solve for coefficient of d alphaj+1/dt and the
                    % remaining part of d alphaj/dt
                
                      r1 = R(j); r2 = R(j+1);
                      cf_succesive = 0;
                     for t = -0.57735 : 2*0.57735 : 0.57735
                        gauss_weight = 1;
                        w = (r2 - ((t*ele_length + r1+r2)/2))/ele_length;
                        N2 = (((t*ele_length + r1+r2)/2)-r1)/ele_length;
                        N1 = w;
                        cf_succesive = cf_succesive + gauss_weight*(w-perturbation)*phi ...
                            *((t*ele_length + r1+r2)/2)*N2*ele_length/2;
                        cf_j = cf_j + gauss_weight*(w-perturbation)*phi*((t*ele_length + r1+r2)/2)...
                            *N1*ele_length/2;
                    end                      
                else

                    cf_succesive=0;
                    cf_j=cf_j+0;
                end
                
                % generate row of the coefficient matrix and the constant
                % vector corresponding to the node currently under
                % consideration. first order accurate backward difference
                % for time discretization is done to reduce computer storage.
                coeff_matrix(j,j-1)=cf_previous/dt-cf1;
                coeff_matrix(j,j)=cf_j/dt+cf1-cf2;
                if (j~=front)
                    coeff_matrix(j,j+1)=cf_succesive/dt+cf2;
                end
                if (j~=front)
                    const_vector(j)=(cf_previous*alpha(i,j-1)/dt)+(cf_j*alpha(i,j)/dt)+(cf_succesive*alpha(i,j+1)/dt);
                else
                    const_vector(j)=(cf_previous*alpha(i,j-1)/dt)+(cf_j*alpha(i,j)/dt); 
                end
            end
            ALPHA=coeff_matrix\const_vector; %final updated degree of cure at the
            % particular time step for the current row of nodes
            
            % place updated degree of cure vector into the existing alpha
            % matrix
            
            alpha(i,1:front)=ALPHA';
            if max(ALPHA) >= alphag
                break;
            end
        end
    for i=nz+2:2*nz+1
        alpha(i,:)=alpha(2*nz+2-i,:);
    end     
    
    % Deduce the viscosity distribution using the current temperature and
    % degree of cure distributions obtained
    
    for i=1:2*nz+1
        for j=1:front
            b = 49.97 - 0.3278*TEMP(i,j) +  0.00055*(TEMP(i,j))^2;
            viscosity_dist(i,j)=Au*exp(Eu/(gasconst*TEMP(i,j)))*(alphag/(alphag-alpha(i,j)))^(a+b*alpha(i,j));
        end
    end
    
    % obtain the thickness averaged viscosity from the through thickness
    % viscosity obtained
    
    for i=1:front
        viscosity(i)=sum(viscosity_dist(:,(i)))/(2*nz+1);
    end
    
    if (front~=nr+1)
        viscosity(front+1)=viscosity(front); % In case the boundary condition
        %changes in the next flow time step
    end   
    
    isstatic = 1;
    current_force = force_computer(P , Vf , fillfactor, Area, front ,...
        R , isstatic);
    if current_force > inj_force
        inj_force = current_force;
    end
    clamping_force = [clamping_force current_force];
    injection_gate_pressure = [injection_gate_pressure Pinj];
    time = [time total_time];
    spot_temp = [spot_temp TEMP(nz+1,20)];
end
% plot(R,P)
% second relaxation phase of 0 second
isstatic = 1;
P = zeros (1,nr+1);
total_time = total_time + 0;
clamping_force = [clamping_force force_computer(P , Vf , fillfactor, Area, front ,...
        R , isstatic)];
injection_gate_pressure = [injection_gate_pressure 0];
time = [time total_time];
%% Compression RTM phase starts
% due to the variable cavity thickness in this phase the fillfactor has an
% areal meaning instead of a volumetric one
front_radius = sqrt((total_resin_volume/(pi*h*phi))+inner_radius^2);
front = 1;
wet_comp_force=0;
while (h > final_h)&&(max(max(alpha)) < alphag) 
    while (front<(nr+1))&&(fillfactor(front)>=0.5)
        front=front+1;
    end
    K=zeros(front,front);
    F=zeros(front,1);
    for i=1:front-1
        n1=mesh(i);
        n2=mesh(i+1);
        u1=viscosity(i);u2=viscosity(i+1); % retrieve nodal viscosity values
        u_avg = (u1 + u2)/2;
        r1 = R(i); r2 = R(i+1);
        kef = (2*pi*h*perm/ele_length)*(r1+r2)/(2*u_avg);
        fef1 = 0; fef2 = 0;        
        for t = -0.57735 : 2*0.57735 : 0.57735
            gauss_weight = 1;
            fef1 = fef1 + gauss_weight*2*pi*h_dot*((t*ele_length + r1+r2)/2)...
                *(r2 - ((t*ele_length + r1+r2)/2))/2;
            fef2 = fef2 + gauss_weight*2*pi*h_dot*((t*ele_length + r1+r2)/2)...
                *(((t*ele_length + r1+r2)/2) - r1)/2;
        end
        Ke=[-kef kef;kef -kef]; % element stiffness matrix
        Fe = [fef1;fef2];
        if(i==front-1)
            Kestored=Ke;
            Festored = Fe;
        end
        K(n1,n1)=K(n1,n1)+Ke(1,1);
        K(n1,n2)=K(n1,n2)+Ke(1,2);
        K(n2,n1)=K(n2,n1)+Ke(2,1);
        K(n2,n2)=K(n2,n2)+Ke(2,2);
        F(n1:n2) = F(n1:n2) - Fe;
    end
    K(front,:)=0;
    K(front,front)=1;
    F(front) = 0;
    pressure=K\F;
    P=[pressure' zeros(1,nr+1-front)]; 
    
    front_discharge=Kestored(2,:)*[P(front-1) 0]'+Festored(2);
    front_velocity = front_discharge/(2*pi*R(front)*h);
    dt=ele_length*phi/(2*front_velocity); 
    if (h - h_dot*dt < final_h)
        dt = (h - final_h)/h_dot;
    end
    wet_compaction_time = wet_compaction_time + dt;
    total_time = total_time+dt;
    
    % determination of the nodal superficial velocities---specific formula
    % derived from mass conservation principles for a circular mould
    for i=1:front
        velocity(i)=h_dot*(R(i)^2-R(1)^2)/(2*h*R(i));
    end
    
    h = h - h_dot*dt;
    Vf = unstressed_Vf*unstressed_h/h;
    phi = (1-Vf);
    perm = A*exp(B*Vf);
    
    front_radius = front_radius + front_velocity*dt/phi;
    for i = 1:nr
        if (R(i) < front_radius) && (R(i+1) <= front_radius)
            fillfactor(i) = 1;
        elseif (R(i) < front_radius) && (R(i+1) > front_radius)
            fillfactor(i) = (front_radius - R(i))/(R(i+1) - R(i));
        else
            fillfactor(i) = 0;
        end
    end   
    
    nodal_spacing = h/(2*nz+1);
    
    % define effective material properties for lumped system
    rhoCp=phi*densityr*Cpr+(1-phi)*densityf*Cpf; 
    keff= Kr*Kf/(phi*Kf+(1-phi)*Kr); 
    
    if (front>2)
            if (alpha(:,front)==zeros(2*nz+1,1))
                travel_time=total_time-time_stored;
                for i=1:nz+1
                    alpha(i,front)=runge_kutta(alpha_stored(i,1),TEMPf,travel_time);
                end
                for i = nz+2 : 2*nz + 1
                    alpha(i,front) = alpha(2*nz +2 -i,front);
                end
                time_stored=total_time;
                alpha_stored=alpha(:,front);
           end
    end 
    
     % apply crank nicolson method for through thickness heat conduction
    
        for j=1:front-1
            coeff_matrix=zeros(nz+2,nz+2);
            const_vector=zeros(nz+2,1);
            coeff_matrix(1,1)=1;const_vector(1)=TEMPu;
            for i=2:nz + 1
                Fo=(keff/rhoCp)*(dt/nodal_spacing^2)*0.5; % half of 
                % fourier's number
                coeff_matrix(i,i-1)=-Fo;
                coeff_matrix(i,i)=1+2*Fo;
                coeff_matrix(i,i+1)=-Fo;
                K1=A1*exp(-E1/(gasconst*(TEMP(i,j))));
                K2=A2*exp(-E2/(gasconst*(TEMP(i,j))));
                R_alpha=(K1+K2*alpha(i,j)^m1)*(1-alpha(i,j))^m2;
                omega=(R_alpha*phi*densityr/rhoCp)*dt*Hr;
                const_vector(i)=TEMP(i-1,j)*Fo+TEMP(i+1,j)*Fo+TEMP(i,j)*(1-2*Fo)+omega;
            end
            coeff_matrix(nz+2,nz+2) = 1;
            coeff_matrix(nz+2,nz) = -1;
            temp = coeff_matrix\const_vector; % updated temperature at given 
            % column due to through thickness heat conduction
            for i = nz+3 : 2*nz + 1
                temp(i) = temp(2*nz + 2 - i);
            end
            
            % place updated temperature vector into the existing TEMP
            % matrix
            
            TEMP(:,j)=temp;
        end    
        
        % apply fourth order Runge Kutta for curing at nodes

        for j=1:front
            for i=1:nz+1
                alpha(i,j)=runge_kutta(alpha(i,j),TEMP(i,j),dt);
            end
            for i = nz+2 : 2*nz+1
                alpha(i,j) = alpha(2*nz+2-i,j);
            end
        end
    
        % apply streamline upwind petrov galerkin method for solving the
        % in-plane convection-diffusion equation
        
        for i=2:nz+1
            coeff_matrix=zeros(front,front);
            const_vector=zeros(front,1);
            coeff_matrix(1,1)=1;const_vector(1)=TEMP(i,1);
            coeff_matrix(front,front)=1;const_vector(front)=TEMPf;
            for j=2:front-1
                              
                % solve for the integrals containing the spacial
                % derivatives
                
                % first solve for coefficient of (Tj-Tj-1)

                % integrate using gaussian quadrature
                
                cf1 = 0;
                r1 = R(j-1); r2 = R(j);
                for t = -0.774597 : 0.774597 : 0.774597 % 3 sample points are chosen as the polynomial is of degree 4.
                    if (t == 0)
                        gauss_weight = 8/9;
                    else
                        gauss_weight = 5/9;
                    end
                    w = (((t*ele_length + r1+r2)/2)-r1)/ele_length;
                    v = (velocity(j-1)*(r2-((t*ele_length + r1+r2)/2))/ele_length)...
                        +(velocity(j)*(((t*ele_length + r1+r2)/2)-r1)/ele_length);
                    cf1 = cf1 + gauss_weight *(w * Cpr * densityr * v ...
                        * ((t*ele_length + r1+r2)/2)^2 + ((t*ele_length + r1+r2)/2)^2*...
                        (1/ele_length)*(densityr*Cpr*v*gama*ele_length/2+keff))/2;
                end
                    
                % next solve for coefficient of (Tj+1-Tj)

                cf2 = 0;
                r1 = R(j) ; r2 = R(j+1);
                for t = -0.774597 : 0.774597 : 0.774597 % 3 sample points are chosen as the polynomial is of degree 4.
                    if (t == 0)
                        gauss_weight = 8/9;
                    else
                        gauss_weight = 5/9;
                    end
                    w = (r2 - ((t*ele_length + r1+r2)/2))/ele_length;
                    v = (velocity(j)*(r2-((t*ele_length + r1+r2)/2))/ele_length)...
                        +(velocity(j+1)*(((t*ele_length + r1+r2)/2)-r1)/ele_length);
                    cf2 = cf2 + gauss_weight *(w * Cpr * densityr * v ...
                        * ((t*ele_length + r1+r2)/2)^2 + ((t*ele_length + r1+r2)/2)^2*...
                        (-1/ele_length)*(densityr*Cpr*v*gama*ele_length/2+keff))/2;
                end
                    
                    
                % solve for integral containing the time derivative
                
                % first solve for coefficient of dTj-1/dt and part of the
                % coefficient of dTj/dt
                
                cf_previous = 0;
                cf_j = 0;
                r1 = R(j-1) ; r2 = R(j);
                for t = -0.774597 : 0.774597 : 0.774597 % 3 sample points are chosen as the polynomial is of degree 4.
                    if (t == 0)
                        gauss_weight = 8/9;
                    else
                        gauss_weight = 5/9;
                    end
                    w = (((t*ele_length + r1+r2)/2)-r1)/ele_length;
                    N2 = w;
                    N1 = (r2 - ((t*ele_length + r1+r2)/2))/ele_length;
                    cf_previous = cf_previous + gauss_weight*(w+perturbation)*rhoCp*...
                        ((t*ele_length + r1+r2)/2)^2*N1*ele_length/2;
                    cf_j = cf_j + gauss_weight*(w+perturbation)*rhoCp*((t*ele_length + r1+r2)/2)^2 ...
                        *N2*ele_length/2;
                end
                
                
                % next solve for the cofficient of dTj+1/dt and remaining
                % part of the coefficient of dTj/dt
                
                cf_succesive = 0;
                r1 = R(j); r2 = R(j+1);
                for t = -0.774597 : 0.774597 : 0.774597 % 3 sample points are chosen as the polynomial is of degree 4.
                    if (t == 0)
                        gauss_weight = 8/9;
                    else
                        gauss_weight = 5/9;
                    end
                    w = (r2 - ((t*ele_length + r1+r2)/2))/ele_length;
                    N2 = (((t*ele_length + r1+r2)/2)-r1)/ele_length;
                    N1 = w;
                    cf_succesive = cf_succesive + gauss_weight*(w-perturbation)*rhoCp*...
                        ((t*ele_length + r1+r2)/2)^2*N2*ele_length/2;
                    cf_j = cf_j + gauss_weight*(w-perturbation)*rhoCp*((t*ele_length + r1+r2)/2)^2 ...
                        *N1*ele_length/2;
                end
                 
                % generate row of the coefficient matrix and the constant
                % vector corresponding to the node currently under
                % consideration. first order accurate backward difference
                % for time discretization is done to reduce computer storage.
                
                coeff_matrix(j,j-1)=cf_previous/dt-cf1;
                coeff_matrix(j,j)=cf_j/dt+cf1-cf2;
                coeff_matrix(j,j+1)=cf_succesive/dt+cf2;
                const_vector(j)=(cf_previous*TEMP(i,j-1)/dt)+(cf_j*TEMP(i,j)/dt)+(cf_succesive*TEMP(i,j+1)/dt);
            end
            temp=coeff_matrix\const_vector; %final updated temperature at the
            % particular time step for the current row of nodes
            
            % place updated temperature vector into the existing TEMP
            % matrix
            
            TEMP(i,1:front)=temp';
        end
        for i=nz+2:2*nz+1
            TEMP(i,:)=TEMP(2*nz+2-i,:);
        end
        
        % apply streamline upwind petrov galerkin method for solution of
        % the species transport equation
        
        for i=1:nz+1
            coeff_matrix=zeros(front,front);
            const_vector=zeros(front,1);
            coeff_matrix(1,1)=1;const_vector(1)=alpha(i,1);
            for j=2:front
                % solve for integrals containing the spatial derivative of
                % the degree of cure
                
                % first solve for the coefficient of alphaj-alphaj-1

                r1 = R(j-1); r2 = R(j);
                cf1 = 0;
                for t = -0.57735 : 2*0.57735 : 0.57735
                    gauss_weight = 1;
                    w=(((t*ele_length + r1+r2)/2)-r1)/ele_length;
                    v = (velocity(j-1)*(r2-((t*ele_length + r1+r2)/2))/ele_length)...
                        +(velocity(j)*(((t*ele_length + r1+r2)/2)-r1)/ele_length);
                    cf1 = cf1+ gauss_weight*(w+perturbation)*((t*ele_length + r1+r2)/2)*v/2;
                end
                
                if (j~=front)
                    % next solve for coefficient of alphaj+1-alphaj
                    
                    r1 = R(j); r2 = R(j+1);
                    cf2 = 0;
                    for t = -0.57735 : 2*0.57735 : 0.57735
                        gauss_weight = 1;
                        w = (r2 - ((t*ele_length + r1+r2)/2))/ele_length;
                        v = (velocity(j)*(r2-((t*ele_length + r1+r2)/2))/ele_length)...
                            +(velocity(j+1)*(((t*ele_length + r1+r2)/2)-r1)/ele_length);
                        cf2 = cf2+ gauss_weight*(w-perturbation)*((t*ele_length + r1+r2)/2)*v/2;
                    end                    
                else
                    cf2=0;
                end
                
                % solve integrals containing time derivative of the degree
                % of cure
                
                % fisrt solve for coefficient of d alphaj-1/dt and part of
                % the coefficient of d alphaj/dt 
             
                r1 = R(j-1); r2 = R(j);
                cf_previous = 0;
                cf_j = 0;
                for t = -0.57735 : 2*0.57735 : 0.57735
                    gauss_weight = 1;
                    w = (((t*ele_length + r1+r2)/2)-r1)/ele_length;
                    N2 = w;
                    N1 = (r2 - ((t*ele_length + r1+r2)/2))/ele_length;
                    cf_previous = cf_previous + gauss_weight*(w+perturbation)*phi ...
                        *((t*ele_length + r1+r2)/2)*N1*ele_length/2;
                    cf_j = cf_j + gauss_weight*(w+perturbation)*phi*((t*ele_length + r1+r2)/2)...
                        *N2*ele_length/2;
                end
                
                if (j~=front)
                    % next solve for coefficient of d alphaj+1/dt and the
                    % remaining part of d alphaj/dt
                
                      r1 = R(j); r2 = R(j+1);
                      cf_succesive = 0;
                     for t = -0.57735 : 2*0.57735 : 0.57735
                        gauss_weight = 1;
                        w = (r2 - ((t*ele_length + r1+r2)/2))/ele_length;
                        N2 = (((t*ele_length + r1+r2)/2)-r1)/ele_length;
                        N1 = w;
                        cf_succesive = cf_succesive + gauss_weight*(w-perturbation)*phi ...
                            *((t*ele_length + r1+r2)/2)*N2*ele_length/2;
                        cf_j = cf_j + gauss_weight*(w-perturbation)*phi*((t*ele_length + r1+r2)/2)...
                            *N1*ele_length/2;
                    end                      
                else

                    cf_succesive=0;
                    cf_j=cf_j+0;
                end
                
                % generate row of the coefficient matrix and the constant
                % vector corresponding to the node currently under
                % consideration. first order accurate backward difference
                % for time discretization is done to reduce computer storage.
                coeff_matrix(j,j-1)=cf_previous/dt-cf1;
                coeff_matrix(j,j)=cf_j/dt+cf1-cf2;
                if (j~=front)
                    coeff_matrix(j,j+1)=cf_succesive/dt+cf2;
                end
                if (j~=front)
                    const_vector(j)=(cf_previous*alpha(i,j-1)/dt)+(cf_j*alpha(i,j)/dt)+(cf_succesive*alpha(i,j+1)/dt);
                else
                    const_vector(j)=(cf_previous*alpha(i,j-1)/dt)+(cf_j*alpha(i,j)/dt); 
                end
            end
            ALPHA=coeff_matrix\const_vector; %final updated degree of cure at the
            % particular time step for the current row of nodes
            
            % place updated degree of cure vector into the existing alpha
            % matrix
            
            alpha(i,1:front)=ALPHA';
        end
    for i=nz+2:2*nz+1
        alpha(i,:)=alpha(2*nz+2-i,:);
    end     
    
    % Deduce the viscosity distribution using the current temperature and
    % degree of cure distributions obtained
    
    for i=1:2*nz+1
        for j=1:front
            b = 49.97 - 0.3278*TEMP(i,j) +  0.00055*(TEMP(i,j))^2;
            viscosity_dist(i,j)=Au*exp(Eu/(gasconst*TEMP(i,j)))*(alphag/(alphag-alpha(i,j)))^(a+b*alpha(i,j));
        end
    end
    
    % obtain the thickness averaged viscosity from the through thickness
    % viscosity obtained
    
    for i=1:front
        viscosity(i)=sum(viscosity_dist(:,(i)))/(2*nz+1);
    end
    
    if (front~=nr+1)
        viscosity(front+1)=viscosity(front); % In case the boundary condition
        %changes in the next flow time step
    end 
    
    isstatic = 0;
    current_force= force_computer(P , Vf , fillfactor, Area, front ,...
        R , isstatic);
    if current_force > wet_comp_force
        wet_comp_force = current_force;
    end
    clamping_force = [clamping_force current_force];
    injection_gate_pressure = [injection_gate_pressure P(1)];
    time = [time total_time];
end
total_time_saved = total_time;
% adjust the degree of cure at the outlet node. fluid particles might not
% reach the outlet 'exactly' due to minor violations of continuity.
% third relaxation phase of 0 second
isstatic = 1;
P = zeros (1,nr+1);
total_time = total_time + 0;
clamping_force = [clamping_force force_computer(P , Vf , fillfactor, Area, front ,...
        R , isstatic)];
injection_gate_pressure = [injection_gate_pressure 0];
time = [time total_time];
total_time = total_time + 2;
clamping_force = [clamping_force force_computer(P , Vf , fillfactor, Area, front ,...
        R , isstatic)];
injection_gate_pressure = [injection_gate_pressure 0];
time = [time total_time];
%% Display the outputs
% clc
if (max(max(alpha)) >= alphag)
    error('the resin could cure prematurely. change process parameters')
elseif (fillfactor(nr) < 0.75)
    error('refine number of elements used.')
else
% % % %     injection_time
% % % %     wet_compaction_time
% % %     figure (01)
% % % %     hold on
% % % %     plot (time, clamping_force);
% % %     title('clamping force vs time')
% % % %     figure (02)
% % % %     hold on
% % % %     plot (R,TEMP(nz+1,:))
% % % %     title('temperature vs position')
% % % %     figure(03)
% % % % %     hold on
% % % %     plot(1:length(spot_temp),spot_temp)
% % % %     figure(03)
% % % %     hold on
% % % %     plot(R,alpha(nz+1,:))
% % % %     title('degree of cure vs position')
% % %     
% % % % Check degree of cure at vent nodes before calling cure simulation.
% % % node_step = round((nr + 1)/4);
% % % cure_profile(:,1) = alpha(1:nz+1,1);
% % % cure_profile(:,2) = alpha(1:nz+1,1 + node_step);
% % % cure_profile(:,3) = alpha(1:nz+1,1 + 2*node_step);
% % % cure_profile(:,4) = alpha(1:nz+1,1 + 3*node_step);
% % % cure_profile(:,5) = alpha(1:nz+1,nr);
% % % cure_profile
% % % TEMP_profile(:,1) = TEMP(1:nz+1,1);
% % % TEMP_profile(:,2) = TEMP(1:nz+1,1 + node_step);
% % % TEMP_profile(:,3) = TEMP(1:nz+1,1 + 2*node_step);
% % % TEMP_profile(:,4) = TEMP(1:nz+1,1 + 3*node_step);
% % % TEMP_profile(:,5) = TEMP(1:nz+1,nr);
% % % TEMP_profile
% % %     if (alpha(:,nr+1)==zeros(2*nz+1,1))
% % %         travel_time=total_time_saved-time_stored;
% % %         for i=1:nz+1
% % %             alpha(i,nr+1)=runge_kutta(alpha_stored(i,1),TEMPf,travel_time);
% % %         end
% % %         for i = nz+2 : 2*nz + 1
% % %             alpha(i,nr+1) = alpha(2*nz + 2 - i, nr+1);
% % %         end
% % %     end
% % %     if(is_it_higher_level == 1) 
% % %         [cure_time, process_improvability] = cure(TEMP,alpha,heat_rate_1,dwell_temperature,dwell_time,heat_rate_2);
% % %         cure_time
% % %         process_improvability
% % %     end    
end
% % % toc
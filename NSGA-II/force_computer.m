function force = force_computer(P , Vf , fillfactor, Area, front , R , isstatic)
% A function that computes the force for given pressure distribution and
% glass fibre compaction. formula used obtained from walbran thesis.
% this function is compatible with circular moulds
a_static_dry=2.97e7;
b_static_dry=-3.6e7;
c_static_dry=1.72e7;
d_static_dry=-3.68e6;
e_static_dry=2.94e5;
a_static_wet=4.55e7;
b_static_wet=-5.80e7;
c_static_wet=2.84e7;
d_static_wet=-6.17e6;
e_static_wet=4.95e5;
a_dynamic_dry=4.68e7;
b_dynamic_dry=-5.14e7;
c_dynamic_dry=2.27e7;
d_dynamic_dry=-4.62e6;
e_dynamic_dry=3.58e5;
a_dynamic_wet=2.96e7;
b_dynamic_wet=-2.65e7;
c_dynamic_wet=9.50e6;
d_dynamic_wet=-1.54e6;
e_dynamic_wet=9.03e4;
forces = zeros ( 1 , length(Area));
ele_length = R(2) - R(1);
if (isstatic == 0)
    % first compute forces arising from fluid pressure in case of dynamic
    % mould
    for i = 1 : front-1
        r1 = R(i); r2 = R(i+1);
        for t = -0.57735 : 2*0.57735 : 0.57735
            gauss_weight = 1;
            forces(i) = forces(i) + gauss_weight*(2*pi/ele_length)*(P(i)*...
                ((t*ele_length + r1+r2)/2)*(r2 - ((t*ele_length + r1+r2)/2)) + ...
                P(i+1)*((t*ele_length + r1+r2)/2)*(((t*ele_length + r1+r2)/2) -r1))...
                *ele_length/2;
        end            
    end
    % add forces arising from fibre compaction
    for i = 1 : length(Area)
        fibre_stress_wet = (a_dynamic_wet*Vf^4 + b_dynamic_wet*Vf^3 +...
            c_dynamic_wet*Vf^2 + d_dynamic_wet*Vf + e_dynamic_wet);
        fibre_stress_dry = (a_dynamic_dry*Vf^4 + b_dynamic_dry*Vf^3 +...
            c_dynamic_dry*Vf^2 + d_dynamic_dry*Vf + e_dynamic_dry);
        forces(i) = forces(i) + fibre_stress_wet*Area(i)*fillfactor(i)+...
            fibre_stress_dry*Area(i)*(1-fillfactor(i));
    end
else
    for i = 1 : front - 1
        r1 = R(i); r2 = R(i+1);
        for t = -0.57735 : 2*0.57735 : 0.57735
            gauss_weight = 1;
            forces(i) = forces(i) + gauss_weight*(2*pi/ele_length)*(P(i)*...
                ((t*ele_length + r1+r2)/2)*(r2 - ((t*ele_length + r1+r2)/2)) + ...
                P(i+1)*((t*ele_length + r1+r2)/2)*(((t*ele_length + r1+r2)/2) -r1))...
                *ele_length/2;
        end             
    end
    % add forces arising from fibre compaction
    for i = 1 : length(Area)
        fibre_stress_wet = (a_static_wet*Vf^4 + b_static_wet*Vf^3 +...
            c_static_wet*Vf^2 + d_static_wet*Vf + e_static_wet);
        fibre_stress_dry = (a_static_dry*Vf^4 + b_static_dry*Vf^3 +...
            c_static_dry*Vf^2 + d_static_dry*Vf + e_static_dry);
        forces(i) = forces(i) + fibre_stress_wet*Area(i)*fillfactor(i)+...
            fibre_stress_dry*Area(i)*(1-fillfactor(i));
    end
end
force = sum(forces);

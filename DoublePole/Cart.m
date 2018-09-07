classdef Cart
  properties
    nPoles;  
    poles;      % objects of pole
    cart_pos;   % cart position
    cart_vel;   % cart velocity
    cart_acc;   % cart acceleration
    gravity;    % gravity, up is positive
    cart_mass;  % mass of the cart
    time;       % time in seconds
    applied_force;
    track_limit; % track limit from the center
    p_failure_angle; % pole failure angle +- 36 degrees from 0
    failed;
    time_step;
    cart_fric;  % cart friction
    p_fric;     % hinge friction
  end
  
  methods
    function object = initialize(object, h_len)  % double pole balancing
      object.nPoles = 2;
      object.cart_pos = 0.0;
      object.cart_vel = 0.0;
      object.cart_acc = 0.0;
      object.gravity = -9.8;
      object.cart_mass = 1.0;
      object.time = 0.0;
      object.applied_force = 0.0;
      object.track_limit = 2.4;
      object.p_failure_angle = degtorad(36);
      object.failed = false;
      object.time_step = 0.01;
      object.cart_fric = 0.0005;
      object.p_fric = 0.000002;
      object.poles = repmat(Pole, object.nPoles, 1);
      
      p_masses = [0.1;h_len/5];
      p_angles = [degtorad(1.0); degtorad(0.0)];
      p_h_lens = [0.5; h_len];
      p_accels = zeros(object.nPoles, 1);
      p_vels = zeros(object.nPoles, 1);
      for i = 1:object.nPoles
        object.poles(i) = initialize(object.poles(i), p_angles(i), ...
          p_vels(i), p_accels(i), p_masses(i), p_h_lens(i));
      end
      
%       p_masses = [0.1;0.03];
%       p_angles = [degtorad(1.0); degtorad(0.0)];
%       p_h_lens = [0.5; 0.15];
%       p_accels = zeros(object.nPoles, 1);
%       p_vels = zeros(object.nPoles, 1);
%       for i = 1:object.nPoles
%         object.poles(i) = initialize(object.poles(i), p_angles(i), ...
%           p_vels(i), p_accels(i), p_masses(i), p_h_lens(i));
%       end
    end
    
    function state = get_state(object)
      % scaled state
      state = zeros(6,1);
      state(1) = object.cart_pos / 2.4;
      state(2) = object.cart_vel / 10.0;
      state(3) = object.poles(1).angle / 0.628329;
      state(4) = object.poles(1).vel / 5.0;
      state(5) = object.poles(2).angle / 0.628329;
      state(6) = object.poles(2).vel / 13.0;
    end
      
    function cart_acc = cart_acc_given_vel(object, cart_vel)
      p_forces = 0.0;
      total_effictive_poles_mass = 0.0;
      for i = 1:object.nPoles
        p = object.poles(i);
        buff = object.p_fric*p.vel/(p.mass*p.h_len) + object.gravity*sin(p.angle);
        force = p.mass*p.h_len*p.vel^2*sin(p.angle) + ...
          3/4*p.mass*cos(p.angle)*buff;
        p_forces  = p_forces + force;
        total_effictive_poles_mass = total_effictive_poles_mass ...
          + p.mass*(1-3/4*cos(p.angle)^2);
      end
      cart_acc = object.applied_force - object.cart_fric*sign(cart_vel) + p_forces;
      cart_acc = cart_acc/(object.cart_mass + total_effictive_poles_mass);
    end
    
    function pole_acc = pole_acc_given_vel(object, pole, pole_vel)
      pole_acc = -3/(4*pole.h_len)*(object.cart_acc*cos(pole.angle) ...
        + object.gravity*sin(pole.angle) + object.p_fric*pole_vel/(pole.mass*pole.h_len));
    end
    
    function object = update_state(object)
      object.time = object.time + object.time_step;
      
      object.cart_vel = runge_kutta4(object, object.cart_vel, ...
        @cart_acc_given_vel, object.time_step);
      for i = 1:object.nPoles
        object.poles(i).vel = runge_kutta4_for_poles(object, object.poles(i), ...
          @pole_acc_given_vel, object.time_step);
      end
      
      object.cart_acc = cart_acc_given_vel(object, object.cart_vel);
      object.cart_pos = object.cart_pos + object.time_step*object.cart_vel;
      for i = 1:object.nPoles
        object.poles(i).angle = object.poles(i).angle ...
          + object.poles(i).vel*object.time_step;
        object.poles(i).acc = pole_acc_given_vel(object, object.poles(i), object.poles(i).vel);
      end
      
      if object.cart_pos > object.track_limit || object.cart_pos < - object.track_limit
        object.failed = true;
      end
      for i = 1:object.nPoles
        if object.poles(i).angle > object.p_failure_angle ...
            || object.poles(i).angle < -object.p_failure_angle
          object.failed = true;
          if object.poles(i).angle >= pi/2 
            object.poles(i).angle = pi/2;
            object.poles(i).vel = -0.9*object.poles(i).vel;
          elseif object.poles(i).angle <= -pi/2
            object.poles(i).angle = -pi/2;
            object.poles(i).vel = -0.9*object.poles(i).vel;
          end
        end
      end
    end
    
    function a = runge_kutta4(object, current_y, y_deriv_func, step)
      k1 = y_deriv_func(object, current_y);
      y1 = current_y + k1 * step / 2;
      k2 = y_deriv_func(object, y1);
      y2 = current_y + k2 * step / 2;
      k3 = y_deriv_func(object, y2);
      y3 = current_y + k3 * step;
      k4 = y_deriv_func(object, y3);
      
      a = current_y + (k1 + 2*k2 + 2*k3 + k4)*step/6;
    end
    
    function a = runge_kutta4_for_poles(object, pole, y_deriv_func, step)
      current_y = pole.vel;
      k1 = y_deriv_func(object, pole, current_y);
      y1 = current_y + k1 * step / 2;
      k2 = y_deriv_func(object, pole, y1);
      y2 = current_y + k2 * step / 2;
      k3 = y_deriv_func(object, pole, y2);
      y3 = current_y + k3 * step;
      k4 = y_deriv_func(object, pole, y3);
      
      a = current_y + (k1 + 2*k2 + 2*k3 + k4)*step/6;
    end
    
  end
end





















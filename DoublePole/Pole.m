classdef Pole
  properties
    angle;  
    vel;    % angle veolcity
    acc;    % acceleration
    mass;   % effective mass
    h_len;  % half length
  end
  
  methods
    function object = initialize(object, angle, velocity, acceleration, mass, half_length)
      object.angle = angle;
      object.vel = velocity;
      object.acc = acceleration;
      object.mass = mass;
      object.h_len = half_length;
    end
  end
end

function z=atan2mvr(y,x)
    z   = atan2(y,x);
    z += 2*pi*(z<0);;  % atan2 ->-pi/pi, now 0..2*pi
end

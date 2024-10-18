 function [var] = trapez(wzm,T,T_max)

if (T >= 0) && (T < T_max/5)
    var = wzm;
elseif (T >= T_max/5) && (T < T_max*4/5)
    var = 0.0001;
elseif (T >= T_max*4/5) && (T <= T_max)
    var = -wzm;
end


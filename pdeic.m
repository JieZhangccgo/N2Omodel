%% Code Initial Conditions
% The initial condition is applied at the first time value and provides the value of u(x,t0) for any value of x.
% The number of initial conditions must equal the number of equations, so for this problem there are 12 initial conditions.
% vector icMesh contains the sample points within the depth interval, and matrix icData contains the corresponding values, 
% each row indicating each component. x is the query points.
% interp1 is used to return interpolated values for multiple sets of data. Use matrix icData' as one input, columns are vectors for components.
% method 'pchip': shape-preserving piecewise cubic interpolation.

%% write a function that returns the initial condition
function u0 = pdeic(x,msInfo) 
icMesh = msInfo.icMesh;
icData = msInfo.icData;
u0 = interp1(icMesh,icData',x,'pchip')';

% a_down = 7;
% a_up = 10;
% if (x>icMesh(a_down-1)) && (x < icMesh(a_down))
%     u0(2)=icData(2,a_down-1); % update initial 2.DOC
%     u0(5)=icData(5,a_down-1); % update initial 5.NH4
%     u0(7)=icData(7,a_down-1);   % update initial 7.NO3
% elseif (x>icMesh(a_up)) && (x < icMesh(a_up+1))
%     u0(2)=icData(2,a_up+1);   % update initial 2.DOC
%     u0(5)=icData(5,a_up+1);   % update initial 5.NH4
%     u0(7)=icData(7,a_up+1);   % update initial 7.NO3
% end

end
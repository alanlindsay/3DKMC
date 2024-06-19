function vout = rotxVec(eta,v)
%   Rotate a list of vectors v around x by ETA (in radians)
%
%	eta  = column vector of angles
%   v    = 3-by-n input of n column vectors to be rotated.
%   vout = 3-by-n output of n column vectors that have been rotated 
% 
vx=v(1,:)';
vy=v(2,:)';
vz=v(3,:)';

vout = [vx';(vy.*cos(eta) -vz.*sin(eta))' ;(vy.*sin(eta) +vz.*cos(eta))']; 

end
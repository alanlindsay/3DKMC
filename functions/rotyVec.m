function vout = rotyVec(eta,v)
%   Rotate a list of vectors v around Y by ETA (in radians)
%
%	eta  = column vector of angles
%   v    = 3-by-n input of n column vectors to be rotated.
%   vout = 3-by-n output of n column vectors that have been rotated 
% 
vx=v(1,:)';
vy=v(2,:)';
vz=v(3,:)';

vout = [(vx.*cos(eta) +vz.*sin(eta))';vy' ;(-vx.*sin(eta)+vz.*cos(eta))']; 

end
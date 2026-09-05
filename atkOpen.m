function lhs = atkOpen(IP, PORT)

if nargin == 2 
    lhs = mexATKConnect('atkOpen', IP, PORT);
elseif nargin == 0
	  lhs = mexATKConnect('atkOpen');
else
	  error('incorrect number of inputs/outputs');
end
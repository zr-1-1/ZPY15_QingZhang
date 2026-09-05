function atkClose(rhs)

if nargin == 1
    mexATKConnect('atkClose', rhs);
elseif nargin == 0
	mexATKConnect('atkClose');
else
	  error('incorrect number of inputs/outputs');
end
function rtnData = atkConnect(conId, command, path, cmdParamString)

if (nargin == 4)
	rtnData = mexATKConnect('atkConnect', conId, command, [path ' ' cmdParamString]);
elseif (nargin == 3)
	rtnData = mexATKConnect('atkConnect', conId, command, path);
else
	error('incorrect number of inputs');
end
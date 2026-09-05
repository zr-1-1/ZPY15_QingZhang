function rtnData = atkExecuteCommand(conId, cmdType, command, path, cmdParamString)

if (nargin == 5)
	rtnData = mexATKConnect('atkExecuteCommand', conId, cmdType, command, [path ' ' cmdParamString]);
elseif (nargin == 4)
	rtnData = mexATKConnect('atkExecuteCommand', conId, cmdType, command, path);
else
	error('incorrect number of inputs');
end
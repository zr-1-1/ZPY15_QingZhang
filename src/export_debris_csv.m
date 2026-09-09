root = pwd;
doc = xmlread(fullfile(root, 'ZPY15.xml'));
satellites = doc.getDocumentElement().getElementsByTagName('Satellite');
fields = {'PositionX', 'PositionY', 'PositionZ', ...
          'VelocityX', 'VelocityY', 'VelocityZ'};
state = zeros(0, 6);
debris_id = zeros(0, 1);
epoch_text = cell(0, 1);
for k = 0:satellites.getLength()-1
    sat = satellites.item(k);
    token = regexp(char(sat.getAttribute('Name')), '^Debris(\d+)$', 'tokens', 'once');
    if isempty(token), continue; end
    orbit = sat.getElementsByTagName('Orbit').item(0);
    row = numel(debris_id) + 1;
    debris_id(row, 1) = str2double(token{1});
    epoch_text{row, 1} = char(orbit.getElementsByTagName('OrbEpoch').item(0).getTextContent());
    for j = 1:6
        node = orbit.getElementsByTagName(fields{j}).item(0);
        state(row, j) = str2double(char(node.getTextContent()));
    end
end
[debris_id, order] = sort(debris_id);
r0 = state(order, 1:3); % m
v0 = state(order, 4:6); % m/s
epoch = datetime(epoch_text(order), 'InputFormat', 'yyyy-MM-dd HH:mm:ss', ...
                 'TimeZone', 'UTC');
save(fullfile(root,'data','debris_initial_state.mat'), 'debris_id', 'r0', 'v0', 'epoch');
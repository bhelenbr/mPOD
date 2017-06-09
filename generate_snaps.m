function generate_snaps(nsnapshots,nelements)

% Random data
% for mode_counter=1:nsnapshots
%     mode = rand([nelements,1]);
%     save(['snapshot' num2str(mode_counter)],'mode');
% end

% Propagating perodic data over 1 period
x = (1/nelements:1/nelements:1)';
for mode_counter=1:nsnapshots
    data = exp(sin(2*pi*(x-mode_counter/nsnapshots)));
    %save(['snapshot' num2str(mode_counter)],'data');
    save(['snapshot' num2str(mode_counter) '.txt'],'data','-ascii');
end




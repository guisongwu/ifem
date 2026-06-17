figure(3);
clf;

plot_m3 = 3;
plot_n3 = 2;

%% u
% subplot(plot_m3, plot_n3, 1);
% trisurf(elem, node(:,1), node(:,2), uh(IUxNode), ...
%         'FaceColor', 'interp', 'EdgeColor', 'interp');
hold on;
patch('faces',elem, 'vertices', node, 'facevertexcdata', uh(IUxNode), ...
      'facecolor', 'interp', 'edgecolor', 'none');
quiver(node(:,1), node(:,2), uh(IUxNode), uh(IUyNode), 0.4, 'color', 'k', 'LineWidth', 1);
axis equal; axis tight; colorbar; view(2);
title('u', 'FontSize', 16);
hold off;
print('uofstokes','-depsc','-painters');


%% v
% subplot(plot_m3, plot_n3, plot_n3+1);
trisurf(elem, node(:,1), node(:,2), uh(IUyNode), ...
        'FaceColor', 'interp', 'EdgeColor', 'interp');
axis equal;
axis tight;
colorbar;
title('v', 'FontSize', 16);
view(2);
print('vofstokes','-depsc','-painters');


%% p
% subplot(plot_m3, plot_n3, 2*plot_n3+1);
trisurf(elem, node(:,1), node(:,2), ph, ...
        'FaceColor', 'interp', 'EdgeColor', 'interp');
axis equal;
axis tight;
colorbar;
title('p', 'FontSize', 16);
view(2);
print('pofstokes','-depsc','-painters');

%% bot
% subplot(plot_m3, plot_n3, 2);
% plot([sft(xbot); sft(xbot)+1], [sft(uh(IUxBot)); sft(uh(IUxBot))], '-');
% title('Bottom u', 'FontSize', 14);
% 
% subplot(plot_m3, plot_n3, 4);
% plot([sft(xbot); sft(xbot)+1], [sft(uh(IUyBot)); sft(uh(IUyBot))], '-');
% title('Bottom v', 'FontSize', 14);




figure(1);



gamma_stab = 1e-11;

stokes_info.nx = nx;
stokes_info.ny = ny;
stokes_info.xtop = xtop;
stokes_info.xbot = xbot;
stokes_info.IUxBot = IUxBot;
stokes_info.IUyBot = IUyBot;
stokes_info.IUxTop = IUxTop;
stokes_info.Mstab = Mstab;
stokes_info.gamma_stab = gamma_stab;

d2Xidm2_func = zeros(Nm);
option.use_newton = false;
for i = 1:Nm
    d2Xidm2_func(:,i) = stokes_hessian(node, elem, stokes_info, bdFlag, m, um, us1, extend_mid(EI(:,i)), option);
end

% dm_NW =  (d2Xidm2_NW + gamma_stab * Mstab) \ (dXidm(:) - gamma_stab * Mstab * m(1:n1));
% dm_GN =  (d2Xidm2_GN + gamma_stab * Mstab) \ (dXidm(:) - gamma_stab * Mstab * m(1:n1));
% dm_FD =  (d2Xidm2_FD + gamma_stab * Mstab) \ (dXidm_FD(:) - gamma_stab * Mstab * m(1:n1));

dm_NW =  (d2Xidm2_NW + gamma_stab * Mstab) \ (dXidm(:));
dm_GN =  (d2Xidm2_GN + gamma_stab * Mstab) \ (dXidm(:));
dm_FD =  (d2Xidm2_FD + gamma_stab * Mstab) \ (dXidm_FD(:));


% dm_NW =  (d2Xidm2_NW + gamma_stab * Mstab) \ (dXidm(:));
% dm_GN =  (d2Xidm2_GN + gamma_stab * Mstab) \ (dXidm(:));
% dm_FD =  (d2Xidm2_FD + gamma_stab * Mstab) \ (dXidm_FD(:));

% dm_NW =  (d2Xidm2_NW + gamma_stab * Mstab) \ (dXidm(:) + gamma_stab * Mstab * m(1:n1));
% dm_GN =  (d2Xidm2_GN + gamma_stab * Mstab) \ (dXidm(:) + gamma_stab * Mstab * m(1:n1));
% dm_FD =  (d2Xidm2_FD + gamma_stab * Mstab) \ (dXidm_FD(:) + gamma_stab * Mstab * m(1:n1));

% dm_func = d2Xidm2_func \ (dXidm(:) + gamma_stab * Mstab * m(1:n1));

dm_func = d2Xidm2_func \ (dXidm(:));

% dm_func =  (d2Xidm2_func + ) \ dXidm(:); 

% d2Xidm2
% dXidm(:)
% error('debug')

figure(2);
clf;

if 1
    % ALL
    plot(sft(xbot), sft_ext(dm_NW), 'o-', ...
         sft(xbot), sft_ext(dm_GN), 'o-', ...
         sft(xbot), sft_ext(dm_func), 'x-', ...
         sft(xbot), sft_ext(dm_FD), 's-', ...
         sft(xbot), sft(m - m0), '*--', ...
         'LineWidth', 2, 'MarkerSize', 10);
    legend('dm Newton', 'dm Gs-NW', 'dm func', 'dm FD', 'm - m0');
         % sft(xbot), sft(m - 2), '*--', ...
         % sft(xbot), sft(m0 - 2), 's-', ...
    %, 'm-mean', 'var'
end
if 0
    % ALL
    plot(sft(xbot), sft_ext(dm_NW) - dm_NW(1), 'o-', ...
         sft(xbot), (sft_ext(dm_GN) - dm_GN(1)) / 2, 'o-', ...
         sft(xbot), sft_ext(dm_func) - dm_func(1), 'x-', ...
         sft(xbot), sft_ext(dm_FD) - dm_FD(1), 's-', ...
         sft(xbot), sft(m - m0) - (m(1)-m0(1)), '*--', ...
         'LineWidth', 2, 'MarkerSize', 10);
    legend('dm Newton', 'dm Gs-NW', 'dm func', 'dm FD', 'm - m0');
end


if 0
    % Without FD
    plot(sft(xbot), sft_ext(dm_NW), 'o-', ...
         sft(xbot), sft_ext(dm_GN), 'o-', ...
         sft(xbot), sft_ext(dm_func), 'x-', ...
         sft(xbot), sft(m - m0), '*--', ...
         'LineWidth', 2, 'MarkerSize', 10);
    legend('dm Newton', 'dm Gs-NW', 'dm func', 'm - m0');
end



if 0
    plot(sft(xbot), sft_ext(dm), 'o-', ...
         sft(xbot), sft_ext(dm_func), 'x-', ...
         sft(xbot), sft_ext(dm_FD), 's-', ...
         sft(xbot), sft (m - m0), '*--', ...
         'LineWidth', 2, 'MarkerSize', 10);
    legend('dm adj', 'dm func', 'dm FD', 'm-m0');
elseif 0
    plot(xbot, dm, 'o-', ...
         xbot, dm_FD, 's-', ...
         xbot, m - m0, '*--', ...
         'LineWidth', 2, 'MarkerSize', 10);
    legend('dm adj', 'dm FD', 'm - m0');
elseif 0
    plot(sft(xbot), sft_ext(dm), 'o-', ...
         sft(xbot), sft(m - m0), '*--', ...
         sft(xbot), 0*xbot, '-', ...
         'LineWidth', 2, 'MarkerSize', 10);
    legend('dm adj', 'm - m0', 'zero');
end







if 0
    plot([sft(xbot); sft(xbot)+1], repmat(sft_ext(dm_NW),2,1), '-', ...
         [sft(xbot); sft(xbot)+1], repmat(sft_ext(dm_GN),2,1), '-', ...
         [sft(xbot); sft(xbot)+1], repmat(sft_ext(dm_FD),2,1), '-', ...
         [sft(xbot); sft(xbot)+1], repmat(sft(m - m0),2,1), 'o--', ...
         'LineWidth', 2, 'MarkerSize', 10);
    legend('dm Newton', 'dm Gs-Nw', 'dm FD', 'm-m0');
    %legend('dm Gs-Nw', 'dm FD', 'm-m0');
end

if 0
    plot(sft(xbot), sft_ext(dm_FD), 's-', ...
         sft(xbot), sft(m - m0), '*--', ...
         'LineWidth', 2, 'MarkerSize', 10);
    legend('dm FD', 'm-m0');
end


if 0
    plot(sft(xbot), sft(dm_NW), 'o-', ...
         sft(xbot), sft(dm_FD), 's-', ...
         'LineWidth', 2, 'MarkerSize', 10);
    legend('dm Newton', 'dm FD');
end


if 0
    plot(sft(xbot), sft(m - m0), '*--', ...
         sft(xbot), sft_ext(dm_GN), 'o-', ...
         sft(xbot), sft_ext(dm_func), 'o-', ...
         'LineWidth', 2, 'MarkerSize', 10);
    legend('m-m0', 'dm GN', 'dm Func');
end


if 0
    %dm_GN
    plot(sft(xbot), sft(m - m0), '*--', ...
         sft(xbot), sft_ext(dm_GN), 'o-', ...
         'LineWidth', 2, 'MarkerSize', 10);
    legend('m-m0', 'dm GN');
end

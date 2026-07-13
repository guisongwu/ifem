if 0
    % 
    % 
    % sum(Mstab) != 0, unstable
    % 
    % 
    Mstab = zeros(Nm);
    for ii = 1:Nm
        for jj = 1:Nm
            % Mstab(ii, jj) = integral_robin_P2(node, elem, bdFlag, ...
            %                                 extend_mid(EI(:,ii)), ...
            %                                 extend_mid(EI(:,jj)), [], option);
            Mstab(ii, jj) = integral_robin_penal(node, elem, bdFlag, ...
                                                 extend_mid(EI(:,ii)), ...
                                                 extend_mid(EI(:,jj)), option);
        end
    end
elseif 1
    % 
    % 
    % TODO: change to 1D gradient
    % 
    % 
    e = ones(Nm,1) / h;
    Mstab = spdiags([-e -e 2*e -e -e], [-(Nm-1) -1 0 1 (Nm-1)], Nm, Nm);
elseif 0
    Mstab = eye(Nm);
end

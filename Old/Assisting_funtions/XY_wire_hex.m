function [rc,zc]=XY_wire_hex(nv,nl,R,do,dr,dz)
    % calculation of design dimentions
%     nt = nv * nl;% total number of turns
    %% Position of each wire
    a = do/2;
    for i = 1:nl
        for j = 1:nv
    %             pos(i,j,1) = R + sin(60/360*2*pi)*(i-1)*do;       % Miguel
    %             pos(i,j,2) = (j-1)*(2*a + dz) + mod(i,2)*a;       % Miguel
                pos(i,j,1) = R+(2*i-1)*a-(i-1)*dr;
                pos(i,j,2) = (j-1/2)*(2*a+dz);
        end
    end

    % Remove the next for loop if using Miguel's code above
    for i = 2:2:N_l
        pos(i,:,2) = pos(i,:,2)+a+ksi_z/2;
    end

    posvect = [reshape( pos(:,:,1).' ,1,numel(pos(:,:,1))) ; reshape( pos(:,:,2).' ,1,numel(pos(:,:,2)))]';
    rc = posvect(:,1)';
    zc = posvect(:,2)';
end
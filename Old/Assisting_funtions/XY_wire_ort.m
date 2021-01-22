function [rc,zc]=XY_wire_ort(nv,nl,R,do,dr,dz) 
%% Position of each wire
a = do/2;
for i = 1:nl
    for j = 1:nv
%             pos(i,j,1) = R+(2*i-1)*a+(i-1)*dr;    % Miguel
%             pos(i,j,2) = (j-1)*(2*a + dz);        % Miguel
            pos(i,j,1) = R+(2*i-1)*a-(i-1)*dr;
            pos(i,j,2) = (j-1/2)*(2*a+dz);
    end
end

posvect = [reshape( pos(:,:,1).' ,1,numel(pos(:,:,1))) ; reshape( pos(:,:,2).' ,1,numel(pos(:,:,2)))]';
rc = posvect(:,1)';
zc = posvect(:,2)';
end
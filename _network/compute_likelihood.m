function gz_vals= compute_likelihood(model,z,X,x_off, y_off)

% compute likelihood vector g= [ log_g(z|x_1), ... , log_g(z|x_M) ] -
% this is for bearings and range case with additive Gaussian noise

M= size(X,2);
P= X([1 3],:);
   % Adjust for x and y offsets
    P(1, :) = P(1, :) + x_off; % Subtract x offset
    P(2, :) = P(2, :) + y_off; % Subtract y offset
%     % Adjust z for offsets
%     z_cartesian_x = z(2) * cos(z(1)) - x_off; % Convert z to Cartesian, adjust x
%     z_cartesian_y = z(2) * sin(z(1)) - y_off; % Convert z to Cartesian, adjust y
% 
%     % Convert adjusted z back to polar coordinates
%     z_adj = zeros(size(z));
%     z_adj(1) = atan2(z_cartesian_y, z_cartesian_x); % Adjusted bearing
%     z_adj(2) = sqrt(z_cartesian_x^2 + z_cartesian_y^2); % Adjusted range
% z=z_adj;
Phi= zeros(2,M);
Phi(1,:)= atan2(P(1,:),P(2,:));
Phi(2,:)= sqrt(sum(P.^2));
e_sq= sum( (diag(1./diag(model.D))*(repmat(z,[1 M])- Phi)).^2 );
gz_vals= exp(-e_sq/2 - log(2*pi*prod(diag(model.D))));

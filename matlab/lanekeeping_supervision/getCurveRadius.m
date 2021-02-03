function radius = getCurveRadius(track)
% Get signed radius of road curvature.
%
% Inputs:
%   track:  set of 2D-points defining the center of the track as a curve.
%
% Outputs:
%   radius: radius of road curvature at each point.

radius = zeros(length(track),1);
for i = 1:length(track)
    % Use three points to calculate the inscribed circle radius:
    if (i==1)
        A = track(:,end);
    else
        A = track(:,i-1);
    end
    B = track(:,i);
    if (i==length(track))
        C = track(:,1);
    else
        C = track(:,i+1);
    end
    
    angle = cross([A-B;0],[A-C;0]);
    
    if (abs(angle) < 1e-15)
        radius(i) = Inf;
    else
        AB = norm(A-B);
        BC = norm(B-C);
        CA = norm(C-A);
        
        R = AB*BC*CA/sqrt( (AB+BC+CA)*(BC+CA-AB)*(CA+AB-BC)*(AB+BC-CA) );
        s = sign(angle);
        radius(i) = s(3)*R;
    end
    
    %     % Use the signed curvature formula by calculating the first and second
    %     % order derivatives along the curve.
    %     if (i==1)
    %         x0 = track(1,end);
    %         y0 = track(2,end);
    %     else
    %         x0 = track(1,i-1);
    %         y0 = track(2,i-1);
    %     end
    %     x1 = track(1,i);
    %     y1 = track(2,i);
    %     if (i==length(track))
    %         x2 = track(1,1);
    %         y2 = track(2,1);
    %     else
    %         x2 = track(1,i+1);
    %         y2 = track(2,i+1);
    %     end
    %
    %     % Three-point formula to approximate first derivative along the curve.
    %     firstDer =  firstDerApprox([x0 y0],[x1 y1],[x2 y2]);
    %
    % 	% Use the above to approximate the second derivative.
    %     hneg =
    %     secondDer =
    
    
end
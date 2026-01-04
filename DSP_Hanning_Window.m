function [Y_win, H_win] = DSP_Hanning_Window(Y, H)
    % Hanning Window applied in DAFT domain
    % Coefficient: 0.5, -0.25, -0.25
    
    % Energy Scaling Factor
    scale = 1 / sqrt(0.375); 
    
    Y = Y(:);
    % Windowing Y: Convolution in Index
    % y[k] = 0.5*y[k] - 0.25*y[k-1] - 0.25*y[k+1]
    Y_win = (0.5*Y - 0.25*circshift(Y,1) - 0.25*circshift(Y,-1)) * scale;

    if nargin > 1
        % Windowing H: Apply to Columns (Impulse Responses)
        H_prev = circshift(H, 1, 1);
        H_next = circshift(H, -1, 1);
        H_win = (0.5*H - 0.25*H_prev - 0.25*H_next) * scale;
    else
        H_win = [];
    end
end
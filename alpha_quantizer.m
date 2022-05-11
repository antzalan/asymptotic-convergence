function quantized_x = alpha_quantizer(x, step, quant_channel, eta_lambda)
% alpha_quantizer: Function that quantizes a matrix in a per-element quantization
% approach. 
% 
% INPUTS:        x: Matrix to quantize
%             step: Quantization step for the quantization 
%    quant_channel: There are two possibilities. Either the quantization is
%                   applied to the input (channel matrix) or to the ouput (precoder) of the system.
%       eta_lambda: If quant_channel==True,  to avoid possible issues with singularities 
%                   and due to the fact that we use matrix inversions, we avoid 
%                   zero values by substituting them with a small enough random 
%                   complex number. The magnitude of this random number is given by 
%                   eta_lambda, which is considered so far as 1/10 of the quantization step.
%
% OUTPUT: Matrix of quantized values

    [nrow, ncol] = size(x);
    
    if nargin == 3
        eta_lambda = step/10;  
    else
        error('Wrong number of input parameter')
    end        
    
    if quant_channel % If quantization is done at input (channel matrix)

        quantized_x = step*floor(x/step + 1/2 + 1/2*1i);
        
        %%% Generating random eta values to avoid zero-NaN issues
        random_sign = floor(rand(2*nrow, 2*ncol) + 0.5);    
        random_sign(random_sign == 0) = -1; % Now it is either -1 or 1
        random_sign_complex = eta_lambda*(random_sign(1:2:2*nrow,(1:2:2*ncol)) + 1i*random_sign(2:2:2*nrow,2:2:2*ncol));

        quantized_x(quantized_x==0) = random_sign_complex(quantized_x==0); 

    else  % If quantization is done at output (precoding matrix)
        quantized_x = step*fix(x/step);        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Antonio Bazco-Nogueras  
% Date: 2022/05/11
% Contact: antonio.bazco@imdea.org
% License: This file can be distributed, remixed, adapted, and other work can be
% built upon it, as long as appropiate credit is included for the original creation. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
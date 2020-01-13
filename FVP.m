% Thanks to Terbe Dániel for this script (a new version formated in the standard format will follow shortly.)

function H = FVP(frames, Fs)
    %FVP Full Video Pulse extraction
    %   frames: sequence of RGB images (N, h, w, c)
    %   Fs: Frame rate
    
    N = size(frames, 1);
    K = 6; 
    L1 = Fs; 
    u0 = 1; 
    L2 = 128; 
    B = [6, 24];            % pulse band (calculated for L2!)
    Jt = zeros(4*K, N, 3);
    
    Pt = zeros(4 * K, N);
    Zt = zeros(4*K, N); 
    H = zeros(1, N);
    
    for i=1:N
        % Generate weighting masks
        Id = imresize(squeeze(frames(i,:,:,:)), [20, 20], 'box'); % downsize
        D = reshape(bsxfun(@rdivide, Id, sum(Id, 3)), [size(Id, 1)*size(Id, 2), 3]);
        
        A = pdist2(D, D, 'euclidean');
        
        [u, ~, ~] = svds(A, K);
        
        % correct the eigenvector signs
        u = bsxfun(@times, u, sign(sum(u0.*u, 1))); 
        u0 = u; 
        
        w = [u(:,1:K), -u(:, 1 : K)]'; 
        w = bsxfun(@minus, w, min(w, [], 2));
        w = bsxfun(@rdivide, w, sum(w, 2));
        
        % Weight and combine spatial pixel values
        J = bsxfun(@times, reshape(Id, [1, size(Id, 1) * size(Id, 2), 3]), w);
        Jt(:, i, :) = [mean(J, 2); var(J, [], 2)];
        
        % Extract rPPg-signals (with POS)
        if(i >= L1)
            C = Jt(:, i - L1 + 1 : i, :); 
            Cn = bsxfun(@rdivide, C, mean(C, 2)) - 1;
            X = Cn(:, :, 2) - Cn(:, :, 3); Y = Cn(:, :, 2) + Cn(:, :, 3) - 2 * Cn(:, :, 1);
            P = X + bsxfun(@times, std(X, [], 2)./std(Y, [], 2), Y); 
            Z = Cn(:, :, 1) + Cn(:, :, 2) + Cn(:, :, 3);
            
            Pt(:, i - L1 + 1 : i) = Pt(:, i - L1 + 1 : i) + bsxfun(@rdivide, bsxfun(@minus, P, mean(P, 2)), std(P, [], 2));
            Zt(:, i - L1 + 1 : i) = Zt(:, i - L1 + 1 : i) + bsxfun(@rdivide, bsxfun(@minus, Z, mean(Z, 2)), std(Z, [], 2));
        end
        
        % Combine rPPG signals into pulse
        if(i >= L2)
            P = Pt(:, i - L2 + 1 : i); 
            Z = Zt(:, i - L2 + 1 : i);
            
            Fp = fft(bsxfun(@rdivide, bsxfun(@minus, P, mean(P, 2)), std(P, [], 2)), [], 2)/L2;
            Fz = fft(bsxfun(@rdivide, bsxfun(@minus, Z, mean(Z, 2)), std(Z, [], 2)), [], 2)/L2;
            
            W = abs(Fp.* conj(Fp))./(1 + abs(Fz.* conj(Fz)));
            
            W(:, 1 : B(1) - 1) = 0;
            W(:, B(2) + 1 : end) = 0;
            h = real(ifft(sum(W.* Fp, 1), [], 2));

            H(:, i - L2 + 1 : i) = H(:, i - L2 + 1 : i) + (h - mean(h))/std(h);
        end
    end
end
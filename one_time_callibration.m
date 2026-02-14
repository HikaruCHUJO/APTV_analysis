clear; clc;

%% ==== ユーザ設定 ====
rootFolder          = '';   
coordCsv = '';%座標データ
subFolders          = 0:30;                                  
output_pre_folder   = '';    % 前処理画像の保存先
output_csv_folder   = '';

show_fit_overlay    = false;  
assert_trace_sigma  = false;  

window_half_size    = 16;


%% ==== 主処理 ====
Traw = readtable(coordCsv);

reqVars = {'id','x_abs','y_abs'};
missing = setdiff(reqVars, Traw.Properties.VariableNames);
Tcoord = table();

Tcoord.label = string(Traw.id);   

% x, y: ガウスフィットの中心を置く座標
Tcoord.x = Traw.x_abs;
Tcoord.y = Traw.y_abs;

% --- 解析パラメータ ---
a = 40; b = 40;                  % 楕円マスク半径（px）
window_half_size = 16;           % ROI 半径（px）
header = {'z','L1','L2','sigmaMajor','sigmaMinor','theta_deg_cont', ...
          'phi_used_deg','delta_deg','x','y','r'};

%% ==== 粒子ごとに 0~24 mm の全 z & 全フレームを一括解析 ====
for k = 1:height(Tcoord)

    % ---- 粒子ラベルと固定座標----
    label_k = Tcoord.label(k);   % 例: "1","2","3",...
    x0 = Tcoord.x(k);          
    y0 = Tcoord.y(k);            

    % 出力先: 
    out_label_folder = fullfile(output_csv_folder, char(label_k));
    if ~exist(out_label_folder, 'dir'); mkdir(out_label_folder); end

    fprintf('[LABEL %s] を解析開始（全 z）\n', label_k);

    for num = subFolders

        % ==== この z フォルダ内の画像一覧を取得 ====
        inputFolder = fullfile(rootFolder, num2str(num));

        imageFiles = dir(fullfile(inputFolder, '*.Bmp'));
        % ファイル名の「数字部分」でソート（0001.bmp, 0002.bmp, ...）
        fileNumbers = nan(numel(imageFiles),1);
        [~, ord]   = sort(fileNumbers);
        sortedFiles = imageFiles(ord);
        numImages = numel(sortedFiles);
        data_rows = [];

        % ==== 各フレーム処理====
        for idx = 1:numImages
            % --- 画像読み込み ---
            thisFile = sortedFiles(idx).name;
            [~, name, ext] = fileparts(thisFile);
            img0 = imread(fullfile(inputFolder, thisFile));
            img0 = im2gray(img0);         
            img  = double(img0);
            
            % 画像サイズ
            [rows, cols] = size(img);
            
            % =====================================================
            %  (1) 画像全体へのマスク処理
            % =====================================================
            [xGrid, yGrid] = meshgrid(1:cols, 1:rows);
            mask_outside = ((xGrid - x0).^2 / a^2) + ((yGrid - y0).^2 / b^2) >= 1;
    
            img_masked = img;             
            img_masked(mask_outside) = 0; 
    
            % ===========================================
            %  (2) マスク済み画像に前処理
            % ===========================================
            imgnormal0   = uint8(img_masked);
            imgnormal    = double(imgnormal0);
            %imgfilt_full = medfilt2(imgnormal, [1,1]);
            imgfilt_full = imgnormal;
            imgfilt_full = imgaussfilt(imgfilt_full,1);
    
            % ===========================================
            %  (3) マスク内の最大輝度位置を ROI の中心にする
            % ===========================================
            [max_val, max_idx] = max(img_masked(:));
            if max_val <= 0
                warning('label=%s, z=%d, frame=%s: マスク内に明るい画素がありません', ...
                        label_k, num, thisFile);
                continue;
            end
            [py_center, px_center] = ind2sub(size(img_masked), max_idx);
    
            % ===========================================
            %  (4) 前処理済み画像から ROI を切り出す（中心は最大輝度）
            % ===========================================
            xmin = max(1, round(px_center) - window_half_size);
            xmax = min(cols, round(px_center) + window_half_size);
            ymin = max(1, round(py_center) - window_half_size);
            ymax = min(rows, round(py_center) + window_half_size);
    
            if (xmax - xmin + 1) < 5 || (ymax - ymin + 1) < 5
                warning('ROI too small: label=%s, z=%d, frame=%s', label_k, num, thisFile);
                continue;
            end
    
            subimg = imgfilt_full(ymin:ymax, xmin:xmax);
            imgfilt = subimg;

            % ----- 回転なし 2D ガウス・フィット -----
            try
                [p, ~, exitflag] = fit2DGauss_subpixel_axisaligned(imgfilt);
                if exitflag <= 0
                    warning('lsqcurvefit 収束せず: label=%s, z=%d, frame=%s', ...
                            label_k, num, thisFile);
                    continue;
                end
            catch ME
                warning('fit 失敗: label=%s, z=%d, frame=%s (%s)', ...
                        label_k, num, thisFile, ME.message);
                continue;
            end

            % ----- ROI 内相対座標 → 画像全体の絶対座標 -----
            rel_cx = p(2);  rel_cy = p(3);     
            abs_cx = xmin + rel_cx - 1;        
            abs_cy = ymin + rel_cy - 1;        

            sigx = max(p(4), eps);  
            sigy = max(p(5), eps);  
            sigmaMajor     = max(sigx, sigy);
            sigmaMinor     = min(sigx, sigy);
            theta_deg_cont = 0;
            phi_used_deg   = 90;
            Delta_deg      = 0;
            L1 = sigy;   % 縦（Y）
            L2 = sigx;   % 横（X）

        
            if assert_trace_sigma
                Sigma = diag([sigx^2, sigy^2]);
                assert(abs((L1^2 + L2^2) - trace(Sigma)) < 1e-6*max(1,trace(Sigma)), ...
                       'L1^2+L2^2 と trace(Sigma) の不一致');
            end

            % 画像中心からの距離 r
            center_imgX = cols/2;
            center_imgY = rows/2;
            r_imgcenter = hypot(abs_cx - center_imgX, abs_cy - center_imgY);
            data_row = [ num, L1, L2, sigmaMajor, sigmaMinor, ...
                         theta_deg_cont, phi_used_deg, Delta_deg, ...
                         abs_cx, abs_cy, r_imgcenter ];
            data_rows = [data_rows; data_row];

            % ---- 可視化 ----
            if show_fit_overlay
                figure(1); clf;
                imshow(uint8(imgfilt_full), []); hold on;

                % FWHM 楕円（画像全体座標）
                fwhmX = 2*sqrt(2*log(2))*sigx;
                fwhmY = 2*sqrt(2*log(2))*sigy;
                t = linspace(0, 2*pi, 200);

                x_e = fwhmX/2*cos(t) + abs_cx; 
                y_e = fwhmY/2*sin(t) + abs_cy;   
                plot(x_e, y_e, 'r-', 'LineWidth', 1.5);
                plot(abs_cx, abs_cy, 'go', 'MarkerSize', 5, 'LineWidth', 1.5);

      
                rectangle('Position', [xmin, ymin, xmax-xmin+1, ymax-ymin+1],'EdgeColor', 'y', 'LineWidth', 1.5);

                title(sprintf('full image (masked+filtered) | z=%d | frame=%s', num, thisFile));
                hold off; drawnow;
            end

        end 

        % ==== abel に対する全フレーム分を CSV に出力 ====
        if ~isempty(data_rows)
            out_csv = fullfile(out_label_folder, sprintf('%d.csv', num));
            T = array2table(data_rows, 'VariableNames', header);
            writetable(T, out_csv);
        end

    end 
end 

% ========================================================================
% 関数
% ========================================================================
function [p, fitimg, exitflag] = fit2DGauss_subpixel_axisaligned(img)
    [Ny, Nx] = size(img);
    [X, Y] = meshgrid(1:Nx, 1:Ny);
    Z = double(img);
    z = Z(:);
    mass = sum(z);
    x0 = sum(X(:).*z) / max(mass, eps);
    y0 = sum(Y(:).*z) / max(mass, eps);
    B0 = prctile(z, 10);
    W  = max(Z - B0, 0);     
    M  = sum(W(:)) + eps;
    mx = sum(sum(W .* X)) / M;
    my = sum(sum(W .* Y)) / M;
    dx = X - mx; dy = Y - my;
    sxx = sum(sum(W .* dx .* dx)) / M;
    syy = sum(sum(W .* dy .* dy)) / M;
    sigX0 = sqrt(max(sxx, 0.5));
    sigY0 = sqrt(max(syy, 0.5));

    % 振幅初期値は「最大−背景」
    A0 = max(max(Z,[],'all') - B0, 1);
    p0 = [A0, x0, y0, sigX0, sigY0, B0];

    % --- モデル関数（回転なし）
    g = @(p,xy) p(1) * exp(-0.5 * ( ((xy(:,1)-p(2))./max(p(4),eps)).^2 + ...
                                     ((xy(:,2)-p(3))./max(p(5),eps)).^2 )) + p(6);

    % --- 下限/上限：中心は画像範囲、sigma は下限 >0、背景は自由
    lb = [0,   1,  1,   0.5, 0.5, -Inf];
    ub = [Inf, Nx, Ny,  Nx,  Ny,   Inf];

    % --- 最適化オプション
    opts = optimoptions('lsqcurvefit','Display','off', ...
        'MaxIter',300,'MaxFunctionEvaluations',5000, ...
        'FunctionTolerance',1e-8,'StepTolerance',1e-8);

    xy = [X(:), Y(:)];
    [p,~,~,exitflag] = lsqcurvefit(g, p0, xy, z, lb, ub, opts);

    % --- 生成画像（確認・デバッグ用）
    fitimg = reshape(g(p, xy), Ny, Nx);
end

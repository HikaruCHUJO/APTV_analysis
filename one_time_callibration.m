% ========================================================================
% APTV: 回転なし（軸平行）2Dガウス・フィット版
%  - 画像群（.Bmp）を読み、回転角θを推定しない軸平行ガウスでフィット
%  - L1/L2 は「縦枝（Y方向）= L1」「横枝（X方向）= L2」として定義
%    [z, L1, L2, sigmaMajor, sigmaMinor, theta_deg_cont, phi_used_deg, delta_deg, x, y, r]
% ========================================================================
clear; clc;

%% ==== ユーザ設定（パスや閾値など） ====込みます
rootFolder          = 'E:\chujo\00_experiment\callibration\2025_12_22\R\3_3';   % z掃引フォルダの
coordCsv = 'E:\chujo\01_code\13_on_time_calibration\2025_12_22\R\3_3\R3_3.csv';%座標データ
subFolders          = 0:30;                                   % 例: 0〜25
output_pre_folder   = 'E:\chujo\callibration\confirm\pre';    % 前処理画像の保存先
output_csv_folder   = 'E:\chujo\02_analysis\stock_data\12_22\R3_3\1_0';

% 表示/ログ
show_fit_overlay    = false;  % 軸平行楕円の重ね描き（必要なら true）
assert_trace_sigma  = false;  % L1^2+L2^2 ≒ trace(Sigma) を検証（数値ラウンド誤差に注意）

% 2Dガウスの切り出し窓（半径）—中心から ±window_half_size の正方ROIを切る
window_half_size    = 16;

% 出力先の作成（存在しない場合作成）
if ~exist(output_pre_folder, 'dir'); mkdir(output_pre_folder); end
if ~exist(output_csv_folder, 'dir'); mkdir(output_csv_folder); end

%% ==== 主処理（各サブフォルダ = 異なる z） ====
%% ==== 主処理：粒子ごとに 0→24 を一気に解析（非クリック・座標CSV使用） ====

%% ==== 座標CSV（Python で処理した left1_final.csv）を読み込む ====
if ~exist(coordCsv, 'file')
    error('粒子座標CSVが見つかりません: %s', coordCsv);
end

% ---- 生のCSVを読み込み（列: id, x_abs, y_abs, x_cv, y_cv, ...） ----
Traw = readtable(coordCsv);

% 必要な列が存在するかチェック（保険）
reqVars = {'id','x_abs','y_abs'};
missing = setdiff(reqVars, Traw.Properties.VariableNames);
if ~isempty(missing)
    error('座標CSVに必要な列がありません: %s', strjoin(missing, ', '));
end

Tcoord = table();

% label: フォルダ名に使うラベル
% ここでは id をそのまま "1","2","3",... のような文字列ラベルにする
Tcoord.label = string(Traw.id);   % → 出力フォルダ名: 1, 2, 3, ...

% x, y: ガウスフィットの中心を置く座標（x_abs, y_abs を使用）
Tcoord.x = Traw.x_abs;
Tcoord.y = Traw.y_abs;

% --- 解析パラメータ（従来どおり） ---
a = 40; b = 40;                  % 楕円マスク半径（px）
window_half_size = 16;           % ROI 半径（px）
header = {'z','L1','L2','sigmaMajor','sigmaMinor','theta_deg_cont', ...
          'phi_used_deg','delta_deg','x','y','r'};

%% ==== 粒子ごとに 0→24 の全 z & 全フレームを一括解析 ====
for k = 1:height(Tcoord)

    % ---- 粒子ラベルと固定座標（CSV 由来） ----
    label_k = Tcoord.label(k);   % 例: "1","2","3",...
    x0 = Tcoord.x(k);            % この粒子の x 座標（全 z で共通）
    y0 = Tcoord.y(k);            % この粒子の y 座標（全 z で共通）

    % 出力先: <output_csv_folder>/<label>/
    out_label_folder = fullfile(output_csv_folder, char(label_k));
    if ~exist(out_label_folder, 'dir'); mkdir(out_label_folder); end

    fprintf('[LABEL %s] を解析開始（全 z）\n', label_k);

    % ---- z を 0→24 まで走査 ----
    for num = subFolders

        % ==== この z フォルダ内の画像一覧を取得 ====
        inputFolder = fullfile(rootFolder, num2str(num));

        % Bmp / bmp の両方を検索
        imageFiles = dir(fullfile(inputFolder, '*.Bmp'));
        if isempty(imageFiles)
            imageFiles = dir(fullfile(inputFolder, '*.bmp'));
        end
        if isempty(imageFiles)
            warning('画像なし: %s', inputFolder);
            continue;
        end

        % ファイル名の「数字部分」でソート（0001.bmp, 0002.bmp, ...）
        fileNumbers = nan(numel(imageFiles),1);
        for ii = 1:numel(imageFiles)
            [~, nm] = fileparts(imageFiles(ii).name);
            fileNumbers(ii) = str2double(nm);  % 数字でない場合は NaN → そのままでも大きな問題はない
        end
        [~, ord]   = sort(fileNumbers);
        sortedFiles = imageFiles(ord);

        % ★ フレーム数
        numImages = numel(sortedFiles);

        % この z・粒子 label に対する全フレーム分の結果をためる配列
        data_rows = [];

        % ==== 各フレーム処理（クリック不要・座標CSVベース）====
        for idx = 1:numImages

            % --- 画像読み込み ---
            thisFile = sortedFiles(idx).name;
            [~, name, ext] = fileparts(thisFile);
            img0 = imread(fullfile(inputFolder, thisFile));
            img0 = im2gray(img0);          % カラーでもグレー化
            img  = double(img0);
            
            % 画像サイズ
            [rows, cols] = size(img);
            
            % =====================================================
            % =====================================================
            %  (1) 画像全体へのマスク処理（楕円の外側を 0 にする）
            % =====================================================
            [xGrid, yGrid] = meshgrid(1:cols, 1:rows);
            mask_outside = ((xGrid - x0).^2 / a^2) + ((yGrid - y0).^2 / b^2) >= 1;
    
            img_masked = img;              % 元画像をコピー
            img_masked(mask_outside) = 0;  % 楕円の外側を 0 にする
    
            % ===========================================
            %  (2) マスク済み画像に前処理（全体に適用）
            % ===========================================
            imgnormal0   = uint8(img_masked);
            imgnormal    = double(imgnormal0);
            %imgfilt_full = medfilt2(imgnormal, [1,1]);
            imgfilt_full = imgnormal;
            imgfilt_full = imgaussfilt(imgfilt_full,1);
    
            % ===========================================
            %  (3) マスク内の最大輝度位置を ROI の中心にする
            % ===========================================
            %   img_masked は「マスク外＝0」なので，
            %   max(img_masked(:)) は必ず楕円内（粒子があるはずの領域）になる
            [max_val, max_idx] = max(img_masked(:));
            if max_val <= 0
                warning('label=%s, z=%d, frame=%s: マスク内に明るい画素がありません', ...
                        label_k, num, thisFile);
                continue;
            end
            [py_center, px_center] = ind2sub(size(img_masked), max_idx);
            % ここで (px_center, py_center) が「実際の粒子中心（推定値）」
    
            % ===========================================
            %  (4) 前処理済み画像から ROI を切り出す（中心は最大輝度）
            % ===========================================
            xmin = max(1, round(px_center) - window_half_size);
            xmax = min(cols, round(px_center) + window_half_size);
            ymin = max(1, round(py_center) - window_half_size);
            ymax = min(rows, round(py_center) + window_half_size);
    
            if (xmax - xmin + 1) < 5 || (ymax - ymin + 1) < 5
                warning('ROI too small: label=%s, z=%d, frame=%s', ...
                        label_k, num, thisFile);
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
            rel_cx = p(2);  rel_cy = p(3);     % ROI 内の中心
            abs_cx = xmin + rel_cx - 1;        % 画像全体での中心X
            abs_cy = ymin + rel_cy - 1;        % 画像全体での中心Y

            sigx = max(p(4), eps);   % 横方向 σ = L2
            sigy = max(p(5), eps);   % 縦方向 σ = L1

            % ---- 派生量（従来と同じ定義） ----
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

            % ---- 1フレーム分の結果を 1 行として追加 ----
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

                x_e = fwhmX/2*cos(t) + abs_cx;   % ★ abs_cx
                y_e = fwhmY/2*sin(t) + abs_cy;   % ★ abs_cy
                plot(x_e, y_e, 'r-', 'LineWidth', 1.5);
                plot(abs_cx, abs_cy, 'go', 'MarkerSize', 5, 'LineWidth', 1.5);

      
                rectangle('Position', [xmin, ymin, xmax-xmin+1, ymax-ymin+1], ...
                          'EdgeColor', 'y', 'LineWidth', 1.5);

                title(sprintf('full image (masked+filtered) | z=%d | frame=%s', ...
                      num, thisFile));
                hold off; drawnow;
            end

        end % for idx (frames)

        % ==== abel に対する全フレーム分を CSV に出力 ====
        if ~isempty(data_rows)
            out_csv = fullfile(out_label_folder, sprintf('%d.csv', num));
            T = array2table(data_rows, 'VariableNames', header);
            writetable(T, out_csv);
        end

    end % for num (z)
end % for k (label)

% ========================================================================
% 関数群（回転なし 2D ガウス・フィット）
% ========================================================================
function [p, fitimg, exitflag] = fit2DGauss_subpixel_axisaligned(img)
% 回転なし（軸平行）2Dガウス
% p = [A, x0, y0, sigmaX, sigmaY, B]
%
% 画像座標 (X:列, Y:行) に沿って楕円等高線が軸平行になるように仮定。


    [Ny, Nx] = size(img);
    [X, Y] = meshgrid(1:Nx, 1:Ny);
    Z = double(img);
    z = Z(:);

    % --- 初期値（重心 & 重み付き二次モーメント）
    mass = sum(z);
    x0 = sum(X(:).*z) / max(mass, eps);
    y0 = sum(Y(:).*z) / max(mass, eps);

    % 背景の初期値は下位パーセンタイル（影響を受けにくい）
    B0 = prctile(z, 10);
    W  = max(Z - B0, 0);     % 背景を引いた非負重み
    M  = sum(W(:)) + eps;
    mx = sum(sum(W .* X)) / M;
    my = sum(sum(W .* Y)) / M;

    % 二次モーメント（軸平行を仮定するので共分散のオフ対角は使わない）
    dx = X - mx; dy = Y - my;
    sxx = sum(sum(W .* dx .* dx)) / M;
    syy = sum(sum(W .* dy .* dy)) / M;

    % sigma 初期値は二次モーメントの平方根（最小値を設定して病的ケースを回避）
    sigX0 = sqrt(max(sxx, 0.5));
    sigY0 = sqrt(max(syy, 0.5));

    % 振幅初期値は「最大−背景」
    A0 = max(max(Z,[],'all') - B0, 1);
    p0 = [A0, x0, y0, sigX0, sigY0, B0];

    % --- モデル関数（回転なし）
    % g(x,y) = A * exp( -0.5 * ( (x-x0)^2/sigX^2 + (y-y0)^2/sigY^2 ) ) + B
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

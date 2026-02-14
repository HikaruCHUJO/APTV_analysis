clear; close all; clc;

%% ---------------- ユーザ設定 ----------------
coord_map_csv     = '';
curves_master_csv = '';
imageFolder       = '';
save_name         = '';

monitor_full = false;     
monitor_stride = 1;       

% ==== 静的マスク（ベアリング・気泡） ====
use_static_mask    = true;
n_static_masks     = 4;
static_mask_mode   = 'manual';  
static_mask_centers = [];
static_mask_radii   = [];

% ==== ROI楕円マスク（円盤領域の外を落とす） ====
use_roi_mask = true;            
a_roi = 520;                    
b_roi = 440;

% ROI中心の決め方（=マスク中心）
mask_mode = 'fixed';            
xc_fixed  = 780;
yc_fixed  = 530;

% ==== L1/L2 定義 ====
useFWHM = false;

% ==== IDW（XY近傍） ====
idw_power = 2;
topK      = 10;

% ==== z探索 ====
z_grid_N  = 1201;
refine    = true;
refine_w  = 0.5;

% ==== しきい値 ====
use_threshold    = true;
dist_threshold   = 30;
threshold_metric = 'wrms';  

% ==== 短命バースト除去====
A_min        = 8;
Lmax         = 2;
eps_return   = 4;
burst_medwin = 7;

% ==== 極座標原点 ====
polar_origin_mode = 'image_center';
ox_manual = 645; oy_manual = 505;

% fps
fps = 167;

% ================================
% 自動追跡
% ================================
do_init_click_particle = true; 


tmpl_half   = 10;
search_half = 40;
fit_half    = 16;

min_corr_peak = 0.85;
gate_jump_px  = 55;

max_missing_consec      = 12;
search_expand_per_miss  = 0.25;

update_template       = false;
template_update_corr  = 0.75;

kf_dt      = 1;
kf_sigma_a = 2.0;
kf_sigma_m = 1.5;

% ================================
% モニタ（任意）
% ================================
show_monitor = true;
monitor_show_processed = true;  

%% ---------------- 1) CSV 読み込み ----------------
cache_mat = fullfile(fileparts(curves_master_csv), 'calib_cache_xy_curves.mat');
[labels, curves, xy_map] = loadOrBuildCalibCache(coord_map_csv, curves_master_csv, cache_mat);

%% ---------------- 2) 画像一覧 & 平均画像 ----------------
D = dir(imageFolder);
D = D(~[D.isdir]);
names = {D.name}';
isNumericBmp = ~cellfun('isempty', regexp(names, '^(?i)\d+\.bmp$', 'once'));
D = D(isNumericBmp);
nums = cellfun(@(s) str2double(regexp(s, '^\d+', 'match', 'once')), {D.name}');
[~, ord] = sort(nums, 'ascend');
imageFiles = D(ord);
assert(~isempty(imageFiles), '数値名の BMP が見つかりません: %s', imageFolder);
N = numel(imageFiles);

Sacc = 0;
for i = 1:N
    Iu = im2gray(imread(fullfile(imageFolder, imageFiles(i).name)));
    Sacc = Sacc + double(Iu);
end
mean_img = Sacc / N;

%% ---- 静的マスクを作成 ----
static_mask = false(size(mean_img));
if use_static_mask
    [Hmean, Wmean] = size(mean_img);
    [Xsm, Ysm] = meshgrid(1:Wmean, 1:Hmean);

    switch lower(static_mask_mode)
        case 'click_multi'
            for k = 1:n_static_masks
                fprintf('静的マスク %d 個目中心をクリック\n', k);
                [mx,my,stopMask] = pickPointRedWithFinish(mean_img);
                if stopMask || isnan(mx) || isnan(my)
                    warning('静的マスクキャンセル → 無効化');
                    use_static_mask = false;
                    static_mask = false(size(mean_img));
                    break;
                end
                R = static_mask_radii(min(k, numel(static_mask_radii)));
                static_mask = static_mask | ((Xsm-mx).^2+(Ysm-my).^2 <= R^2);
            end

        case 'manual'
            for k = 1:n_static_masks
                mx = static_mask_centers(k,1);
                my = static_mask_centers(k,2);
                R  = static_mask_radii(min(k, numel(static_mask_radii)));
                static_mask = static_mask | ((Xsm-mx).^2+(Ysm-my).^2 <= R^2);
            end

        case 'none'
            use_static_mask = false;
            static_mask = false(size(mean_img));

        otherwise
            error('static_mask_mode 不正');
    end
end

%% ---------------- ROI中心（xc,yc）決定 ----------------
xc = NaN; yc = NaN;
switch lower(mask_mode)
    case 'fixed'
        xc = xc_fixed; yc = yc_fixed;

    case 'once'
        Iclick = mean_img;
        if use_static_mask, Iclick(static_mask)=0; end
        [xc, yc, stopNow] = pickPointRedWithFinish(Iclick);
        if stopNow
            fprintf('ROI中心決定前に Finish → 終了\n');
            return;
        end

    case 'per_frame'
        xc = NaN; yc = NaN;

    otherwise
        error('mask_mode 不正');
end

%% ---------------- 出力テーブル ----------------
res = table('Size',[N 22], 'VariableTypes', ...
    {'double','double','double','double','double','double','double','double','double','double', 'string','double','double','double','double','double','double','double','double','double', 'double','double'}, ...
    'VariableNames', ...
    {'frame','X','Y','sigmaX','sigmaY','FWHMX','FWHMY', 'sigmaMaj','sigmaMin','theta_deg_raw','labels_used','wrms_L','K_used', 'z_euclid_idw','z_euclid_raw','dist_metric','theta_used_deg','aspect_meas','L1','L2', 'rho_px','time_s'});

trkXY = nan(N,2);
kF = 2*sqrt(2*log(2));

%% ---------------- tracking state ----------------
haveTrack = false;
missCount = 0;
tmpl = [];
kf = [];
last_valid_xy = [NaN; NaN];

%% ---------------- 観測----------------
mon = [];
if show_monitor
    mon = initMonitorFigure(mean_img, use_roi_mask, xc, yc, a_roi, b_roi, ...
                            use_static_mask, static_mask_centers, static_mask_radii);

    mon.hist_max = 300;


    mon.full   = monitor_full;
    mon.stride = max(1, round(monitor_stride));

    if ~monitor_full
        set(mon.pt_pred,'Visible','off');
        set(mon.pt_cand,'Visible','off');
        set(mon.pt_rej ,'Visible','off');
        set(mon.traj   ,'Visible','off');
        set(mon.rect_search,'Visible','off');
        set(mon.rect_fit   ,'Visible','off');
        set(mon.roi,'Visible','off');     
        mon.draw_windows = false;
    else
        mon.draw_windows = true;
    end
end
%% ---------------- 3) main loop ----------------
iProc = 0;

for i = 1:N

    Iraw_u8 = im2gray(imread(fullfile(imageFolder, imageFiles(i).name)));
    Iraw = double(Iraw_u8);
    [H,W] = size(Iraw);

    if strcmpi(mask_mode,'per_frame')
        Iclick = Iraw;
        if use_static_mask, Iclick(static_mask)=0; end
        [xc,yc,stopNow] = pickPointRedWithFinish(Iclick);
        if stopNow, break; end
    end
    IclickP = Iraw;
    if use_static_mask, IclickP(static_mask)=0; end

    Iproc = Iraw;
    if use_static_mask, Iproc(static_mask)=0; end

    Iproc = imgaussfilt(medfilt2(Iproc,[1,1]), 1);

    if use_roi_mask
        [Xg,Yg] = meshgrid(1:W, 1:H);
        roi_out = ((Xg-xc).^2/a_roi^2 + (Yg-yc).^2/b_roi^2) >= 1;
        Iproc(roi_out) = 0;
    end

    if ~haveTrack
        if do_init_click_particle
            [x0,y0,stopNow] = pickPointRedWithFinish(IclickP);
            if stopNow
                fprintf('初期粒子クリック前に Finish → 終了\n');
                break;
            end
            cx0=x0; cy0=y0;
        else
            [~,idxMax] = max(Iproc(:));
            [cy0,cx0] = ind2sub(size(Iproc), idxMax);
        end

        tmpl0 = getPatchFixed(Iproc, cx0, cy0, tmpl_half);
        if isempty(tmpl0) || std(double(tmpl0(:))) < 1e-6
            error('初期テンプレが作れません（端orフラット）。クリック位置/ROIサイズを見直してください。');
        end
        tmpl = double(tmpl0);

        kf = initKalmanCV2D([cx0;cy0], kf_dt, kf_sigma_a, kf_sigma_m);
        haveTrack = true;
        missCount = 0;
        last_valid_xy = [cx0;cy0];
    end

    kf = kf_predict(kf);
    x_pred = kf.x(1);
    y_pred = kf.x(2);

    expand = 1 + search_expand_per_miss * missCount;
    sh = round(search_half * expand);

    if isempty(tmpl) || std(tmpl(:)) < 1e-6
        tmplTry = getPatchFixed(Iproc, last_valid_xy(1), last_valid_xy(2), tmpl_half);
        if ~isempty(tmplTry) && std(double(tmplTry(:))) >= 1e-6
            tmpl = double(tmplTry);
        end
    end

    searchPatch = getPatchFixed(Iproc, x_pred, y_pred, sh);
    peakCorr = -Inf;
    x_cand = NaN; y_cand = NaN;

    if ~isempty(searchPatch) && ~isempty(tmpl) && all(size(searchPatch) > size(tmpl)) ...
            && std(tmpl(:)) >= 1e-6
        c = normxcorr2(tmpl, searchPatch);
        [peakCorr, idxPk] = max(c(:));
        [yPk, xPk] = ind2sub(size(c), idxPk);

        yOff = yPk - size(tmpl,1);
        xOff = xPk - size(tmpl,2);

        xMinS = round(x_pred) - sh;
        yMinS = round(y_pred) - sh;

        x_cand = xMinS + xOff + floor(size(tmpl,2)/2);
        y_cand = yMinS + yOff + floor(size(tmpl,1)/2);
    end

    isValid = isfinite(peakCorr) && (peakCorr >= min_corr_peak) ...
        && isfinite(x_pred) && isfinite(y_pred) ...
        && isfinite(x_cand) && isfinite(y_cand) ...
        && (hypot(x_cand-x_pred, y_cand-y_pred) <= gate_jump_px);

    cx=NaN; cy=NaN; sX=NaN; sY=NaN;

    if ~isValid
        missCount = missCount + 1;
        if missCount > max_missing_consec
            fprintf('missing連続 %d → 終了\n', missCount);
            break;
        end
    else
        subFit = getPatchFixed(Iproc, x_cand, y_cand, fit_half);
        if isempty(subFit)
            missCount = missCount + 1;
        else
            [pfit, ~, exitflag] = fit2DGauss_subpixel_axisaligned(subFit);
            if exitflag <= 0
                missCount = missCount + 1;
            else
                xMinF = round(x_cand) - fit_half;
                yMinF = round(y_cand) - fit_half;

                cx = xMinF + pfit(2) - 1;
                cy = yMinF + pfit(3) - 1;
                sX = max(pfit(4), eps);
                sY = max(pfit(5), eps);

                kf = kf_update(kf, [cx;cy]);
                last_valid_xy = [cx;cy];
                missCount = 0;

                if update_template && (peakCorr >= template_update_corr)
                    tmplNew = getPatchFixed(Iproc, cx, cy, tmpl_half);
                    if ~isempty(tmplNew) && std(double(tmplNew(:))) >= 1e-6
                        tmpl = double(tmplNew);
                    end
                end
            end
        end

        if missCount > max_missing_consec
            fprintf('fit/NCC missing連続 → 終了\n');
            break;
        end
    end

    if show_monitor && isvalid(mon.fig)
    
        imgShow = uint8(min(max(Iproc,0),255));
    
        updateMonitor(mon, imgShow, i, missCount, x_pred, y_pred, x_cand, y_cand, cx, cy, peakCorr, isValid, sh, fit_half);
    
        drawnow limitrate;
    
        if isappdata(mon.fig,'stopNow') && getappdata(mon.fig,'stopNow')
            fprintf('Finish pressed -> stop at frame %d\n', i);
            break;
        end
    end

    if isfinite(sX) && isfinite(sY)
        if useFWHM
            L1_meas = kF*sY; L2_meas = kF*sX;
        else
            L1_meas = sY;    L2_meas = sX;
        end
        asp_meas = L1_meas / max(L2_meas, eps);

        FWHMX = kF*sX; FWHMY = kF*sY;
        sigmaMaj = max(sX,sY);
        sigmaMin = min(sX,sY);

        z_each    = nan(numel(labels),1);
        dmin_each = nan(numel(labels),1);
        has_curve = false(numel(labels),1);

        for j = 1:numel(labels)
            lbch = char(labels(j));
            if ~isKey(curves, lbch), continue; end
            S = curves(lbch);

            zq  = linspace(S.zmin, S.zmax, z_grid_N);
            L1q = ppval(S.pp_L1, zq);
            L2q = ppval(S.pp_L2, zq);

            d2 = (L2_meas - L2q).^2 + (L1_meas - L1q).^2;
            [~,kmin] = min(d2);
            z_star = zq(kmin);

            if refine && z_grid_N >= 5
                dz = max(diff(zq));
                aL = max(S.zmin, z_star - max(refine_w, dz));
                aU = min(S.zmax, z_star + max(refine_w, dz));
                fun = @(z) (L2_meas-ppval(S.pp_L2,z)).^2 + (L1_meas-ppval(S.pp_L1,z)).^2;
                z_star = fminbnd(fun, aL, aU);
                d2_min = fun(z_star);
            else
                d2_min = d2(kmin);
            end

            z_each(j)    = z_star;
            dmin_each(j) = sqrt(d2_min);
            has_curve(j) = true;
        end

        use_idx = find(has_curve);
        z_raw = NaN; z_final = NaN; wrms = NaN; Kused = 0; lbl_used = string("");
        dmin_nearest = NaN;

        if ~isempty(use_idx)
            dXY = nan(numel(use_idx),1);
            for t = 1:numel(use_idx)
                lbch = char(labels(use_idx(t)));
                xy_lab = xy_map(lbch);
                dXY(t) = hypot(cx-xy_lab(1), cy-xy_lab(2));
            end
            dXY(dXY==0)=eps;

            if ~isempty(topK) && isfinite(topK) && topK < numel(use_idx)
                [~,ord2] = sort(dXY,'ascend');
                keep = ord2(1:topK);
            else
                keep = 1:numel(use_idx);
            end

            w = 1./(dXY(keep).^idw_power);
            w = w/sum(w);

            z_raw  = sum(w .* z_each(use_idx(keep)));
            wrms   = sqrt(sum(w .* (dmin_each(use_idx(keep)).^2)));
            dmin_nearest = min(dmin_each(use_idx(keep)));
            Kused  = numel(keep);
            lbl_used = strjoin(string(labels(use_idx(keep))), ',');
        end

        switch lower(threshold_metric)
            case 'wrms', dsel = wrms;
            case 'min',  dsel = dmin_nearest;
            otherwise,   dsel = wrms;
        end

        z_final = z_raw;
        if use_threshold && ~(isfinite(dsel) && dsel <= dist_threshold)
            z_final = NaN;
        end
    else
        [z_raw,z_final,wrms,Kused,lbl_used,dsel, ...
         L1_meas,L2_meas,asp_meas,FWHMX,FWHMY,sigmaMaj,sigmaMin] = ...
         deal(NaN,NaN,NaN,0,string(""),NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN);
    end
    iProc = iProc + 1;

    res.frame(iProc)        = i;
    res.X(iProc)            = cx;
    res.Y(iProc)            = cy;
    res.sigmaX(iProc)       = sX;
    res.sigmaY(iProc)       = sY;
    res.FWHMX(iProc)        = FWHMX;
    res.FWHMY(iProc)        = FWHMY;
    res.sigmaMaj(iProc)     = sigmaMaj;
    res.sigmaMin(iProc)     = sigmaMin;
    res.theta_deg_raw(iProc)= 0;
    res.theta_used_deg(iProc)= 90;

    res.aspect_meas(iProc)  = asp_meas;
    res.L1(iProc)           = L1_meas;
    res.L2(iProc)           = L2_meas;

    res.z_euclid_raw(iProc) = z_raw;
    res.z_euclid_idw(iProc) = z_final;
    res.labels_used(iProc)  = lbl_used;
    res.wrms_L(iProc)       = wrms;
    res.K_used(iProc)       = Kused;
    res.dist_metric(iProc)  = dsel;

    trkXY(iProc,:) = [cx,cy];

    if isfinite(fps), res.time_s(iProc)=(i-1)/fps; else, res.time_s(iProc)=NaN; end
end

res   = res(1:iProc,:);
trkXY = trkXY(1:iProc,:);

if iProc==0
    warning('解析フレーム 0');
    return;
end

%% ================= 4) スパイク形状外れ値除去 =================
z0 = res.z_euclid_idw(:);
z_clean = z0;
flag_burst_rm = false(numel(z0),1);

m = movmedian(z0, burst_medwin, 'omitnan');
dev = isfinite(z0) & isfinite(m) & (abs(z0-m) > A_min);
idx = find(dev);

if ~isempty(idx)
    brk = [true; diff(idx)>1; true];
    bgn = idx(brk(1:end-1));
    edd = idx(brk(2:end));

    for k = 1:numel(bgn)
        i0=bgn(k); i1=edd(k); L=i1-i0+1;
        if L>Lmax, continue; end

        pre = find(isfinite(z0(1:i0-1)) & isfinite(m(1:i0-1)), 1, 'last');
        if isempty(pre), continue; end
        m_pre = m(pre);

        post = find(isfinite(z0(i1+1:end)) & isfinite(m(i1+1:end)), 1, 'first');
        if isempty(post), continue; end
        post = post + i1;

        if abs(z0(post)-m_pre) <= eps_return
            z_clean(i0:i1)=NaN;
            flag_burst_rm(i0:i1)=true;
        end
    end
end

res.z_after_burst_clean = z_clean;
res.flag_burst_removed  = flag_burst_rm;

%% ================= 5) 極座標（θ,ρ） =================
[Hmean,Wmean] = size(mean_img);
switch lower(polar_origin_mode)
    case 'manual'
        ox=ox_manual; oy=oy_manual;
    case 'image_center'
        ox=Wmean/2; oy=Hmean/2;
    case 'click_once'
        [ox,oy] = pickPointRedWithFinish(mean_img);
    otherwise
        error('polar_origin_mode 不正');
end

dx = res.X-ox; dy = res.Y-oy;
theta_rad = atan2(-(dy), dx);
res.theta_used_deg = mod(rad2deg(theta_rad), 360);
res.rho_px = hypot(dx,dy);

%% ================= 6) CSV 保存 =================
traj_dir = fullfile(imageFolder, save_name);
if ~exist(traj_dir,'dir'), mkdir(traj_dir); end
traj_csv = fullfile(traj_dir, sprintf('traj_xyztheta_%s.csv', datestr(now,'yyyymmdd_HHMMSS')));

Ttraj = table( ...
    res.frame, res.time_s, res.X, res.Y, res.rho_px, res.theta_used_deg, ...
    res.z_euclid_raw, res.z_after_burst_clean, res.wrms_L, res.K_used, res.labels_used, res.L1, res.L2, 'VariableNames', {'frame','time_s','x_px','y_px','rho_px','theta_deg', 'z_raw_mm','z_clean_mm','wrms_L','K_used','labels_used','L1','L2'});

writetable(Ttraj, traj_csv);
fprintf('[traj] Saved: %s (rows=%d)\n', traj_csv, height(Ttraj));

%% ================= 7) 流跡線 =================
figure('Color','w');
imshow(repmat(mat2gray(mean_img),[1 1 3]), []); hold on;
plot(trkXY(:,1), trkXY(:,2), 'r-', 'LineWidth',1.6);
title('Trajectory (background: mean image)');

%% ================= 8) z vs theta =================
figure('Color','w'); hold on; grid on; box on;
ok = isfinite(res.theta_used_deg) & isfinite(res.z_after_burst_clean);
theta = res.theta_used_deg(ok); zplot = res.z_after_burst_clean(ok);
plot(theta, zplot, 'r.-', 'LineWidth',1.2, 'MarkerSize',10);
xlim([0 360]); xticks(0:45:360);
xlabel('\theta [deg]'); ylabel('z [mm]'); title('z vs \theta');

%% ================= 9) 3D =================
figure('Color','w'); hold on; grid on; box on;
ok3 = isfinite(res.X) & isfinite(res.Y) & isfinite(res.z_after_burst_clean);
X3=res.X(ok3); Y3=res.Y(ok3); Z3=res.z_after_burst_clean(ok3);
plot3(X3,Y3,Z3,'-','Color',[0.4 0.4 0.4],'LineWidth',1.2);
scatter3(X3,Y3,Z3,25,Z3,'filled');
xlabel('x [px]'); ylabel('y [px]'); zlabel('z [mm]');
title('3D trajectory (color=z)');
axis equal; view(45,25); colorbar;

%% ========================================================================
%  関数う
% ========================================================================

function mon = initMonitorFigure(img, use_roi_mask, xc, yc, a_roi, b_roi, use_static_mask, static_mask_centers, static_mask_radii)

    mon.fig = figure('Color','w','Name','Auto tracking monitor', 'NumberTitle','off','MenuBar','none','ToolBar','none');

    mon.ax  = axes('Parent',mon.fig);
    mon.im  = imshow(uint8(min(max(img,0),255)), 'Parent', mon.ax);
    hold(mon.ax,'on');   
    axis(mon.ax,'image'); axis(mon.ax,'on');

    mon.traj = plot(mon.ax, NaN, NaN, 'r-', 'LineWidth', 1.2);  

    mon.pt_meas = plot(mon.ax, NaN, NaN, 'ro', 'MarkerFaceColor','r', 'MarkerSize',3);
    mon.pt_pred = plot(mon.ax, NaN, NaN, 'c+', 'LineWidth',1.6, 'MarkerSize',10);     
    mon.pt_cand = plot(mon.ax, NaN, NaN, 'yo', 'LineWidth',1.2, 'MarkerSize',8);      
    mon.pt_rej  = plot(mon.ax, NaN, NaN, 'mx', 'LineWidth',1.6, 'MarkerSize',9);      

    mon.rect_search = rectangle(mon.ax,'Position',[1 1 1 1], ...
        'EdgeColor','y','LineWidth',1.2,'LineStyle','-','Visible','off');
    mon.rect_fit = rectangle(mon.ax,'Position',[1 1 1 1], ...
        'EdgeColor','g','LineWidth',1.2,'LineStyle','-','Visible','off');

    mon.roi = plot(mon.ax, NaN, NaN, 'w--', 'LineWidth', 1.2);
    if use_roi_mask && isfinite(xc) && isfinite(yc) && isfinite(a_roi) && isfinite(b_roi)
        tt = linspace(0,2*pi,360);
        xr = xc + a_roi*cos(tt);
        yr = yc + b_roi*sin(tt);
        set(mon.roi,'XData',xr,'YData',yr);
    end

    if use_static_mask && ~isempty(static_mask_centers)
        for k = 1:size(static_mask_centers,1)
            cxm = static_mask_centers(k,1);
            cym = static_mask_centers(k,2);
            R   = static_mask_radii(min(k,numel(static_mask_radii)));
            tt = linspace(0,2*pi,180);
            plot(mon.ax, cxm + R*cos(tt), cym + R*sin(tt), 'w-', 'LineWidth',1.0);
        end
    end

    mon.txt = text(mon.ax, 10, 20, '', 'Color','y','FontSize',12,'FontWeight','bold');

    title(mon.ax,'Finish: button or [F]','FontSize',12);
    setappdata(mon.fig,'stopNow',false);
    uicontrol('Style','pushbutton','String','Finish (F)',...
        'Units','normalized','Position',[0.80 0.93 0.18 0.06],...
        'FontSize',11,'Callback',@(~,~)setappdata(mon.fig,'stopNow',true));
    set(mon.fig,'KeyPressFcn',@(~,evt)onKey(evt));
    function onKey(evt)
        switch lower(evt.Key)
            case {'f','escape'}
                setappdata(mon.fig,'stopNow',true);
        end
    end

    mon.hist_max = 300;
    mon.draw_windows = true;
end


function updateMonitor(mon, img, frameIdx, missCount, x_pred, y_pred, x_cand, y_cand, cx, cy, corrPeak, isValid, sh, fit_half)
    mon.im.CData = img;
    
    if isfinite(cx) && isfinite(cy)
        mon.pt_meas.XData = cx;
        mon.pt_meas.YData = cy;
    else
        mon.pt_meas.XData = NaN;
        mon.pt_meas.YData = NaN;
    end

    if isfield(mon,'full') && ~mon.full
        mon.txt.String = sprintf('frame=%d  miss=%d', frameIdx, missCount);
        return;
    end

    if isempty(mon) || ~isvalid(mon.fig), return; end

    if isfield(mon,'full') && mon.full
        mon.im.CData = img;
    end

    if isfinite(cx) && isfinite(cy)
        mon.pt_meas.XData = cx;
        mon.pt_meas.YData = cy;
    else
        mon.pt_meas.XData = NaN;
        mon.pt_meas.YData = NaN;
    end

    if isfield(mon,'full') && ~mon.full
        mon.txt.String = sprintf('frame=%d  miss=%d', frameIdx, missCount);
        return;
    end

    if isfinite(x_pred) && isfinite(y_pred)
        mon.pt_pred.XData = x_pred;
        mon.pt_pred.YData = y_pred;
    else
        mon.pt_pred.XData = NaN; mon.pt_pred.YData = NaN;
    end

    if isfinite(x_cand) && isfinite(y_cand)
        mon.pt_cand.XData = x_cand;
        mon.pt_cand.YData = y_cand;
    else
        mon.pt_cand.XData = NaN; mon.pt_cand.YData = NaN;
    end

    if isValid && isfinite(cx) && isfinite(cy)
        mon.pt_meas.XData = cx;
        mon.pt_meas.YData = cy;

        xd = mon.traj.XData; yd = mon.traj.YData;
        xd = [xd, cx]; yd = [yd, cy];

        if numel(xd) > mon.hist_max
            xd = xd(end-mon.hist_max+1:end);
            yd = yd(end-mon.hist_max+1:end);
        end
        mon.traj.XData = xd;
        mon.traj.YData = yd;

        mon.pt_rej.XData = NaN; mon.pt_rej.YData = NaN;

    else
        mon.pt_meas.XData = NaN; mon.pt_meas.YData = NaN;
        if isfinite(x_cand) && isfinite(y_cand)
            mon.pt_rej.XData = x_cand;
            mon.pt_rej.YData = y_cand;
        else
            mon.pt_rej.XData = NaN; mon.pt_rej.YData = NaN;
        end
    end
    if mon.draw_windows
        if isfinite(x_pred) && isfinite(y_pred) && isfinite(sh)
            mon.rect_search.Position = [x_pred-sh, y_pred-sh, 2*sh, 2*sh];
            mon.rect_search.Visible = 'on';
        else
            mon.rect_search.Visible = 'off';
        end
    
        if isfinite(x_cand) && isfinite(y_cand) && isfinite(fit_half)
            mon.rect_fit.Position = [x_cand-fit_half, y_cand-fit_half, 2*fit_half, 2*fit_half];
            mon.rect_fit.Visible = 'on';
        else
            mon.rect_fit.Visible = 'off';
        end
    end

    if isfinite(corrPeak)
        mon.txt.String = sprintf('frame=%d  miss=%d  corr=%.3f  valid=%d', ...
            frameIdx, missCount, corrPeak, isValid);
    else
        mon.txt.String = sprintf('frame=%d  miss=%d', frameIdx, missCount);
    end
end

function patch = getPatchFixed(I,cx,cy,halfSize)
    [H,W]=size(I);
    if ~isfinite(cx)||~isfinite(cy), patch=[]; return; end
    cx=round(cx); cy=round(cy);
    xMin=cx-halfSize; xMax=cx+halfSize;
    yMin=cy-halfSize; yMax=cy+halfSize;
    if xMin<1||yMin<1||xMax>W||yMax>H, patch=[]; return; end
    patch = I(yMin:yMax, xMin:xMax);
end

function kf = initKalmanCV2D(xy0,dt,sigma_a,sigma_m)
    kf.dt=dt;
    kf.x=[xy0(1);xy0(2);0;0];
    kf.F=[1 0 dt 0; 0 1 0 dt; 0 0 1 0; 0 0 0 1];
    kf.H=[1 0 0 0; 0 1 0 0];
    kf.P=diag([25,25,100,100]);
    q=sigma_a^2;
    kf.Q=q*[dt^4/4 0 dt^3/2 0; 0 dt^4/4 0 dt^3/2; dt^3/2 0 dt^2 0; 0 dt^3/2 0 dt^2];
    r=sigma_m^2;
    kf.R=r*eye(2);
    kf.I=eye(4);
end

function kf = kf_predict(kf)
    kf.x = kf.F*kf.x;
    kf.P = kf.F*kf.P*kf.F' + kf.Q;
end

function kf = kf_update(kf,z)
    y = z - kf.H*kf.x;
    S = kf.H*kf.P*kf.H' + kf.R;
    K = kf.P*kf.H'/S;
    kf.x = kf.x + K*y;
    kf.P = (kf.I - K*kf.H)*kf.P;
end

function [p,fitimg,exitflag] = fit2DGauss_subpixel_axisaligned(img)
    [Ny,Nx]=size(img);
    [X,Y]=meshgrid(1:Nx,1:Ny);
    Z=double(img); z=Z(:);

    mass=sum(z);
    x0=sum(X(:).*z)/max(mass,eps);
    y0=sum(Y(:).*z)/max(mass,eps);

    B0=prctile(z,10);
    W=max(Z-B0,0);
    M=sum(W(:))+eps;
    mx=sum(sum(W.*X))/M;
    my=sum(sum(W.*Y))/M;

    dx=X-mx; dy=Y-my;
    sxx=sum(sum(W.*dx.*dx))/M;
    syy=sum(sum(W.*dy.*dy))/M;

    sigX0=sqrt(max(sxx,0.5));
    sigY0=sqrt(max(syy,0.5));
    A0=max(max(Z,[],'all')-B0,1);
    p0=[A0,x0,y0,sigX0,sigY0,B0];

    g=@(p,xy) p(1)*exp(-0.5*(((xy(:,1)-p(2))./max(p(4),eps)).^2 + ((xy(:,2)-p(3))./max(p(5),eps)).^2)) + p(6);

    lb=[0,1,1,0.5,0.5,-Inf];
    ub=[Inf,Nx,Ny,Nx,Ny,Inf];

    opts=optimoptions('lsqcurvefit','Display','off','MaxIter',300,'MaxFunctionEvaluations',5000, 'FunctionTolerance',1e-8,'StepTolerance',1e-8);

    xy=[X(:),Y(:)];
    [p,~,~,exitflag]=lsqcurvefit(g,p0,xy,z,lb,ub,opts);
    fitimg=reshape(g(p,xy),Ny,Nx);
end

function s2 = normalize_labels(s)
    s=string(s); s=strip(s);
    zenkaku=["０","１","２","３","４","５","６","７","８","９"];
    hankaku=["0","1","2","3","4","5","6","7","8","9"];
    for k=1:10, s=replace(s,zenkaku(k),hankaku(k)); end
    s=replace(s,"＿","_"); s=replace(s,"－","-");
    s2=lower(s);
end

function [x,y,stopNow] = pickPointRedWithFinish(img)
    stopNow=false;
    hFig=figure('Color','w','Name','クリックで選択 / Finishで終了', 'NumberTitle','off','MenuBar','none','ToolBar','none');
    hAx=axes('Parent',hFig);

    imgd = double(img);
    
    valid = isfinite(imgd) & (imgd > 0);
    if nnz(valid) < 100
        valid = isfinite(imgd);
    end
    
    lo = prctile(imgd(valid), 1);
    hi = prctile(imgd(valid), 99.9);   
    if ~isfinite(lo) || ~isfinite(hi) || hi <= lo
        lo = min(imgd(valid));
        hi = max(imgd(valid));
    end
    
    imshow(img, [lo hi], 'Parent', hAx);

    title(hAx,'クリックで確定 / 終了は [Finish] または [F]','FontSize',12);

    [H,W]=size(img);
    hHx=line(hAx,[1 W],[NaN NaN],'Color','r','LineWidth',1.5);
    hHy=line(hAx,[NaN NaN],[1 H],'Color','r','LineWidth',1.5);

    uicontrol('Style','pushbutton','String','Finish (F)',...
        'Units','normalized','Position',[0.76 0.93 0.22 0.06], 'FontSize',11,'Callback',@(~,~)finishNow());

    set(hFig,'WindowButtonMotionFcn',@onMove);
    set(hFig,'WindowButtonDownFcn',@onClick);
    set(hFig,'KeyPressFcn',@onKey);

    uiwait(hFig);
    pos=getappdata(hFig,'clickPos');
    if isempty(pos) && ~stopNow
        cp=get(hAx,'CurrentPoint');
        pos=cp(1,1:2);
    end
    if isvalid(hFig), close(hFig); end

    if stopNow, x=NaN; y=NaN;
    else, x=pos(1); y=pos(2);
    end

    function onMove(~,~)
        if ~isvalid(hAx), return; end
        cp=get(hAx,'CurrentPoint');
        xcp=cp(1,1); ycp=cp(1,2);
        if xcp>=1 && xcp<=W && ycp>=1 && ycp<=H
            set(hHx,'YData',[ycp ycp]);
            set(hHy,'XData',[xcp xcp]);
        end
    end
    function onClick(~,~)
        if ~isvalid(hAx), return; end
        cp=get(hAx,'CurrentPoint');
        setappdata(hFig,'clickPos',cp(1,1:2));
        uiresume(hFig);
    end
    function onKey(~,evt)
        switch lower(evt.Key)
            case {'f','escape'}
                finishNow();
            case {'return','enter'}
                onClick([],[]);
        end
    end
    function finishNow()
        stopNow=true;
        uiresume(hFig);
    end
end

function [labels,curves,xy_map] = loadOrBuildCalibCache(coord_map_csv, curves_master_csv, cache_mat)
    if exist(cache_mat,'file')
        dCache=dir(cache_mat); dM=dir(coord_map_csv); dC=dir(curves_master_csv);
        if dCache.datenum >= max(dM.datenum,dC.datenum)
            S=load(cache_mat,'labels','curves','xy_map');
            labels=S.labels; curves=S.curves; xy_map=S.xy_map;
            fprintf('[CACHE] Loaded: %s\n', cache_mat);
            return;
        end
    end

    optsM=detectImportOptions(coord_map_csv);
    optsM.SelectedVariableNames={'label','x','y'};
    optsM=setvartype(optsM,'label','string');
    M=readtable(coord_map_csv,optsM);

    optsC=detectImportOptions(curves_master_csv);
    optsC.SelectedVariableNames={'label','z','L1','L2'};
    optsC=setvartype(optsC,'label','string');
    C=readtable(curves_master_csv,optsC);

    M.label=normalize_labels(M.label);
    C.label=normalize_labels(C.label);

    labelsM=unique(M.label,'stable');
    labelsC=unique(C.label,'stable');
    labels=labelsM(ismember(labelsM,labelsC));

    xy_map=containers.Map('KeyType','char','ValueType','any');
    for i=1:height(M)
        xy_map(char(M.label(i))) = [M.x(i), M.y(i)];
    end

    curves=containers.Map('KeyType','char','ValueType','any');
    for i=1:numel(labels)
        lb=labels(i); lbch=char(lb);
        rows=(C.label==lb);
        if ~any(rows), continue; end
        zi=C.z(rows); L1i=C.L1(rows); L2i=C.L2(rows);
        [zuniq,ia]=unique(zi);
        L1u=L1i(ia); L2u=L2i(ia);
        curves(lbch)=struct('zmin',min(zuniq),'zmax',max(zuniq),...
            'pp_L1',pchip(zuniq,L1u),'pp_L2',pchip(zuniq,L2u));
    end

    save(cache_mat,'labels','curves','xy_map','-v7.3');
    fprintf('[CACHE] Saved: %s\n', cache_mat);
end

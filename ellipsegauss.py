import cv2
import numpy as np
from scipy.optimize import curve_fit
import csv
from pathlib import Path

# --- 2D Gaussian function ---
def twoD_Gaussian(coords, amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
    x, y = coords
    xo = float(xo); yo = float(yo)
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    g = offset + amplitude * np.exp(-(a*(x-xo)**2 + 2*b*(x-xo)*(y-yo) + c*(y-yo)**2))
    return g.ravel()

# --- ★ここで id の開始番号（オフセット）を設定する ---
# 例: 1枚目の校正板: 0, 
ID_OFFSET = 0  # ←ここを実験ごとに変える

# --- Load image ---
img_path = r"E:\chujo\01_code\13_on_time_calibration\12_13\l4_6\0-12142025032618-0.Bmp"
img = cv2.imread(img_path, cv2.IMREAD_GRAYSCALE)
if img is None:
    raise FileNotFoundError(f"Could not read the image: {img_path}")

in_path = Path(img_path)
out_stem = in_path.with_suffix('')

# --- Parameters ---
threshold_value = 97
min_contour_area = 5
max_contour_area = 100

# 長方形マスクの設定（解析する範囲）
use_rect_mask = True
mask_x0, mask_y0 = 530, 70
mask_x1, mask_y1 = 1300, 850

# --- Thresholding ---
_, thresh = cv2.threshold(img, threshold_value, 255, cv2.THRESH_BINARY)

# ★ 長方形マスクの適用 -----------------------------
if use_rect_mask:
    mask = np.zeros_like(thresh, dtype=np.uint8)
    cv2.rectangle(mask, (mask_x0, mask_y0), (mask_x1, mask_y1), 255, thickness=-1)
    thresh_masked = cv2.bitwise_and(thresh, mask)
else:
    thresh_masked = thresh

# --- Find contours ---
contours, _ = cv2.findContours(
    thresh_masked, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE
)
filtered_contours = [
    c for c in contours
    if min_contour_area <= cv2.contourArea(c) <= max_contour_area
]

# --- Prepare output image ---
output = cv2.cvtColor(img, cv2.COLOR_GRAY2BGR)
if use_rect_mask:
    cv2.rectangle(output, (mask_x0, mask_y0), (mask_x1, mask_y1), (255, 0, 255), 2)

rows = []

# --- Loop over each contour ---
for local_idx, cnt in enumerate(filtered_contours, start=1):
    if len(cnt) < 5:
        continue

    # ★ グローバル id = オフセット + ローカル番号
    gid = ID_OFFSET + local_idx

    # Fit ellipse
    ellipse = cv2.fitEllipse(cnt)
    (x_center, y_center), (x_len, y_len), angle_deg = ellipse

    # Draw ellipse and center
    cv2.ellipse(output, ellipse, (0, 0, 255), 2)
    cv2.circle(output, (int(x_center), int(y_center)), 3, (0, 255, 0), -1)
    cv2.putText(
        output, str(gid),  # ★ここも gid を表示
        (int(x_center - x_len//2), int(y_center - y_len//2)),
        cv2.FONT_HERSHEY_SIMPLEX, 0.6, (255, 0, 0), 2
    )

    # --- Crop ROI around the ellipse for Gaussian fit ---
    x0 = max(0, int(x_center - x_len))
    x1 = min(img.shape[1], int(x_center + x_len))
    y0 = max(0, int(y_center - y_len))
    y1 = min(img.shape[0], int(y_center + y_len))
    roi = img[y0:y1, x0:x1]
    if roi.size == 0:
        continue

    x = np.linspace(0, roi.shape[1]-1, roi.shape[1])
    y = np.linspace(0, roi.shape[0]-1, roi.shape[0])
    x, y = np.meshgrid(x, y)

    # ROI 内での中心座標に cv2.fitEllipse の中心を使う
    xo0 = x_center - x0
    yo0 = y_center - y0

    amp0    = float(roi.max() - roi.min())
    offset0 = float(roi.min())

    initial_guess = (
        amp0,
        xo0,
        yo0,
        max(x_len/4, 1.0),
        max(y_len/4, 1.0),
        0.0,
        offset0
    )


    try:
        #popt, _ = curve_fit(twoD_Gaussian, (x, y), roi.ravel(), p0=initial_guess)
        h, w = roi.shape

        lower_bounds = [0,     0,      0,    0.5, 0.5, -np.pi/2, 0]
        upper_bounds = [255*2, w,      h,    w,   h,   np.pi/2,  255*2]

        popt, _ = curve_fit(
            twoD_Gaussian,
            (x, y),
            roi.ravel(),
            p0=initial_guess,
            bounds=(lower_bounds, upper_bounds),
            maxfev=8000  # ← 試行回数を少し増やす
        )

    except RuntimeError:
        print(f"Could not fit Gaussian for ellipse gid={gid}")
        continue

    amplitude, xo, yo, sigma_x, sigma_y, theta, offset = popt

    x_abs = x0 + xo
    y_abs = y0 + yo

    print(
        f"Ellipse gid={gid}: cv_center=({x_center:.1f},{y_center:.1f}), "
        f"gauss_center=({x_abs:.2f},{y_abs:.2f}), "
        f"sigma_x={sigma_x:.2f}, sigma_y={sigma_y:.2f}, "
        f"theta(deg)={np.rad2deg(theta):.1f}"
    )

    rows.append({
        "id":       int(gid),       # ★ここも gid
        "x_abs":    float(x_abs),
        "y_abs":    float(y_abs),
        "x_cv":     float(x_center),
        "y_cv":     float(y_center),
        "sigma_x":  float(sigma_x),
        "sigma_y":  float(sigma_y),
        "theta_deg":float(np.rad2deg(theta)),
    })

# --- CSV出力 ---
out_csv = out_stem.parent / f"{out_stem.name}_particles.csv"
with open(out_csv, "w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(
        f,
        fieldnames=["id","x_abs","y_abs","x_cv","y_cv","sigma_x","sigma_y","theta_deg"]
    )
    writer.writeheader()
    writer.writerows(rows)
print(f"[info] Saved CSV: {out_csv} (n={len(rows)})")

# --- 画像保存 ---
out_img = out_stem.parent / f"{out_stem.name}_center_annotated.png"
ok = cv2.imwrite(str(out_img), output)
if ok:
    print(f"[info] Saved annotated image: {out_img}")
else:
    print("[warn] Failed to save annotated image.")

cv2.imshow("Gaussian Fitted Ellipses (masked)", output)
cv2.waitKey(0)
cv2.destroyAllWindows()

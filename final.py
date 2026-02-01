import csv
from pathlib import Path
import math  # 今回は使っていないが，将来的に距離計算などで使うかもしれないので残しておく

# ============================================================
#  粒子中心のズレが小さいデータだけを残して CSV を作り直すスクリプト
#
#  条件:
#    dx = |x_abs - x_cv| <= threshold_px
#    dy = |y_abs - y_cv| <= threshold_px
#
#  さらに改良:
#    ・条件を満たした行      → *_final.csv に保存
#    ・条件を満たさなかった行 → *_fail.csv に保存（除外データ一覧）
#
#  使い方:
#    1) in_csv_path を自分の CSV のパスに変える
#    2) 実行すると、同じフォルダに
#         - 元ファイル名 + "_final.csv"（採用された粒子）
#         - 元ファイル名 + "_fail.csv" （除外された粒子）
#       が出力される
# ============================================================

# --- (1) 入力CSVのパスを指定 --------------------------------
in_csv_path = r"E:\chujo\01_code\13_on_time_calibration\12_13\l1_1\0-12132025234334-0_particles.csv"
in_path = Path(in_csv_path)

# 出力CSVのパス（元ファイル名の末尾に "_final" / "_fail" を付ける）
out_csv_ok_path   = in_path.with_name(in_path.stem + "_final.csv")  # 条件を満たした行
out_csv_fail_path = in_path.with_name(in_path.stem + "_fail.csv")   # 条件を満たさなかった行

# --- (2) 閾値（ピクセル単位） --------------------------------
# ここを変更することで，どれくらいのズレまで許容するかを調整できる
threshold_px = 3.0  # 例: 「3ピクセル以内」を採用

# --- (3) CSVを読み込み（全行を rows_in に格納）--------------
rows_in = []
with open(in_csv_path, "r", encoding="utf-8") as f:
    reader = csv.DictReader(f)
    # ヘッダ（列名）はここで取得しておく
    fieldnames = reader.fieldnames
    if fieldnames is None:
        raise ValueError("入力CSVにヘッダ(列名)がありません。")

    # すべての行を一旦メモリに読み込む
    for row in reader:
        rows_in.append(row)

# --- (4) 出力用の行リストを作成 ------------------------------
rows_ok   = []  # 条件を満たした行（採用）
rows_fail = []  # 条件を満たさなかった行（除外）

# 追加で dx, dy も CSV に入れたい場合は True
add_dx_dy_columns = True

for row in rows_in:
    try:
        # 文字列 → float に変換
        x_abs = float(row["x_abs"])
        y_abs = float(row["y_abs"])
        x_cv  = float(row["x_cv"])
        y_cv  = float(row["y_cv"])
    except (KeyError, ValueError) as e:
        # 列がない / 数値に変換できない行はスキップ
        # 必要なら print をコメントアウト解除してデバッグに使う
        print(f"[warn] この行はスキップしました: {row} ({e})")
        continue

    # --- dx, dy を計算 ---
    dx = abs(x_abs - x_cv)
    dy = abs(y_abs - y_cv)

    # デバッグ・解析用に dx, dy を行データに追加
    if add_dx_dy_columns:
        row["dx"] = f"{dx:.3f}"
        row["dy"] = f"{dy:.3f}"

    # --- 条件で振り分ける ----------------------------------
    # 条件: dx <= threshold_px かつ dy <= threshold_px のものだけ「採用」
    if dx <= threshold_px and dy <= threshold_px:
        rows_ok.append(row)      # 採用データ
    else:
        rows_fail.append(row)    # 除外データ

# --- (5) 出力用のヘッダを準備 -------------------------------
# 元のヘッダをコピーして，必要なら dx, dy を追加
out_fieldnames = fieldnames.copy()
if add_dx_dy_columns:
    if "dx" not in out_fieldnames:
        out_fieldnames.append("dx")
    if "dy" not in out_fieldnames:
        out_fieldnames.append("dy")

# --- (6) 採用データ用 CSV (*_final.csv) を書き出し ----------
with open(out_csv_ok_path, "w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(f, fieldnames=out_fieldnames)
    writer.writeheader()
    writer.writerows(rows_ok)

# --- (7) 除外データ用 CSV (*_fail.csv) を書き出し -----------
with open(out_csv_fail_path, "w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(f, fieldnames=out_fieldnames)
    writer.writeheader()
    writer.writerows(rows_fail)

# --- (8) 結果を表示 -----------------------------------------
print(f"[info] 入力行数: {len(rows_in)}")
print(f"[info] 採用行数 (dx,dy <= {threshold_px}): {len(rows_ok)}")
print(f"[info] 除外行数 (dx,dy >  {threshold_px}): {len(rows_fail)}")
print(f"[info] 採用データ保存先: {out_csv_ok_path}")
print(f"[info] 除外データ保存先: {out_csv_fail_path}")

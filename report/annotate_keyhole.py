"""
annotate_keyhole.py
-------------------
Draws dimension annotations (depth and width) on keyhole.png
at the positions defined by the pink measurement box dots.
Saves the result as results/keyhole_annotated.png.
"""

from PIL import Image, ImageDraw, ImageFont

IMAGE_PATH  = "/home/mzomoro1/bin/lbf3/results/keyhole.png"
OUTPUT_PATH = "/home/mzomoro1/bin/lbf3/results/keyhole_annotated.png"

# Measured values
results = [
    {"label": "316L",        "depth": 340, "width": 57,
     # pink dot bounding box: x_min, x_max, y_min, y_max (from measure_keyhole.py)
     "x_min": 227.8, "x_max": 246.5, "y_min": 232.8, "y_max": 343.8},
    {"label": "316L + PLC",  "depth": 324, "width": 52,
     "x_min": 746.0, "x_max": 762.8, "y_min": 221.0, "y_max": 326.5},
]

img = Image.open(IMAGE_PATH).convert("RGB")
draw = ImageDraw.Draw(img)

# Try to load a small font; fall back to default
try:
    font = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 11)
    font_sm = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 10)
except:
    font = ImageFont.load_default()
    font_sm = font

WHITE  = (255, 255, 255)
YELLOW = (255, 220,  50)
TICK   = 4   # half-length of end ticks in px

for r in results:
    x0, x1 = r["x_min"], r["x_max"]
    y0, y1 = r["y_min"], r["y_max"]
    cx = (x0 + x1) / 2
    cy = (y0 + y1) / 2

    # ── depth annotation (vertical arrow on the left side of the box) ──
    ax = x0 - 8          # place it just left of the box
    # vertical line
    draw.line([(ax, y0), (ax, y1)], fill=YELLOW, width=1)
    # end ticks
    draw.line([(ax - TICK, y0), (ax + TICK, y0)], fill=YELLOW, width=1)
    draw.line([(ax - TICK, y1), (ax + TICK, y1)], fill=YELLOW, width=1)
    # label — rotated text not available in basic PIL, so place beside midpoint
    depth_text = f"{r['depth']} µm"
    bbox = draw.textbbox((0, 0), depth_text, font=font_sm)
    tw = bbox[2] - bbox[0]
    th = bbox[3] - bbox[1]
    draw.text((ax - tw - 3, cy - th / 2), depth_text, fill=YELLOW, font=font_sm)

    # ── width annotation (horizontal arrow below the box) ──
    ay = y1 + 8          # place it just below the box
    # horizontal line
    draw.line([(x0, ay), (x1, ay)], fill=YELLOW, width=1)
    # end ticks
    draw.line([(x0, ay - TICK), (x0, ay + TICK)], fill=YELLOW, width=1)
    draw.line([(x1, ay - TICK), (x1, ay + TICK)], fill=YELLOW, width=1)
    # label centred below the line
    width_text = f"{r['width']} µm"
    bbox = draw.textbbox((0, 0), width_text, font=font_sm)
    tw = bbox[2] - bbox[0]
    th = bbox[3] - bbox[1]
    draw.text((cx - tw / 2, ay + 3), width_text, fill=YELLOW, font=font_sm)

img.save(OUTPUT_PATH)
print(f"Saved: {OUTPUT_PATH}")

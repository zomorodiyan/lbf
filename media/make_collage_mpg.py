#!/usr/bin/env python3
"""
Compose 4-panel collage video:
  Left  column: testrun19_hires  top (with info overlay) stacked on bot
  Right column: testrun22_hires  top (with info overlay) stacked on bot
Outputs: collage_tr19_tr22.mp4  (25/3 fps ≈ 1 min for 501 frames)
"""

import os
import subprocess
import sys
from PIL import Image, ImageDraw, ImageFont

TR19 = "/home/mzomoro1/bin/lbf3/tutorials/laserbeamFoam/testrun19_hires/gif"
TR22 = "/home/mzomoro1/bin/lbf3/tutorials/laserbeamFoam/testrun22_hires/gif"
OUTDIR = "/home/mzomoro1/bin/lbf3/tutorials/laserbeamFoam"
TMPDIR = "/tmp/collage_frames"
N_FRAMES = 501

FONT_B = "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf"
FONT_R = "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"

# Info-box anchor (top-left corner of the empty grey area marked by user)
BOX_X, BOX_Y = 45, 175
LINE_H = 40

# Case metadata from transportProperties
CASES = {
    "tr19": {
        "title":    "316L",
        "subtitle": "(testrun19)",
        "sigma":    "σ [N/m] = 1.8648",
        "dsigma1":  "dσ/dT [N/(m·K)] =",
        "dsigma2":  "  −2.28×10⁻⁴",
    },
    "tr22": {
        "title":    "316L + 1% PLC",
        "subtitle": "(testrun21)",
        "sigma":    "σ [N/m] = 1.7675",
        "dsigma1":  "dσ/dT [N/(m·K)] =",
        "dsigma2":  "  +3.41×10⁻⁴",
    },
}


def load_fonts(size_b=27, size_r=24):
    try:
        fb = ImageFont.truetype(FONT_B, size_b)
        fr = ImageFont.truetype(FONT_R, size_r)
    except Exception:
        fb = fr = ImageFont.load_default()
    return fb, fr


def draw_info(img: Image.Image, meta: dict, fb, fr) -> Image.Image:
    """Draw info text directly onto img (no background box) and return RGB result."""
    img = img.convert("RGB")
    draw = ImageDraw.Draw(img)

    lines = [
        (meta["title"],    fb, "white"),
        (meta["subtitle"], fr, "white"),
        (meta["sigma"],    fr, "white"),
        (meta["dsigma1"],  fr, "white"),
        (meta["dsigma2"],  fr, "white"),
        (None, None, None),                        # spacer
        ("ARROW  velocity",  fr, "white"),         # white arrow + white label
        ("SQUARE  pressure", fr, "white"),         # solid black square + white label
    ]

    y = BOX_Y
    for text, font, color in lines:
        if text is None:
            y += LINE_H // 2
            continue

        if text.startswith("ARROW"):
            # White arrow, black label
            ax, ay = BOX_X + 2, y + LINE_H // 2
            draw.line([(ax, ay), (ax + 18, ay)], fill="white", width=3)
            draw.polygon([(ax+14, ay-5), (ax+24, ay), (ax+14, ay+5)], fill="white")
            label = text.replace("ARROW  ", "")
            draw.text((BOX_X + 32, y), label, fill="white", font=font)

        elif text.startswith("SQUARE"):
            # Solid black square + black outline, white label
            sq = 14
            sx, sy = BOX_X + 5, y + (LINE_H - sq) // 2
            draw.rectangle([sx, sy, sx + sq, sy + sq], fill="black", outline="black")
            label = text.replace("SQUARE  ", "")
            draw.text((BOX_X + 32, y), label, fill="white", font=font)

        else:
            draw.text((BOX_X, y), text, fill=color, font=font)

        y += LINE_H

    return img


def compose_frame(idx: int, fb, fr) -> Image.Image:
    s = f"{idx:04d}"
    top19 = Image.open(f"{TR19}/top_tr19h.{s}.png").convert("RGB")
    bot19 = Image.open(f"{TR19}/bot_tr19h.{s}.png").convert("RGB")
    top22 = Image.open(f"{TR22}/top_tr22h.{s}.png").convert("RGB")
    bot22 = Image.open(f"{TR22}/bot_tr22h.{s}.png").convert("RGB")

    top19 = draw_info(top19, CASES["tr19"], fb, fr)
    top22 = draw_info(top22, CASES["tr22"], fb, fr)

    w = top19.width          # 1000
    h_top = top19.height     # 911
    h_bot = bot19.height     # 400
    BORDER_TOP = 75          # black strip at top
    BORDER_BOT = 112         # black strip at bottom (1.5 × top)
    h_total = h_bot + h_top + BORDER_TOP + BORDER_BOT
    if h_total % 2:
        h_total += 1

    col = Image.new("RGB", (w * 2, h_total), color=(0, 0, 0))
    # bot_ on top, top_ on bottom
    col.paste(bot19, (0, BORDER_TOP))
    col.paste(top19, (0, BORDER_TOP + h_bot))
    col.paste(bot22, (w, BORDER_TOP))
    col.paste(top22, (w, BORDER_TOP + h_bot))

    # Vertical separator between the two cases
    sep_draw = ImageDraw.Draw(col)
    sep_draw.line([(w, 0), (w, h_total)], fill="white", width=3)

    return col


def main():
    preview_only = "--preview" in sys.argv
    os.makedirs(TMPDIR, exist_ok=True)

    fb, fr = load_fonts()

    if preview_only:
        frame = compose_frame(0, fb, fr)
        out = f"{OUTDIR}/collage_preview.png"
        frame.save(out)
        print(f"Preview saved: {out}")
        return

    print(f"Compositing {N_FRAMES} frames …")
    for i in range(N_FRAMES):
        compose_frame(i, fb, fr).save(f"{TMPDIR}/frame_{i:04d}.png")
        if i % 50 == 0:
            print(f"  {i}/{N_FRAMES-1}")

    out_mp4 = f"{OUTDIR}/collage_tr19_tr22.mp4"
    cmd = [
        "ffmpeg", "-y",
        "-framerate", "25/3",
        "-i", f"{TMPDIR}/frame_%04d.png",
        "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "18",
        "-vf", "crop=trunc(iw/2)*2:trunc(ih/2)*2",
        out_mp4,
    ]
    print("Encoding …")
    subprocess.run(cmd, check=True)
    print(f"Done: {out_mp4}")


if __name__ == "__main__":
    main()

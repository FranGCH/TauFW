import os
from PIL import Image

# Directory containing your PNGs
IMG_DIR = "./output_plots/"          # change to your path
OUT_DIR = "./combined_pre_post/" # choose your output folder
os.makedirs(OUT_DIR, exist_ok=True)

files = [f for f in os.listdir(IMG_DIR) if f.endswith(".png")]

# Create dictionaries for fast lookup
prefit = {f.replace("prefit_", ""): f for f in files if f.startswith("prefit_")}
postfit = {f.replace("postfit_", ""): f for f in files if f.startswith("postfit_")}

for key in sorted(prefit.keys()):
    if key in postfit:
        pre_path = os.path.join(IMG_DIR, prefit[key])
        post_path = os.path.join(IMG_DIR, postfit[key])

        pre_img = Image.open(pre_path)
        post_img = Image.open(post_path)

        # Match heights (resize prefit to postfit height)
        h = max(pre_img.height, post_img.height)
        pre_img = pre_img.resize((pre_img.width, h))
        post_img = post_img.resize((post_img.width, h))

        # Create combined image
        combined = Image.new("RGB", (pre_img.width + post_img.width, h))
        combined.paste(pre_img, (0, 0))
        combined.paste(post_img, (pre_img.width, 0))

        # Output name
        out_name = f"combined_{key}"
        combined.save(os.path.join(OUT_DIR, out_name))

        print(f"Created {out_name}")

print("Done!")

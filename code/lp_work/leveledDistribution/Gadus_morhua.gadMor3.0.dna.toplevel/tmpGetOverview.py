from pathlib import Path
import shutil

# Directory containing m=... directories
base_dir = Path(".")

# Destination directory
overview_dir = base_dir / "overview"
overview_dir.mkdir(exist_ok=True)

for m in range(10, 20):
    source_dir = base_dir / f"m={m},{m},{m}_t=32,64,128"
    source_file = source_dir / "cumulated_p_M_per_level.png"

    # Rename the copied file using the value of m
    destination_file = overview_dir / f"{m}.png"

    if source_file.exists():
        shutil.copy2(source_file, destination_file)
        print(f"Copied: {source_file} -> {destination_file}")
    else:
        print(f"File not found: {source_file}")
import os
import subprocess
import sys

input_file = sys.argv[1]
output_dir = sys.argv[2]

name = os.path.basename(input_file)

if name.startswith("adapter_content"):
    subprocess.run([
        "python3",
        "scripts/plot_adapter.py",
        input_file,
        output_dir,
    ])

elif name.startswith("per_cycle"):
    subprocess.run([
        "python3",
        "scripts/plot_quality.py",
        input_file,
        output_dir,
    ])

elif name.startswith("quality_distribution"):
    subprocess.run([
        "python3",
        "scripts/plot_distribution.py",
        input_file,
        output_dir,
    ])

else:
    print("Unknown TSV:", name)
    sys.exit(1)
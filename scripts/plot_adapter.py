#!/usr/bin/env python3

import csv
import os
import sys

import matplotlib.pyplot as plt


def read_tsv(filename):
    data = {}

    with open(filename, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")

        for field in reader.fieldnames:
            data[field] = []

        for row in reader:
            for field in reader.fieldnames:
                data[field].append(float(row[field]))

    return data


def main():

    if len(sys.argv) != 3:
        print("Usage:")
        print("python plot_adapter.py <adapter_content.tsv> <output_dir>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_dir = sys.argv[2]

    if not os.path.exists(input_file):
        print(f"File not found: {input_file}")
        sys.exit(1)

    data = read_tsv(input_file)

    if "pos" not in data:
        print("Invalid TSV: missing 'pos' column")
        sys.exit(1)

    x = data["pos"]

    colors = {
        "Universal": "#1f77b4",
        "SmallRNA3'": "#ff7f0e",
        "SmallRNA5'": "#2ca02c",
        "Nextera": "#d62728",
    }

    plt.figure(figsize=(10, 6))

    plotted = False

    for name in data:

        if name == "pos":
            continue

        y = data[name]

        #
        # Не рисуем полностью пустые линии
        #
        if max(y) == 0:
            continue

        plotted = True

        plt.plot(
            x,
            y,
            linewidth=2.0,
            color=colors.get(name, None),
            label=name,
        )

    #
    # Если адаптеров нет вообще
    #
    if not plotted:

        plt.text(
            max(x) / 2,
            50,
            "No adapters detected",
            fontsize=16,
            ha="center",
            va="center",
            color="gray",
        )

    plt.title(
        "Adapter Content per Position",
        fontsize=18,
        weight="bold",
    )

    plt.xlabel(
        "Position in read (bp)",
        fontsize=14,
    )

    plt.ylabel(
        "Adapter Content (%)",
        fontsize=14,
    )

    plt.xlim(1, max(x))
    plt.ylim(1, 100)

    plt.xticks(fontsize=11)
    plt.yticks(fontsize=11)

    plt.grid(
        linestyle="--",
        linewidth=0.7,
        alpha=0.5,
    )

    if plotted:
        plt.legend(
            loc="upper right",
            framealpha=0.95,
            fontsize=11,
        )

    plt.tight_layout()

    os.makedirs(output_dir, exist_ok=True)

    base = os.path.splitext(os.path.basename(input_file))[0]

    outfile = os.path.join(
        output_dir,
        f"adapters_plot_{base}.png",
    )

    plt.savefig(outfile, dpi=200)

    print(f"Saved: {outfile}")


if __name__ == "__main__":
    main()
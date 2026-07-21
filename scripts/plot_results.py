import matplotlib.pyplot as plt
import sys

if len(sys.argv) != 2:
    print("Usage: python plot_results.py <results_file>")
    sys.exit(1)

results_file = sys.argv[1]
print(f"Reading results from: {results_file}")

x = []
y = []

with open(results_file, 'r') as file:
    for line in file:
        row = line.split()
        if row:
            x.append(float(row[0]))
            y.append(float(row[1]))


plt.plot(x, y)

plt.ylim(0, 1)


plt.xlabel('Position in read (bp)')
plt.ylabel('% Adapter')
plt.title('Adapter Content')
plt.grid()

plt.savefig('../plots/adapters_plot_' + results_file[25:-4] + '.png')
print(f"Plot saved as: ../plots/adapters_plot_{results_file[25:-4]}.png")
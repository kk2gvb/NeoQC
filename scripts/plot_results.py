import matplotlib.pyplot as plt

x = []
y = []

with open('../results/adapter_stats.txt', 'r') as file:
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

plt.savefig('../plots/adapters_plot.png')

#TODO: Надо сделать так, чтобы файл с результатами анализатора передавался в скрипт как аргумент, а не был захардкожен.
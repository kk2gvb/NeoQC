import sys
import os
import csv
import matplotlib.pyplot as plt

def main():
    # 1. Проверка аргументов (соответствует вызову из PlotRunner.cpp)
    if len(sys.argv) != 3:
        print("Usage: python plot_results.py <input_tsv> <output_dir>")
        sys.exit(1)

    input_tsv = sys.argv[1]
    output_dir = sys.argv[2]

    if not os.path.exists(input_tsv):
        print(f"Error: Input file not found: {input_tsv}")
        sys.exit(1)

    # 2. Чтение данных (динамически, через csv модуль)
    data = {}
    headers = []
    
    try:
        with open(input_tsv, 'r', encoding='utf-8') as f:
            reader = csv.reader(f, delimiter='\t')
            headers = next(reader)  # Первая строка — заголовки
            
            # Инициализируем списки для каждого столбца
            for h in headers:
                data[h] = []
            
            # Читаем строки
            for row in reader:
                for i, val in enumerate(row):
                    if i < len(headers):
                        try:
                            data[headers[i]].append(float(val))
                        except ValueError:
                            data[headers[i]].append(0.0)
    except Exception as e:
        print(f"Error reading TSV file: {e}")
        sys.exit(1)

    if not headers or 'pos' not in headers:
        print("Error: Invalid TSV format. 'pos' column not found.")
        sys.exit(1)

    # 3. Построение графика
    plt.figure(figsize=(10, 6))
    
    x = data['pos']
    
    # Рисуем линии для всех адаптеров, кроме столбца 'pos'
    for header in headers:
        if header != 'pos':
            plt.plot(x, data[header], label=header, linewidth=1.5, alpha=0.8)

    # 4. Настройка осей и оформления
    plt.xlabel('Position in read (bp)', fontsize=12)
    plt.ylabel('Adapter Content (%)', fontsize=12)  # Исправлено: теперь это проценты
    plt.title('Adapter Content per Position', fontsize=14)
    
    plt.ylim(-0.01, 1)
    plt.xlim(0, max(x) if x else 10)
    
    plt.legend(loc='upper right', fontsize='small', framealpha=0.9)
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.tight_layout()

    # 5. Сохранение (безопасное формирование пути)
    os.makedirs(output_dir, exist_ok=True)
    
    # Берём только имя файла без расширения, например: "adapter_content_R1"
    base_name = os.path.basename(input_tsv).replace('.tsv', '')
    output_path = os.path.join(output_dir, f"adapters_plot_{base_name}.png")
    
    plt.savefig(output_path, dpi=150)
    print(f"Success: Plot saved to {output_path}")

if __name__ == "__main__":
    main()
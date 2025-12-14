import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import make_interp_spline

# Чтение данных из файлов
data1 = np.loadtxt('signal_data.txt')
x = data1[:, 0]
y_original = data1[:, 1]
y_filtered = data1[:, 2]

# Интерполяция сплайнами
x_smooth = np.linspace(x.min(), x.max(), 300)
spline_original = make_interp_spline(x, y_original, k=3)
spline_filtered = make_interp_spline(x, y_filtered, k=3)
y_original_smooth = spline_original(x_smooth)
y_filtered_smooth = spline_filtered(x_smooth)

# Построение графика
plt.figure(figsize=(12, 6))
plt.plot(x_smooth, y_original_smooth, 'b-', linewidth=2, alpha=0.7, label='Исходный сигнал')
plt.plot(x_smooth, y_filtered_smooth, 'r-', linewidth=2, label='После фильтрации')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Фильтрация зашумленного сигнала')
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('filtering_result.png', dpi=150)
plt.show()
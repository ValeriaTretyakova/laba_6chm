#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <chrono>
#include <algorithm>
#include <fstream>
#include <direct.h>
#include <iomanip> 

using namespace std;
using Complex = complex<double>;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#define PI M_PI
#endif

// Реализация DFT и IDFT
vector<Complex> DFT(const vector<Complex>& input) {
    int N = input.size();
    vector<Complex> output(N);
    for (int k = 0; k < N; ++k) {
        Complex sum(0, 0);
        for (int n = 0; n < N; ++n) {
            double angle = -2 * PI * k * n / N;
            sum += input[n] * Complex(cos(angle), sin(angle));
        }
        output[k] = sum;
    }
    return output;
}

vector<Complex> IDFT(const vector<Complex>& input) {
    int N = input.size();
    vector<Complex> output(N);
    for (int n = 0; n < N; ++n) {
        Complex sum(0, 0);
        for (int k = 0; k < N; ++k) {
            double angle = 2 * PI * k * n / N;
            sum += input[k] * Complex(cos(angle), sin(angle));
        }
        output[n] = sum / static_cast<double>(N);
    }
    return output;
}

// Реализация FFT и IFFT (in-place)
void FFT(vector<Complex>& a) {
    int N = a.size();
    if (N <= 1) return;
    vector<Complex> even(N / 2), odd(N / 2);
    for (int i = 0; i < N / 2; ++i) {
        even[i] = a[i * 2];
        odd[i] = a[i * 2 + 1];
    }
    FFT(even);
    FFT(odd);
    for (int k = 0; k < N / 2; ++k) {
        Complex t = polar(1.0, -2 * PI * k / N) * odd[k];
        a[k] = even[k] + t;
        a[k + N / 2] = even[k] - t;
    }
}

void IFFT(vector<Complex>& a) {
    int N = a.size();
    for (auto& x : a) x = conj(x);
    FFT(a);
    for (auto& x : a) x = conj(x) / static_cast<double>(N);
}

// Генерация зашумленного сигнала
vector<Complex> generateSignal(int N, double A, double omega1, double phi, double B, double omega2) {
    vector<Complex> signal(N);
    for (int j = 0; j < N; ++j) {
        double value = A * cos(2.0 * PI * omega1 * j / N + phi) + B * cos(2.0 * PI * omega2 * j / N);
        signal[j] = Complex(value, 0);
    }
    return signal;
}

// Генерация сигнала с разрывами
vector<Complex> generateSignal2(int N, double A, double B, double omega2) {
    vector<Complex> z(N);
    int N4 = N / 4;
    int N2 = N / 2;
    int N34 = 3 * N / 4;

    for (int j = 0; j < N; ++j) {
        double value = 0.0;
        if (j >= N4 && j <= N2) {
            value = A + B * cos(2.0 * PI * omega2 * j / N);
        }
        else if (j > N34) {
            value = A + B * cos(2.0 * PI * omega2 * j / N);
        }
        z[j] = Complex(value, 0);
    }
    return z;
}

// Вывод таблицы спектров
void printSpectrumTable(const vector<Complex>& z_original, const vector<Complex>& Z, const string& title) {
    int N = static_cast<int>(Z.size());
    cout << "\n" << title << ":\n";
    cout << setw(3) << "m" << " | "
        << setw(12) << "Re(z)" << " | "
        << setw(12) << "Re(Z)^" << " | "
        << setw(12) << "Im(Z)^" << " | "
        << setw(15) << "Амплитуда" << " | "
        << setw(12) << "Фаза" << endl;

    cout << scientific << setprecision(6);

    for (int m = 0; m < N; ++m) {
        double amp = abs(Z[m]);
        if (amp > 1e-10) {
            double phase = arg(Z[m]);

            cout << setw(3) << m << " | "
                << setw(12) << z_original[m].real() << " | "
                << setw(12) << Z[m].real() << " | "
                << setw(12) << Z[m].imag() << " | "
                << setw(15) << amp << " | "
                << setw(12) << phase << endl;
        }
    }
    cout << defaultfloat;
}

// Обнуление высокочастотных шумовых компонентов
vector<Complex> removeHighFrequencyNoise(const vector<Complex>& Z) {
    int N = static_cast<int>(Z.size());
    vector<Complex> Z_filtered = Z;

    int m_noise1 = 197;
    int m_noise2 = N - 197; 

    Z_filtered[m_noise1] = Complex(0, 0);
    Z_filtered[m_noise2] = Complex(0, 0);

    int bandwidth = 2; // ±2 отсчета
    for (int delta = -bandwidth; delta <= bandwidth; ++delta) {
        int m1 = m_noise1 + delta;
        int m2 = m_noise2 + delta;

        if (m1 >= 0 && m1 < N) Z_filtered[m1] = Complex(0, 0);
        if (m2 >= 0 && m2 < N) Z_filtered[m2] = Complex(0, 0);
    }

    return Z_filtered;
}

// Сохранение данных для построения графиков
void saveToFile(const string& filename, const vector<Complex>& signal1, const vector<Complex>& signal2 = vector<Complex>()) {
    ofstream file(filename);
    if (!file.is_open()) {
        cerr << "Ошибка открытия файла: " << filename << endl;
        return;
    }

    int N = static_cast<int>(signal1.size());
    for (int i = 0; i < N; ++i) {
        file << i << " " << signal1[i].real();
        if (!signal2.empty() && static_cast<int>(signal2.size()) == N) {
            file << " " << signal2[i].real();
        }
        file << endl;
    }
    file.close();
}

// Основная функция
int main() {
    setlocale(0, "");

    int n = 9;
    int N = static_cast<int>(pow(2, n));

    double A = 2.94;
    double omega1 = 2.0;
    double phi = PI / 2;
    double B = 0.27;
    double omega2 = 197.0;

    // Путь для сохранения файлов 
    string desktopPath = "C:\\Users\\DELL\\Desktop\\график\\";
    _mkdir(desktopPath.c_str());

    // Генерация сигнала
    vector<Complex> z = generateSignal(N, A, omega1, phi, B, omega2);

    // DFT и FFT с замером времени
    cout << "   DFT... ";
    auto start = chrono::high_resolution_clock::now();
    vector<Complex> Z_dft = DFT(z);
    auto end = chrono::high_resolution_clock::now();
    auto dft_time = chrono::duration_cast<chrono::microseconds>(end - start).count();
    cout << dft_time << " мкс\n";

    cout << "   FFT... ";
    start = chrono::high_resolution_clock::now();
    vector<Complex> Z_fft_copy = z;
    FFT(Z_fft_copy);
    end = chrono::high_resolution_clock::now();
    auto fft_time = chrono::duration_cast<chrono::microseconds>(end - start).count();
    cout << fft_time << " мкс\n";

    cout << "\n   Сравнение производительности:\n";
    cout << "   DFT: " << dft_time << " микросекунд\n";
    cout << "   FFT: " << fft_time << " микросекунд\n";
    if (fft_time > 0) {
        double speedup = static_cast<double>(dft_time) / fft_time;
        cout << "   Ускорение FFT: " << speedup << " раз\n";
    }

    // Таблицы спектров
    printSpectrumTable(z, Z_dft, "Спектр DFT");
    printSpectrumTable(z, Z_fft_copy, "Спектр FFT");

    // Проверка, что DFT и FFT дают одинаковые результаты
    double max_diff = 0.0;
    for (int i = 0; i < N; ++i) {
        double diff = abs(Z_dft[i] - Z_fft_copy[i]);
        if (diff > max_diff) max_diff = diff;
    }
    cout << "Максимальная разница между DFT и FFT: " << max_diff << endl;
    if (max_diff < 1e-10) {
    }

    // Обнуление конкретных шумовых компонентов
    vector<Complex> Z_filtered = removeHighFrequencyNoise(Z_fft_copy);

    // Табличка спектра после фильтрации
    printSpectrumTable(z, Z_filtered, "Спектр после фильтрации");

    vector<Complex> z_filtered_copy = Z_filtered;
    IFFT(z_filtered_copy);

    // Сохранение спектров для сравнения
    saveToFile(desktopPath + "spectrum_before.txt", Z_fft_copy);
    saveToFile(desktopPath + "spectrum_after.txt", Z_filtered);

    // Сохранение сигналов
    string file1 = desktopPath + "signal_data.txt";
    saveToFile(file1, z, z_filtered_copy);

    // Генерация сигнала с разрывами
    vector<Complex> z2 = generateSignal2(N, A, B, omega2);

    // Сохранение данных для второй задачи
    string file2 = desktopPath + "signal2_data.txt";
    saveToFile(file2, z2);

    return 0;
}

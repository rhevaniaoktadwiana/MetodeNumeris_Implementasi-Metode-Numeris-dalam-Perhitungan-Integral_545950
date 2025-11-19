#include <iostream>
#include <cmath>
#include <iomanip>
#include <vector>

// --- Fungsi yang Akan Diintegrasikan: f(x) = x^2 + 1 ---
double f_integral(double x) {
    return x * x + 1.0;
}

// =======================================================
//                   1. METODE TRAPEZOID
// =======================================================
double trapezoid_rule(double (*f)(double), double a, double b, int n) {
    double h = (b - a) / n;
    double integral = 0.5 * (f(a) + f(b));
    
    for (int i = 1; i < n; ++i) {
        double x_i = a + i * h;
        integral += f(x_i);
    }
    
    integral *= h;
    return integral;
}

// =======================================================
//         2. METODE RECURSIVE TRAPEZOID (ADAPTIF)
// =======================================================
double adaptive_trapezoid(double (*f)(double), double a, double b, double tolerance, int max_depth, int depth = 1) {
    // T1: Perkiraan Trapezoid pada interval [a, b]
    double T1 = (b - a) * (f(a) + f(b)) / 2.0;
    
    // T2: Perkiraan Trapezoid pada 2 sub-interval
    double c = (a + b) / 2.0;
    double h = (b - a) / 2.0;
    double T2 = 0.5 * T1 + h * f(c);
    
    // Perkiraan kesalahan (Error Estimate)
    double error_estimate = std::abs(T2 - T1) / 3.0;
    
    if (error_estimate < tolerance || depth >= max_depth) {
        return T2;
    } else {
        // Rekursi untuk sub-interval
        double left_integral = adaptive_trapezoid(f, a, c, tolerance, max_depth, depth + 1);
        double right_integral = adaptive_trapezoid(f, c, b, tolerance, max_depth, depth + 1);
        
        return left_integral + right_integral;
    }
}

// =======================================================
//                   3. METODE ROMBERG
// =======================================================
double romberg_integration(double (*f)(double), double a, double b, int k_max) {
    // R[i][j] adalah tabel Romberg
    std::vector<std::vector<double>> R(k_max, std::vector<double>(k_max));

    for (int i = 0; i < k_max; ++i) {
        int n = 1 << i; // n = 2^i
        
        // Kolom Pertama (Recursive Trapezoid)
        if (i == 0) {
            R[i][0] = (b - a) * (f(a) + f(b)) / 2.0;
        } else {
            double h = (b - a) / n;
            double sum_midpoints = 0.0;
            // Hanya hitung titik tengah baru (indeks ganjil)
            for (int j = 1; j < n; j += 2) {
                sum_midpoints += f(a + j * h);
            }
            R[i][0] = 0.5 * R[i-1][0] + h * sum_midpoints;
        }

        // Ekstrapolasi Richardson
        for (int j = 1; j <= i; ++j) {
            double power_of_4 = std::pow(4.0, j);
            // Rumus: R(i, j) = R(i, j-1) + [R(i, j-1) - R(i-1, j-1)] / (4^j - 1)
            R[i][j] = R[i][j-1] + (R[i][j-1] - R[i-1][j-1]) / (power_of_4 - 1.0);
        }
    }
    
    return R[k_max - 1][k_max - 1];
}

// =======================================================
//         4. METODE GAUSSIAN QUADRATURE (3-Titik)
// =======================================================
double gaussian_quadrature_3pt(double (*f)(double), double a, double b) {
    // 1. Simpul (Nodes) t_i di interval [-1, 1]
    // Ini adalah akar dari Polinomial Legendre P_3(t) = 0
    const double t1 = -std::sqrt(3.0 / 5.0);
    const double t2 = 0.0;
    const double t3 = std::sqrt(3.0 / 5.0);

    // 2. Bobot (Weights) w_i
    const double w1 = 5.0 / 9.0;
    const double w2 = 8.0 / 9.0;
    const double w3 = 5.0 / 9.0;

    // 3. Transformasi Batas (Scaling dan Shifting factors)
    double c1 = (b - a) / 2.0; // Faktor skala
    double c2 = (b + a) / 2.0; // Faktor pergeseran
    
    // 4. Transformasi Simpul t_i ke Simpul x_i di [a, b]
    double x1 = c1 * t1 + c2;
    double x2 = c1 * t2 + c2;
    double x3 = c1 * t3 + c2;
    
    // 5. Rumus Kuadratur Umum (Termasuk faktor skala c1)
    // Integral = c1 * [w1*f(x1) + w2*f(x2) + w3*f(x3)]
    double integral = c1 * (w1 * f(x1) + w2 * f(x2) + w3 * f(x3));
    
    return integral;
}

// =======================================================
//                        MAIN PROGRAM
// =======================================================
int main() {
    // --- Pengaturan Umum ---
    double a = 0.0;
    double b = 2.0;
    
    // Solusi Eksak: 14/3
    double exact_solution = 14.0 / 3.0;

    std::cout << std::fixed << std::setprecision(10);
    std::cout << "=================================================" << std::endl;
    std::cout << "Integral dari f(x) = x^2 + 1 dari 0 sampai 2" << std::endl;
    std::cout << "Solusi Eksak: " << exact_solution << std::endl;
    std::cout << "=================================================" << std::endl;
    
    // --- 1. Jalankan Metode Trapezoid ---
    int n = 8;
    double result_trap = trapezoid_rule(f_integral, a, b, n);
    std::cout << "1. Trapezoid Rule (n=" << n << "): " << result_trap << std::endl;

    // --- 2. Jalankan Metode Recursive Trapezoid ---
    double tolerance = 1e-9;
    int max_depth = 10;
    double result_recursive = adaptive_trapezoid(f_integral, a, b, tolerance, max_depth);
    std::cout << "2. Recursive Trapezoid (max_depth=" << max_depth << "): " << result_recursive << std::endl;

    // --- 3. Jalankan Metode Romberg ---
    int k_max = 5;
    double result_romb = romberg_integration(f_integral, a, b, k_max);
    std::cout << "3. Romberg Integration (R(4,4)):" << result_romb << std::endl;

    // --- 4. Jalankan Metode Gaussian Quadrature ---
    double result_gauss = gaussian_quadrature_3pt(f_integral, a, b);
    std::cout << "4. Gaussian Quadrature (3-point): " << result_gauss << std::endl;
    
    std::cout << "=================================================" << std::endl;

    return 0;
}
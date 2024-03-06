#include <iostream>
#include <vector>
#include <complex>
#include <numbers>
#include <iomanip>
#include <algorithm>
#include <random>
#include <chrono>

using comp = std::complex<double>;
namespace dft {
    constexpr double pi = std::numbers::pi;

    // Only for testing. DO NOT use this o(n^2) algorithm!!!
    // sgn=-1: dft, sgn=1: idft
    void _dft_test(std::vector<comp>& x, int sgn = -1) {
        int n = x.size();
        const comp root = std::polar(1.0, 2 * pi * sgn / n);
        std::vector<comp> temp(n);

        for (auto [i, base] = std::tuple{0, comp(1)}; i < n; ++i, base *= root) {
            for (auto [j, angle] = std::tuple{0, comp(1)}; j < n; ++j, angle *= base) {
                temp[i] += x[j] * angle;
            }
        }
        std::swap(temp, x);
    }

    // Only for testing. DO NOT use this o(n^2) algorithm!!!
    void _idft_test(std::vector<comp>& x) {
        // conj version
        // std::transform(x.begin(), x.end(), x.begin(), [](const auto &c) {
        //     return std::conj(c); 
        // });
        // _dft_test(x);
        // const auto norm = 1.0 / x.size();
        // std::transform(x.begin(), x.end(), x.begin(), [&](const auto &c) {
        //     return norm * std::conj(c); 
        // });
        _dft_test(x, 1);
        const auto norm = 1.0 / x.size();
        std::transform(x.begin(), x.end(), x.begin(), [&](const auto &c){
            return norm * c;
        });
    }

    void fft(std::vector<comp>& x) {

    }
};

int main() {
    std::mt19937 rng(std::chrono::steady_clock::now().time_since_epoch().count());
    std::uniform_real_distribution<double> dist(0, 10000);

    int N = 10000;
    std::vector<double> A(N);
    std::generate(A.begin(), A.end(), [&](){ return dist(rng); });

    std::vector<comp> B(N);
    for (int i = 0; i < N; ++i) {
        B[i] = A[i];
    }
    auto t0 = std::chrono::high_resolution_clock::now();
    dft::_dft_test(B);
    dft::_idft_test(B);
    auto t1 = std::chrono::high_resolution_clock::now();
    auto time = std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0);
    std::cout << "Operation takes " << time.count() << " seconds.\n";

    for (int i = 0; i < N; ++i) {
        if (std::real(B[i]) - A[i] >= 1e-4 || std::imag(B[i]) >= 1e-4) {
            std::cout << "Test fail!!!\n";
            exit(-1);
        }
    }
    
    std::cout << "Test pass!\n";
    return 0;
}
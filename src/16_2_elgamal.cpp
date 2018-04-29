#define DO_PRINT

#include <math/math.hpp>

using namespace math;

template <typename T, typename U, typename V>
T sm(T b, U e, V m) {
  return sqm<true>(b, e, m);
}

void elgamal(unsigned _p, unsigned _a, unsigned _d, unsigned _i, unsigned _x) {
  trace("Verschlüsseln: \\\\\n");
  const auto al = aligned();

  variable p{_p, "p"}, a{_a, "\\alpha"}, d{_d, "d"}, i{_i, "i"}, x{_x, "x"};

  const auto kpub_b = rename(sm(a, d, p), "k_{pub, B}");
  al.b();
  const auto kpub_a = rename(sm(a, i, p), "k_{pub, A}");
  al.b();
  const auto k_m = rename(sm(kpub_b, i, p), "k_M");
  al.b();
  const auto y = rename((x * k_m) % p, "y");
  al.close();

  trace("Entschlüsseln: \\\\\n");
  al.open();
  const auto k_m2 = rename(sm(kpub_a, d, p), "k_M");
  al.close();

  const auto [_, _t, __] = eea<true, std::int64_t>(p.v(), k_m2.v());

  al.open();
  auto t = variable{_t, "k_M^{-1}"};
  if (t < 0) {
    t = rename((t + p) % p, "k_M^{-1}");
  }
  const auto x_2 = rename((y * t) % p, "x");
}

int main() {
  trace("\\item\n");
  elgamal(383, 20, 175, 82, 137);
  trace("\\item\n");
  elgamal(383, 20, 175, 311, 137);
}

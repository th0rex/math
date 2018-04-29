#define DO_PRINT

#include <math/math.hpp>

using namespace math;

void elgamal(unsigned _p, unsigned _a, unsigned _d, unsigned _i, unsigned _x) {
  trace("Verschlüsseln: \\\\\n");
  const auto al = aligned();

  variable p{_p, "p"}, a{_a, "alpha"}, d{_d, "d"}, i{_i, "i"}, x{_x, "d"};

  const auto kpub_b = rename(z_exp(a, d, p), "k_{pub, B}");
  al.b();
  const auto kpub_a = rename(z_exp(a, i, p), "k_{pub, A}");
  al.b();
  const auto k_m = rename(z_exp(kpub_b, i, p), "k_M");
  al.b();
  const auto y = rename((x * k_m) % p, "y");
  al.close();

  trace("Entschlüsseln: \\\\\n");
  al.open();
  const auto k_m2 = rename(z_exp(kpub_a, d, p), "k_M");
  al.close();

  auto [_, _t, __] = eea<True>(p.v(), k_m2.v());
  auto t = rename(_t, "k_M^{-1}");
  if (t < 0) {
    t = rename((t + p) % p, "k_M^{-1}");
  }
  al.open();
  const auto x_2 = (y * t) % p;
}

int main() {
  trace("\\item\n");
  elgamal(383, 20, 175, 82, 137);
  trace("\\item\n");
  elgamal(383, 20, 175, 311, 137);
}

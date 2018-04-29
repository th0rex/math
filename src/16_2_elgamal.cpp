#define DO_PRINT

#include <math/math.hpp>

using namespace math;

template <typename T, typename U, typename V>
T sm(T b, U e, V m) {
  return sqm<true>(b, e, m);
}

// TODO: put this into som generic thing

template <typename T>
auto var(T t, const char* s) {
  return variable<T, format_expr>{t, s};
}

template <typename T>
auto var(variable<T, format_expr> const& v, const char* s) {
  return variable{v, s};
}

template <op o, template <typename, bool> typename F, typename L, typename R>
auto var(expression<o, F, L, R> e, const char* s) {
  return variable{std::move(e), s};
}

void elgamal(unsigned _p, unsigned _a, unsigned _d, unsigned _i, unsigned _x) {
  trace("Verschlüsseln: \\\\\n");
  const auto al = aligned();

  const auto p = var(_p, "P");
  const auto a = var(_a, "\\alpha");
  const auto d = var(_d, "d");
  const auto i = var(_i, "i");
  const auto x = var(_x, "x");

  const auto kpub_b = var(sm(a, d, p), "k_{pub, B}");
  al.b();
  const auto kpub_a = var(sm(a, i, p), "k_{pub, A}");
  al.b();
  const auto k_m = var(sm(kpub_b, i, p), "k_M");
  al.b();
  const auto y = var((x * k_m) % p, "y");
  al.close();

  trace("Entschlüsseln: \\\\\n");
  al.open();
  const auto k_m2 = var(sm(kpub_a, d, p), "k_M");
  al.close();

  const auto [_, _t, __] = eea<true, std::int64_t>(p.v(), k_m2.v());

  al.open();
  auto t = var(_t, "k_M^{-1}");
  if (t < 0) {
    t = var((t + p) % p, "k_M^{-1}");
  }
  const auto x_2 = var((y * t) % p, "x");
}

int main() {
  trace("\\item\n");
  elgamal(383, 20, 175, 82, 137);
  trace("\\item\n");
  elgamal(383, 20, 175, 311, 137);
}

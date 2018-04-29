#define DO_PRINT

#include <math/math.hpp>

using namespace math;

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

int main() {
  {
    const auto _ = aligned{};
    auto v = var(0, "i");
    v += 1;
    v = v * v + v * 5;
    v %= 3;

    auto x = v;
    v *= 2;
    v = x + 5;
  }

  std::cout << "\n\n";

  const auto x = sqm(5u, 10u, 13u);
  std::cout << "Square multiply: \\\\\n";
  {
    const auto _ = aligned{};
    const auto y =
        var(sqm(var(5u, "alpha"), var(10u, "d"), var(13u, "p")), "k_{pub, B}");
  }

  std::cout << baby_step_giant_step(17, 101, 167) << "\n";

  return 0;
}

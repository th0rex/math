#define DO_PRINT

#include <math/math.hpp>

using namespace math;

auto var = make_var_factory<format_expr>();

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

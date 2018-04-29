#define DO_PRINT

#include <math/math.hpp>

using namespace math;

int main() {
  {
    const auto _ = aligned{};
    variable v{0, "i"};
    v += 1;
    v = v * v + v * 5;
    v %= 3;

    variable x = v;
    v *= 2;
    v = x + 5;
  }

  std::cout << "\n\n";

  const auto x = sqm(5u, 10u, 13u);
  std::cout << "Square multiply: \\\\\n";
  {
    const auto _ = aligned{};
    const auto y = rename(
        sqm(variable{5u, "alpha"}, variable{10u, "d"}, variable{13u, "p"}),
        "k_{pub, B}");
  }

  std::cout << baby_step_giant_step(17, 101, 167) << "\n";

  return 0;
}

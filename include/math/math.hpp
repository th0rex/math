#pragma once

#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <optional>
#include <type_traits>
#include <unordered_map>
#include <utility>

namespace math {
template <typename T>
std::ostream &operator<<(std::ostream &os, std::optional<T> const &x) {
  if (x) {
    return os << *x;
  }
  return os << "(none)";
}

unsigned indent = 0;
bool is_aligned = false;

#ifdef DO_PRINT
std::ostream *out = &std::cout;

template <bool First = false, typename... T>
void trace(T &&... ts) {
  if constexpr (First) {
    for (auto i = 0; i < indent * 2; ++i) {
      std::cout << ' ';
    }
  }
  (*out << ... << ts);
}
#else
template <typename... T>
void trace(T &&... ts) {}
#endif

template <typename T, typename = std::void_t<>>
struct get_real_type {
  using type = T;
};

template <typename T>
struct get_real_type<T, std::void_t<typename T::real_type>> {
  using type = typename T::real_type;
};

template <typename T>
using real_type_t = typename get_real_type<T>::type;

template <typename T>
struct variable;

enum class op { ADD, SUB, MUL, DIV, MOD };

template <op o, typename L, typename R>
struct expression;

template <op o, typename L, typename R>
auto make_expression(L l, R r) {
  return expression<o, L, R>{std::move(l), std::move(r)};
}

template <op o, typename L, typename R>
struct expression {
  using real_type = std::common_type_t<real_type_t<L>, real_type_t<R>>;

  L _lhs;
  R _rhs;

  expression(L lhs, R rhs) : _lhs{std::move(lhs)}, _rhs{std::move(rhs)} {}

  expression(expression const &) = delete;
  expression(expression &&) = default;

  template <typename B>
  auto operator+(B rhs) {
    return make_expression<op::ADD>(std::move(*this), std::move(rhs));
  }
  template <typename B>
  auto operator-(B rhs) {
    return make_expression<op::SUB>(std::move(*this), std::move(rhs));
  }
  template <typename B>
  auto operator*(B rhs) {
    return make_expression<op::MUL>(std::move(*this), std::move(rhs));
  }
  template <typename B>
  auto operator/(B rhs) {
    return make_expression<op::DIV>(std::move(*this), std::move(rhs));
  }
  template <typename B>
  auto operator%(B rhs) {
    return make_expression<op::MOD>(std::move(*this), std::move(rhs));
  }
};

template <typename T>
struct get_value {
  constexpr static auto value(T t) { return t; }
};

template <typename L, typename R>
struct get_value<expression<op::ADD, L, R>> {
  constexpr static auto value(expression<op::ADD, L, R> e) {
    return get_value<L>::value(std::move(e._lhs)) +
           get_value<R>::value(std::move(e._rhs));
  }
};
template <typename L, typename R>
struct get_value<expression<op::SUB, L, R>> {
  constexpr static auto value(expression<op::SUB, L, R> e) {
    return get_value<L>::value(std::move(e._lhs)) -
           get_value<R>::value(std::move(e._rhs));
  }
};
template <typename L, typename R>
struct get_value<expression<op::MUL, L, R>> {
  constexpr static auto value(expression<op::MUL, L, R> e) {
    return get_value<L>::value(std::move(e._lhs)) *
           get_value<R>::value(std::move(e._rhs));
  }
};
template <typename L, typename R>
struct get_value<expression<op::DIV, L, R>> {
  constexpr static auto value(expression<op::DIV, L, R> e) {
    return get_value<L>::value(std::move(e._lhs)) /
           get_value<R>::value(std::move(e._rhs));
  }
};
template <typename L, typename R>
struct get_value<expression<op::MOD, L, R>> {
  constexpr static auto value(expression<op::MOD, L, R> e) {
    return get_value<L>::value(std::move(e._lhs)) %
           get_value<R>::value(std::move(e._rhs));
  }
};

template <typename T, typename = void>
struct format_expr;

template <typename T>
struct format_expr<
    T, std::enable_if_t<std::is_integral_v<T> || std::is_floating_point_v<T>>> {
  constexpr static void format(T const &t) { trace(t); }
};

template <typename T, typename = void>
struct format_delimited {
  constexpr static void do_format(T const &t, const char *l, const char *r) {
    trace(l);
    format_expr<T>::format(t);
    trace(r);
  }
};

template <typename T>
struct format_delimited<
    T, std::enable_if_t<std::is_integral_v<T> || std::is_floating_point_v<T>>> {
  constexpr static void do_format(T const &t, const char *, const char *) {
    format_expr<T>::format(t);
  }
};

template <typename T>
struct format_delimited<variable<T>, void> {
  constexpr static void do_format(variable<T> const &t, const char *,
                                  const char *) {
    format_expr<variable<T>>::format(t);
  }
};

template <typename T>
void format_paren(T const &t) {
  format_delimited<T>::do_format(t, "\\left(", "\\right)");
}

template <typename L, typename R>
struct format_expr<expression<op::ADD, L, R>, void> {
  constexpr static void format(expression<op::ADD, L, R> const &e) {
    format_paren(e._lhs);
    trace(" + ");
    format_paren(e._rhs);
  }
};
template <typename L, typename R>
struct format_expr<expression<op::SUB, L, R>, void> {
  constexpr static void format(expression<op::SUB, L, R> const &e) {
    format_paren(e._lhs);
    trace("-");
    format_paren(e._rhs);
  }
};
template <typename L, typename R>
struct format_expr<expression<op::MUL, L, R>, void> {
  constexpr static void format(expression<op::MUL, L, R> const &e) {
    format_paren(e._lhs);
    trace("*");
    format_paren(e._rhs);
  }
};
template <typename L, typename R>
struct format_expr<expression<op::DIV, L, R>, void> {
  constexpr static void format(expression<op::DIV, L, R> const &e) {
    format_paren(e._lhs);
    trace("/");
    format_paren(e._rhs);
  }
};
template <typename L, typename R>
struct format_expr<expression<op::MOD, L, R>, void> {
  constexpr static void format(expression<op::MOD, L, R> const &e) {
    format_paren(e._lhs);
    trace(" \\bmod ");
    format_paren(e._rhs);
  }
};

template <typename T>
struct get_value<variable<T>> {
  constexpr static typename variable<T>::real_type value(variable<T> const &v);
};

template <typename T>
struct variable {
  using real_type = T;

  T _value;
  const char *_name;
  std::uint64_t _counter;

  void format() {
    trace<true>(_name, is_aligned ? " &= " : " = ", _value, "\\\\\n");
  }

  variable(T &&value)
      : _value{std::forward<T>(value)}, _name{nullptr}, _counter{0} {}

  variable(T &&value, const char *name)
      : _value{std::forward<T>(value)}, _name{name}, _counter{0} {
    format();
  }

  variable(T value, const char *name)
      : _value{std::move(value)}, _name{name}, _counter{0} {
    format();
  }

  variable(variable const &o, const char *name)
      : _value{o._value}, _name{name}, _counter{0} {
    format();
  }

  variable(variable const &) = default;
  variable(variable &&) = default;

  template <typename A>
  auto operator+(A a) {
    return make_expression<op::ADD>(*this, std::move(a));
  }
  template <typename A>
  auto operator-(A a) {
    return make_expression<op::SUB>(*this, std::move(a));
  }
  template <typename A>
  auto operator*(A a) {
    return make_expression<op::MUL>(*this, std::move(a));
  }
  template <typename A>
  auto operator/(A a) {
    return make_expression<op::DIV>(*this, std::move(a));
  }
  template <typename A>
  auto operator%(A a) {
    return make_expression<op::MOD>(*this, std::move(a));
  }

  template <typename A>
  using condition = std::enable_if_t<std::is_same_v<real_type_t<A>, T>>;

  template <typename A, typename = condition<A>>
  variable &operator=(A a) {
    assert(_name && "can't assign to temporary variables");
    _counter++;
    trace<true>(_name, "_{", _counter, is_aligned ? "} &= " : "} = ");
    format_expr<A>::format(a);
    _value = get_value<A>::value(std::move(a));
    trace(" = ", _value, "\\\\\n");
    return *this;
  }

  template <typename A, typename = condition<A>>
  variable &operator+=(A a) {
    return *this = make_expression<op::ADD>(*this, std::move(a));
  }

  template <typename A, typename = condition<A>>
  variable &operator-=(A a) {
    return *this = make_expression<op::SUB>(*this, std::move(a));
  }

  template <typename A, typename = condition<A>>
  variable &operator*=(A a) {
    return *this = make_expression<op::MUL>(*this, std::move(a));
  }

  template <typename A, typename = condition<A>>
  variable &operator/=(A a) {
    return *this = make_expression<op::DIV>(*this, std::move(a));
  }

  template <typename A, typename = condition<A>>
  variable &operator%=(A a) {
    return *this = make_expression<op::MOD>(*this, std::move(a));
  }

  bool operator==(const T &o) const { return _value == o; }

  bool operator!=(const T &o) const { return !(*this == o); }

  bool operator<(const T &o) const { return _value < o; }
  // TODO: more operators

  explicit operator real_type() const { return _value; }

  real_type v() const { return _value; }
};

template <typename T>
variable(T, const char *)->variable<T>;

template <typename T>
struct format_expr<variable<T>> {
  constexpr static void format(variable<T> const &v) {
    if (v._counter == 0) {
      trace(v._name);
    } else {
      trace(v._name, "_{", v._counter, "}");
    }
  }
};

template <typename T>
constexpr typename variable<T>::real_type get_value<variable<T>>::value(
    variable<T> const &v) {
  return v._value;
}

template <typename T>
struct is_variable : std::false_type {};

template <typename T>
struct is_variable<variable<T>> : std::true_type {};

template <typename T>
constexpr bool is_variable_v = is_variable<std::remove_reference_t<T>>::value;

template <typename T>
auto rename(T &&value, const char *name) {
  if constexpr (is_variable_v<T>) {
    return variable{value, name};
  }
  return std::forward<T>(value);
}

struct aligned {
  aligned() { open(); }

  void open() const {
    trace("$\n\\begin{aligned}[t]\\\\\n");
    indent += 1;
    is_aligned = true;
  }

  // break an aligned block
  void b() const {
    if (is_aligned) {
      close();
      open();
    }
  }

  void close() const {
    if (is_aligned) {
      indent -= 1;
      trace("\\end{aligned}\n$\n");
      is_aligned = false;
    }
  }

  aligned(aligned &&) = default;

  ~aligned() { close(); }
};

// TODO: implement Z_{n}

template <typename T>
T z_exp(T b, T e, T m) {
  using real_type = real_type_t<T>;
  static_assert(std::is_integral_v<real_type>,
                "z_exp only works for integrals");
  constexpr auto bit_size = sizeof(real_type) * 8;

  if (e == 0) {
    return {1};
  }

  T acc = rename(b, "tmp");

  unsigned i = __builtin_clz((real_type)e) + 1;
  while (i < bit_size) {
    acc = (acc * acc) % m;
    if ((real_type)e & (1 << (bit_size - 1 - i++))) {
      acc = (acc * b) % m;
    }
  }

  return acc;
}

template <typename T>
struct eea_result {
  T s;
  T t;
  T gcd;
};

template <typename T>
struct h {
  T &a0;
  T &a1;
  T &a2;

  h(T *arr, std::uint64_t offset)
      : a0{arr[offset % 3]},
        a1{arr[(offset + 2) % 3]},
        a2{arr[(offset + 1) % 3]} {}
};

template <typename T>
h(T *, std::uint64_t)->h<T>;

template <bool Trace = false, typename T>
eea_result<T> eea(T r0, T r1) noexcept {
  static_assert(std::is_integral_v<T>, "eea only works for integrals");

  if (r0 <= r1) {
    std::swap(r0, r1);
  }

  T r[3] = {r0, r1, 0};
  T s[3] = {1, 0, 0};
  T t[3] = {0, 1, 0};

  if constexpr (Trace) {
    trace(
        "\\begin{tabular}{|c|c|c|c|c|}\n\\hline\ni & r & q & s & t \\\\ "
        "\\hline \\hline\n");
    for (auto i = 0; i < 2; ++i) {
      trace(i, " & ", r[i], " & & ", s[i], " & ", t[i], "\\\\ \\hline \n");
    }
  }

  std::uint64_t i = 2;
  for (; r[(i + 2) % 3] != 0; ++i) {
    // This might be oddly named.
    // x0 = x_i, x1 = x_{i-1}, x2 = x_{i-2}
    auto [r0, r1, r2] = h{r, i};
    auto [s0, s1, s2] = h{s, i};
    auto [t0, t1, t2] = h{t, i};

    r0 = r2 % r1;

    const auto q = (r2 - r0) / r1;

    s0 = s2 - q * s1;
    t0 = t2 - q * t1;

    if constexpr (Trace) {
      trace(i, " & ", r0, " & ", q, " & ", s0, " & ", t0, "\\\\ \\hline \n");
    }
  }

  if constexpr (Trace) {
    trace("\\end{tabular}\\\\\n");
  }

  return {s[(i + 1) % 3], t[(i + 1) % 3], r[(i + 1) % 3]};
}

template <typename T>
// returns dlog_a(b) mod p
std::optional<T> baby_step_giant_step(T a, T b, T p) {
  static_assert(std::is_integral_v<T>,
                "baby_step_giant_step only works for integrals");

  const auto m = (T)ceil(sqrt(p));
  auto [_, a_inv, __] = eea(a, p);
  a_inv = (a_inv + p) % p;
  const auto a_inv_m = z_exp(a_inv, m, p);

  // a^{x_b} to x_b
  std::unordered_map<T, T> lookup;
  for (T i = 0; i < m; ++i) {
    lookup[z_exp(a, i, p)] = i;
  }

  for (T i = 0; i < m; ++i) {
    const auto res = (z_exp(a_inv_m, i, p) * b) % p;
    if (const auto r = lookup.find(res); r != lookup.end()) {
      return {i * m + r->second};
    }
  }

  return {};
}

}  // namespace math

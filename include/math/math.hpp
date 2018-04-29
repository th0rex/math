#pragma once

#include <cassert>
#include <cstdint>
#include <iostream>
#include <type_traits>
#include <utility>

namespace math {
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

  explicit operator real_type() const { return _value; }
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
  aligned() {
    trace("$\n\\begin{aligned}[t]\\\\\n");
    indent += 1;
    is_aligned = true;
  }

  aligned(aligned &&) = default;

  ~aligned() {
    indent -= 1;
    trace("\\end{aligned}\n$\n");
    is_aligned = false;
  }
};

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

}  // namespace math

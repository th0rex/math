#define DO_PRINT

#include <math/math.hpp>

using namespace math;

template <typename T, typename U, typename V>
T sm(T b, U e, V m) {
  return sqm<true>(b, e, m);
}

auto var = make_var_factory<format_expr_with_value>();

template <typename T, typename U>
struct cert {
  T x;
  U sig;
};

template <typename T>
struct elgamal_sig {
  T r;
  T s;
};

template <typename E, typename T>
auto elgamal_sign(E& env, T x, T k_e, T d, T a, T p) {
  auto r = var(sm(a, k_e, p), "r");

  env.close();

  const auto p_1 = var(p - 1, "p''");
  const auto [_, _k_e_inv, __] = eea<true, std::int32_t>(k_e.v(), p_1.v());

  env.open();
  auto k_e_inv = var(_k_e_inv, "k_E^{-1}");
  if (k_e_inv < 0) {
    k_e_inv = var((k_e_inv + (p - 1)) % (p - 1), "k_E^{-1}");
  }

  auto s = var(((x - d * r) * k_e_inv) % (p - 1), "s");
  if (s < 0) {
    s = var((s + (p - 1)) % (p - 1), "s");
  }
  return elgamal_sig<T>{r, s};
}

template <typename E, typename T>
auto make_cert(E& env, T dh_a, T dh_priv, T dh_p, T id, T k_e, T p, T d, T a) {
  auto k_pub = var(sm(dh_a, dh_priv, dh_p), "k_{pub,i}");
  auto x = var(k_pub + id, "x_i");
  return cert<T, elgamal_sig<T>>{x, elgamal_sign(env, x, k_e, d, a, p)};
}

template <typename E, typename T>
bool elgamal_verify(E& env, T x, elgamal_sig<T> sig, T a, T b, T p) {
  const auto t_1 = var(sm(b, sig.r, p), "t_1");
  const auto t_2 = var(sm(sig.r, sig.s, p), "t_2");
  auto t = var((t_1 * t_2) % p, "t");
  env.b();
  auto u = var(sm(a, x, p), "u");
  return t == u.v();
}

template <typename T, typename U>
void print_cert(const char* name, cert<T, U> const& cert) {
  trace("Cert_{", name, "}\\,\\text{:}\\\\");
  trace("x &= ", cert.x.v(), "\\\\");
  trace("sig &= [", cert.sig.r.v(), "\\,\\,\\text{,}\\,\\,", cert.sig.s.v(),
        "]\\\\");
}

template <typename E, typename C, typename T>
void verify_cert(E& env, C const& cert, T a, T b, T p) {
  if (elgamal_verify(env, cert.x, cert.sig, a, b, p)) {
    trace("\\text{t == u => Zertifikat ist valide}");
  } else {
    trace("\\text{t != u => Zertifikat ist nicht valide}");
  }
  env.b();
}

int main() {
  trace("\\item\n");
  const auto al = aligned();

  const auto dh_p = var(107, "p");
  const auto dh_alpha = var(28, "\\alpha");

  const auto elgamal_p = var(467, "p'");
  const auto elgamal_alpha = var(45, "\\alpha'");
  const auto elgamal_d = var(378, "d'");

  const auto cert_alice =
      make_cert(al, dh_alpha, var(12, "a"), dh_p, var(65, "id"), var(29, "k_E"),
                elgamal_p, elgamal_d, elgamal_alpha);
  print_cert("A", cert_alice);

  al.b();

  const auto cert_bob =
      make_cert(al, dh_alpha, var(34, "b"), dh_p, var(66, "id"), var(71, "k_E"),
                elgamal_p, elgamal_d, elgamal_alpha);

  print_cert("B", cert_bob);
  al.b();

  const auto cert_charley =
      make_cert(al, dh_alpha, var(56, "c"), dh_p, var(67, "id"),
                var(113, "k_E"), elgamal_p, elgamal_d, elgamal_alpha);

  print_cert("C", cert_charley);
  al.close();

  trace("\\item\n");
  al.open();

  const auto elgamal_beta =
      var(sm(elgamal_alpha, elgamal_d, elgamal_p), "\\beta");
  al.b();

  trace("\\text{Alice:}\\\\");
  verify_cert(al, cert_alice, elgamal_alpha, elgamal_beta, elgamal_p);
  al.b();
  trace("\\text{Bob:}\\\\");
  verify_cert(al, cert_bob, elgamal_alpha, elgamal_beta, elgamal_p);
  al.b();
  trace("\\text{Charley:}\\\\");
  verify_cert(al, cert_charley, elgamal_alpha, elgamal_beta, elgamal_p);
  al.b();
}

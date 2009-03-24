#include <iostream>
#include <fstream>
#include <jlt/math.hpp>
#include <jlt/stlio.hpp>
#include <jlt/exceptions.hpp>
#include <jlt/reciprocal_polynomial.hpp>

#define DOUBLECHECK_SPECTRAL_RADIUS

#ifdef DOUBLECHECK_SPECTRAL_RADIUS
#include <jlt/mathmatrix.hpp>
#include <jlt/eigensystem.hpp>
#endif

namespace jlt {

template<class S, int g, class T>
T findroot(const reciprocal_polynomial<S,g>& p, const T x0, const T tol);

#ifdef DOUBLECHECK_SPECTRAL_RADIUS
// Spectral radius (magnitude of largest eigenvalue) of polynomial.
template<class S, int g, class T>
void spectral_radius(const reciprocal_polynomial<S,g>& p, T& sr);
#endif

template<class T>
bool is_candidate(T const l, T const lmin, T const lmax);

template<class T>
bool increment_vector(std::vector<T>& a,
		      const std::vector<T>& amin,
		      const std::vector<T>& amax);

template<class S, int g>
bool traces_to_reciprocal_poly(const std::vector<S>& Tr,
			       jlt::reciprocal_polynomial<S,g>& p);

template<class T>
bool reject(const std::vector<T>& a);

} // namespace jlt


int main()
{
  using jlt::reciprocal_polynomial;
  using jlt::operator<<;
  using jlt::is_candidate;
  using jlt::Abs;
  using std::cout;
  using std::cerr;
  using std::endl;

  typedef double T;

  const int g = 5;
  long long int N = 0;
  long long int Np = 0; // Number of candidate polynomials found.
  long long int N_found_positive_root = 0;
  long long int N_found_negative_root = 0;
#ifdef DOUBLECHECK_SPECTRAL_RADIUS
  long long int N_discard_positive_by_spectral_radius = 0;
  long long int N_discard_negative_by_spectral_radius = 0;
#endif
  long long int N_failed_to_converge = 0;
  long long int N_fractional_coefficients = 0;
  const T tol = 1e-4;
  T lambdamax;
  // Use bound from below as well: Penner's 2^1/(12g-12).
  const T lambdamin = jlt::Pow(2.0,1.0/(12*g-12));

  std::vector<int> Trmin(g), Trmax(g), Tr(g);

  if (g == 2)
    {
      lambdamax = 1.72208380573904;
    }
  else if (g == 3)
    {
      lambdamax = 1.72208380573904;
    }
  else if (g == 4)
    {
      lambdamax = 1.34371999561122;
    }
  else if (g == 5)
    {
      lambdamax = 1.27247673631011;
    }
  else if (g == 6)
    {
      lambdamax = 1.22571747523471;
    }
  else if (g == 7)
    {
      lambdamax = 1.19266542682899;
    }
  else
    {
      cerr << "What should the maximum root be?\n";
      exit(1);
    }

  for (int k = 0; k < g; ++k)
    {
      Trmax[k] = g * (jlt::Pow(lambdamax,k+1) + 1);
      Trmin[k] = -Trmax[k];
    }
  cerr << "Maximum traces:            " << Trmax << endl;
  long int to_try = Trmax[0]+1;
  for (int k = 1; k < g; ++k) to_try *= 2*Trmax[k]+1;
  cerr << "Cases to try:                  " << to_try << endl;

  reciprocal_polynomial<int,g> p;

  T x0 = 1.2*lambdamax;

  // Initial value for traces: the first ones begin at 0, to take
  // advantage of the sp(p(x))=sp(p(-x)) symmetry of the spectral
  // radius.
  Tr[0] = 0;
  for (int m = 1; m < g; ++m) Tr[m] = Trmin[m];

  std::ofstream ostr("poly.m");
  ostr << "{";

  do
    {
      ++N;

      // traces-to_poly returns false if the current list of traces do
      // not give a polynomial over Z.
      if (!traces_to_reciprocal_poly(Tr,p))
	{
	  ++N_fractional_coefficients;
	  continue;
	}

      bool cand = false;
      T lambda = 0;
      try
	{
	  lambda = jlt::findroot(p,x0,tol);
	  // If lambda < 0, converged to a negative root from a
	  // positive initial guess, so can't possibly have a dominant
	  // positive real root.
	  if (is_candidate(lambda,lambdamin,lambdamax) && lambda > 0)
	    { cand = true; ++N_found_positive_root; }
	}
      catch (jlt::failed_to_converge<T>&) { ++N_failed_to_converge; }

      if (!cand)
	{
	  // If we didn't get a candidate, try to start from negative
	  // x, in case the dominant real eigenvalue is negative.
	  try
	    {
	      lambda = jlt::findroot(p,-x0,tol);
	      // If lambda > 0, converged to a positive root from a
	      // negative initial guess, so can't possibly have a
	      // dominant negative real root.
	      if (is_candidate(lambda,lambdamin,lambdamax) && lambda < 0)
		{ cand = true; ++N_found_negative_root; }
	    }
	  catch (jlt::failed_to_converge<T>&) { ++N_failed_to_converge; }
	}

#ifdef DOUBLECHECK_SPECTRAL_RADIUS
      // The most "wasteful" part is by far dominated by polynomials
      // that have a real root with |lambda|<lambdamax, but a complex
      // root |lambdac|>|lambda|.  Thus, it's worth explicitly finding
      // the spectral radius of the polynomial in this case.
      if (cand)
	{
	  T lambda_sr;
	  jlt::spectral_radius(p,lambda_sr);
	  if (!is_candidate(lambda_sr,lambdamin,lambdamax))
	    {
	      // Looks like the spectral radius is too large after
	      // all: drop this candidate.
	      if (lambda > 0) ++N_discard_positive_by_spectral_radius;
	      if (lambda < 0) ++N_discard_negative_by_spectral_radius;
	      cand = false;
	    }
	  else
	    lambda = lambda_sr;
	}
#endif // DOUBLECHECK_SPECTRAL_RADIUS

      if (cand)
	{
	  if (Np++ != 0) ostr << ",\n";
	  ostr << p.to_polynomial() << std::flush;
	}
    } while (jlt::increment_vector(Tr,Trmin,Trmax));

  ostr << "}\n";
  cerr << "Discarded fractional poly:     ";
  cerr << N_fractional_coefficients;
  cerr << "\nFound positive candidate root: ";
  cerr << N_found_positive_root;
#ifdef DOUBLECHECK_SPECTRAL_RADIUS
  cerr << "\t(" << N_discard_positive_by_spectral_radius << " discarded)";
#endif
  cerr << endl;
  cerr << "Found negative candidate root: ";
  cerr << N_found_negative_root;
#ifdef DOUBLECHECK_SPECTRAL_RADIUS
  cerr << "\t(" << N_discard_negative_by_spectral_radius << " discarded)";
#endif
  cerr << endl;
  cerr << "Failed to converge:            ";
  cerr << N_failed_to_converge << endl;
  cerr << "Found " << Np << "/" << N << " candidates.\n";
}


namespace jlt {

template<class S, int g, class T>
inline T findroot(const reciprocal_polynomial<S,g>& p,
		  const T x0, const T tol)
{
  T px(p(x0)), x(x0);
  int i = 0;
  const int itmax = 10;

  while (Abs(px) > tol && i++ < itmax)
    {
      x = x - px / p.derivative_at(x);
      px = p(x);
    }

  if (i == itmax)
    throw
      failed_to_converge<T>
      ("Failed to converge to specified accuracy.\n",Abs(px));
  else
    return x;
}


#ifdef DOUBLECHECK_SPECTRAL_RADIUS
// Spectral radius (magnitude of largest eigenvalue) of polynomial.
template<class S, int g, class T>
void spectral_radius(const reciprocal_polynomial<S,g>& p, T& sr)
{
  const int n = p.degree();
  mathmatrix<T> M(n,n);

  for (int i = 0; i < n-1; ++i) M(i,i+1) = 1;
  // Copy the coefficients of the polynomial to the matrix.
  for (int m = 0; m < n; ++m) M(n-1,m) = -(T)p[m]/p[n];

  sr = spectral_radius(M);
}
#endif // DOUBLECHECK_SPECTRAL_RADIUS


template<class T>
inline bool is_candidate(T const l, T const lmin, T const lmax)
{
  return (Abs(l) >= lmin && Abs(l) < lmax);
}


template<class T>
inline bool increment_vector(std::vector<T>& a,
			     const std::vector<T>& amin,
			     const std::vector<T>& amax)
{
  //
  // Increment elements of the vector a (from the last element),
  // resetting each position a[m] to amin[m] if it passes amax[m].
  //
  const int n = a.size();
  bool incr = false;
  for (int m = n-1; m >= 0; --m)
    {
      // The sequence that coefficients cycle through is
      // amin[m]...amax[m] for a[m], 0 <= m < n
      if (a[m] < amax[m])
	{
	  // Increment the coefficient.
	  incr = true;
	  ++a[m];
	  break;
	}
      else
	{
	  // Otherwise reset coefficient at that position, and let
	  // the loop move on to the next position.
	  a[m] = amin[m];
	}
    }

  return incr;
}



template<class S, int g>
bool traces_to_reciprocal_poly(const std::vector<S>& Tr,
			       jlt::reciprocal_polynomial<S,g>& p)
{
  // See Silva, J. Math. Phys. 39, 6206 (1998), Theorem 1.
  const int n = 2*g;

  for (int k = 1; k <= g; ++k)
    {
      p[n-k] = -Tr[k-1];
      for (int m = 1; m <= k-1; ++m)
	{
	  p[n-k] -= Tr[k-m-1]*p[m];
	}
      if (p[n-k] % k) return false;
      p[n-k] /= k;
    }

  return true;
}


template<class T>
inline bool reject(const std::vector<T>& a)
{
  // If the first nonzero coefficient of an odd x power is
  // negative then we can skip this case using the
  // sp(p(x))=sp(p(-x)) symmetry of the spectral radius.
  const int n = a.size();
  bool skip = false;
  for (int m = 0; m < n; m += 2)
    {
      if (a[m] > 0) break;
      if (a[m] < 0) { skip = true; break; }
    }
  return skip;
}

} // namespace jlt

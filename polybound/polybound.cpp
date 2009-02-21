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
T findroot(const reciprocal_polynomial<S,g,T>& p, const T x0, const T tol);

#ifdef DOUBLECHECK_SPECTRAL_RADIUS
// Spectral radius (magnitude of largest eigenvalue) of polynomial.
template<class S, int g, class T>
T spectral_radius(const reciprocal_polynomial<S,g,T>& p);
#endif

template<class T>
bool is_candidate(T const l, T const lmin, T const lmax);

template<class T>
bool increment_vector(std::vector<T>& a,
		      const std::vector<T>& amin,
		      const std::vector<T>& amax);

template<class T>
bool reject(const std::vector<T>& a);

} // namespace jlt


int main()
{
  using jlt::reciprocal_polynomial;
  using jlt::Abs;
  using jlt::operator<<;
  using jlt::is_candidate;
  using std::cout;
  using std::cerr;
  using std::endl;

  typedef double T;

  const int g = 4;
  long long int N = 0;
  long long int Np = 0; // Number of candidate polynomials found.
  long long int N_found_positive_root = 0;
  long long int N_found_negative_root = 0;
#ifdef DOUBLECHECK_SPECTRAL_RADIUS
  long long int N_discard_positive_by_spectral_radius = 0;
  long long int N_discard_negative_by_spectral_radius = 0;
#endif
  long long int N_failed_to_converge = 0;
  const T tol = 1e-4;
  T lambdamax;
  // Use bound from below as well: Penner's 2^1/(12g-12).
  const T lambdamin = jlt::Pow(2.0,1.0/(12*g-12));

  std::vector<int> amin(g), amax(g), a(g);

  if (g == 3)
    {
      // 12,765 cases
      amax[0] = 6; amin[0] = 0;
      amax[1] = 18; amin[1] = -amax[1];
      amax[2] = 26; amin[2] = -amax[2];
      lambdamax = 1.72208380573904;
    }
  else if (g == 4)
    {
      // 9,889,930 cases
      amax[0] = 8; amin[0] = 0;
      amax[1] = 30; amin[1] = -amax[1];
      amax[2] = 61; amin[2] = -amax[2];
      amax[3] = 77; amin[3] = -amax[3];
      lambdamax = 1.34371999561122;
    }
  else if (g == 5)
    {
#if 0
      // 54,873,202,455 cases
      amax[0] = 10; amin[0] = 0;
      amax[1] = 46; amin[1] = -amax[1];
      amax[2] = 123; amin[2] = -amax[2];
      amax[3] = 217; amin[3] = -amax[3];
      amax[4] = 261; amin[4] = -amax[4];
      lambdamax = 1.17628081825992;
#endif
      // 63,523,102,800 cases
      amax[0] = 10; amin[0] = 0;
      amax[1] = 47; amin[1] = -amax[1];
      amax[2] = 128; amin[2] = -amax[2];
      amax[3] = 226; amin[3] = -amax[3];
      amax[4] = 273; amin[4] = -amax[4];
      lambdamax = 1.27247673631011;
    }
  else
    {
      cerr << "What should the coefficients be?\n";
      exit(1);
    }

  reciprocal_polynomial<int,g,T> p;

  T x0 = 1.2*lambdamax;

  // Initial value for coefficients: the odd ones begin at 0, to take
  // advantage of the sp(p(x))=sp(p(-x)) symmetry of the spectral
  // radius.  If amin[m]<0, they eventually become negative as they
  // increment past amax[m] and their value resets to amin[m].
  // /* What happens if both amin and amax are negative? */
  for (int m = 0; m < g; m += 2) a[m] = 0;
  for (int m = 1; m < g; m += 2) a[m] = amin[m];

  std::ofstream ostr("poly.m");
  ostr << "{";

  do
    {
      // If the first nonzero coefficient of an odd x power is
      // negative then we can skip this case using the
      // sp(p(x))=sp(p(-x)) symmetry of the spectral radius.
      if (jlt::reject(a)) continue;
      ++N;

      // Print current status to stderr once in a while.
      if (g > 3)
	{
	  bool prnt = true;
	  for (int i = g-3; i < g; ++i)
	    if (a[i] != 0) { prnt = false; break; }
	  if (prnt)
	    {
	      for (int i = 0; i < g-3; ++i)
		{
		  cerr << "a[" << i << "]=" << a[i] << "\t";
		}
	      cerr << "\tN=" << N << endl;
	    }
	}

      // Copy the coefficients to the polynomial object.
      for (int m = 0; m < g; ++m) p[m+1] = a[m];

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
	  T lambda_sr = jlt::spectral_radius(p);
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
    } while (jlt::increment_vector(a,amin,amax));

  ostr << "}\n";
  cerr << "Found " << Np << "/" << N << " candidates.\n";
  cerr << "Found positive candidate root: ";
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
}


namespace jlt {

template<class S, int g, class T>
inline T findroot(const reciprocal_polynomial<S,g,T>& p,
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
T spectral_radius(const reciprocal_polynomial<S,g,T>& p)
{
  const int n = p.degree();
  mathmatrix<T> M(n,n);

  for (int i = 0; i < n-1; ++i) M(i,i+1) = 1;
  // Copy the coefficients of the polynomial to the matrix.
  for (int m = 0; m < n; ++m) M(n-1,m) = -(T)p[m]/p[n];

  return spectral_radius(M);
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

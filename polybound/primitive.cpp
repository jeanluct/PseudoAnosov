#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include <jlt/mathmatrix.hpp>
#include <jlt/polynomial.hpp>

template<class T>
bool increment_vector(std::vector<T>& a,
		      const std::vector<T>& amax);

#ifdef ROWSUM
template<class T>
bool increment_upper_triangle(std::vector<T>& a, const int n,
			      T& norm,
			      const T maxnorm,
			      std::vector<T>& colsum,
			      const T maxmincolsum);
#else
template<class T>
bool increment_upper_triangle(std::vector<T>& a, const int n,
			      T& norm, const T maxnorm);
#endif

template<class T>
T rowsum_bound(jlt::mathmatrix<T>& A);

bool skipstring(std::ifstream& strm, const std::string& s);

int main()
{
  using jlt::operator<<;
  using std::cout;
  using std::cerr;
  using std::endl;

  typedef jlt::mathmatrix<int>			Mat;
  typedef std::vector<int>			Vec;
  typedef Vec::iterator				Vecit;
  typedef Vec::const_iterator			Veccit;
  typedef jlt::polynomial<int>			Poly;
  typedef std::vector<Poly>			PVec;
  typedef PVec::iterator			PVecit;
  typedef PVec::const_iterator			PVeccit;
  typedef long long int				llint;

  std::ifstream indata;
  indata.open("polycoeffs_n=6.m");

  PVec pl;

#if 0
  const int n = 6;
  double lambdamax = 1.8311;

  while (skipstring(indata,"{"))
    {
      Poly p;

      for (int i = n; i >= 0; --i)
	{
	  indata >> p[i];
	  if (i > 0) skipstring(indata,",");
	}
      skipstring(indata,"}");
      pl.push_back(p);
    }
  cout << "Checking " << pl.size() << " polynomials:\n";
  for (PVeccit pi = pl.begin(); pi != pl.end(); ++pi) { cout << *pi << endl;} 

  Poly p = pl.front();
#else
  const int n = 4;
  double lambdamax = 1.73;

  Poly p;
  p[0] = 1;
  p[n] = 1;
  p[n/2] = -1;
  p[n/2+1] = -1;
  p[n/2-1] = -1;
#endif

  int tr = -p[n-1];

  // Create lower-triangle of matrices.
  std::list<Mat> Alow;
  // Vector of positions of the "1" entry.
  // 0 <= al[k] <= k+1.  al[k]=k+1 corresponds to no 1 at all.
  Vec al(n), almax(n);
  for (int k = 0; k < n; ++k) almax[k] = k+1;
  do
    {
      Mat A(n,n);
      int trA = 0;
      for (int k = 0; k < n; ++k)
	{
	  if (al[k] < k+1) A(k,al[k]) = 1;
	  if (al[k] == k) ++trA;
	}
      if (trA == tr) Alow.push_back(A);
    }
  while(increment_vector(al,almax));
  // Now figure out which matrices lead to a zero determinant.
  /* */

  cout << Alow.size() << endl;

  //
  // Upper-triangle of matrices.
  //
  std::list<Mat> Aup;
  int Nup = (n*(n-1))/2, aupmax = 2;
  Vec a(Nup), amax(Nup);
  for (int k = 0; k < Nup; ++k) { amax[k] = aupmax; }
  // Make index pair for eack k.
  Vec rowidx(Nup), colidx(Nup);
  int row = 0, col = 1;
  for (int k = 0; k < Nup; ++k)
    {
      rowidx[k] = row;
      colidx[k] = col++;
      if (col == n) { ++row; col = row+1; }
    }

  llint N = 0;
  llint rowsumexceeded = 0;	// Total times exceeded matrix row sum?

  int maxnorm = std::ceil(std::pow(lambdamax,(double)n)) + n - 1;
  cout << "maxnorm = " << maxnorm << endl;

  //  for (std::list<Mat>::const_iterator li = Alow.begin(); li != Alow.end(); ++li)
    std::list<Mat>::const_iterator li = Alow.begin();
    {
      // Calculate the norm for the lower matrix.
      int Alownorm = 0;
#ifdef ROWSUM
      Vec rowsum(n);
#endif
      for (int i = 0; i < n; ++i)
	{
	  for (int j = 0; j < n; ++j)
	    {
#ifdef ROWSUM
	      rowsum[i] += li->operator()(j,i);
#endif
	      Alownorm += li->operator()(i,j);
	    }
	}
      cout << "Norm of lower = " << Alownorm << endl;
#ifdef ROWSUM
      cout << "Row sums of lower = " << rowsum << endl;
#endif
      int Aupnorm = 0;

      do
	{
	  // Form matrix.
	  Mat A(*li);
	  for (int k = 0; k < Nup; ++k) A(rowidx[k],colidx[k]) = a[k];

	  ++N;
	  if (!(N % 200000))
	    { cerr << "a = " << a << endl; }

	  if (rowsum_bound(A) > lambdamax)
	    {
	      // cerr << "row/column sum bound exceeded.\n";
	      ++rowsumexceeded;
	    }

	  // Compute characteristic polynomial.
#if 0
	  Poly cpoly(A.charpoly());

	  if (cpoly == p)
	    {
	      cout << "Got it!\n";
	    }
#endif
	}
#ifdef ROWSUM
      while(increment_upper_triangle(a,n,Aupnorm,maxnorm-Alownorm,rowsum,1));
#else
      while(increment_upper_triangle(a,n,Aupnorm,maxnorm-Alownorm));
#endif
    }

    cout << rowsumexceeded << endl;
    cout << N << endl;
}


template<class T>
inline bool increment_vector(std::vector<T>& a,
			     const std::vector<T>& amax)
{
  //
  // Increment elements of the vector a (from the last element),
  // resetting each position a[m] to 0 if it passes amax[m].
  //
  const int n = a.size();
  bool incr = false;
  for (int m = n-1; m >= 0; --m)
    {
      // The sequence that coefficients cycle through is 0...amax,
      // rolling over from amax to 0 if necessary.
      ++a[m];
      if (a[m] > amax[m]) { a[m] = 0; continue; }
      incr = true;
      break;
    }
  return incr;
}


#ifdef ROWSUM
template<class T>
bool increment_upper_triangle(std::vector<T>& a, const int n,
			      T& norm,
			      const T maxnorm,
			      std::vector<T>& rowsum,
			      const T maxminrowsum)
#else
template<class T>
bool increment_upper_triangle(std::vector<T>& a, const int n,
			      T& norm, const T maxnorm)
#endif
{
  int m0 = 0;
  for (int k = 0; k < n; ++k)
    {
      for (int l = k+1; l < n; ++l)
	{
	  int m = m0 + (l-(k+1));
	  ++a[m];
	  ++norm;
#ifdef ROWSUM
	  ++rowsum[l];
#endif
	  if (norm > maxnorm)
	    {
	      norm -= a[m];
#ifdef ROWSUM
	      rowsum[l] -= a[m];
#endif
	      a[m] = 0; continue;
	    }
#ifdef ROWSUM
	  typename std::vector<T>::iterator
	    minrowsum = min_element(rowsum.begin(),rowsum.end());
	  using jlt::operator<<;
	  // std::cerr << "rowsum = " << rowsum << std::endl;
	  // std::cerr << "minrowsum = " << *minrowsum << std::endl;
	  if (*minrowsum > maxminrowsum)
	    {
	      // std::cout << "Row sum exceeded.\n";
	      norm -= a[m];
	      rowsum[l] -= a[m];
	      a[m] = 0;
	      continue;
	    }
#endif

#if 0
	  // Sanity check.
	  int checknorm = 0;
	  for (int i = 0; i < (n*(n-1))/2; ++i) checknorm += a[i];
	  if (norm != checknorm)
	    {
	      std::cerr << "Norms don't match...\n";
	      exit(-1);
	    }
#endif

	  return true;
	}
      m0 += (n-1-k);
    }
  return false;
}


template<class T>
inline T rowsum_bound(jlt::mathmatrix<T>& A)
{
  T rowsummin = -1;
  int n = A.dim();
  for (int j = 0; j < n; ++j)
    {
      int rowsum = 0;
      for (int i = 0; i < n; ++i)
	{
	  rowsum += A(i,j);
	}
      if (rowsum < rowsummin || rowsummin == -1) rowsummin = rowsum;
    }
  return rowsummin;
}


// Skip a string s in the stream strm.
//  * Check that the string matches;
//  * Completely ignore whitespace;
//  * ? is a wildcard.
bool skipstring(std::ifstream& strm, const std::string& s)
{
  char c;
  std::string s2(s);

  // Remove whitespace from s
  s2.erase(std::remove(s2.begin(),s2.end(),' '),s2.end());

  // Now skip the string, ignoring whitespace characters
  for (unsigned int i = 0; i < s2.length(); ++i)
    {
      // Read in, do not count whitespace
      do
	{
	  strm >> c;
	  if (strm.eof()) return false;
          if (c != ' ' && c != s2[i] && s2[i] != '?')
	    {
	      std::cerr << "Can't skip: string " << s;
	      std::cerr << " doesn't match.\n";
	      exit(1);
	    }
	}
      while (c == ' ');
    }

  return true;
}

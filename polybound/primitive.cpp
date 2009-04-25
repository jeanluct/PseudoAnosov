#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include <jlt/mathmatrix.hpp>
#include <jlt/polynomial.hpp>

template<class T>
bool increment_vector(std::vector<T>& a,
		      const std::vector<T>& amax);

template<class T>
bool increment_row_vector(std::vector<T>& a, const int n,
			  int& norm, const int maxnorm);

template<class T>
T rowcol_bound(jlt::mathmatrix<T>& A);

bool skipstring(std::ifstream& strm, const std::string& s);

int main()
{
  static const int debug = 0;

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
  llint normexceeded = 0;	// Total times exceeded matrix norm?
  llint rowcolsumexceeded = 0;	// Total times exceeded matrix column sum?
  llint colsumexceeded = 0;	// Total times exceeded matrix column sum?
  llint rowsumexceeded = 0;	// Total times exceeded matrix row sum?

  int maxnorm = std::ceil(std::pow(lambdamax,(double)n)) + n - 1;
  cout << "maxnorm = " << maxnorm << endl;

  //  for (std::list<Mat>::const_iterator li = Alow.begin(); li != Alow.end(); ++li)
    std::list<Mat>::const_iterator li = Alow.begin();
    {
      // Calculate the norm for the lower matrix.
      int Alownorm = 0;
      for (int i = 0; i < n; ++i)
	{
	  for (int j = 0; j < n; ++j)
	    {
	      Alownorm += li->operator()(i,j);
	    }
	}
      cout << "Norm of lower = " << Alownorm << endl;
      int Aupnorm = 0;

      do
	{
	  // Form matrix.
	  Mat A(*li);
	  for (int k = 0; k < Nup; ++k) A(rowidx[k],colidx[k]) = a[k];

	  ++N;
	  if (!(N % 100000))
	    { cerr << "a = " << a << endl; }

	  if (rowcol_bound(A) > lambdamax)
	    {
	      // cerr << "row/column sum bound exceeded.\n";
	      ++rowcolsumexceeded;
	    }
#if 0
	  int Mnorm = 0, colsummin = 0, rowsummin = 0;
	  for (int i = 0; i < n; ++i)
	    {
	      int colsum = 0, rowsum = 0;
	      for (int j = 0; j < n; ++j)
		{
		  colsum += A(j,i);
		  rowsum += A(i,j);
		}
	      Mnorm += colsum;
	      if (Mnorm > maxnorm)
		{
		  if (debug)
		    {
		      std::cerr << "Exceeded matrix norm ";
		    }
		  ++normexceeded;
		  continue;

		}
	      if (colsum < colsummin || colsummin == 0) colsummin = colsum;
	      if (rowsum < rowsummin || rowsummin == 0) rowsummin = rowsum;
	    }

	  // The minimum column sum is a lower bound on the spectral radius.
	  if (colsummin > lambdamax)
	    {
	      if (debug)
		{
		  std::cerr << "Exceeded column sum at pathlength ";
		}
	      ++colsumexceeded;
	      continue;
	    }
	  // The minimum row sum is also a lower bound on the spectral radius.
	  if (rowsummin > lambdamax)
	    {
	      if (debug)
		{
		  std::cerr << "Exceeded row sum at pathlength ";
		}
	      ++rowsumexceeded;
	      continue;
	    }
#endif
	  /*
	  cout << "norm = " << Mnorm;
	  cout << " colsum = " << colsummin;
	  cout << " rowsum = " << rowsummin << endl;
	  */

	  // Compute characteristic polynomial.
#if 0
	  Poly cpoly(A.charpoly());

	  if (cpoly == p)
	    {
	      cout << "Got it!\n";
	    }
#endif
	}
      while(increment_row_vector(a,n,Aupnorm,maxnorm-Alownorm));
    }

    cout << rowcolsumexceeded << endl;
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


template<class T>
inline bool increment_row_vector(std::vector<T>& a, const int n,
				 int& norm, const int maxnorm)
{
  int m0 = 0;
  for (int k = 0; k < n; ++k)
    {
      for (int l = k+1; l < n; ++l)
	{
	  int m = m0 + (l-(k+1));
	  ++a[m];
	  ++norm;
	  if (norm > maxnorm) { norm -= a[m]; a[m] = 0; continue; }

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
inline T rowcol_bound(jlt::mathmatrix<T>& A)
{
  T colsummin = 0, rowsummin = 0;
  int n = A.dim();
  for (int i = 0; i < n; ++i)
    {
      int colsum = 0, rowsum = 0;
      for (int j = 0; j < n; ++j)
	{
	  colsum += A(j,i);
	  rowsum += A(i,j);
	}
      if (colsum < colsummin || colsummin == 0) colsummin = colsum;
      if (rowsum < rowsummin || rowsummin == 0) rowsummin = rowsum;
    }
  return std::max(rowsummin,colsummin);
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

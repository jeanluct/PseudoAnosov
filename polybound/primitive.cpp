#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include <jlt/mathmatrix.hpp>
#include <jlt/polynomial.hpp>

template<class T>
bool increment_vector(std::vector<T>& a,
		      const std::vector<T>& amin,
		      const std::vector<T>& amax,
		      const std::vector<T>& astart,
		      const std::vector<T>& aend);

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

  const int n = 6;
  std::ifstream indata;
  indata.open("polycoeffs_n=6.m");

  PVec pl;

#if 1
  double lambdamax = 1.840;

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
  Vec al(n), almin(n), almax(n);
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
  while(increment_vector(al,almin,almax,almin,almax));
  // Now figure out which matrices lead to a zero determinant.
  /* */

  cout << Alow.size() << endl;

  //
  // Upper-triangle of matrices.
  //
  std::list<Mat> Aup;
  int Nup = (n*(n-1))/2, aupmax = 2;
  Vec a(Nup), amin(Nup), amax(Nup);
  for (int k = 0; k < Nup; ++k) amax[k] = aupmax;
  // Make index pair for eack k.
  Vec rowidx(Nup), colidx(Nup);
  int row = 0, col = 1;
  for (int k = 0; k < Nup; ++k)
    {
      rowidx[k] = row;
      colidx[k] = col++;
      if (col == n) { ++row; col = row+1; }
    }

  llint normexceeded = 0;	// Total times exceeded matrix norm?
  llint colsumexceeded = 0;	// Total times exceeded matrix column sum?
  llint rowsumexceeded = 0;	// Total times exceeded matrix row sum?

  int maxnorm = std::ceil(std::pow(lambdamax,(double)n)) + n - 1;
  cout << "maxnorm = " << maxnorm << endl;

  //  for (std::list<Mat>::const_iterator li = Alow.begin(); li != Alow.end(); ++li)
  std::list<Mat>::const_iterator li = Alow.begin();
    {
      do
	{
	  // Form matrix.
	  Mat A(*li);
	  for (int k = 0; k < Nup; ++k) A(rowidx[k],colidx[k]) = a[k];

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

	  /*
	  cout << "norm = " << Mnorm;
	  cout << " colsum = " << colsummin;
	  cout << " rowsum = " << rowsummin << endl;
	  */

	  // Compute characteristic polynomial.
	  Poly cpoly(A.charpoly());

	  if (cpoly == p)
	    {
	      cout << "Got it!\n";
	    }
	}
      while(increment_vector(a,amin,amax,amin,amax));
    }
}


template<class T>
inline bool increment_vector(std::vector<T>& a,
			     const std::vector<T>& amin,
			     const std::vector<T>& amax,
			     const std::vector<T>& astart,
			     const std::vector<T>& aend)
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
      // astart[m]...aend[m] for a[m], 0 <= m < n, rolling over from
      // amax[m] to amin[m] if necessary.
      ++a[m];
      // aend[m]=0 is the same as aend[m]=-1, since we skip 0.
      if (a[m] == aend[m]+1 || (a[m] == 0 && aend[m] == 0))
	{ a[m] = astart[m]; continue; }
      if (a[m] > amax[m]) { a[m] = amin[m]; }
      // Skip 0.
      if (a[m] == 0) ++a[m];
      incr = true;
      break;
    }

  return incr;
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

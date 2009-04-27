#include <iostream>
#include <fstream>
#include <vector>
#include <jlt/mathmatrix.hpp>
#include <jlt/polynomial.hpp>

template<class T>
bool increment_vector(std::vector<T>& a,
		      const std::vector<T>& amax);

template<class T>
bool increment_upper_triangle(std::vector<T>& a,
			      T& norm, const T maxnorm);

template<class T>
T mincolsum(jlt::mathmatrix<T>& A);

template<class T>
T matrix_norm(const jlt::mathmatrix<T>& A);

bool skipstring(std::ifstream& strm, const std::string& s);

int main()
{
  using jlt::operator<<;
  using std::cout;
  using std::cerr;
  using std::endl;
  using std::setw;

  typedef jlt::mathmatrix<int>			Mat;
  typedef std::vector<int>			Vec;
  typedef Vec::iterator				Vecit;
  typedef Vec::const_iterator			Veccit;
  typedef jlt::polynomial<int>			Poly;
  typedef std::vector<Poly>			PVec;
  typedef PVec::iterator			PVecit;
  typedef PVec::const_iterator			PVeccit;
  typedef std::vector<Mat>			MVec;
  typedef MVec::iterator			MVecit;
  typedef MVec::const_iterator			MVeccit;
  typedef long long int				llint;

  std::ifstream indata;
  indata.open("polycoeffs_n=6.m");

  PVec pl;

#if 1
  const int n = 6;

  // Read in Mathematica file of polynomials.
  skipstring(indata,
	     "(* Created by Wolfram Mathematica ??? : www.wolfram.com *)");
  skipstring(indata,"{");

  do
    {
      skipstring(indata,"{");

      Poly p;

      for (int i = n; i >= 0; --i)
	{
	  indata >> p[i];
	  if (i > 0) skipstring(indata,",");
	}
      skipstring(indata,"}");
      pl.push_back(p);
    }
  while (skipstring(indata,","));
  cerr << "Checking " << pl.size() << " polynomials:\n";
  for (PVeccit pi = pl.begin(); pi != pl.end(); ++pi) { cerr << *pi << endl;} 

#if 1
  // x^6 - x^4 - x^3 - x^2 + 1
  // 1710 matrices: all of type (0^28,1^8) or (0^27,1^9)
  Poly p = pl[0];
  double lambdamax = 1.40127;
#endif
#if 0
  // x^6 - x^5 + x^4 - 3 x^3 + x^2 - x + 1
  // No matrices found!  But it can clearly be realized with 7x7.
  Poly p = pl[1];
  double lambdamax = 1.46557;
#endif
#if 0
  // x^6 - x^5 - x^3 - x + 1
  // 1129 matrices: type (0^28,1^8), (0^27,1^9), (0^26,1^10), (0^26,1^9,2^1)
  Poly p = pl[2];
  double lambdamax = 1.50614;
#endif
#if 0
  // x^6 - x^5 - x^4 + x^3 - x^2 - x + 1
  Poly p = pl[3];
  double lambdamax = 1.55603;
#endif
#if 0
  // x^6 - 2 x^5 + x^3 - 2 x - x + 1
  Poly p = pl[17];
  double lambdamax = 1.83108;
#endif
#if 0
  // x^6 - x^4 - 4 x^3 - x^2 + 1
  Poly p = pl.back();
  double lambdamax = 1.83929;
#endif

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
  MVec Alow;
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

  cerr << Alow.size() << " lower-triangular forms\n";
  int todo = Alow.size();

  //
  // Upper-triangle of matrices: find valid patterns
  //
  int Nup = (n*(n-1))/2;
  Vec aup(Nup), aupmax(Nup);
  for (int k = 0; k < Nup; ++k) { aupmax[k] = 1; }
  // Make index pair for eack k.
  Vec rowidx(Nup), colidx(Nup);
  {
    int row = 0, col = 1;
    for (int k = 0; k < Nup; ++k)
      {
	rowidx[k] = row;
	colidx[k] = col++;
	if (col == n) { ++row; col = row+1; }
      }
  }

  llint N = 0;
  llint colsumexceeded = 0;	// Total times exceeded matrix row sum?
  llint maxnormexceeded = 0;
  llint validpatterns = 0;
  llint reduciblepatterns = 0;

  int maxnorm = (int)std::ceil(std::pow(lambdamax,(double)n)) + n - 1;
  cerr << "maxnorm = " << maxnorm << endl;

  cout << "{\n";
  bool thefirst = true;

  for (int li = 0; li < (int)todo; ++li)
    {
      // Calculate the norm for the lower matrix.
      int Alownorm = matrix_norm(Alow[li]);
      cerr << "Lower-triangular form " << setw(4) << li+1;
      cerr << " (norm " << Alownorm << "): ";

      // Vector of allowable patterns.
      std::vector<Vec> pattern;
      int n1min = -1, n1max = -1;

      do
	{
	  // Form matrix.
	  Mat A(Alow[li]);
	  for (int k = 0; k < Nup; ++k) A(rowidx[k],colidx[k]) = aup[k];

	  ++N;

	  if (mincolsum(A) <= lambdamax)
	    {
	      if (matrix_norm(A) <= maxnorm)
		{
		  if (!A.isReducible())
		    {
		      ++validpatterns;
		      pattern.push_back(aup);
		      int n1 = std::count(aup.begin(),aup.end(),1);
		      if (n1 < n1min || n1min == -1) n1min = n1;
		      if (n1 > n1max || n1max == -1) n1max = n1;
		    }
		  else
		    {
		      ++reduciblepatterns;
		    }
		}
	      else
		{
		  ++maxnormexceeded;
		}
	    }
	  else
	    {
	      ++colsumexceeded;
	    }
	}
      while(increment_vector(aup,aupmax));

      cerr << setw(5) << pattern.size() << " valid patterns";

      // Skip lower-triangular forms with no allowable patterns.
      if (pattern.empty()) { cerr << endl; continue; }

      cerr << " (size ";
      cerr << setw(2) << n1min+Alownorm << " to ";
      cerr << setw(2) << n1max+Alownorm << ")\n";

      for (int pa = 0; pa < (int)pattern.size(); ++pa)
	{
	  int Npat = std::count(pattern[pa].begin(),pattern[pa].end(),1);
	  // Make index pair for eack k.
	  Vec patrowidx(Npat), patcolidx(Npat);
	  {
	    int row = 0, col = 1, m = 0;
	    for (int k = 0; k < (n*(n-1))/2; ++k)
	      {
		if (pattern[pa][k] == 1)
		  {
		    patrowidx[m] = row;
		    patcolidx[m] = col;
		    ++m;
		  }
		if (++col == n) { ++row; col = row+1; }
	      }
	  }

	  // Start with the pattern equal to all ones.
	  Vec apat(Npat,1);
	  int Aupnorm = apat.size();

	  // Now loop over vales of the entries of the pattern
	  do
	    {
	      // Form matrix.
	      Mat A(Alow[li]);
	      for (int k = 0; k < Npat; ++k)
		A(patrowidx[k],patcolidx[k]) = apat[k];

	      ++N;
	      /*
		if (!(N % 1000000))
		{ cerr << "a = " << apat << endl; }
	      */

	      // Compute characteristic polynomial.
	      Poly cpoly(A.charpoly());

	      if (cpoly == p)
		{
		  cerr << "Got it!\n";
		  if (!thefirst) cout << "," << endl; else thefirst = false;
		  A.printMathematicaForm(cout);
		}
	    }
	  while(increment_upper_triangle(apat,Aupnorm,maxnorm-Alownorm));
	}
    }


  cout << "\n}\n";

  cerr << validpatterns << endl;
  cerr << reduciblepatterns << endl;
  cerr << colsumexceeded << endl;
  cerr << maxnormexceeded << endl;
  cerr << N << endl;
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
bool increment_upper_triangle(std::vector<T>& a,
			      T& norm, const T maxnorm)
{
  const int n = a.size();
  for (int m = n-1; m >= 0; --m)
    {
      ++a[m];
      ++norm;
      if (norm > maxnorm) { norm -= (a[m]-1); a[m] = 1; continue; }

#if 0
      // Sanity check.
      int checknorm = 0;
      for (int i = 0; i < n; ++i) checknorm += a[i];
      if (norm != checknorm)
	{
	  std::cerr << "Norms don't match...\n";
	  exit(-1);
	}
#endif

      return true;
    }
  return false;
}


template<class T>
inline T mincolsum(jlt::mathmatrix<T>& A)
{
  T colsummin = -1;
  int n = A.dim();
  for (int j = 0; j < n; ++j)
    {
      int colsum = 0;
      for (int i = 0; i < n; ++i)
	{
	  colsum += A(i,j);
	}
      if (colsum < colsummin || colsummin == -1) colsummin = colsum;
    }
  return colsummin;
}


template<class T>
inline T matrix_norm(const jlt::mathmatrix<T>& A)
{
  T norm = 0;
  int n = A.dim();
  for (int j = 0; j < n; ++j)
    {
      for (int i = 0; i < n; ++i)
	{
	  norm += A(i,j);
	}
    }
  return norm;
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
#if 0
	      std::cerr << "Can't skip: string " << s;
	      std::cerr << " doesn't match.\n";
#endif
	      return false;
	    }
	}
      while (c == ' ');
    }

  return true;
}

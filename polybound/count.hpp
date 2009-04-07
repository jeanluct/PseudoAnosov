#ifndef JLT_COUNT_HPP
#define JLT_COUNT_HPP

/* Make this into a class to automatically increase a counter */

//
// count.hpp
//

#include <vector>

class vector_iterator<T = long long int>
{
private:
  typedef typename std::vector<T>			Vec;
  typedef typename std::vector<T>::size_type		size_type;

  Vec a, amin, amax;

public:
  vector_iterator(Vec& _a, Vec& _amin,  Vec& _amax)
    : a(_a), amin(_amin), amax(_amax) { set_to(_a); }

  vector_iterator(Vec& _amin,
		  Vec& _amax) : a(_amin), amin(_amin), amax(_amax) {}

  vector_iterator(Vec& _amax)
    : a(_amax.size()), amin(_amax.size()), amax(_amax) {}

  vector_iterator(size_type _n, T _amin, T _amax)
    : a(_n,_amin), amin(_n,_amin), amax(_n,_amax) {}

  void set_to(const Vec& _a)
  {
    // Check to make sure that a has the right size and is in the
    // right range.
    /* */

    a = _a;
  }

  Vec& operator++()
  {
  }
}

} // namespace jlt

#endif // JLT_COUNT_HPP

/* Copyright (C) 2005-2011 M. T. Homer Reid
 *
 * This file is part of SCUFF-EM.
 *
 * SCUFF-EM is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * SCUFF-EM is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

/* libVec.h, contributed by johannes feist */
/*
  This library consists only of the include file libVec.h, which
  defines a template class Vec describing a vector in 3-dimensional
  Euclidean space (and generalizations depending on the template
  type).  Through operator overloading, these vectors can be
  manipulated very intuitively.  All functions are inline, so a
  (smart) compiler should be able to convert most expressions to
  direct evaluation without creating too many temporaries or actually
  calling functions.
  
  In principle, this class just provides syntactic sugar for handling
  T [3] arrays (where T is the type in the template argument), which
  is all the data in the class. In order to use the class with other
  codes expecting raw T[3] arrays, we provide a constructor of the
  class from T*, so a T* can be passed to a function expecting a
  dVec. In addition, we provide an operator T*, which can be used to
  pass an object of type Vec<T> to a function expecting a T* array.
  It is the user's responsibility to ensure that any usage as T* only
  refers to T[3] arrays!

  If T is a complex type (i.e. T=std::complex<U>), 
  some template metaprogramming tricks are used to allow specialization
  of the functions.

  Two specializations are provided as typedefs: dVec for Vec<double>
  and zVec for Vec<std::complex<double> >. 
 */

#ifndef LIBVEC_H
#define LIBVEC_H

#include <cmath>   // for sqrt and acos
#include <cstring> // for memset
#include <iosfwd>
#include <complex>

// some helper templates that allow for partial specialization of the routines below
// depending on whether the type is a normal type or an extension of <complex>
template <bool val> struct bool_type { typedef bool_type<val> type; };
typedef bool_type<false> false_type;
typedef bool_type<true>   true_type;

template <typename T> struct is_complex                   : false_type { typedef T norm_type; };
template <typename U> struct is_complex<std::complex<U> > :  true_type { typedef U norm_type; };

template <typename T> class Vec {
private:
  T v[3];
public:
  inline Vec() { zero(); }
  inline Vec(T a0, T a1, T a2) { v[0]=a0; v[1]=a1; v[2]=a2; }
  inline Vec(const T* const a) { v[0]=a[0]; v[1]=a[1]; v[2]=a[2]; }
  template <typename U> inline Vec(Vec<U> uv) { v[0]=uv[0]; v[1]=uv[1]; v[2]=uv[2]; }

  inline operator T*() { return v;} // this allows to use a Vec anywhere a T* is needed
  T      & operator[] (int ii)       { return v[ii]; }
  T const& operator[] (int ii) const { return v[ii]; }

  inline static Vec xHat() { return Vec(1,0,0); }
  inline static Vec yHat() { return Vec(0,1,0); }
  inline static Vec zHat() { return Vec(0,0,1); }

  inline Vec& operator*=(T alpha)       { v[0]*=alpha; v[1]*=alpha; v[2]*=alpha; return *this;}
  inline Vec& operator/=(T alpha)       { v[0]/=alpha; v[1]/=alpha; v[2]/=alpha; return *this;}
  inline Vec& operator+=(Vec const& rhs) { v[0]+=rhs[0]; v[1]+=rhs[1]; v[2]+=rhs[2]; return *this;}
  inline Vec& operator-=(Vec const& rhs) { v[0]-=rhs[0]; v[1]-=rhs[1]; v[2]-=rhs[2]; return *this;}
  inline Vec const operator* (T alpha) const  { return Vec(alpha*v[0], alpha*v[1], alpha*v[2]); }
  friend inline Vec const operator* (T alpha, Vec const& rhs) { return rhs*alpha; }
  inline Vec const operator/ (T alpha) const  { T inva = (T)1/alpha; return *this * inva; }
  inline Vec const operator+ (Vec const& rhs) const { return Vec(v[0]+rhs[0],v[1]+rhs[1],v[2]+rhs[2]);}
  inline Vec const operator- (Vec const& rhs) const { return Vec(v[0]-rhs[0],v[1]-rhs[1],v[2]-rhs[2]);}
  // unary minus
  inline Vec const operator- () const { return Vec(-v[0],-v[1],-v[2]);}
  friend inline std::ostream& operator<<(std::ostream& os, Vec const& vv) { return os << vv[0] << " " << vv[1] << " " << vv[2]; }
  friend inline std::istream& operator>>(std::istream& is, Vec&       vv) { return is >> vv[0] >> vv[1] >> vv[2]; }

  inline void zero() { memset(v,0,3*sizeof(T)); } // caution!! this only works with "primitive" types.

  // functions where we need to distinguish between complex / not complex
  typedef typename is_complex<T>::type      is_complex_type;
  typedef typename is_complex<T>::norm_type norm_type;

  inline norm_type norm2() const { return norm2_( is_complex_type() ); }
  inline norm_type norm()  const { return sqrt(norm2()); }
  // explicitly use double here because SWIG gets confused otherwise.
  // this should not hurt because numerical types can all be converted from/to double
  inline norm_type normalize(double R = 1.) { norm_type d=norm(); *this *= R/d; return d; };
  friend inline norm_type norm2(Vec const& dv) { return dv.norm2(); }
  friend inline norm_type norm (Vec const& dv) { return dv.norm (); }

  friend inline norm_type  dist(Vec const& lhs, Vec const& rhs) { return (lhs-rhs).norm(); }
  friend inline T           dot(Vec const& lhs, Vec const& rhs) { return lhs[0]*rhs[0] + lhs[1]*rhs[1] + lhs[2]*rhs[2]; }
  // dot product with complex conjugation - only makes sense for complex Vecs!
  friend inline T          dotC(Vec const& lhs, Vec const& rhs) { return dot(conj(lhs),rhs); }

  // angles only make sense for non-complex Vecs
  friend inline T          cosangle(Vec const& lhs, Vec const& rhs) { return dot(lhs,rhs)/(lhs.norm() * rhs.norm()); }
  friend inline T             angle(Vec const& lhs, Vec const& rhs) { return acos(cosangle(lhs,rhs)); }
  
  friend inline Vec const     cross(Vec const& lhs, Vec const& rhs) { return Vec(lhs[1]*rhs[2] - lhs[2]*rhs[1],
										 lhs[2]*rhs[0] - lhs[0]*rhs[2],
										 lhs[0]*rhs[1] - lhs[1]*rhs[0]); }

  // these only make sense for complex Vecs
  friend inline Vec const      conj(Vec const& dv) { return Vec(conj(dv[0]),conj(dv[1]),conj(dv[2])); }

  inline Vec<norm_type> real() const { return Vec<norm_type>(v[0].real(),v[1].real(),v[2].real()); }
  inline Vec<norm_type> imag() const { return Vec<norm_type>(v[0].imag(),v[1].imag(),v[2].imag()); }
  
  inline std::ostream& printreim(std::ostream& os) { return os << v[0].real() << " " << v[0].imag() << " " 
							       << v[1].real() << " " << v[1].imag() << " " 
							       << v[2].real() << " " << v[2].imag(); }


private:
  inline norm_type norm2_( false_type ) const { return dot(*this,*this); }
  inline norm_type norm2_(  true_type ) const { return std::norm(v[0])+std::norm(v[1])+std::norm(v[2]); }
};

typedef Vec<double> dVec;
typedef Vec<std::complex<double> > zVec;

// cross product with complex conjugation of first argument
template<typename T> inline Vec<T> const crossC(Vec<T> const& lhs, Vec<T> const& rhs) { return cross(conj(lhs),rhs); }

template<typename T> inline Vec<T> const GetOrthogonalNormVec(Vec<T> const& vin) {
  Vec<T> vout = cross(Vec<T>(1,0,0),vin);
  // if test vector was (nearly) collinear, we take an orthogonal one
  if (norm(vout) < 1.e-5 * norm(vin))
    vout = cross(Vec<T>(0,1,0),vin);
  vout.normalize();
  return vout;
}

/* this calculates the local cylindrical coordinates r and z of a point. input:
   Xrel: cartesian vector of the point relative to the origin of the system
   nHat: unit vector along the z-axis of the coordinate system

   we use the formulas from 
   http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html, 
   where x1 = (0,0,0) and x2 = nHat (|nHat|==1!), and
   formula (13) from http://mathworld.wolfram.com/Point-PlaneDistance.html,
   where Xi = (0,0,0) (any point on the plane) and n = nHat
*/
template<typename T> inline void GetLocalCylinderCoordinates(Vec<T> const& Xrel, Vec<T> const& nHat, T& r, T& z) {
  r = norm(cross(nHat, Xrel));
  z = dot(nHat, Xrel);
}

// specialization to forbid this operation for complex vectors, where it does not make sense.
template <typename U> inline void GetLocalCylinderCoordinates(Vec<std::complex<U> > const& Xrel, Vec<std::complex<U> > const& nHat, std::complex<U>& r, std::complex<U>& z) {
  // do an operation that is not allowed (call a complex number like a function)
  r(5);
}


#endif // #ifndef LIBVEC_H

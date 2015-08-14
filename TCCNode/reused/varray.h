//
//  varray.h
//  VTest2
//
//  Created by Denis Kovacs on 6/28/12.
//  Copyright 2012 __MyCompanyName__. All rights reserved.
//

#ifndef VARRAY_H
#define VARRAY_H

#include <valarray>
#define NDEBUG
#include <assert.h>
#include <vector>
#include <set>

namespace std {
    
    enum V_All_Type { _ };
    enum V_ColRow_Type { v_col, v_row };
    
    template<class T> struct v_array;
    
    
    template<class T> struct vindirect_array;
    template<class T> struct vmask_array;
    template<class T> struct vslice_array;
    
    
    template<class T>
    struct v_array
    {
        size_t dim[2];
        std::valarray<T> v;
        
        v_array(size_t d0=0, size_t d1=0, T c=T()): v(c, d0*d1) { this->dim[0]=d0; this->dim[1]=d1; }
        v_array(size_t d0, size_t d1, const std::valarray<T> &va):v(va) { dim[0] = d0; dim[1] = d1; }
        v_array(T *va, size_t n):v(va, n) { dim[0] = n; dim[1] = 1; }
        
        v_array(const vslice_array<T> &vs)   :v((*vs.src).v[vs.sl]) { dim[0] = (vs.cr==v_col) ? vs.src->dim[0]: 1 ; dim[1] = (vs.cr==v_col) ? 1 : vs.src->dim[1]; }
        v_array(const vmask_array<T> &vm)    :v((*vm.src).v[vm.mask.v]) { dim[0] = v.size(); dim[1] = 1; }
        v_array(const vindirect_array<T> &vm):v(T(), vm.ind.dim[0]*vm.ind.dim[1]) { dim[0] = vm.ind.dim[0]; dim[1]=vm.ind.dim[1]; v=(*vm.src).v[vm.ind.v]; }
        v_array(const std::valarray<T> &va)  :v(va) { dim[0] = v.size(); dim[1] = 1; }
        v_array(const v_array<T> &va)        :v(va.v) { 
            dim[0] = va.dim[0]; 
            dim[1] = va.dim[1]; 
        }
        
        T operator[](size_t idx) const { assert(idx<v.size()); return v[idx]; }
        T operator()(size_t idx) const { assert(idx<v.size()); return v[idx]; }
        T operator()(size_t r, size_t c) const { assert(r<dim[0] && c<dim[1]); return v[c*this->dim[0] + r]; }

        T &operator[](size_t idx)  { assert(idx<v.size()); return v[idx]; }
        T &operator()(size_t idx)  { assert(idx<v.size()); return v[idx]; }
        T &operator()(size_t r, size_t c)  { assert(r<dim[0] && c<dim[1]); return v[c*this->dim[0] + r]; }
        
        std::vindirect_array<T> operator ()(const v_array<size_t> &c);
        std::vmask_array<T>     operator ()(const v_array<bool> &c);
        
        std::vindirect_array<T> operator ()(V_All_Type d,  v_array<size_t> c);
        std::vmask_array<T>     operator ()(V_All_Type d,  v_array<bool> c);
        std::vindirect_array<T> operator()( v_array<size_t> r, V_All_Type d);
        std::vmask_array<T>     operator()( v_array<bool> r, V_All_Type d);
        
        std::vslice_array<T>    operator()(V_All_Type d, size_t c);
        std::vslice_array<T>    operator()(size_t r, V_All_Type d);
        
        void reshape(size_t d1, size_t d2);
        void resize(size_t d1, size_t d2, T c=T());
        
        void operator =(T val) { v = val; }
        void operator =( std::v_array<T> o) { dim[0] = o.dim[0]; dim[1] = o.dim[1]; v.resize(o.v.size()); v = o.v; }

        void operator +=(T val) { v += val; }
        void operator +=( std::v_array<T> o) { assert(dim[0]==o.dim[0] &&dim[1]==o.dim[1]); v += o.v; }

        void operator -=(T val) { v -= val; }
        void operator -=( std::v_array<T> o) { assert(dim[0]==o.dim[0] &&dim[1]==o.dim[1]); v -= o.v; }
        
        void operator *=(T val) { v *= val; }
        void operator *=( std::v_array<T> o) { assert(dim[0]==o.dim[0] &&dim[1]==o.dim[1]); v *= o.v; }
        
        void operator /=(T val) { v /= val; }
        void operator /=( std::v_array<T> o) { assert(dim[0]==o.dim[0] &&dim[1]==o.dim[1]); v /= o.v; }
        
        size_t nR() const { return dim[0]; }
        size_t nC() const { return dim[1]; }
    };
    

#define VINDIRECT_ARRAY_ASSIGNMENT(_Op) \
 void operator _Op(T val)                 { (*src).v[ind.v] _Op val; } \
 void operator _Op( v_array<T> o)         { (*src).v[ind.v] _Op o.v; }\
 void operator _Op( vindirect_array<T> o) { (*src).v[ind.v] _Op (*o.src).v[o.ind.v]; }\

    
    template<class T>
    struct vindirect_array
    {
        v_array<T> *src;
        v_array<size_t> ind; // indirection
        
        vindirect_array(v_array<T> *vs, const v_array<size_t> &vi): src(vs), ind(vi) { }
        
        VINDIRECT_ARRAY_ASSIGNMENT(=)
        VINDIRECT_ARRAY_ASSIGNMENT(+=)
        VINDIRECT_ARRAY_ASSIGNMENT(-=)
        VINDIRECT_ARRAY_ASSIGNMENT(*=)
        VINDIRECT_ARRAY_ASSIGNMENT(/=)
    };

#define VMASK_ARRAY_ASSIGNMENT(_Op) \
 void operator _Op(T val) { (*src).v[mask.v] _Op val; } \
 void operator _Op( v_array<T> o)     { (*src).v[mask.v] _Op o.v; } \
 void operator _Op( vmask_array<T> o) { (*src).v[mask.v] _Op (*o.src).v[mask.v]; }
    
    template<class T>
    struct vmask_array
    {
        v_array<T> *src;
        v_array<bool> mask;
        
        vmask_array(v_array<T> *vs, const v_array<bool> &m): src(vs), mask(m) { }
        
        VMASK_ARRAY_ASSIGNMENT(=)
        VMASK_ARRAY_ASSIGNMENT(+=)
        VMASK_ARRAY_ASSIGNMENT(-=)
        VMASK_ARRAY_ASSIGNMENT(*=)
        VMASK_ARRAY_ASSIGNMENT(/=)
    };
    
#define VSLICE_ARRAY_ASSIGNMENT(_Op) \
 void operator _Op(T val)              { (*src).v[sl] _Op val; } \
 void operator _Op( v_array<T> o)      { (*src).v[sl] _Op o.v; } \
 void operator _Op( vslice_array<T> o) { (*src).v[sl] _Op (*o.src).v[o.sl]; }

    
    template<class T>
    struct vslice_array
    {
        v_array<T> *src;
        std::slice sl;
        V_ColRow_Type cr;
        
        vslice_array(v_array<T> *vs, std::slice slice, V_ColRow_Type t): src(vs), sl(slice), cr(t) {  }

        VSLICE_ARRAY_ASSIGNMENT(=)
        VSLICE_ARRAY_ASSIGNMENT(+=)
        VSLICE_ARRAY_ASSIGNMENT(-=)
        VSLICE_ARRAY_ASSIGNMENT(*=)
        VSLICE_ARRAY_ASSIGNMENT(/=)
    };
    
    
    template<class T> 
    ostream& operator <<(ostream& out,  v_array<T> v);
    
    template<class T> 
    size_t size(const std::v_array<T> &a, V_ColRow_Type cr);
    
    template<class T> T max(v_array<T> v);
    template<class T> T min(v_array<T> v);
    
    v_array<size_t> v_range(size_t k1, size_t k2);
    
    template<class T> 
    v_array<T> v_cat(size_t dim, v_array<T> v1, v_array<T> v2);
    
    
    // implementation 
    template<class T>
    std::vindirect_array<T> v_array<T>::operator ()(const v_array<size_t> &c)
    {
        return vindirect_array<T>(this, c);
    }
    
    
    
    template<class T>
    inline void v_array<T>::resize(size_t d1, size_t d2, T c)
    {
        dim[0]=d1; dim[1]=d2;
        size_t n = d1*d2;
        
        if (this->v.size()!=n)
        {
            // vector style behavior
            std::valarray<T> temp(this->v);
            this->v.resize(n, c);
            std::copy( &temp[0], &temp[n < temp.size() ? n : temp.size()], &(this->v[0]));
        }
    }
    
    template<class T>
    inline void v_array<T>::reshape(size_t d1, size_t d2) { 
        assert(d1*d2 == v.size()); dim[0]=d1; dim[1]=d2; 
    }
    
    // array mask
    template<class T>
    std::vmask_array<T> v_array<T>::operator ()(const v_array<bool> &c)
    {
        size_t nnz = 0;
        for (size_t k=0; k<c.v.size(); k++) if (c.v[k]) nnz++; 
        
        return vmask_array<T>(this, c);
    }
    // select one column)
    template<class T> 
    inline std::vslice_array<T> v_array<T>::operator()(V_All_Type d, size_t c)
    { 
        assert((c <= this->dim[1]));
        return vslice_array<T>(this, std::slice(c*this->dim[0], this->dim[0], 1), v_col);
    }
    
    // select one row
    template<class T> 
    inline std::vslice_array<T> v_array<T>::operator()(size_t r, V_All_Type d)
    { 
        assert((r <= this->dim[0]));
        return vslice_array<T>(this, std::slice(r, this->dim[1], this->dim[0]), v_row);
    }
    // select columns (indirection)
    template<class T> 
    inline std::vindirect_array<T> v_array<T>::operator()(V_All_Type d,  v_array<size_t> c) 
    { 
        assert((c.v.max() <= this->dim[1]));
        std::v_array<size_t> indirect(this->dim[0], c.v.size());
        for (size_t k=0; k<c.v.size(); k++)
        {
            for (size_t kR=0; kR<this->dim[0]; kR++)
            {
                indirect(kR, k) = c(k)*this->dim[0] + kR;
            }
        }
        return vindirect_array<T>(this, indirect);
    }
    // select columns (mask)
    template<class T> 
    inline std::vmask_array<T> v_array<T>::operator()(V_All_Type d,  v_array<bool> c)
    {
        assert(c.v.size() == this->dim[1]);
        std::v_array<bool> mask(this->dim[0], this->dim[1], false);
        size_t nC = 0;
        for (size_t k=0; k<this->dim[1]; k++)
        {
            if (c(k)) { mask(_ , k) = true; nC++; }
        }
        
        return vmask_array<T>(this, mask);
    }
    // select rows (indirection)
    template<class T> 
    inline std::vindirect_array<T> v_array<T>::operator()( v_array<size_t> r, V_All_Type d)
    {
        assert( (r.v.max() <= this->dim[0]) );
        std::v_array<size_t> indirect(r.v.size(), this->dim[1]);
        for (size_t k=0; k<r.v.size(); k++)
        {
            for (size_t kC=0; kC<this->dim[1]; kC++)
            {
                indirect(k, kC) = kC*this->dim[0] + r(k);
            }
        }
        
        return vindirect_array<T>(this, indirect);
    }
    // select rows (mask)
    template<class T> 
    inline std::vmask_array<T> v_array<T>::operator()( v_array<bool> r, V_All_Type d)
    {
        assert(r.v.size() == this->dim[0]);
        std::v_array<bool> mask(this->dim[0], this->dim[1], false);
        size_t nR = 0;
        for (size_t k=0; k<this->dim[0]; k++)
        {
            if (r(k)) { mask(k, _) = true; nR++; }
        }
        return vmask_array<T>(this, mask);
    }

    
    typedef v_array<double> v_double;
    typedef v_array<size_t> v_size_t;
    typedef v_array<bool> v_bool;
    typedef v_array<uint64_t> v_uint64;
    

#define V_ARRAY_CONSTANTS(_Name, _Type) \
inline v_array<_Type> cv_##_Name(_Type v0) { v_array<_Type> v(1, 1); v(0) = v0; return v; } \
inline v_array<_Type> cv_##_Name(_Type v0, _Type v1) { v_array<_Type> v(1, 2); v(0) = v0; v(1) = v1; return v; } \
inline v_array<_Type> cv_##_Name(_Type v0, _Type v1, _Type v2) { v_array<_Type> v(1, 3); v(0) = v0; v(1) = v1; v(2) = v2; return v; } \
inline v_array<_Type> cv_##_Name(_Type v0, _Type v1, _Type v2, _Type v3) { v_array<_Type> v(1, 4); v(0) = v0; v(1) = v1; v(2) = v2; v(3) = v3; return v; } \
inline v_array<_Type> cv_##_Name(_Type v0, _Type v1, _Type v2, _Type v3, _Type v4) { v_array<_Type> v(1, 5); v(0) = v0; v(1) = v1; v(2) = v2; v(3) = v3; v(4) = v4; return v; } 
    
    V_ARRAY_CONSTANTS(bool, bool);
    V_ARRAY_CONSTANTS(double, double);
    V_ARRAY_CONSTANTS(size_t, size_t);
    V_ARRAY_CONSTANTS(int, int);
    V_ARRAY_CONSTANTS(uint64_t, uint64_t);
    
    
#define _V_ARRAY_BINARY_ASSERT assert(v1.dim[0]==v2.dim[0] && v1.dim[1]==v2.dim[1])
#define _V_ARRAY_DEFINE_BINARY_OPERATOR(_Op, _Type, _RetType) \
inline v_array<_RetType> operator _Op(const v_array<_Type> &v1, const v_array<_Type> &v2) { _V_ARRAY_BINARY_ASSERT; return v_array<_RetType>(v1.dim[0], v1.dim[1], v1.v _Op v2.v); }   \
\
inline v_array<_RetType> operator _Op(const _Type s1, const v_array<_Type> &v1) { \
std::valarray<_RetType> bl(v1.v.size()); \
for (size_t k=0; k<v1.v.size(); k++) bl[k]=s1 _Op v1.v[k];\
return v_array<_RetType>(v1.dim[0], v1.dim[1], bl); \
}  \
\
inline v_array<_RetType> operator _Op(const v_array<_Type> &v1, const _Type s1) { \
std::valarray<_RetType> bl(v1.v.size()); \
for (size_t k=0; k<v1.v.size(); k++) bl[k]=v1.v[k] _Op s1;\
return v_array<_RetType>(v1.dim[0], v1.dim[1], bl); \
} 
    
    _V_ARRAY_DEFINE_BINARY_OPERATOR(+, bool, bool)
    _V_ARRAY_DEFINE_BINARY_OPERATOR(-, bool, bool)
    _V_ARRAY_DEFINE_BINARY_OPERATOR(*, bool, bool)
    _V_ARRAY_DEFINE_BINARY_OPERATOR(/, bool, bool)
    _V_ARRAY_DEFINE_BINARY_OPERATOR(&, bool, bool)
    _V_ARRAY_DEFINE_BINARY_OPERATOR(|, bool, bool)
    _V_ARRAY_DEFINE_BINARY_OPERATOR(%, bool, bool)
    
    _V_ARRAY_DEFINE_BINARY_OPERATOR(+, size_t, size_t)
    _V_ARRAY_DEFINE_BINARY_OPERATOR(-, size_t, size_t)
    _V_ARRAY_DEFINE_BINARY_OPERATOR(*, size_t, size_t)
    _V_ARRAY_DEFINE_BINARY_OPERATOR(/, size_t, size_t)
    _V_ARRAY_DEFINE_BINARY_OPERATOR(&, size_t, size_t)
    _V_ARRAY_DEFINE_BINARY_OPERATOR(|, size_t, size_t)
    _V_ARRAY_DEFINE_BINARY_OPERATOR(%, size_t, size_t)
    _V_ARRAY_DEFINE_BINARY_OPERATOR(==, size_t, bool)
    _V_ARRAY_DEFINE_BINARY_OPERATOR(!=, size_t, bool)
    _V_ARRAY_DEFINE_BINARY_OPERATOR(<, size_t, bool)
    _V_ARRAY_DEFINE_BINARY_OPERATOR(>, size_t, bool)
    _V_ARRAY_DEFINE_BINARY_OPERATOR(<=, size_t, bool)
    _V_ARRAY_DEFINE_BINARY_OPERATOR(>=, size_t, bool)
    
    _V_ARRAY_DEFINE_BINARY_OPERATOR(+, uint64_t, uint64_t)
    _V_ARRAY_DEFINE_BINARY_OPERATOR(-, uint64_t, uint64_t)
    _V_ARRAY_DEFINE_BINARY_OPERATOR(*, uint64_t, uint64_t)
    _V_ARRAY_DEFINE_BINARY_OPERATOR(/, uint64_t, uint64_t)
    _V_ARRAY_DEFINE_BINARY_OPERATOR(&, uint64_t, uint64_t)
    _V_ARRAY_DEFINE_BINARY_OPERATOR(|, uint64_t, uint64_t)
    _V_ARRAY_DEFINE_BINARY_OPERATOR(%, uint64_t, uint64_t)
    _V_ARRAY_DEFINE_BINARY_OPERATOR(==, uint64_t, bool)
    _V_ARRAY_DEFINE_BINARY_OPERATOR(!=, uint64_t, bool)
    _V_ARRAY_DEFINE_BINARY_OPERATOR(<, uint64_t, bool)
    _V_ARRAY_DEFINE_BINARY_OPERATOR(>, uint64_t, bool)
    _V_ARRAY_DEFINE_BINARY_OPERATOR(<=, uint64_t, bool)
    _V_ARRAY_DEFINE_BINARY_OPERATOR(>=, uint64_t, bool)

    _V_ARRAY_DEFINE_BINARY_OPERATOR(+, double, double)
    _V_ARRAY_DEFINE_BINARY_OPERATOR(-, double, double)
    _V_ARRAY_DEFINE_BINARY_OPERATOR(*, double, double)
    _V_ARRAY_DEFINE_BINARY_OPERATOR(/, double, double)
    _V_ARRAY_DEFINE_BINARY_OPERATOR(==, double, bool)
    _V_ARRAY_DEFINE_BINARY_OPERATOR(!=, double, bool)
    _V_ARRAY_DEFINE_BINARY_OPERATOR(<, double, bool)
    _V_ARRAY_DEFINE_BINARY_OPERATOR(>, double, bool)
    _V_ARRAY_DEFINE_BINARY_OPERATOR(<=, double, bool)
    _V_ARRAY_DEFINE_BINARY_OPERATOR(>=, double, bool)
    
    
    inline v_array<bool> operator !(v_array<bool> v1) { return v_array<bool>(v1.dim[0], v1.dim[1], !v1.v); }
    
    template<class T> inline T max(v_array<T> v) { return v.v.max(); }
    template<class T> inline T min(v_array<T> v) { return v.v.min(); }
    
    
    template<class T> inline size_t size(const std::v_array<T> &a, V_ColRow_Type cr) { return a.dim[cr==v_row?0:1]; }
    
    inline v_size_t v_range(size_t k1, size_t k2)
    {
        v_size_t va(1, k2-k1+1);
        for (size_t k=k1; k<=k2; k++) va(k-k1) = k;
        return va;
    }
    
    template<class T> 
    inline size_t length(std::v_array<T> &a)
    {
        return a.v.size();
    }

    template<class T> 
    inline size_t length(std::vector<T> &a)
    {
        return a.size();
    }
    
    template<class T> inline void v_transpose_1D(v_array<T> &v1)
    {
        v1.reshape(v1.dim[1], v1.dim[0]);
    }
    
    
    template<class T> 
    inline v_array<T> v_cat(V_ColRow_Type cr, v_array<T> v1, v_array<T> v2)
    {
        v_array<T> va;
        if (cr==v_row) // concatenation along rows 
        {
            assert(v1.dim[1] == v2.dim[1]);
            va.resize(v1.dim[0] + v2.dim[0], v1.dim[1]);
            
            size_t n = 0;
            for (size_t k1 = 0; k1<v1.dim[1]; k1++)
            {
                for (size_t k0 = 0; k0<v1.dim[0]; k0++) { va(n++) = v1(k0, k1); }
                for (size_t k0 = 0; k0<v2.dim[0]; k0++) { va(n++) = v2(k0, k1); }
            }
        } 
        else  // concatenation along columns
        {
            assert(v1.dim[0] == v2.dim[0]);
            va.resize(v1.dim[0], v1.dim[1] + v2.dim[1]);
            
            size_t n = 0;
            for (size_t k1 = 0; k1<v1.dim[0]*v1.dim[1]; k1++) { va(n++) = v1(k1); }
            for (size_t k1 = 0; k1<v2.dim[0]*v2.dim[1]; k1++) { va(n++) = v2(k1); }
        }
        return va;
    }
    
    template<class T>
    ostream& operator <<(ostream& out,  v_array<T> v)
    {
        out<<"["<<v.dim[0]<<"x"<<v.dim[1]<<"] "<<endl;
        for (size_t kR=0; kR<v.dim[0]; kR++)
        {
            out<<"  ";
            for (size_t kC=0; kC<v.dim[1]; kC++)
            {
                if (kC<v.dim[1]-1)
                    out<<v(kR, kC)<<",";
                else
                    out<<v(kR, kC)<<std::endl;
                
            }
        }
        return out;
    }
    
    template<class T>
    inline v_array<T> cshift(v_array<T> va, int n, V_ColRow_Type cr)
    {
        if (cr==v_col)
        {
            v_size_t idx(1, va.dim[1]);
            for (size_t k=0; k<va.dim[1]; k++) idx[k] = k;
            idx.v=idx.v.cshift(n);
            
            return va(_, idx);
        } 
        else // if (cr==v_row)    
        {
            v_size_t idx(va.dim[0], 1);
            for (size_t k=0; k<va.dim[0]; k++) idx[k] = k;
            idx.v=idx.v.cshift(n);
            
            return va(idx, _);
        }
    }
    
    template<class T>
    inline v_array<T> circshift(v_array<T> va, int n, V_ColRow_Type cr) { return cshift(va, -n, cr); } // MATLAB convention opposite to C++ convention.
    
    template<class T>
    inline v_array<T> sum(v_array<T> va, V_ColRow_Type cr)
    {
        if (cr==v_row)
        {
            v_array<T> res(1, va.dim[0], 0);
            for (size_t k=0; k<va.dim[0]; k++)
            {
                res += va(k,_);
            }
            return res;
        } 
        else // if (cr==v_col)
        {
            v_array<T> res(va.dim[0], 1, 0);
            for (size_t k=0; k<va.dim[1]; k++)
            {
                res += va(_, k);
            }
            return res;
        }
    }

    inline v_size_t bsum(v_bool va, V_ColRow_Type cr)
    {
        if (cr==v_row)
        {
            v_size_t res(1, va.dim[0], 0);
            for (size_t k=0; k<va.dim[0]; k++)
            {
                res(va(k,_)) = res(va(k,_)) + 1;
            }
            return res;
        } 
        else // if (cr==v_col)
        {
            v_size_t res(va.dim[0], 1, 0);
            for (size_t k=0; k<va.dim[1]; k++)
            {
                res(va(_, k)) = res(va(_, k)) + 1 ;
            }
            return res;
        }
    }
    
    
    template<class T>
    inline T sum(const v_array<T> &va)
    {
        T res=0;
        for (size_t k=0; k<va.v.size(); k++) res += va.v[k];
        return res;
    }

    inline bool all(v_bool v1)
    {
        bool b=true;
        for (size_t k=0; k<v1.v.size(); k++) { b &=v1.v[k]; }
        return b;
    }
    
    inline bool any(v_bool v1)
    {
        bool b=false;
        for (size_t k=0; k<v1.v.size(); k++) { b |=v1.v[k]; }
        return b;
    }
    
    template<class T>
    inline bool isempty(const v_array<T> &va)
    {
        return va.dim[0] == 0 || va.dim[1] == 0;
    }

    inline v_size_t find(const std::v_bool &v)
    {
        std::vector<size_t> idx;
        for (size_t k=0; k<v.v.size(); k++) if (v[k]) idx.push_back(k);
        
        v_size_t res(&idx[0], idx.size());
        if (size(v, v_row) == 1) res.reshape(1, idx.size());
        return res;
    }
    
    inline v_size_t find_first(const std::v_bool &v)
    {
        std::vector<size_t> idx;
        for (size_t k=0; k<v.v.size(); k++) 
        {
            for (size_t k=0; k<v.v.size(); k++) if (v[k]) { idx.push_back(k); break; }
        }
        
        v_size_t res(&idx[0], idx.size());
        if (size<bool>(v, v_row) == 1) res.reshape(1, idx.size());
        return res;
    }

    template<class T>
    inline v_array<T> unique(const v_array<T> &va)
    {
        std::set<T> uniqueSet(&va.v[0], &va.v[0]+va.v.size());
        
        v_array<T> res(1,uniqueSet.size());
        size_t k=0;
        for(typename std::set<T>::iterator it=uniqueSet.begin(); it!=uniqueSet.end(); it++, k++)
        {
            res[k] = *it;            
        }
        
        return res;
    }

}

#endif
#ifndef KVC_H
#define KVC_H



#include <string>
#include <valarray>
#include <iostream>
#include <vector>
#include <stdint.h>
#include <typeinfo>
#include <complex>
#include <memory>

#include "msockets.h"
#include "varray.h"

struct MatlabVar;
template<class C> inline MatVarType get_MatVarType();

struct KVCObject
{
    std::auto_ptr<MatlabVar> mv;
    
    virtual ~KVCObject() { } 
    virtual MatVarType getMatVarType() { return MAT_EMPTY; };
    virtual bool isValid() { return false; }
};

struct KVCArray: public KVCObject
{
    size_t n;
    std::vector<size_t> dims; // dims is never empty!

    KVCArray(): n(0), dims(2, 0) { }
    virtual bool isMutable() { return false; };
    
    virtual void setDims(std::vector<size_t> d)
    {
        dims.clear(); n=0; if (d.empty()) return;
        n=1;
        dims.reserve(d.size());
        for (size_t k=0; k<d.size(); k++) { n*=d[k]; dims.push_back(d[k]); }
    }
    
    virtual void setDims(size_t d1, size_t d2 = 1)
    {
        dims.clear(); dims.reserve(2);
        n = d1*d2;
        dims.push_back(d1); dims.push_back(d2);
    }
    
    void setToRowVector() 
    {
        dims.clear(); dims.push_back(n); dims.push_back(1);
    }
    
    void setToColumnVector()
    {
        dims.clear(); dims.push_back(1); dims.push_back(n);
    }
    
    virtual bool isValid()
    {
        if (dims.size()==0) return false;
        return true;
    }
};

template<class T>
struct KVCImmutableArray: public KVCArray
{
    T *pr, *pi;

    KVCImmutableArray(): pr(0), pi(0) {}
    virtual MatVarType getMatVarType() { return get_MatVarType<T>(); };
    
    virtual bool isValid() { if (!KVCArray::isValid()) return false; size_t sum=dims.size()==0?0:1; for (size_t k=0; k<dims.size(); k++) sum*=dims[k]; return (this->n==sum); }
};

template<class T>
struct KVCRealArray: public KVCImmutableArray<T>
{
    T &operator[](size_t idx) { assert(idx<this->n); return( ((T *)this->pr)[idx]); }
};

template<class T>
struct KVCComplexArray: public KVCImmutableArray<T>
{
    std::complex<T> operator[](size_t idx) { assert(idx<this->n); return( std::complex<T>(((T *)this->pr)[idx], ((T *)this->pi)[idx]) ); }
    void set(size_t idx, std::complex<T> v) { assert(idx<this->n); ((T *)this->pr)[idx] = v.real(); ((T *)this->pi)[idx] = v.imag(); }
};

template<class T>
struct KVCIndexArray: public KVCImmutableArray<T>
{
    const T operator[](size_t idx) { assert(idx<this->n); return(((T *)this->pr)[idx]-1); }
    void set(size_t idx, T val) { assert(idx<this->n); ((T *)this->pr)[idx] = val+1; }
};


template<class T>
struct KVCMutableArray: public KVCArray
{
    void fixDims() // used after push_back
    {
        if ((this->dims.size())!=2)
        {
            this->dims.clear(); this->dims.push_back(0); this->dims.push_back(0);
        }
        if (this->dims[0]>1) 
        {
            this->dims[0]=this->n; this->dims[1]=1;
        } else {
            this->dims[0]=1; this->dims[1]=this->n;
        }
    }        
    
    
    virtual MatVarType getMatVarType() { return get_MatVarType<T>(); };
    virtual bool isMutable() { return true; };
};






template<class MT, class T>
struct KVCMutableVarray: public KVCMutableArray<MT>, public std::v_array<T>
{
    void setDimsFromVArray()
    {
        KVCArray::setDims(this->dim[0], this->dim[1]);
    }
    
    virtual void setDims(std::vector<size_t> d)
    {
        KVCArray::setDims(d);
        std::v_array<T>::resize(d[0], d[1]);
    }

    virtual void setDims(size_t d1, size_t d2 = 1)
    {
        KVCArray::setDims(d1, d2);
        std::v_array<T>::resize(d1, d2);
//        this->nR = this->dims[0]; this->nC = this->dims[1];
    }

    virtual bool isValid() 
    { 
        size_t sum=this->dims.size()==0?0:1; 
        for (size_t k=0; k<this->dims.size(); k++) sum*=this->dims[k]; 
        return ((this->n==sum)&&(sum==this->v.size())); 
    }
    
    void resizeDim(size_t d1, size_t d2, T c=T()) 
    { 
        this->dims.clear(); this->dims.push_back(d1); this->dims.push_back(d2);
        this->n=d1*d2;
        std::v_array<T>::resize(d1, d2, c);
    }
    
    void operator =(const std::v_array<T> &o) { 
        this->dims[0] = o.dim[0]; 
        this->dims[1] = o.dim[1]; 
        std::v_array<T>::operator =(o);
        this->n = o.v.size();
    }
};

template<class T> struct KVCMutableRealArray: public KVCMutableVarray<T, T> 
{ 
    void operator =(const std::v_array<T> &o) { KVCMutableVarray<T, T>::operator =(o); }
};
template<class T> struct KVCMutableComplexArray: public KVCMutableVarray<T, std::complex<T> > 
{ 
    void operator =(const std::v_array<T> &o) { KVCMutableVarray<T, std::complex<T> >::operator =(o); }
};
template<class T> struct KVCMutableIndexArray: public KVCMutableVarray<T, size_t> 
{ 
    void operator =(const std::v_array<T> &o) { KVCMutableVarray<T, size_t>::operator =(o); }
    void operator =(const KVCMutableIndexArray<T> &o) { KVCMutableVarray<T, size_t>::operator =(o); }
};


#define INIT_KVC(k, v) KVC(k, v, sizeof(v)/sizeof(KVCObject **))

struct KVCStruct: public KVCObject
{
    size_t        KVC_n;
    std::string  *KVC_keys;
    KVCObject   **KVC_vals;
    
    KVCObject *valueForKey(std::string s)
    {
        for (int k=0; k<KVC_n; k++) if (KVC_keys[k]==s) return KVC_vals[k];
        return 0;
    }
    
    virtual MatVarType getMatVarType() { return MAT_STRUCT; };
    
    template<class T>
    void setValueForKey(std::string s, T *v2)
    {
        for (int k=0; k<KVC_n; k++) if (KVC_keys[k]==s) { T *v = dynamic_cast<T *>(KVC_vals[k]); if (v) *v = *v2; }
    }
    
    void KVC(std::string *keys, KVCObject **vals, size_t n)
    {
        KVC_keys = new std::string[n];
        KVC_vals = new KVCObject *[n];
        for (size_t k=0; k<n; k++)
        {
            KVC_keys[k]=keys[k];
            KVC_vals[k]=vals[k];
        }
        KVC_n = n;
    }
    
    virtual bool isValid() { for (size_t k=0; k<KVC_n; k++) if (!KVC_vals[k]->isValid()) return false; return true; }
    
    virtual ~KVCStruct() { 
        delete [] KVC_vals; 
        delete [] KVC_keys; 
    }
};

typedef KVCStruct KVCMutableStruct;

struct KVCCellBase: public KVCArray
{
    std::vector<KVCObject *> children;
    
    void clear()
    {
        for (size_t k=0; k<children.size(); k++) delete children[k]; 
        children.clear();
    }
    
    virtual void allocateElements() = 0;
    virtual ~KVCCellBase() { 
        clear();
    }
    
    virtual bool isValid(){ 
        if (children.size()!=n) return false; 
        for (size_t k=0; k<children.size(); k++) 
            if (!children[k]->isValid()) return false;
        return true;
    }
};

template<class C>
struct KVCCell: public KVCCellBase
{
    C *newChild() { children.push_back(new C()); this->n++; return (C *)children.back(); } // we are the owner of the child objects now!
    C *getChild(size_t idx) { return (C *)children[idx]; }
    
    virtual void allocateElements()
    {
        for (size_t k=0; k<children.size(); k++) delete children[k]; 
        children.clear();
        for (size_t k=0; k<n; k++) children.push_back(new C());
    }
};





template<class C> inline MatVarType get_MatVarType()    { return MAT_EMPTY; };
template<> inline MatVarType get_MatVarType<bool>()     { return MAT_LOGICAL; }
template<> inline MatVarType get_MatVarType<char>()     { return MAT_CHAR; }
template<> inline MatVarType get_MatVarType<uint8_t>()  { return MAT_UINT8; }
template<> inline MatVarType get_MatVarType<int8_t>()   { return MAT_INT8; }
template<> inline MatVarType get_MatVarType<uint16_t>() { return MAT_UINT16; }
template<> inline MatVarType get_MatVarType<int16_t>()  { return MAT_INT16; }
template<> inline MatVarType get_MatVarType<uint32_t>() { return MAT_UINT32; }
template<> inline MatVarType get_MatVarType<int32_t>()  { return MAT_INT32; }
template<> inline MatVarType get_MatVarType<uint64_t>() { return MAT_UINT64; }
template<> inline MatVarType get_MatVarType<int64_t>()  { return MAT_INT64; }
template<> inline MatVarType get_MatVarType<float>()    { return MAT_FLOAT32; }
template<> inline MatVarType get_MatVarType<double>()   { return MAT_DOUBLE; }

#endif
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
//    size_t n;
//    std::vector<size_t> dims; // dims is never empty!

    KVCArray() {}  //: n(0), dims(2, 0) { }
    virtual bool isMutable() { return false; };

    virtual size_t nDim() = 0;
    virtual size_t &dim(size_t nD) = 0;
    virtual void setDims(std::vector<size_t> d) = 0;
    virtual void setDims(size_t d1, size_t d2 = 1) = 0;
    virtual size_t size() = 0;

    virtual bool isValid()
    {
        if (this->size()==0) return false;
        return true;
    }
};


template<class MT, class T>
struct KVCVarray: public KVCArray, public std::v_array<T>
{
    virtual size_t nDim() { return 2; }
    
    virtual size_t &dim(size_t nD)
    {
        return std::v_array<T>::dim[nD];
    }
    
    virtual void setDims(std::vector<size_t> d)
    {
        std::v_array<T>::resize(d[0], d[1]);
    }

    virtual void setDims(size_t d1, size_t d2 = 1)
    {
        std::v_array<T>::resize(d1, d2);
    }
    
    virtual size_t size() { return std::v_array<T>::v.size(); }

    virtual bool isValid() 
    { 
        return (this->size()==this->v.size());
    }
    
    void operator =(const std::v_array<T> &o) {
        this->setDims(o.nR(), o.nC());
        std::v_array<T>::operator =(o);
    }
    
    virtual MatVarType getMatVarType() { return get_MatVarType<MT>(); };
};

template<class T> struct KVCRealArray: public KVCVarray<T, T>
{ 
    void operator =(const std::v_array<T> &o) { KVCVarray<T, T>::operator =(o); }
};
template<class T> struct KVCComplexArray: public KVCVarray<T, std::complex<T> >
{ 
    void operator =(const std::v_array<T> &o) { KVCVarray<T, std::complex<T> >::operator =(o); }
};
template<class T> struct KVCIndexArray: public KVCVarray<T, size_t>
{ 
    void operator =(const std::v_array<T> &o) { KVCVarray<T, size_t>::operator =(o); }
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
    
    virtual size_t nDim() { return 1; }
    
    virtual size_t &dim(size_t nD)
    {
        return KVC_n;
    }
    
    virtual void setDims(std::vector<size_t> d)
    {
        std::cerr<<"setDims not implemented for KVCStruct!"<<std::endl;
    }
    virtual void setDims(size_t d1, size_t d2 = 1)
    {
        std::cerr<<"setDims not implemented for KVCStruct!"<<std::endl;
    }
    
    virtual size_t size() { return KVC_n; }
    
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
    std::vector<size_t> dimensions;
    
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

    virtual size_t nDim() { return dimensions.size(); }
    virtual size_t &dim(size_t nD)
    {
        return dimensions[nD];
    }
    
    virtual void setDims(std::vector<size_t> d)
    {
        dimensions = d;
    }
    
    virtual void setDims(size_t d1, size_t d2 = 1)
    {
        dimensions.clear(); dimensions.push_back(d1); dimensions.push_back(d2);
    }

    virtual size_t size()
    {
        size_t N = (this->nDim()==0) ? 0 : 1;
        for (size_t k=0; k<this->nDim(); k++) N*=this->dim(k);
        return N;
    }

    virtual bool isValid()
    {
        if (children.size()!=this->size()) return false;

        for (size_t k=0; k<children.size(); k++)
        {
            if (!children[k]->isValid()) return false;
        }
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
        
        size_t sum = (dimensions.size()==0) ? 0 : 1;
        for (size_t k=0; k<this->dims().size(); k++) sum*=this->dim(k);

        for (size_t k=0; k<sum; k++) children.push_back(new C());
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
#include "kvcMatlabVar.h"
#include "msockets.h"

#include <memory>
using namespace std;

typedef bool (*MatVarBuilderFct_t)(MatlabVar *parent,  KVCArray *srcAr, KVCArray *tgtAr);

typedef bool (*KVCObject_from_MatlabVar_Handler)(KVCObject *kvco, MatlabVar *mv);

typedef bool (*MatlabVar_from_KVCArray_Handler)(MatVarType mt, KVCArray *ko);

static bool lookups_initialized;
static std::map<MatVarType, KVCObject_from_MatlabVar_Handler> lookup_init_KVCObject_from_MatlabVar;
static std::map<MatVarType, MatlabVar_from_KVCArray_Handler> lookup_create_MatlabVar_from_KVCArray;

static std::map<std::string, MatVarType> lookup_MatlabVar_from_KVC_eltype;

void setup_KVC_lookups();


bool init_Nothing_from_MatlabVar(KVCObject *kvco, MatlabVar *mv)
{
    return true;
}

bool init_KVCArray_from_MatlabVar(KVCArray *ka, MatlabVar *mv)
{
    std::vector<size_t> d;
    for (size_t k=0; k<*mv->ndim; k++) d.push_back(mv->dims[k]);   
    ka->setDims(d);

    return true;
}


template<class C>
bool init_KVCArray_from_MatlabVar(KVCObject *kvco, MatlabVar *mv)
{
    C *pr = (C *)mv->pr, *pi = (C *)mv->pi;
    
    KVCImmutableArray<C> *kia = dynamic_cast<KVCImmutableArray<C> *>(kvco);
    if (kia)
    {
        init_KVCArray_from_MatlabVar(kia, mv);
        kia->pr = (C *)mv->pr;
        kia->pi = (C *)mv->pi;
        
        kia->mv.reset(mv);
        
        return true;
    }
    
    if (*mv->varComplexity==MAT_COMPLEX)
    {
        KVCMutableComplexArray<C> *ka = dynamic_cast<KVCMutableComplexArray<C> *>(kvco); if (!ka) return false;
        init_KVCArray_from_MatlabVar(ka, mv);
        ka->v.resize(ka->n);
        for (size_t k=0; k<ka->n; k++)
        {
            (*ka)[k] = std::complex<C>(pr[k], pi[k]);
        }
        return true;
    } 
    
    KVCMutableRealArray<C> *ka = dynamic_cast<KVCMutableRealArray<C> *>(kvco);
    if (ka)
    {
        init_KVCArray_from_MatlabVar(ka, mv);
        ka->v.resize(ka->n);
        for (size_t k=0; k<ka->n; k++)
        {
            (*ka)[k] = pr[k];
        }
        return true;
    }
    
    KVCMutableIndexArray<C> *kmia = dynamic_cast<KVCMutableIndexArray<C> *>(kvco);
    if (kmia)
    {
        init_KVCArray_from_MatlabVar(kmia, mv);
        kmia->v.resize(kmia->n);
        for (size_t k=0; k<kmia->n; k++)
        {
            (*kmia)[k] = size_t(pr[k] - 1); // index base 1 to index base 0
        }
        return true;
    }
    
    std::cout<<"Type mismatch - could not initialize KVCObject of type "<<MatVarStrings[kvco->getMatVarType()]<<" with a MatlabVar of type "<<MatVarStrings[*mv->varType]<<std::endl;
    return false;
}


bool init_KVCStruct_from_MatlabVar(KVCObject *kvco, MatlabVar *mv)
{
    if (!lookups_initialized) setup_KVC_lookups();
    KVCStruct *ks = dynamic_cast<KVCStruct *>(kvco); if (!ks) return false;
    
    for (int k=0; k<mv->fieldNames.size(); k++)
    {
        KVCObject *child_kvco = ks->valueForKey( mv->fieldNames[k] );
        if (!child_kvco) continue;
        
        KVCObject_from_MatlabVar_Handler fct = lookup_init_KVCObject_from_MatlabVar[*mv->children[k]->varType];
        
        if (!fct)
        {
            std::cerr<<"Don't know how to handle "<<mv->fieldNames[k]<<"... Ignoring."<<std::endl;
            continue;
        }
        
        if (fct == &init_Nothing_from_MatlabVar)
        {
            std::cout<<"No handler for "<<mv->fieldNames[k]<<"... Ignoring. "<<std::endl;
            continue;
        }
        if (!fct(child_kvco, mv->children[k])) 
        {
            std::cout<<"Error initializing child "<<mv->fieldNames[k]<<std::endl;
            return false;
        }
    }
    return true;
}


bool init_KVCCell_from_MatlabVar(KVCObject *kvco, MatlabVar *mv)
{
    KVCCellBase *kc = dynamic_cast<KVCCellBase *>(kvco); 
    if (!kc) return false;

    init_KVCArray_from_MatlabVar(kc, mv);

    kc->allocateElements();
    
    for (int k=0; k<kc->n; k++)
    {
        KVCObject *child_kvco = kc->children[k];
        if (!child_kvco) continue;
        
        KVCObject_from_MatlabVar_Handler fct = lookup_init_KVCObject_from_MatlabVar[*mv->children[k]->varType];
        if (!fct)
        {
            std::cerr<<"Don't know how to handle for child"<<mv->fieldNames[k]<<"... Ignoring."<<std::endl;
            continue;
        }
        
        if (fct == &init_Nothing_from_MatlabVar)
        {
            std::cout<<"No handler for "<<mv->fieldNames[k]<<"... Ignoring. "<<std::endl;
            continue;
        }
        if (!fct(child_kvco, mv->children[k])) return false;
    }
    
    return true;
}


auto_ptr<mwSize> buildMWSizeDims(std::vector<size_t> &d)
{
    if (d.empty()) return auto_ptr<mwSize>(0);
    mwSize *dims = new mwSize[d.size()];
    for (size_t k=0; k<d.size(); k++) 
    {
        dims[k] = (mwSize)d[k];
    }
    return auto_ptr<mwSize>(dims);
}


template<class C>
inline bool create_MatlabVar_from_KVCRealArray(MatVarType mt, KVCRealArray<C> *ka)
{
    ka->mv.reset(new MatlabVar(mt, MAT_REAL, (mwSize) ka->dims.size(), buildMWSizeDims(ka->dims).get()));
    ka->pr = (C *)(ka->mv->pr);
    
    return true;
}

template<class C>
inline bool create_MatlabVar_from_KVCIndexArray(MatVarType mt, KVCIndexArray<C> *ka)
{
    ka->mv.reset(new MatlabVar(mt, MAT_REAL, (mwSize) ka->dims.size(), buildMWSizeDims(ka->dims).get()));
    ka->pr = (C *)ka->mv->pr;
    
    return true;
}

template<class C>
inline bool create_MatlabVar_from_KVCComplexArray(MatVarType mt, KVCComplexArray<C> *ka)
{
    ka->mv.reset(new MatlabVar(mt, MAT_COMPLEX, (mwSize) ka->dims.size(), buildMWSizeDims(ka->dims).get()));
    ka->pr = (C *)ka->mv->pr;
    ka->pi = (C *)ka->mv->pi;
    
    return true;
}

template<class C>
inline bool create_MatlabVar_from_KVCMutableRealArray(MatVarType mt,  KVCMutableRealArray<C> *ka)
{
    std::vector<mwSize> dims; for (size_t k=0; k<ka->dims.size(); k++) dims.push_back((mwSize)ka->dims[k]); // convert to mwSize
    
    MatlabVar *mv = new MatlabVar(mt, MAT_REAL, (mwSize)dims.size(), &dims[0]);
    ka->mv.reset(mv);
    
    C *pr = (C *)mv->pr;
    
    for (size_t k=0; k<ka->n; k++)
    {
        pr[k] = (*ka)[k];
    }
    
    return true;
}

template<class C>
inline bool create_MatlabVar_from_KVCMutableIndexArray(MatVarType mt,  KVCMutableIndexArray<C> *ka)
{
    std::vector<mwSize> dims; for (size_t k=0; k<ka->dims.size(); k++) dims.push_back((mwSize)ka->dims[k]);
    
    MatlabVar *mv = new MatlabVar(mt, MAT_REAL, (mwSize)dims.size(), &dims[0]);
    ka->mv.reset(mv);

    C *pr = (C *)mv->pr;

    for (size_t k=0; k<ka->n; k++)
    {
        pr[k] = C((*ka)[k]) + 1;  // index base 0 to index base 1;
    }

    return true;
}

template<class C>
inline bool create_MatlabVar_from_KVCMutableComplexArray(MatVarType mt,  KVCMutableComplexArray<C> *ka)
{
    std::vector<mwSize> dims; for (size_t k=0; k<ka->dims.size(); k++) dims.push_back((mwSize)ka->dims[k]);

    MatlabVar *mv = new MatlabVar(mt, MAT_COMPLEX, (mwSize)dims.size(), &dims[0]);
    ka->mv.reset(mv);

    C *pr = (C *)mv->pr, *pi = (C *)mv->pi;

    for (size_t k=0; k<ka->n; k++)
    {
        pr[k] = (*ka)[k].real();
        pi[k] = (*ka)[k].imag();
    }
    
    return true;
}

template<class C>
bool create_MatlabVar_from_KVCArray(MatVarType mt, KVCArray *ka)
{
    KVCMutableRealArray<C> *kmra = dynamic_cast<KVCMutableRealArray<C> *>(ka);
    if (kmra) return create_MatlabVar_from_KVCMutableRealArray(mt, kmra); 
    
    KVCMutableIndexArray<C> *kmia = dynamic_cast<KVCMutableIndexArray<C> *>(ka);
    if (kmia) return create_MatlabVar_from_KVCMutableIndexArray(mt, kmia);
    
    KVCRealArray<C> *kra = dynamic_cast<KVCRealArray<C> *>(ka);
    if (kra) return create_MatlabVar_from_KVCRealArray(mt, kra);
    
    KVCIndexArray<C> *kia = dynamic_cast<KVCIndexArray<C> *>(ka);
    if (kia) return create_MatlabVar_from_KVCIndexArray(mt, kia);
    
    KVCMutableComplexArray<C> *kmca = dynamic_cast<KVCMutableComplexArray<C> *>(ka);
    if (kmca) return create_MatlabVar_from_KVCMutableComplexArray(mt, kmca);

    KVCComplexArray<C> *kca = dynamic_cast<KVCComplexArray<C> *>(ka);
    if (kca) return create_MatlabVar_from_KVCComplexArray(mt, kca);
    
    return false;
}


bool create_MatlabVar_from_KVCCell(KVCCellBase *kc);


bool create_MatlabVar_from_KVCStruct(KVCStruct *ks)
{
    if (!lookups_initialized) setup_KVC_lookups();
    
    char * *cfNames = new char *[ks->KVC_n]; 
    for (mwSize k=0; k<ks->KVC_n; k++) cfNames[k] = (char *)ks->KVC_keys[k].c_str();
    MatlabVar *mv = new MatlabVar(MAT_STRUCT, (mwSize)ks->KVC_n, (const char **)cfNames); 
    delete [] cfNames;
    
    for (int k=0; k<ks->KVC_n; k++)
    {
        KVCArray *ka = dynamic_cast<KVCArray *>(ks->KVC_vals[k]);
        if (ka)
        {
            MatVarType mt = ka->getMatVarType();
            
            if (mt==MAT_EMPTY)
            {
                KVCCellBase *kcb = dynamic_cast<KVCCellBase *> (ks->KVC_vals[k]);
                if (kcb)
                {
                    create_MatlabVar_from_KVCCell(kcb);
                    mv->addChild( kcb->mv );
                    continue;
                }
                
                std::cerr<<"Don't know how to handle "<<ks->KVC_keys[k]<<"... Ignoring."<<std::endl;
                mv->addChild(auto_ptr<MatlabVar>(new MatlabVar()));
                continue;
            }
            
            MatlabVar_from_KVCArray_Handler fct = lookup_create_MatlabVar_from_KVCArray[mt];
                
            if (!fct)
            {
                std::cout<<"No handler for "<<ks->KVC_keys[k]<<"... Ignoring. "<<std::endl;
                mv->addChild(auto_ptr<MatlabVar>(new MatlabVar()));
                continue;
            }
                
            if (!fct(mt, ka)) return false;
            mv->addChild(ka->mv);
            
            continue;
        } 
        
        KVCStruct *kcs = dynamic_cast<KVCStruct *>(ks->KVC_vals[k]);
        if (kcs)
        {
            create_MatlabVar_from_KVCStruct(kcs);
            mv->addChild(kcs->mv);
            continue;
        }
        
        std::cerr<<"Couldn't find type for cell #"<<k<<"... Ignoring. "<<std::endl;
        mv->addChild(auto_ptr<MatlabVar>(new MatlabVar()));
    }

    ks->mv.reset(mv);
    
    return true;
}

bool create_MatlabVar_from_KVCCell(KVCCellBase *kc)
{
    if (!lookups_initialized) setup_KVC_lookups();

    MatlabVar *mv = new MatlabVar(MAT_CELL, (mwSize)kc->dims.size(), buildMWSizeDims(kc->dims).get());
    
    for (int k=0; k<kc->children.size(); k++)
    {
        KVCArray *ka = dynamic_cast<KVCArray *>(kc->children[k]);
        if (ka)
        {
            MatVarType mt = ka->getMatVarType();
            
            if (mt==MAT_EMPTY)
            {
                KVCCellBase *kcb = dynamic_cast<KVCCellBase *> (kc->children[k]);
                if (kcb)
                {
                    create_MatlabVar_from_KVCCell(kcb);
                    mv->addChild(kcb->mv);
                    continue;
                }
                
                std::cerr<<"Don't know how to handle cell #"<<k<<"... Ignoring."<<std::endl;
                mv->addChild(auto_ptr<MatlabVar>(new MatlabVar()));
                continue;
            }
            
            MatlabVar_from_KVCArray_Handler fct = lookup_create_MatlabVar_from_KVCArray[mt];
            
            if (!fct)
            {
                std::cout<<"No handler for cell #"<<k<<"... Ignoring. "<<std::endl;
                mv->addChild(auto_ptr<MatlabVar>(new MatlabVar()));
                continue;
            }

            if (!fct(mt, ka)) return false;
            mv->addChild(ka->mv);
            continue;
        } 
        
        KVCStruct *kcs = dynamic_cast<KVCStruct *>(kc->children[k]);
        if (kcs)
        {
            create_MatlabVar_from_KVCStruct(kcs);
            mv->addChild(kcs->mv);
            continue;
        }
        
        std::cerr<<"Couldn't find type for cell #"<<k<<"... Ignoring. "<<std::endl;
        mv->addChild(auto_ptr<MatlabVar>(new MatlabVar()));
    }
    
    kc->mv.reset(mv);
    
    return true;
}

void setup_KVC_lookups()
{
    lookup_init_KVCObject_from_MatlabVar[MAT_EMPTY]   = &init_Nothing_from_MatlabVar;
    lookup_init_KVCObject_from_MatlabVar[MAT_LOGICAL] = &init_KVCArray_from_MatlabVar<bool>;
    lookup_init_KVCObject_from_MatlabVar[MAT_CHAR]    = &init_KVCArray_from_MatlabVar<char>;
    lookup_init_KVCObject_from_MatlabVar[MAT_UINT8]   = &init_KVCArray_from_MatlabVar<uint8_t>;
    lookup_init_KVCObject_from_MatlabVar[MAT_INT8]    = &init_KVCArray_from_MatlabVar<int8_t>;
    lookup_init_KVCObject_from_MatlabVar[MAT_UINT16]  = &init_KVCArray_from_MatlabVar<uint16_t>;
    lookup_init_KVCObject_from_MatlabVar[MAT_INT16]   = &init_KVCArray_from_MatlabVar<int16_t>;
    lookup_init_KVCObject_from_MatlabVar[MAT_UINT32]  = &init_KVCArray_from_MatlabVar<uint32_t>;
    lookup_init_KVCObject_from_MatlabVar[MAT_INT32]   = &init_KVCArray_from_MatlabVar<int32_t>;
    lookup_init_KVCObject_from_MatlabVar[MAT_UINT64]  = &init_KVCArray_from_MatlabVar<uint64_t>;
    lookup_init_KVCObject_from_MatlabVar[MAT_INT64]   = &init_KVCArray_from_MatlabVar<int64_t>;
    lookup_init_KVCObject_from_MatlabVar[MAT_FLOAT32] = &init_KVCArray_from_MatlabVar<float>;
    lookup_init_KVCObject_from_MatlabVar[MAT_SINGLE]  = &init_KVCArray_from_MatlabVar<float>;
    lookup_init_KVCObject_from_MatlabVar[MAT_FLOAT64] = &init_KVCArray_from_MatlabVar<double>;
    lookup_init_KVCObject_from_MatlabVar[MAT_DOUBLE]  = &init_KVCArray_from_MatlabVar<double>;
    lookup_init_KVCObject_from_MatlabVar[MAT_STRUCT]  = &init_KVCStruct_from_MatlabVar;
    lookup_init_KVCObject_from_MatlabVar[MAT_CELL]    = &init_KVCCell_from_MatlabVar;
        
    lookup_create_MatlabVar_from_KVCArray[MAT_LOGICAL] = &create_MatlabVar_from_KVCArray<bool>;
    lookup_create_MatlabVar_from_KVCArray[MAT_CHAR]    = &create_MatlabVar_from_KVCArray<char>;
    lookup_create_MatlabVar_from_KVCArray[MAT_UINT8]   = &create_MatlabVar_from_KVCArray<uint8_t>;
    lookup_create_MatlabVar_from_KVCArray[MAT_INT8]    = &create_MatlabVar_from_KVCArray<int8_t>;
    lookup_create_MatlabVar_from_KVCArray[MAT_UINT16]  = &create_MatlabVar_from_KVCArray<uint16_t>;
    lookup_create_MatlabVar_from_KVCArray[MAT_INT16]   = &create_MatlabVar_from_KVCArray<int16_t>;
    lookup_create_MatlabVar_from_KVCArray[MAT_UINT32]  = &create_MatlabVar_from_KVCArray<uint32_t>;
    lookup_create_MatlabVar_from_KVCArray[MAT_INT32]   = &create_MatlabVar_from_KVCArray<int32_t>;
    lookup_create_MatlabVar_from_KVCArray[MAT_UINT64]  = &create_MatlabVar_from_KVCArray<uint64_t>;
    lookup_create_MatlabVar_from_KVCArray[MAT_INT64]   = &create_MatlabVar_from_KVCArray<int64_t>;
    lookup_create_MatlabVar_from_KVCArray[MAT_FLOAT32] = &create_MatlabVar_from_KVCArray<float>;
    lookup_create_MatlabVar_from_KVCArray[MAT_SINGLE]  = &create_MatlabVar_from_KVCArray<float>;
    lookup_create_MatlabVar_from_KVCArray[MAT_FLOAT64] = &create_MatlabVar_from_KVCArray<double>;
    lookup_create_MatlabVar_from_KVCArray[MAT_DOUBLE]  = &create_MatlabVar_from_KVCArray<double>;
    
    lookups_initialized = true;
}




bool init_KVCObject_from_MatlabVar(KVCObject *kvco, auto_ptr<MatlabVar> mv)
{
    if (!lookups_initialized) setup_KVC_lookups();
    
    KVCObject_from_MatlabVar_Handler fct = lookup_init_KVCObject_from_MatlabVar[*mv->varType];
    
    if ( (!fct) || (fct == &init_Nothing_from_MatlabVar) )
    {
        std::cerr<<"Don't know how to handle MatlabVar..."<<std::endl;
        return false;
    }
    
    return fct(kvco, mv.release());
}


bool create_MatlabVar_from_KVCObject(KVCObject *kvco)
{
    if (!lookups_initialized) setup_KVC_lookups();
    
    if (!kvco->isValid()) 
    {
        std::cerr<<"KVCObject is invalid!"<<std::endl;
        return false;
    }
    
    KVCArray *ka = dynamic_cast<KVCArray *>(kvco);
    if (ka)
    {
        MatVarType mt = ka->getMatVarType();
        
        
        if (mt==MAT_EMPTY)
        {
            KVCCellBase *kcb = dynamic_cast<KVCCellBase *> (kvco);
            if (kcb)
            {
                return create_MatlabVar_from_KVCCell(kcb);
            }

            std::cerr<<"Don't know how to handle KVCArray..."<<std::endl;
            return false;
        }
        
        MatlabVar_from_KVCArray_Handler fct = lookup_create_MatlabVar_from_KVCArray[mt];
        
        if (!fct)
        {
            std::cout<<"No handler for KVCArray..."<<std::endl;
            return false;
        }
        return fct(mt, ka);
    } 
    
    KVCStruct *kcs = dynamic_cast<KVCStruct *>(kvco);
    if (kcs)
    {
        return create_MatlabVar_from_KVCStruct(kcs);
    }
    
    return false; // should not happen
}


bool init_KVCObject_from_Socket(KVCObject *kvco, int remote_socket)
{
    auto_ptr<MatlabVar> mv(new MatlabVar());
    
    if (!socketRecv(remote_socket, mv.get())) return false;
    
    return init_KVCObject_from_MatlabVar(kvco, mv);
}


bool send_KVCObject_to_Socket(KVCObject *kvco, int remote_socket)
{
    if (!kvco->mv.get()) // does a MatlabVar backend already exist?
    {
        if (!create_MatlabVar_from_KVCObject(kvco)) return false;
    }
    
    return socketSend(remote_socket, kvco->mv.get());
}
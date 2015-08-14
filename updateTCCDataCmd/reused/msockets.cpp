/*
 *  msockets.cpp
 *  SocketServer
 *
 *  Created by Denis Kovacs on 10/6/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include "msockets.h"

#include <assert.h>
#include <pthread.h>

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <sys/signal.h>
#include <sys/un.h>
#include <unistd.h>
#include <netinet/in.h>
#include <unistd.h>
#include <string.h>
#include <netdb.h>
#include <arpa/inet.h>


#define SOCK_PATH "echo_socket"

using namespace std;

std::string MatVarStrings[] = {
	"MAT_EMPTY",
	"MAT_LOGICAL",
	"MAT_CHAR",
	"MAT_UINT8",
	"MAT_INT8",
	"MAT_UINT16",
	"MAT_INT16",
	"MAT_UINT32",
	"MAT_INT32",
	"MAT_UINT64",
	"MAT_INT64",
	"MAT_SINGLE",
	"MAT_DOUBLE",
	"MAT_STRUCT",
	"MAT_CELL",
};


void 
socketDumpSocketBlock(struct SocketBlock d)
{
    int k;
    
    printf("Size: %d, Content :\n", (int) d.block_size);
    for (k=0; k<d.block_size; k++)
    {
        printf("%03d ", (unsigned int)d.block_data[k]);
        if (k%4==3) printf(",");
        if (k%90==89) printf("\n");
    }
    printf("\n");
};



#pragma mark Socket Server


void socketAllocateBlock(SocketBlock &b, size_t size)
{
    b.block_data = new unsigned char[ size ];
    if (b.block_data) b.block_reserved = size;
}

void socketDeallocateBlock(SocketBlock &b)
{
    if (b.block_reserved) delete [] b.block_data;
    b.block_data = 0;
    b.block_reserved = b.block_size = 0; 
}


int
socketReceiveBlock(int remote_socket, struct SocketBlock &d)
{
    // Receive the count
    int cnt = 0, ret, recvlen;
    while (cnt < (int) sizeof(int)) {
        ret = (int) recv(remote_socket,((char *)&recvlen)+cnt, sizeof(int)-cnt,0);
        
        if (ret <= 0) {
            perror("recv");
            return 0;
        }
        cnt += ret;
    }
    
    if(recvlen <= 0) {
        return 0;
    }
    
    if (d.block_reserved < recvlen)
    {
        if (d.block_reserved) socketDeallocateBlock(d);
        socketAllocateBlock(d, recvlen);
        if (!d.block_reserved) return 0;
    } 
    d.block_size = recvlen;
    
    cnt = 0;
    while(cnt < recvlen)
    {
        ret = (int) recv(remote_socket, (d.block_data+cnt), recvlen-cnt,0);
        if (ret <= 0) {
            perror("recv");
            return 0;
        }
        cnt += ret;
    }
    
    return recvlen + (int)sizeof(int);
}   



int
socketSendBlock(int remote_socket, struct SocketBlock &d)
{
    int ret;
    
    ret = (int) send(remote_socket,(const char *)&d.block_size,sizeof(int),0);
    if(ret == -1) {
        perror("send");
        return ret;
    }
    
    int cnt = 0;
    while(cnt < d.block_size) {
        ret = (int) send(remote_socket,d.block_data+cnt,d.block_size-cnt,0);
        if(ret == -1) {
            perror("send");
            return -1;
        }
        cnt += ret;
    }
    
    return cnt + (int) sizeof(int);
}


int
socketSendMultiBlock(int remote_socket, size_t nBlocks, struct SocketBlock **d)
{
    // compute total data size
    int datasize = 0, k;
    for (k=0; k<nBlocks; k++)
        datasize += (int) d[k]->block_size;
    
    int ret, nSent = 0;
    
    ret = (int) send(remote_socket,(const char *)&datasize,sizeof(int),0);
    if(ret == -1) {
        perror("send");
        return ret;
    }
    
    nSent += (int) sizeof(int);
    
    for (k=0; k<nBlocks; k++)
    {
        int cnt = 0;
        while(cnt < d[k]->block_size) {
            ret = (int) send(remote_socket,d[k]->block_data + cnt, d[k]->block_size - cnt,0);
            if(ret == -1) {
                perror("send");
                return ret;
            }
            cnt += ret;
        }
        
        nSent += cnt;
    }
    
    return nSent;
}



void
socketWaitForUpdate(struct SocketServerInfo *srvInfo)
{
    pthread_cond_wait(&srvInfo->data_update_cv, &srvInfo->data_update_mutex);
}



int
socketServerAccept(int local_socket)
{
    struct sockaddr_in remote_addr;
    socklen_t t = sizeof(remote_addr);
    return accept(local_socket, (struct sockaddr *)&remote_addr, &t);
}



int
socketConnect(int port, char *hostname)
{
    struct sockaddr_in pin;
    struct hostent *hp;
    int sock;
    
	memset(&pin,0,sizeof(pin));
	pin.sin_family = AF_INET;
	pin.sin_port = htons(port);
	if((hp = gethostbyname(hostname))!=0) 
		pin.sin_addr.s_addr = 
        ((struct in_addr *)hp->h_addr)->s_addr;
	else 
		pin.sin_addr.s_addr = inet_addr(hostname);
    
    sock = (int)socket(AF_INET,SOCK_STREAM,0);
	
	if(sock == -1) {
		perror("socket");
		return -1;
	}
	
	if(connect(sock,(const struct sockaddr *)&pin,sizeof(pin))) {
		perror("connect");
		close(sock);
		return -1;
	}
    
    return sock;
}

void* 
socketThreadMain(void *data)
{
    struct SocketServerInfo *srvInfo = (struct SocketServerInfo *)data;
    
    int ret=0;
    
    char str[101];
    
    str[100]=0;
    
    for(;;)
    {
        for (;;)
        {
            printf("Waiting for connection...\n");
            if (srvInfo->isServer)
            {
                srvInfo->socket = socketServerAccept(srvInfo->server_socket);
            }
            else 
            {
                srvInfo->socket = socketConnect(srvInfo->port, srvInfo->hostname);
            }
            
            
            if (srvInfo->socket == -1)
            {
                perror("accept/connect");
                sleep(1);
            }
            else 
            {
                break;
            }
        }
        printf("Connected.\n");
        
        if (srvInfo->comInit)
        {
            ret = srvInfo->comInit(srvInfo);
        }
        
        if (srvInfo->comUpdate)
        {
            while (ret==0)
            {
                ret = srvInfo->comUpdate(srvInfo);
            }       
        }
        
        printf("%d, closing connection\n", ret);
        
        close(srvInfo->socket);
    }
    
    // never reaches this
    return NULL;
}


#pragma mark Socket Server Setup


int 
socketServerListen(int port)
{
    int one = 1;
    int local_socket;
    
    signal(SIGPIPE, SIG_IGN);
    
    if((local_socket = socket(AF_INET,SOCK_STREAM,0)) == -1) {
        perror("socket");
        return -1;
    }
    
    setsockopt(local_socket,SOL_SOCKET,SO_REUSEADDR, (const char *)&one,sizeof(int));
    
    struct sockaddr_in local_addr;
    memset(&local_addr,0,sizeof(local_addr));
    local_addr.sin_family = AF_INET;
    local_addr.sin_addr.s_addr = INADDR_ANY;
    local_addr.sin_port = htons(port);
    
    if(bind(local_socket,(struct sockaddr *)&local_addr,sizeof(local_addr)) == -1) {
        perror("bind");
        close(local_socket);
        return -1;
    }
    
    if(listen(local_socket,100)== -1) {
        perror("listen");
        close(local_socket);
        return -1;
    }	
    
    return local_socket;
}



int 
socketClose(int socket)
{
    close(socket); // close server-side socket
    
    return 0;
}



void 
socketStopThread(struct SocketServerInfo *srvInfo)
{
    int retVal;
    
    socketClose(srvInfo->socket);
    
    if (srvInfo->isServer) socketClose(srvInfo->server_socket);
    
    retVal = pthread_attr_destroy(&srvInfo->attr);
    
    pthread_mutex_destroy(&srvInfo->data_update_mutex);
    pthread_cond_destroy(&srvInfo->data_update_cv);
    
    free(srvInfo);
    
    assert(!retVal);
}



struct SocketServerInfo *
socketLaunchClientThread(const char *hostname, int server_port, int (*comInitF)(struct SocketServerInfo *), int (*comUpdateF)(struct SocketServerInfo *))
{
    // Create the thread using POSIX routines.
    int retVal;
    
    struct SocketServerInfo *srvInfo = (struct SocketServerInfo *) malloc(sizeof(struct SocketServerInfo));
    srvInfo->isServer = 0;
    srvInfo->hostname = (char *)hostname;
    srvInfo->port = server_port;
    srvInfo->comInit   = comInitF;
    srvInfo->comUpdate = comUpdateF;
    
    // set mutex
    pthread_mutex_init(&srvInfo->data_update_mutex, NULL);
    pthread_cond_init (&srvInfo->data_update_cv, NULL);
    
    retVal = pthread_attr_init(&srvInfo->attr);
    assert(!retVal);
    retVal = pthread_attr_setdetachstate(&srvInfo->attr, PTHREAD_CREATE_DETACHED);
    assert(!retVal);
    
    
    int threadError = pthread_create(&srvInfo->posixThreadID, &srvInfo->attr, &socketThreadMain, srvInfo);
    if (threadError != 0)
    {
        printf("Problem creating thread\n");
    }
    
    return srvInfo;
}



struct SocketServerInfo *
socketLaunchServerThread(int server_port, int (*comInitF)(struct SocketServerInfo *), int (*comUpdateF)(struct SocketServerInfo *))
{
    // Create the thread using POSIX routines.
    int retVal;
    
    struct SocketServerInfo *srvInfo = (struct SocketServerInfo *) malloc(sizeof(struct SocketServerInfo));
    srvInfo->isServer = 1;
    srvInfo->port = server_port;
    srvInfo->comInit   = comInitF;
    srvInfo->comUpdate = comUpdateF;
    
    srvInfo->server_socket = socketServerListen(srvInfo->port);
    assert(srvInfo->server_socket != -1);
    
    // set mutex
    pthread_mutex_init(&srvInfo->data_update_mutex, NULL);
    pthread_cond_init (&srvInfo->data_update_cv, NULL);
    
    retVal = pthread_attr_init(&srvInfo->attr);
    assert(!retVal);
    retVal = pthread_attr_setdetachstate(&srvInfo->attr, PTHREAD_CREATE_DETACHED);
    assert(!retVal);
    
    
    int threadError = pthread_create(&srvInfo->posixThreadID, &srvInfo->attr, &socketThreadMain, srvInfo);
    if (threadError != 0)
    {
        printf("Problem creating thread\n");
    }
    
    return srvInfo;
}



void
socketLockDataMutex(struct SocketServerInfo *srvInfo)
{
    pthread_mutex_lock(&srvInfo->data_update_mutex);
}



void 
socketUnlockDataMutex(struct SocketServerInfo *srvInfo)
{
    pthread_mutex_unlock(&srvInfo->data_update_mutex);
}



void
socketUpdate(struct SocketServerInfo *srvInfo)
{
    pthread_cond_signal(&srvInfo->data_update_cv);
}




// ------ Matlab-specific -----------------------------------------------------------

#include <cstring>
#include <vector>
#include <assert.h>

static const int varSize[] = 
{
    0,
    1, //sizeof(mxLogical),
    2, //sizeof(mxChar),
    1,1,2,2,4,4,8,8,4,8, sizeof(char *), sizeof(char *)
};


// this is only called if MatlabVar is constructed from scratch
void
compute_blocksize(MatlabVar &m, MatVarType matType, MatVarComplexity matCmplx, mwSize nDimV, mwSize *dimV)
{
    switch (matType)
    {
        case MAT_EMPTY:
            m.block_size = sizeof(MatVarType);
            break;
            
        case MAT_STRUCT:
            m.block_size  = sizeof(MatVarType);
            m.block_size += sizeof(int) + sizeof(int) * nDimV;
            
            m.block_size += sizeof(int);
            for (size_t k=0; k<m.fieldNames.size(); k++)
            {
                m.block_size += strlen(m.fieldNames[k])+1;
            }
            
            break;
            
        case MAT_CELL:
            m.block_size  = sizeof(MatVarType);
            m.block_size += sizeof(int) + sizeof(int) * nDimV;
            
            break;
            
        default:
            mwSize elements = 1;
            for(mwSize i=0; i<nDimV; i++) elements *= dimV[i];
            
            m.block_size  = sizeof(MatVarType);
            m.block_size += sizeof(MatVarComplexity);
            m.block_size += sizeof(int) + sizeof(int) * nDimV;
            m.block_size += ((matCmplx == MAT_COMPLEX) ? 2 : 1 ) * elements * varSize[matType];
            break;
    }
    
    
}



mwSize MatlabVar::getNumElements()
{
    mwSize elements = 1;
    for(mwSize i=0; i<*ndim; i++) elements *= dims[i];
    
    return elements;
}


void MatlabVar::deleteChildren()
{
    for (mwSize k=0; k<children.size(); k++)
    {
        delete children[k];
    }
    nchildren = 0; 
    
    children.clear(); 
}


// init Pointers from serialized block
void initPointers(MatlabVar &m, unsigned char *cdata)
{
    m.block_data = cdata;
    m.deleteChildren();
    m.fieldNames.clear();
    
    m.varType = (MatVarType *)cdata; cdata += sizeof(MatVarType);
    m.varComplexity = 0; // only number arrays have this set.
    
    switch (*m.varType)
    {
        case MAT_EMPTY:
            break;
            
        case MAT_STRUCT:
            m.ndim = (mwSize *)cdata; cdata += sizeof(mwSize);
            m.dims = (mwSize *)cdata; cdata += sizeof(mwSize) * (*m.ndim);
            
            m.nFields   = (int *)cdata; cdata += sizeof(int);
            m.nchildren = m.getNumElements() * (mwSize)(*m.nFields);
            
            for(int k=0; k<*m.nFields; k++)
            {
                m.fieldNames.push_back((char *)cdata);
                cdata += strlen((char *)cdata) + 1;
            }
            break;
            
        case MAT_CELL:
            m.ndim = (mwSize *)cdata; cdata += sizeof(mwSize);
            m.dims = (mwSize *)cdata; cdata += sizeof(mwSize) * (*m.ndim);
            
            m.nchildren = m.getNumElements();
            
            break;
            
        default:
            m.varComplexity = (MatVarComplexity *)cdata; cdata += sizeof(MatVarComplexity);
            
            m.ndim = (mwSize *)cdata; cdata += sizeof(mwSize);
            m.dims = (mwSize *)cdata; cdata += sizeof(mwSize) * (*m.ndim);
            
            mwSize elements = m.getNumElements();
            
            m.pr = (unsigned char *)cdata; cdata += varSize[*m.varType] * elements;
            if (*m.varComplexity == MAT_COMPLEX)
            {
                m.pi = (unsigned char *)cdata; cdata += varSize[*m.varType] * elements;
            }
            
            break;
            
    }
    
    m.block_size = (mwSize) (cdata - m.block_data);
    
}



// this is called when we construct a MatlabVar from scratch

void initPointers(MatlabVar &m, MatVarType matType, MatVarComplexity matCmplx, mwSize nDimV, mwSize *dimV)
{
    unsigned char *cdata = m.block_data;
    
    m.deleteChildren();
    // field names already have to be set!
    
    m.varType = (MatVarType *)cdata; *m.varType = matType; cdata += sizeof(MatVarType);
    
    switch (matType)
    {
        case MAT_EMPTY:
            break;
            
        case MAT_STRUCT:
            m.ndim = (mwSize *)cdata; *m.ndim = nDimV; cdata += sizeof(mwSize);
            m.dims = (mwSize *)cdata;                    cdata += sizeof(mwSize) * (*m.ndim);
            for (mwSize k=0; k<*m.ndim; k++)
            {
                m.dims[k] = dimV[k]; 
            }
            
            m.nFields = (int *)cdata; *m.nFields = (int) m.fieldNames.size(); cdata += sizeof(int);
            for (size_t k=0; k<*m.nFields; k++)
            {
                size_t len = strlen(m.fieldNames[k])+1;
                memcpy(cdata, m.fieldNames[k], len );
                m.fieldNames[k] = (char *)cdata; cdata += len;
            }
            
            m.nchildren = m.getNumElements() * (*m.nFields);
            
            break;
            
        case MAT_CELL:
            m.ndim = (mwSize *)cdata; *m.ndim = nDimV; cdata += sizeof(mwSize);
            m.dims = (mwSize *)cdata;                    cdata += sizeof(mwSize) * (*m.ndim);
            for (int k=0; k<*m.ndim; k++)
            {
                m.dims[k] = dimV[k]; 
            }
            m.nchildren = m.getNumElements();
            
            break;
            
        default:
            m.varComplexity = (MatVarComplexity *)cdata; *m.varComplexity = matCmplx; cdata += sizeof(MatVarComplexity);
            
            m.ndim = (mwSize *)cdata; *m.ndim = nDimV; cdata += sizeof(mwSize);
            m.dims = (mwSize *)cdata;                    cdata += sizeof(mwSize) * (*m.ndim); 
            
            for (int k=0; k<*m.ndim; k++)
            {
                m.dims[k] = dimV[k];
            }
            
            mwSize elements = m.getNumElements();
            
            m.pr = (unsigned char *)cdata; cdata += varSize[matType] * elements;
            memset((void *)m.pr, 0, varSize[matType] * elements);
            if (matCmplx == MAT_COMPLEX)
            {
                m.pi = (unsigned char *)cdata; cdata += varSize[matType] * elements;
                memset((void *)m.pi, 0, varSize[matType] * elements);
            }
            break;
    }
}



void
allocateMemory(MatlabVar &m, MatVarType matType, MatVarComplexity matCmplx, mwSize nDimV, mwSize *dimV)
{
    compute_blocksize(m, matType, matCmplx, nDimV, dimV);
    
    m.block_data = new unsigned char[ m.block_size ];
    m.block_reserved = m.block_size;
    
    initPointers(m, matType, matCmplx, nDimV, dimV);
}


void
MatlabVar::addChild(auto_ptr<MatlabVar> child_ap)
{
    MatlabVar *child = child_ap.release();
    
    assert(*varType == MAT_CELL || *varType == MAT_STRUCT);
    
    children.push_back(child);
}

// constructors for Arrays

MatlabVar::MatlabVar(MatVarType matType, MatVarComplexity matCmplx, mwSize nDimV, mwSize *dimV)
{
    assert(matType != MAT_CELL  && matType != MAT_STRUCT);
    
    allocateMemory(*this, matType, matCmplx, nDimV, dimV);
}

MatlabVar::MatlabVar(MatVarType matType, MatVarComplexity matCmplx, mwSize d1, mwSize d2)
{
    assert(matType != MAT_CELL  && matType != MAT_STRUCT);
    
    mwSize dimV[] = {d1, d2};
    mwSize nDimV = d2==0 ? 1 : 2;
    
    allocateMemory(*this, matType, matCmplx, nDimV, dimV);
}

// constructors for Structs

MatlabVar::MatlabVar(MatVarType matType, mwSize nF, const char **cFNames)
{
    assert(matType == MAT_STRUCT);
    for (mwSize k=0; k<nF; k++) fieldNames.push_back(cFNames[k]);
    
    mwSize dimV[] = {1, 1};
    
    allocateMemory(*this, MAT_STRUCT, MAT_REAL, 2, dimV);
}

/*
MatlabVar::MatlabVar(MatVarType matType, std::vector<std::string> cFNames)
{
    assert(matType == MAT_STRUCT);
    for (mwSize k=0; k<cFNames.size(); k++) fieldNames.push_back(cFNames[k].c_str());
    
    mwSize dimV[] = {1, 1};
    
    allocateMemory(*this, MAT_STRUCT, MAT_REAL, 2, dimV);
}
*/


MatlabVar::MatlabVar(MatVarType matType, mwSize nF, const char **cFNames, mwSize nDimV, mwSize *dimV )
{
    assert(matType == MAT_STRUCT);
    for (mwSize k=0; k<nF; k++) fieldNames.push_back(cFNames[k]);
    
    allocateMemory(*this, MAT_STRUCT, MAT_REAL, nDimV, dimV);
}

// constructors for Cells

MatlabVar::MatlabVar(MatVarType matType, mwSize d1, mwSize d2)
{
    assert(matType == MAT_CELL);
    
    mwSize dimV[] = {d1, d2};
    mwSize nDimV = d2==0 ? 1 : 2;
    
    allocateMemory(*this, MAT_CELL, MAT_REAL, nDimV, dimV);
    
}

MatlabVar::MatlabVar(MatVarType matType, mwSize nDimV, mwSize *dimV)
{
    assert(matType == MAT_CELL);
    
    allocateMemory(*this, MAT_CELL, MAT_REAL, nDimV, dimV);
    
}

bool 
MatlabVar::isValid()
{
    switch(*varType)
    {
        case MAT_CELL:
        case MAT_STRUCT:
            if (nchildren != children.size()) return false;
            for (size_t k=0; k<nchildren; k++)
            {
                if (!children[k]->isValid()) return false;
            }
            break;
            
        default:
            break;
    }
    
    return true;
}


void MatlabVar::deallocate()
{
    //    printf("Deallocating.\n");
    if (block_reserved) delete [] block_data; 
    block_reserved = block_size = 0; block_data = 0;
    
    deleteChildren();
}

MatlabVar::~MatlabVar()
{
    deallocate();
}

/*
MatlabVar &MatlabVar::operator=(const MatlabVar& v)
{
    deallocate();
    
    switch (*v.varType)
    {
        case MAT_EMPTY:
            allocateMemory(*this, *v.varType, MAT_REAL, 0, 0);
            break;
            
        case MAT_STRUCT:
            fieldNames = v.fieldNames;
            allocateMemory(*this, *v.varType, *v.varComplexity, *v.ndim, v.dims);
            break;
            
        case MAT_CELL:
            allocateMemory(*this, *v.varType, *v.varComplexity, *v.ndim, v.dims);
            break;
            
        default:
            allocateMemory(*this, *v.varType, *v.varComplexity, *v.ndim, v.dims);
            memcpy(
            break;
            
    }
    return *this;
}
*/


MatlabVar *
MatlabVar::duplicate()
{
    unsigned char *newblock = new unsigned char[ block_size ];
    memcpy(newblock, block_data, block_size);

    MatlabVar *mv = new MatlabVar(newblock);
    mv->block_reserved = block_size;
    
    return mv;
}

mwSize 
getTotalBlockSize(MatlabVar &m)
{
    mwSize bs = (mwSize) m.block_size;
    for (mwSize k=0; k<m.nchildren; k++)
    {
        bs += getTotalBlockSize(*m.children[k]);
    }
    
    return bs;
}


void 
deserialize(MatlabVar &m, unsigned char *cdata)
{
    if (cdata == 0) return;
    
    initPointers(m, (unsigned char *)cdata);
    
    cdata += m.block_size;
    
    for (mwSize k=0; k<m.nchildren; k++)
    {
        m.children.push_back(new MatlabVar(cdata));
        cdata += getTotalBlockSize(*m.children.back());
    }
}

// constructor from block data

MatlabVar::MatlabVar(unsigned char *cdata)
{
    deserialize(*this, cdata);
}



mwSize 
getNumSubBlocks(MatlabVar &m)
{
    mwSize nc = 1;                                  // this block...
    
    for (mwSize k=0; k<m.nchildren; k++)
    {
        nc += getNumSubBlocks(*m.children[k]);        // ... and all its children.
    }
    
    return nc;
}



SocketBlock **
getSubBlocks(MatlabVar &m, SocketBlock **blocks)
{
    *blocks = &m; blocks++;                       // this block ...
    
    for (mwSize k=0; k<m.nchildren; k++)
    {
        blocks = getSubBlocks(*m.children[k], blocks); // ... and all its children
    }
    
    return blocks;
}



int 
socketSend(int remote_socket, MatlabVar *m)
{
    mwSize allSubBlocks = getNumSubBlocks(*m);
    
    SocketBlock **blocks = new SocketBlock *[allSubBlocks];
    
    getSubBlocks(*m, blocks);
    
    int retVal = socketSendMultiBlock(remote_socket, allSubBlocks, blocks);
    
    delete [] blocks;
    
    return retVal;
}



int
socketRecv(int remote_socket, MatlabVar *m)
{
    int ret = socketReceiveBlock(remote_socket, *m);
    
    if (!ret) return ret;
    
    deserialize(*m, m->block_data);
    
    return ret;
}






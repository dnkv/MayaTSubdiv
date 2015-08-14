/*
 *  msockets.h
 *  SocketServer
 *
 *  Created by Denis Kovacs on 10/6/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef MSOCKETS_H
#define MSOCKETS_H


#include <string.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <sys/signal.h>
#include <sys/un.h>
#include <unistd.h>
#include <netinet/in.h>

#include <stdint.h>

struct SocketBlock
{
    size_t         block_size;     // active block size used
    size_t         block_reserved; // allocated block size
    unsigned char *block_data;
    
    SocketBlock(): block_size(0), block_reserved(0), block_data(0) {}
};



// MATLAB-specific --------------------------------

#include <vector>
#include <map>
#include <string>
#include <assert.h>

typedef enum {
	MAT_EMPTY = 0,
	MAT_LOGICAL = 1,
	MAT_CHAR = 2,
	MAT_UINT8 = 3,
	MAT_INT8 = 4,
	MAT_UINT16 = 5,
	MAT_INT16 = 6,
	MAT_UINT32 = 7,
	MAT_INT32 = 8,
	MAT_UINT64 = 9,
	MAT_INT64 = 10,
	MAT_FLOAT32 = 11,
	MAT_SINGLE = 11,
	MAT_FLOAT64 = 12,
	MAT_DOUBLE = 12,
	MAT_STRUCT = 13,
	MAT_CELL = 14,
	MAT_VAREND = 15
} MatVarType;

extern std::string MatVarStrings[];


typedef enum {
	MAT_REAL = 0,
	MAT_COMPLEX = 1
} MatVarComplexity;


typedef unsigned int mwSize;

struct MatlabVar: public SocketBlock
{
    MatlabVar() { block_data=0, block_size=0; block_reserved=0; nchildren = 0; }
    MatlabVar(unsigned char *cdata); // build from from serialized data
    
    MatlabVar(MatVarType matType, MatVarComplexity matCmplx, mwSize nDimV, mwSize *dimV);
    MatlabVar(MatVarType matType, MatVarComplexity matCmplx, mwSize d1, mwSize d2 = 0);
    
    MatlabVar(MatVarType matType, mwSize nF, const char **cFNames);
    MatlabVar(MatVarType matType, mwSize nF, const char **cFNames, mwSize nDimV, mwSize *dimV );
    
    MatlabVar(MatVarType matType, mwSize d1, mwSize d2=0);
    MatlabVar(MatVarType matType, mwSize nDimV, mwSize *dimV );
    
    void deallocate();
    void deleteChildren();
    ~MatlabVar();
    
    MatlabVar *duplicate();
    void addChild(std::auto_ptr<MatlabVar> mv);
    bool isValid();
    mwSize getNumElements(); // returns number of elements
    
    MatVarType       *varType;
    MatVarComplexity *varComplexity;
    mwSize           *ndim;
    mwSize           *dims;
    
    void             *pr;
    void             *pi;
    
    int                       *nFields;
    std::vector<const char *>  fieldNames;
    
    mwSize                     nchildren;
    std::vector<MatlabVar *>   children;
};


int // returns number of bytes sent
socketSend(int remote_socket, MatlabVar *v);

int // returns number of bytes received
socketRecv(int remote_socket, MatlabVar *v = 0);


// communication - specific

struct SocketServerInfo
{
    pthread_attr_t  attr;
    pthread_t       posixThreadID;
    
    pthread_mutex_t data_update_mutex;
    pthread_cond_t  data_update_cv;
    
    int             server_socket;
    int             socket;
    char           *hostname;
    
    char            isServer;
    int             port;
    
    int (*comInit)(struct SocketServerInfo *);
    int (*comUpdate)(struct SocketServerInfo *);
};



void 
socketDumpSocketBlock(struct SocketBlock d);


// basic I/O

int 
socketServerListen(int port);

int
socketServerAccept(int local_socket);

int
socketConnect(int port, char *hostname);

int 
socketClose(int socket);

// mutex

void
socketLockDataMutex(struct SocketServerInfo *srvInfo);

void 
socketUnlockDataMutex(struct SocketServerInfo *srvInfo);

void
socketUpdate(struct SocketServerInfo *srvInfo);

// communication thread

void 
socketWaitForUpdate(struct SocketServerInfo *srvInfo);

// main thread

struct SocketServerInfo *
socketLaunchClientThread(const char *hostname, int server_port, int (*comInitF)(struct SocketServerInfo *), int (*comUpdateF)(struct SocketServerInfo *));

struct SocketServerInfo *
socketLaunchServerThread(int server_port, int (*comInitF)(struct SocketServerInfo *), int (*comUpdateF)(struct SocketServerInfo *));

void 
socketStopThread(struct SocketServerInfo *srvInfo);

#endif
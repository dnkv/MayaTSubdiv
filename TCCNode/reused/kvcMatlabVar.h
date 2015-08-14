#ifndef KVCMATLABVAR_H
#define KVCMATLABVAR_H

#include "msockets.h"
#include "kvc.h"

bool create_MatlabVar_from_KVCObject(KVCObject *kvco);

bool init_KVCObject_from_MatlabVar(KVCObject *kvco, std::auto_ptr<MatlabVar> mv);
bool init_KVCObject_from_Socket(KVCObject *kvco, int remote_socket);
bool send_KVCObject_to_Socket(KVCObject *kvco, int remote_socket);

#endif
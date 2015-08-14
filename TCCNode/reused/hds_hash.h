//
//  hds_hash.h
//  TCCNode
//
//  Created by Denis Kovacs on 9/25/12.
//  Copyright 2012 __MyCompanyName__. All rights reserved.
//

#ifndef HDS_HASH 
#define HDS_HASH

// modified from http://stackoverflow.com/questions/7110301/

namespace std {
    template <class T>
    struct hash<std::pair<T, T> >
    {
        size_t
        operator()(std::pair<T, T> const& tt) const
        {
            size_t seed = 0;
            seed ^= hash<T>()(tt.first) + 0x9e3779b9 + (seed<<6) + (seed>>2);
            seed ^= hash<T>()(tt.second) + 0x9e3779b9 + (seed<<6) + (seed>>2);
            
            return seed;
        }
        
    };
}

#endif //HDS_HASH
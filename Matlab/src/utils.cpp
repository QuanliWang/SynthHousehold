#include "mex.h"
#include <vector>
#include <cmath>

#define HEAD 1
#define SPOUSE 2
#define CHILD 3
#define CHILDINLAW 4
#define PARENT 5
#define PARENTINLAW 6
#define SIBLING 7
#define SIBLINGINLAW 8
#define GRANDCHILD 9


//1 = head/householder, 2 = spouse, 3 = child, 4 = child-in-law, 5 = parent, 6 = parent-in- law, 7 = sibling, 8 = sibling-in-law,
//9 = grandchild,
//10 = other relatives,
//11 = partner, friend, visitor,
//12 = other non-relatives
inline bool IsHead(double relate, double age) {
    return (relate == HEAD && age >=17);
}

inline bool MoreThanOneHead(double *record, int hhsize) {
    int nhead = 0;
    for (int i = 1; i <= hhsize; i++) {
        if (record[2*hhsize+i]==HEAD) {
            nhead++;
        }
    }
    return (nhead >1);
}

inline int GetValidSpouse(double *record, int hhsize) {
    int nspouse = 0;
    int spouse = -1;
    for (int i = 1; i <= hhsize; i++) {
        if (record[2*hhsize+i]==SPOUSE) {
            nspouse++;
            spouse = i;
        }
    }
    if (nspouse > 1) {return 0;} //too many spouse
    if (nspouse == 0) { return -1;} //no spouse
    if (record[hhsize+spouse]<17) {return 0;} //spouse is under-age
    return spouse;
}

inline bool IsValidCouple(double *record, int hh_size, int spouse) {
    if (spouse ==0) { //bad spouse or too many spouses
        return false;
    } else { //valid spouse or no spouse
        if (spouse>0) {//the only spouse, so check sex, and age difference
            if (record[1] == record[spouse]) {return false;}
            if (std::abs(record[hh_size + 1] - record[hh_size + spouse]) > 52) {return false;}
        }
    }
    return true;
}

//return -1 if no child
//return the record index of the oldest child otherwise
inline int GetOldestChild(double *record, int hhsize) {
    double age = -1;
    int child = -1;  //no childen
    for (int i = 1; i <= hhsize; i++) {
        if (record[2*hhsize+i]==CHILD) {
            if (record[hhsize+i] > age) {
                age = record[hhsize+i];
                child = i;
            }
        }
    }
    return child;
}

inline bool IsValidChild(double *record, int hh_size, int child) {
    if (child>0) {//get a child, check age difference
        if (record[hh_size + 1] - record[hh_size + child] <12) {return false;}
    }
    return true;
}


//return -1 if no child
//return the record index of the oldest child otherwise
inline int GetOldestChildInLaw(double *record, int hhsize) {
    double age = -1;
    int child = -1;  //no childen
    for (int i = 1; i <= hhsize; i++) {
        if (record[2*hhsize+i]==CHILDINLAW) {
            if (record[hhsize+i] > age) {
                age = record[hhsize+i];
                child = i;
            }
        }
    }
    return child;
}

inline bool IsValidChildInLaw(double *record, int hh_size, int child) {
    if (child>0) {//get a child, check age difference
        if (record[hh_size + 1] - record[hh_size + child] <10) {return false;}
    }
    return true;
}

//return -1 if no parent
//return the record index of the youngest parent otherwise
inline int GetYoungestParent(double *record, int hhsize) {
    double age = 1000;
    int parent = -1;  //no childen
    for (int i = 1; i <= hhsize; i++) {
        if (record[2*hhsize+i]==PARENT) {
            if (record[hhsize+i] < age) {
                age = record[hhsize+i];
                parent = i;
            }
        }
    }
    return parent;
}

inline bool IsValidParent(double *record, int hh_size, int parent) {
    if (parent>0) {//get a child, check age difference
        if (record[hh_size + parent] -record[hh_size + 1] <13) {return false;}
    }
    return true;
}

inline int GetYoungestParentInLaw(double *record, int hhsize) {
    double age = 1000;
    int parent = -1;  //no childen
    for (int i = 1; i <= hhsize; i++) {
        if (record[2*hhsize+i]==PARENTINLAW) {
            if (record[hhsize+i] < age) {
                age = record[hhsize+i];
                parent = i;
            }
        }
    }
    return parent;
}

inline bool IsValidParentInLaw(double *record, int hh_size, int parent) {
    if (parent>0) {//get a child, check age difference
        if (record[hh_size + parent] -record[hh_size + 1] <9) {return false;}
    }
    return true;
}

inline bool IsValidSiblingOrSiblingInLaw(double *record, int hhsize) {
    for (int i = 1; i <= hhsize; i++) {
        if (record[2*hhsize+i]== SIBLING || record[2*hhsize+i] == SIBLINGINLAW) {
            if (std::abs(record[hhsize + i] - record[hhsize + 1]) >33) {return false;}
        }
    }
    return true;
}

inline bool IsValidGrandChild(double *record, int hhsize, int spouse) {
    for (int i = 1; i <= hhsize; i++) {
        if (record[2*hhsize+i]== GRANDCHILD) {
            if (record[hhsize + 1] < 33) {return false;} //too young to be grand parent for the HEAD
            if (spouse > 0) { //make sure the spouse(if any) is not too young
                if (record[hhsize + spouse] < 33) {return false;}
            }
            if (record[hhsize + 1] - record[hhsize + i] <30 ) {return false;}
        }
    }
    return true;
}


double isValid(double *datah, int hh_size) {
    
    //mexPrintf("check head 0\n");
    if (!IsHead(datah[2 * hh_size + 1], datah[hh_size + 1])) {return 0;}
    //mexPrintf("check head 1\n");
    if (MoreThanOneHead(datah,hh_size)) {return 0;}
    
    //mexPrintf("check head\n");
    int spouse = GetValidSpouse(datah,hh_size);
    if (!IsValidCouple(datah,hh_size,spouse)) {return 0;}
    
    //mexPrintf("check spouse\n");
    
    int oldestChild = GetOldestChild(datah,hh_size);
    if (!IsValidChild(datah,hh_size,oldestChild)) {return 0;}
    
    //mexPrintf("%d\n",oldestChild);
    int oldestChildInLaw = GetOldestChildInLaw(datah,hh_size);
    if (!IsValidChildInLaw(datah,hh_size,oldestChildInLaw)) {return 0;}
    
    //mexPrintf("2\n");
    int youngestParent = GetYoungestParent(datah,hh_size);
    if (!IsValidParent(datah,hh_size,youngestParent)) {return 0;}
    
    //mexPrintf("3\n");
    int youngestParentInLaw = GetYoungestParentInLaw(datah,hh_size);
    if (!IsValidParentInLaw(datah,hh_size,youngestParentInLaw)) {return 0;}
    
    //mexPrintf("4\n");
    if (!IsValidSiblingOrSiblingInLaw(datah,hh_size)) {return 0;}
    
    //mexPrintf("5\n");
    if (!IsValidGrandChild(datah,hh_size,spouse)) {return 0;}
    
    return 1;
    
}

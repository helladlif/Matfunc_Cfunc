//
// Created by lidon on 2023/12/13.
//

#include "findc.h"

///找到大于0的位置
MatrixXi find_loc(const MatrixXf &x){
    MatrixXi loc=MatrixXi ::Zero(x.rows(),1);
    int j=0;
    for(int i=0;i<x.rows();i++){
        if(x(i,0)>0){
            loc(j++,0)=i;
        }
    }
    return loc.topRows(j);
}
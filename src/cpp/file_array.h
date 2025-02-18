//
//  file_array.h
//  Template functions must defined where declared, cannot be in separate .cpp file
//
//  Created by Bruce Zhou on 10/23/18.
//

#ifndef file_array_h
#define file_array_h
#include <fstream>
#include <vector>
#include <string.h>
#include <iostream>


template <typename T> bool xy_txtfile(const std::string &ofilename,
                                      const std::vector<T> &time,
                                      const std::vector<T> &timevec,
                                      const int &nTime)
{
    std::ofstream ofile(ofilename.c_str());
    
    for(unsigned int i=0; i<nTime; i++)
    {
        ofile <<time[i]<<" "<<timevec[i]<< std::endl;
    }
    ofile.close();
    return true;
}

template <typename T> bool x_txtfile(const std::string &ofilename,
                                     const std::vector<T> &data)
{
    std::ofstream ofile(ofilename.c_str());
    
    for(unsigned int i=0; i<data.size(); i++)
    {
        ofile << data[i]<< std::endl;
    }
    ofile.close();
    return true;
}
#endif /* file_array_h */


#ifndef _ARRAYDOC_H_
#define _ARRAYDOC_H_

#include <vector>


namespace SXSI 
{

class ArrayDoc {
public:
    ArrayDoc(BlockArray *input)
        : data(input)
    {

    }
    ArrayDoc(FILE *fp)
        : data(0)
    {
        data = new BlockArray(fp);
    }

    ~ArrayDoc()
    {
        delete data;
    }

    void save(FILE *fp)
    {
        data->Save(fp);
    }
    
    inline uint access(uint i)
    {
        return (*data)[i];
    }
    inline std::vector<int> accessAll(uint i, uint j)
    {
        std::vector<int> res;
        res.reserve(j-i+1);

        for (; i <= j; ++i)
            res.push_back((*data)[i]);

        return res;
    }
    
    std::vector<int> access(uint i, uint j, uint min, uint max)
    {
        std::vector<int> res;
        res.reserve(j-i+1);

        for (; i <= j; ++i)
            if ((*data)[i] >= min && (*data)[i] <= max)
                res.push_back((*data)[i]);

        return res;
    }
    

    uint count(uint i, uint j, uint min, uint max)
    {
        uint c = 0;
        for (; i <= j; ++i)
            if ((*data)[i] >= min && (*data)[i] <= max)
                ++c;
        return c;
    }
    
    
private:
    BlockArray *data;
};

};

#endif

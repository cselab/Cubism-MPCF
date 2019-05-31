/*
 *  SerializerIO.h
 *  Cubism
 *
 *  Created by Diego Rossinelli on 5/28/08.
 *  Copyright 2008 CSE Lab, ETH Zurich. All rights reserved.
 *
 */
#pragma once
#include <fstream>
#include <vector>
#include <string>

#include "Common.h"

CUBISM_NAMESPACE_BEGIN

template<typename GridType, typename TStreamer>
class SerializerIO
{
public:

    // Typedefs
    typedef typename GridType::BlockType TBlock;
    typedef typename TBlock::ElementType TElement;

    // Virtual methods
    virtual void Write(GridType & inputGrid, std::string fileName)
    {
        std::ofstream output(fileName.c_str(),  ios::out);

        output << inputGrid;

        const std::vector<BlockInfo> vInfo = inputGrid.getBlocksInfo();
        for(std::vector<BlockInfo>::const_iterator it = vInfo.begin(); it!= vInfo.end(); ++it)
            ((TBlock*)(it->ptrBlock))->template Write<TStreamer>(output);
    }

    virtual void Read(GridType & inputGrid, std::string fileName)
    {
        std::ifstream input(fileName.c_str(), ios::in);

        input >> inputGrid;

        const std::vector<BlockInfo> vInfo = inputGrid.getBlocksInfo();
        for(std::vector<BlockInfo>::const_iterator it = vInfo.begin(); it!= vInfo.end(); ++it)
            ((TBlock*)(it->ptrBlock))->template Read<TStreamer>(input);
    }
};

CUBISM_NAMESPACE_END

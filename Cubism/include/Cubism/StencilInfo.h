/*
 *  StencilInfo.h
 *  Cubism
 *
 *  Created by Diego Rossinelli on 11/17/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */
#pragma once

#include <vector>
#include <cstdarg>
#include <cassert>
#include <iostream>
#include <cstdlib>

#include "Common.h"

CUBISM_NAMESPACE_BEGIN

struct StencilInfo
{
    int sx, sy, sz, ex, ey, ez;
    std::vector<int> selcomponents;

    bool tensorial;

    StencilInfo() {}

    // Deprecated, use std::vector constructor instead.
    StencilInfo(int _sx, int _sy, int _sz, int _ex, int _ey, int _ez, bool _tensorial, const int ncomponents, ...):
    sx(_sx), sy(_sy), sz(_sz), ex(_ex), ey(_ey), ez(_ez), selcomponents(), tensorial(_tensorial)
    {
        assert(ncomponents > 0);

        va_list ev_list;

        va_start(ev_list, ncomponents);
        for(int i=0; i<ncomponents; i++)
            selcomponents.push_back(va_arg(ev_list, int));
        va_end(ev_list);

        if (!isvalid())
        {
            std::cout << "Stencilinfo instance not valid. Aborting\n";
            abort();
        }
    }

  StencilInfo(int _sx, int _sy, int _sz, int _ex, int _ey, int _ez, bool _tensorial, const std::vector<int> &components) :
  sx(_sx), sy(_sy), sz(_sz), ex(_ex), ey(_ey), ez(_ez), selcomponents(components), tensorial(_tensorial)
  {
    assert(selcomponents.size() > 0);

        if (!isvalid())
        {
            std::cout << "Stencilinfo instance not valid. Aborting\n";
            abort();
        }
    }

    StencilInfo(const StencilInfo& c):
        sx(c.sx), sy(c.sy), sz(c.sz),
        ex(c.ex), ey(c.ey), ez(c.ez),
        selcomponents(c.selcomponents), tensorial(c.tensorial)
    {
    }

    std::vector<int> _all() const
    {
        int extra[] = {sx, sy, sz, ex, ey, ez, (int)tensorial};
        std::vector<int> all(selcomponents);
        all.insert(all.end(), extra, extra + sizeof(extra)/sizeof(int));

        return all;
    }

    bool operator<(StencilInfo s) const
    {
        std::vector<int> me = _all(), you = s._all();

        const int N = std::min(me.size(), you.size());

        for(int i=0; i<N; ++i)
            if (me[i] < you[i])
                return true;
            else
                if (me[i] > you[i])
                    return false;

        return me.size() < you.size();
    }

    bool isvalid() const
    {
        const bool not0 = selcomponents.size() == 0;
        const bool not1 = sx>0 || ex<=0 || sx>ex;
        const bool not2 = sy>0 || ey<=0 || sy>ey;
        const bool not3 = sz>0 || ez<=0 || sz>ez;

        return !(not0 || not1 || not2 || not3);
    }
};

CUBISM_NAMESPACE_END

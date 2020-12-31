/*
 * TdZdd: a Top-down/Breadth-first Decision Diagram Manipulation Framework
 * by Hiroaki Iwashita <iwashita@erato.ist.hokudai.ac.jp>
 * Copyright (c) 2014 ERATO MINATO Project
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 */

#pragma once

#include <cassert>
#include <cstring>
#include <stdint.h>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>

#include <tdzdd/DdSpec.hpp>

#include "Graph.hpp"

/*struct ComponentWeightSpecCount {
    int16_t comp; ///< uncolored edge component counter.

    ComponentWeightSpecCount()
            : comp(0) {
    }

    ComponentWeightSpecCount(int16_t uncoloredEdgeComponents)
            : comp(uncoloredEdgeComponents) {
    }

    size_t hash() const {
        return comp;
    }

    bool operator==(ComponentWeightSpecCount const& o) const {
        return comp == o.comp;
    }

    friend std::ostream& operator<<(std::ostream& os,
            ComponentWeightSpecCount const& o) {
        return os << o.comp;
    }
};*/

class ComponentWeightSpecMate {
public:
    typedef int32_t Offset;
    typedef uint32_t uOffset;

private:
    Offset hoc; ///< offset to head or vertex-weight.
    Offset nxt; ///< offset to next connected vertex.

public:
    ComponentWeightSpecMate(Offset hoc = 0)
            : hoc(hoc), nxt(0) {
    }

    bool operator==(ComponentWeightSpecMate const& o) const {
        return this == &o;
    }

    bool operator!=(ComponentWeightSpecMate const& o) const {
        return this != &o;
    }

    void clear() {
        hoc = 0;
        nxt = 0;
    }

    bool isHead() const {
        return hoc >= 0;
    }

    bool isTail() const {
        return nxt == 0;
    }

    bool isIsolated() const {
        return isHead() && isTail();
    }

    ComponentWeightSpecMate& head() {
        return isHead() ? *this : *(this + hoc);
    }

    ComponentWeightSpecMate const& head() const {
        return isHead() ? *this : *(this + hoc);
    }

    ComponentWeightSpecMate& next() {
        return *(this + nxt);
    }

    ComponentWeightSpecMate const& next() const {
        return *(this + nxt);
    }

    int getWeight() const {
        assert(head().hoc >= 0);
        return head().hoc;
    }

    void print() const {
        std::cerr << "(" << hoc << "," << nxt << ")";
    }

    void print(std::ostream& ost) const {
        ost << "(" << hoc << "," << nxt << ")";
    }

    void print(FILE* fp) const {
        fprintf(fp, "(%d,%d)", hoc, nxt);
    }

    void mergeLists(ComponentWeightSpecMate& o1, ComponentWeightSpecMate& o2,
                    ComponentWeightSpecMate* mate) {
        ComponentWeightSpecMate* p1 = &o1.head();
        ComponentWeightSpecMate* p2 = &o2.head();
        if (p1 == p2) return;
        if (p1 > p2) std::swap(p1, p2);

        p1->hoc += p2->hoc;

        for (ComponentWeightSpecMate* q = p2;; q += q->nxt) {
            q->hoc = p1 - q;
            if (q->nxt == 0) break;
        }

        ComponentWeightSpecMate* p = p1;
        ComponentWeightSpecMate* q = p2;

        while (true) {
            assert(p != q);
            ComponentWeightSpecMate* pp = p + p->nxt;
            assert(p <= pp && pp != q);

            while (p < pp && pp < q) {
                p = pp;
                pp += pp->nxt;
                assert(p <= pp && pp != q);
            }

            assert(p == pp || q < pp);
            p->nxt = q - p;
            if (p == pp) break;
            p = q, q = pp;
        }
    }

    void replaceHeadWith(ComponentWeightSpecMate& newHead, ComponentWeightSpecMate* mate) const {
        ComponentWeightSpecMate const* p = &head();
        ComponentWeightSpecMate* q = &newHead;

        Offset v = p->hoc;
        assert(v >= 0);

        q->hoc = v;

        while (q->nxt > 0) {
            q += q->nxt;
            q->hoc = &newHead - q;
        }
    }

    void removeFromList(ComponentWeightSpecMate const& o) {
        if (o.nxt == 0) {
            for (ComponentWeightSpecMate* p = this; p <= &o; ++p) {
                if (p + p->nxt == &o) p->nxt = 0;
            }
        }
        else {
            for (ComponentWeightSpecMate* p = this; p <= &o; ++p) {
                if (p + p->nxt == &o) p->nxt += o.nxt;
            }
        }
    }

    friend std::ostream& operator<<(std::ostream& os,
            ComponentWeightSpecMate const& o) {
        return os << "<" << o.hoc << "," << o.nxt << ">";
    }
};

class ComponentWeightSpec: public tdzdd::PodArrayDdSpec<ComponentWeightSpec,
        ComponentWeightSpecMate,2> {
    typedef ComponentWeightSpecMate Mate;

    Graph const& graph;
    int const m;
    int const n;
    int const lower;
    int const upper;
    int const mateSize;
    std::vector<Mate> initialMate;
    bool const noLoop;
    bool const lookahead;

    int takable(Mate const* mate, Graph::EdgeInfo const& e) const {
        Mate const& w1 = mate[e.v1 - e.v0];
        Mate const& w2 = mate[e.v2 - e.v0];

        // don't connect again
        //if (noLoop && w1.head() == w2.head()) return false;

        if (w1.head() != w2.head() && w1.getWeight() + w2.getWeight() > upper) {
            return false;
        }

        if (e.v1final && e.v2final) {
            if (w1.isIsolated() && w2.isIsolated()) { // new component leaves immediately
                if (w1.getWeight() + w2.getWeight() < lower) {
                    return false;
                }
            }
            else if (w1.isHead() && w2 == w1.next() && w2.isTail()) { // existing component leaves
                if (w1.getWeight() < lower) {
                    return false;
                }
            }
        }
        return true;
    }

    bool doTake(Mate* mate, Graph::EdgeInfo const& e) const {
        if (!takable(mate, e)) return false;

        mate[0].mergeLists(mate[e.v1 - e.v0], mate[e.v2 - e.v0], mate);
        assert(mate[0].getWeight() <= upper);
        return true;
    }

    bool doNotTake(Mate* mate, Graph::EdgeInfo const& e) const {
        Mate& w1 = mate[e.v1 - e.v0];
        Mate& w2 = mate[e.v2 - e.v0];

        if (e.v1final && w1.isIsolated()) {
            if (w1.getWeight() < lower) {
                return false;
            }
        }

        if (e.v2final && w2.isIsolated()) {
            if (w2.getWeight() < lower) {
                return false;
            }
        }

        if (e.v1final && e.v2final && w1.isHead() && w2 == w1.next()
                && w2.isTail()) { // existing component leaves
            if (w1.getWeight() < lower) {
                return false;
            }
        }
        return true;
    }

    void update(Mate* mate, Graph::EdgeInfo const& e,
            Graph::EdgeInfo const& ee) const {
        int const d = ee.v0 - e.v0;
        assert(d >= 0);
        Mate* p1 = &mate[e.v1 - e.v0];
        Mate* p2 = &mate[e.v2 - e.v0];
        Mate* pd = p1 + d;

        for (Mate* q = p1; q < pd; ++q) {
            Mate* qq = &q->next();
            if (qq >= pd) {
                q->replaceHeadWith(*qq, mate);
            }
        }

        if (e.v2final) {
            mate[0].removeFromList(*p2);
            p2->clear();
        }

        if (e.v1final) {
            mate[0].removeFromList(*p1);
            p1->clear();
        }

        if (d > 0) {
            std::memmove(p1, pd, (mateSize - d) * sizeof(*mate));
            for (int i = mateSize - d; i < mateSize; ++i) {
                p1[i] = initialMate[ee.v0 + i];
            }
        }
    }

public:
    // weight_list: [weight of vertex 1, weight of vertex 2,..., weight of vertex m]
    ComponentWeightSpec(Graph const& graph, const std::vector<uint32_t>& weight_list,
                        uint32_t lower, uint32_t upper,
            bool noLoop = false, bool lookahead = true)
            : graph(graph), m(graph.vertexSize()), n(graph.edgeSize()),
              lower(lower), upper(upper),
              mateSize(graph.maxFrontierSize()), initialMate(1 + m + mateSize),
              noLoop(noLoop), lookahead(lookahead) {
        this->setArraySize(mateSize);

        for (int v = 1; v <= m; ++v) {
            std::ostringstream oss;
            oss << v;
            int u = graph.getVertex(oss.str());
            initialMate[u] = Mate(ComponentWeightSpecMate::
                                  Offset(weight_list[v - 1]));
        }
    }

    int getRoot(Mate* mate) const {
        for (int v = 1; v <= m; ++v) {
            if (initialMate[v].getWeight() > upper) {
                return 0;
            }
        }

        int const v0 = graph.edgeInfo(0).v0;

        for (int i = 0; i < mateSize; ++i) {
            mate[i] = initialMate[v0 + i];
        }

        return n;
    }

    int getChild(Mate* mate, int level, int take) const {
        assert(1 <= level && level <= n);
        int i = n - level;
        Graph::EdgeInfo const* e = &graph.edgeInfo(i);

        //std::cerr << "level = " << level << ", take = " << take;
        //for (int j = 0; j < mateSize; ++j) {
        //    mate[j].print(std::cerr);
        //}
        //std::cerr << std::endl;

        if (take) {
            if (!doTake(mate, *e)) return 0;
        }
        else {
            if (!doNotTake(mate, *e)) return 0;
        }

        if (++i == n) return -1;

        Graph::EdgeInfo const* ee = &graph.edgeInfo(i);
        update(mate, *e, *ee);

        while (lookahead) {
            e = ee;
            if (takable(mate, *e)) break;
            if (!doNotTake(mate, *e)) return 0;

            if (++i == n) return -1;

            ee = &graph.edgeInfo(i);
            update(mate, *e, *ee);
        }

        assert(i < n);
        return n - i;
    }

    //size_t hashCode(Count const& count) const {
    //    return count.hash();
    //}
};

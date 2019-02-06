/*//////////////////////////////////////////////////////////////////////////////
//
//  IMPORTANT: READ BEFORE DOWNLOADING, COPYING, INSTALLING OR USING.
//
//  By downloading, copying, installing or using the software you agree to this
//  license. If you do not agree to this license, do not download, install,
//  copy or use the software.
//
// This file originates from the openFABMAP project:
// [http://code.google.com/p/openfabmap/] -or-
// [https://github.com/arrenglover/openfabmap]
//
// For published work which uses all or part of OpenFABMAP, please cite:
// [http://ieeexplore.ieee.org/xpl/articleDetails.jsp?arnumber=6224843]
//
// Original Algorithm by Mark Cummins and Paul Newman:
// [http://ijr.sagepub.com/content/27/6/647.short]
// [http://ieeexplore.ieee.org/xpl/articleDetails.jsp?arnumber=5613942]
// [http://ijr.sagepub.com/content/30/9/1100.abstract]
//
//                           License Agreement
//
// Copyright (C) 2012 Arren Glover [aj.glover@qut.edu.au] and
//                    Will Maddern [w.maddern@qut.edu.au], all rights reserved.
//
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//  * Redistribution's of source code must retain the above copyright notice,
//    this list of conditions and the following disclaimer.
//
//  * Redistribution's in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other Materials provided with the distribution.
//
//  * The name of the copyright holders may not be used to endorse or promote
//    products derived from this software without specific prior written
///   permission.
//
// This software is provided by the copyright holders and contributors "as is"
// and any express or implied warranties, including, but not limited to, the
// implied warranties of merchantability and fitness for a particular purpose
// are disclaimed. In no event shall the Intel Corporation or contributors be
// liable for any direct, indirect, incidental, special, exemplary, or
// consequential damages (including, but not limited to, procurement of
// substitute goods or services; loss of use, data, or profits; or business
// interruption) however caused and on any theory of liability, whether in
// contract, strict liability,or tort (including negligence or otherwise)
// arising in any way out of the use of this software, even if advised of the
// possibility of such damage.
//////////////////////////////////////////////////////////////////////////////*/

#include "inference.hpp"

namespace of2 {

int InferBase::pq(int q) {
    return (int)clTree.at<double>(0,q);
}

double InferBase::Pzq(int q, bool zq) {
    return (zq) ? clTree.at<double>(1,q) : 1 - clTree.at<double>(1,q);
}

double InferBase::PzqGzpq(int q, bool zq, bool zpq) {
    if (zpq) {
        return (zq) ? clTree.at<double>(2,q) : 1 - clTree.at<double>(2,q);
    } else {
        return (zq) ? clTree.at<double>(3,q) : 1 - clTree.at<double>(3,q);
    }
}

double InferBinary::PzqGeq(bool zq, bool eq) {
    if (eq) {
        return (zq) ? PzGe : 1 - PzGe;
    } else {
        return (zq) ? PzGNe : 1 - PzGNe;
    }
}

double InferBinary::PeqGLzq(int q, bool Lzq, bool eq) {
    double alpha, beta;
    alpha = PzqGeq(Lzq, true) * Pzq(q, true);
    beta = PzqGeq(Lzq, false) * Pzq(q, false);

    if (eq) {
        return alpha / (alpha + beta);
    } else {
        return 1 - (alpha / (alpha + beta));
    }
}

double InferBinary::PzqGL(int q, bool zq, bool /*zpq*/, bool Lzq,
                              const bool & newPlace /*= false*/)
{
    double p = (newPlace ? Pzq(q, false) : PeqGLzq(q, Lzq, false)) * PzqGeq(zq, false) +
            (newPlace ? Pzq(q, true) : PeqGLzq(q, Lzq, true)) * PzqGeq(zq, true);

    return p;
}

double InferBinary::PzqGzpqL(int q, bool zq, bool zpq, bool Lzq,
                                 const bool & newPlace /*= false*/) {
    double p;
    double alpha, beta;

    alpha = Pzq(q,  zq) * PzqGeq(!zq, false) * PzqGzpq(q, !zq, zpq);
    beta  = Pzq(q, !zq) * PzqGeq( zq, false) * PzqGzpq(q,  zq, zpq);
    p = (newPlace ? Pzq(q, false) : PeqGLzq(q, Lzq, false))
            * beta / (alpha + beta);

    alpha = Pzq(q,  zq) * PzqGeq(!zq, true) * PzqGzpq(q, !zq, zpq);
    beta  = Pzq(q, !zq) * PzqGeq( zq, true) * PzqGzpq(q,  zq, zpq);
    p += (newPlace ? Pzq(q, true) : PeqGLzq(q, Lzq, true))
            * beta / (alpha + beta);

    return p;
}

} // namespace of2

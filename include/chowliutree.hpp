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
//    and/or other materials provided with the distribution.
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

#ifndef CHOWLIUTREE_H_
#define CHOWLIUTREE_H_

#include <opencv2/core/core.hpp>

#include <vector>
#include <list>

namespace of2 {

///
/// \brief A Chow-Liu tree implementation designed for FAB-MAP.
///
/// A Chow-Liu tree is required by FAB-MAP. The Chow-Liu tree provides an
/// estimate of the full distribution of visual words using a minimum spanning
/// tree. The tree is generated from training data.
///
class CV_EXPORTS ChowLiuTree {
public:
    ChowLiuTree();
    virtual ~ChowLiuTree();

    //@{
    ///
    /// \brief You add data to the chow-liu tree before calling make.
    /// \param imgDescriptor A \#imgs x \#words bag of words descriptor.
    ///
    void add(const cv::Mat& imgDescriptor);
    ///
    /// \brief You add data to the chow-liu tree before calling make.
    /// \param imgDescriptors A vector of \#imgs x \#words bag of words descriptors.
    ///
    void add(const std::vector<cv::Mat>& imgDescriptors);
    //@}

    const std::vector<cv::Mat>& getImgDescriptors() const;

    ///
    /// \brief Builds the Chow Liu tree from the descriptors that have been added.
    /// \param infoThreshold Ignores word pairs whose mutual information is below this threshold.
    /// \return A Mat containing the 4 x |v| Chow Liu tree,
    /// where (0,q) is parent (p) index, (1,q) is P(q), (2,q) is P(q|p), (3,q) is P(q|~p)
    ///
    cv::Mat make(double infoThreshold = 0.0);

private:
    std::vector<cv::Mat> imgDescriptors;
    cv::Mat mergedImgDescriptors;

    typedef struct info {
        float score;
        short word1;
        short word2;
    } info;

    //probabilities extracted from mergedImgDescriptors
    double P(int a, bool za);
    double JP(int a, bool za, int b, bool zb); //a & b
    double CP(int a, bool za, int b, bool zb); // a | b

    //calculating mutual inforMation of all edges
    void createBaseEdges(std::list<info>& edges, double infoThreshold);
    double calcMutInfo(int word1, int word2);
    static bool sortInfoScores(const info& first, const info& second);

    //selecting minimum spanning egdges with maximum inforMation
    bool reduceEdgesToMinSpan(std::list<info>& edges);

    //building the tree sctructure
    cv::Mat buildTree(int root_word, std::list<info> &edges);
    void recAddToTree(cv::Mat &cltree, int q, int pq,
                      std::list<info> &remaining_edges);
    std::vector<int> extractChildren(std::list<info> &remaining_edges, int q);

};

} // namespace of2

#endif /* CHOWLIUTREE_H_ */

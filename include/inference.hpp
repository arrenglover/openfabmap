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

#ifndef INFERENCE_H_
#define INFERENCE_H_

#include <opencv2/core/core.hpp>

namespace of2 {

///
/// \brief This class contains some base functionality for FabMap inference
///
class CV_EXPORTS InferBase
{
public:
    InferBase (cv::Ptr<cv::Mat> _clTreePtr,
               const double & _PzGe,
               const double & _PzGNe,
               const bool & _naiveBayes)
        : PzGe(_PzGe), PzGNe(_PzGNe),
          naiveBayes(_naiveBayes),
          clTreePtr(_clTreePtr), clTree(*clTreePtr)
    {
        //check for a valid Chow-Liu tree
        CV_Assert(clTree.type() == CV_64FC1);
        cv::checkRange(clTree.row(0), false, NULL, 0, clTree.cols);
        cv::checkRange(clTree.row(1), false, NULL, DBL_MIN, 1);
        cv::checkRange(clTree.row(2), false, NULL, DBL_MIN, 1);
        cv::checkRange(clTree.row(3), false, NULL, DBL_MIN, 1);
    }

    //@{Chow-Liu Tree
    ///
    /// \brief Probability of a word's measurement being true (from training).
    /// \param q The index of the word.
    /// \return The probability.
    ///
    int pq(int q);
    ///
    /// \brief The probability of a particular word measurement (from training).
    /// \param q The index of the word.
    /// \param zq The measurement value.
    /// \return The probability.
    ///
    double Pzq(int q, bool zq);
    ///
    /// \brief The probability of a particular word measurement given the CLT parent measurement.
    /// \param q The index of the word.
    /// \param zq The measurement value.
    /// \param zpq The CLT parent word's measurement value.
    /// \return The probability.
    ///
    double PzqGzpq(int q, bool zq, bool zpq);
    //@}

    ///
    /// \brief Get the Chow-Liu tree
    /// \return A const reference to the Chow-Liu tree
    ///
    const cv::Mat & getClTree()
    {
        return clTree;
    }

    ///
    /// \brief Get the vocabulary size as specified by the size of the Chow-Liu Tree
    /// \return The vocabulary size
    ///
    const int & vocabSize()
    {
        return clTree.cols;
    }

protected:

    //@{
    double PzGe;    ///< Probability of measurement given the generator is present
    double PzGNe;   ///< Probability of measurement given the generator is absent
    //@}

    /// Are we using naiveBayes?
    bool naiveBayes;
    /// A pointer to the Chow-Liu tree stored in the FabMap class
    cv::Ptr<cv::Mat> clTreePtr;
    /// A reference to the pointer, for legacy access (remove this eventually)
    cv::Mat & clTree;
};

///
/// \brief This class implements the binary inference for FabMap.
///
class CV_EXPORTS InferBinary : public InferBase
{
public:
    InferBinary(cv::Ptr<cv::Mat> _clTree,
                const double & _PzGe,
                const double & _PzGNe,
                const bool & _naiveBayes)
        : InferBase(_clTree, _PzGe, _PzGNe, _naiveBayes)
    {
    }

    //@{ FAB-MAP Core
    ///
    /// \brief Calculates the measurement probability given the generator state.
    /// \param zq If the word was seen in the image.
    /// \param eq If the generator is present in the scene.
    /// \return The probability of the observation given the generator.
    ///
    double PzqGeq(bool zq, bool eq);
    ///
    /// \brief Calculates the generator probability given the past measurement (singular) at the location.
    /// TODO: Replace this with a saved state, can update with multiple measurements.
    /// \param q The index of the word in the vocabulary.
    /// \param Lzq If the word previously observed at this location.
    /// \param eq If the generator is present in the scene.
    /// \return The probability that the generator is present at the location
    ///
    double PeqGLzq(int q, bool Lzq, bool eq);
    ///
    /// \brief Returns the generator probability given the past measurements at the location.
    /// \param q The index of the word in the vocabulary.
    /// \param Lm The index of the location.
    /// \param eq If the generator is present in the scene.
    /// \return The probability that the generator is present at the location
    ///
    double PeqGL(const int & q, const int & Lm, bool eq);
    ///
    /// \brief Calculates the Naive-Bayes measurement probability given past measurement.
    /// The parent measurement is unused.
    /// \param q The index of the word in the vocabulary.
    /// \param zq If the word was present in the image.
    /// \param zpq If the parent word in the Chow-Liu tree was present in the image (unused).
    /// \param Lzq If the word previously observed at this location.
    /// \param newPlace If we'd like the mean-field probability of being in a new place.
    /// \return The measurement probability given the past measurement.
    ///
    double PzqGL(int q, bool zq, bool zpq, bool Lzq,
                         const bool & newPlace = false);
    ///
    /// \brief Calculates the measurement probability given past and parent measurements.
    /// \param q The index of the word in the vocabulary.
    /// \param zq If the word was present in the image.
    /// \param zpq If the parent word in the Chow-Liu tree was present in the image.
    /// \param Lzq If the word previously observed at this location.
    /// \param newPlace If we'd like the mean-field probability of being in a new place.
    /// \return The measurement probability given past and parent measurements.
    ///
    double PzqGzpqL(int q, bool zq, bool zpq, bool Lzq,
                            const bool & newPlace = false);
    ///
    /// \brief Pointer to the function that will calculate the measurement probability (Naive-Bayes or Chow-Liu Tree).
    /// \param q The index of the word in the vocabulary.
    /// \param zq If the word was present in the image.
    /// \param zpq If the parent word in the Chow-Liu tree was present in the image.
    /// \param Lzq If the word previously observed at this location.
    /// \return The measurement probability given past (and parent measurements for CLT).
    ///
    double PzGL(int q, bool zq, bool zpq, bool Lzq,
                              const bool & newPlace)
    {
        return naiveBayes ? PzqGL(q, zq, zpq, Lzq, newPlace)
                          : PzqGzpqL(q, zq, zpq, Lzq, newPlace);
    }
    //@}
};

} // namespace of2

#endif /* INFERENCE_H_ */

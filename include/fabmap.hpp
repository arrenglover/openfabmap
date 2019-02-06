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

#ifndef FABMAP_H_
#define FABMAP_H_

#include "inference.hpp"

#include <opencv2/core/core.hpp>

#include <vector>
#include <list>
#include <map>
#include <set>

#ifdef __GNUC__
#define DEPRECATED(func) func __attribute__ ((deprecated))
#elif defined(_MSC_VER)
#define DEPRECATED(func) __declspec(deprecated) func
#else
#pragma message("WARNING: You need to implement DEPRECATED for this compiler")
#define DEPRECATED(func) func
#endif

namespace of2 {

///
/// \brief The return data format for FabMap comparisons and localizations.
///
struct CV_EXPORTS IMatch {

    /// Default constructor.
    IMatch() :
        queryIdx(-1), imgIdx(-1),
        groupSize(0),
        likelihood(-DBL_MAX), match(-DBL_MAX)
    {
    }
    ///
    /// \brief Initialization constructor
    /// \param _queryIdx Sets IMatch::queryIdx.
    /// \param _imgIdx Sets IMatch::imgIdx.
    /// \param _likelihood Sets IMatch::likelihood.
    /// \param _match Sets IMatch::match.
    /// \param _groupSize Sets IMatch::groupSize.
    ///
    IMatch(int _queryIdx, int _imgIdx,
           double _likelihood, double _match,
           unsigned int _groupSize = 0) :
        queryIdx(_queryIdx), imgIdx(_imgIdx),
        groupSize(_groupSize),
        likelihood(_likelihood), match(_match)
    {
    }

    int queryIdx;    ///< Query descriptor index (descriptor being compared).
    int imgIdx;      ///< Test descriptor index (reference descriptor being compared to).
    int groupSize;   ///< Size of image groups being compared (for PiraMap)
    int placeId;     ///< The index of the location (L) for the query after loop closure.

    double likelihood;  ///< The likelihood of the descriptors coming from the same location.
    double match;       ///< The normalized probability that the descriptors come from the same location

    /**
     * @brief IMatch operator < for match probabilities.
     * @param m The RHS IMatch object
     * @return If match probability for LHS < RHS
     */
    bool operator<(const IMatch& m) const {
        return match < m.match;
    }

};



///
/// \brief This class defines the common functionality for the FabMap derivatives.
///
class CV_EXPORTS FabMap {
public:

    /// FabMap flag options
    enum {
        MEAN_FIELD = 1,
        SAMPLED = 2,
        NAIVE_BAYES = 4,
        CHOW_LIU = 8,
        MOTION_MODEL = 16
    };

    ///
    /// \brief Base FabMap constructor
    /// \param clTree Chow Liu tree from training
    /// \param PzGe Measurement probability given the generator is present
    /// \param PzGNe Measurement probability given the generator is absent
    /// \param flags Flag for options
    /// \param numSamples Number of samples to use for new place sampling
    ///
    FabMap(const cv::Mat& clTree, double PzGe, double PzGNe, int flags,
           int numSamples = 0);
    virtual ~FabMap();

    //@{
    /// Method to add training data for sampling method.
    virtual void addTraining(const cv::Mat& queryImgDescriptor);
    virtual void addTraining(const std::vector<cv::Mat>& queryImgDescriptors);
    //@}

    //@{
    /// Method to add to the test data (the map).
    virtual void add(const cv::Mat& queryImgDescriptor);
    virtual void add(const std::vector<cv::Mat>& queryImgDescriptors);
    //@}

    //@{
    /// Access the descriptors with a read-only reference
    const std::vector<cv::Mat>& getTrainingImgDescriptors() const;
    const std::vector<cv::Mat>& getTestImgDescriptors() const;
    //@}

    //@{ Image comparisons
    ///
    /// \brief FabMap image comparison (not full localization).
    /// \param queryImgDescriptor The descriptor for the query image
    /// \param testImgDescriptors The descriptors for the test images
    /// \param[out] matches Contains the match probabilities for the comparisons
    /// \param mask Not currently implemented?
    ///
    void compare(const cv::Mat& queryImgDescriptor,
                 const cv::Mat& testImgDescriptors, std::vector<IMatch>& matches,
                 const cv::Mat& mask = cv::Mat());
    ///
    /// \brief FabMap image comparison (not full localization).
    /// \param queryImgDescriptor The descriptor for the query image
    /// \param testImgDescriptors The descriptors for the test images
    /// \param[out] matches Contains the match probabilities for the comparisons
    /// \param mask Not currently implemented?
    ///
    void compare(const cv::Mat& queryImgDescriptor,
                 const std::vector<cv::Mat>& testImgDescriptors,
                 std::vector<IMatch>& matches, const cv::Mat& mask = cv::Mat());
    ///
    /// \brief FabMap image comparison (not full localization).
    /// \param queryImgDescriptors The descriptors for the query images
    /// \param testImgDescriptors The descriptors for the test images
    /// \param[out] matches Contains the match probabilities for the comparisons
    /// \param mask Not currently implemented?
    ///
    void compare(const std::vector<cv::Mat>& queryImgDescriptors,
                 const std::vector<cv::Mat>& testImgDescriptors,
                 std::vector<IMatch>& matches, const cv::Mat& mask = cv::Mat());
    //@}

    //@{ Localization against map
    ///
    /// \brief FabMap localization against the map (deprecated, see FabMap::localize).
    /// \param queryImgDescriptors The descriptors for the query images
    /// \param[out] matches Contains the match probabilities for the comparisons
    /// \param addQuery Add the query to the test descriptors when done
    /// \param mask Not currently implemented?
    /// \deprecated This function has been deprecated, use FabMap::localize instead
    ///
    DEPRECATED (void compare(const std::vector<cv::Mat>& queryImgDescriptors, std::vector<
                             IMatch>& matches, bool addQuery = false, const cv::Mat& mask =
            cv::Mat()));
    ///
    /// \brief FabMap localization against the map.
    /// \param queryImgDescriptors The descriptors for the query images
    /// \param[out] matches Contains the match probabilities for the comparisons
    /// \param addQuery Add the query to the test descriptors when done
    /// \param mask Not currently implemented?
    /// TODO: Add addThreshold that replaces addQuery, and addIfLatest that doesn't add to the latest place.
    ///
    virtual void localize(const std::vector<cv::Mat>& queryImgDescriptors, std::vector<
                          IMatch>& matches, bool addQuery = false, const cv::Mat& mask =
            cv::Mat());
    ///
    /// \brief FabMap localization against the map (deprecated, see FabMap::localize).
    /// \param queryImgDescriptor The descriptor for the query image
    /// \param[out] matches Contains the match probabilities for the comparisons
    /// \param addQuery Add the query to the test descriptors when done
    /// \param mask Not currently implemented?
    /// \deprecated This function has been deprecated, use FabMap::localize instead
    ///
    DEPRECATED (void compare(const cv::Mat& queryImgDescriptor,
                             std::vector<IMatch>& matches, bool addQuery = false,
                             const cv::Mat& mask = cv::Mat()));
    ///
    /// \brief FabMap localization against the map.
    /// \param queryImgDescriptor The descriptor for the query image
    /// \param[out] matches Contains the match probabilities for the comparisons
    /// \param addQuery Add the query to the test descriptors when done
    /// \param mask Not currently implemented?
    /// TODO: Add addThreshold that replaces addQuery, and addIfLatest that doesn't add to the latest place.
    ///
    virtual void localize(const cv::Mat& queryImgDescriptor,
                          std::vector<IMatch>& matches, bool addQuery = false,
                          const cv::Mat& mask = cv::Mat());
    //@}

protected:

    //@{ Base Image descriptor operations
    ///
    /// \brief Compares a query to test descriptors.
    /// It calculates the likelihood that the two descriptors came from the same location.
    /// \param queryImgDescriptor The descriptor for the query image
    /// \param queryIndex The index of the query image (for query index in matches).
    /// \param testImgDescriptors The image descriptors to compare against.
    /// \param[out] matches The result of the descriptor comparisons in IMatch::likelihood.
    ///
    void compareImgDescriptor(const cv::Mat& queryImgDescriptor,
                              int queryIndex, const std::vector<cv::Mat>& testImgDescriptors,
                              std::vector<IMatch>& matches);
    ///
    /// \brief Adds an image descriptor to the test descriptors (the map).
    /// Does no comparisons to the test data.
    /// \param queryImgDescriptor The query image descriptor that will be added.
    ///
    void addImgDescriptor(const cv::Mat& queryImgDescriptor);
    //@}

    //@{ Likelihood functions, these change with each type of FabMap
    ///
    /// \brief The getLikelihoods method is overwritten for each different FabMap
    /// \param queryImgDescriptor The descriptor to be compared.
    /// \param testImgDescriptors The descriptors to compare against.
    /// \param[out] matches Result of the comparison of query to test.
    ///
    virtual void getLikelihoods(const cv::Mat& queryImgDescriptor,
                                const std::vector<cv::Mat>& testImgDescriptors,
                                std::vector<IMatch>& matches) = 0;
    ///
    /// \brief Calculates the likelihood the query comes from a new place.
    /// \param queryImgDescriptor The descriptor from the query image.
    /// \return The log-likelihood the query descriptor came from a new place.
    ///
    virtual double getNewPlaceLikelihood(const cv::Mat& queryImgDescriptor);
    //@}

    ///
    /// \brief Turns measurement likelihoods into location probabilities.
    /// Also applies the motion model if in use (specified in the options flag).
    /// \param[in,out] matches Contains the input likelihoods, and output probabilities.
    ///
    void normaliseDistribution(std::vector<IMatch>& matches);

    //@{ Data
    /// Image descriptors seen in training (for sampled new location probability)
    std::vector<cv::Mat> trainingImgDescriptors;
    /// Image descriptors seen so far
    std::vector<cv::Mat> testImgDescriptors;
    /// Prior match probabilities for motion model p(L^k|Z^k-1)
    std::vector<IMatch> priormatches;
    /// Generator states for all locations
    std::vector<cv::Mat> peGL;
    //@}

    //@{ Parameters
    double Pnew;    ///< Prior probability of entering a new location (motion model)

    double mBias;   ///< Forward motion bias, 1 all forward, 0 all backward (motion model)
    double sFactor; ///< Smoothing factor for matches, applied at the end (see papers)

    int flags;      ///< Flags for different modes
    int numSamples; ///< Number of samples to use for sampled new location probability
    //@}

    //@{
    /// Chow Liu Tree
    cv::Ptr<cv::Mat> clTree;
    /// Inference object (uses clTree)
    cv::Ptr<InferBinary> infer;
    //@}

};

///
/// \brief The original FAB-MAP algorithm, developed based on:
/// http://ijr.sagepub.com/content/27/6/647.short.
///
/// See the FabMap base class for more inforMation.
///
/// Note: It does not currently associate multiple measurements with a location (apply
/// the loop closure). This is a 'todo'.
///
class CV_EXPORTS FabMap1: public FabMap {
public:
    FabMap1(const cv::Mat& clTree, double PzGe, double PzGNe, int flags,
            int numSamples = 0);
    virtual ~FabMap1();
protected:

    ///
    /// \brief FabMap1 implementation of likelihood comparison
    /// \param queryImgDescriptor The query descriptor to be compared.
    /// \param testImgDescriptors The reference descriptors to compare against.
    /// \param[out] matches The results of the comparisons in IMatch::likelihood.
    ///
    void getLikelihoods(const cv::Mat& queryImgDescriptor,
                        const std::vector<cv::Mat>& testImgDescriptors,
                        std::vector<IMatch>& matches);
};

///
/// \brief A computationally faster Look-Up-Table version of the original FAB-MAP algorithm.
///
/// See the FabMap base class for more inforMation.
///
/// A look-up-table is used to precompute many of the reoccuring calculations.
///
class CV_EXPORTS FabMapLUT: public FabMap {
public:
    FabMapLUT(const cv::Mat& clTree, double PzGe, double PzGNe,
              int flags, int numSamples = 0, int precision = 6);
    virtual ~FabMapLUT();
protected:

    //FabMap look-up-table implementation of the likelihood comparison
    void getLikelihoods(const cv::Mat& queryImgDescriptor, const std::vector<
                        cv::Mat>& testImgDescriptors, std::vector<IMatch>& matches);

    //precomputed data
    int (*table)[8];

    //data precision
    int precision;
};

///
/// \brief The Accelerated FAB-MAP algorithm, developed based on:
/// http://ieeexplore.ieee.org/xpl/articleDetails.jsp?arnumber=5613942.
///
/// See the FabMap base class for more inforMation.
///
class CV_EXPORTS FabMapFBO: public FabMap {
public:
    FabMapFBO(const cv::Mat& clTree, double PzGe, double PzGNe, int flags,
              int numSamples = 0, double rejectionThreshold = 1e-8, double PsGd =
            1e-8, int bisectionStart = 512, int bisectionIts = 9);
    virtual ~FabMapFBO();

protected:

    //FabMap Fast Bail-out implementation of the likelihood comparison
    void getLikelihoods(const cv::Mat& queryImgDescriptor, const std::vector<
                        cv::Mat>& testImgDescriptors, std::vector<IMatch>& matches);

    ///
    /// \brief Stucture used to determine word comparison order
    ///
    struct WordStats {
        WordStats() :
            q(0), info(0), V(0), M(0) {
        }

        WordStats(int _q, double _info) :
            q(_q), info(_info), V(0), M(0) {
        }

        int q;
        double info;
        mutable double V;
        mutable double M;

        bool operator<(const WordStats& w) const {
            return info < w.info;
        }

    };

    //private fast bail-out necessary functions
    void setWordStatistics(const cv::Mat& queryImgDescriptor, std::multiset<WordStats>& wordData);
    double limitbisection(double v, double m);
    double bennettInequality(double v, double m, double delta);
    static bool compInfo(const WordStats& first, const WordStats& second);

    //parameters
    double PsGd;
    double rejectionThreshold;
    int bisectionStart;
    int bisectionIts;
};

///
/// \brief The FAB-MAP2.0 algorithm, developed based on:
/// http://ijr.sagepub.com/content/30/9/1100.abstract
///
/// See the FabMap base class for more inforMation.
///
class CV_EXPORTS FabMap2: public FabMap {
public:

    FabMap2(const cv::Mat& clTree, double PzGe, double PzGNe, int flags);
    virtual ~FabMap2();

    //FabMap2 builds the inverted index and requires an additional training/test
    //add function
    void addTraining(const cv::Mat& queryImgDescriptors) {
        FabMap::addTraining(queryImgDescriptors);
    }
    void addTraining(const std::vector<cv::Mat>& queryImgDescriptors);

    void add(const cv::Mat& queryImgDescriptors) {
        FabMap::add(queryImgDescriptors);
    }
    void add(const std::vector<cv::Mat>& queryImgDescriptors);

protected:

    //FabMap2 implementation of the likelihood comparison
    void getLikelihoods(const cv::Mat& queryImgDescriptor, const std::vector<
                        cv::Mat>& testImgDescriptors, std::vector<IMatch>& matches);
    double getNewPlaceLikelihood(const cv::Mat& queryImgDescriptor);

    //the likelihood function using the inverted index
    void getIndexLikelihoods(const cv::Mat& queryImgDescriptor, std::vector<
                             double>& defaults, std::map<int, std::vector<int> >& invertedMap,
                             std::vector<IMatch>& matches);
    void addToIndex(const cv::Mat& queryImgDescriptor,
                    std::vector<double>& defaults,
                    std::map<int, std::vector<int> >& invertedMap);

    //data
    std::vector<double> d1, d2, d3, d4;
    std::vector<std::vector<int> > children;

    // TODO: inverted map a std::vector?

    std::vector<double> trainingDefaults;
    std::map<int, std::vector<int> > trainingInvertedMap;

    std::vector<double> testDefaults;
    std::map<int, std::vector<int> > testInvertedMap;

};

} // namespace of2

#endif /* FABMAP_H_ */

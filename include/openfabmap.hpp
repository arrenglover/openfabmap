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

#ifndef OPENFABMAP_H_
#define OPENFABMAP_H_


#include <opencv2/opencv.hpp>
//#include "opencv2/core/core.hpp"
//#include "opencv2/features2d/features2d.hpp"

#include <vector>
#include <list>
#include <map>
#include <set>
#include <valarray>

#ifdef __GNUC__
#define DEPRECATED(func) func __attribute__ ((deprecated))
#elif defined(_MSC_VER)
#define DEPRECATED(func) __declspec(deprecated) func
#else
#pragma message("WARNING: You need to implement DEPRECATED for this compiler")
#define DEPRECATED(func) func
#endif

using cv::Mat;

namespace of2 {

using std::list;
using std::map;
using std::multiset;

///
/// \brief The return data format for FabMap comparisons and localizations.
///
struct CV_EXPORTS IMatch {

    /// Default constructor.
    IMatch() :
        queryIdx(-1), imgIdx(-1),
        likelihood(-DBL_MAX), match(-DBL_MAX),
        groupSize(0)
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
        likelihood(_likelihood), match(_match),
        groupSize(_groupSize)
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
/// \brief Stores a location in FabMap's topological map.
///
struct CV_EXPORTS Place
{
public:
    typedef std::pair<Place*,Place*> placePtrPair_t;
    /// The generator state p(e)
    Mat generators;
    /// The indices for the observation descriptors at this location
    std::list<int> descriptorIdx;
    ///
    /// \brief Get a reference to the previous and next places
    /// \return The previous and next places
    ///
    const std::list<placePtrPair_t> & getPlacesPrevNext()
    {
        return placesPrevNext;
    }
    ///
    /// \brief Get a const reference to the other linked places
    /// \return The other linked places.
    ///
    const std::list<Place*> & getPlacesConnected()
    {
        return placesConnected;
    }

private:
    /// \brief Pointers to other locations connected to this one in sequence.
    /// The .first is the previous location in sequence, the .second is the next in sequence.
    /// This size of this list is <= the size of the observation list (normally ==).
    std::list<placePtrPair_t> placesPrevNext;
    /// \brief Pointers to other locations connected to this one, but not merged.
    std::list<Place*> placesConnected;
    /// Pointers can only be managed by the PlaceGraph class
    friend class PlaceGraph;
};

///
/// \brief Stores the graph of all the places we've been in FabMap's topological map.
///
class CV_EXPORTS PlaceGraph
{
public:

    /// Default Constructor
    PlaceGraph()
    {}

    ///
    /// \brief Add a new place to the end of the visit list.
    /// \return The new place that has been created and visited.
    ///
    Place * addNew()
    {
        // Add the new place
        places.push_back(Place());
        // Point forward
        CV_DbgAssert(lastPlace->placesPrevNext.empty() == false);
        CV_DbgAssert(lastPlace->placesPrevNext.front().second == NULL);
        lastPlace->placesPrevNext.front().second = &places.back();
        // Point backward
        places.back().placesPrevNext.clear();
        places.back().placesPrevNext.push_back(Place::placePtrPair_t(lastPlace, NULL));
        // Update last place
        lastPlace = &places.back();
    }

    ///
    /// \brief Visit a previously visited place (connect it to previous and next places)
    /// \param place The place to be visited.
    ///
    void visitPlace(Place * place)
    {
        // Point forward
        CV_DbgAssert(lastPlace->placesPrevNext.empty() == false);
        CV_DbgAssert(lastPlace->placesPrevNext.front().second == NULL);
        lastPlace->placesPrevNext.front().second = place;
        // Point backward
        place->placesPrevNext.push_back(Place::placePtrPair_t(lastPlace, NULL));
        // Update last place
        lastPlace = place;
    }

    /// \brief Remove a place from the map.
    /// This location will no longer be compared against, and all references are removed.
    /// Not currently implemented
    void removePlace(Place * place);

    /// \brief Add connection
    void connectPlaces(Place * placeA, Place * placeB)
    {
        placeA->placesConnected.push_back(placeB);
        placeB->placesConnected.push_back(placeA);
    }

    /// Accessor for a read-only reference to the places
    const std::list<Place> & getPlaces() const
    {
        return places;
    }

protected:

    /// All of the places we've been, ordered by first visit
    std::list<Place> places;
    /// The most recently visited place.
    Place* lastPlace;
};

class CV_EXPORTS InferBase
{
public:
    InferBase (cv::Ptr<Mat> _clTreePtr,
               const double & _PzGe,
               const double & _PzGNe,
               const bool & _naiveBayes)
        : PzGe(_PzGe), PzGNe(_PzGNe),
          clTreePtr(_clTreePtr), clTree(*clTreePtr), naiveBayes(_naiveBayes)
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
    const Mat & getClTree()
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

    /// Are we using naiveBayes?
    bool naiveBayes;
    /// A pointer to the Chow-Liu tree stored in the FabMap class
    cv::Ptr<Mat> clTreePtr;
    /// A reference to the pointer, for legacy access (remove this eventually)
    Mat & clTree;
    /// A graph of all the places we've been, and the generators present
    PlaceGraph placeGraph;

    //@{
    double PzGe;    ///< Probability of measurement given the generator is present
    double PzGNe;   ///< Probability of measurement given the generator is absent
    //@}
};

class CV_EXPORTS InferBinary : public InferBase
{
public:
    InferBinary(cv::Ptr<Mat> _clTree,
                const double & _PzGe,
                const double & _PzGNe,
                const bool & _naiveBayes)
        : InferBase(_clTree, _PzGe, _PzGNe, _naiveBayes)
    {
        if (naiveBayes) {
            PzGL = &InferBinary::PzqGL;
        } else {
            PzGL = &InferBinary::PzqGzpqL;
        }
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
    double (InferBinary::*PzGL)(int q, bool zq, bool zpq, bool Lzq,
                              const bool & newPlace);
    //@}
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
    FabMap(const Mat& clTree, double PzGe, double PzGNe, int flags,
           int numSamples = 0);
    virtual ~FabMap();

    //@{
    /// Method to add training data for sampling method.
    virtual void addTraining(const Mat& queryImgDescriptor);
    virtual void addTraining(const std::vector<Mat>& queryImgDescriptors);
    //@}

    //@{
    /// Method to add to the test data (the map).
    virtual void add(const Mat& queryImgDescriptor);
    virtual void add(const std::vector<Mat>& queryImgDescriptors);
    //@}

    //@{
    /// Access the descriptors with a read-only reference
    const std::vector<Mat>& getTrainingImgDescriptors() const;
    const std::vector<Mat>& getTestImgDescriptors() const;
    //@}

    //@{ Image comparisons
    ///
    /// \brief FabMap image comparison (not full localization).
    /// \param queryImgDescriptor The descriptor for the query image
    /// \param testImgDescriptors The descriptors for the test images
    /// \param[out] matches Contains the match probabilities for the comparisons
    /// \param mask Not currently implemented?
    ///
    void compare(const Mat& queryImgDescriptor,
                 const Mat& testImgDescriptors, std::vector<IMatch>& matches,
                 const Mat& mask = Mat());
    ///
    /// \brief FabMap image comparison (not full localization).
    /// \param queryImgDescriptor The descriptor for the query image
    /// \param testImgDescriptors The descriptors for the test images
    /// \param[out] matches Contains the match probabilities for the comparisons
    /// \param mask Not currently implemented?
    ///
    void compare(const Mat& queryImgDescriptor,
                 const std::vector<Mat>& testImgDescriptors,
                 std::vector<IMatch>& matches, const Mat& mask = Mat());
    ///
    /// \brief FabMap image comparison (not full localization).
    /// \param queryImgDescriptors The descriptors for the query images
    /// \param testImgDescriptors The descriptors for the test images
    /// \param[out] matches Contains the match probabilities for the comparisons
    /// \param mask Not currently implemented?
    ///
    void compare(const std::vector<Mat>& queryImgDescriptors,
                 const std::vector<Mat>& testImgDescriptors,
                 std::vector<IMatch>& matches, const Mat& mask = Mat());
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
    DEPRECATED (void compare(const std::vector<Mat>& queryImgDescriptors, std::vector<
                             IMatch>& matches, bool addQuery = false, const Mat& mask =
            Mat()));
    ///
    /// \brief FabMap localization against the map.
    /// \param queryImgDescriptors The descriptors for the query images
    /// \param[out] matches Contains the match probabilities for the comparisons
    /// \param addQuery Add the query to the test descriptors when done
    /// \param mask Not currently implemented?
    /// TODO: Add addThreshold that replaces addQuery, and addIfLatest that doesn't add to the latest place.
    ///
    virtual void localize(const std::vector<Mat>& queryImgDescriptors, std::vector<
                          IMatch>& matches, bool addQuery = false, const Mat& mask =
            Mat());
    ///
    /// \brief FabMap localization against the map (deprecated, see FabMap::localize).
    /// \param queryImgDescriptor The descriptor for the query image
    /// \param[out] matches Contains the match probabilities for the comparisons
    /// \param addQuery Add the query to the test descriptors when done
    /// \param mask Not currently implemented?
    /// \deprecated This function has been deprecated, use FabMap::localize instead
    ///
    DEPRECATED (void compare(const Mat& queryImgDescriptor,
                             std::vector<IMatch>& matches, bool addQuery = false,
                             const Mat& mask = Mat()));
    ///
    /// \brief FabMap localization against the map.
    /// \param queryImgDescriptor The descriptor for the query image
    /// \param[out] matches Contains the match probabilities for the comparisons
    /// \param addQuery Add the query to the test descriptors when done
    /// \param mask Not currently implemented?
    /// TODO: Add addThreshold that replaces addQuery, and addIfLatest that doesn't add to the latest place.
    ///
    virtual void localize(const Mat& queryImgDescriptor,
                          std::vector<IMatch>& matches, bool addQuery = false,
                          const Mat& mask = Mat());
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
    void compareImgDescriptor(const Mat& queryImgDescriptor,
                              int queryIndex, const std::vector<Mat>& testImgDescriptors,
                              std::vector<IMatch>& matches);
    ///
    /// \brief Adds an image descriptor to the test descriptors (the map).
    /// Does no comparisons to the test data.
    /// \param queryImgDescriptor The query image descriptor that will be added.
    ///
    void addImgDescriptor(const Mat& queryImgDescriptor);
    //@}

    //@{ Likelihood functions, these change with each type of FabMap
    ///
    /// \brief The getLikelihoods method is overwritten for each different FabMap
    /// \param queryImgDescriptor The descriptor to be compared.
    /// \param testImgDescriptors The descriptors to compare against.
    /// \param[out] matches Result of the comparison of query to test.
    ///
    virtual void getLikelihoods(const Mat& queryImgDescriptor,
                                const std::vector<Mat>& testImgDescriptors,
                                std::vector<IMatch>& matches) = 0;
    ///
    /// \brief Calculates the likelihood the query comes from a new place.
    /// \param queryImgDescriptor The descriptor from the query image.
    /// \return The log-likelihood the query descriptor came from a new place.
    ///
    virtual double getNewPlaceLikelihood(const Mat& queryImgDescriptor);
    //@}

    ///
    /// \brief Turns measurement likelihoods into location probabilities.
    /// Also applies the motion model if in use (specified in the options flag).
    /// \param[in,out] matches Contains the input likelihoods, and output probabilities.
    ///
    void normaliseDistribution(std::vector<IMatch>& matches);

    //@{ Data
    /// Image descriptors seen in training (for sampled new location probability)
    std::vector<Mat> trainingImgDescriptors;
    /// Image descriptors seen so far
    std::vector<Mat> testImgDescriptors;
    /// Prior match probabilities for motion model p(L^k|Z^k-1)
    std::vector<IMatch> priormatches;
    /// Generator states for all locations
    std::vector<Mat> peGL;
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
    cv::Ptr<Mat> clTree;
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
    FabMap1(const Mat& clTree, double PzGe, double PzGNe, int flags,
            int numSamples = 0);
    virtual ~FabMap1();
protected:

    ///
    /// \brief FabMap1 implementation of likelihood comparison
    /// \param queryImgDescriptor The query descriptor to be compared.
    /// \param testImgDescriptors The reference descriptors to compare against.
    /// \param[out] matches The results of the comparisons in IMatch::likelihood.
    ///
    void getLikelihoods(const Mat& queryImgDescriptor,
                        const std::vector<Mat>& testImgDescriptors,
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
    FabMapLUT(const Mat& clTree, double PzGe, double PzGNe,
              int flags, int numSamples = 0, int precision = 6);
    virtual ~FabMapLUT();
protected:

    //FabMap look-up-table implementation of the likelihood comparison
    void getLikelihoods(const Mat& queryImgDescriptor, const std::vector<
                        Mat>& testImgDescriptors, std::vector<IMatch>& matches);

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
    FabMapFBO(const Mat& clTree, double PzGe, double PzGNe, int flags,
              int numSamples = 0, double rejectionThreshold = 1e-8, double PsGd =
            1e-8, int bisectionStart = 512, int bisectionIts = 9);
    virtual ~FabMapFBO();

protected:

    //FabMap Fast Bail-out implementation of the likelihood comparison
    void getLikelihoods(const Mat& queryImgDescriptor, const std::vector<
                        Mat>& testImgDescriptors, std::vector<IMatch>& matches);

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
    void setWordStatistics(const Mat& queryImgDescriptor, multiset<WordStats>& wordData);
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

    FabMap2(const Mat& clTree, double PzGe, double PzGNe, int flags);
    virtual ~FabMap2();

    //FabMap2 builds the inverted index and requires an additional training/test
    //add function
    void addTraining(const Mat& queryImgDescriptors) {
        FabMap::addTraining(queryImgDescriptors);
    }
    void addTraining(const std::vector<Mat>& queryImgDescriptors);

    void add(const Mat& queryImgDescriptors) {
        FabMap::add(queryImgDescriptors);
    }
    void add(const std::vector<Mat>& queryImgDescriptors);

protected:

    //FabMap2 implementation of the likelihood comparison
    void getLikelihoods(const Mat& queryImgDescriptor, const std::vector<
                        Mat>& testImgDescriptors, std::vector<IMatch>& matches);
    double getNewPlaceLikelihood(const Mat& queryImgDescriptor);

    //the likelihood function using the inverted index
    void getIndexLikelihoods(const Mat& queryImgDescriptor, std::vector<
                             double>& defaults, map<int, std::vector<int> >& invertedMap,
                             std::vector<IMatch>& matches);
    void addToIndex(const Mat& queryImgDescriptor,
                    std::vector<double>& defaults,
                    map<int, std::vector<int> >& invertedMap);

    //data
    std::vector<double> d1, d2, d3, d4;
    std::vector<std::vector<int> > children;

    // TODO: inverted map a std::vector?

    std::vector<double> trainingDefaults;
    map<int, std::vector<int> > trainingInvertedMap;

    std::vector<double> testDefaults;
    map<int, std::vector<int> > testInvertedMap;

};

///
/// \brief A Chow-Liu tree implementation designed for FAB-MAP.
///
/// A Chow-Liu tree is required by FAB-MAP. The Chow-Liu tree provides an estiMate of the
/// full distribution of visual words using a minimum spanning tree. The tree is
/// generated through training data.
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
    void add(const Mat& imgDescriptor);
    ///
    /// \brief You add data to the chow-liu tree before calling make.
    /// \param imgDescriptors A std::vector of \#imgs x \#words bag of words descriptors.
    ///
    void add(const std::vector<Mat>& imgDescriptors);
    //@}

    const std::vector<Mat>& getImgDescriptors() const;

    Mat make(double infoThreshold = 0.0);

private:
    std::vector<Mat> imgDescriptors;
    Mat mergedImgDescriptors;

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
    void createBaseEdges(list<info>& edges, double infoThreshold);
    double calcMutInfo(int word1, int word2);
    static bool sortInfoScores(const info& first, const info& second);

    //selecting minimum spanning egdges with maximum inforMation
    bool reduceEdgesToMinSpan(list<info>& edges);

    //building the tree sctructure
    Mat buildTree(int root_word, list<info> &edges);
    void recAddToTree(Mat &cltree, int q, int pq,
                      list<info> &remaining_edges);
    std::vector<int> extractChildren(list<info> &remaining_edges, int q);

};

///
/// \brief A custom vocabulary training method based on:
/// http://www.springerlink.com/content/d1h6j8x552532003/.
///
class CV_EXPORTS BOWMSCTrainer: public cv::BOWTrainer {
public:
    BOWMSCTrainer(double clusterSize = 0.4);
    virtual ~BOWMSCTrainer();

    // Returns trained vocabulary (i.e. cluster centers).
    virtual Mat cluster() const;
    virtual Mat cluster(const Mat& descriptors) const;

protected:

    double clusterSize;

};

} // namespace of2

#endif /* OPENFABMAP_H_ */

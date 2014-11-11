// Main Page Doxygen Documentation

///
/// \file
///
/// \brief TODO: Move this to README.md and translate to markdown
///

namespace cv
{

namespace of2
{

/*! \mainpage OpenFABMAP
 *
 *
 * This is an open and modifiable code-source which implements the Fast Appearance-based
 * Mapping algorithm (FAB-MAP) originally developed by Mark Cummins and Paul Newman.
 * OpenFABMAP was designed from published FAB-MAP theory and is for personal and research
 * use.
 *
 * FAB-MAP is a Simultaneous Localisation and Mapping algorithm which operates in
 * appearance space only. FAB-MAP performs location matching between places that have
 * been visited within the world as well as providing a measure of the probability of
 * being at a new, previously unvisited location. Camera images form the sole input to
 * the system, from which bag-of-words models are formed through the extraction of
 * appearance-based (e.g. SURF) features.
 *
 * The code has implementations of
 * * Feature Detection and Extraction and Bag-of-words models using OpenCV
 * * BOWMSCTrainer, a Bag-of-Words vocabulary trainer
 *   (<a href="http://www.springerlink.com/content/d1h6j8x552532003/">Teynor & Burkhardt 2007</a>)
 * * ChowLiuTree, a Chow-Liu tree implementation
 *   (<a href="http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=1054142">Chow & Liu 1968</a>)
 * * FabMap1, the original FabMap algorithm
 *   (<a href="http://ijr.sagepub.com/content/27/6/647.short">Cummins & Newman 2008</a>)
 * * FabMapLUT which uses a look-up table to precompute commonly used calculations for FabMap
 * * FabMapFBO which uses the fast bail-out speed-up for FabMap
 *   (<a href="http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=5613942">Cummins & Newman 2010</a>)
 * * FabMap2 which is able to handle much larger environments than the earlier FabMap algorithms
 *   (<a href="http://ijr.sagepub.com/content/30/9/1100.short">Cummins & Newman 2010</a>)
 *
 * For an overview of OpenFABMAP see
 *   (<a href="http://eprints.qut.edu.au/50317/1/glover_ICRA2012_final.pdf">Glover et al. 2012</a>).
 * OpenFABMAP was first used in
 *   (<a href="http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=5509547&tag=1">Glover et al. 2010</a>).
 *
 * As of the latest version, openFABMAP is dependent solely on <a
 * href="http://code.opencv.org/projects/opencv/wiki#OpenCV-Wiki">OpenCV2.3</a> or
 * higher. The project is designed to integrate with OpenCV 2.3 2D feature-based methods
 * and storage methods. The project has a <a href="http://www.cmake.org/">CMake</a> build
 * environment for general use on both Linux and Windows systems. See the README file for
 * more information on compiling the code.
 *
 * OpenFABMAP is also designed to integrate with <a href="http://wiki.ros.org/">Robot
 * Operating System (ROS)</a>. See the <a
 * href="https://wiki.qut.edu.au/display/cyphy/cyphy+ROS+wiki+page">CyPhy-ROS</a> page
 * for a package that has implemented openFABMAP as a ROS node.
 *
 * For questions on how to modify the source to your specific implementation, bug
 * reporting, comments and suggestions, or if you would like to become involved in
 * developing the openFABMAP project beyond the current implementation, contact via
 * <a href="https://github.com/arrenglover/openfabmap">github</a>.
 *
 * Citations <a href="openFABMAP.enw" target="_blank"><b>Endnote</b></a>
 *   <a href="openFABMAP.bib" target="_blank"><b>Bibtex</b></a>
 *
 * \image html surf_small.jpg "SURF features used by FAB-MAP"
 *
 */

}

}

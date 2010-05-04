/*
 * HashClustering.h
 *
 *  Created on: 26.04.2010
 *      Author: steffen
 */

#ifndef HASHCLUSTERING_H_
#define HASHCLUSTERING_H_

#include <OpenMS/DATASTRUCTURES/HashGrid.h>
#include <OpenMS/DATASTRUCTURES/DataSubset.h>
#include <OpenMS/COMPARISON/CLUSTERING/ClusteringMethod.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>

namespace OpenMS
{

/**
		@brief Hierarchical clustering based on geometric hashing. Only elemens in the surrounding are clustered.
		So it is not more necessary to calculate all distances which leads to a linear runtime instead of quadratic.

		@see ClusterHierarchical
		@ingroup SpectraClustering
	*/

class HashClustering : public ProgressLogger{
private:
	/**
	 * @brief current minimal distance
	 */
	DoubleReal min_distance;
	/**
	 * @brief two DataSubsets with current minimal distance
	 */
	std::pair<DataSubset*,DataSubset*> min_distance_subsets;
	/**
	 * @brief set of distances
	 */
	DistanceSet distances;
	/**
	 * @brief method which is used to calculate the distances
	 */
	ClusteringMethod* method;
	/**
	 * @brief the grid for geometric hashing
	 */
	HashGrid grid;
	/**
	 * @brief Calculates initial distances
	 */
	void init();
	/**
	 * @brief Merges two DataSubsets
	 */
	void merge();
	/**
	 * @brief Finds the two DataSubsets with minimal distance
	 * @param subset1 first DataSubset
	 * @param subset2 second DataSubset
	 */
	void updateMinElements();
	/**
	* @brief Calculates the distance of two DataSubsets using <i>getDistance</i> of the clustering method
	* @param subset1 first DataSubset
	 * @param subset2 second DataSubset
	 */
	DoubleReal getDistance(DataSubset& subset1,DataSubset& subset2);
	/**
	 * @brief Calculates the distance of two DataPoints using <i>getDistance</i> of the clustering method
	 * @param point1 first DataPoint
	 * @param point2 second DataPoint
	 */
	DoubleReal getDistance(DataPoint& point1,DataPoint& point2);
public:
	/**
	 * @brief Detailed constructor
	 * @param data this data points will be clustered
	 * @param rt_threshold height of the grid cells
	 * @param mz_threshold width of the grid cells
	 * @param method_ method to use for calculating distances
	 */
	HashClustering(std::vector<DataPoint>& data, int rt_threshold, int mz_threshold, ClusteringMethod& method_);
	/**
	 * @brief Calculates the silhouette values for any possible cluster number
	 * @param tree hierarchical clustering tree
	 */
	std::vector< Real > averageSilhouetteWidth(std::vector<BinaryTreeNode>& tree);
	/**
			@brief Method to calculate a partition resulting from a certain step in clustering given by the number of clusters

			@param cluster_quantity Size giving the number of clusters
			@param tree vector of BinaryTreeNodes representing the clustering
			@param clusters vector of vectors holding the clusters
			@see BinaryTreeNode

			after call of this method the argument clusters is filled corresponding to the given @p cluster_quantity with the indices of the elements clustered
	 */
	void cut(int cluster_quantity, std::vector< std::vector<DataPoint*> >& clusters, std::vector<BinaryTreeNode>& tree);

	/**
	 * @brief Starts the clustering and returns a vector of subtrees when finished
	 * @param subtrees vector of subtrees, which will be filled after the clustering process
	 */
	void performClustering(std::vector<std::vector<BinaryTreeNode > >& subtrees );
};
}



#endif /* HASHCLUSTERING_H_ */

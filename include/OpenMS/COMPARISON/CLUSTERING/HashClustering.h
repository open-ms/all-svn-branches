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

class HashClustering : public ProgressLogger{
private:
	ClusteringMethod* method;
	HashGrid grid;
	DoubleReal min_distance;
	std::pair<DataSubset*,DataSubset*> min_distance_subsets;
	void init();
	void merge(DataSubset& subset1,DataSubset& subset2);
	void updateMinElements();
public:
	HashClustering(std::vector<DataPoint>& data, int rt_threshold, int mz_threshold, ClusteringMethod& method_);
	std::vector< Real > averageSilhouetteWidth(std::vector<BinaryTreeNode>& tree);
	void cut(int cluster_quantity, std::vector< std::vector<DataPoint*> >& clusters, std::vector<BinaryTreeNode>& tree);
	void cut(int cluster_quantity, std::vector< std::vector<BinaryTreeNode> >& subtrees, std::vector<BinaryTreeNode>& tree);
	DoubleReal getDistance(DataSubset& subset1,DataSubset& subset2);
	DoubleReal getDistance(DataPoint& point1,DataPoint& point2);
	std::vector<std::vector<BinaryTreeNode> > performClustering();
};
}



#endif /* HASHCLUSTERING_H_ */
